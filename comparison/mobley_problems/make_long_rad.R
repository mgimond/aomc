# Parse Monte Carlo radiance output into a long (tidy) table
# Script pulled from the raomc package

make_long_rad <- function(
    workdir = ".",
    rad_file = "rad.out",
    lambbot_file = "lambbot.inp",
    verbose = FALSE
) {
  rad_path  <- file.path(workdir, rad_file)
  lamb_path <- file.path(workdir, lambbot_file)
  
  # --- Read wavelength table (lambbot.inp) ---
  if (!file.exists(lamb_path)) stop("lambbot file not found: ", lamb_path)
  
  n_rows <- suppressWarnings(as.integer(readLines(lamb_path, n = 1)))
  if (is.na(n_rows) || n_rows <= 0) stop("First line of lambbot.inp must be a positive integer row count.")
  
  dat <- tryCatch({
    read.table(
      lamb_path,
      skip       = 1,
      nrows      = n_rows,
      col.names  = c("wavelength","frac1","frac2","frac3","frac4"),
      colClasses = rep("numeric", 5)
    )
  }, error = function(e) {
    stop("Failed to read lambbot.inp numeric block: ", conditionMessage(e))
  })
  
  wavelengths <- dat$wavelength
  
  if (verbose) {
    message("Read ", n_rows, " wavelength rows from lambbot.inp.")
  }
  
  # --- Read rad.out ---
  if (!file.exists(rad_path)) stop("rad.out not found: ", rad_path)
  lines <- readLines(rad_path)
  
  # Find the start of each wavelength block (Angle header line)
  angle_header_idx <- grep("^\\s*Angle", lines, ignore.case = TRUE)
  if (length(angle_header_idx) == 0) {
    warning("No 'Angle' header lines found in rad.out. Check header spelling/spacing.")
    return(data.frame(Wavelength = numeric(), depth = numeric(), alpha = numeric(), phi = numeric(), radiance = numeric()))
  }
  
  n_blocks <- length(angle_header_idx)
  n_use <- min(n_blocks, length(wavelengths))
  
  # Helpers ---------------------------------------------------------------
  
  # Extract all numeric tokens from a string
  extract_nums <- function(s) {
    as.numeric(unlist(regmatches(s, gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", s))))
  }
  
  # Resolve phi centers given header angles and data column count
  resolve_phi_centers <- function(phi_hdr, ncols, verbose, block_id) {
    n_hdr <- length(phi_hdr)
    
    if (ncols == n_hdr) {
      # Header lists centers
      list(phi_centers = phi_hdr)
    } else if (ncols == n_hdr - 1) {
      # Header lists edges -> compute centers
      phi_centers <- (phi_hdr[-n_hdr] + phi_hdr[-1]) / 2
      list(phi_centers = phi_centers)
    } else if (n_hdr >= (ncols + 1) && ncols == floor((n_hdr - 1) / 2)) {
      # Full-plane edges in header, half-plane bins in data
      edges_half  <- phi_hdr[seq_len(ncols + 1)]
      phi_centers <- (edges_half[-length(edges_half)] + edges_half[-1]) / 2
      list(phi_centers = phi_centers)
    } else if (ncols == ceiling(n_hdr / 2)) {
      # Header lists centers full-plane, output half-plane -> take first half
      phi_centers <- phi_hdr[seq_len(ncols)]
      list(phi_centers = phi_centers)
    } else {
      # Fallback: truncate
      list(phi_centers = phi_hdr[seq_len(min(n_hdr, ncols))])
    }
  }
  
  # ----------------------------------------------------------------------
  
  result <- vector("list", n_use)
  
  for (w in seq_len(n_use)) {
    start_idx <- angle_header_idx[w]
    end_idx   <- if (w < length(angle_header_idx)) angle_header_idx[w + 1] - 1 else length(lines)
    block     <- lines[start_idx:end_idx]
    
    if (length(block) < 2) next
    
    header_line <- block[1]
    
    # --- φ angles from header line ---
    # We parse the numbers exactly as they appear (degrees or radians)
    phi_hdr <- extract_nums(header_line)
    
    if (length(phi_hdr) == 0) {
      if (verbose) message(sprintf("[Block %d] No numeric φ angles found in header.", w))
      next
    }
    
    block <- block[-1] # Drop header
    
    # --- Depth lines ---
    # Regex matches lines ending in "number m", handling "Depth:" prefix
    depth_line_idx <- grep("[-+]?[0-9]*\\.?[0-9]+\\s*m\\s*$", block)
    if (length(depth_line_idx) == 0) next
    
    rows_out <- list()
    
    for (d in seq_along(depth_line_idx)) {
      d_start <- depth_line_idx[d]
      d_end   <- if (d < length(depth_line_idx)) depth_line_idx[d + 1] - 1 else length(block)
      
      depth_line <- block[d_start]
      depth_val  <- as.numeric(sub(".*?([-+]?[0-9]*\\.?[0-9]+)\\s*m\\s*$", "\\1", depth_line))
      
      # Sentinel: First depth is -1 (often "in air" or "just above surface")
      if (d == 1) depth_val <- -1
      
      depth_block <- block[(d_start + 1):d_end]
      if (length(depth_block) == 0) next
      
      # Alpha rows start with a number
      alpha_rows_idx <- grep("^\\s*[-+]?\\d", depth_block)
      if (length(alpha_rows_idx) == 0) next
      alpha_rows <- depth_block[alpha_rows_idx]
      
      # Determine φ centers using the first data row
      if (d == 1 || !exists("phi_centers")) {
        first_nums <- extract_nums(alpha_rows[1])
        if (length(first_nums) < 2) next
        n_phi_cols <- length(first_nums) - 1
        
        phi_res <- resolve_phi_centers(phi_hdr, n_phi_cols, verbose, w)
        phi_centers <- phi_res$phi_centers
      }
      
      # Parse each alpha row
      for (row in alpha_rows) {
        nums <- extract_nums(row)
        if (length(nums) < 2) next
        
        alpha_val <- nums[1] # Preserves unit from file
        radiances <- nums[-1]
        
        # Align radiances
        if (length(radiances) < length(phi_centers)) {
          radiances <- c(radiances, rep(NA_real_, length(phi_centers) - length(radiances)))
        } else if (length(radiances) > length(phi_centers)) {
          radiances <- radiances[seq_along(phi_centers)]
        }
        
        rows_out[[length(rows_out) + 1]] <- data.frame(
          Wavelength = rep(wavelengths[w], length(phi_centers)),
          depth      = rep(depth_val,     length(phi_centers)),
          alpha      = rep(alpha_val,     length(phi_centers)),
          phi        = phi_centers,
          radiance   = radiances,
          row.names  = NULL
        )
      }
    }
    
    result[[w]] <- if (length(rows_out)) do.call(rbind, rows_out) else NULL
  }
  
  out <- if (length(result)) do.call(rbind, result) else data.frame(
    Wavelength = numeric(), depth = numeric(), alpha = numeric(), phi = numeric(), radiance = numeric()
  )
  
  out <- out[!is.na(out$radiance), , drop = FALSE]
  
  if (verbose) {
    message("Parsed rows: ", nrow(out))
  }
  out
}
