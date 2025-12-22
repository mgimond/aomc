#  Read AOP output file and return long-form data
#  Function pulled from the raomc package

make_long_aop <- function(workdir = ".", aop_file = "aop.out", lambbot_file = "lambbot.inp") {
  aop_path <- file.path(workdir, aop_file)
  lamb_path <- file.path(workdir, lambbot_file)
  
  if (!file.exists(aop_path)) stop("AOP file not found: ", aop_path)
  if (!file.exists(lamb_path)) stop("lambbot.inp file not found: ", lamb_path)
  
  # --- Read wavelengths from lambbot.inp ---
  lamb_lines <- readLines(lamb_path, warn = FALSE)
  lamb_lines <- lamb_lines[nzchar(trimws(lamb_lines))]
  if (length(lamb_lines) < 2) stop("lambbot.inp does not contain wavelength data.")
  # Drop first line (record count), parse remaining
  num_lines <- lamb_lines[-1]
  num_lines <- num_lines[grepl("^[0-9]", trimws(num_lines))]
  dat <- read.table(text = paste(num_lines, collapse = "\n"), header = FALSE)
  wavelengths <- as.numeric(dat[[1]])
  
  
  # --- Read aop.out ---
  lines <- readLines(aop_path, warn = FALSE)
  trim <- function(x) sub("^\\s+|\\s+$", "", x)
  split_ws <- function(s) strsplit(trim(s), "\\s+")[[1]]
  
  hdr_idx <- grep("^\\s*depth\\b", lines, ignore.case = TRUE)
  if (length(hdr_idx) == 0) stop("No 'depth (m)' header lines found in aop.out")
  
  block_ranges <- lapply(seq_along(hdr_idx), function(i) {
    start <- hdr_idx[i]
    end <- if (i < length(hdr_idx)) hdr_idx[i + 1] - 1 else length(lines)
    c(start = start, stop = end)
  })
  
  if (length(wavelengths) < length(block_ranges)) {
    warning("Fewer wavelengths than AOP blocks; recycling.")
    wavelengths <- rep(wavelengths, length.out = length(block_ranges))
  }
  
  out_list <- vector("list", length(block_ranges))
  
  for (b in seq_along(block_ranges)) {
    rng <- block_ranges[[b]]
    header_line <- lines[rng["start"]]
    hdr_tokens <- split_ws(header_line)
    hdr_tokens <- hdr_tokens[hdr_tokens != "(m)"]
    hdr_tokens[1] <- "depth"
    hdr_tokens <- gsub("^R\\*$", "R", hdr_tokens)
    
    rows <- lines[(rng["start"] + 1):rng["stop"]]
    rows <- rows[!grepl("^\\s*\\*", rows)]
    rows <- rows[nzchar(trim(rows))]
    if (!length(rows)) next
    
    tc <- textConnection(paste(rows, collapse = "\n"))
    df_block <- read.table(tc, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    close(tc)
    
    K <- length(hdr_tokens)
    if (ncol(df_block) != K) {
      if (ncol(df_block) > K) {
        df_block <- df_block[seq_len(K)]
      } else {
        for (j in (ncol(df_block) + 1):K) df_block[[j]] <- NA_real_
      }
    }
    names(df_block) <- hdr_tokens
    df_block[] <- lapply(df_block, function(x) suppressWarnings(as.numeric(x)))
    if (nrow(df_block) >= 1) df_block$depth[1] <- -1
    
    aop_cols <- setdiff(names(df_block), "depth")
    long <- reshape(df_block,
                    varying = aop_cols,
                    v.names = "value",
                    timevar = "aop",
                    times = aop_cols,
                    direction = "long")
    rownames(long) <- NULL
    long$wavelength <- wavelengths[b]
    long <- long[, c("wavelength", "depth", "aop", "value")]
    out_list[[b]] <- long
  }
  
  out <- do.call(rbind, out_list)
  out <- out[order(out$wavelength, out$depth, out$aop), ]
  rownames(out) <- NULL
  out
}
