# Tabulated Sky Radiance Input (`skydist.inp`)

This document describes the format for the `skydist.inp` file, which allows the user to provide a custom, tabulated sky radiance distribution for the model's light source.

**Important**: This file is only read by the model if parameter **2) Light source from model(0) or file(1)** in `amc.inp` is set to `1`.

## File Format

The `skydist.inp` file has a multi-part structure consisting of a header, a list of wavelengths, and a data block of radiance values. The file is read in free-format, meaning all numerical values should be separated by spaces or commas.

### 1. Header Line

The very first line of the file must contain **4 integer values** that define the dimensions of the sky data:

1.  `numskyelem`: The total number of individual sky radiance data points that will be listed in the data block section of the file.
2.  `sky_zenith`: The number of zenith angle bins used in the sky distribution (i.e., the number of divisions from the horizon to the point directly overhead).
3.  `sky_azimuth`: The number of azimuthal angle bins. Note that the azimuth index in the data block will range from `0` to this value.
4.  `numskywave`: The number of distinct wavelengths for which sky radiance data is provided.

### 2. Wavelengths Line

The second line of the file must list the actual wavelength values. The number of values on this line **must** match the `numskywave` parameter from the header line.

### 3. Sky Radiance Data Block

Following the first two lines, the remainder of the file consists of the sky radiance data points. The number of lines in this section **must** match the `numskyelem` parameter from the header line.

Each line in this block represents a single sky "patch" or element and must contain **4 values**:

1.  `Wavelength Index`: An integer (1-based) that points to a wavelength in the Wavelengths Line. For example, a value of `1` corresponds to the first wavelength listed.
2.  `Zenith Index`: An integer (1-based) specifying the zenith angle bin, from `1` to `sky_zenith`.
3.  `Azimuth Index`: An integer (0-based) specifying the azimuth angle bin, from `0` to `sky_azimuth`.
4.  `Radiance Value`: A floating-point number representing the radiance intensity for that specific sky patch.

This format allows for a sparse representation of the sky; you only need to provide entries for the sky patches that have non-zero radiance.

### Internal Processing

After reading this file, the AOMC model processes the data for each wavelength by:
1.  Summing all radiance values for that wavelength.
2.  Normalizing the values so they represent a probability distribution.
3.  Multiplying this distribution by the total number of iterations (from `amc.inp`) to determine how many photons should be initiated from each specific sky patch.

### Wavelength Matching Logic

It is not necessary for the wavelengths listed in this file to exactly match the simulation wavebands defined in `lambbot.inp`. The model **does not interpolate** the sky data.

Instead, for each simulation waveband, the model performs a **nearest-neighbor search**. It identifies the single, numerically closest wavelength provided in `skydist.inp` and uses the complete sky radiance distribution from that single entry for the entire simulation run of that waveband.

For example, if `skydist.inp` contains data for 500nm and 600nm, and the simulation is running for wavebands at 520nm, 540nm, and 580nm:
*   The 520nm run will use the sky data for 500nm.
*   The 540nm run will also use the sky data for 500nm (since 540 is closer to 500 than to 600).
*   The 580nm run will use the sky data for 600nm.

### Note on Comments

The AOMC model expects input files to contain only the required numerical data in the correct order. 

## Example

Below is a valid `skydist.inp` file, followed by an explanation of its contents.

**File Content:**
```
6 2 2 2
550.0 650.0
1 1 0 10.5
1 1 1 12.0
1 2 1 20.0
2 1 0 8.2
2 1 1 9.5
2 2 1 15.3
```

**Explanation of the Example:**

*   **Header Line (`6 2 2 2`):**
    *   `6`: There are 6 radiance data lines to read in the data block.
    *   `2`: The sky is divided into 2 zenith bins (indexed 1 and 2).
    *   `2`: The sky is divided into 3 azimuth bins (indexed 0, 1, and 2).
    *   `2`: The file provides data for 2 different wavelengths.

*   **Wavelengths Line (`550.0 650.0`):**
    *   The first wavelength (`index=1`) is `550.0`.
    *   The second wavelength (`index=2`) is `650.0`.

*   **Sky Radiance Data Block:**
    *   The line `1 1 0 10.5` corresponds to: Wavelength Index 1, Zenith Index 1, Azimuth Index 0, with a radiance value of 10.5.
    *   The line `2 2 1 15.3` corresponds to: Wavelength Index 2, Zenith Index 2, Azimuth Index 1, with a radiance value of 15.3, and so on for all 6 lines.
