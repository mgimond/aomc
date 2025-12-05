# Specific Scattering Coefficient Input Parameters (`scat.inp`)

This document describes the format for the `scat.inp` file, which defines the wavelength-dependent specific scattering coefficients for the various water quality constituents simulated by the AOMC model.

## File Format

The structure of the `scat.inp` file is strict and must conform to the following rules for the model to parse it correctly:

1.  **Number of Wavebands**: The very first value in the file must be an integer representing the total number of wavebands. This number **must** match the value specified for parameter **18) Number of wavebands** in the `amc.inp` file.

2.  **Scattering Coefficients**: Following the number of wavebands, the file must contain a matrix of specific scattering coefficient values.
    *   **Rows**: The number of rows in this matrix must match the number of wavebands declared on the first line.
    *   **Columns**: The number of columns must match the number of system constituents specified in parameter **3) number of system constituents** in the `amc.inp` file. Each column represents a constituent (e.g., pure water, phytoplankton, non-algal particles, CDOM).
    *   **Data Format**: The numerical values are read in free-format (list-directed). This means the numbers should be separated by spaces or commas. They do not need to be aligned in fixed columns, as long as each line contains the correct number of values.

3.  **Additional Information**: Any text or data placed after the specific scattering coefficient matrix is ignored by the model. This space can be used for comments, notes, or metadata, as shown in the example below.

## Example

Here is a sample `scat.inp` file configured for **12 wavebands** and **3 system constituents**.

```
12
0.020   0.50   0.10
0.019   0.48   0.10
0.018   0.45   0.10
0.015   0.40   0.10
0.012   0.35   0.09
0.009   0.30   0.08
0.007   0.25   0.07
0.005   0.20   0.06
0.003   0.15   0.05
0.002   0.10   0.04
0.001   0.08   0.03
0.001   0.06   0.02

--- NOTES ---
This section is not read by the AOMC model.
The columns represent the following constituents:
1. Pure Water
2. Phytoplankton
3. Non-Algal Particles
Data sourced from internal lab measurements, 2023.
```
