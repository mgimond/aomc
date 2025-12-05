# Specific Absorption Coefficient Input Parameters (`abs.inp`)

This document describes the format for the `abs.inp` file, which defines the wavelength-dependent specific absorption coefficients for the various water quality constituents simulated by the AOMC model.

## File Format

The structure of the `abs.inp` file is strict and must conform to the following rules for the model to parse it correctly:

1.  **Number of Wavebands**: The very first value in the file must be an integer representing the total number of wavebands. This number **must** match the value specified for parameter **18) Number of wavebands** in the `amc.inp` file.

2.  **Absorption Coefficients**: Following the number of wavebands, the file must contain a matrix of absorption coefficient values.
    *   **Rows**: The number of rows in this matrix must match the number of wavebands declared on the first line.
    *   **Columns**: The number of columns must match the number of system constituents specified in parameter **3) number of system constituents** in the `amc.inp` file. Each column represents a constituent (e.g., pure water, phytoplankton, non-algal particles, CDOM).
    *   **Data Format**: The numerical values are read in free-format (list-directed). This means the numbers should be separated by spaces or commas. They do not need to be aligned in fixed columns, as long as each line contains the correct number of values.

3.  **Additional Information**: Any text or data placed after the absorption coefficient matrix is ignored by the model. This space can be used for comments, notes, or metadata, as shown in the example below.

## Example

Here is a sample `abs.inp` file configured for **12 wavebands** and **3 system constituents**.

```
12
0.0049  0.032  0.004
0.0051  0.033  0.004
0.0058  0.035  0.004
0.0125  0.040  0.004
0.0210  0.046  0.004
0.0380  0.052  0.004
0.0510  0.058  0.004
0.0620  0.061  0.004
0.0810  0.063  0.004
0.1090  0.062  0.004
0.1550  0.059  0.004
0.2220  0.051  0.004

--- NOTES ---
This section is not read by the AOMC model.
The columns represent the following constituents:
1. Pure Water
2. Phytoplankton
3. Non-Algal Particles
Data sourced from internal lab measurements, 2023.
```
