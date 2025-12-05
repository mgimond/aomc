# AMC Input Parameters (`amc.inp`)

This document describes the input parameters for the AOMC model, as defined in the `amc.inp` file. The `mc.f90` program reads this file, and is very specific about the data format and range of values for each parameter.

**Important Note for Float Parameters:** When providing input for parameters with a "Float" data type, you must include a decimal point, even if the value is a whole number (e.g., `100.0` instead of `100`).

The following table details each parameter, its data type, valid range, and any additional notes or restrictions.

| Input Parameter                               | Data Type | Range of Values          | Note                                                                                                                              |
| --------------------------------------------- | --------- | ------------------------ | --------------------------------------------------------------------------------------------------------------------------------- |
| 1) Iterations                                 | Integer   | `[0, 99999999]`          | Number of photons (iterations) to simulate.                                                                                       |
| 2) Light source from model(0) or file(1)      | Integer   | `[0, 1]`                 | Use boundary light condition from a lookup table (1) or let the model generate it (0).                                            |
| 2a) direct azimuth angle(deg)                 | Integer   | `[0, 359]`               | Azimuth angle of specular light source (if using model's skylight distribution).                                                  |
| 2b) direct zenith angle(deg)                  | Integer   | `[0, 89]`                | Zenith angle of specular light source. 0 deg is the poleward angle.                                                               |
| 3) number of system constituents              | Integer   | `[1, 20]`                | Number of water constituents (including water). Minimum value is 1.                                                              |
| 4) depth (M)                                  | Float     | `]0, 10000[`             | Depth of bottom boundary in meters.                                                                                               |
| 5) side boundaries (M)                        | Float     | `]0, 100000[`            | Width of side boundary (same for both X and Y) in meters.                                                                         |
| 6) Number of depth recording layers           | Integer   | `[2, 500]`               | Number of layers in which apparent optical properties are to be recorded.                                                         |
| 6a) Thickness of each recording layer       | Float     | `[-1] U ]0, +inf[`       | Thickness of recording layers in meters. If set to -1.0, layers will span the entire water column.                                  |
| 7) # of alpha intervals to bin (MAX=180)      | Integer   | `[1, 180]`               | Number of equal latitude bins for photon vector angle. Must be an even number.                                                     |
| 8) # of phi intervals to bin (MAX=180)        | Integer   | `[1, 360]`               | Number of equal longitude bins for photon vector angle.                                                                           |
| 8a) Use polar cap as grid?(yes=1, no=0)       | Integer   | `[0, 1]`                 | A value of 1 will create a single polar bin. A value of 0 will break down the polar cap into PHI number of bins.                   |
| 8b) Interval type (ang = 1,cos = 0)           | Integer   | `[0, 1]`                 | Defines interval type.                                                                                                            |
| 9) X,Y,Z for absorption (yes=1,no=0)          | Integer   | `[0, 1]`                 | Log all instances of absorption within the water column. Only for single wavelength mode.                                         |
| 10) BOT.-output(yes=1,no=0)                   | Integer   | `[0, 1]`                 | Log all instances of photons absorbed by the bottom boundary. Only for single wavelength mode.                                    |
| 11) SIDE-output(yes=1,no=0)                   | Integer   | `[0, 1]`                 | Log all instances of photons absorbed by the side boundary. Only for single wavelength mode.                                      |
| 12) SPF from file? (yes=1,no=0)               | Integer   | `[0, 1]`                 | Determine if the cumulative SPF is to be read from a file (1) or randomly determined from an isotropic function (0).               |
| 12a) Number of angles                       | Integer   | `[0, 360]`               | Number of scattering phase function angles to read in.                                                                            |
| 13) # of concentration dependent layers       | Integer   | `[1, 100]`               | Number of layers with different optical characteristics.                                                                          |
| 14) Number of cells in X-Y direction(0=NOGRID)| Integer   | `[0, 101]`               | Number of rows and columns for grid output. If 0, no grid file is created. Must be an odd number.                                   |
| 14a) Cell size (in unit meters)             | Float     | `]0, inf[`               | Size of each cell (same for dX and dY) in meters.                                                                                 |
| 14b) Output to Surfer format (yes=1,no=0)     | Integer   | `[0, 1]`                 | If 1, output is in Surfer GRID ASCII format.                                                                                      |
| 15) Light source: point[1],circle[2],rect.[3] | Integer   | `[1, 2, 3]`              | Type of light source.                                                                                                             |
| 15a) If circle, define diameter             | Float     | `]0, +inf[`              | Diameter of the collimated light source if circular.                                                                              |
| 15b) If rectangle, define dX                | Float     | `]0, +inf[`              | X dimension if rectangular light source.                                                                                          |
| 15c) and define dY                          | Float     | `]0, +inf[`              | Y dimension if rectangular light source.                                                                                          |
| 16) Bottom target: NO[0],YES[1]               | Integer   | `[0, 1]`                 | Specify a bottom target.                                                                                                          |
| 16a) define dX                              | Float     | `]0, +inf[`              | X dimension of the target.                                                                                                        |
| 16b) define dY                              | Float     | `]0, +inf[`              | Y dimension of the target.                                                                                                        |
| 17) Normalize output to air(0) or water(1)    | Integer   | `[0, 1]`                 | Normalize optical properties to incident radiation on water surface (0) or just below (1).                                        |
| 18) Number of wavebands (minimum = 1)         | Integer   | `[1, 500]`               | Number of wavebands to simulate.                                                                                                  |
