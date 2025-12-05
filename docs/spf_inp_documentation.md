# Scattering Phase Function Input (`spf.inp`)

This document describes the format for the `spf.inp` file, which provides the cumulative scattering phase function (SPF) for each water constituent. This data determines the new direction of a photon after a scattering event.

## File Format

The structure of the `spf.inp` file is as follows:

1.  **Number of Angles**: The very first value must be an integer representing the number of angles in the phase function table. This number **must** match the value for parameter **12a) Number of angles** in `amc.inp`.

2.  **Cumulative SPF Table**: Following the first line, the file must contain one line for each angle. Each line is read in free-format (space- or comma-separated) and contains `N+1` values, where `N` is the number of constituents (from `amc.inp` parameter 3):
    *   **Column 1**: The scattering angle in **radians**.
    *   **Column 2 to N+1**: The cumulative scattering probability for each constituent, from 0.0 to 1.0. Each column corresponds to a constituent in the same order they are defined elsewhere in the model.

## How the Model Uses This File

When a photon is scattered, the model uses the data from `spf.inp` to determine its new direction in a two-step process:

1.  **Select the Scattering Particle**: First, the model randomly determines which type of particle (constituent) caused the scattering. This is based on the relative scattering contribution of each constituent at that wavelength.

2.  **Determine the New Angle**: Once a particle is chosen, the model uses its corresponding cumulative SPF column as a probability lookup table to find the new scattering angle.
    *   A random number between 0 and 1 is generated.
    *   The model finds where this random number would fit within the chosen constituent's cumulative SPF column, identifying the two consecutive probability values it falls between.
    *   **By default, the model then assigns the discrete angle associated with the *upper* boundary of this probability interval.** For example, if the random number falls between the probability for 0.1 radians and 0.2 radians, the new scattering angle will be set to 0.2 radians.

**Note on Interpolation**: The `vsf.f90` subroutine contains a commented-out (disabled) code block that would perform linear interpolation to calculate a more precise angle within the interval. However, this is not the default behavior. As a result, the slope values calculated when reading `spf.inp` are not used in the default configuration.

## Example

Here is a sample `spf.inp` file configured for **5 angles** and **2 system constituents**.

```
5
0.000   0.000   0.000
0.785   0.600   0.500
1.571   0.800   0.750
2.356   0.950   0.900
3.142   1.000   1.000

--- NOTES ---
This section is not read by the AOMC model.
Format: [Angle in Radians] [Constituent 1 Cum. SPF] [Constituent 2 Cum. SPF]
```
