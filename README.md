# Aquatic Optics Monte Carlo (AOMC) Model

[![Documentation](https://img.shields.io/website?label=documentation&style=for-the-badge&up_message=online&url=https%3A%2F%2Fmgimond.github.io%2Faomc)](https://mgimond.github.io/aomc)

AOMC is an open-source model written in FORTRAN 90 to simulate the propagation of light in an optically shallow, vertically heterogeneous aquatic medium.

This repository contains the source code for the model. For complete documentation, including theory, detailed compilation instructions, and input file formats, please visit the **[full documentation website](https://mgimond.github.io/aomc)**.

## Quick Start

1.  **Clone the repository:**
    ```sh
    git clone https://github.com/mgimond/aomc.git
    cd aomc
    ```

2.  **Compile the model:**
    Ensure you have a Fortran compiler like `gfortran` installed, then run the build script. This will create an executable named `aomc` in the project's root directory.
    ```sh
    ./build.sh
    ```

3.  **Run a simulation:**
    The default input files are located in the `input/` directory. To run the model with these files:
    ```sh
    cd input/
    ../aomc
    ```
    Output files (`*.out`) will be generated in the `input/` directory.

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.