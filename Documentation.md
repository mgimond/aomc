## Description of each FORTRAN files:

+ `geom2.f90`  Subroutine computes the new X,Y and Z position of the photon after a scattering event.
+ `global.f90` Subroutine contains all the global variables used in the AOMC model.
+ `interface.f90` Determines the outcome of the photon's interaction with the air-water interface.
+ `light_internal.f90` Subroutine computes the incident angle of the new photon on the water surface.
+ `logbin.f90` Subroutine logs the photon's trajectory.
+ `mc.f90` (main program) This is the main program.
+ `modules.f90` Contains the random number generator function "rand", and the angle summation routine "anglesum".
+ `rand_global.f90` Module is contained in subroutine "modules.f90" and computes the random number "rand".
+ `vsf.f90` Subroutine determines the photon's angle of scatter using the scattering phase function.
+ `water.f90` Subroutine traces the path of a photon from the moment it enters the water to the moment when it lost from the system.

## The input files:
+ `amc.inp`:  Main input file. Defines general environmental parameters such as bottom depth, vertical heterogenaty of water body, and solar zenith angle. Also defines model parameters such as sky model to use, number of photons to run and number of bands to simulate.
+ `conc.inp`: Defines the concentration for each parameters simulated at the different depths (if a vertically heterogenous water body is simulated). This file also allows the user to define the refractive index of each layer.
+ `lamba.inp`
+ `lambs.inp`
+ `lambda.inp`
+ `spf.inp`

## The output files:
The AOMC model will output 4  " *.out "  files:
+ `aop.out`  A listing by depth and wavelength of all of the apparent optical properties except radiance.
+ `iop.out`   A listing of the inherent optical properties.
+ `rad.out`  A listing by depth and wavelength of radiance at various angles
+ `wave.out` A listing by wavelength of the "just below the surface (0-)" values for several apparent optical properties (this is a subset of the data presented in the aop.out file).

## `amc.inp` input file parameters

1) Number of photons (iterations) to simulate. <br>
   UNITS  = n/a <br>
   FORMAT = I12 <br>
   LIMIT  = [0,99999999] <br>


2) Use boundary light condition from a lookup table (1) 
   or let the model generate it from a direct or uniformly
   diffuse source (1). If data are read from file, skylight
   distribution must wavelength dependent. <br>
   UNITS  = n/a <br>
   FORMAT = I1 <br>
   LIMIT  = [0;1] <br>
   
   a) Azimuth angle of specular light source (if using
      model's skylight distribution module) <br>
      UNITS  = deg <br>
      FORMAT = I8 <br>
      LIMIT  = [0,359] <br>

   b) Zenith angle of specular light source. 0deg is the
      polward angle. (if using model's skylight 
	distribution module) <br>
      UNITS  = deg <br>
      FORMAT = I8 <br>
      LIMIT  = [0,90[ <br>

3) Number of water constituents (including water).
   The minimum value should be 1 to include water. <br>
   UNITS  = n/a <br>
   FORMAT = I8 <br>
   LIMIT  = [1,20] <br>

4) Depth of bottom boundary.  <br>
   UNITS  = meter <br>
   FORMAT = F12.2 <br>
   LIMIT  = ]0,10000[ <br>

5) Width of side boundary (same for both X and Y). If
   a photon interacts with the side boundary, it is lost 
   from the system. <br>
   UNITS  = meter <br>
   FORMAT = F12.2 <br>
   LIMIT  = ]0,100000[ <br>

6) Number of layers in which apparent optical 
   properties are to be recorded.  <br>
   UNITS  = n/a <br>
   FORMAT = I8 <br>
   LIMIT  = [2,500] <br>
   
   a) Specify the thickness of the recording layers. If the
      number of layers are to span the entire water column,
	set this value to -1.0 ( The actual depth of these
      recordings will then be calculated from the ratio 
      (depth)/(depth intervals)) <br>
	    UNITS  = meter <br>
      FORMAT = F12.4 <br>
      LIMIT  = [-1] U ]0,+inf[ <br>

7) When the model records a photon's transfer across a layer 
   interval, it bins the angle of the photon's vector. The
   bins are divided into equal degrees of latitude. <br>
   UNITS  = n/a <br>
   FORMAT = I8 <br>
   LIMIT  = [1,180] <br>

8) When the model records a photon's transfer across a depth 
   interval, it bins the angle of the photon's vector. The
   bins are divided into equal degrees of longitude.
   UNITS  = n/a
   FORMAT = I8
   LIMIT  = [1,360]
    
   a)The user may opt to have the polward bin concatenated into
     a single polar bin (polar cap). A value of 1 will create a
     single polar bin. A value of 0 will break down the polar cap
     into PHI number of bins. <br>
     UNITS  = n/a <br>
     FORMAT = I1 <br>
     LIMIT  = [0;1] <br>

9) The user may request to log all instances of absorption
   within the water column. The log records the position
   X, Y, Z of the absorbing event, the total optical pathlength
   of the photon and it's last angle of propagation.
   This file is only applicable to single wavelength mode.
   (1=yes 0=no).   <br>
   UNITS  = n/a <br>
   FORMAT = I1 <br>
   LIMIT  = [0;1]     <br>

10) The user may request to log all instances of photons
    absorbed by the bottom boundary. The log records the 
    position X, Y, Z of the absorbing event, the total optical 
    pathlength of the photon and it's last angle of propagation.
    This file is only applicable to single wavelength mode.
    (1=yes 0=no).   <br>
    UNITS  = n/a <br>
    FORMAT = I1 <br>
    LIMIT  = [0;1] <br>


11) The user may request to log all instances of photons
    absorbed by the side boundary. The log records the 
    position X, Y, Z of the absorbing event, the total optical 
    pathlength of the photon and it's last angle of propagation.
    This file is only applicable to single wavelength mode.
    (1=yes 0=no).   <br>
    UNITS  = n/a <br>
    FORMAT = I1 <br>
    LIMIT  = [0;1] <br>

13) Determine if the cumulative SPF is to be read from a file (1) 
    or randomly determined from an isotropic function. <br>
    UNITS  = n/a <br>
    FORMAT = I1 <br>
    LIMIT  = [0;1] <br>
 
14) Number of layers with different optical characteristics
    (varying concentrations and refractive index) <br>
    UNITS  = n/a <br>
    FORMAT = I8 <br>
    LIMIT  = [1,100] <br>

15) Grid file output: select the number of rows and columns for 
    the grid output. If a value of 0 is selected, no grid file 
    will be created. <br>
    UNITS  = n/a <br>
    FORMAT = I8 <br>
    LIMIT  = [0,101] (must be odd number) <br>
    a) Size of each cell (same for dX and dY) <br>
       UNITS  = meters <br>
       FORMAT = F12.6 <br>
       LIMIT  = }0,inf[ <br>
    b) If file output is to be in SURFER (Golden Software) format, 
       only the upwelling radiance from the air/water surface is
       reported. The SURFER file format is a SURFER GRID ASCII
       file. <br>
       UNITS  = n/a <br>
       FORMAT = I8 <br>
       LIMIT  = [0;1] <br>

16) The light source type can be chosen to be a point (value=1),
   a circle (value=2) or a rectangle (value=3). <br>
   UNITS  = n/a <br>
   FORMAT = I8 <br>
   LIMIT  = [1;2;3] <br>
   
   + a) If a circular light source is chosen, enter the diameter of
     the collimated light source. <br>
     UNITS  = m <br>
     FORMAT = F12.6 <br>
     LIMIT  = ]0,+inf[ <br>
   + b) If a rectangular light source is chosen, enter the X
     value of its dimensions. <br>
     UNITS  = m <br>
     FORMAT = F12.6 <br>
     LIMIT  = ]0,+inf[ <br>
   + c) If a rectangular light source is chosen, enter the Y
     value of its dimensions. <br>
     UNITS  = m <br>
     FORMAT = F12.6 <br>
     LIMIT  = ]0,+inf[ <br>

16) A bottom target can be specified by the user <br>
   UNITS  = n/a <br>
   FORMAT = I8 <br>
   LIMIT  = [1,2] <br>
   a)The X and Y dimensions of the target need to be specified. <br>
     UNITS  = meter <br>
     FORMAT = F12.6 <br>
     LIMIT  = ]0,+inf[ <br>

17) Determine if optical properties are to be normalized to the
    incident radiation impinging on the water surface (0) or to that
    just below the air-water interface (1). <br>
     UNITS  = meter <br>
     FORMAT = F12.6 <br>
     LIMIT  = [0;1] <br>
