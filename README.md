## What is the AOMC model? 

The AOMC model is an Open Source Aquatic Optics Monte Carlo model. It is written in FORTRAN 90. It simulates the propagation of light in an optically shallow, vertically heterogeneous aquatic medium. The model is in continuous development. Some of the routines have been extensively tested while others are still experimental. 

## What will the model simulate? 

The model simulates radiometric quantities (as measured by radiometers)  such as:
* Upwelling and downwelling scalar and planar irradiances (Euo, Edo, Eu and Ed)
* Radiance
* Mean cosines
* Shape factors (as used in the two-flow or two-stream model)
* Surface Albedo 


## What will the model NOT simulate? 

Currently, the model does not:
* simulate the influence of surface waves (however, this may be implemented in subsequent versions),
* allow for a horizontally variable bottom depth,
* allow for internal light source,
* simulate inelastic scattering (such as Raman scattering).

## What hardware is required to run this model? 

The model will run on any platform that supports a FORTRAN 90 compiler. Intel provides a somewhat free F90 compiler for Linux based PC  systems. This model was written using strict F90 standards. It has been successfully ported (1) to Windows 98 and Windows 2000 using Lahey and Compaq compilers, (2) to Linux OS using Intel’s non-commercial F90 compiler, (3) and SGI using MIPS F90 compiler. 

## Is output graphical? 

No. It is all console (command line) based. However, the output can easily be imported into a spreadsheet or a favorite graphing program.  

## Is the model free? 

Yes! It follows the Open Source Definition (see https://opensource.org/osd.html). 

## Can I modify the program? 

Of course! The development and improvement of the model is highly encouraged. All that is asked is that proper credit(s) be given where credit is due. 

## Is there support?

Yes and No. I haven't touched this project since 2003 so support is very limited.