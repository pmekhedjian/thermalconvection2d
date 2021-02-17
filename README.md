# thermalconvection2d
Thermal Convection Code (in 2D)

All code rights and credit belong to Gary Glatzmaier at UC Santa Cruz
Website: https://websites.pmc.ucsc.edu/~glatz/

Author of:  Introduction to Modeling Convection in Planets and Stars: Magnetic Field, Density Stratification, Rotation
https://press.princeton.edu/books/hardcover/9780691141725/introduction-to-modeling-convection-in-planets-and-stars

Physics 227: Advanced Fluid Mechanics
https://websites.pmc.ucsc.edu/~glatz/phys227.html

-----------------------------------------------------------------------

Steps to run Fortran fluid code and produce IDL video:

0. Make whatever adjustments you need to in order to define your x,z grid. [See: parameter (nz=101,nn=50,nx=201) - nn is the number of modes in the perturbation, nz/nx is the level of granularity in the grid]

1. Compile either gfd_lin.f or gfd_nonlin.f using ifort or gfortran. The most interesting thermal convection phenomena occur in the nonlinear case.

2. Use IDL or GDL to run ka2.pro. You will need to .compile the included .pro files as they are pre-requisites to the procedures in ka2.pro. 

3. If you choose to produce a video from the resulting PPMs, make sure you have ffmpeg or mpeg_encode installed and in your PATH.
