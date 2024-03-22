# glmt-hb
For computing scattering cross sections for gaussian beams incident on homogeneous spheres of different radii

Contents:

=== Programs ===
- init_glmt.o - initialises the job
- init_glmt.f90 - source code for the above
- spher_f_mono_HB.o - Mie code written by Mischenko - adapted to loop over multiple radii (radii contained in text file input.txt, which is generated automatically)
- spher_f_mono_HB.f - source code for the above
- glmt.o - generalised Lorenz-Mie thoery by Jianqi Shen et al. (see references in source code) - adapted to loop over multiple radii
- Light_scattering_calculation_of_an_Elliptical_Gaussian_Beam_HB.cpp - source code for the above
- final_glmt.0 - finalises the job and outputs scattering cross sections to file CEXT_GLMT
- final-glmt.f90 - source code for the above

=== Input Files ===
- Parameters for calculations.txt
You will likely not need to change this file very often. The code is very sensitive to poor formatting. 

The first line must be a single integer (0, 1, or 2) which indicates how the beam shape coefficients should be calculated. 0 is angular decomposition spectra, 1 is 2D quadrature of GLMT, and 2 is the localised approximation. Option 2 is MUCH faster than options 0 and 1 for large beam dimensions but the accuracy of this method has been shown to be unreliable in the 2018 paper cited within the glmt source code (see files list above). It should be mentioned, that the beam shape coefficients only need calculating once (described in further detail below), so using methods 0 and 1 may be preferred if the results appear to fit Mie data better for small radii.

The second line contains 4 integers (0 or 1) separated by tabs, which indicate the computations that should be performed. The first integer indicates whether the beam shape coefficients should be recalculated. If using angular spectra decomposition or 2D quadrature methods, you should set this to 0 if you have already calculated the beam shape coefficients, which speeds up the program runtime massively. In general, the second line should really only be set to 1 of 2 possible combinations:  
i) 1 0 1 0 - Recalculate beam shape coefficients and then compute far field scattering  
ii) 0 0 1 0 - Reuse pre-existing beam shape coefficients and then compute far field scattering  

- Parameters of incident beam.txt
This is where the beam shape and wavelength are defined. Dimensions given in microns.

The first line contains the beam waist in x direction, beam waist in y direction, and wavlength. These are all separated by tabs. 

The second line contains the beam offset, which has not been tested for any values other than 0, 0, 0. Note that, while it may be tempting to set the beam waist in one direction to a huge number, this will probably lead to extremely long beam shape coefficient computation times. Although as described above, beam shape coefficients only need to be calculated once.

- Parameters for particle.txt
This is where the particle sizes and refractive index are defined.

First line contains the start radius, end radius, and number of radius increments (same format as Matlab function linspace()). All separated by tabs. Dimensions are in microns.

Second line contains refractive index. Real part followed by imaginary part, separated by tab.

== Output Files ===
All output files are placed into a new directory named job#, where # is an integer that is determined automatically. Some other output files are jumped in the main directory because I was too lazy to put the somewhere else. The job directory contains the Mie scattering cross sections, the glmt phase function data (theta only), and the scattering cross sections from the glmt method.

**Note: After uploading to the cluster, you may receive an "access denied" error. I am not sure what causes this. A workaround is to recompile the source code:**
```
module load intel; ifort init_glmt.f90 -o init_glmt.o; ifort final_glmt.f90 -o final_glmt.o; ifort spher_f_mono_HB.f -o spher_f_mono_HB.o; g++ Light_scattering_calculation_of_an_Elliptical_Gaussian_Beam_HB.cpp -o glmt.o
```

