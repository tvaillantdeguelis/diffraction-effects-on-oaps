# diffraction-effects-on-oaps

This script give a function ("compute_diffraction") which computes the 
diffraction image from an opaque planar object illuminated by a
monochromatic plane wave.

The function "plot_diffraction_image" allows to plot the diffracted
image and the binary one which would result from the 2D-S measurement.

The main program gives an example of the use of these functions for a
rectangular shape.

Inputs
----------
	- Z: distance to object plane (µm)
	- l: wavelength (µm)
	- pixel_size: (µm)
	- M: transmittance mask (particle shape) with size (y.size, x.size)
	- filter_fmax: (µm-1) low freq filter on FFT to remove noise (you can turn to False in the function if not needed)
  
Output
----------
	- I: intensity of the diffracted light on a screen orthogonal to the incident direction at a distance Z from the object.
  
References
----------
T. Vaillant de Guélis, V. Shcherbakov, A. Schwarzenböck, 2019: Diffraction patterns from opaque planar objects 
    simulated with Maggi–Rubinowicz method and angular spectrum theory, Opt. Express, 27, 9372–9381, 
    [doi:10.1364/OE.27.009372](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-27-7-9372).

T. Vaillant de Guélis, A. Schwarzenböck, V. Shcherbakov, C. Gourbeyre, B. Laurent, R. Dupuy, 
    P. Coutris, and C. Duroure, 2019: Study of the diffraction pattern of cloud particles and respective response 
    of Optical Array Probes, Atmos. Meas. Tech., 12, 2513–2529, 
    [doi:10.5194/amt-12-2513-2019](https://www.atmos-meas-tech.net/12/2513/2019/amt-12-2513-2019.html).
