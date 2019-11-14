# diffraction-effects-on-oaps

These scripts give a function ("compute_diffraction") which computes the 
diffraction image from an opaque planar object illuminated by a
monochromatic plane wave.

The function "plot_diffraction_image" allows to plot the diffracted
image and the binary one which would result from the 2D-S measurement. Return *diffraction_image.png*.

Instrument configuration can be set in config.py

Particle shape and distance to the object plane are set in 
particle_shape_and_distance.py

Inputs
----------
In particle_shape_and_distance.py:
	- Z: distance to object plane (in cm)
	- M: transmittance mask (particle shape) with size (y.size, x.size)
	
In config.py:
	- l: wavelength of the monochromatic incident light (in nm)
	- pixel_size: image resolution on which diffraction pattern is computed (in µm)
	- x_2DSpixel_size: vertical axis (photodiode array axis) resolution (in µm)
	- y_2DSpixel_size: horizontal axis (aircraft displacement axis) resolution (in µm)
	
Output
----------
	- I: intensity of the diffracted light on a screen orthogonal to the incident direction at a distance Z from the object
	- I_2DS: theoretical binary record by the 2D-S
  
References
----------
T. Vaillant de Guélis, V. Shcherbakov, A. Schwarzenböck, 2019: Diffraction patterns from opaque planar objects 
    simulated with Maggi–Rubinowicz method and angular spectrum theory, Opt. Express, 27, 9372–9381, 
    [doi:10.1364/OE.27.009372](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-27-7-9372).

T. Vaillant de Guélis, A. Schwarzenböck, V. Shcherbakov, C. Gourbeyre, B. Laurent, R. Dupuy, 
    P. Coutris, and C. Duroure, 2019: Study of the diffraction pattern of cloud particles and respective response 
    of Optical Array Probes, Atmos. Meas. Tech., 12, 2513–2529, 
    [doi:10.5194/amt-12-2513-2019](https://www.atmos-meas-tech.net/12/2513/2019/amt-12-2513-2019.html).
