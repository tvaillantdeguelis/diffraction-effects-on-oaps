#!/usr/bin/env python
# coding: utf8
"""
These scripts give a function ("compute_diffraction") which computes the 
diffraction image from an opaque planar object illuminated by a
monochromatic plane wave.

The function "plot_diffraction_image" allows to plot the diffracted
image and the binary one which would result from the 2D-S measurement.

Instrument configuration can be set in config.py

Particle shape and distance to the object plane are set in 
particle_shape_and_distance.py

Author: Thibault Vaillant de Guélis
Date: 2019/11/14

References:
- Vaillant de Guélis et al., Opt. Express, 27, 9372–9381, 
  doi:10.1364/OE.27.009372, 2019.
- Vaillant de Guélis et al., Atmos. Meas. Tech., 12, 2513–2529,
  doi:10.5194/amt-12-2513-2019, 2019.
"""


def load_config():

	l = 783 # (nm) wavelength
	pixel_size = 1.1 # (µm) image resolution on which diffraction pattern is computed
	x_2DSpixel_size = 11.4 # (µm) vertical axis (photodiode array axis)
	y_2DSpixel_size = 10.0 # (µm) horizontal axis (aircraft displacement axis)
	# max one decimals (ex: 11.4) for 2DS pixel sizes

	return l, pixel_size, x_2DSpixel_size, y_2DSpixel_size
