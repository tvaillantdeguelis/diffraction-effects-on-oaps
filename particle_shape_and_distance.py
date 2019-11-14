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

import numpy as np


def load_shape_and_distance(x, y, X, Y):
	"""
	Set the distance Z of the particle from the object plane.
	Set the shape of the particle in the array M which has the same size
	and resolution than the image on which diffraction pattern will be 
	computed (0 where particle stops light, 1 where particle lets light 
	by).
	"""

	############################
	# DISTANCE Z FROM OBJECT PLANE
	Z = 2.35 # (cm) distance Z


	########################
	# PARTICLE SHAPE ARRAY M

	# Example with rectangle:

	# Parameters
	rect_w = 200 # (µm) rectangle width 
	rect_l = 100 # (µm) rectangle length
	
	# Create mask M of opaque rectangle
	M = np.ones((y.size, x.size))
	A = (np.abs(Y) <= (rect_w)/2.)
	B = (np.abs(X) <= (rect_l)/2.)
	M[A&B] = 0

	return Z, M