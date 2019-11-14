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
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MultipleLocator, FixedLocator
import seaborn as sns
import matplotlib.gridspec as gridspec

from config import load_config
from particle_shape_and_distance import load_shape_and_distance


def compute_diffraction(Z, l, pixel_size, x, y, X, Y, M, filter_fmax,
						lowpass_filter=True):
	"""
	Arguments needed for the function are:
	- Z: distance to object plane (µm)
	- l: wavelength (µm)
	- pixel_size: (µm)
	- M: transmittance mask (particle shape) with size (y.size, x.size)
	- filter_fmax: (µm-1) low freq filter on FFT to remove noise (you 
	  can turn to False in the function if not needed)
	- x, y, X, Y: obtained with:
		x_min = -1024*pixel_size # (µm)
		x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
		y_min = -1024*pixel_size # (µm)
		y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
		x = np.arange(x_min, x_max+.001, pixel_size)
		y = np.arange(y_min, y_max+.001, pixel_size)
		X, Y = np.meshgrid(x, y)
	Return: 
	- I: intensity of the diffracted light on a screen orthogonal to the 
	incident direction at a distance Z from the object.
	"""

	# Parameters
	k = 2*np.pi/l
	
	delta_fx = 1./(x.size*pixel_size) # In order to be able to use FFT, 
	delta_fy = 1./(y.size*pixel_size) # see Eq. (4.4-12) of Ersoy (2006)
	fx_max = x.size/2*delta_fx # µm-1
	fy_max = y.size/2*delta_fy # µm-1
	fx = np.arange(-fx_max, fx_max-delta_fx+.0000001, delta_fx)
	fy = np.arange(-fy_max, fy_max-delta_fy+.0000001, delta_fy)
	# freq after fft = index until N ÷ total size (in distance)
	# x.size/2 and y.size/2 used to match with fftshift

	FX, FY = np.meshgrid(fx, fy)

	# x and y boundaries
	x_b = np.append(x - (x[1]-x[0])/2., x[-1] + (x[1]-x[0])/2.)
	y_b = np.append(y - (y[1]-y[0])/2., y[-1] + (y[1]-y[0])/2.)

	# Field at Z = 0
	U0 = np.copy(M)

	# Angular spectrum at Z = 0 (Fourier tranform), Eq. (4.3-2) of Ersoy (2006) 
	A0 = np.fft.fftshift(np.fft.fft2(U0))
	# fft2: Compute the 2-dimensional discrete Fourier Transform
	# fftshift: Shift the zero-frequency component to the center of the spectrum
	if False:
		plot_image(fy, fx, np.abs(A0))

	# Transfer function, Eq. (4.3-14) of Ersoy (2006)
	non_homogeneous_waves = 4*np.pi**2*(FX**2+FY**2) > k**2 # Eq. (4.3-13) of Ersoy (2006)
	H = np.exp( 1j*Z*np.abs(np.sqrt(k**2-4*np.pi**2*(FX**2+FY**2))) )
	H[non_homogeneous_waves] = 0

	# Angular spectrum at Z = Z, Eq. (4.3-11) of Ersoy (2006)
	AZ = A0*H

	# Lowpass filter to remove noise
	if lowpass_filter:
		filtre = np.zeros((y.size, x.size))
		A = (np.abs(FY) <= filter_fmax)
		B = (np.abs(FX) <= filter_fmax)
		filtre[A&B] = 1
		AZ = AZ*filtre

	# Field at Z = Z, Eq. (4.3-12) of Ersoy (2006)
	UZ = np.fft.ifft2(np.fft.fftshift(AZ))

	# Intensity at Z = Z
	I = np.abs(UZ)**2

	# Phase at Z = Z
	P = np.arctan(np.imag(UZ)/np.real(UZ))

	return I


def degrade_image(I, x, y, pixel_size, x_2DSpixel_size, y_2DSpixel_size):              
	"""
	Degrade image to 2D-S resolution
	"""

	# Initialization
	I_2DS = np.ones((129,129))*-9999. # 129 to allow centering
	y_2DS = np.ones(129)*-9999.
	x_2DS = np.ones(129)*-9999.
	x_center = np.where((x<0.1) & (x>-0.1))[0][0]
	y_center = np.where((y<0.1) & (y>-0.1))[0][0]

	# Pixel of I should be (1.1 x 1.1 µm)
	# Decompose each pixel of I in 100 (0.1 x 0.1 µm) => easier averaging after
	pixel_len_decomp = 10
	I_0_1 = np.repeat(I, pixel_len_decomp, axis=0)
	I_0_1 = np.repeat(I_0_1, pixel_len_decomp, axis=1)
	x_center_0_1 = np.int(x_center*pixel_len_decomp + pixel_len_decomp/2)
	y_center_0_1 = np.int(y_center*pixel_len_decomp + pixel_len_decomp/2)

	# Average at 2D-S resolution
	for i_2DS in np.arange(129):
		i = i_2DS - 64 # to center
		for j_2DS in np.arange(129):
			j = j_2DS - 64
			min_y = y_center_0_1-np.int(y_2DSpixel_size*pixel_len_decomp/2) +\
					i*np.int(y_2DSpixel_size*pixel_len_decomp)
			max_y = y_center_0_1-np.int(y_2DSpixel_size*pixel_len_decomp/2) +\
					(i+1)*np.int(y_2DSpixel_size*pixel_len_decomp)
			min_x = x_center_0_1-np.int(x_2DSpixel_size*pixel_len_decomp/2) +\
					j*np.int(x_2DSpixel_size*pixel_len_decomp)
			max_x = x_center_0_1-np.int(x_2DSpixel_size*pixel_len_decomp/2) +\
					(j+1)*np.int(x_2DSpixel_size*pixel_len_decomp)
			I_2DS[i_2DS, j_2DS] = np.mean(I_0_1[min_y:max_y, min_x:max_x])
													 # x=vertical, y=horizontal
		
	y_2DS = np.arange(-64*y_2DSpixel_size, 64*y_2DSpixel_size+.0001, 
					  y_2DSpixel_size)

	x_2DS = np.arange(-64*x_2DSpixel_size, 64*x_2DSpixel_size+.0001, 
					  x_2DSpixel_size)

	return I_2DS, x_2DS, y_2DS


def plot_diffraction_image(Z, I, x, y, I_2DS, x_2DS, y_2DS, y_2DSpixel_size):
	"""
	Return the diffracted image and the binary one which would result 
	from the 2D-S measurement.
	Note: x = vertical, y = horizotal
	"""
													
	# Colormap
	print('\tI min = %.f' % np.min(I.T*100) + ' %')
	print('\tI max = %.f' % np.max(I.T*100) + ' %')
	vmax = 140
	my_cmap = cm.Greys_r

	# x and y boundaries
	x_b = np.append( x - (x[1]-x[0])/2., x[-1] + (x[1]-x[0])/2.)
	y_b = np.append( y - (y[1]-y[0])/2., y[-1] + (y[1]-y[0])/2.)

	with sns.axes_style("whitegrid", {'axes.edgecolor':'0.', 
									  'grid.color':'1.',
									  'xtick.color':'0.', 
									  'xtick.direction':'out', 
									  'xtick.major.size':5, 
									  'xtick.minor.size':2,
									  'ytick.color':'0.', 
									  'ytick.direction':'out', 
									  'ytick.major.size':5, 
									  'ytick.minor.size':2}):
		
		# Create figure
		fig = plt.figure(figsize=(19/2.54, 6/2.54))

		# Create subplots
		gs0 = gridspec.GridSpec(1, 4, width_ratios=[6, 0.5, 4, 6], wspace=0.2)
		
		#-----------------------------
		# Plot diffraction image
		ax = plt.subplot(gs0[0])
		ax.set_aspect('equal')
		pc = plt.pcolormesh(y_b, x_b, I.T*100, cmap=my_cmap)
		plt.clim(0, vmax) # same as in the other 
		plt.xlabel('$y\ (\mu m)$')
		plt.ylabel('$x\ (\mu m)$')
		ax.xaxis.set_major_locator(MultipleLocator(500))
		ax.xaxis.set_minor_locator(MultipleLocator(100))
		ax.yaxis.set_major_locator(MultipleLocator(500))
		ax.yaxis.set_minor_locator(MultipleLocator(100))
		ax.xaxis.set_ticks_position('both')
		ax.yaxis.set_ticks_position('both')
		plt.xlim(-y_2DSpixel_size*128/2, y_2DSpixel_size*128/2)
		plt.ylim(-y_2DSpixel_size*128/2, y_2DSpixel_size*128/2)
														# 2DS = 128 photodiodes
		plt.title('Theoretical diffracted image', fontsize=12, y=1.1)

		#-----------------------------
		# Colorbar
		ax = plt.subplot(gs0[1])
		cb = plt.colorbar(pc, cax=ax, label='Light intensity (%)', extend='max')
		cb.set_ticks(np.array((0, 50, 100, 140)))
		ax.yaxis.set_minor_locator(MultipleLocator(10))
		# cb.ax.yaxis.set_minor_locator(FixedLocator(np.arange(0, vmax+.1, 10.)/vmax))
		
		#-----------------------------
		# Z distance
		plt.text(0.52, -0.07, 'Z = %.2f cm' % (Z/1e4), ha='center', va='bottom', 
				color='royalblue', fontsize=14, transform=fig.transFigure)

		#-----------------------------
		# 2D-S image
		ax = plt.subplot(gs0[3])
		ax.set_aspect('equal')

		# I with 50 % threshold
		I_2DS[I_2DS<=0.5] = 0.
		I_2DS[I_2DS>0.5] = 1.

		my_cmap = cm.Greys_r

		# x and y boundaries
		x_2DS_b = np.append(x_2DS     - (x_2DS[1]-x_2DS[0])/2., 
							x_2DS[-1] + (x_2DS[1]-x_2DS[0])/2.)
		y_2DS_b = np.append(y_2DS     - (y_2DS[1]-y_2DS[0])/2.,
							y_2DS[-1] + (y_2DS[1]-y_2DS[0])/2.)

		# If no pixel detection
		if np.min(I_2DS) > 0.5:
			ax.text(0.5, 0.5, u'No detection', ha='center', va='center', 
					fontsize=10, fontweight='bold', color='r', 
					transform=ax.transAxes)
			my_cmap = cm.Greys # to have a white image (not black)
			
		# Plot 2D-S image
		pc = plt.pcolormesh(y_2DS_b, x_2DS_b, I_2DS.T, cmap=my_cmap)

		plt.xlabel('$y\ (\mu m)$')
		plt.ylabel('$x\ (\mu m)$')
		ax.xaxis.set_major_locator(MultipleLocator(500))
		ax.xaxis.set_minor_locator(MultipleLocator(100))
		ax.yaxis.set_major_locator(MultipleLocator(500))
		ax.yaxis.set_minor_locator(MultipleLocator(100))
		ax.xaxis.set_ticks_position('both')
		ax.yaxis.set_ticks_position('both')
		plt.xlim(-y_2DSpixel_size*128/2, y_2DSpixel_size*128/2) 
		plt.ylim(-y_2DSpixel_size*128/2, y_2DSpixel_size*128/2)
														# 2DS = 128 photodiodes
		plt.title('Theoretical 2D-S record', fontsize=12, y=1.1)

	return


# **********************************************************************
# MAIN PROGRAM
# **********************************************************************

# Configuration
l, pixel_size, x_2DSpixel_size, y_2DSpixel_size = load_config()
l = l/1000. # nm to µm
x_min = -1024*pixel_size # (µm)
x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
y_min = -1024*pixel_size # (µm)
y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
x = np.arange(x_min, x_max+.001, pixel_size)
y = np.arange(y_min, y_max+.001, pixel_size)
X, Y = np.meshgrid(x, y)
filter_fmax = 0.06 # (µm-1) low freq filter on FFT to remove noise

# Load particle shape and distance
Z, M = load_shape_and_distance(x, y, X, Y)
Z = Z*10000 # cm to µm

# Compute diffracted image I
I = compute_diffraction(Z, l, pixel_size, x, y, X, Y, M, filter_fmax)

# Degrade image to 2D-S resolution
I_2DS, x_2DS, y_2DS = degrade_image(I, x, y, pixel_size, x_2DSpixel_size, 
									y_2DSpixel_size)

# Plot
plot_diffraction_image(Z, I, x, y, I_2DS, x_2DS, y_2DS, y_2DSpixel_size)
fig_filename_2DS = "diffraction_image.png"
plt.savefig(fig_filename_2DS, bbox_inches='tight', dpi=500)
print("%s saved" % fig_filename_2DS)