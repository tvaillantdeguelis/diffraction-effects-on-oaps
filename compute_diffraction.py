# Arguments needed for the following function are:
# - Z: distance to object plane (µm)
# - l: wavelength (µm)
# - pixel_size: (µm)
# - M: transmittance mask (particle shape) with size (y.size, x.size)
# - filter_fmax: (µm-1) low freq filter on FFT to remove noise (you can turn to False in the function if not needed)
# - x, y, X, Y: obtained with:
# x_min = -1024*pixel_size # (µm)
# x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
# y_min = -1024*pixel_size # (µm)
# y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
# x = np.arange(x_min, x_max+.001, pixel_size)
# y = np.arange(y_min, y_max+.001, pixel_size)
# X, Y = np.meshgrid(x, y)

def compute_diffraction(Z, l, pixel_size, x, y, X, Y, M, filter_fmax):

	# Parameters
	k = 2*np.pi/l
	
	delta_fx = 1./(x.size*pixel_size) # In order to be able to use FFT, Eq. (4.4-12) of Ersoy (2006)
	delta_fy = 1./(y.size*pixel_size)
	fx_max = x.size/2*delta_fx # µm-1
	fy_max = y.size/2*delta_fy # µm-1
	fx = np.arange(-fx_max, fx_max-delta_fx+.0000001, delta_fx)
	fy = np.arange(-fy_max, fy_max-delta_fy+.0000001, delta_fy)
	# freq après fft = indices jusqu'à N ÷ taille totale (en distance)
	# x.size/2 et y.size/2 utilisés pour coller avec fftshift

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
	if True:
		filtre = np.zeros((y.size, x.size))
		A = (np.abs(FY) <= filter_fmax)
		B = (np.abs(FX) <= filter_fmax)
		filtre[A&B] = 1
		if False:
			plot_image(y_b, x_b, filtre)
		AZ = AZ*filtre

	# Field at Z = Z, Eq. (4.3-12) of Ersoy (2006)
	UZ = np.fft.ifft2(np.fft.fftshift(AZ))

	# Intensity at Z = Z
	I = np.abs(UZ)**2

	# Phase at Z = Z
	P = np.arctan(np.imag(UZ)/np.real(UZ))

	return I