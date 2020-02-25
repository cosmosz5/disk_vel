import numpy as np
import astropy.io.fits as pyfits
import pdb


def elliptical_gauss(map_sz, scale_mas, FWHM, PA, incl):
    sz = map_sz
    if sz % 2 == 0: ### Define the grid (x, and y are inverted !!!!!)
        x, y = np.mgrid[-np.floor(sz/2-1):np.floor(sz/2):sz*1j, np.floor(sz/2-1):-np.floor(sz/2):sz* 1j]
    else:
        x, y = np.mgrid[-np.floor(sz/2):np.floor(sz/2):sz*1j, np.floor(sz/2):-np.floor(sz/2):sz* 1j]

    x = x.flatten()
    y = y.flatten()
    xy = np.vstack((x, y))
    c, s = np.cos(np.deg2rad(PA)), np.sin(np.deg2rad(PA))
    R_x = np.matrix([[c, -s], [s, c]])
    xy_prime = np.dot(xy.T, R_x)


    sigma = (FWHM / scale_mas) / 2.3548
    #g = np.exp(-1.*(np.array(xy_prime[:, 0])**2. + xy_prime[:, 1]**2.))
    g = np.exp(-1.*(np.array(xy_prime[:,0])**2. + np.array(xy_prime[:,1])**2. * (1. / np.cos(np.deg2rad(incl))**2.)) / (2 * sigma**2.))

    g = np.reshape(g, (sz, sz))
    return g
