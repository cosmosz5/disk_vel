import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pdb

def matrix_index(map_sz):
    sz = map_sz
    if sz % 2 == 0:  ### Define the grid (x, and y are inverted !!!!!)
        x, y = np.mgrid[-np.floor(sz / 2 - 1):np.floor(sz / 2):sz * 1j, np.floor(sz / 2 - 1):-np.floor(sz / 2):sz * 1j]
    else:
        x, y = np.mgrid[-np.floor(sz / 2):np.floor(sz / 2):sz * 1j, np.floor(sz / 2):-np.floor(sz / 2):sz * 1j]

    return x, y