import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import astropy.constants as const
import pdb


def keplerian_vel(m_star, r):
    ### Mstar should be in meters!
    ### r: radius should be in meters
    m_star = m_star * const.M_sun
    r = r * const.au
    V_sigma = np.sqrt((const.iau2015.G * m_star) / r)
    return V_sigma


def vel_proj(V_sigmax, V_sigmay, PA, incl): ### Keplerian Velocity
    #vproj = (V_sigma * np.sin(np.deg2rad(PA))) * np.sin(np.deg2rad(incl))
    vproj = - (-V_sigmay * np.sin(np.deg2rad(PA)) - V_sigmax * np.cos(np.deg2rad(PA))) * np.sin(np.deg2rad(incl))
    return vproj

def keplerian_vel2D(m_star, map_sz, scale_mas, dist):
    ### Mstar should be in meters!
    ### r: radius should be in meters
    sz = map_sz
    if sz % 2 == 0:  ### Define the grid (x, and y are inverted !!!!!)
        x, y = np.mgrid[-np.floor(sz / 2 - 1):np.floor(sz / 2):sz * 1j, np.floor(sz / 2 - 1):-np.floor(sz / 2):sz * 1j]
    else:
        x, y = np.mgrid[-np.floor(sz / 2):np.floor(sz / 2):sz * 1j, np.floor(sz / 2):-np.floor(sz / 2):sz * 1j]

    ind_x = np.where(x == 0)
    ind_y = np.where(y == 0)

    x[ind_x] = 1e-8
    y[ind_y] = 1e-8

    m_star = m_star * const.M_sun
    rr = np.sqrt(np.array(x)**2 + np.array(y)**2)
    au_mas = (dist) / 1000.0 ## To convert the distance to kilo-parsecs
    rr = rr * scale_mas * au_mas * const.au ### To convert coordinates to au
    rr_x = x * scale_mas * au_mas * const.au ### To convert coordinates to au
    rr_y = y * scale_mas * au_mas * const.au ### To convert coordinates to au
    V_sigma = np.sqrt((const.iau2015.G * m_star) / rr) / 1000.0 ### To convert m -> km
    Vx_sigma = V_sigma * np.sin(np.arctan2(x,y))
    Vy_sigma = V_sigma * np.cos(np.arctan2(x,y))

    return np.array(V_sigma), np.array(Vx_sigma), np.array(Vy_sigma)