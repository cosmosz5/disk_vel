import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pdb
import importlib
import intensity_profiles as iprof
import kinematic_disk as kd
import astropy.constants as const
import scipy.signal as signal
import oi_tools
importlib.reload(kd)
importlib.reload(iprof)
importlib.reload(oi_tools)

##### Star parameters ##########
m_star = 4  ##In m_sun units
r = np.linspace(1., 50., 50)### In AU

##### Image parameters #####
imsize = 129 ###Pixels
imscale = 0.001 ### Milliarcsec/pixel

##### Disk parameters #####
FWHM = 0.1 ## In milliarcseconds
PA = 50.0 ## In degrees
incl = 60.0 ## In degrees
dist =  150.0 ### Parsecs

##### Wavelength range #####
eff_wave = np.linspace(2.16e-6, 2.17e-6, 10)
eff_band = eff_wave[1] - eff_wave[0]
lambda0 = 2.1661e-6
vline = np.array((eff_wave - lambda0) / lambda0 * const.c / 1000.0)

###################################################
##### AUTOMARIC FROM HERE #########################

##### Create the Intensity profile of the disk #####
Iprof = iprof.elliptical_gauss(imsize, imscale, FWHM, PA, incl)
Iprof = Iprof / np.sum(Iprof)
pyfits.writeto('test.fits', Iprof, overwrite=True)

##### Compute the 2D projected velocity map #########
vel, velx, vely  = kd.keplerian_vel2D(m_star, imsize, imscale, dist)

pyfits.writeto('test_vel.fits', vel, overwrite=True)
pyfits.writeto('test_vel_x.fits', velx, overwrite=True)
pyfits.writeto('test_vel_y.fits', vely, overwrite=True)

vel_proj = kd.vel_proj(velx, vely, PA, incl)
pyfits.writeto('test_vel_proj.fits', vel_proj, overwrite=True)

##### Compute the iso-velocity region ###############

RR = np.zeros([len(eff_wave), imsize, imsize])
Doppler_vel = np.array(const.c / 1000.0 * (eff_wave - lambda0) / lambda0)
inc = np.array(const.c / 1000.0 * (eff_band) / lambda0)
for i in range(len(eff_wave)):
    iso_map = np.zeros(vel_proj.shape)
    ind = np.where( (vel_proj >= Doppler_vel[i] - inc/2 ) & (vel_proj <= Doppler_vel[i] + inc/2 ) )
    iso_map[ind] = 1.0

    RR[i, :,:] = iso_map * Iprof
pyfits.writeto('iso_vel.fits', RR, overwrite=True)
pdb.set_trace()

