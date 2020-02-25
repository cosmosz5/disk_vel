import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import pdb

oi_data = pyfits.open('COMB_SCI_BrG_HD45677.fits')
n_continuum = 4 ## Number of channels on each side of the wavelength range to be used for the continuum


###### Read the data from the file ########
eff_wave = oi_data['OI_WAVELENGTH'].data['EFF_WAVE']
eff_band = oi_data['OI_WAVELENGTH'].data['EFF_BAND']
oi_vis = oi_data['OI_VIS'].data['VISAMP']
oi_vis_err = oi_data['OI_VIS'].data['VISAMPERR']
oi_phi = oi_data['OI_VIS'].data['VISPHI']
oi_phi_err = oi_data['OI_VIS'].data['VISPHIERR']
ucoord = oi_data['OI_VIS'].data['UCOORD']
vcoord = oi_data['OI_VIS'].data['VCOORD']
oi_flux = oi_data['OI_FLUX'].data['FLUX']
oi_flux_err = oi_data['OI_FLUX'].data['FLUXERR']


fig1, (ax1) = plt.subplots(1,1)
#### Fit the line flux ########
mean_flux = np.mean(oi_flux, axis = 0)
mean_flux_err = np.sqrt(np.sum(oi_flux_err**2, axis=0))
## Normalize by the continuum #####
cont1 = mean_flux[0:4]
cont2 = mean_flux[4:]
eff_wave_cont1 = eff_wave[0:4]
eff_wave_cont2 = eff_wave[4:]

ax1.plot(cont1, 'or')
ax1.errorbar(range(len(mean_flux)), mean_flux, yerr=mean_flux_err)

plt.show()
pdb.set_trace()