import numpy as np
from scipy.interpolate import interp1d

from phoenix_model import get_phoenix_model
from radmc3d_wavelength_setup import default_wavelengths
from natconst import um, Lsun, sig_boltz


# the wavelengths write to this file must contain the same wavelengths
# get stellar model that is converted to cgs units
# interpolate star model on to default radmc wavelengths then extrapolate
# write star file


def star_setup(TEFF, LBOL, logg, metalicity = 0, dts = 1):
    
    # wavelengths for wavelength_micron.inp file
    radmc_lam_um = default_wavelengths()
    
    # calculate star radius
    rstar = np.sqrt(LBOL * Lsun / (4 * np.pi * sig_boltz * TEFF**4))

    # get star phoenix model [lambda: um, flux density: erg/(sec*cm^2*Hz)]
    star_lam, star_fnu = get_phoenix_model(TEFF, LBOL, logg, metalicity, dts)
    
    # smoothing: interpolate and extrapolate onto radmc wavelenghts
    log_lam, log_fnu = np.log10([star_lam, star_fnu])
    log_fnu_interp   = interp1d(log_lam, log_fnu, fill_value="extrapolate")
    log_fnu_new      = log_fnu_interp(np.log10(radmc_lam_um))
    radmc_fnu_star   = 10**log_fnu_new


    with open('stars.inp', 'w+') as f:
        # Write formating section
        f.write('2\n')                # iformat
        f.write(f'1 {radmc_lam_um.size}\n')             # number of stars and wavelenths
        f.write(f'{rstar} 0 0 0 0\n') # Stellar information: centered

        # Write wavelengths [um] first then flux [cgs]
        for lam in radmc_lam_um:
            f.write(f'{lam}\n')
        for fnu in radmc_fnu_star:
            f.write(f'{fnu}\n')

# For testing
# star_setup(5800, 1, 3.5)

