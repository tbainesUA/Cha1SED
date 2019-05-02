import numpy as np

"""
Tyler Baines
March 29,2019

This script is used to generate the wavelength input file in microns that wil be read radmc.

"""
def default_wavelengths():
    #default wavelengths for setup in microns
    nlam    = 200 
    lam_min = 0.1
    lam_max = 5000
    return np.geomspace(lam_min, lam_max, nlam)

def wavelength_setup():
    # writes the wavelength points for continuum radi-trans calculation
    with open('wavelength_micron.inp', 'w') as f:
        lam = default_wavelengths()
        f.write(f'{len(lam)}' + '\n')
        for l in lam:
            f.write(f'{l}' + '\n')
    
def wavelength_camera_setup(lam):
    # required input wavelengths i.e. flux observation at some wavelengths
    # writes the wavelength points to be used to produce spectra, SED, or images
    with open(' camera wavelength micron.inp', 'w') as f:
        f.write(f'{len(lam)}' + '\n')
        for l in lam:
            f.write(f'{l}' + '\n')
    

"""if __name__ == "__main__":
    wavelength_setup()
    wavelength_camera_setup(np.linspace(0,10,11))"""
