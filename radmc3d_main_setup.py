

"""
Tyler Baines
March 31, 2019
This scipt writes the main input settings file for to run radmc3d
"""

def radmc_main_setup(nphot = 1e7):
    # may need corrections
    with open('radmc3d.inp', 'w+') as f:
        f.write(f'nphot = {nphot}\n')
        f.write('istar_sphere = 1\n')
        f.write('scattering_mode_max = 3\n')