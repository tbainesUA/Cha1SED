import numpy as np
import os
from natconst import sig_boltz, pc, Lsun , c_light

"""
Tyler Baines
March 31,2019

Script that loads in a stellar atmosphere phoenix model and converts model wavelength and flux 
normalize to some distance (Default is 1 pc). Function returns model in cgs units
"""

def get_phoenix_model(TEFF, LBOL, logg, metalicity = 0, dts = 1):
    
    um = 1e-4 # conversion factor from angstroms to microns

    # Define the range of model temperatures & logg that exist in stellar library.
    TS1         = np.arange(1200, 2100, 100)                    # [K]
    TS2         = np.arange(2050, 2450, 50 )                    # [K]
    TS3         = np.arange(2500, 7000, 100)                    # [K]
    Temperature = np.concatenate([TS1, TS2, TS3]).astype(float) # [K]
    
    Logg = np.arange(2.5, 5.5, 0.5)

    # parent directory
    parentDIR = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    
    
    # condition aims to find the closest respectable model given an effective temperature/logg value
    # and return the value
    available_values    = [Temperature, Logg]
    given_values        = [TEFF, logg]
    get_nearest_index   = lambda available, given : np.argmin(np.abs(available - given))
    T_index, logg_index = list(map(get_nearest_index, available_values, given_values))
    
    
    Tstar    = Temperature[T_index]
    logg_mod = Logg[logg_index]

    # format file name string that will be loaded in
    star_file = f'lte0{Tstar/100}-{logg_mod}-{metalicity}.0a+0.0.BT-Settl.spec.7.txt'

    star_mod_path = os.path.join(parentDIR, 'PhoenixModel Library', star_file)

    # load model data [lambda angstrom, flambda per ang] 
    model = np.genfromtxt(star_mod_path) 

    # convert model wavelength and flux denisty values to micron dependencies
    lam_ang = model[:,0]   # angstroms
    lam_um  = lam_ang * um # microns


    flam_per_ang = model[:,1]        # ergs/(sec*cm^2*A) 
    flam_per_um  = flam_per_ang / um # ergs/(sec*cm^2*um)

    
    # converting flambda to fnu that is normalize to some source distance
    # dts -> distance to source 
    dts          = dts * pc                                          # [cm]
    Lstar        = LBOL * Lsun                                       # [erg/sec]
    T4           = Tstar**4                                          # [K4]
    surface_area = Lstar / (sig_boltz * T4)                          # [cm2]
    flam_at_dts  = flam_per_um * surface_area / (4 * np.pi * dts**2) # erg/(sec*cm^2*um) rescaled
    fnu_at_dts   = flam_at_dts * lam_um**2 / c_light / um            # erg/(sec*cm^2*Hz)

    return lam_um, fnu_at_dts


