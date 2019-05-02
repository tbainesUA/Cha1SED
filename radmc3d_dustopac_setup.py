import numpy as np
import os
import shutil

"""
Tyler Baines
March 29, 2019

This script generates the dust opacity input file that will be read by radmc.
Dust Opacity Library used contains DSharp opacities with mie scattering matrices
at various grain sizes. file name to be 

"""

def dustopac_setup(logagrain):
    
    
    # grain sizes used to make dust opac files in logaritmic space
    # Fine gridding up to 4 coarse grid to 5
    # this reason is because of time to make dust opac model bhmie code
    available_logagrain_S1 = np.linspace(-2,5,1000)                             # orginal models to made
    available_logagrain_S1 = available_logagrain_S1[available_logagrain_S1 < 4] # make correction here
    available_logagrain_S2 = np.linspace(4,5,25)                                # coarse grid of large grains
    available_logagrain_S2 = np.append(available_logagrain_S2, 4.033)
    
    # combine both fine and coarse grids and round values to 4 to match dust opacity filename format
    available_logagrain = np.round(np.concatenate([available_logagrain_S1, available_logagrain_S2]) ,4)
    
    
    parentDIR = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    

    # loop through grain sizes that are in dust opacity library and copy them over
    for agrain in logagrain:
        
        # condition to find closest respectable value of grain size model file
        agrain_close = available_logagrain[(np.abs(available_logagrain - agrain)).argmin()]
        
        
        dust_file    = f'dustkapscatmat_DSHARPmix_{agrain_close}_-1_5_500_181_0.05_20_5.0_True.inp'
        
        source       = os.path.join(parentDIR, 'DustKappaScatmat Library', dust_file)
        destination  = dust_file[:-34] + '.inp'
        
        
        shutil.copy2(source, destination)
        # end of loop
    
    # get all dust species name that were copied to be written to file
    dust_species = [i for i in os.listdir() if i.startswith('dustkapscatmat_DSHARPmix_')]
    nspecies = len(dust_species)
    
    with open('dustopac.inp', 'w+') as f:
        # Write formating section
        f.write('2\n')                             # iformat
        f.write(f'{nspecies}\n')                   # number species to be used
        f.write('-----------------------------\n') # Stellar information: centered

        # Write species sections 
        for i, dopac in enumerate(dust_species):
            f.write(f'inputstyle[{i}]\n')
            f.write('iquantum[0]\n')
            f.write(f'{dopac}\n')
            f.write('-----------------------------\n')

    return


#Testing purposes
#dustopac_setup(np.array([-2.01, 0.003]))
#dustopac_setup([-2.01])
#$dustopac_setup(-2.01, 0.003)
