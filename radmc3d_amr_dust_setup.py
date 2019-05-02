import numpy as np
from natconst import au, Msun

""" 
Tyler Baines
March 28, 2019

This script requires two input that are lists of model parameters that will be use
to generate Radmc-3d input files. A spherical coordinate system is used for all modelling
purposes. CGS units are used through out all computations as well as fine (f) and coarse 
(c) grids. 

Input list must have the following format:

    _Disk parameters_                        _Grid parameters_ ('#' = number)
    rin_au     : Inner Radius            |   nr_cells_inner    : # of inner cells     
    rc_au      : Characteristic Radius   |   nr_cells_outer    : # of outter cells
    rout_au    : Outer Radius            |   r_scale           : scale factor applied to r_in 
    gamma      : Spatial Index           |   nphi_cells        : # of polar angles
    H100au     : Scale Height at 100 AU  |   ntheta_cells_low  : # altitude angle near midplane
    beta       : Flaring Index           |   ntheta_cells_high : # altitude angle near z-axis
    Mdust_Msun : Dust mass in Sol Mass   |   theta_min         : minimum angle 
                                         |   theta_split       : angle grid style change
"""

def amr_dust_setup(rin_au, rc_au, rout_au, gamma, H100au, beta, Mdust_Msun,
                   nr_cells_inner, nr_cells_outer, r_scale, nphi_cells,
                   ntheta_cells_low, ntheta_cells_high, theta_min, 
                   theta_split):

        
    # constants used and coverted to CGS units if required
    r_in     = rin_au * au
    r_out    = rout_au * au
    r_change = r_scale * r_in
    rc       = rc_au * au 
    R0       = 100 * au
    Mdust    = Mdust_Msun * Msun
    H100     = H100au * au
    twoPI    = 2 * np.pi
    Sigma0   = (
        (Mdust * (gamma - 2)) / 
        (twoPI * rc**2) / 
        (np.exp(-(rout_au/rc_au)**(2 - gamma)) - np.exp(-(rin_au/rc_au)**(2 - gamma)))
    )

    
    # ----- Setup radial cell walls ----- #

    # Note:
    # to avoid duplication, endpoint is set to FALSE where rchange 
    # will not be accounted for in inner grid but instead recovered in outer
    
    r_inner = np.linspace(r_in, r_change, nr_cells_inner + 1, False) # inner FINE grid style
    r_outer = np.geomspace(r_change, r_out, nr_cells_outer + 1)      # outer COARSE grid style
    r_grid  = np.concatenate([r_inner, r_outer])                     # complete grid of radial cell walls

    # ----- Setup polar angle theta cell walls ----- #
    theta_high     = np.linspace(theta_min, theta_split, ntheta_cells_high + 1, False) # COARSE grid spacing
    theta_low      = np.linspace(theta_split, np.pi/2, ntheta_cells_low + 1)           # FINE grid spacing near the mid-plane

    # define cell walls above midplane to be reflected (inverted) for each region to be uniform

    theta_above_midplane = np.concatenate([theta_high, theta_low])            # 0 to pi/2 
    theta_below_midplane = np.pi - theta_above_midplane[::-1][1:]             # pi/2 to pi
    theta_grid = np.concatenate([theta_above_midplane, theta_below_midplane]) # completed grid of theta cell walls 

    # ----- Setup polar angle phi cell walls ----- #
    phi_grid = np.linspace(0, twoPI, nphi_cells + 1)


    # ----- Find coords of cell walls center ----- #
    r_center_grid     = 0.5 * (r_grid[:-1] + r_grid[1:])
    phi_center_grid   = 0.5 * (phi_grid[:-1] + phi_grid[1:])
    theta_center_grid = 0.5 * (theta_grid[:-1] + theta_grid[1:])

    nr_cen, nphi_cen, ntheta_cen = list(map(lambda x: x.size, [r_center_grid, phi_center_grid, theta_center_grid]))
    print(f'nr = {nr_cen}, {len(r_grid)}')
    print(f'nt = {ntheta_cen},{len(theta_grid)}')
    print(f'np = {nphi_cen},{len(phi_grid)}')

    # ----- Calculate dust density for each cell ----- #
    
    # Method to avoid looping, although phi array not important because dust model equations do not 
    # depend on the phi.
    phi, theta, r_sph = np.meshgrid(phi_center_grid, theta_center_grid, r_center_grid, indexing = 'ij')
    



    # spherical to cylindrical coord transform
    r_cyln = r_sph * np.sin(theta)
    z      = r_sph * np.cos(theta) 

    # compute scale heights
    H = H100 * (r_cyln / R0)**beta
    
    # compute dust surface density
    Sigma_dust = (
        Sigma0 * 
        (r_cyln / rc)**(-gamma) * 
        np.exp(-(r_cyln / rc)**(2 - gamma))
    )

    # compute volume dust density 
    rho_dust_density = (
        Sigma_dust / (np.sqrt(twoPI) * H) *
        np.exp(-0.5 * (z / H)**2)
    )
    
    # ----- write files ----- #

    with open('amr_grid.inp','w') as f:
        # Write formating section
        f.write('1\n')                              # iformat
        f.write('0\n')                              # AMR grid style  (0=regular grid, no AMR)
        f.write('100\n')                            # Coordinate system
        f.write('0\n')                              # grid info
        f.write('1 1 1\n')                          # include coordinate: r, phi, theta   #
        f.write(f'{nr_cen} {nphi_cen} {ntheta_cen}' + '\n') 

        # Write grid values on one line
        # line #: a_1 a_2 a_3 ... a_n+1
        for r in r_grid:
            f.write(f'{r}' + ' ')
        f.write('\n')
        for phi in phi_grid:
            f.write(f'{phi}' + ' ')
        f.write('\n')
        for theta in theta_grid:
            f.write(f'{theta}' + ' ')
    
    with open('dust_density.inp', 'w') as f:
        f.write('1\n')                     # iformat
        f.write('{}\n'.format(rho_dust_density.size)) # nrcells
        f.write('1\n')                     # dust speicies
        
        for rho in rho_dust_density.flatten():
            f.write(f'{rho}' + '\n') 
    
    return 

if __name__ == "__main__":
    rin_au = 1
    rc_au = 5
    rout_au = 20
    gamma = 1
    H100au = 10
    beta = 2
    Mdust_Msun = 1e-4

    ppdisk = [rin_au, rc_au, rout_au, gamma, H100au, beta, Mdust_Msun]

    nr_cells_inner    = 5
    nr_cells_outer    = 10
    r_scale           = 2
    ntheta_cells_low  = 5
    ntheta_cells_high = 5
    theta_split       = np.pi/8
    theta_min         = 0
    nphi_cells        = 3
    spatial_grid = [nr_cells_inner, nr_cells_outer, r_scale, nphi_cells,
                   ntheta_cells_low, ntheta_cells_high, theta_min, 
                   theta_split]

    amr_dust_setup(*ppdisk, *spatial_grid)