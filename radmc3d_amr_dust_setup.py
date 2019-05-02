import numpy as np
from natconst import au, Msun

''' 
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
'''


def dust_grain_weights(logamin_um, logamax_um, na, alpha, matdens):
    '''
    This function computes the weightings for a powerlaw distribution: n(a)~a**alpha, a grain size distribution.
    However, this distribution is to be converted to a grain mass distribution who bin-walls are centered at location
    sqrt(m_i * m_i+1)
    '''
    
    # convert from grain size distribution to grain mass distribution
    dustGrainSizeGrid     = 10.**(np.linspace(logamin_um, logamax_um, na))*1e-4 # [cm]
    dustMassGrainGrid     = (4.*np.pi/3)* matdens*grainSizeGrid**3              # [g] 
    massDistributionGrid  = dustMassGrainGrid**((alpha - 2.)/3.)
    
    # initialize bin-walls of mass distribution and populate with center of walls (geometric mean)
    dustMassGrainCellWalls       = np.zeros(na+1)
    dustMassGrainCellWalls[0]    = dustMassGrainGrid[0]
    dustMassGrainCellWalls[-1]   = dustMassGrainGrid[-1]
    dustMassGrainCellWalls[1:-1] = np.sqrt(dustMassGrainGrid[1:]*dustMassGrainGrid[:-1])
    
    deltaMassGrain = dustMassGrainGrid[1:]-dustMassGrainGrid[:-1]
    dummy          = (massDistributionGrid*dustMassGrainGrid*deltaMassGrain).sum()
    
    # Normalize distribtion
    massDistributionGrid = massDistributionGrid/dummy
    
    # Calculate the weights
    weights = massDistributionGrid*dustMassGrainGrid*deltaMassGrain
    return weights 

def vertical_dust_settling(Hg, matdens, gasdens, agrain, tau):
    ''''''
    return Hg / np.sqrt(matdens * agrain / (gasdens * Hg * tau))


def amr_dust_setup(rin_au, rc_au, rout_au, gamma, H100au, beta, Mdust_Msun, 
                   logamin_um, logamax_um, alpha, na, matdens,
                   nr_cells_inner, nr_cells_outer, r_scale, nphi_cells,
                   ntheta_cells_low, ntheta_cells_high, theta_min, 
                   theta_split):
    '''
    '''
    
    # constants used and coverted to CGS units if required
    r_in     = rin_au * au
    r_out    = rout_au * au
    r_change = r_scale * r_in
    rc       = rc_au * au 
    R0       = 100. * au
    Mdust    = Mdust_Msun * Msun
    H100     = H100au * au
    twoPI    = 2. * np.pi
    
    '''# ----- Setup radial cell walls ----- #'''
    
    # Note:
    # to avoid duplication, endpoint is set to FALSE where rchange 
    # will not be accounted for in inner grid but instead recovered in outer
    
    radiusCellWallsInner = np.linspace(r_in, r_change, nr_cells_inner + 1, False)       # inner FINE grid style
    radiusCellWallsOuter = np.geomspace(r_change, r_out, nr_cells_outer + 1)            # outer COARSE grid style
    radiusCellWalls      = np.concatenate([radiusCellWallsInner, radiusCellWallsOuter]) # complete grid of radial cell walls

    ''' # ----- Setup polar angle theta cell walls ----- #'''
    thetaCellWallsHigh = np.linspace(theta_min, theta_split, ntheta_cells_high + 1, False) # COARSE grid spacing
    thetaCellWallsLow  = np.linspace(theta_split, np.pi/2., ntheta_cells_low + 1)       # FINE grid spacing near the mid-plane

    # define cell walls above midplane to be reflected (inverted) for each region to be uniform

    thetaAboveMidplane = np.concatenate([thetaCellWallsHigh, thetaCellWallsLow])  # 0 to pi/2 
    thetaBelowMidplane = np.pi - thetaAboveMidplane[::-1][1:]                     # pi/2 to pi
    thetaCellWalls     = np.concatenate([thetaAboveMidplane, thetaBelowMidplane]) # completed grid of theta cell walls 

    '''# ----- Setup polar angle phi cell walls ----- #'''
    phiCellsWalls = np.linspace(0, twoPI, nphi_cells + 1)
    
    '''# ----- Find grid of cell walls center and number of cell walls ----- #'''
    sphericalCoordCellWalls  = [radiusCellWalls, thetaCellWalls, phiCellsWalls]
    sphericalCellWallCenters = list(map(lambda x: (0.5 * (x[:-1] + x[1:]), x.size - 1), sphericalCoordCellWalls))    
    
    radiusCellCenters , nrCells     = sphericalCellWallCenters[0]
    thetaCellCenters  , nthetaCells = sphericalCellWallCenters[1]
    phiCellCenters    , nphiCells   = sphericalCellWallCenters[2]
    
    '''# ----- Calculate dust grain mass distribution weights ----- #'''
    ''' INSERT DUST STUFF HERE '''
    dustGrainMassWeights = dust_grain_weights(logamin_um, logamax_um, alpha, na, matdens)
    
    
    '''# ----- Calculate dust volume density for each cell for a single polar angle ----- #'''
    # Note: for loops = BORING, Numpy = :D
    
    # Dust density is independent of phi which means we can evaluate a single polar angle then tile
    # Dust density is calculated on a grid using cylindrical coordinates ~ rho(r,z) thus from convert
    # spherical coords
    
    radiusSphericalGrid   , thetaGrid = np.meshgrid(radiusCellCenters, thetaCellCenters)
    radiusCylindricalGrid , zGrid     = map(lambda r, theta: (r * np.sin(theta), r * np.cos(theta)))
    
    '''### INSERT CHECK CONDITI0NS HERE ###'''
    
    #Scale the radial components to the disk characteristic size
    radiusScalingGrid = radiusCylindricalGrid / rc
    
    # Function used to value at disk inner and outer bounds
    viscosityProifile = np.exp(-radiusScalingGrid**(2 - gamma)) 
    
    # surface density scaled from the difference of viscosity at the inner region and outer
    # while the entire grid requires all values
    
    initialSurfaceDensity = Mdust * (gamma - 2) / (twoPI * rc**2) / (viscosityProifile[-1] - viscosityProifile[0])
    surfaceDensityGrid    = initialSurfaceDensity * radiusScalingGrid * viscosityProifile
    scaleHeightGrid       = H100 * (radiusCylindricalGrid / R0)**beta
    volumeDensityGrid     = surfaceDensityGrid / (scaleHeightGrid * twoPI**0.5) * np.exp(-0.5 * (zGrid/scaleHeightGrid)**2)
    '''### INSERT CHECK HERE ###'''
    
    '''# ----- Make the full disk density grid ----- #'''
    diskVolumeDensityGrid = np.tile(volumeDensityGrid.flatten(), nphiCells)
    
    '''# ----- write input files ----- #'''
    
    with open('amr_grid.inp','w') as f:
        # Write formating section
        f.write('1\n')                                    # iformat
        f.write('0\n')                                    # AMR grid style  (0=regular grid, no AMR)
        f.write('100\n')                                  # Coordinate system
        f.write('0\n')                                    # grid info
        f.write('1 1 1\n')                                # include coordinate: r, phi, theta   #
        f.write(f'{nrCells} {nthetaCells} {nphiCells}\n') # number of grid cells for each coordinate

        # Write grid values on one line
        # line #: a_1 a_2 a_3 ... a_n+1
        for r in radiusCellWalls:
            f.write(f'{r}' + ' ')
            
        f.write('\n')
        
        for phi in phiCellsWalls:
            f.write(f'{phi}' + ' ')
            
        f.write('\n')
        
        for theta in thetaCellWalls:
            f.write(f'{theta}' + ' ')
    
    
    with open('dust_density.inp', 'w') as f:
        f.write('1\n')                             # iformat
        f.write(f'{diskVolumeDensityGrid.size}\n') # nrcells
        f.write('1\n')                             # dust speicies
        
        for denisty in diskVolumeDensityGrid:
            f.write(f'{density}\n') 
    
    
    
