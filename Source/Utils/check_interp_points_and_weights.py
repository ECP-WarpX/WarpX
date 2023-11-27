# Copyright 2020 Edoardo Zoni
#
# This file is part of WarpX
#
# License: BSD-3-Clause-LBNL

#-------------------------------------------------------------------------------
# Compute interpolation points and weights for coarsening and refinement in IO
# and MR applications in 1D (extensions to 2D and 3D are trivial). Weights are
# computed in order to guarantee total charge conservation for both cell-centered
# data (equal weights) and nodal data (weights depend on distance between points
# on fine and coarse grids).
#
# Notation:
# - index i  refers to points on coarse grid
# - index ii refers to points on fine   grid
# - sc denotes the staggering of the coarse data
# - sf denotes the staggering of the fine   data
# - cr denotes the coarsening ratio (must be cr=1,2,4 only)
#
# For MR applications only the cases sc=sf=0 and sc=sf=1 are considered. Terms
# multiplied by (1-sf)*(1-sc) are ON for cell-centered data and OFF for nodal data,
# while terms multiplied by sf*sc are ON for nodal data and OFF for cell-centered
# data. C++ implementation in Source/ablastr/coarsen/average.(H/.cpp) and
# Source/ablastr/coarsen/sample.(H/.cpp)
#-------------------------------------------------------------------------------

import sys

import numpy as np

# # Refinement for MR: interpolation points and weights
# def refinement_points_and_weights( ii, sc, sf, cr ):
#     if   ( cr==1 ):
#         numpts = 1
#         idxmin = ii
#     elif ( cr>=2 ):
#         if   ( ii%cr==0 ):
#            numpts = (1-sf)*(1-sc)+sf*sc
#         elif ( ii%cr!=0 ):
#            numpts = (1-sf)*(1-sc)+2*sf*sc
#         idxmin = (ii//cr)*(1-sf)*(1-sc)+(ii//cr)*sf*sc
#     weights = np.zeros( numpts )
#     for ir in range( numpts ):
#         i = idxmin+ir
#         if   ( ii==iimin or ii==iimax ):
#            weights[ir] = (1-sf)*(1-sc)+((abs(cr-abs(ii-i*cr)))/(cr)+(cr/2-0.5))*sf*sc
#         else:
#            weights[ir] = (1-sf)*(1-sc)+((abs(cr-abs(ii-i*cr)))/(cr))*sf*sc
#     return [ numpts, idxmin, weights ]

def coarsening_fine_grid_limits( sc, sf, cr ):
    if   ( sf == 0 ): # cell-centered
        iimin = 1
        iimax = 4*cr
    elif ( sf == 1 ): # nodal
        iimin = 0
        iimax = 4*cr
    return [ iimin, iimax ]

def coarsening_coarse_grid_limits( sc, sf, cr, iimin, iimax ):
    # imin = int( iimin/cr ) # original
    # imax = int( iimax/cr )-(1-sc)*sf+(1-sf)*sc # original
    imin = 10000
    imax = -10000
    for i in range(-100,101):
        numpts, idxmin, weights = coarsening_points_and_weights(i, sc, sf, cr)
        if iimin <= idxmin + numpts - 1:
            imin = min(imin,i)
        if iimax >= idxmin:
            imax = max(imax,i)
    return [ imin, imax ]

def coarsening_points_and_weights( i, sc, sf, cr ):
    twoImin = -cr*sc + sf - 1
    if (twoImin % 2 == 0):
        idxmin = i*cr + twoImin//2
        numpts = cr+1
        weights = np.zeros( numpts )
        weights[0] = 1.0/(2*cr)
        weights[numpts-1] = 1.0/(2*cr)
        for ir in range( 1, numpts-1 ):
            weights[ir] = 1.0/cr
    else:
        idxmin = i*cr + (twoImin+1)//2
        numpts = cr
        weights = np.zeros( numpts )
        for ir in range( 0, numpts ):
            weights[ir] = 1.0/cr

    return [ numpts, idxmin, weights ]

# def refinement_fine_grid_limits( sc, sf, cr ):
#     if   ( sf == 0 ): # cell-centered
#         iimin = 0 # original
#         iimax = 7 # original
#     elif ( sf == 1 ): # nodal
#         iimin = 0 # original
#         iimax = 8 # original
#     return [ iimin, iimax ]

# def refinement_coarse_grid_limits( sc, sf, cr, iimin, iimax ):
#     # imin = int( iimin/cr ) # original
#     # imax = int( iimax/cr )-(1-sc)*sf+(1-sf)*sc # original
#     imin = 10000
#     imax = -10000
#     for i in range(-100,101):
#         numpts, idxmin, weights = refinement_points_and_weights(i, sc, sf, cr)
#         if iimin <= idxmin + numpts - 1:
#             imin = min(imin,i)
#         if iimax >= idxmin:
#             imax = max(imax,i)
#     return [ imin, imax ]

# # Refinement for MR: interpolation points and weights
# def refinement_points_and_weights( ii, sc, sf, cr ):
#     if   ( cr==1 ):
#         numpts = 1
#         idxmin = ii
#     elif ( cr>=2 ):
#         if   ( ii%cr==0 ):
#             numpts = (1-sf)*(1-sc)+sf*sc
#         elif ( ii%cr!=0 ):
#             numpts = (1-sf)*(1-sc)+2*sf*sc
#         idxmin = (ii//cr)*(1-sf)*(1-sc)+(ii//cr)*sf*sc
#     weights = np.zeros( numpts )
#     for ir in range( numpts ):
#         i = idxmin+ir
#         if   ( ii==iimin or ii==iimax ):
#             weights[ir] = (1-sf)*(1-sc)+((abs(cr-abs(ii-i*cr)))/(cr)+(cr/2-0.5))*sf*sc
#         else:
#             weights[ir] = (1-sf)*(1-sc)+((abs(cr-abs(ii-i*cr)))/(cr))*sf*sc
#     return [ numpts, idxmin, weights ]

## TODO Coarsening for IO: interpolation points and weights
#def coarsening_points_and_weights_for_IO( i, sf, sc, cr ):
#    if   ( cr==1 ):
#        numpts = 1+abs(sf-sc)
#        idxmin = i-sc*(1-sf)
#    elif ( cr>=2 ):
#        numpts = 2-sf
#        idxmin = i*cr+cr//2*(1-sc)-(1-sf)
#    weights = np.zeros( numpts )
#    for ir in range( numpts ):
#        weights[ir] = (1/numpts)*(1-sf)*(1-sc)+(1/numpts)*sf*sc
#    return [ numpts, idxmin, weights ]

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Input coarsening ratio
cr = int( input( "\n Select coarsening ratio (cr=1,2,4): cr=" ) )
# if ( cr!=1 and cr!=2 and cr!=4 ):
#     print()
#     sys.exit( 'coarsening ratio cr={} is not valid'.format( cr ) )

# Loop over possible staggering of coarse and fine grid (cell-centered or nodal)
for sc in [0,1]:
    for sf in [0,1]:

        print( '\n **************************************************' )
        print( ' * Staggering of coarse grid: sc={}'.format( sc ), end='' )
        if ( sc == 0 ):
            print( ' cell-centered  *' )
        elif ( sc == 1 ):
            print( ' nodal          *' )
        print( ' * Staggering of fine   grid: sf={}'.format( sf ), end='' )
        if ( sf == 0 ):
            print( ' cell-centered  *' )
        elif ( sf == 1 ):
            print( ' nodal          *' )
        print( ' **************************************************' )

        # if ( sf!=sc ):
        #     print( '\n WARNING: sc={} not equal to sf={}, not implemented for MR, continue ...'.format( sc, sf ) )
        #     continue

        print( '\n Coarsening for MR: check interpolation points and weights' )
        print( ' ---------------------------------------------------------' )

        iimin,iimax = coarsening_fine_grid_limits( sc, sf, cr )
        imin, imax  = coarsening_coarse_grid_limits( sc, sf, cr, iimin, iimax )

        print( '\n Min and max index on coarse grid:  imin={}  imax={}'.format( imin, imax ) )
        print(   ' Min and max index on fine   grid: iimin={} iimax={}'.format( iimin, iimax ) )

        # Number of grid points
        # nc = imax-imin+1
        # nf = iimax-iimin+1
        # print( '\n Number of points on coarse grid: nc={}'.format( nc ) )
        # print(   ' Number of points on fine   grid: nf={}'.format( nf ) )

        # Coarsening for MR: interpolation points and weights
        for i in range( imin, imax+1 ): # index on coarse grid
            numpts,idxmin,weights = coarsening_points_and_weights( i, sc, sf, cr )
            print( '\n Find value at i={} by interpolating over the following points and weights:'.format( i ) )
            wtotal = 0.0
            for ir in range( numpts ): # interpolation points and weights
                ii = idxmin+ir
                wtotal += weights[ir]
                print( ' (ii={},w={:.3f})'.format( ii, weights[ir] ), end='' )
                if not ( ir == numpts-1 ):
                    print( ' ', end='' )
            print()
            if (abs(wtotal - 1.0) > 1E-9):
                print ('\n ERROR: total weight wtotal={} should be 1 for coarse index i={}'.format(wtotal,i))

        # Coarsening for MR: check conservation properties
        for ii in range( iimin, iimax+1 ): # index on fine grid
            ws = 0.0
            for i in range( imin, imax+1 ): # index on coarse grid
                numpts,idxmin,weights = coarsening_points_and_weights( i, sc, sf, cr )
                for ir in range( numpts ): # interpolation points and weights
                    jj = idxmin+ir
                    if ( jj==ii ): # interpolation point matches point on fine grid
                       ws += weights[ir]
            if (abs(ws - 1.0/cr) > 1E-9):
                print( '\n ERROR: sum of weights ws={} should be 1/cr={} for ii={}'.format( ws, 1.0/cr, ii ) )

        # print( '\n Refinement for MR: check interpolation points and weights' )
        # print( ' ---------------------------------------------------------' )

        # iimin,iimax = refinement_fine_grid_limits( sc, sf, cr )
        # imin ,imax  = refinement_coarse_grid_limits( sc, sf, cr, iimin, iimax )
        # # Number of grid points
        # nc = imax-imin+1
        # nf = iimax-iimin+1

        # print( '\n Min and max index on coarse grid:  imin={}  imax={}'.format( imin, imax ) )
        # print(   ' Min and max index on fine   grid: iimin={} iimax={}'.format( iimin, iimax ) )
        # print( '\n Number of points on coarse grid: nc={}'.format( nc ) )
        # print(   ' Number of points on fine   grid: nf={}'.format( nf ) )

        # # Refinement for MR: interpolation points and weights
        # for ii in range ( nf ): # index on fine grid
        #     numpts,idxmin,weights = refinement_points_and_weights( ii, sc, sf, cr )
        #     print( '\n Find value at ii={} by interpolating over the following points and weights:'.format( ii ) )
        #     for ir in range( numpts ): # interpolation points and weights
        #         i = idxmin+ir
        #         sign = '+'
        #         if (i == 0):
        #             sign = ' '
        #         if (i < 0):
        #             sign = '-'
        #         print( ' (i='+sign+'{},w={:.2f})'.format( abs(i), weights[ir] ), end='' )
        #         if not ( ir == numpts-1 ):
        #             print( ' ', end='' )
        #     print()

        # # Refinement for MR: check conservation properties
        # for i in range( nc ): # index on coarse grid
        #     ws = 0.0
        #     for ii in range( nf ): # index on fine grid
        #         numpts,idxmin,weights = refinement_points_and_weights( ii, sc, sf, cr )
        #         for ir in range( numpts ): # interpolation points and weights
        #             jj = idxmin+ir
        #             if ( jj==i ): # interpolation point matches point on coarse grid
        #                 ws += weights[ir]
        #     if ( abs(ws - cr) > 1E-9 ):
        #         print( '\n ERROR: sum of weights ws={:.2f} should be cr={} for i={}'.format( ws, cr, i ) )
