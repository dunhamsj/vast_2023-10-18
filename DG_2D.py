#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

def Lagrange( x, xq, i ):

    L = 1.0

    for j in range( xq.shape[0] ):

        if j != i:

            L *= ( x - xq[j] ) / ( xq[i] - xq[j] )

    return L

def rhoh( x, y, rho_q, xq ):

    N = xq.shape[0]

    rho = 0.0
    for i in range( N ):
        for j in range( N ):
            rho += rho_q[i,j] * Lagrange( x, xq, i ) * Lagrange( y, xq, j )

    return rho

def GetQuadrature( N ):

    if N == 1:

        wq = np.array( [ +1.0 ], np.float64 )
        xq = np.array( [ +0.0 ], np.float64 )

    elif N == 2:

        wq = np.array( [ +1.0, +1.0 ], np.float64 ) / 2.0
        xq = np.array( [ -1.0, +1.0 ], np.float64 ) / np.sqrt( 12.0 )

    elif N == 3:

        wq = np.array( [ +5.0, +8.0, +5.0 ], np.float64 ) / 18.0
        xq = np.array( [ -1.0, +0.0, +1.0 ], np.float64 ) \
               * np.sqrt( 3.0 / 20.0 )

    else:

      print( 'Invalid choice of N: {:}'.format( N ) )
      exit( 'Exiting...' )

    return wq, xq

def CreateElement():

    fig, ax = plt.subplots( 1, 1 )
    ax.axis( 'off' )

    xM = [ -0.5, +0.5 ]
    yM = [ -0.5, +0.5 ]
    ax.plot( [ xM[0], xM[1] ], [ yM[0], yM[0] ], 'k-' )
    ax.plot( [ xM[1], xM[1] ], [ yM[0], yM[1] ], 'k-' )
    ax.plot( [ xM[1], xM[0] ], [ yM[1], yM[1] ], 'k-' )
    ax.plot( [ xM[0], xM[0] ], [ yM[1], yM[0] ], 'k-' )

    ax.set_aspect( 'equal' )

    return fig, ax

def AddQuadraturePoints( ax, N ):

    wq, xq = GetQuadrature( N )

    for i in range( N ):
        for j in range( N ):
            ax.plot( xq[i], xq[j], 'r.' )

    return

def AddFacePoints( ax, N ):

    wq, xq = GetQuadrature( N )

    for i in range( N ):
        ax.plot( -0.5, xq[i], 'rs', mfc = 'None' )
        ax.plot( +0.5, xq[i], 'rs', mfc = 'None' )
        ax.plot( xq[i], -0.5, 'rs', mfc = 'None' )
        ax.plot( xq[i], +0.5, 'rs', mfc = 'None' )

def rhoExact( x, y ):
    return 0.25 + x * y + x**2 * y + y**2 * x + y**4

def PlotExact( fig, ax ):

    x = np.linspace( -0.5, +0.5, 100 )
    y = np.linspace( -0.5, +0.5, 100 )

    rho = np.empty( (x.shape[0],y.shape[0]), np.float64 )
    for i in range( rho.shape[0] ):
        for j in range( rho.shape[0] ):
            rho[i,j] = rhoExact( x[i], y[j] )

    vmin = rho.min()
    vmax = rho.max()

    eps = 0.001
    im = ax.imshow( rho.T, \
                    extent = [ -0.5+eps, +0.5-eps, -0.5+eps, +0.5-eps ], \
                    origin = 'lower', \
                    vmin = vmin, vmax = vmax )

    cbar = fig.colorbar( im )

    return vmin, vmax

def PlotDensity( N, fig, ax, vmin, vmax ):

    wq, xq = GetQuadrature( N )

    rho_q = np.empty( (N,N), np.float64 )

    for i in range( N ):
        for j in range( N ):
            rho_q[i,j] = rhoExact( xq[i], xq[j] )

    x = np.linspace( -0.5, +0.5, 100 )
    y = np.linspace( -0.5, +0.5, 100 )

    rho = np.empty( (x.shape[0],y.shape[0]), np.float64 )
    for i in range( rho.shape[0] ):
        for j in range( rho.shape[0] ):
            rho[i,j] = rhoh( x[i], y[j], rho_q, xq )

    eps = 0.001
    im = ax.imshow( rho.T, \
                    extent = [ -0.5+eps, +0.5-eps, -0.5+eps, +0.5-eps ], \
                    origin = 'lower', \
                    vmin = vmin, vmax = vmax )

    return

if __name__ == '__main__':

    fs = 17

    fig, ax = CreateElement()
    ax.set_title( 'Exact', fontsize = fs )
    vmin, vmax = PlotExact( fig, ax )
    plt.savefig( 'fig.DG_2D_Exact.png', dpi = 300 )
#    plt.show()
    plt.close()

    fig, ax = CreateElement()
    ax.set_title( r'$N=2$', fontsize = fs )
    N = 2
    PlotDensity( N, fig, ax, vmin, vmax )
    AddQuadraturePoints( ax, N )
    AddFacePoints      ( ax, N )
    plt.savefig( 'fig.DG_2D_N{:d}.png'.format( N ), dpi = 300 )
#    plt.show()
    plt.close()

    fig, ax = CreateElement()
    ax.set_title( r'$N=3$', fontsize = fs )
    N = 3
    PlotDensity( N, fig, ax, vmin, vmax )
    AddQuadraturePoints( ax, N )
    AddFacePoints      ( ax, N )
    plt.savefig( 'fig.DG_2D_N{:d}.png'.format( N ), dpi = 300 )
#    plt.show()
    plt.close()
