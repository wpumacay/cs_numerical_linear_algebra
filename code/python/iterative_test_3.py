import numpy as np

MIN = -3.0
STEP = 1.0
MAX = 3.0
N = int( ( MAX - MIN ) / STEP )

x1_set = [ MIN + STEP * q for q in range( N ) ]
x2_set = [ MIN + STEP * q for q in range( N ) ]
x3_set = [ MIN + STEP * q for q in range( N ) ]
x4_set = [ MIN + STEP * q for q in range( N ) ]
x5_set = [ MIN + STEP * q for q in range( N ) ]
x6_set = [ MIN + STEP * q for q in range( N ) ]


def getSpectralRadius_Gauss( x1, x2, x3, x4, x5, x6 ) :
    coeffs = [1, 
              x1 * x3 * x5 - x1 * x4 - x2 * x5 - x3 * x6,
              x2 * x4 * x6]

    leigs = np.roots( coeffs )
    maxnorm = np.abs( leigs[0] )
    for q in range( 1, len( leigs ) ) :
        _norm = np.absolute( leigs[q] )
        if ( _norm > maxnorm ) :
            maxnorm = _norm

    return maxnorm

def getSpectralRadius_Jacobi( x1, x2, x3, x4, x5, x6 ) :
    coeffs = [1,
              0, 
              -( x1 * x4 + x2 * x5 + x3 * x6 ),
              x1 * x3 * x5 + x2 * x4 * x6]

    leigs = np.roots( coeffs )
    maxnorm = np.abs( leigs[0] )
    for q in range( 1, len( leigs ) ) :
        _norm = np.absolute( leigs[q] )
        if ( _norm > maxnorm ) :
            maxnorm = _norm

    return maxnorm

cases = []

for x1 in x1_set :
    for x2 in x2_set :
        for x3 in x3_set :
            for x4 in x4_set :
                for x5 in x5_set :
                    for x6 in x6_set :
                        _sr_j = getSpectralRadius_Jacobi( x1, x2, x3, x4, x5, x6 )
                        _sr_g = getSpectralRadius_Gauss( x1, x2, x3, x4, x5, x6 )
                        
                        if ( _sr_j < 1 and _sr_g > 1 ) :
                            print '_sr_j: ' , _sr_j
                            print '_sr_g: ' , _sr_g
                            cases.append( [x1,x2,x3,x4,x5,x6] )

for case in cases :
    print 'case: ', case
