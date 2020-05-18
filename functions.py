import numpy as np

##
## Weigthed boxcar smoothing algorithm
##
def bcsmooth(spectrum):
  ## smooth to 3 pixels
  N = len(spectrum)
  z = [0. for i in range(N)]
  w = [1., 3., 1.]
  for i in range(1,N-2):
    c=0.
    for j in range(0,3):
      c += w[j]*spectrum[i-1+j]
    c /= 5.
    z[i] = c
  return z


##
## This code is a modified version of
## Algorithm D in Morhac, M. 2009, NIMA, 600, 478
##
## A pseudo-continuum is estimated via a
## Statistics-sensitive Nonlinear Iterative
## Peak-clipping (SNIP) algorithm.
## In the reference above, the method was
## designed to isolate 'peaks' from a
## 'background'. Below, an 'absorption line'
## is separated from the 'pseudo-continuum'
## in a stellar spectrum.
##
## INPUTS:
##    spectrum -> flux values
##    m        -> half the 'feature' size to smooth over
##    w        -> simultaneous smoothing parameter
##    BOXCAR   -> boxcar smooth spectrum before (true/false)
##
## OUTPUTS:
##    bkgr     -> pseudo-continuum values

def bkgrnd(spectrum, m=42, w=1, BOXCAR=False):
    N = len(spectrum)
    if BOXCAR:
        y = bcsmooth(spectrum)
    else:
        y = [spectrum[i] for i in range(N)]
    bkgr = [spectrum[i] for i in range(N)]
##    for p in range(m,0,-1):
    for p in range(m,w,-1):
        for i in range(p,N-p):
            a = y[i]
            b = (y[i-p] + y[i+p])/2.
            c = 0.
            for j in range(i-w,i+w+1):
                c += y[j]
            c /= (2.*w+1.)
            bkgr[i] = c if a>b else b
        for i in range(p,N-p):
            y[i] = bkgr[i]
    return bkgr

##
## Cosmic Ray Rejection Algorithm
##
def crreject( flux, level=10., niter=1, BOX=5, VERBOSE=False ):
    N = flux.size
    for j in range(niter):
        count = 0
        for i in range(BOX, N-(2*BOX+1)):
            faveL = np.mean(flux[(i-BOX):i])
            faveH = np.mean(flux[(i+BOX):(i+2*BOX)])
            fstdL = np.std(flux[(i-BOX):i])
            fstdH = np.std(flux[(i+BOX):(i+2*BOX)])
            fave = 0.5*(faveL+faveH)
            fstd = 0.5*(fstdL+fstdH)
            if( flux[i] > fave + level*fstd ):
                flux[i] = ( (faveH-faveL)/(2*BOX-1)*(BOX/2.) + faveL )
                count+=1
        if( VERBOSE ): print('Replaced %d cosmic rays.'%count)
