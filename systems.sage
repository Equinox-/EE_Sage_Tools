bandPassMinFreq=1
bandPassMaxFreq=1e15

def tfWrap(H):
    return lambda x: phasorMag(evalFn(H, I*x));

def bandPassFreq(H, fMin=bandPassMinFreq, fMax=bandPassMaxFreq):
    return find_local_maximum(tfWrap(H), fMin, fMax, maxfun=500)[1]

def bandStopFreq(H, fMin=bandPassMinFreq, fMax=bandPassMaxFreq):
    return find_local_minimum(tfWrap(H), fMin, fMax, maxfun=500)[1]

def bandPassInfo(H, fMin=bandPassMinFreq, fMax=bandPassMaxFreq, center=None):
    omega = var('omega')
    if center == None:
        center = bandPassFreq(H, fMin, fMax)
    fmag=phasorMag(evalFn(H, I*center))/sqrt(2)
    lkfn = (lambda omega: phasorMag(evalFn(H, I*omega)) - fmag)
    left=find_root(lkfn, fMin, center - 1)
    right=find_root(lkfn, center + 1, fMax)
    return [left, center, right]

def bandStopInfo(H, fMin=bandPassMinFreq, fMax=bandPassMaxFreq, center=None):
    omega = var('omega')
    if center == None:
        center = bandStopFreq(H, fMin, fMax)
    fmag=phasorMag(evalFn(H, I*infinity))/sqrt(2)
    lkfn = (lambda omega: phasorMag(evalFn(H, I*omega)) - fmag)
    left=find_root(lkfn, fMin, center - 1)
    right=find_root(lkfn, center + 1, fMax)
    return [left, center, right]