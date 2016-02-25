def bandPassFreq(H, fMin=1, fMax=1e12):
    omega = var('omega')
    return find_local_maximum(lambda omega: phasorMag(evalFn(H, I*omega)), fMin, fMax, maxfun=500)[1]

def bandStopFreq(H, fMin=1, fMax=1e12):
    omega = var('omega')
    return find_local_minimum(lambda omega: phasorMag(evalFn(H, I*omega)), fMin, fMax, maxfun=500)[1]
