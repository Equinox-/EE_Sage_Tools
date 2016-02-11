def fromPhasor(mag, phase):
    rm = sqrt(mag*mag / (1 + tan(phase)^2))
    return rm + I*tan(phase)*rm

def phasorLatex(v, precision):
    return SI_numeric(phasorMag(v), precision).str(skip_zeroes=1) + "\phase{" + SI_numeric(phasorPhase(v)*180/pi, precision).str(skip_zeroes=1) + "^{\circ}}"

def phasorPhase(v):
    return atan2(imag(v),real(v))

def phasorMag(v):
    return sqrt(real(v)^2 + imag(v)^2)