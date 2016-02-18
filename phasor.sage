def fromPhasor(mag, phase):
    rm = sqrt(mag*mag / (1 + tan(phase)^2))
    return rm + I*tan(phase)*rm

def phasorLatex(v, precision):
    return SI_numeric(phasorMag(v), precision).str(skip_zeroes=1) + "\phase{" + SI_numeric(phasorPhase(v)*180/pi, precision).str(skip_zeroes=1) + "^{\circ}}"

def phasorPhase(v):
    return atan2(imag(v),real(v))

def phasorMag(v):
    return abs(v)

def phasorToTime(phasor, omega, var):
    return phasorMag(phasor)*cos(omega * var + phasorPhase(phasor))

def phasorToTimeLatex(phasor, omega, precision):
    phas = phasorPhase(phasor)
    return SI_numeric(phasorMag(phasor), precision).str(skip_zeroes=1) + "\\cos\\left(" + SI_numeric(omega, precision).str(skip_zeroes=1) + "t " + ("+" if phas > 0 else "-") + " " + SI_numeric(phas, precision).str(skip_zeroes=1) + "\\right)"
