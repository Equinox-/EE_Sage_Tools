def fromPhasor(mag, phase):
    return mag*exp(CC(0,1)*phase)

def fromPhasorDeg(mag, phase):
    return fromPhasor(mag, phase*pi/180)

def phasorLatex(v, precision):
    return SI_str(SI_numeric(phasorMag(v), precision)) + "\phase{" + SI_str(SI_numeric(phasorPhase(v)*180/pi, precision)) + "^{\circ}}"

def phasorPhase(v):
    return atan2(imag(v),real(v))

def phasorPhase2(v):
    return ln(v / phasorMag(v))

def phasorMag(v):
    return abs(v)

def phasorMag2(v):
    return sqrt(real(v)^2 + imag(v)^2)

def phasorToTime(phasor, omega, var):
    return phasorMag(phasor)*cos(omega * var + phasorPhase(phasor))

def phasorToTimeLatex(phasor, omega, precision):
    phas = phasorPhase(phasor)
    return SI_str(SI_numeric(phasorMag(phasor), precision)) + "\\cos\\left(" + SI_str(SI_numeric(omega, precision)) + "t " + ("+" if phas > 0 else "-") + " " + SI_str(abs(SI_numeric(phas, precision))) + "\\right)"


def phasorToTimeLatex2(phasor, precision):
    phas = phasorPhase(phasor)
    return SI_str(SI_numeric(phasorMag(phasor), precision)) + "\\cos\\left(\\omega t " + ("+" if phas > 0 else "-") + " " + SI_str(abs(SI_numeric(phas, precision))) + "\\right)"