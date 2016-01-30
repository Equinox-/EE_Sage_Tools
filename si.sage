def undoEscape(s):
    return s.replace("\v", "\\v").replace("\t", "\\t").replace("\a", "\\a").replace("\p", "\\p");

def SI_numeric(val, precision):
    if not(hasattr(val, 'n')):
        val = RR(val)
    return val.n(digits=precision)

def SI(val, unit, precision):
    va = SI_numeric(val, precision)
    if (imag(va) != 0):
        if (real(va) == 0):
            fix = "-j" if (imag(va) < 0) else "j"
            return fix + SI(abs(imag(va)), unit, precision)
        else:
            sig = "-" if (imag(va) < 0) else "+"
            return SI(real(va), "", precision) + " " + sig + " j" + SI(abs(imag(va)), unit, precision)
    return r"\SI{%s}{%s}" % (va.str(skip_zeroes=1),undoEscape(unit))