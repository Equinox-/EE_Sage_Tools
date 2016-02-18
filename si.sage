def undoEscape(s):
    return s.replace("\v", "\\v").replace("\t", "\\t").replace("\a", "\\a").replace("\p", "\\p").replace("\f", "\\f").replace("\n", "\\n");

SI_prefixes = [ [1e-12, "\\pico"], [1e-9, "\\nano"], [1e-6, "\\micro"], [1e-3, "\\milli"], [1, ""], [1e3, "\\kilo"], [1e6, "\\mega"], [1e9, "\\giga"]]

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
    return r"\SI{%s}{%s}" % (va.str(),undoEscape(unit))

def SI_smart(val, unit, precision):
    va = SI_numeric(val, precision)
    prev = [1, ""]
    for prefix in SI_prefixes:
        if abs(va) < prefix[0]:
            return SI(val / prev[0], prev[1] + unit, precision)
        prev = prefix
    return SI(val, unit, precision)
