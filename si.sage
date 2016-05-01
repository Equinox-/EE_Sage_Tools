ESC_UNDO = [["\a", "\\a"], ["\b", "\\b"], ["\c", "\\c"], ["\d", "\\d"], ["\e", "\\e"], ["\f", "\\f"], ["\g", "\\g"], ["\h", "\\h"], ["\i", "\\i"], ["\j", "\\j"], ["\k", "\\k"], ["\l", "\\l"], ["\m", "\\m"], ["\n", "\\n"], ["\o", "\\o"], ["\p", "\\p"], ["\q", "\\q"], ["\r", "\\r"], ["\s", "\\s"], ["\t", "\\t"], ["\w", "\\w"], ["\v", "\\v"], ["\y", "\\y"], ["\z", "\\z"]]

def undoEscape(s):
    for u in ESC_UNDO:
        s = s.replace(u[0], u[1])
    return s

SI_prefixes = [ [1e-18, "\\atto"], [1e-15, "\\femto"], [1e-12, "\\pico"], [1e-9, "\\nano"], [1e-6, "\\micro"], [1e-3, "\\milli"], [1, ""], [1e3, "\\kilo"], [1e6, "\\mega"], [1e9, "\\giga"], [1e12, "\\tera"], [1e15, "\\peta"], [1e18, "\\exa"]]

def SI_numeric(val, precision):
    if not(hasattr(val, 'n')):
        val = RR(val)
    return val.n(digits=precision)

def SI_str(va):
    return cleanStR(va.str(skip_zeroes=1).rstrip("."))

def SI(val, unit, precision):
    if hasattr(val, '__iter__'):
        st="\\left[\\begin{tabular}{" + "c"*len(val) + "}"
        flag=False
        for k in val:
            if flag:
                st+="&"
            st+=SI(k, "", precision)
            flag=True
        return st+"\\end{tabular}\\right]" + (r"\SI{}{%s}" % (undoEscape(unit)))
    va = SI_numeric(val, precision)
    if (imag(va) != 0):
        if (real(va) == 0):
            fix = "-j" if (imag(va) < 0) else "j"
            return fix + SI(abs(imag(va)), unit, precision)
        else:
            sig = "-" if (imag(va) < 0) else "+"
            return SI(real(va), "", precision) + " " + sig + " j" + SI(abs(imag(va)), unit, precision)
    if numpy.isfinite(va):
        return r"\SI{%s}{%s}" % (SI_str(va),undoEscape(unit))
    else:
        return r"%s \SI{}{%s}" % (SI_str(va),undoEscape(unit))

def SI_smart(val, unit, precision):
    va = SI_numeric(val, precision)
    if numpy.isfinite(va):
        prev = [1, ""]
        for prefix in SI_prefixes:
            if abs(va) < prefix[0]:
                return SI(val / prev[0], prev[1] + unit, precision)
            prev = prefix
    return SI(val, unit, precision)