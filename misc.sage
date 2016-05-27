def numpy_to_sage(v):
    if v in ZZ:
        return ZZ(v)
    elif v in QQ:
        return QQ(v)
    elif isinstance(v, sage.symbolic.expression.Expression):
        return v
    elif v in RR:
        return RR(v)
    else:
        return CC(v)

def constIneq(s):
    return isinstance(s, bool) or isinstance(s, numpy.bool_) or (s.lhs().is_constant() and s.rhs().is_constant())

def evalFn(fn, t, var=None):
    if isinstance(fn, sage.symbolic.function.Function):
        return fn(t)
    elif hasattr(fn, '__call__'):
        try:
            if var == None:
                var = fn.default_variable();
            return fn({var:t})
        except:
            return fn(t)
    else:
        return fn

import re

def cleanTeX(v):
    return cleanStR(latex(v))

def cleanStR(v):
    return re.sub(r'\.[0]+([^0-9])', r'\1', re.sub(r'(\.[0-9]+?)[0]+([^0-9])', r'\1\2', v + " ")).rstrip(" ")

def __replRoundDecimal(match):
    begin,cut,suffix = match.groups()
    try:
        mod=begin[-1]
        decider=cut[0]
        if int(decider) >= 5:
            # round required
            i = -1
            copy=list(begin)
            while i>=-len(copy):
                ch=copy[i]
                if ch != '.':
                    va=int(ch)+1
                    if va <= 9:
                        copy[i]=str(va)
                        break
                    else:
                        copy[i]='0'
                i=i-1
            begin=''.join(copy)
            if i<-len(copy):
                begin = "1" + begin
    except e:
        print(e)
        pass
    return begin+suffix

# Count >= 0
def roundDecimal(ss, count):
    return re.sub(r'([^0-9])1(|\.[0]+)[ ]+\\times (10\^)', r'\1\3', re.sub(r'([0-9]*\.[0-9]{'+str(count)+'})([0-9]+)([^0-9])', __replRoundDecimal, " " + ss + " ")).rstrip(" ").lstrip(" ")


def numApproxSR(v, perc):
    if hasattr(v, "is_constant") and v.is_constant() and (not(hasattr(v, "is_relational")) or not(v.is_relational())):
        if not(v.is_integer()) or abs(v) > 9999:
            if not(v in RR):
                v = CC(numApproxSR(real(v), perc), numApproxSR(imag(v), perc))
            else:
                v = v.n(perc)
    elif isinstance(v, sage.symbolic.expression.Expression) and v.operator() is not None:
        operands=v.operands()[:]
        for j in xrange(0, len(operands)):
            operands[j] = numApproxSR(operands[j], perc)
        v=v.operator()(*operands)
    return v

def stringToTable(string, delim):
    tb=string.replace("\r", "").split("\n")
    csv=[]
    for i in xrange(0, len(tb)):
        if len(tb[i].strip()) > 0:
            csv.append(tb[i].split(delim))
    return csv

def readFileFully(path):
    content = ""
    with open(path, 'r') as content_file:
        content = content_file.read()
    return content

def readCSV(path):
    table = stringToTable(readFileFully(path), ",")
    for k in table:
        for j in xrange(0, len(k)):
            k[j] = RR(k[j])
    return table