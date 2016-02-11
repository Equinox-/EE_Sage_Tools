def numpy_to_sage(v):
    if v in ZZ:
        return ZZ(v)
    elif v in QQ:
        return QQ(v)
    elif v in CC:
        return CC(v)
    else:
        return RR(v)

def evalFn(fn, t):
    if isinstance(fn, sage.symbolic.function.Function):
        return fn(t)
    elif hasattr(fn, '__call__'):
        return fn({fn.default_variable():t})
    else:
        return fn