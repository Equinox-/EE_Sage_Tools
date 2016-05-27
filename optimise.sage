def findRootMultidim(func, vrs, guess, iters=10):
    jdef=[]
    for pt in xrange(0, len(func)):
        rw=[]
        for v in vrs:
            rw.append(func[pt].diff(v))
        jdef.append(rw)
    J=matrix(jdef)
    J_inv=J.I
    func = vector(func)
    R = vector(guess).n(32)
    for rrr in xrange(0, iters):
        dct=dict(zip(vrs, R))
        dp = J_inv.subs(dct).n(32) * -(func.subs(dct).n(32))
        R = R + dp
    return dict(zip(vrs,R))