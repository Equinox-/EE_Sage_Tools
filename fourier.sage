def fourier_series_symbolic(a, N, var):
    ee = 0
    for k in xrange(0, len(a)):
        ee += a[k] * exp(I*k*2*pi/N*var)
    return ee

def idt_fourier_transform(ft, wvar, nvar, recurse=True, algs=['sympy', 'maxima']):
    for alg in algs:
        try:
            return 1/(2*pi) * integral(ft * exp(I*wvar*nvar), wvar, -pi, pi, algorithm=alg)
        except:
            # Ignore
            print('alg ' + alg + ' failed')
    if not(recurse):
        raise Exception("Broken")
    # Expand and retry
    H = ft.expand()
    if H.operator() == sage.symbolic.operators.add_vararg:
        ops = H.operands()
        out = 0
        for op in ops:
            print("Recurse " + str(op))
            out += idt_fourier_transform(op, wvar, nvar, False)
        return out
    else:
        raise Exception("Yeah nope")


def fourier_transform_of_fourier_series(a, N, svar):
    ss = 0
    for k in xrange(0, len(a)):
        ss += 2*pi*numpy_to_sage(a[k])*dirac_delta(svar - k*2*pi/N)
    return ss

# Computes periodic sum over one period
def periodicSum(fn, svar, N, span=1):
    out = 0
    for k in xrange(-span, span+1):
        out += evalFn(fn, svar + N * k)
    return out

def cconv_sym(a, b, svar, N, span=1):
    var('tau')
    AN = periodicSum(a, svar, N, span)
    BN = periodicSum(b, svar, N, span)
    return integrate(evalFn(AN, tau) * evalFn(BN, svar - tau), tau, -N/2, N/2, algorithm='sympy')