def fourier_series_symbolic(a, N, var):
    ee = 0
    for k in xrange(0, len(a)):
        ee += a[k] * exp(I*k*2*pi/N*var)
    return ee