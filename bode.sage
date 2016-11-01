# Encapsulate matplotlib figure inside figure
class MatPlotFigure:
    def __init__(self, mpl):
        self.mpl = mpl

    def save(self, filename=None, **kwds):
        dpi = 100 if not('dpi' in kwds) else round(kwds['dpi'])
        transparent = ('transparent' in kwds) and bool(kwds['transparent'])
        fig_tight = not('fig_tight' in kwds) or bool(kwds['fig_tight'])

        if filename is None:
            from sage.misc.superseded import deprecation
            deprecation(17234,'the filename argument is now mandatory')
            from sage.misc.temporary_file import graphics_filename
            filename = graphics_filename()
        ext = os.path.splitext(filename)[1].lower()

        from matplotlib import rcParams
        rc_backup = (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'], rcParams['text.usetex'])
        mympl = self.mpl
        mympl.tight_layout()

        opts = dict(dpi=dpi, transparent=transparent)
        if fig_tight is True:
            opts['bbox_inches'] = 'tight'
        mympl.savefig(filename, **opts)
        (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'], rcParams['text.usetex']) = rc_backup

# Encapsulate tachyon string inside figure
class TachyonFigure:
    def __init__(self, mpl):
        self.mpl = mpl

    def save(self, filename=None, **kwds):
        tachyon_rt(self.mpl, outfile=filename)

# Workaround for noframe figure
class Figure3D:
    def __init__(self, mpl):
        self.mpl = mpl

    def save(self, filename=None, **kwds):
        self.mpl.save_image(filename, **kwds)

# Creates a matplotlib instance of the given transfer function.
def bodePlot_data(freq, mag, phase, **kwargs):
    radians = ('radians' in kwargs) and bool(kwargs['radians'])
    decibel = ('decibel' in kwargs) and bool(kwargs['decibel'])
    square = ('square' in kwargs) and bool(kwargs['square'])
    pdb = ('pdb' in kwargs) and bool(kwargs['pdb'])
    angfreq = ('angfreq' in kwargs) and bool(kwargs['angfreq'])
    decibel |= pdb
    square |= pdb

    multi = isinstance(mag[0], list)
    channels = len(mag) if multi else 1

    bodeAmp = []
    bodePhase = []
    lastPhase = []

    for l in xrange(0, channels):
        bodeAmp.append([])
        bodePhase.append([])
        lastPhase.append([])
        magf = mag[l] if multi else mag
        phasef = phase[l] if multi else phase
        for j in xrange(0, len(magf)):
            m = magf[j]
            if decibel:
                if square:
                    m = m * m
                else:
                    m = abs(m)
                m = 10 * log(m) / log(10)
            p = phasef[j]
            while j > 0 and abs(p - lastPhase[l]) > pi:
                if p < lastPhase[l]:
                    p += 2 * pi
                else:
                    p -= 2 * pi
            lastPhase[l] = p
            if not(radians):
                p *= 180 / pi
            f = freq[j]
            if not(angfreq):
                f /= 2 * pi
            bodeAmp[l].append((f.n(), m.n()))
            bodePhase[l].append((f.n(), p.n()))
    import numpy
    import matplotlib
    import matplotlib.pyplot as plt
    fig = plt.figure()
    # Create the first Y-axis on the left by default and
    # plot the first curve associated to the left Y-axis
    ax1 = fig.add_subplot(111)
    ax1.set_xscale('log')
    if angfreq:
        ax1.set_xlabel('$\omega$')
    else:
        ax1.set_xlabel('$f$ (Hz)')
    ax1.set_ylabel('$|H(j\omega)|$')
    if decibel:
        ax1.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: '%g dB'%x))

    for l in xrange(0, channels):
        bodeAmp_x = zip(*(bodeAmp[l]))[0]
        bodeAmp_y = zip(*(bodeAmp[l]))[1]
        ax1.plot(bodeAmp_x, bodeAmp_y, '-', color='red')

    # Paint the tick labels on the left Y-axis in red 
    # to match the color of a curve
    for tl in ax1.get_yticklabels():
        tl.set_color('red')

    # Create the second Y-axis on the right and
    # plot the second curve associated to the right Y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('$\\angle H(j\omega)$')

    limMin = 1e10
    limMax = -1e10
    for l in xrange(0, channels):
        bodePhase_x = zip(*(bodePhase[l]))[0]
        bodePhase_y = zip(*(bodePhase[l]))[1]
        ax2.plot(bodePhase_x, bodePhase_y, '--', color='blue')
        if not(radians):
            limMin = min(limMin, floor(min(bodePhase_y) - .001))
            limMax = max(limMax, ceil(max(bodePhase_y) + .001))
    if not(radians):
        ax2.set_ylim(limMin, limMax)
        ax2.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: '$%.0f^\circ$'%x))
        count = max(len(ax1.get_ygridlines()) - 1, 1)
        step = ceil((limMax - limMin) / count + .001)
        ax2.yaxis.set_ticks(numpy.arange(limMin, limMax + step, step))

    # Paint the tick labels on the left Y-axis in blue 
    # to match the color of a curve
    for tl in ax2.get_yticklabels():
        tl.set_color('blue')

    # Create the vertical and horizontal grid lines
    ax1.xaxis.grid(color='grey', linestyle='--', linewidth=0.5)
    ax1.yaxis.grid(color='grey', linestyle='--', linewidth=0.5)

    # Save the figure (to see it in Sage)
    return plt

def bodePlot_func(h, fMin, fMax, **kwargs):
    angfreq = ('angfreq' in kwargs) and bool(kwargs['angfreq'])
    step = 1.1 if not('step' in kwargs) else float(kwargs['step'])
    if step <= 1:
        raise Exception("Bad step size")
    multi = isinstance(h, list)
    channels = len(h) if multi else 1

    freq = []
    mag = []
    phase = []
    x = fMin
    while x < fMax:
        f = x
        if not(angfreq):
            f *= 2 * pi.n()
        freq.append(f)
        x *= step

    for j in xrange(0, channels):
        mag.append([])
        phase.append([])
        for i in xrange(0, len(freq)):
            v = evalFn(h[j] if multi else h, I*freq[i])
            m = phasorMag(v)
            p = phasorPhase(v)
            mag[j].append(m)
            phase[j].append(p)
    return bodePlot_data(freq, mag, phase, **kwargs)
