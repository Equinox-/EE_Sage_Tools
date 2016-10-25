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

    bodeAmp = []
    bodePhase = []
    exPhase = 0
    lastPhase = 0
    for j in xrange(0, len(mag)):
        m = mag[j]
        if decibel:
            if square:
                m = m * m
            else:
                m = abs(m)
            m = 10 * log(m) / log(10)
        p = phase[j] + exPhase
        while j > 0 and abs(p - lastPhase) > pi:
            if p < lastPhase:
                p += 2 * pi
            else:
                p -= 2 * pi
        lastPhase = p
        if not(radians):
            p *= 180 / pi
        f = freq[j]
        if not(angfreq):
            f /= 2 * pi
        bodeAmp.append((f.n(), m.n()))
        bodePhase.append((f.n(), p.n()))

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
    bodeAmp_x = zip(*bodeAmp)[0]
    bodeAmp_y = zip(*bodeAmp)[1]
    if decibel:
        ax1.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: '%g dB'%x))
    ax1.plot(bodeAmp_x, bodeAmp_y, '-', color='red')

    # Paint the tick labels on the left Y-axis in red 
    # to match the color of a curve
    for tl in ax1.get_yticklabels():
        tl.set_color('red')

    # Create the second Y-axis on the right and
    # plot the second curve associated to the right Y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('$\\angle H(j\omega)$')
    bodePhase_x = zip(*bodePhase)[0]
    bodePhase_y = zip(*bodePhase)[1]
    if not(radians):
        limMin = floor(min(bodePhase_y) - .001)
        limMax = ceil(max(bodePhase_y) + .001)
        ax2.set_ylim(limMin, limMax)
        ax2.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: '$%.0f^\circ$'%x))
        count = max(len(ax1.get_ygridlines()) - 1, 1)
        step = ceil((limMax - limMin) / count + .001)
        ax2.yaxis.set_ticks(numpy.arange(limMin, limMax + step, step))
    ax2.plot(bodePhase_x, bodePhase_y, '--', color='blue')

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
    freq = []
    mag = []
    phase = []
    x = fMin
    while x < fMax:
        f = x
        if not(angfreq):
            f *= 2 * pi.n()
        v = evalFn(h, I*f)
        m = phasorMag(v)
        p = phasorPhase(v)
        freq.append(f)
        mag.append(m)
        phase.append(p)
        x *= step
    return bodePlot_data(freq, mag, phase, **kwargs)
