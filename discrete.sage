import numpy

def dstFindNonZero(fn, center=0):
    if fn != 0: # Doesn't simplify to zero
        for n in xrange(100): # Scan region
            if evalFn(fn, center+n) != 0:
                return center+n
            elif evalFn(fn, center-n) != 0:
                return center-n
    return 0

def dstFitTime(fn, step, t=0, needy=25):
    pv = evalFn(fn, t)
    sn = 0
    while sn < needy:
        t = t + step
        cv = evalFn(fn, t)
        if cv==pv:
            sn = sn + 1
        else:
            sn = 0
        pv = cv
    return t - step * needy

class Discrete(object):
    @staticmethod
    def fromFunction(fn, left = None, right = None, periodic = False):
        center = None
        if left == None:
            if center == None:
                center = dstFindNonZero(fn)
            left = dstFitTime(fn, -1, center)
        if right == None:
            if center == None:
                center = dstFindNonZero(fn)
            right = dstFitTime(fn, 1, center) + 1
        return Discrete(numpy.array([evalFn(fn, t) for t in xrange(left, right)]), left, periodic)

    def __init__(self, vals, left, periodic = False):
        self.vals = numpy.array(vals)
        self.left = left
        self.periodic = periodic

    def right(self):
        return self.left + len(self.vals)

    def xrange(self):
        return xrange(self.left, self.right())

    def copy(self):
        return Discrete(numpy.copy(self.vals), self.left)

    def __len__(self):
        return len(self.vals)

    def __getitem__(self, key):
        if not(self.periodic):
            if key < self.left:
                return self.vals[0]
            if key - self.left >= len(self.vals):
                return self.vals[len(self.vals) - 1]
        return self.vals[(key - self.left) % len(self)]

    def __setitem__(self, key, value):
        self.vals[key - self.left] = value

    def __repr__(self):
        return "Discrete[left=" + str(self.left) + ", " + str(self.vals) + "]"

    def __elementBinaryDD(self, other, fn):
        nl = min(self.left, other.left)
        nr = max(self.right(), other.right())
        nv = Discrete(numpy.zeros(nr - nl), nl)
        for v in xrange(nl, nr):
            nv[v] = fn(self[v], other[v])
        return nv

    def reverse(self):
        return Discrete(self.vals[::-1], -self.right())

    def __abs__(self):
        return Discrete(abs(self.vals), self.left)

    def __neg__(self):
        return Discrete(-self.vals, self.left)

    def __pos__(self):
        return self.copy()

    def __lshift__(self, other):
        return Discrete(numpy.copy(self.vals), self.left - other)

    def __rshift__(self, other):
        return Discrete(numpy.copy(self.vals), self.left + other)

    def __mul__(self, other):
        if isinstance(other, Discrete):
            return self.__elementBinaryDD(other, lambda x,y: x*y)
        else:
            return Discrete(self.vals * other, self.left)

    def __pow__(self, other):
        if isinstance(other, Discrete):
            return self.__elementBinaryDD(other, lambda x,y: pow(x,y))
        else:
            return Discrete(pow(self.vals, other), self.left)

    def __div__(self, other):
        if isinstance(other, Discrete):
            return self.__elementBinaryDD(other, lambda x,y: x/y)
        else:
            return Discrete(self.vals / other, self.left)

    def __add__(self, other):
        if isinstance(other, Discrete):
            return self.__elementBinaryDD(other, lambda x,y: x+y)
        else:
            return Discrete(self.vals + other, self.left)

    def __sub__(self, other):
        if isinstance(other, Discrete):
            return self.__elementBinaryDD(other, lambda x,y: x-y)
        else:
            return Discrete(self.vals - other, self.left)

    def impulses(self, left = None, right = None, **kwargs):
        g = Graphics()
        start = true
        expl = True
        expl = True
        if left == None:
            left = self.left
            expl = False
        if right == None:
            right = self.right()
            expl = False
        for i in xrange(left, right):
            if abs(self[i]) > 1e-10 or expl or (i+1 < self.right() and abs(self[i+1]) > 1e-10) or (i > self.left and abs(self[i-1]) > 1e-10):
                if start:
                    g += arrow2d((i, 0), (i, self[i]), kwargs)
                    start = false
                else:
                    g += arrow2d((i, 0), (i, self[i]))
        return g

    def points(self):
        pts = []
        for i in self.xrange():
            pts.append((i, self[i]))
        return pts

    def convolve(self, other):
        left = self.left + other.left
        right = self.right() + other.right()
        out = Discrete(numpy.zeros(right - left), left)
        for t in xrange(left, right):
            for k in xrange(self.left - len(self), self.right() + len(self)):
                out[t] += self[k] * other[t - k]
        return out

    def nMax(self):
        bi = 0
        for i in self.xrange():
            if self[i] >= self[bi]:
                bi = i
        return bi

    def nMin(self):
        bi = 0
        for i in self.xrange():
            if self[i] <= self[bi]:
                bi = i
        return bi

    def find(self, val, center=0, step=0):
        for i in xrange(0, 100):
            if step == 0:
                if self[center + i] == val:
                    return center + i
                elif self[center - i] == val:
                    return center - i
            elif self[center + step * i] == val:
                return center + step * i
        return None

    def symbolic_steps(self, v):
        eps = 1e-10
        fn = self.vals[0]
        ln = self.vals[0]
        for t in self.xrange():
            cn = self[t]
            d = cn - ln
            if abs(d) > eps:
                if abs(d - 1) < eps:
                    fn += unit_step(v - t)
                elif abs(d + 1) < eps:
                    fn -= unit_step(v - t)
                else:
                    fn += (cn - ln) * unit_step(v - t)
                ln = cn
        return fn

    def fourier_series(self):
        N = numpy_to_sage(len(self))
        out = [ZZ(0)]*N
        left = numpy_to_sage(self.left)
        for k in xrange(0, N):
            kt = numpy_to_sage(k)
            for n in xrange(0, N):
                nt = numpy_to_sage(n)
                out[k] += numpy_to_sage(self.vals[n]) * exp(-I*kt*2*pi*nt/N)
            out[k] *= (1/N) * exp(-I*kt*2*pi*left/N)
        return out

    def fourier_series_symbolic(self, var):
        return fourier_series_symbolic(self.fourier_series(), len(self), var)
