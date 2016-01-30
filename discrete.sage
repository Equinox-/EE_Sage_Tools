import numpy

def dstEval(fn, t):
    if isinstance(fn, sage.symbolic.function.Function):
        return fn(t)
    elif hasattr(fn, '__call__'):
        return fn({fn.default_variable():t})
    else:
        return fn

def dstFindNonZero(fn, center=0):
    if fn != 0: # Doesn't simplify to zero
        for n in xrange(100): # Scan region
            if dstEval(fn, center+n) != 0:
                return center+n
            elif dstEval(fn, center-n) != 0:
                return center-n
    return 0

def dstFitTime(fn, step, t=0, needy=25):
    pv = dstEval(fn, t)
    sn = 0
    while sn < needy:
        t = t + step
        cv = dstEval(fn, t)
        if cv==pv:
            sn = sn + 1
        else:
            sn = 0
        pv = cv
    return t - step * needy

class Discrete(object):
    @staticmethod
    def fromFunction(fn, left = None, right = None):
        center = None
        if left == None:
            if center == None:
                center = dstFindNonZero(fn)
            left = dstFitTime(fn, -1, center)
        if right == None:
            if center == None:
                center = dstFindNonZero(fn)
            right = dstFitTime(fn, 1, center) + 1
        return Discrete(numpy.array([dstEval(fn, t) for t in xrange(left, right)]), left)

    def __init__(self, vals, left):
        self.vals = numpy.array(vals)
        self.left = left

    def right(self):
        return self.left + len(self.vals)

    def xrange(self):
        return xrange(self.left, self.left + len(self.vals))

    def copy(self):
        return Discrete(numpy.copy(self.vals), self.left)

    def __len__(self):
        return len(self.vals)

    def __getitem__(self, key):
        if key < self.left:
            return self.vals[0]
        if key - self.left >= len(self.vals):
            return self.vals[len(self.vals) - 1]
        return self.vals[key - self.left]

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

    def impulses(self, **kwargs):
        g = Graphics()
        start = true
        for i in self.xrange():
            if abs(self[i]) > 1e-10:
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
            for k in self.xrange():
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