def controlVars():
    G=globals()
    G['t']=var('t', domain='positive')
    G['s']=var('s')
    G['K']=var('K')
    G['B']=var('B')
    for i in range(1, 10):
        # Constants
        G['K_'+str(i)]=var('K_'+str(i))
        G['B_'+str(i)]=var('B_'+str(i))
        G['M_'+str(i)]=var('M_'+str(i))
        G['J_'+str(i)]=var('J_'+str(i))
        # Spatial Position
        G['x_'+str(i)]=xr=function('x_'+str(i))(t)
        G['dx_'+str(i)]=xr.diff(1)
        G['ddx_'+str(i)]=xr.diff(2)
        G['X_'+str(i)]=function('X_'+str(i))(s)
        # Spatial Force
        G['f_'+str(i)]=function('f_'+str(i))(t)
        G['F_'+str(i)]=function('F_'+str(i))(t)
        # Radial Position
        G['theta_'+str(i)]=tr=function('theta_'+str(i))(t)
        G['dtheta_'+str(i)]=tr.diff(1)
        G['ddtheta_'+str(i)]=tr.diff(2)
        G['Theta_'+str(i)]=function('Theta_'+str(i))(s)
        # Radial Force
        G['tau_'+str(i)]=function('tau_'+str(i), latex_name='\\tau_{' + str(i) + '}')(t)
        G['Tau_'+str(i)]=function('Tau_'+str(i), latex_name='\\tau_{' + str(i) + '}')(s)

def eom1Transfer(eom, xf, ff):
    tt = eom.right().collect(ff) / ff
    bb = eom.left().collect(xf) / xf
    return (tt, bb)

def eom1NTransfer(eom, xf, ff, ss):
    tt=eom1Transfer(eom, xf, ff)
    return normalizeTransfer(tt[0], tt[1], ss)

def normalizeTransfer(top, bottom, v):
    f=bottom.coefficients(v, sparse=False)[-1]
    ring=PolynomialRing(SR, v)
    return ring((top/f).coefficients(s, sparse=False)) / ring((bottom/f).coefficients(s, sparse=False))

def zerosInternal(f, v, **kwargs):
    numeric = kwargs['numeric'] if 'numeric' in kwargs else False
    if hasattr(f, 'roots'):
        try:
            rrr=f.roots()
            out = list()
            for k,v in rrr:
                for j in xrange(0, v):
                    if numeric:
                        k=n(k)
                    out.append(k)
            return out
        except:
            pass
    solus = solve(f==0, v, solution_dict=true)
    out = list()
    for s in solus:
        if numeric:
            out.append(n(s[v]))
        else:
            out.append(s[v])
    return out

def tsPoles(ts, vs=None):
    if vs == None:
        vs = ts.default_variable()
    return zerosInternal(ts.denominator(), vs)

def tsZeros(ts, vs=None):
    if vs == None:
        vs = ts.default_variable()
    return zerosInternal(ts.numerator(), vs)

def tsPoleZero(ts, **kwargs):
    label = kwargs.get('label', None)
    poles = tsPoles(ts)
    zeros = tsZeros(ts)
    plt = sage.plot.point.point2d(complexToXY(poles), marker='x', axes_labels=['$\Re$', '$\Im$'], **kwargs)
    plt += sage.plot.point.point2d(complexToXY(zeros), marker='o', **kwargs)
    if label != None:
        for p in poles:
            plt += sage.plot.text.text(label, (real(p), imag(p)), horizontal_alignment='left', vertical_alignment='top', **kwargs)
        for p in zeros:
            plt += sage.plot.text.text(label, (real(p), imag(p)), horizontal_alignment='left', vertical_alignment='top', **kwargs)
    return plt

def desolve_SS(tvar, U, A, B, C, D, **kwargs):
    steps = kwargs['step'] if 'step' in kwargs else 200
    end = kwargs['end'] if 'end' in kwargs else 1
    start = kwargs['start'] if 'start' in kwargs else 0
    ics = kwargs['ics'] if 'ics' in kwargs else vector([0] * len(A.rows()))
    split = kwargs['split'] if 'split' in kwargs else False

    step = (end - start) / steps
    xv = ics
    tv = start
    tout = []
    out = []
    while True:
        dt = (end - tv) if ((tv+step) > end) else step
        Uv = U.subs({tvar: tv})
        out.append(C * xv + D * Uv)
        tout.append(tv)
        if tv >= end:
            break

        dxdt = A * xv + B * Uv
        xv += dxdt * dt
        tv += dt
    if split:
        outs = []
        for k in range(0, len(C.rows())):
            pts = []
            for i in range(0, len(out)):
                pts.append((tout[i], out[i][k]))
            outs.append(pts)
        return outs
    return tout, out

def evans(expr, kMax=None, accuracy=1):
    import itertools

    # limits
    smax=.002/accuracy
    smin=smax/3
    stepMin = 1e-15/accuracy
    iterMax=20000

    va = expr.default_variable()
    n = expr.numerator()
    d = expr.denominator()
    if kMax == None:
        kMax = round(500 * vector(SR(d).coefficients(s, sparse=False)).norm() / vector(SR(n).coefficients(s, sparse=False)).norm())

    vzeros = zerosInternal(n, va)
    vpoles = zerosInternal(d, va)
    if len(vzeros) + len(vpoles) == 0:
        raise ValueError("No poles or zeros")

    nrm = 2.0 * max([abs(v) for v in vzeros] + [abs(v) for v in vpoles])
    md = d.degree(va)

    kv = 1e-4
    step = 1
    itr = 0

    routeCount = 0
    poleRoutes = {}
    # initialize poleRoutes @ k=0
    extent = 0

    kPoles = zerosInternal(d, va, numeric=True)
    for v in kPoles:
        extent = max(extent, v)
        poleRoutes[routeCount] = [v]
        routeCount += 1

    while kv <= kMax and itr < iterMax:
        kPoles = zerosInternal(d+kv*n, va, numeric=True)
        # map poles to routes
        kRoutes = None
        kRouteErr = None
        # brute force the best match every time
        for perm in itertools.permutations(range(0, len(kPoles))):
            kTmpRoutes = {}
            kTmpRouteErr = 0
            for i in xrange(0, len(perm)):
                j = perm[i]
                # pole i maps to route j
                if j in poleRoutes:
                    # raw distance
                    kMyErr = abs(poleRoutes[j][-1] - kPoles[i])

                    # dist of midpoint from line connecting self and self-2
                    if len(poleRoutes[j]) >= 2:
                        pt1 = kPoles[i]
                        pt2 = poleRoutes[j][-2]
                        test = poleRoutes[j][-1]
                        dx=real(pt2)-real(pt1)
                        dy=imag(pt2)-imag(pt1)
                        
                        dxt=real(pt1)-real(test)
                        dyt=imag(pt1)-imag(test)

                        midLineDistance = abs(dx * dyt - dy * dxt) / kMyErr
                        kMyErr = (midLineDistance * 1e3) + (kMyErr / 3e1)
                    kTmpRouteErr = max(kTmpRouteErr, kMyErr)
                kTmpRoutes[j] = kPoles[i]
            if kRoutes == None or kTmpRouteErr < kRouteErr:
                kRouteErr = kTmpRouteErr
                kRoutes = kTmpRoutes
        kRouteErr /= nrm

        #print(str(kv) + "/" + str(kMax) + "\t\t" + str(kRouteErr) + "\t" + str(step) + "\t" + str([poleRoutes[k][-1] for k in poleRoutes]))
        # check err
        if kRouteErr > smax and step > stepMin:
            step = max(stepMin, step / 2e0)
        else:
            if kRouteErr < smin:
                step *= 2e0
            for k in kRoutes:
                v=kRoutes[k]
                extent = max(extent, abs(v))
                routeCount = max(routeCount, k+1)
                if k in poleRoutes:
                    poleRoutes[k].append(v)
                else:
                    poleRoutes[k] = [v]
            kv += step
        itr = itr + 1

    g = tsPoleZero(n/d)
    for k in poleRoutes:
        g += line(complexToXY(poleRoutes[k]))

    # draw asmyptotes
    asymLen = extent * 1.4
    count = len(vpoles) - len(vzeros)
    if count > 0:
        sigma = (sum(vpoles) - sum(vzeros)) / count
        for pole in xrange(0, count):
            theta = pi * (1 + 2*pole) / count
            g += line([(real(sigma), imag(sigma)), (real(sigma) + cos(theta) * asymLen, imag(sigma) + sin(theta) * asymLen)], linestyle='--', color='red')

    return g

class ControlSystem:
    def __init__(self, tv, sv, *args):
        self.tv = tv
        self.sv = sv
        self.things=list()
        # ground
        self.things.append({'m':infinity, 'l':0, 'x':0, 'k':list(), 'b':list(), 'f':list()})
        for count, thing in enumerate(args):
            d = dict()
            d['m'] = thing[2]
            d['l'] = thing[1]
            d['x'] = thing[0]
            d['k'] = list()
            d['b'] = list()
            d['f'] = list()
            self.things.append(d)

    def k(self, v, a, b):
        self.things[a]['k'].append((v, b))
        self.things[b]['k'].append((v, a))
        return self

    def b(self, v, a, b):
        self.things[a]['b'].append((v, b))
        self.things[b]['b'].append((v, a))
        return self

    def f(self, a, vV, vL):
        self.things[a]['f'].append((vV, vL))
        return self

    def eom(self, i):
        tt = self.things[i]
        expr = tt['m'] * diff(tt['x'], 2)
        for b in tt['b']:
            expr += b[0] * (diff(tt['x'], 1) - diff(self.things[b[1]]['x'], 1))
        for k in tt['k']:
            expr += k[0] * (tt['x'] - self.things[k[1]]['x'])

        ff = 0
        for f in tt['f']:
            ff += f[0]
        if isinstance(ff, sage.symbolic.expression.Expression):
            ff = ff.simplify_full()
        if isinstance(expr, sage.symbolic.expression.Expression):
            expr = expr.simplify_full()
        return expr == ff

    def eom_LT_Vec(self, i):
        row = [0] * (len(self.things) + 1)
        s = self.sv
        tt = self.things[i]
        row[i] += tt['m'] * s^2
        for b in tt['b']:
            row[i] += b[0] * s
            row[b[1]] -= b[0] * s
        for k in tt['k']:
            row[i] += k[0]
            row[k[1]] -= k[0]
        for f in tt['f']:
            row[len(self.things)] += f[1]
        for j in range(0, len(row)):
            if isinstance(row[j], sage.symbolic.expression.Expression):
                row[j] = row[j].simplify_full()
        return (vector(row[1:len(self.things)]), row[len(self.things)])

    def eom_LT(self, i):
        vec, force = self.eom_LT_Vec(i)
        expr = 0
        for j in range(1, len(self.things)):
            expr += vec[j-1] * self.things[j]['l']
        return expr == force

    def eom_LT_Mat(self):
        rows = list()
        ans = list()
        for j in range(1, len(self.things)):
            vec, force = self.eom_LT_Vec(j)
            rows.append(vec)
            ans.append(force)
        return matrix(rows), vector(ans)

    def transfer_Mat(self):
        mat, ans = self.eom_LT_Mat()
        imat = mat.inverse()
        # mat * VV * imat = ans * imat
        return imat, ans

    def transfer_Vec(self, i):
        imat, ans = self.transfer_Mat()
        return vector(imat[i-1]), ans

    def transfer_Expr(self, i):
        ivec, ans = self.transfer_Vec(i)
        return ivec * ans

    def transfer(self, i, force):
        return (self.transfer_Expr(i) / force).simplify_full()

    def transfer_Normalized(self, i, force):
        tf = self.transfer(i, force)
        return normalizeTransfer(tf.numerator(), tf.denominator(), self.sv)

    def eom(self, i):
        tt = self.things[i]
        expr = tt['m'] * diff(tt['x'], 2)
        for b in tt['b']:
            expr += b[0] * (diff(tt['x'], 1) - diff(self.things[b[1]]['x'], 1))
        for k in tt['k']:
            expr += k[0] * (tt['x'] - self.things[k[1]]['x'])

        ff = 0
        for f in tt['f']:
            ff += f[0]
        if isinstance(ff, sage.symbolic.expression.Expression):
            ff = ff.simplify_full()
        if isinstance(expr, sage.symbolic.expression.Expression):
            expr = expr.simplify_full()
        return expr == ff

    def eom_SS(self):
        X = []
        for j in range(1, len(self.things)):
            X.append(self.things[j]['x'])
            X.append(diff(self.things[j]['x'], 1))
        FV = dict()
        for j in range(1, len(self.things)):
            tt = self.things[j]
            for f in tt['f']:
                if not(f in FV):
                    FV[f[0]] = len(FV)

        mat_A = []
        mat_B = []
        for j in range(1, len(self.things)):
            tt = self.things[j]
            row = [0] * len(X)
            ff = [0] * len(FV)
            m = tt['m']
            for b in tt['b']:
                row[1+(j-1)*2] -= b[0] / m
                if b[1] > 0:
                    row[1+(b[1]-1)*2] += b[0] / m
            for k in tt['k']:
                row[(j-1)*2] -= k[0] / m
                if k[1] > 0:
                    row[(k[1]-1)*2] += k[0] / m
            for f in tt['f']:
                ff[FV[f[0]]] += (1/m)
            row2=[0] * len(X)
            row2[1+(j-1)*2]=1
            mat_A.append(row2)
            mat_A.append(row)
            mat_B.append([0] * len(FV))
            mat_B.append(ff)
        F = [0] * len(FV)
        for k, v in FV.iteritems():
            F[v] = k
        return (matrix(mat_A), vector(X), matrix(mat_B), vector(F))