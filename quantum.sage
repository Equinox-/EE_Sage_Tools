## Quantum Mechanics Stuff

def quantumVars():
    G = globals()
    
    G['len_xyz'] = G['len_x'], G['len_y'], G['len_z'] = var('L_x L_y L_z', domain='positive')
    G['wv_xyz'] = G['wv_x'], G['wv_y'], G['wv_z'] = var('k_x k_y k_z', domain='real')
    G['n_xyz'] = G['n_x'], G['n_y'], G['n_z'] = var('n_x n_y n_z', domain='integer')
    G['m_xyz'] = G['m_x'], G['m_y'], G['m_z'] = var('n_x n_y n_z', domain='integer')


    G['planck_bar'] = var('hbar', domain='positive', latex_name='\hbar')
    G['planck'] = var('h', domain='positive', latex_name='h')
    G['electron_charge'] = G['e_0'] = G['q'] = var('e_0', domain='positive')
    G['boltzmann'] = var('k_b', domain='positive')
    G['x'], G['y'], G['z'] = var('x y z', domain='real')
    G['mass'] = var('m', domain='positive')
    G['mass_electron'] = var('m_e', domain='positive')
    G['mass_conduction'] = var('m_c', domain='positive')
    G['mass_valence'] = var('m_v', domain='positive')
    G['consta'] = {
        planck_bar: 1.0545718e-34,
        electron_charge: 1.60218e-19,
        boltzmann: 1.3806e-23,
        mass_electron: 9.1e-31
    }
    consta[planck] = consta[planck_bar] * 2 * pi

class NumericalHamiltonian:
    def __init__(self, minX, maxX, qV, subst):
        self.minX = minX
        self.maxX = maxX
        self.qV = qV
        self.subst = subst
        self.evs = None
        self.resolution(10)

    def resolution(self, incs):
        self.increments = incs
        self.dX = (self.maxX - self.minX) / self.increments
        self.evs = None

    def computeNumericalEigenfunctions(self):
        t = -planck_bar^2 / (2 * mass * self.dX * self.dX)
        H = [[]] * self.increments
        for r in xrange(0, self.increments):
            H[r] = [0] * self.increments
            for c in xrange(0, self.increments):
                if r == c:
                    H[r][c] = -2*t + qV(self.minX + r * self.dX)
                elif c == r-1 or c == r+1:
                    H[r][c] = t
                else:
                    H[r][c] = 0
        self.hamilt = matrix(H)
        self.evs = self.hamilt.subs(self.subst).change_ring(RDF).eigenvectors_left()
        return self.evs

    def sortAndSplitEigenfunctions(self):
        if self.evs == None:
            self.computeNumericalEigenfunctions()
        # unpack eigenvectors
        evsu = []
        for val,vecs,order in evs:
           for vec in vecs:
               vvec = vector(vec)
               vvec = vvec / sqrt(vvec.dot_product(vvec) * dX) # normalize
               vlvec = [((minX + (x*dX)).n(53), vvec[x].n(53)) for x in xrange(0, vvec.length())]
               evsu.append({"E": val, "chi": vvec, "chi_line": vlvec})
        evsu.sort(lambda x,y: int(-1) if x[0] < y[0] else int(1))
        return evsu