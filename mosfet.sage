vtn(vt0, gamma, phif, vsb)=vt0+gamma*(sqrt(vsb+phif)-sqrt(phif))

nmos_cutoff(vtnn, kn, lam, vgs, vds)=0
def nmos_is_cutoff(vtnn, kn, lam, vgs, vds):
    return vgs <= vtnn

nmos_triode(vtnn, kn, lam, vgs, vds)=kn*(vgs-vtnn-vds/2)*vds
def nmos_is_triode(vtnn, kn, lam, vgs, vds):
    if not(constIneq(vgs <= vtnn)) or not(constIneq((vgs - vtnn) > vds)):
        raise Exception("Can't simplify " + str(vgs) + " and " + str(vtnn) + " and " + str(vds))
    return (vgs > vtnn) and ((vgs - vtnn) > vds)

nmos_saturation(vtnn, kn, lam, vgs, vds)=(kn/2)*(vgs-vtnn)^2*(1+lam*vds)
def nmos_is_saturation(vtnn, kn, lam, vgs, vds):
    if not(constIneq(vgs <= vtnn)) or not(constIneq((vgs - vtnn) > vds)):
        raise Exception("Can't simplify: " + str(vgs) + " and " + str(vtnn) + " and " + str(vds))
    return (vgs > vtnn) and ((vgs - vtnn) <= vds)

def nmos(vtn, kn, lam, vgs, vds):
    if nmos_is_cutoff(vtn, kn, lam, vgs, vds):
        return nmos_cutoff(vtn, kn, lam, vgs, vds)
    if nmos_is_triode(vtn, kn, lam, vgs, vds):
        return nmos_triode(vtn, kn, lam, vgs, vds)
    else:
        return nmos_saturation(vtn, kn, lam, vgs, vds)

vtp(vt0, gamma, phif, vbs)=vt0-gamma*(sqrt(vbs+phif)-sqrt(phif))

pmos_cutoff(vtpp, kp, lam, vgs, vds)=0
def pmos_is_cutoff(vtpp, kp, lam, vgs, vds):
    return vgs >= vtpp
    
pmos_triode(vtpp, kp, lam, vgs, vds)=kp*(vgs-vtpp-vds/2)*vds
def pmos_is_triode(vtpp, kp, lam, vgs, vds):
    if not(constIneq(vgs >= vtpp)):
        raise Exception("Can't simplify " + str(vgs) + " and " + str(vtpp) + " and " + str(vds))
    if not(constIneq(abs(vgs-vtpp) >= abs(vds))):
        raise Exception("Can't simplify " + str(vgs) + " and " + str(vtpp) + " and " + str(vds))
    return (vgs < vtpp) and abs(vgs-vtpp) >= abs(vds)
    
pmos_saturation(vtpp, kp, lam, vgs, vds)=(kp/2)*(vgs-vtpp)^2*(1+lam*abs(vds))
def pmos_is_saturation(vtpp, kp, lam, vgs, vds):
    if not(constIneq(vgs >= vtpp)):
        raise Exception("Can't simplify " + str(vgs) + " and " + str(vtpp) + " and " + str(vds))
    if not(constIneq(abs(vgs-vtpp) >= abs(vds))):
        raise Exception("Can't simplify " + str(vgs) + " and " + str(vtpp) + " and " + str(vds))
    return (vgs < vtpp) and abs(vgs-vtpp) < abs(vds)

def pmos(vtp, kp, lam, vgs, vds):
    if pmos_is_cutoff(vtp, kp, lam, vgs, vds):
        return pmos_cutoff(vtp, kp, lam, vgs, vds)
    if pmos_is_triode(vtp, kp, lam, vgs, vds):
        return pmos_triode(vtp, kp, lam, vgs, vds)
    else:
        return pmos_saturation(vtp, kp, lam, vgs, vds)