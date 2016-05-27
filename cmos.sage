def cmos_solve_internal(kp, kn, lam_p, lam_n, vt_p, vt_n, vgs_p, vgs_n, vds_p, vds_n, freevar, space, sourced_current=0, regions=10):
    # Solve the 8 states:
    p_cut=pmos_cutoff(vt_p, kp, lam_p, vgs_p, vds_p)
    p_tri=pmos_triode(vt_p, kp, lam_p, vgs_p, vds_p)
    p_sat=pmos_saturation(vt_p, kp, lam_p, vgs_p, vds_p)
    n_cut=nmos_cutoff(vt_n, kn, lam_n, vgs_n, vds_n)
    n_tri=nmos_triode(vt_n, kn, lam_n, vgs_n, vds_n)
    n_sat=nmos_saturation(vt_n, kn, lam_n, vgs_n, vds_n)
    p_ops=[p_cut, p_sat, p_tri]
    p_assert=[pmos_is_cutoff, pmos_is_saturation, pmos_is_triode]
    n_ops=[n_cut, n_sat, n_tri]
    n_assert=[nmos_is_cutoff, nmos_is_saturation, nmos_is_triode]
    mmodes=["Cutoff", "Saturation", "Triode"]
    
    # vds_n - vds_p = vdd - vss
    for i in xrange(0, 3):
        i_p=p_ops[i]
        for j in xrange(0, 3):
            i_n=n_ops[j]
            regions=10
            step=(space[1] - space[0])/regions
            expr=sourced_current + i_n - i_p
            print(mmodes[i] + " and " + mmodes[j] + " find roots for " + str(expr) + " on " + str(space))
            for r in xrange(0, regions):
                try:
                    left=space[0] + step*r
                    right= left + step
                    s = {freevar:find_root(expr, left, right)}

                    p_good=p_assert[i](vt_p.subs(s), kp, lam_p, vgs_p.subs(s), vds_p.subs(s))
                    n_good=n_assert[j](vt_n.subs(s), kn, lam_n, vgs_n.subs(s), vds_n.subs(s))
                    state={"pmos":{"vds":vds_p.subs(s), "i":i_p.subs(s), "vgs":vgs_p.subs(s), "mode":mmodes[i]}, "nmos":{"vds":vds_n.subs(s), "i":i_n.subs(s), "vgs":vgs_n.subs(s), "mode":mmodes[j]}}
                    #print(str(state) + " is " + str("Yes" if p_good else "No") + " and " + str("Yes" if n_good else "No"))
                    if p_good and n_good:
                        return state
                except Exception as err:
                    pass
    return {"pmos":{"vds":0,"i":0,"vgs":0, "mode":"BAD"}, "nmos":{"vds":0,"i":0,"vgs":0, "mode":"BAD"}}    

def cmos_solve(kp, kn, lam_p, lam_n, vt_p, vt_n, vdd, vss, v_in, sourced_current=0, regions=10):
    vgs_p = v_in - vdd
    vgs_n = v_in - vss
    var('vds_p', domain='real')
    vds_n = vdd - vss + vds_p

    return cmos_solve_internal(kp, kn, lam_p, lam_n, vt_p, vt_n, vgs_p, vgs_n, vds_p, vds_n, vds_p, (-vdd, -vss), sourced_current, regions)

def cmos_switch_threshold(kp, kn, lam_p, lam_n, vt_p, vt_n, vdd, vss, sourced_current=0, regions=10):
    var('v_i', domain='real')
    vgs_n = vds_n = v_i - vss
    vgs_p = vds_p = v_i - vdd + vss

    return cmos_solve_internal(kp, kn, lam_p, lam_n, vt_p, vt_n, vgs_p, vgs_n, vds_p, vds_n, v_i, (vss, vdd), sourced_current, regions)

def find_root_safe(f, v_i, a, b, steps=5):
    if (steps == 0):
        return NaN
    try:
        cl=f.subs(v_i=a)
    except:
        cl=NaN
    try:
        cr=f.subs(v_i=b)
    except:
        cr=NaN
    if not(cl in RR) or not(cr in RR) or math.isnan(cl) or math.isnan(cr):
        mid=(a+b)/2
        left=find_root_safe(f, v_i, a, mid, steps-1)
        if (left in RR) and not(math.isnan(left)):
            return left
        right=find_root_safe(f, v_i, mid, b, steps-1)
        return right
    try:
        return find_root(f, a, b)
    except:
        return NaN

def cmos_crit_internal(sols, v_i, vss, vdd):
    # Scan through these, looking for the one with a negative slope at the middle:
    for s in sols:
        try:
            xx = find_root_safe(s[v_o].diff(v_i) + 1, v_i, vss, vdd)
            if (xx in RR) and not(math.isnan(xx)):
                return (xx, s[v_o](v_i=xx))
        except e:
            pass
    return (None, None)

def cmos_critical_points(kp, kn, lam_p, lam_n, vt_p, vt_n, vdd, vss):
    # V_IL, V_OH -> N saturated, P Triode  (Q1)
    # V_IH, V_OL -> P saturated, N Triode  (Q2)
    
    var('v_i', domain='real')
    var('v_o', domain='real')

    vgs_p = v_i - vdd
    vgs_n = v_i - vss
    vds_p = v_o - vdd
    vds_n = v_o - vss

    n_args=(vt_n, kn, lam_n, vgs_n, vds_n)
    p_args=(vt_p, kp, lam_p, vgs_p, vds_p)

    n_tri=nmos_triode(*n_args)
    n_sat=nmos_saturation(*n_args)
    p_tri=pmos_triode(*p_args)
    p_sat=pmos_saturation(*p_args)

    bound_low = vss+vt_n
    bound_high = vdd+vt_p

    Q1sol=solve(n_sat == p_tri, v_o, solution_dict=true)
    Q2sol=solve(p_sat == n_tri, v_o, solution_dict=true)

    Q1=cmos_crit_internal(Q1sol, v_i, bound_low, bound_high)
    Q2=cmos_crit_internal(Q2sol, v_i, bound_low, bound_high)
    return {"vil":Q1[0], "voh":Q1[1], "vih":Q2[0], "vol":Q2[1]}