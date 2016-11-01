def semiconductorVars():
	G = globals()
	G['out'] = var('out')
	G['planck_bar'] = 1.0545718e-34
	G['planck'] = planck_bar*2*pi
	G['electron_charge'] = -1.60218e-19
	G['electron_volt'] = -electron_charge
	G['electron_mass'] = 9.1e-31
	G['speed_of_light'] = 299792458
	G['boltzmann'] = 1.38064852e-23
	G['T'] = var('T')
	G['E_F'] = var('E_F')
	G['E'] = var('E')
	G['E_v'] = var('E_v')
	G['E_c'] = var('E_c')
	G['E_g'] = var('E_g')
	G['m_c'] = var('m_c')
	G['m_v'] = var('m_v')
	G['fermi_function'] = 1/(1 + exp((E - E_F) / (boltzmann * T)))
	G['DOS_c'] = m_c * sqrt(2 * m_c * (E - E_c)) / (pi^2 * planck_bar^3)
	G['DOS_v'] = m_v * sqrt(2 * m_v * (E_v - E)) / (pi^2 * planck_bar^3)
	G['prob_v'] = DOS_v * (1 - fermi_function)
	G['prob_c'] = DOS_c * fermi_function
	G['semiconductor_assume'] = [E_v + E_g == E_c]