def semiconductorVars(dimen = 1):
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
	G['N_D'] = G['num_donors'] = var('N_D')
	G['N_A'] = G['num_acceptors'] = var('N_A')
	G['n_i'] = var('n_i')
	G['fermi_function'] = 1/(1 + exp((E - E_F) / (boltzmann * T)))
	G['DOS_c'] = m_c * sqrt(2 * m_c * (E - E_c)) / (pi^2 * planck_bar^3)
	G['DOS_v'] = m_v * sqrt(2 * m_v * (E_v - E)) / (pi^2 * planck_bar^3)
	G['prob_v'] = DOS_v * (1 - fermi_function)
	G['prob_c'] = DOS_c * fermi_function
	G['semiconductor_assume'] = [E_v + E_g == E_c]

	posvar = [0] * dimen
	if dimen <= 3:
		posvar[0] = G['r_x'] = G['position_x'] = var('r_x')
		if dimen > 1:
			posvar[1] = G['r_y'] = G['position_y'] = var('r_y')
		if dimen > 2:
			posvar[2] = G['r_z'] = G['position_z'] = var('r_z')
	else:
		for j in xrange(1, dimen + 1):
			posvar[j-1] = G['e_' + str(j)] = G['position_' + str(j)] = var('e_' + str(j))
	if dimen == 1:
		G['grad'] = lambda x: diff(x, posvar[0])
	else:
		G['grad'] = lambda x: [diff(x, p) for p in posvar]

	G['t'] = G['time'] = var('t')
	G['n'] = G['electrons'] = function('n')(t, *posvar)
	G['p'] = G['holes'] = function('p')(t, *posvar)

	G['mu_n'] = G['mobility_electrons'] = var('mu_n')
	G['mu_p'] = G['mobility_holes'] = var('mu_p')

	G['D_n_expr'] = boltzmann * T / abs(electron_charge) * mu_n
	G['D_p_expr'] = boltzmann * T / abs(electron_charge) * mu_p

	G['D_n'] = G['diffusion_constant_electrons'] = var('D_n')
	G['D_p'] = G['diffusion_constant_holes'] = var('D_p')
	G['Jdiff_n'] = G['diffusion_current_electrons'] = electron_charge * diffusion_constant_electrons * -grad(electrons)
	G['Jdiff_p'] = G['diffusion_current_holes'] = -electron_charge * diffusion_constant_holes * -grad(holes)

	G['A'] = G['area'] = var('A')
	G['L_n'] = G['diffusion_length_n'] = var('L_n')
	G['L_p'] = G['diffusion_length_p'] = var('L_p')
	G['V_a'] = G['voltage_bias'] = var('V_a')
	G['J_ideal_diode'] = abs(electron_charge) * n_i^2 * (D_n / (N_A*L_n) + D_p / (N_D * L_p)) * (exp(abs(electron_charge) * V_a / (boltzmann*T)) - 1)
	G['I_ideal_diode'] = J_ideal_diode * A