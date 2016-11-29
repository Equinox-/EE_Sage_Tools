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
	G['epsilon_naut'] = G['epsilon_0'] = G['vacuum_permittivity'] = 8.854187817e-12
	G['epsilon'] = var('epsilon')
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
	G['EE'] = G['electric_field'] = var('EE', latex_name='\\mathcal{E}')
	G['Jdrift_n'] = G['drift_current_electrons'] = abs(electron_charge) * mobility_electrons * electrons * electric_field
	G['Jdrift_p'] = G['drift_current_holes'] = abs(electron_charge) * mobility_holes * holes * electric_field
	G['Jdiff_n'] = G['diffusion_current_electrons'] = electron_charge * diffusion_constant_electrons * -grad(electrons)
	G['Jdiff_p'] = G['diffusion_current_holes'] = -electron_charge * diffusion_constant_holes * -grad(holes)
	G['Jcurrent'] = G['total_current'] = Jdrift_n + Jdrift_p + Jdiff_n + Jdiff_p

	G['A'] = G['area'] = var('A')
	G['L_n'] = G['diffusion_length_n'] = var('L_n')
	G['L_p'] = G['diffusion_length_p'] = var('L_p')
	G['V_a'] = G['voltage_bias'] = var('V_a')
	G['J_ideal_diode'] = abs(electron_charge) * n_i^2 * (D_n / (N_A*L_n) + D_p / (N_D * L_p)) * (exp(abs(electron_charge) * V_a / (boltzmann*T)) - 1)
	G['I_ideal_diode'] = J_ideal_diode * A

	G['v_bi'] = G['built_in_voltage'] = var('V_bi')
	G['depletion_n_width'] = G['depletion_electron_width'] = sqrt(2*epsilon/ abs(electron_charge) * (N_A / (N_D * (N_A+N_D))) * v_bi)
	G['depletion_p_width'] = G['depletion_hole_width'] = sqrt(2*epsilon/ abs(electron_charge) * (N_D / (N_A * (N_A+N_D))) * v_bi)
	G['depletion_width'] = (depletion_electron_width + depletion_hole_width).simplify_full()
	G['v_bi_expr'] = G['built_in_voltage_expr'] = boltzmann * T / abs(electron_charge) * ln(N_A * N_D / (n_i ^ 2))


	for sfk in ['m', 'i', 's', 'o']:
		G['work_function_volt_' + sfk] = var('phi_' + sfk)
		G['work_function_' + sfk] = var('qphi_' + sfk)
		G['work_function_' + sfk + '_expr'] = abs(electron_charge) * G['work_function_volt_' + sfk]
		G['electron_affinity_volt_' + sfk] = var('chi_' + sfk)
		G['electron_affinity_' + sfk] = var('qchi_' + sfk)
		G['electron_affinity_' + sfk + '_expr'] = abs(electron_charge) * G['electron_affinity_volt_' + sfk]

	schottky_prefactor = ((m_c * abs(electron_charge) * boltzmann^2) / (2*pi^2 * planck_bar^3)) * T^2
	G['phi_b'] = G['schottky_barrier'] = work_function_volt_m - electron_affinity_volt_s
	G['qphi_b'] = G['electron_schottky_barrier'] = work_function_m - electron_affinity_s
	G['J_schottky_diode'] = schottky_prefactor * (exp(abs(electron_charge) * V_a / (boltzmann*T)) - 1) * exp(-electron_schottky_barrier / (boltzmann*T))
	G['I_schottky_diode'] = area * J_schottky_diode

	
	G['silicon_vals'] = { n_i: 10^10, epsilon: epsilon_naut * 11.68 }