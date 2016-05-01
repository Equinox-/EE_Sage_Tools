COMPLEX_PROPAGATION_CONSTANT(R,L,G,C, omega) = sqrt((R+L*i*omega) * (G+C*i*omega)) # henceforth known as cpc (small gamma)
PHASE_VELOCITY(cpc, omega) = omega / imag(cpc) # henceforth known as pv
ATTENUATION_CONST(cpc) = real(cpc) # alpha
PHASE_CONSTANT(cpc) = imag(cpc) # beta
CHARACTERISTIC_IMPEDANCE(R,L,G,C,omega) = sqrt((R+L*i*omega)/(G+C*i*omega)) # known as z0
# Assuming lossless line
REFLECTION_CONSTANT(Zo, Zl) = (Zl-Zo)/(Zl+Zo) # Known as rc (big gamma)
ZIN(Zo, Gamma, gamma, z) = Zo*((1+Gamma*exp(I*z*PHASE_CONSTANT(gamma)))/(1-Gamma*exp(I*z*PHASE_CONSTANT(gamma))))
STANDING_WAVE_RATIO(Gamma) = abs(Gamma+1)/abs(1-Gamma)