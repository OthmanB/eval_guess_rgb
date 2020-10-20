import numpy as np
from bump_DP import make_synthetic_asymptotic_star
from function_rot import  amplitude_ratio

# Calculate the model for a given degree
# Elle tient compte des multiplets (splitting d'ordre 1)
# Elle tient compte de la force centrifuge (ordre 2) par le biais de a2=asym[0]= eta * Dnl ~ eta * 0.7 [no unit]
# Elle tient compte de l'effet de rotation latitudinal (ordre 3) par le biais de a3=asym[1] [microHz]
# A noter que a2 ~ a1^2/(G rho_sun) * Dnu_sun/Dnu [no unit]. Typiquement, rho=[0.1, 1.5] g/cm^3 pour une solar-like. G=6.67d-8 cm^3/g/s^2
def build_l_mode(x_l,Amp_l,f_c_l,a1,a2,a3, asym,gamma_l,l,V):

result=0

#ymax=max(Amp_l*V)
#plot,x_l,x_l,/nodata,background=fsc_color('white'),yr=[0,ymax],color=fsc_color('black'), $
#	xr=[f_c_l - 2*l*a1 - gamma_l, f_c_l + 2*l*a1 + gamma_l]
for m in range(-l,l+1):
	if l != 0 then:
		Qlm=1.*(l*(l+1.) - 3.*m^2)/((2.*l-1.)*(2.*l+3.)) #  coeficient for centrifugal force
		if l == 1:
			clm=m #   latitudinal coeficient for l=1
		if l == 2:
			clm=1.*(5*m^3 - 17*m)/3. #  latitudinal coeficient for l=2
		if l == 3:
			clm=0 # a3=asym[1] IS NOT IMPLEMENTED FOR l=3
		profile=2*( x_l- f_c_l*( 1d + a2*Qlm) + m*a1 + clm*a3) / gamma_l 
	else:
		profile=2*( x_l-f_c_l)/gamma_l

  profile=profile^2
  asymetry=(1d + asym*(x_l/f_c_l - 1e0))^2 + (0.5*gamma_l*asym/f_c_l)^2 # term of asymetry for each lorentzian

  #if abs(m) eq 0 then col=fsc_color('blue')
  #if abs(m) eq 1 then col=fsc_color('red')
  #if abs(m) eq 2 then col=fsc_color('Orange')
  #oplot, x_l,Amp_l*V[m+l]/(1+profile),color=col;,linestyle=2
  #plots, [f_c_l+m*a1, f_c_l+m*a1], [0,ymax], color=fsc_color('black'), linestyle=2
  result=result + asymetry*Amp_l*V[m+l]/(1e0+profile)

endfor

return result

def model_asymptotic(numax_star, Dnu_star, epsilon_star, beta_p_star=0., gamma_max_l0=1, alpha_star=0., rot_core=0.1, rot_envelope=0.1, inclination=45, dfmin=6, dfmax=6):
	'''
		numnax_star: numax of the star
		Dnu_star: Large separation Dnu of the star
		epsilon_star: epsilon of the star
		beta_p_star: curvature of the modes 
		gamma_max_l0: Width at numax for the l=0 modes
		alpha_star: Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core. Close to 0 if no nuclear reactions
		rot_core: rotation rate of the core in microHz
		rot_envelope: rotation rate of the envelope in microHz
		dfmin: Integer specifying how many radial orders below numax will be created
		dfmax: Integer specifying how many radial orders above numax will be created	
	'''
	
	# ---- Constants ----
	el=1 # Degree of the simulated mixed mode
	Vl=[1,1.5,0.6, 0.07]
	delta0l_percent=1
	# -------------------

	# Parameters for p modes that follow exactly the asymptotic relation of p modes
	delta0l_star=-el*(el + 1) * delta0l_percent / 100.

	## Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	#numax_star=numax_from_stello2009(Dnu_star)
	fmin=numax_star - dfmin*Dnu_star
	fmax=numax_star + dfmax*Dnu_star

	nmax_star=numax_star/Dnu_star - epsilon_star
	alpha_p_star=beta_p_star/nmax_star

	noise_params=[ 1., -2. , 0. , 1. ,-1. , 0. , 2.  ,1.]

	# ----- Setup to not touch ----
	Teff_star=-1
	rot_ratio=-1
	H0_spread=0
	output_file_rot='test.rot'
	#filetemplate="Configurations/templates/Sun.template"
	filetemplate="templates/11771760.template"
	# ----------------

	nu_l0, nu_p_l1, nu_g_l1, nu_m_l1, nu_l2, nu_l3, width_l0, width_m_l1, width_l2, width_l3, height_l0, height_l1, height_l2, height_l3, a1_l1, a1_l2, a1_l3=make_synthetic_asymptotic_star(Teff_star, numax_star, Dnu_star, epsilon_star, delta0l_percent, alpha_p_star, nmax_star, DP1_star, alpha_star, q_star, fmin, fmax, maxHNR_l0=maxHNR_l0, noise_params=noise_params, Gamma_max_l0=Gamma_max_l0, rot_env_input=rot_envelope, rot_ratio_input=rot_ratio, rot_core_input=rot_core, output_file_rot=output_file_rot, Vl=Vl, H0_spread=H0_spread, filetemplate=filetemplate)
	model=model_core_alt(x, nu_l0=nu_l0, nu_l1=nu_m_l1, nu_l2=nu_l2, nu_l3=nu_l3, 
				   height_l0=height_l0, height_l1=height_l1, height_l2=height_l2, height_l3=height_l3, 
				   width_l0=width_l0, width_l1=width_l1, width_l2=width_l2, width_l3=width_l3, 
				   a1_l1=a1_l1, a1_l2=a1_l1, a1_l3=a1_l1, inc=45)

def model_core_alt(x, nu_l0=[], nu_l1=[], nu_l2=[], nu_l3=[], 
				   height_l0=[], height_l1=[], height_l2=[], height_l3=[], 
				   width_l0=[], width_l1=[], width_l2=[], width_l3=[], 
				   a1_l1=[], a1_l2=[], a1_l3=[], inc=45):
	''' 
		This construct the model with very limited approximation consideration. 
		a2, a3 and the mode asymetry are set to 0. 
		Lorentzian parameters are regrouped by degree up to l=3, as per visible from the parameters above
	'''	
	a2=0
	a3=0
	asym=0

	model_l0=0
	model_l1=0
	model_l2=0
	model_l3=0

	el=0
	for en in range(len(nu_l0)):
		V=np.zeros(1) + 1.
		model_l0=model_l0 + build_l_mode(x, height_l0[en], nu_l0[en], 0, 0, 0, asym, width_l0[en], el, V) # H, f, a1, a2, a3, asym, Width, l, V
	el=1
	V=amplitude_ratio(el, inc) # l, inc
	for en in range(len(nu_l1)):
		model_l1=model_l1 + build_l_mode(x, height_l1[en], nu_l1[en], a1_l1[en], a2, a3, asym, width_l1[en], el, V) # H, f, a1, a2, a3, asym, Width, l, V
	el=2
	V=amplitude_ratio(el, inc) # l, inc
	for en in range(len(nu_l2)):
		model_l2=model_l2 + build_l_mode(x, height_l2[en], nu_l2[en], a1_l2[en], a2, a3, asym, width_l2[en], el, V) # H, f, a1, a2, a3, asym, Width, l, V
	el=3
	V=amplitude_ratio(el, inc) # l, inc
	for en in range(len(nu_l3)):
		model_l3=model_l3 + build_l_mode(x, height_l3[en], nu_l3[en], a1_l3[en], a2, a3, asym, width_l3[en], el, V) # H, f, a1, a2, a3, asym, Width, l, V

	# ----- degug ----
	plt.plot(x, model_l0)
	plt.plot(x, model_l1)
	plt.plot(x, model_l2)
	plt.plot(x, model_l3)
	plt.show()
	return model_l0+model_l1+model_l2+model_l3

def model_core(params, x):
	''' 
		This construct the model without any approximation consideration. It requires as many parameters as 
		Lorentzian. Thus params is a succession of N_Lparams=6 for N_Lorentzian 
		Each Lorentzian must have parameters in this order: l, Height, Freq, Width, a1, a2,a3, inclination, aymmetry
	'''

	N_Lparams=9 # CONSTANT DEFINING THE NUMBER OF PARAMETERS FOR EACH LORENTZIAN

	N_Lorentzian=len(params)/N_Lparams # Must be an integer
	i0=0
	model=0.
	for i in range(N_Lorentzian):
		p=params[i0:i0+N_Lparams]
		V=amplitude_ratio(p[0], p[8]) # l, inc
		model=model + build_l_mode(x, p[1], p[2], p[4], p[5], p[6], p[7], p[3], p[0], V) # H, f, a1, a2, a3, asym, Width, l, V
		i0=i0+N_Lparams+1
	return model
