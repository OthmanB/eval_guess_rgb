'''
	Main procedures to make the grids
'''
from io_tune import *
from model import model_asymptotic
import itertools
from scipy.io import readsav 

def save_cfg(rootdir, mcmc_file, s1_file, KIC, Dnu, epsilon, variables_name, variables, cfg_file='config.cfg'):
	'''
		A function that saves the information about the grid structure and processing 
		into a txt file inside the rootdir directory
	'''
	f=open(os.path.join(rootdir , cfg_file), 'w')
		f.write('MCMC_file=' + mcmc_file)
		f.write('s1_file=' +  s1_file)
		f.write('KIC=' + KIC)
		f.write('Dnu=' + str(Dnu))
		f.write('epsilon=' + str(epsilon))
		for i in range(len(variables_name)):
			f.write(v + ' = ' + variables_name[i])
			f.write('  ', str(variables[i]))
	f.close()


def make_specplot():
	'''
		Function that generate an image file with the plot of the fit on the top of the data
	'''
	return 0
def search_grid():
	'''
		A function that perform a search for the best chi22 into the grid
	'''

	#target_subdir_lintree(rootdir, subdir_seq, verbose=False)
	return 0


def likelihood(likelihood_name='chi22'):
	'''
		Function that handles the likelihood functions
	'''

	return l

def set_inputs(variables_name, variables):
	'''
		A function that regroups all the initial setups and handles the possibility
		that some critical input are not provided by setting some by defaut.
		Also properly format them as numpy arrays
	'''
	try:
		DP1_var=np.array(variables[variables_name.index('DP1')], dtype=np.float)
	except:
		## ---- Evaluation of DP if not a parameter ----
		# Super rough estimate derived by visual inspection of the Mosser+2015, Fig.1
		c=36.8222
		b=2.63897
		a=0.0168202
		DP1_0=a*Dnu**2 + b*Dnu + c
		DP1_var=np.array(np.linspace(0.8, 1.4, 10)*DP1_0, dtype=np.float) # 10 Values of DP1 by default

	# If gamma_max_l0 is specified use it
	try:
		gamma_max_l0_var=np.array(variables[variables_name.index('gamma_max_l0')], dtype=np.float)
	except: # Otherwise, empirically set some mode width using a combination of Dnu and numax
		coef_Dnu=0.01
		coef_numax=0.001
		gamma_max_l0_var=np.array([np.median([coef_Dnu*Dnu, coef_numax*numax])], dtype=np.float)
	try: # If Hsig is provided, use it, relative to the fiduciary value 3*Dnu
		Hsig_var=np.array(variables[variables_name.index('Hsig')], dtype=np.float)
		Hsig_var=3*Dnu*Hsig_var
	except: # Otherwise, just use the fiduciary value 3*Dnu
		Hsig_var=np.array([3*Dnu], dtype=np.float) # 1 value of Hsig by defautl
	
	Hmax_factor=np.array([1],dtype=np.float)

	try:
		rot_core_var=np.array(variables[variables_name.index('rot_core')], dtype=np.float)
	except:
		rot_core_var=np.array([0.1, 0.2, 0.4], dtype=np.float)

	try:
		rot_core_var=np.array(variables[variables_name.index('rot_envelope')], dtype=np.float)
	except:
		rot_envelope_var=np.array([0.1, 0.2, 0.4], dtype=np.float)

	try:
		inclination_var=np.array(variables[variables_name.index('inclination')], dtype=np.float)
	except:
		inclination_var=np.array([0, 45, 90], dtype=np.float) # Default is 3 inclination values

	return DP1_var, gamma_max_l0_var, Hmax_factor, Hsig_var, rot_core_var, rot_envelope_var, inclination_var

def grid_maker(rootdir, variables, variables_name, s1_file, mcmc_file):
	'''
		Core function that generates a grid into a specified directory
		rootdir: The root directory in which the subdirectories specifiying the 
			     variables of the model will be created
		vars: Values of the variables that need to be explored. Must be a list of list
			  The dimension of the main list must be equal to the dimension of vars_name
			  The sublists must contain the 'node values' at which the grid is calculated
			  eg. if you have vars_name=['p1', 'p2'] then vars=[[1,2,3,4], [-1,-2,-3]]
		vars_name: Name of the variables. Must be a 1D list.
	'''
	# Variable definition in function of the user choice
	DP1_var, gamma_max_l0_var, Hmax_factor, Hsig_var, rot_core_var, rot_envelope_var, inclination_var=set_inputs(variables_name, variables)

	# Read the MCMC file
	KIC, Dnu, C0, n_range, extra_priors, type_vals, eigen_vals, noise_guess, noise_fit=read_MCMCfile(file)
	epsilon=C0/Dnu-np.fix(C0/Dnu)
	fmin=n_range[0]*Dnu - C0
	fmax=n_range[1]*Dnu + C0

	# Read the s1 file
	s1_data=readsav(s1_file)
	signumax=s1_data['output'][-1]
	numax=s1_data['output'][-2]
	Hnumax=s1_data['output'][-3]
	
	# --- Defining spectrum data ---
	freq=output['freq']
	pos_keep=np.where(np.bitwise_and(freq >= fmin, freq <= fmax))
	freq=freq[pos_keep] # Limit the range for the data to where there is modes
	spec=output['spec_reg']
	spec=spec[pos_keep]/noise_harvey(freq, s1_data['output'][0:-4]) # Spectrum normalized by the noise... so that we do not need to worry about the noise background so much

	# Write the configuration used for the grid in a file before starting
	save_cfg(rootdir, mcmc_file, s1_file, KIC, Dnu, epsilon, variables_name, variables, cfg_file='config.cfg')

	# The main loop making the grid
	# NEED TO GET A WAY TO SCAN ALL THE COMBINATION HERE
	# THIS MIGHT IMPLIE MAKING A FULL LIST OF ALL POSSIBLE COMBINATIONS AND THEN LOOP OVER THEN USING A SINGLE LOOP INDEX
	for g in grid:
		param_vect=[]
		param_vect.append() # NEED TO GET A WAY TO SCAN ALL THE COMBINATION HERE
		make_subdir_lintree(rootdir, param_vect, verbose=False)
		
		maxHNR_l0=hmax_from_envelope(hnumax, Dnu, gamma_max_l0) # build the maxHNR from gamma_max_l0_var as those are closely related to the NRJ of the modes
		maxHNR_l0=Hmax_factor*maxHNR_l0

		model, nu_l0, nu_p_l1, nu_g_l1, nu_m_l1, nu_l2, nu_l3, 
		width_l0, width_m_l1, width_l2, width_l3, 
		height_l0, height_l1, height_l2, height_l3, 
		a1_l1, a1_l2, a1_l3=model_asymptotic(numax, Dnu, epsilon, DP1, beta_p_star=0., gamma_max_l0=gamma_max_l0, maxHNR_l0=maxHNR_l0, alpha_star=0., 
								rot_core=rot_core, rot_envelope=rot_envelope, inclination=inc, dfmin=6, dfmax=6, output_file_rot='test.rot')
		#l=likelihood(spec, likelihood_name='chi22')
		#make_specplot()
	return 0

def test_grid_maker(kic='2141436'):
	rootdir='/Volumes/home/2020/ML-Siddarth/tune-rgb/test/'
	s1_dir='/Volumes/home/2020/ML-Siddarth/s1/'
	mcmc_file_dir='/Volumes/home/2020/ML-Siddarth/setups/'

	# NEED TO DEFINE HERE THE GRID
	#variables_name=['DP1', 'gamma_max_l0', 'Hsig', 'rot_core', 'rot_envelope', 'inclination']
	variables_name=[] # ALL WILL BE DEFAULT VALUES
	variables=[] 
	
	# Execute the grid_maker
	s1_file= os.path.join(s1_dir , kic)
	mcmc_file= os.path.join(mcmc_file_dir , kic) #format_ID(kic)
	grid_maker(rootdir, variables, variables_name, s1_file, mcmc_file)
