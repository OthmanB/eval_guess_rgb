'''
	Main procedures to make the grids
'''
from io_tune import *
from model import *
import itertools
from scipy.io import readsav 
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import os

def floatvect2str(vect, precision='{:6f}', list_format=True):
	'''
		Convert a vector/list of floats into a fixed precision string
	'''
	s=""
	for v in vect:
		s=s + " " + precision.format(v)
	if list_format == True:
		s=s.split()
	return s

def get_combi(lists_of_list): #, dtype=np.float):
	'''
		Give all the possible combinations between the N lists inside the list
		e.g. lists_of_list=[[1,2], [0,10], [100,200]] returns
			 r=[(1, 0, 100), (1, 0, 200), (1, 10, 100), (1, 10, 200), (2, 0, 100), (2, 0, 200), (2, 10, 100), (2, 10, 200)]
	'''
	r=list(itertools.product(*lists_of_list))
	return r #np.array(r, dtype=dtype)


def save_cfg(rootdir, mcmc_file, s1_file, KIC, Dnu, epsilon, variables_name, variables, cfg_file='config.cfg'):
	'''
		A function that saves the information about the grid structure and processing 
		into a txt file inside the rootdir directory
	'''
	f=open(os.path.join(rootdir , cfg_file), 'w')
	f.write('MCMC_file=' + mcmc_file + '\n')
	f.write('s1_file=' +  s1_file + '\n')
	f.write('KIC=' + KIC + '\n')
	f.write('Dnu=' + str(Dnu) + '\n')
	f.write('epsilon=' + str(epsilon) + '\n')
	for i in range(len(variables_name)):
		f.write(variables_name[i] + ' = ' +  str(variables[i]) + '\n')
	f.close()

def make_specplot(freq, spec, model, file_out, scoef=[], fmin=-1, fmax=-1):
	'''
		Function that generate an image file with the plot of the fit on the top of the data
	'''
	if fmin <0:
		fmin=min(freq)
	if fmax <0:
		fmax=max(freq)
	
	col=['Blue', 'Orange', 'Cyan', 'Olive' ]
	if len(scoef) > 4:
		print('Warning: requesting too many smoothed curved. The maximum is 4.')
		print('         We will truncate the scoef vector for compatibility.')
		print('         Pursuing...')
		scoef=scoef[:4]

	ymin=0 # To change only if plots are in ylog
	ymax=max(spec[np.where(np.bitwise_and(freq < fmax, freq > fmin))])

	resol=freq[10] - freq[9]

	plt.figure()
	plt.plot(freq, spec, color='Black')
	if len(scoef) > 0:
		i=0
		ymax=[]
		for s in scoef:
			sfactor=np.int(s/resol) # Smooth over the specified width in microHz
			y_smooth=gaussian_filter(spec, s, mode='mirror')
			plt.plot(freq, y_smooth, color=col[i])
			ymax.append(max(y_smooth[np.where(np.bitwise_and(freq < fmax, freq > fmin))]))
			i=i+1
		ymax=max(ymax)
	plt.plot(freq, model, color='Red')

	plt.xlabel('Frequency (microHz)')
	plt.ylabel('Power (ppm^2/microHz)')
	plt.xlim(fmin, fmax)
	plt.ylim(ymin, ymax)
	plt.savefig(file_out, dpi=300)

	return 0

def make_specplot_v2(freq, spec, model_l0, model_l1, model_l2, model_l3, Delta, file_out_core, scoef=[], fmin=-1, fmax=-1, col_spec=['lightgrey', 'gray', 'darkgrey','black'], likelihood=None):
	'''
		Function that generate an image file with the plot of the fit on the top of the data
		More elaborated than make_specplot()
	'''
	if fmin <0:
		fmin=min(freq)
	if fmax <0:
		fmax=max(freq)

	fmax_limit=fmax
		
	if len(scoef) > 3:
		print('Warning: requesting too many smoothed curved. The maximum is 3.')
		print('         We will truncate the scoef vector for compatibility.')
		print('         Pursuing...')
		scoef=scoef[:4]

	ymin=0 # To change only if plots are in ylog
	ymax=max(spec[np.where(np.bitwise_and(freq < fmax, freq > fmin))])

	resol=freq[10] - freq[9]
	print('resol =', resol)

	col_model=['blue', 'red', 'orange', 'olive']


	# ------ Global view ------
	plt.figure()
	if len(scoef) > 0:
		i=0
		ymax=[]
		for s in scoef:
			sfactor=np.int(s/resol) # Smooth over the specified width in microHz
			y_smooth=gaussian_filter(spec, s, mode='mirror')
			plt.plot(freq, y_smooth, color=col_spec[i])
			ymax.append(max(y_smooth[np.where(np.bitwise_and(freq < fmax, freq > fmin))]))
			i=i+1
		ymax=1.1*max(ymax)
	
	plt.plot(freq, model_l0, color=col_model[0])
	plt.plot(freq, model_l1, color=col_model[1])
	plt.plot(freq, model_l2, color=col_model[2])
	plt.plot(freq, model_l3, color=col_model[3])
	plt.xlabel('Frequency (microHz)')
	plt.ylabel('Power (ppm^2/microHz)')
	plt.xlim(fmin, fmax)
	plt.ylim(ymin, ymax)
	plt.savefig(file_out_core + '.png', dpi=300)


	# ------ Local views -----	
	Delta_overlap=0.05
	slice_n=1
	fmax=fmin + Delta + Delta_overlap*Delta
	while fmax < fmax_limit:
	#for slice_n in range(Nslices):
		if fmax > fmax_limit:
			fmax=fmax_limit
		#print('fmin/fmax:', fmin, '  /  ', fmax)

		plt.figure()
		if len(scoef) > 0:
			i=0
			ymax=[]
			for s in scoef:
				sfactor=np.int(s/resol) # Smooth over the specified width in microHz
				y_smooth=gaussian_filter(spec, s, mode='mirror')
				plt.plot(freq, y_smooth, color=col_spec[i])
				ymax.append(max(y_smooth[np.where(np.bitwise_and(freq < fmax, freq > fmin))]))
				i=i+1
			ymax=1.1*max(ymax)
	
		plt.plot(freq, model_l0, color=col_model[0])
		plt.plot(freq, model_l1, color=col_model[1])
		plt.plot(freq, model_l2, color=col_model[2])
		plt.plot(freq, model_l3, color=col_model[3])

		plt.xlabel('Frequency (microHz)')
		plt.ylabel('Power (ppm^2/microHz)')
		plt.xlim(fmin, fmax)
		plt.ylim(ymin, ymax)
		plt.savefig(file_out_core + '_' + str(slice_n) + '.png', dpi=300)
		slice_n=slice_n+1

		fmin=fmin + Delta - Delta_overlap*Delta
		fmax=fmin + Delta + Delta_overlap*Delta

	return 0
def search_grid():
	'''
		A function that perform a search for the best chi22 into the grid
	'''

	#target_subdir_lintree(rootdir, subdir_seq, verbose=False)
	return 0


def likelihood(y, model, likelihood_name='chi22'):
	'''
		Function that handles the likelihood functions
	'''
	if likelihood_name == 'chi22':
		l = -np.sum(y/model)-np.sum(np.log(model))
	else:
		print('Error: The program only supports likelihood_name = "chi22"')
	return l

def set_inputs(variables_name, variables, Dnu, numax):
	'''
		A function that regroups all the initial setups and handles the possibility
		that some critical input are not provided by setting some by defaut.
		Also properly format them as numpy arrays
	'''
	vnames=[]

	v='DP1'
	status=variables_name.count(v)
	if status == True:
		DP1_var=np.array(variables[variables_name.index(v)], dtype=np.float)
		vnames.append(variables_name[variables_name.index(v)]) # This in order to organise in specific order
	else:
		## ---- Evaluation of DP if not a parameter ----
		# Super rough estimate derived by visual inspection of the Mosser+2015, Fig.1
		c=36.8222
		b=2.63897
		a=0.0168202
		DP1_0=a*Dnu**2 + b*Dnu + c
		DP1_var=np.array(np.linspace(0.8, 1.4, 10)*DP1_0, dtype=np.float) # 10 Values of DP1 by default
		vnames.append(v)

	# If gamma_max_l0 is specified use it
	v='gamma_max_l0'
	status=variables_name.count(v)
	if status == True:
		gamma_max_l0_var=np.array(variables[variables_name.index(v)], dtype=np.float)
		vnames.append(variables_name[variables_name.index(v)])
	else: # Otherwise, empirically set some mode width using a combination of Dnu and numax
		coef_Dnu=0.01
		coef_numax=0.001
		gamma_max_l0_var=np.array([np.median([coef_Dnu*Dnu, coef_numax*numax])], dtype=np.float)
		vnames.append(v)
		
	''' # NOT USED AT THE MOMENT
	v='Hsig'
	status=variables_name.count(v)
	if status == True: # If Hsig is provided, use it, relative to the fiduciary value 3*Dnu
		Hsig_var=np.array(variables[variables_name.index(v)], dtype=np.float)
		Hsig_var=3*Dnu*Hsig_var
		vnames.append(variables_name[variables_name.index(v)])
	else: # Otherwise, just use the fiduciary value 3*Dnu
		Hsig_var=np.array([3*Dnu], dtype=np.float) # 1 value of Hsig by default
		vnames.append(v)
	'''


	#Hmax_factor=np.array([1],dtype=np.float)
	#vnames.append('Hmax_factor')
	v='Hmax_factor'
	status=variables_name.count(v)
	if status == True:
		Hmax_factor=np.array(variables[variables_name.index(v)], dtype=np.float)
		vnames.append(variables_name[variables_name.index(v)])
	else:
		Hmax_factor=np.array([1], dtype=np.float)
		vnames.append(v)
		
	v='rot_core'
	status=variables_name.count(v)
	if status == True:
		rot_core_var=np.array(variables[variables_name.index(v)], dtype=np.float)
		vnames.append(variables_name[variables_name.index(v)])
	else:
		rot_core_var=np.array([0.1, 0.2, 0.4], dtype=np.float)
		vnames.append(v)
	
	v='rot_envelope'
	status=variables_name.count(v)
	if status == True:
		rot_envelope_var=np.array(variables[variables_name.index(v)], dtype=np.float)
		vnames.append(variables_name[variables_name.index(v)])
	else:
		rot_envelope_var=np.array([0.1, 0.2, 0.4], dtype=np.float)
		vnames.append(v)

	v='inclination'
	status=variables_name.count(v)
	if status == True:
		inclination_var=np.array(variables[variables_name.index(v)], dtype=np.float)
		vnames.append(variables_name[variables_name.index(v)])
	else:
		inclination_var=np.array([0, 45, 90], dtype=np.float) # Default is 3 inclination values
		vnames.append(v)

	v='d0l_percent'
	status=variables_name.count(v)
	if status == True:
		d0l_var=np.array(variables[variables_name.index(v)], dtype=np.float)
		vnames.append(variables_name[variables_name.index(v)])
	else:
		d0l_var=np.array([1.5], dtype=np.float) # Default is 3 inclination values
		vnames.append(v)


	#all_vectors=[DP1_var, gamma_max_l0_var, Hmax_factor, Hsig_var, rot_core_var, rot_envelope_var, inclination_var]
	all_vectors=[DP1_var, gamma_max_l0_var, Hmax_factor, rot_core_var, rot_envelope_var, inclination_var,  d0l_var]
	return vnames, all_vectors

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
	# Read the MCMC file
	print('  [a] Reading the MCMC file...')
	KIC, Dnu, C0, n_range, extra_priors, type_vals, eigen_vals, noise_guess, noise_fit=read_MCMCfile(mcmc_file)
	epsilon=C0/Dnu-np.fix(C0/Dnu)
	fmin=n_range[0]*Dnu - C0
	fmax=n_range[1]*Dnu + C0

	# Read the s1 file
	print('  [b] Reading the s1 file...')
	s1_data=readsav(s1_file)
	signumax=s1_data['output'][-1]
	numax=s1_data['output'][-2]
	Hnumax=s1_data['output'][-3]
	Hnumax=Hnumax/noise_harvey(np.zeros(1) + numax, s1_data['output'][0:-4]) # We need to remove the noise background

	# Variable definition in function of the user choice. Also sorting to specific order
	print('  [c] Defining/Sorting variables...')
	variables_name, variables=set_inputs(variables_name, variables, Dnu, numax)
	
	# --- Defining spectrum data ---
	print('  [d] Setting/Normalising the spectrum...')
	freq=s1_data['freq']
	pos_keep=np.where(np.bitwise_and(freq >= fmin, freq <= fmax))
	freq=freq[pos_keep] # Limit the range for the data to where there is modes
	spec=s1_data['spec_reg']
	spec=spec[pos_keep]/noise_harvey(freq, s1_data['output'][0:-4]) # Spectrum normalized by the noise... so that we do not need to worry about the noise background so much

	# Write the configuration used for the grid in a file before starting
	print('  [e] Saving the overall configuration...')
	save_cfg(rootdir, mcmc_file, s1_file, KIC, Dnu, epsilon, variables_name, variables, cfg_file='config.cfg')

	# Define all possible combinations
	print('  [f] Defining all the requested combinations...')
	combi=get_combi(variables)
	print('			Mapping of the parameter space is the following:')
	for i in range(len(variables)):
		print('			', variables_name[i], ' :', variables[i])
	print('			Total number of combinations:', len(combi))
	print('----')
	# The main loop making the grid
	print('  [g] Computing the grid...')
	Ncombi=len(combi)
	for i in range(Ncombi):
		param_vect_list=floatvect2str(combi[i])
		print('			[', i+1, '/', Ncombi, '] ')
		print('          ', variables_name)
		print('          ', param_vect_list)
		make_subdir_lintree(rootdir, param_vect_list, verbose=False) # Creates the sequence of directory whenever necessary
		dir_out=target_subdir_lintree(rootdir, param_vect_list, list_content=False, verbose=False) # Points toward the directory of the current node in the grid 
		
		DP1=float(param_vect_list[variables_name.index('DP1')])
		Hmax_factor=float(param_vect_list[variables_name.index('Hmax_factor')])
		gamma_max_l0=float(param_vect_list[variables_name.index('gamma_max_l0')])
		rot_core=float(param_vect_list[variables_name.index('rot_core')])
		rot_envelope=float(param_vect_list[variables_name.index('rot_envelope')])
		inc=float(param_vect_list[variables_name.index('inclination')])
		delta_0l_percent=float(param_vect_list[variables_name.index('d0l_percent')])

		maxHNR_l0=hmax_from_envelope(Hnumax, Dnu, gamma_max_l0) # build the maxHNR from gamma_max_l0_var as those are closely related to the NRJ of the modes
		maxHNR_l0=Hmax_factor*maxHNR_l0
		#maxHNR_l0=100
		scoef=gamma_max_l0*2
		maxHNR_l0=Hmax_factor*max(gaussian_filter(spec, scoef, mode='mirror'))

		print('		dir_out=', dir_out)
		nu_l0, nu_p_l1, nu_g_l1, nu_m_l1, nu_l2, nu_l3, width_l0, width_m_l1, width_l2, width_l3, height_l0, height_l1, height_l2, height_l3, a1_l1, a1_l2, a1_l3=model_asymptotic(numax, Dnu, epsilon, DP1, beta_p_star=0., delta0l_percent=delta_0l_percent, gamma_max_l0=gamma_max_l0, maxHNR_l0=maxHNR_l0, alpha_star=0., 
								rot_core=rot_core, rot_envelope=rot_envelope, inclination=inc, dfmin=6, dfmax=6, output_file_rot=os.path.join(dir_out,'info.rot'))
		model_l0, model_l1, model_l2, model_l3=model_core_alt(freq, nu_l0=nu_l0, nu_l1=nu_m_l1, nu_l2=nu_l2, nu_l3=nu_l3, 
				   height_l0=height_l0, height_l1=height_l1, height_l2=height_l2, height_l3=height_l3, 
				   width_l0=width_l0, width_l1=width_m_l1, width_l2=width_l2, width_l3=width_l3, 
				   a1_l1=a1_l1, a1_l2=a1_l1, a1_l3=a1_l1, inc=inc)

		noise_background=1 # This because we removed early in this routine the noise background... ==> Only flat noise background of mean=1 now
		l=likelihood(spec, model_l0 + model_l1 + model_l2 + model_l3 + noise_background, likelihood_name='chi22') # compute the likelihood for the model over the specified range
		file_out_l=os.path.join(dir_out, 'likelihood.txt')
		fo=open(file_out_l, 'w')
		fo.write("{:16f}".format(l))
		fo.close()
		file_out_core=os.path.join(dir_out , 'spectrum')
		make_specplot_v2(freq, spec, model_l0 + noise_background, model_l1 + noise_background, model_l2 + noise_background, model_l3 + noise_background, Dnu, file_out_core, scoef=[0.1, 0.4, 0.8], fmin=-1, fmax=-1, col_spec=['lightgrey', 'darkgray', 'black'], likelihood=l)
		#exit()
	return 0

def main_grid_maker(kic='3426673'):
	rootdir='/Volumes/home/2020/ML-Siddarth/tune-rgb/grids/'
	s1_dir='/Volumes/home/2020/ML-Siddarth/s1/'
	mcmc_file_dir='/Volumes/home/2020/ML-Siddarth/setups/'

	# NEED TO DEFINE HERE THE GRID
	#variables_name=['DP1', 'gamma_max_l0', 'Hsig', 'rot_core', 'rot_envelope', 'inclination']
	variables_name=['DP1', 'Hmax_factor', 'rot_core', 'rot_envelope',  'inclination', 'd0l_percent'] # ALL WILL BE DEFAULT VALUES
	# VALUES FROM SIDDARTH TABLE: 
	# KIC      Dnu  Published_Dp  Published_core_rot(mu_hz)  published_observed_rotational_components  Predicted_Dp Predicted_q Predicted_core_rotation(mu_hz)  Predicted_inclination
	# 3426673 13.56 83.10         0.38                       2 										   82.90 		0.13 	    0.89 							43.50
	DP1=[80, 81, 82, 83, 84]
	Hmax_factor=[0.7, 0.8, 0.9]
	#gamma_max_l0=[] # We let the default value
	rot_core=[0.3, 0.4, 0.5, 0.8, 0.9, 1.0]
	rot_envelope=[0.1, 0.2]
	inclination=[30, 45, 60, 90]
	d0l_percent=[2]
	variables=[DP1, Hmax_factor, rot_core, rot_envelope, inclination, d0l_percent] 
	
	# Execute the grid_maker
	s1_file= os.path.join(s1_dir , kic) + '.sav'
	mcmc_file= os.path.join(mcmc_file_dir , kic) + '.MCMC' #format_ID(kic)
	print('[1] Starting the grid maker...')

	# Ensuring that the rootdir exists
	make_subdir_lintree(rootdir, [kic], verbose=True)
	grid_maker(os.path.join(rootdir,kic), variables, variables_name, s1_file, mcmc_file)

	print('All Done')

main_grid_maker()