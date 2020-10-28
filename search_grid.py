'''
	Main procedures to extract information from a grid
'''

import numpy as np
import glob
import os
import itertools
from io_tune import *
from model import *
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
#from collections import deque

def rotate_list(l, x):
	'''
		A basic function that will shift by x the list l and push the last element of l to the begining
	'''
	return l[-x:] + l[:-x]

def get_likelihood_file(directory, file, extension='txt', filter_type=None, filter_params=None):
	'''
		filter_type : None ==> No filter
		filter_type : 'index' ==> search all files that starts by file and return the index specified in filter_params
	'''
	file_out=file
	if filter_type != None: # If no filter, return the file(s) without discrimination
		if filter_type == 'index':
			file_out=glob.glob(os.path.join(directory, file + '*.' + extension)) 
			try:
				file_out=file_out[filter_params]
			except:
				print('directory: ', directory)
				print('file     : ', file)
				print('file_out : ', file_out)
				
			file_out=file_out.split('/')[-1]
	return file_out

def read_likelihood_file(dir_file, likelihood_file='likelihood.txt'):
	'''
		Read a file that contains the likelihood into a directory
		dir_file: directory into which to look for the likelihood file
		likelihood_file: name of the file to be read
	'''
	f=open(os.path.join(dir_file , likelihood_file), 'r')
	l=np.float(f.readline())
	f.close()
	return l

def get_depth_dir(dir):
	'''
		Find the depth of a directory. e.g. for '/home/hello/how/are/you'
		The depth is 5
	'''
	return np.float(len(dir.split('/'))-1)

def search_grid(dir_grid, Nbest=5, combi_file=None, minimize=False, likelihood_filename='likelihood.txt', verbose=False, file_filter_type=None, file_filter_params=None):
	'''
		A function that perform a search for the best value of the statistical
		criteria calculated over the grid
		dir_grid: root directory where all the grid content is
		Nbest: Maximum number of best solutions to be retrieved
		combi_file: [NOT YET IMPLEMENTED] If set, the program will read a 
			        file to determine the directories that have to be searched
			        If the grid is large, this is more efficient than scanning the directory to get
			        all its topology (default strategy)
		minimize: If set to True, will look for the lowest values of the statistical criteria
				  Useful if e.g. the -logLikelihood was calculated (instead of the +logLikelihood)
				  If set to False will look the highest values of the statistical criteria
		likelihood_filename: Name of the likelihood file to look for into each directories of the grid

	'''
	if combi_file == None:
		print('[1] Scanning the topology of ', dir_grid, '...')
		out, out_type=list_files(dir_grid, verbose=False)
		#dirs=out[np.where(out_type =='Directory')] # Retrieve only directories
		indices = [i for i, s in enumerate(out_type) if 'Directory' in s]
		all_dirs = [out[i] for i in indices]
		# All directories are not what we want. We want only those that have the maximum depth (e.g. 5 parameters must lead to a directory depth of 5)
		depths=np.zeros(len(all_dirs))
		i=0
		for  d in all_dirs:
			depths[i]=get_depth_dir(d)
			i=i+1
		try:
			keep=np.where(depths == np.max(depths))[0] # We keep only the longest directories
		except:
			print(max(depths))
			print(np.where(depths == np.max(depths)))
			exit()
		dirs=[]
		for i in range(len(keep)):
			dirs.append(all_dirs[int(keep[i])])
		print('		A total of ', len(dirs), 'Combinations found')
	else:
		#print('[1] Reading the combinatory file ', os.path.join(dir_grid, combi_file), '...')
		# Read a combinatory file
		#subdir_seq=read_combi()
		# Compose the names of the directories to be explored	
		#target_subdir_lintree(dir_grid, subdir_seq, verbose=False)
		#print('		A total of ', len(dirs), 'Combinations found')
		print('Error: The search-on-grid approach based on a combinatory file ')
		print('       Is not yet implemented. Please use the default approach (directory scanning)')
		print('       by setting combi_file=None in search_grid.py::search_grid()')
		print('       The program will exit now')
		exit()

	print('[2] Searching for best likelihood...')
	# ------ Initialise the tables --------
	failed_list=[] # If we failed to get a solution, we put it here
	best_sols=np.zeros(Nbest, dtype=np.float) + np.nan
	dir_best_sols=list(itertools.repeat("", Nbest))
	for i in range(Nbest):
		file=get_likelihood_file(dirs[i], likelihood_filename, extension='txt', filter_type=file_filter_type, filter_params=file_filter_params)
		#l=read_likelihood_file(d, likelihood_file=likelihood_filename)
		try:
			l=read_likelihood_file(dirs[i], likelihood_file=file)	
			best_sols[i]=l
			dir_best_sols[i]=dirs[i]
		except:
			failed_list.append(dirs[i])
	s=np.argsort(best_sols)
	tmp=dir_best_sols
	for i in range(len(s)):
		dir_best_sols[i]=tmp[s[i]]
	if minimize == False:
		best_sols=best_sols[::-1]
		dir_best_sols=dir_best_sols[::-1]
	# ----- Main loop -----
	for d in dirs[Nbest:]:
		file=get_likelihood_file(d, likelihood_filename, extension='txt', filter_type=file_filter_type, filter_params=file_filter_params)
		#l=read_likelihood_file(d, likelihood_file=likelihood_filename)
		try:
			l=read_likelihood_file(d, likelihood_file=file)
			if minimize == False:
				if best_sols[0] <= l:
					best_sols=shift(best_sols,1)
					#dir_best_sols=dir_best_sols.rotate(1)		
					dir_best_sols=rotate_list(dir_best_sols,1)
					best_sols[0]=l
					dir_best_sols[0]=d
					if verbose == True:
						print('likelihood=', l)
						print('dir=', d)
						print('-----')
			else:
				if best_sols[0] >=l:
					best_sols=shift(best_sols,1)
					#dir_best_sols=dir_best_sols.rotate(1)
					dir_best_sols=rotate_list(dir_best_sols,1)
					best_sols[0]=l
					dir_best_sols[0]=d
					if verbose == True:
						print('likelihood=', l)
						print('dir=', d)
						print('-----')
		except:
			print('Failed to find likelihood files in :', d)
			print('				Pursuing...')
			failed_list.append(d)
	return best_sols, dir_best_sols

# KIC      Dnu  Published_Dp  Published_core_rot(mu_hz)  published_observed_rotational_components  Predicted_Dp Predicted_q Predicted_core_rotation(mu_hz)  Predicted_inclination
# 3426673  13.56 83.10         0.38                       2 										   82.90 		0.13 	    0.89 							43.50
	
#dir_grid='/Volumes/home/2020/ML-Siddarth/eval_guess_rgb/grids/3426673/'
dir_grid='/Volumes/home/2020/ML-Siddarth/eval_guess_rgb/grids/fine/3426673/'
#dir_grid='/Volumes/home/2020/ML-Siddarth/eval_guess_rgb/grids/3426673/80.000000/0.157259/0.700000/' # For test only
#best_sols, dir_best_sols=search_grid(dir_grid, Nbest=5, combi_file=None, minimize=False, likelihood_filename='likelihood.txt')

best_sols, dir_best_sols=search_grid(dir_grid, Nbest=10, combi_file=None, minimize=False, likelihood_filename='likelihood_', file_filter_type="index", file_filter_params=3)
print('[Slice 3] Best solutions found at those locations:')
print("{:10} {:10}".format('likelihood', 'directory'))
for i in range(len(best_sols)):
	if dir_best_sols[i] != '':
		print("{:10}     {}".format(best_sols[i], dir_best_sols[i]))

best_sols, dir_best_sols=search_grid(dir_grid, Nbest=10, combi_file=None, minimize=False, likelihood_filename='likelihood_', file_filter_type="index", file_filter_params=4)
print('[Slice 4] Best solutions found at those locations:')
print("{:10} {:10}".format('likelihood', 'directory'))
for i in range(len(best_sols)):
	if dir_best_sols[i] != '':
		print("{:10}     {}".format(best_sols[i], dir_best_sols[i]))

best_sols, dir_best_sols=search_grid(dir_grid, Nbest=10, combi_file=None, minimize=False, likelihood_filename='likelihood_', file_filter_type="index", file_filter_params=5)
print('[Slice 5] Best solutions found at those locations:')
print("{:10} {:10}".format('likelihood', 'directory'))
for i in range(len(best_sols)):
	if dir_best_sols[i] != '':
		print("{:10}     {}".format(best_sols[i], dir_best_sols[i]))


