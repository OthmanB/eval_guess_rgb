@read_MCMC_input.pro
function build_init_param_a1a2a3, input_data_file


; ***** CONSTANTS *****
G=6.667d-8
Teff_sun= 5777d ; same values as in the function: seismic_vsini
Dnu_sun=135.1d
numax_sun=3150d
;numax_sun=3090d
R_sun=6.96342d5 ; in km
M_sun=1.98855d30 ; in kg
rho_sun=M_sun*1d3/(4d*!pi*(R_sun*1d5)^3d/3d) ; in g.cm-3


;***** Construction du Vecteur de variables Initiales *****
file_guess=file_search(input_data_file)

if n_elements(file_guess) eq 1 then begin
	guesses=Read_MCMC_input(file_guess, Dnu=Dnu, C_l=C_l, mode_range=mode_range)
endif else begin
	print, 'multiplie initial guesses files met the KIC criteria... Emmergency stop'
	stop
endelse

Dnu0=Dnu ; use for the 'Analysis' case

epsilon=C_l/Dnu

DD=1.5d*Dnu/100
n_0=mode_range[0] & n0=n_0

print, 'Ordre n0 ',n_0

init_config=guesses.eigen_param

posl0=where(init_config[0,*] eq 0)
posl1=where(init_config[0,*] eq 1)
posl2=where(init_config[0,*] eq 2)
posl3=where(init_config[0,*] eq 3)

  print, '------------------------------------'
  print, '... Begining of the setup for the case phase_ACTION = Analysis'
  print, ''
  
  lmax=max(init_config[0,*])

  nu_input_A=reform(init_config[1,posl0])
  if posl1[0] ne -1 then nu_input_B=reform(init_config[1,posl1])

  tol_l0=init_config[2:3,posl0]
  if posl1[0] ne -1 then begin
    tol_l1=init_config[2:3,posl1]
  endif

  if lmax ge 2 then begin
    if posl2[0] ne -1 then begin
   	  nu_input_A2=reform(init_config[1,posl2])
  	 tol_l2=init_config[2:3,posl2]
    endif; else begin
      ;nu_input_A2=-1
    ;endelse
  endif
  if lmax ge 3 then begin
  	nu_input_B2=reform(init_config[1,posl3])
  	tol_l3=init_config[2:3,posl3]
  endif

  nb_ordres_l0=n_elements(nu_input_A)
  nb_ordres_l1=n_elements(nu_input_B)
  if nb_ordres_l0 eq nb_ordres_l1 then begin
  	nb_ordres=nb_ordres_l0 
  endif else begin
  	print, 'The number of modes per radial order MUST be the same!'
  	print, 'If lmax=2, Check that Nmodes(l=0) = Nmodes(l=1) = Nmodes(l=2)'
  	print, 'If lmax=3, Check that Nmodes(l=0) = Nmodes(l=1) = Nmodes(l=2) = Nmodes(l=3)'
  	print, 'Emmergency stop'
  	stop
  endelse

  Ampl_list=reform(init_config[5,posl0]) 
  Larg_param=reform(init_config[4,posl0])
  
  if n_elements(guesses.noise_param) gt 1 then begin; Case of Multiple Harvey profile
  	Bruit_param=reform(guesses.noise_param[3:*])
	Bruit_priors=guesses.noise_s2[*,3:*]
  	if 100*bruit_priors[1,3]/bruit_priors[0,3] lt 5 then begin
  		bruit_priors[1,3]=5*bruit_priors[0,3]/100 ; in case of problem with the prior uncertainty
  		bruit_priors[2,3]=5*bruit_priors[0,3]/100
  		print, 'Warning: uncertainty on prior for Height of second harvey is smaller than 5% input ==> uncertainty forced to 5%'
  	endif
  	if 100*bruit_priors[1,4]/bruit_priors[0,4] lt 0.5 then begin
  		 bruit_priors[1,4]=0.5*bruit_priors[0,4]/100 ;  tc1 limit in case of problem
  		 bruit_priors[2,4]=0.5*bruit_priors[0,4]/100
  		 print, 'Warning: uncertainty on prior for tc of second harvey is smaller than 0.5% input ==> uncertainty forced to 0.5%'
  	endif
  	if 100*bruit_priors[1,6]/bruit_priors[0,6] lt 0.05 then begin
  		bruit_priors[1,6]=0.05*bruit_priors[0,6]/100 ; white noise limit in case of problem
  		bruit_priors[2,6]=0.05*bruit_priors[0,6]/100 ; white noise limit in case of problem  
  		print, 'Warning: uncertainty on prior for white noise coefficient of second harvey is smaller than 0.05% input ==> uncertainty forced to 0.05%'
  	endif
  endif else begin ; Case of a simple white noise (For narrow band fit or for simulations)
 	 Bruit_param=guesses.noise_param
 	 Bruit_priors=dblarr(3, 1) ; one single noise parameters, with +/- sigma ==> 3 variables
 	 Bruit_priors[*,0]=guesses.noise_s2
 	 if 100*bruit_priors[1,0]/bruit_priors[0,0] lt 0.05 then begin
  		bruit_priors[1,0]=0.05*bruit_priors[0,0]/100 ; white noise limit in case of problem
  		bruit_priors[2,0]=0.05*bruit_priors[0,0]/100 ; white noise limit in case of problem  
  		print, 'Warning: uncertainty on prior for white noise coefficient of second harvey is smaller than 0.05% input ==> uncertainty forced to 0.05%'
  	endif
  endelse

  ; -------- ? DO WE WISH A FIT AROUND THE ASSYMPTOTIC RELATION ? ------
  use_alaw=0 ; if set to 0, we directly fit the frequencies. Otherwise, fit the departure from expected frequency linear pattern.
  if use_alaw eq 0 then begin
  	Dnu=0. ; input[Nmax+vis]
  	epsilon=0
   endif else begin
  	fit_p=poly_fit(findgen(n_elements(nu_input_A)),nu_input_A, 1, /double, yfit=reg) 
  	Dnu=fit_p[1]
  	epsilon=fit_p[0]/fit_p[1]
  endelse
  ; -----------------------------------------

  ; ---- Define the frequency parameters -----
  corr_nu=dblarr(lmax+1,nb_ordres)
  for i=0, nb_ordres-1 do begin
 
  	corr_nu[0,i]=nu_input_A[i]-(i+epsilon)*Dnu
  	if posl1[0] ne -1 then corr_nu[1,i]=nu_input_B[i]-(i+epsilon+0.5)*Dnu
  	if lmax ge 2 then corr_nu[2,i]=nu_input_A2[i] -(i+epsilon)*Dnu
  	if lmax ge 3 then corr_nu[3,i]=nu_input_B2[i] -(i+epsilon+0.5)*Dnu
  endfor

  ; ------ Si peak bagging pour les guess -----
  f_list=dblarr((lmax+1)*nb_ordres)
  cpt=0
  for i=0, nb_ordres-1 do begin
  	f_list[cpt]=corr_nu[0,i]
  	f_list[cpt+1]=corr_nu[1,i]
  	if lmax ge 2 then f_list[cpt+2]=corr_nu[2,i]
  	if lmax ge 3 then f_list[cpt+3]=corr_nu[3,i]
  	cpt=cpt+lmax+1
  endfor
  ;--------------------------------------------

  freq_param=[Dnu,epsilon,f_list]

  ; ********** SPLITTING *******
  ;pos_nus=where(strtrim(guesses.modes_common_names,2) eq 'val_splitting_a1' or guesses.modes_common_names eq 'val_splitting_a1')
   pos_nus=1
  if pos_nus[0] ne -1 then begin
  	splitting_param=0. ;double(guesses.modes_common[pos_nus,0])
  	nus=0. ;splitting_param
  	prior_nus=[0.,1] ;double(reform(guesses.modes_common[pos_nus, 1:*]))
  endif else begin
  	print, 'Splitting parameters not found in the configuration file!'
  	print, 'Emmergency stop!'
  	stop
  endelse

  high_order_vals=[0.,0., 0, 0, 0] ; default, no centrifugal force and no latitudinal rotation
  rho=(Dnu0/Dnu_sun)^2d * rho_sun
  a1_init=splitting_param; initial rotation rate that is used for evaluating the a2 term... It is is not a critical parameter as soon as we do a MCMC fit.
  Dnl=0.7
  eta=(4./3.)*!pi*Dnl*(a1_init*1d-6)^2/(rho*G)

  print, 'Asphericity of the star taken into account for the splitting assuming initial eta = 4/3 !pi Dnl * a1^2/(G*rho) ~ ' + strtrim(eta, 2)
  high_order_vals[0]=eta
 
  a3=0.00
  print, 'No latitudinal effect taken into account for the splitting'
  high_order_vals[1]=0.
  
  print, 'No asymetry for the modes'
  high_order_vals[4]=0.
  
  splitting_param=[splitting_param, high_order_vals]
  if lmax ge 2 then begin
  	asym_l2=[0] ;[0.,0.,0.,0.]  6 - 5 parameters for high_order_vals - 2 for asym_l1
  	splitting_param=[splitting_param,asym_l2]
  endif
  if lmax ge 3 then begin
  	asym_l3=[0.,0.,0.,0.,0.,0.]
  	splitting_param=[splitting_param,asym_l3]
  endif

  ratio_param=[1.5]
  if lmax ge 2 then ratio_param=[ratio_param,0.53]
  ;if lmax ge 3 then ratio_param=[ratio_param,0.027]
  if lmax ge 3 then ratio_param=[ratio_param,0.08]
  ;if lmax ge 3 then ratio_param=[ratio_param,0.05]
  
  inc=45.
  if inc gt 89.9999 then inc=89.9999
  prior_inc=[0, 1]
  angle=inc

  parameters_length=[n_elements(Ampl_list),n_elements(ratio_param),n_elements(freq_param),n_elements(splitting_param),$
  	n_elements(larg_param),n_elements(bruit_param),1,n_elements(bruit_orbital_param),n_elements(hyper_param)]
  input=[Ampl_list,ratio_param,freq_param,splitting_param,Larg_param,bruit_param,angle];,hyper_param];,bruit_orbital_param];hyper_param]
  
  var_names=[replicate('Height', parameters_length[0]), replicate('Visibility', parameters_length[1]), replicate('Frequency', parameters_length[2]), $
  	     'Splitting a1', 'Splitting eta', 'Splitting a3', 'Splitting (not used)', 'Splitting (not used)', 'Lorentzian asymetry', replicate('Splitting (not used)', n_elements(splitting_param)-6), $
	     replicate('Width', parameters_length[4]), replicate('Noise', parameters_length[5]), 'Inclination']
;--
  
  all_input=input

  longueur=n_elements(input) & relax=dblarr(longueur)
  relax[*]=1
  
  relax[where(var_names eq 'Inclination')]=1 ; choose to fix/free the stellar inclination

  relax[nb_ordres+n_elements(ratio_param)]=0. ; parametre Dnu (droite directrice)
  relax[nb_ordres+n_elements(ratio_param)+1]=0. ; parametre epsilon (droite directrice)

  relax[nb_ordres]=1. ; parametre ratio A1/A0 libre
  if lmax ge 2 then relax[nb_ordres+1]=1. ;parametres ratio A2/A0 libre
  if lmax ge 3 then relax[nb_ordres+2]=1. ;parametres ratio A3/A0 libre

  if n_elements(bruit_param) gt 1 then begin
  	relax[n_elements(input)-(7+n_elements(hyper_param))-1:n_elements(input)-(7+n_elements(hyper_param))-1-1+3]=0 ; le premier profile de Harvey est fix�
   endif else begin ; else, we have only a white noise...
    relax[n_elements(input)-(7+n_elements(hyper_param))-1]=1 ; Free White noise
  endelse
  ; --- Handling the Second order term of the splitting (free only if asked in the .MCMC file) ---
  high_order=[0,0,0,0,0]  ;  by defaut, centrifugal force (element 0) and latitudinal rotation (element 1) and the asymetry (element 4) are not accounted for
   pos_0=n_elements(Ampl_list)+lmax+n_elements(freq_param)+1
  relax[pos_0:pos_0+n_elements(high_order)-1]=high_order
  if lmax le 1 then begin
  	print, ' Problem! The current code requires lmax>=2'
  	print, ' The program will stop now'
  	stop
  endif
  if lmax ge 2 then relax[pos_0+n_elements(high_order):pos_0+5]=0 ; extra slots for asymetry fixed to 0
  if lmax ge 3 then relax[pos_0+n_elements(high_order):pos_0+11]=0 ; extra slots for asymetry fixed to 0

  if n_elements(splitting_param) ne lmax*(lmax+1) + 1 then begin
  	print, 'Incorrect number of parameters for the relax of splitting_param'
  	print, 'Expected number = l(l+1) +1 = ', lmax*(lmax+1) + 1
  	print, 'Nsplitting_param=', n_elements(splitting_param)
  	print, ' These numbers should be the same. Debug required. The program will stop now'
  	stop
  endif


  print,' ***** Input values AND relax list *****'
  for i=0, n_elements(input)-1 do print, var_names[i], input[i], relax[i]
  print,' *****-----------------------------*****'

  nb_param=total(parameters_length)
  nb_free=n_elements(where(relax eq 1))
  print, 'The fit includes '+strtrim(nb_param,1)+' parameters with ' + strtrim(nb_free,2) +' free parameters, including :'
  print, n_elements(Ampl_list), ' Heights'
  print, n_elements(larg_param), ' Widths'
  print, n_elements(freq_param), ' Frequencies'
  ;if lmax eq 2 then  print, n_elements(splitting_param), ' Splitting including ', strtrim(n_elements(where(relax[pos_0:pos_0+5] eq 0)),1), ' fixes'
  ;if lmax eq 3 then  print, n_elements(splitting_param), ' Splitting including ', strtrim(n_elements(where(relax[pos_0:pos_0+11] eq 0)),1), ' fixes'
  print, n_elements(bruit_param), ' Parameters of Noise'
  ;print, n_elements(bruit_orbital_param), ' Parameters of Orbital Noise (CoRoT Stars only)'
  print, n_elements(hyper_param), ' Hyper-parameters'
  print, 'Over a total of '+strtrim(nb_ordres,1)+' Radial orders '+strtrim(lmax+1,1)+' degrees (lmax='+strtrim(lmax,1)+')'

  print, '... Setup finished for the case phase_ACTION = Analysis'
  print, ''


result_init={input:dblarr(nb_param),parameters_length:dblarr(n_elements(parameters_length)), var_names:strarr(nb_param), $
	n0:long(0), fmin:0d, fmax:0d}

result_init.var_names=var_names
result_init.input=input
result_init.parameters_length=parameters_length
result_init.n0=n0

result_init.fmin=(mode_range[0]-1)*Dnu0
result_init.fmax=(mode_range[1]+1)*Dnu0

;print, 'Debug stop'
;stop
return,result_init
end

; N: number of columns
; K: maximum number of lines
function read_Ncolumns, file,N, K, skip=skip

if n_elements(skip) eq 0 then skip=1 ; defaut we skip one line only


openr, 3, file

	param=dblarr(K,N)
	a=''
	i=0d
      while EOF(3) ne 1 do begin
      	if i lt skip then readf,3,format='(q)'
        if i ge skip then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)-1
          	for j=0,N_uu-1 do begin
          		param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          	endfor
			param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		endif
		i=i+1
      endwhile

close,3
param0=param
param=param[where(param[*,1] ne 0),*]
print, 'END read'
;save, param, filename=file+'.sav'
return,param
end


