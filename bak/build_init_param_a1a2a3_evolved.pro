@read_MCMC_input.pro
function build_init_param_a1a2a3_evolved, input_data_file

;print, 'NOT IMPLEMENTED'
;stop

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
  if posl1[0] ne -1 then begin
    nu_input_B=reform(init_config[1,posl1])
    tol_l1=init_config[2:3,posl1]
  endif
  tol_l0=init_config[2:3,posl0]

  if lmax ge 2 then begin
  	if posl2[0] ne -1 then begin
      nu_input_A2=reform(init_config[1,posl2])
  	  tol_l2=init_config[2:3,posl2]
    endif
  endif
  if lmax ge 3 then begin
    if posl3[0] ne -1 then begin
    	nu_input_B2=reform(init_config[1,posl3])
    	tol_l3=init_config[2:3,posl3]
    endif
  endif

  nb_ordres_l0=n_elements(nu_input_A)
  nb_ordres_l1=n_elements(nu_input_B)
  nb_ordres_l2=n_elements(nu_input_A2)
  nb_ordres_l3=n_elements(nu_input_B2)

  ;nb_ordres=nb_ordres_l0 
  nb_ordres=max([nb_ordres_l0, nb_ordres_l1, nb_ordres_l2, nb_ordres_l3])

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

  f_list=-1

  if posl0[0] ne -1 then $
    if f_list[0] eq -1 then begin
      f_list=nu_input_A
    endif else begin
      f_list=[f_list, nu_input_A]
    endelse
  if posl1[0] ne -1 then $
    if f_list[0] eq -1 then begin
      f_list=nu_input_B
    endif else begin
      f_list=[f_list, nu_input_B]
    endelse
  if posl2[0] ne -1 then $
    if f_list[0] eq -1 then begin
      f_list=nu_input_A2
    endif else begin
      f_list=[f_list, nu_input_A2]
    endelse
  if posl3[0] ne -1 then $
    if f_list[0] eq -1 then begin
      f_list=nu_input_B2
    endif else begin
      f_list=[f_list, nu_input_B2]
    endelse

  ;freq_param=[0,0,f_list]
  freq_param=f_list

  ; ********** SPLITTING *******
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
  
  splitting_param=[a1_init, high_order_vals]
  ;if lmax ge 2 then begin
  ;	asym_l2=[0] ;[0.,0.,0.,0.]  6 - 5 parameters for high_order_vals - 2 for asym_l1
  ; 	splitting_param=[splitting_param,asym_l2]
  ; endif
  ;if lmax ge 3 then begin
  ;	asym_l3=[0];[0.,0.,0.,0.,0.,0.]
  ;	splitting_param=[splitting_param,asym_l3]
  ;endif

  ;ratio_param=[1.5]
  ;if lmax ge 2 then ratio_param=[ratio_param,0.53]
  ;if lmax ge 3 then ratio_param=[ratio_param,0.08]
  ratio_param=[1.5,0.53,0.08]

  inc=45.
  if inc gt 89.9999 then inc=89.9999
  prior_inc=[0, 90]
  angle=inc

  parameters_length=[n_elements(Ampl_list),n_elements(ratio_param),nb_ordres_l0, nb_ordres_l1, nb_ordres_l2, nb_ordres_l3,n_elements(splitting_param),$
  	n_elements(larg_param),n_elements(bruit_param),n_elements(inc)]
  input=[Ampl_list,ratio_param,freq_param,splitting_param,Larg_param,bruit_param,angle]
  
  var_names=[replicate('Height', parameters_length[0]), replicate('Visibility', parameters_length[1])]
  if posl0[0] ne -1 then var_names=[var_names, replicate('Frequency l=0', parameters_length[2])]
  if posl1[0] ne -1 then var_names=[var_names, replicate('Frequency l=1', parameters_length[3])]
  if posl2[0] ne -1 then var_names=[var_names, replicate('Frequency l=2', parameters_length[4])]
  if posl3[0] ne -1 then var_names=[var_names, replicate('Frequency l=3', parameters_length[5])]
  var_names=[var_names, 'Splitting a1', 'Splitting eta', 'Splitting a3', 'Splitting (not used)', 'Splitting (not used)', 'Lorentzian asymetry', replicate('Width', parameters_length[7]), $
        replicate('Noise', parameters_length[8]), 'Inclination']

  all_input=input



nb_param=total(parameters_length)
result_init={input:dblarr(nb_param),parameters_length:dblarr(n_elements(parameters_length)), var_names:strarr(nb_param), $
	n0:long(0), fmin:0d, fmax:0d}

result_init.var_names=var_names
result_init.input=input
result_init.parameters_length=parameters_length
result_init.n0=n0

result_init.fmin=(mode_range[0]-1)*Dnu0
result_init.fmax=(mode_range[1]+1)*Dnu0

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


