;*******************************************************************************
;This function is used by the nr_powell minimization routine
;FIT GLOBAL du l=0,1,2,3
;*******************************************************************************
;param : vecteur des parametres
;length_param : nombre de variables par parametres
;x , y : donn�es abscisse et ordonn�es
@function_rot.pro
function global_fit_a1a2a3,param,length_param,x

ndata=n_elements(x)
param2=param


;vecteur length_param :
;[Nmax,1,1,j,k,l,m,1]
; dans l'ordre : n amplitudes ; Ratio A1/A0 ; Ratio A2/A0 ; j param freq ; k param split; l param largeur; m bruits; 1 angle
;[ 0...Nmax-1 ; Nmax ... Nmax+1 ; Nmax+2 ... Nmax+j+1 ; Nmax+j+2 ... Nmax+j+k+1 ; Nmax+j+k+2 ... Nmax+j+k+l+1; Nmax+j+k+l+2 ... Nmax+j+k+m+1 ; Nmax+j+k+m+2]
Nmax=length_param[0] ; nombre d'amplitudes l=0 (mais aussi le nombre de grande separation � ajuster)
vis=length_param[1] ; nombre de visibilit�s consider�s (1 si lmax=1, 2 si lmax=2, etc...)
j=length_param[2] ; 3 x nombres de degrees pour la frequence
k=length_param[3] ; nombre de degrees pour le splitting
l=length_param[4] ; nombre de degrees pour la largeur
m=length_param[5] ; nombre de degrees pour le bruit
orb=length_param[7] ; ; nombre de degrees pour le bruit orbital ou d'hyper_param

Nmax=length_param[0]

Dnu=param2[Nmax+vis]
angle=param2[Nmax+j+k+l+m+vis] ; position n+j+k+m+l+2
param2[Nmax]=abs(param2[Nmax]) ; on traite avec abs du ratio ... le prior se charge de penaliser les modeles ou c'est pas le cas
if vis ge 2 then param2[Nmax+1]=abs(param2[Nmax+1]) ; on traite avec abs du ratio
if vis ge 3 then param2[Nmax+2]=abs(param2[Nmax+2]) ; on traite avec abs du ratio

ratios_l1=amplitude_ratio(1,angle)
ratios_l2=amplitude_ratio(2,angle)
ratios_l3=amplitude_ratio(3,angle)

;**** boucle de creation du spectre total ! ****
model_final=0.

cpt=0.
f_l0_base=dblarr(Nmax) & f_l1_base=f_l0_base & f_l2_base=f_l0_base & f_l3_base=f_l0_base
for n=0., Nmax-1 do begin
	f_l0_base[n]= (n +param2[Nmax+vis+1])*Dnu 
	f_l1_base[n]= (n +param2[Nmax+vis+1]+0.5)*Dnu
	f_l0_base[n]= f_l0_base[n] + param2[Nmax+vis+2+cpt]
	f_l1_base[n]= f_l1_base[n] + param2[Nmax+vis+3+cpt]

	if vis ge 2 then begin
	    f_l2_base[n] = (n +param2[Nmax+vis+1]+ 1.)*Dnu-Dnu
		f_l2_base[n] = f_l2_base[n] + param2[Nmax+vis+4+cpt]; les corrections � la loi asymptotique
	endif
	if vis ge 3 then begin
		f_l3_base[n]= (n +param2[Nmax+vis+1]+3./2)*Dnu-Dnu
		f_l3_base[n]=f_l3_base[n] + param2[Nmax+vis+5+cpt]; les corrections � la loi asymptotique
	endif
cpt=cpt+vis+1
endfor

for n=0., Nmax-1 do begin
;**** conversion log --> base 10 ****
; ****** interpretation des Amplitudes ****
Amp_l0=abs(param2[n])
Amp_l1=param2[Nmax]*Amp_l0 ; n c'est l'indice de Amp_l0, Nmax l'indice du ratio_l1/l0
if vis ge 2 then Amp_l2=param2[Nmax+1]*Amp_l0 ; n c'est l'indice de Amp_l0, Nmax+1 l'indice du ratio_l2/l0
if vis ge 3 then Amp_l3=param2[Nmax+2]*Amp_l0 ; n c'est l'indice de Amp_l0, Nmax+2 l'indice du ratio_l3/l0
; ****** Creation des polynomes ******
; **** les frequences ****
f_l0=0 & f_l1=0 & f_l2=0 & f_l3=0

f_l0=f_l0_base[n]
f_l1=f_l1_base[n]

if vis ge 2 then begin
	f_l2=f_l2_base[n]
endif
if vis ge 3 then begin
	f_l3=f_l3_base[n]
endif

cpt=cpt+vis+1

; **** Largeurs ****
ggMax=l
gamma_l0=0. & gamma_l1=0. & gamma_l2=0. & gamma_l3=0.

gamma_l0=abs(param2[Nmax+j+k+vis+n]) ; une largeur par Dnu et conversion log --> normal !!
gamma_l1=abs(interpol(param2[Nmax+j+k+vis:Nmax+j+k+vis+l-1], f_l0_base, f_l1)) ; it is in principle a better approximation
gamma_l2=abs(interpol(param2[Nmax+j+k+vis:Nmax+j+k+vis+l-1], f_l0_base, f_l2)) ; it is in principle a better approximation
if vis ge 3 then gamma_l3=abs(interpol(param2[Nmax+j+k+vis:Nmax+j+k+vis+l-1], f_l0_base, f_l3)) ; it is in principle a better approximation

; **** splitting ****
ssMax=k
splitting_l1=0. & splitting_l2=0. & splitting_l3=0.
ss=0
splitting_l1=param2[Nmax+j+vis+ss]
splitting_l2=param2[Nmax+j+vis+ss]
splitting_l3=param2[Nmax+j+vis+ss]

asym_l0=param2[Nmax+j+vis+ss+1:Nmax+j+vis+5] ; we pass it here because l=0 uses asym[4]
asym_l1=asym_l0
if vis ge 2 then asym_l2=asym_l0 ; we assume the same asymmetry parameters for all l
if vis ge 3 then asym_l3=asym_l0

splitting_l1=abs(splitting_l1) & splitting_l2=abs(splitting_l2) & splitting_l3=abs(splitting_l3)

;--------- Creation des modeles de fit lorentzien et en frequence ---------------
;------ l=0,1,2 ------
model_l0=optimised_lorentzian_calculation(x,Amp_l0,f_l0,0.,asym_l0,gamma_l0,0,1,ndata)
model_l1=optimised_lorentzian_calculation(x,Amp_l1,f_l1,splitting_l1,asym_l1,gamma_l1,1,ratios_l1,ndata)
if vis ge 2 then model_l2=optimised_lorentzian_calculation(x,Amp_l2,f_l2,splitting_l2,asym_l2,gamma_l2,2,ratios_l2,ndata) else $
	modele_l2=0.
if vis ge 3 then begin
	model_l3=optimised_lorentzian_calculation(x,Amp_l3,f_l3,splitting_l3,asym_l3,gamma_l3,3,ratios_l3,ndata)
	endif else $
		model_l3=0
model_final=model_l0+model_l1+model_l2+model_l3+model_final
endfor

if m eq 1 then begin
	bruit=abs(param2[Nmax+j+k+l+vis])
endif else begin
	; **** vecteur bruit ****
	param2[Nmax+j+k+l+vis]=abs(param2[Nmax+j+k+l+vis]) ; valeurs positives de la cte de tps (le prior se charge des val neg)
	param2[Nmax+j+k+l+vis+1]=abs(param2[Nmax+j+k+l+vis+1]) ; valeurs positives de la cte de tps (le prior se charge des val neg)
	param2[Nmax+j+k+l+vis+2]=abs(param2[Nmax+j+k+l+vis+2])
	param2[Nmax+j+k+l+vis+3]=abs(param2[Nmax+j+k+l+vis+3])
	
	; ---- si model de Harvey ---
	bruit=param2[Nmax+j+k+l+vis]/(1.+(param2[Nmax+j+k+l+vis+1]*(1e-3*x))^param2[Nmax+j+k+l+vis+2]) ;modele veritable de Harvey
	if m gt 4 then begin
		param2[Nmax+j+k+l+vis+4]=abs(param2[Nmax+j+k+l+vis+4])
		param2[Nmax+j+k+l+vis+5]=abs(param2[Nmax+j+k+l+vis+5])
		param2[Nmax+j+k+l+vis+6]=abs(param2[Nmax+j+k+l+vis+6])
		bruit=bruit+param2[Nmax+j+k+l+vis+3]/(1.+(param2[Nmax+j+k+l+vis+4]*(1e-3*x))^param2[Nmax+j+k+l+vis+5]) ;modele veritable de Harvey
		bruit=bruit+param2[Nmax+j+k+l+vis+6] ; Ajout d'un bruit blanc
	endif else $
		bruit=bruit+param2[Nmax+j+k+l+vis+3] ; Ajout d'un bruit blanc
endelse
; ---------------------------

model_final=model_final+bruit

return,model_final
end


; Fonctions qui calcul de maniere optimis� l'amplitude d'un mode l
; donn�e. Elle tient compte des multiplets (splitting d'ordre 1)
function optimised_lorentzian_calculation, x,Amp,f_c,f_s,asym,gamma_mode,l,V,ndata

c=20. ; critere choisi pour definir la zone de calcul

if gamma_mode ge 1. AND f_s ge 1. then begin
if l ne 0 then begin
	pmin=f_c-c*(l*f_s+gamma_mode)
	pmax=f_c+c*(l*f_s+gamma_mode)
	endif else begin
	pmin=f_c-c*gamma_mode*2.2
	pmax=f_c+c*gamma_mode*2.2
	endelse
endif

if gamma_mode le 1. AND f_s ge 1. then begin
	if l ne 0 then begin
	pmin=f_c-c*(l*f_s+1.)
	pmax=f_c+c*(l*f_s+1.)
	endif else begin
	pmin=f_c-c*2.2
	pmax=f_c+c*2.2
	endelse
endif

if gamma_mode ge 1. AND f_s le 1. then begin
	if l ne 0 then begin
	pmin=f_c-c*(l+gamma_mode)
	pmax=f_c+c*(l+gamma_mode)
	endif else begin
	pmin=f_c-c*2.2*gamma_mode
	pmax=f_c+c*2.2*gamma_mode
	endelse
endif

if gamma_mode le 1. AND f_s le 1. then begin
	if l ne 0 then begin
	pmin=f_c-c*(l+1.)
	pmax=f_c+c*(l+1.)
	endif else begin
	pmin=f_c-c*2.2
	pmax=f_c+c*2.2
	endelse
endif

;*** gestion des bornes xmin et xmax ****
pas=x[1]-x[0]
xmin=x[0] & xmax=x[ndata-1]
if (pmax-pas) lt xmin then pmax=xmin+c
if (pmin+pas) ge xmax then pmin=xmax-c
;****************************************
d_p=where(x ge pmin AND x le pmax)

if d_p[0] eq -1 then begin
	print, 'pmin='+strtrim(pmin,1)+'  pmax='+strtrim(pmax,1)
	print, 'f_c='+strtrim(f_c,1)+'  l='+strtrim(l,1)+'  f_s='+strtrim(f_s,1)+'  gamma_mode='+strtrim(gamma_mode,1)
	print, 'probleme avec d_p !!! Ensemble vide !!!'
endif

x_p=x[d_p]

model_l=build_l_mode(x_p,Amp,f_c,f_s,asym,gamma_mode,l,V)
;wait,0.5
;stop
result=dblarr(ndata)
result[d_p]=model_l

return, result
end

; Calculate the model for a given degree
; Elle tient compte des multiplets (splitting d'ordre 1)
; Elle tient compte de la force centrifuge (ordre 2) par le biais de a2=asym[0]= eta * Dnl ~ eta * 0.7 [no unit]
; Elle tient compte de l'effet de rotation latitudinal (ordre 3) par le biais de a3=asym[1] [microHz]
; A noter que a2 ~ a1^2/(G rho_sun) * Dnu_sun/Dnu [no unit]. Typiquement, rho=[0.1, 1.5] g/cm^3 pour une solar-like. G=6.67d-8 cm^3/g/s^2
function build_l_mode, x_l,Amp_l,f_c_l,f_s,asym,gamma_l,l,V

result=0

;window,1
;ymax=max(Amp_l*V)
;plot,x_l,x_l,/nodata,background=fsc_color('white'),yr=[0,ymax],color=fsc_color('black'), $
	xr=[f_c_l - 2*l*f_s - gamma_l, f_c_l + 2*l*f_s + gamma_l]
;i=0
for m=-l,l do begin
	if l ne 0 then begin
		Qlm=1.*(l*(l+1.) - 3.*m^2)/((2.*l-1.)*(2.*l+3.)) ;  a2=asym[0] coeficient for centrifugal force
		if l eq 1 then clm=m ;  a3=asym[1] latitudinal coeficient for l=1
		if l eq 2 then clm=1.*(5*m^3 - 17*m)/3. ;  a3=asym[1] latitudinal coeficient for l=2
		if l eq 3 then clm=0 ; a3=asym[1] IS NOT IMPLEMENTED FOR l=3

		profile=2*( x_l- f_c_l*( 1d + asym[0]*Qlm) + m*f_s + clm*asym[1]) / gamma_l 
	endif else $
		profile=2*( x_l-f_c_l)/gamma_l

  profile=profile^2
  asymetry=(1d + asym[4]*(x_l/f_c_l - 1d))^2 + (0.5*gamma_l*asym[4]/f_c_l)^2 ; term of asymetry for each lorentzian

  ;if abs(m) eq 0 then col=fsc_color('blue')
  ;if abs(m) eq 1 then col=fsc_color('red')
  ;if abs(m) eq 2 then col=fsc_color('Orange')
  ;oplot, x_l,Amp_l*V[m+l]/(1+profile),color=col;,linestyle=2
  ;plots, [f_c_l+m*f_s, f_c_l+m*f_s], [0,ymax], color=fsc_color('black'), linestyle=2
  ;stop
  result=result + asymetry*Amp_l*V[m+l]/(1+profile)

;i=i+2
endfor

;stop
return, result
end

