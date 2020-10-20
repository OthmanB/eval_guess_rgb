; a classical analysis
; updated on 13/08/2015 to include the splitting and inclination (or any  5 extra parameters)
function Read_MCMC_input, file_in_res, Dnu=Dnu, C_l=C_l, mode_range=mode_range

openr, 1, file_in_res


      nl=200.
	  el=4

      n=fltarr(nl)
      param=fltarr(nl,el)
      param_type=strarr(nl)

	  a='' & i=0 & out=0d
      while out lt 3 AND eof(1) eq 0 do begin ; until we didn't reach the next # symbol...
          readf,1,a ; read data
          test=byte(a)
          if test[0] eq 33 AND test[1] ne 33 then begin
          	uu=strsplit(a,' ')
          	Dnu=double(strmid(a, uu[1]))
          endif
          if test[0] eq 33 AND test[1] eq 33 then begin
          	uu=strsplit(a,' ')
          	C_l=double(strmid(a, uu[1]))
          endif
          if test[0] eq 42 then begin
          	uu=strsplit(a,' ')
          	mode_range=double([strmid(a, uu[1], uu[2]), strmid(a, uu[2])])
		  endif
          if test[0] ne 35 AND test[0] ne 33 AND test[0] ne 42 then begin
          	uu=strsplit(a,' ')
          	N_uu=N_elements(uu)-1
          	cpt=0
          	 for j=0,N_uu-1 do begin
          	 	  val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          	 	  val0=byte(val)
          	 	  val=strtrim(val0[where(val0 ne 32)],1)
          	 	  if strtrim(val,1) eq 'p ' OR strtrim(val,1) eq 'g ' OR strtrim(val,1) eq 'co' OR $
          	 	  	strtrim(val,1) eq 'p' OR strtrim(val,1) eq 'g' OR $
          	 	  	strtrim(val,1) eq 'p  ' OR strtrim(val,1) eq 'g  ' OR strtrim(val,1) eq 'co ' then begin
          	 	  	param_type(i)=strtrim(val,1)
          	 	  	cpt=cpt+1
          	 	  endif else param(i,j)=float(val)
          	endfor
		 	param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 endif
		 if test[0] eq 35 then out=out+1
			i=i+1
      endwhile
	indexes=where(param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		param=param[indexes,*]
		param_type=param_type[indexes]
	endif else param=-1

	if indexes[0] ne -1 then param=transpose(param)
	param=param[1:*,*]  ; If you don't understand why this is here... DON'T TOUCH... We simply remove an unuseful column


	i=0d & s=0 & priors=0
	 while out lt 5 AND eof(1) eq 0 do begin ; the priors, until we reach the 5th # symbol
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
				priors=[priors, double(a)]
				s=[s,1]
			endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	priors=priors[where(s ne 0)]

	i=0d & s=0 & eigen_param=fltarr(200,6)
	while out lt 6 AND eof(1) eq 0 do begin ; the priors, until we reach the 6th # symbol
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			eigen_param(i,j)=float(val)
          		endfor
		 		eigen_param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(eigen_param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		eigen_param=eigen_param[indexes,*]
	endif else eigen_param=-1
	if indexes[0] ne -1 then eigen_param=transpose(eigen_param)


	i=0d & s=0 & noise_param=-1d
	while out lt 7 AND eof(1) eq 0 do begin ; the priors, until we reach the 7th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		noise_tmp=fltarr(N_uu+1)
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(j)=float(val)
          		endfor
		 		if N_uu ne 0 then noise_tmp(N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 		if N_uu eq 0 then noise_tmp= float(a)

		 		noise_param=[noise_param, noise_tmp]
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_param[*] ne -1) ; remove the zero at the end of the table
	noise_param=noise_param[indexes]


	i=0d & s=0 & noise_tmp=fltarr(20, 3)-1
	while out lt 8 AND eof(1) eq 0 do begin ; the priors, until we reach the 8th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(i,j)=float(val)
          		endfor
		 		noise_tmp(i,N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))

		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_tmp[*,0] ne -1) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		noise_tmp=noise_tmp[indexes,*]
	endif else noise_tmp=-1
	if indexes[0] ne -1 then noise_s2=transpose(noise_tmp)

; -------- splitting and inclination follows ------
	i=0d & s=0 & modes_common=fltarr(5,3) ; up to 5 parameters
	modes_common_names=strarr(5) ; we name them here (e.g. 'splitting', 'inclination')
	while out lt 9 AND eof(1) eq 0 do begin ; the initial values for splitting/inclination + priors, until we reach the 9th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
			modes_common_names[i]=strmid(a,uu(0),uu(1)-uu(0)-1) ; take the name of the variable
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=1,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			modes_common(i,j-1)=float(val)
          		endfor
		 		modes_common(i,N_uu-1)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))

		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(modes_common[*,0] ne -1) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		modes_common=modes_common[indexes,*]
		modes_common_names=modes_common_names[indexes]
	endif else begin
		modes_common=-1
		modes_common_names=-1
	endelse
	;if indexes[0] ne -1 then noise_s2=transpose(noise_tmp)

close,1

struc={param:fltarr(n_elements(param[*,0]), n_elements(param[0,*])), param_type:strarr(n_elements(param_type)), $
	eigen_param:fltarr(n_elements(eigen_param[*,0]),n_elements(eigen_param[0,*])), $
	priors:dblarr(n_elements(priors)), noise_param:dblarr(n_elements(noise_param)), $
	noise_s2:dblarr(n_elements(noise_s2[*,0]), n_elements(noise_s2[0,*])), $
	modes_common_names:strarr(n_elements(modes_common_names)), $
	modes_common:strarr(n_elements(modes_common[*,0]), n_elements(modes_common[0,*]))}

struc.param=param
struc.param_type=param_type
struc.eigen_param=eigen_param
struc.priors=priors
struc.noise_param=noise_param
struc.noise_s2=noise_s2
struc.modes_common_names=modes_common_names
struc.modes_common=modes_common


return, struc
end

; for a full mode structure analysis (fit with > 500 parameters)
function Read_MCMCv2_input, file_in_res, Dnu=Dnu, C_l=C_l, mode_range=mode_range

openr, 1, file_in_res


      nl=400.
	  el=4

      n=fltarr(nl)
      param=fltarr(nl,el)
      param_type=strarr(nl)

	  a='' & i=0 & out=0d
      while out lt 3 AND eof(1) eq 0 do begin ; until we didn't reach the next # symbol...
          readf,1,a ; read data
          test=byte(a)
          if test[0] eq 33 AND test[1] ne 33 then begin
          	uu=strsplit(a,' ')
          	Dnu=double(strmid(a, uu[1]))
          endif
          if test[0] eq 33 AND test[1] eq 33 then begin
          	uu=strsplit(a,' ')
          	C_l=double(strmid(a, uu[1]))
          endif
          if test[0] eq 42 then begin
          	uu=strsplit(a,' ')
          	mode_range=double([strmid(a, uu[1], uu[2]), strmid(a, uu[2])])
		  endif
          if test[0] ne 35 AND test[0] ne 33 AND test[0] ne 42 then begin
          	uu=strsplit(a,' ')
          	N_uu=N_elements(uu)-1
          	cpt=0
          	 for j=0,N_uu-1 do begin
          	 	  val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          	 	  val0=byte(val)
          	 	  val=strtrim(val0[where(val0 ne 32)],1)
          	 	  if strtrim(val,1) eq 'p ' OR strtrim(val,1) eq 'g ' OR strtrim(val,1) eq 'co' OR $
          	 	  	strtrim(val,1) eq 'p' OR strtrim(val,1) eq 'g' OR $
          	 	  	strtrim(val,1) eq 'p  ' OR strtrim(val,1) eq 'g  ' OR strtrim(val,1) eq 'co ' then begin
          	 	  	param_type(i)=strtrim(val,1)
          	 	  	cpt=cpt+1
          	 	  endif else param(i,j)=float(val)
          	endfor
		 	param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 endif
		  if test[0] eq 35 then out=out+1
			i=i+1
      endwhile
	indexes=where(param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		param=param[indexes,*]
		param_type=param_type[indexes]
	endif else param=-1

	if indexes[0] ne -1 then param=transpose(param)
	param=param[1:*,*]  ; If you don't understand why this is here... DON'T TOUCH... We simply remove an unuseful column


	i=0d & s=0 & priors=0
	 while out lt 5 AND eof(1) eq 0 do begin ; the priors, until we reach the 5th # symbol
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
				priors=[priors, double(a)]
				s=[s,1]
			endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	priors=priors[where(s ne 0)]

	i=0d & s=0 & eigen_param=fltarr(200,6+6) ; last four are hyper priors on splitting / inclination + index number
	while out lt 6 AND eof(1) eq 0 do begin ; the priors, until we reach the 6th # symbol (and last one)
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			eigen_param(i,j)=float(val)
          		endfor
		 		eigen_param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(eigen_param[*,3] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		eigen_param=eigen_param[indexes,*]
	endif else eigen_param=-1
	if indexes[0] ne -1 then eigen_param=transpose(eigen_param)


	i=0d & s=0 & noise_param=-1d
	while out lt 7 AND eof(1) eq 0 do begin ; the priors, until we reach the 7th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		noise_tmp=fltarr(N_uu+1)
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(j)=float(val)
          		endfor
		 		if N_uu ne 0 then noise_tmp(N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 		if N_uu eq 0 then noise_tmp= float(a)

		 		noise_param=[noise_param, noise_tmp]
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_param[*] ne -1) ; remove the zero at the end of the table
	noise_param=noise_param[indexes]


	i=0d & s=0 & noise_tmp=fltarr(20, 3)-1
	while out lt 8 AND eof(1) eq 0 do begin ; the priors, until we reach the 8th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(i,j)=float(val)
          		endfor
		 		noise_tmp(i,N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))

		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_tmp[*,0] ne -1) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		noise_tmp=noise_tmp[indexes,*]
	endif else noise_tmp=-1
	if indexes[0] ne -1 then noise_s2=transpose(noise_tmp)

close,1

struc={param:fltarr(n_elements(param[*,0]), n_elements(param[0,*])), param_type:strarr(n_elements(param_type)), $
	eigen_param:fltarr(n_elements(eigen_param[*,0]),n_elements(eigen_param[0,*])), $
	priors:dblarr(n_elements(priors)), noise_param:dblarr(n_elements(noise_param)), $
	noise_s2:dblarr(n_elements(noise_s2[*,0]), n_elements(noise_s2[0,*]))}

struc.param=param
struc.param_type=param_type
struc.eigen_param=eigen_param
struc.priors=priors
struc.noise_param=noise_param
struc.noise_s2=noise_s2

return, struc
end


function Read_MCMCvslice_input, file_in_res, Dnu=Dnu, C_l=C_l, mode_range=mode_range

openr, 1, file_in_res


      nl=400.
	  el=4

      n=fltarr(nl)
      param=fltarr(nl,el)
      param_type=strarr(nl)

	  a='' & i=0 & out=0d
      while out lt 3 AND eof(1) eq 0 do begin ; until we didn't reach the next # symbol...
          readf,1,a ; read data
          test=byte(a)
          if test[0] eq 33 AND test[1] ne 33 then begin
          	uu=strsplit(a,' ')
          	Dnu=double(strmid(a, uu[1]))
          endif
          if test[0] eq 33 AND test[1] eq 33 then begin
          	uu=strsplit(a,' ')
          	C_l=double(strmid(a, uu[1]))
          endif
          if test[0] eq 42 then begin
          	uu=strsplit(a,' ')
          	mode_range=double([strmid(a, uu[1], uu[2]), strmid(a, uu[2])])
		  endif
          if test[0] ne 35 AND test[0] ne 33 AND test[0] ne 42 then begin
          	uu=strsplit(a,' ')
          	N_uu=N_elements(uu)-1
          	cpt=0
          	 for j=0,N_uu-1 do begin
          	 	  val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          	 	  val0=byte(val)
          	 	  val=strtrim(val0[where(val0 ne 32)],1)
          	 	  if strtrim(val,1) eq 'p ' OR strtrim(val,1) eq 'g ' OR strtrim(val,1) eq 'co' OR $
          	 	  	strtrim(val,1) eq 'p' OR strtrim(val,1) eq 'g' OR $
          	 	  	strtrim(val,1) eq 'p  ' OR strtrim(val,1) eq 'g  ' OR strtrim(val,1) eq 'co ' then begin
          	 	  	param_type(i)=strtrim(val,1)
          	 	  	cpt=cpt+1
          	 	  endif else param(i,j)=float(val)
          	endfor
		 	param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 endif
		  if test[0] eq 35 then out=out+1
			i=i+1
      endwhile
	indexes=where(param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		param=param[indexes,*]
		param_type=param_type[indexes]
	endif else param=-1

	if indexes[0] ne -1 then param=transpose(param)
	param=param[1:*,*]  ; If you don't understand why this is here... DON'T TOUCH... We simply remove an unuseful column


	i=0d & s=0 & priors=0
	 while out lt 5 AND eof(1) eq 0 do begin ; the priors, until we reach the 5th # symbol
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
				priors=[priors, double(a)]
				s=[s,1]
			endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	priors=priors[where(s ne 0)]

	i=0d & s=0 & eigen_param=fltarr(200,6+6) ; last four are hyper priors on splitting / inclination + index number
	while out lt 6 AND eof(1) eq 0 do begin ; the priors, until we reach the 6th # symbol (and last one)
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			eigen_param(i,j)=float(val)
          		endfor
		 		eigen_param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(eigen_param[*,3] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		eigen_param=eigen_param[indexes,*]
	endif else eigen_param=-1
	if indexes[0] ne -1 then eigen_param=transpose(eigen_param)


	i=0d & s=0 & noise_param=-1d
	while out lt 7 AND eof(1) eq 0 do begin ; the priors, until we reach the 7th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		noise_tmp=fltarr(N_uu+1)
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(j)=float(val)
          		endfor
		 		if N_uu ne 0 then noise_tmp(N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 		if N_uu eq 0 then noise_tmp= float(a)

		 		noise_param=[noise_param, noise_tmp]
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_param[*] ne -1) ; remove the zero at the end of the table
	noise_param=noise_param[indexes]


	i=0d & s=0 & noise_tmp=fltarr(20, 3)-1
	while out lt 8 AND eof(1) eq 0 do begin ; the priors, until we reach the 8th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(i,j)=float(val)
          		endfor
		 		noise_tmp(i,N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))

		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_tmp[*,0] ne -1) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		noise_tmp=noise_tmp[indexes,*]
	endif else noise_tmp=-1
	if indexes[0] ne -1 then noise_s2=transpose(noise_tmp)


	i=0d & s=0 & Windows_spec=fltarr(200,2) ; Windows for the fits
	while out lt 9 AND eof(1) eq 0 do begin ; the priors, until we reach the 9th # symbol (and last one)
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			Windows_spec(i,j)=float(val)
          		endfor
		 		Windows_spec(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(Windows_spec[*,0] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		Windows_spec=Windows_spec[indexes,*]
	endif else Windows_spec=-1
	if indexes[0] ne -1 then Windows_spec=transpose(Windows_spec)

close,1

struc={param:fltarr(n_elements(param[*,0]), n_elements(param[0,*])), param_type:strarr(n_elements(param_type)), $
	eigen_param:fltarr(n_elements(eigen_param[*,0]),n_elements(eigen_param[0,*])), $
	Windows_spec:fltarr(n_elements(Windows_spec[*,0]),n_elements(Windows_spec[0,*])), $
	priors:dblarr(n_elements(priors)), noise_param:dblarr(n_elements(noise_param)), $
	noise_s2:dblarr(n_elements(noise_s2[*,0]), n_elements(noise_s2[0,*]))}

struc.param=param
struc.param_type=param_type
struc.eigen_param=eigen_param
struc.Windows_spec=Windows_spec
struc.priors=priors
struc.noise_param=noise_param
struc.noise_s2=noise_s2

return, struc
end



; a classical analysis but we have a tag of which modes should be fitted
function Read_MCMC_inputv3, file_in_res, Dnu=Dnu, C_l=C_l, mode_range=mode_range

openr, 1, file_in_res


      nl=200.
	  el=4

      n=fltarr(nl)
      param=fltarr(nl,el)
      param_type=strarr(nl)

	  a='' & i=0 & out=0d
      while out lt 3 AND eof(1) eq 0 do begin ; until we didn't reach the next # symbol...
          readf,1,a ; read data
          test=byte(a)
          if test[0] eq 33 AND test[1] ne 33 then begin
          	uu=strsplit(a,' ')
          	Dnu=double(strmid(a, uu[1]))
          endif
          if test[0] eq 33 AND test[1] eq 33 then begin
          	uu=strsplit(a,' ')
          	C_l=double(strmid(a, uu[1]))
          endif
          if test[0] eq 42 then begin
          	uu=strsplit(a,' ')
          	mode_range=double([strmid(a, uu[1], uu[2]), strmid(a, uu[2])])
		  endif
          if test[0] ne 35 AND test[0] ne 33 AND test[0] ne 42 then begin
          	uu=strsplit(a,' ')
          	N_uu=N_elements(uu)-1
          	cpt=0
          	 for j=0,N_uu-1 do begin
          	 	  val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          	 	  val0=byte(val)
          	 	  val=strtrim(val0[where(val0 ne 32)],1)
          	 	  if strtrim(val,1) eq 'p ' OR strtrim(val,1) eq 'g ' OR strtrim(val,1) eq 'co' OR $
          	 	  	strtrim(val,1) eq 'p' OR strtrim(val,1) eq 'g' OR $
          	 	  	strtrim(val,1) eq 'p  ' OR strtrim(val,1) eq 'g  ' OR strtrim(val,1) eq 'co ' then begin
          	 	  	param_type(i)=strtrim(val,1)
          	 	  	cpt=cpt+1
          	 	  endif else param(i,j)=float(val)
          	endfor
		 	param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 endif
		 if test[0] eq 35 then out=out+1
			i=i+1
      endwhile
	indexes=where(param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		param=param[indexes,*]
		param_type=param_type[indexes]
	endif else param=-1

	if indexes[0] ne -1 then param=transpose(param)
	param=param[1:*,*]  ; If you don't understand why this is here... DON'T TOUCH... We simply remove an unuseful column


	i=0d & s=0 & priors=0
	 while out lt 5 AND eof(1) eq 0 do begin ; the priors, until we reach the 5th # symbol
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
				priors=[priors, double(a)]
				s=[s,1]
			endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	priors=priors[where(s ne 0)]

	i=0d & s=0 & eigen_param=fltarr(200,7) ; last element is a tag saying if 1/0 we fix the mode
	while out lt 6 AND eof(1) eq 0 do begin ; the eigenparameters, until we reach the 6th # symbol (and last one)
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			eigen_param(i,j)=float(val)
          		endfor
		 		eigen_param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(eigen_param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		eigen_param=eigen_param[indexes,*]
	endif else eigen_param=-1
	if indexes[0] ne -1 then eigen_param=transpose(eigen_param)


	i=0d & s=0 & noise_param=-1d
	while out lt 7 AND eof(1) eq 0 do begin ; the priors, until we reach the 7th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		noise_tmp=fltarr(N_uu+1)
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(j)=float(val)
          		endfor
		 		if N_uu ne 0 then noise_tmp(N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 		if N_uu eq 0 then noise_tmp= float(a)

		 		noise_param=[noise_param, noise_tmp]
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_param[*] ne -1) ; remove the zero at the end of the table
	noise_param=noise_param[indexes]


	i=0d & s=0 & noise_tmp=fltarr(20, 3)-1
	while out lt 8 AND eof(1) eq 0 do begin ; the priors, until we reach the 8th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(i,j)=float(val)
          		endfor
		 		noise_tmp(i,N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))

		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_tmp[*,0] ne -1) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		noise_tmp=noise_tmp[indexes,*]
	endif else noise_tmp=-1
	if indexes[0] ne -1 then noise_s2=transpose(noise_tmp)

close,1

struc={param:fltarr(n_elements(param[*,0]), n_elements(param[0,*])), param_type:strarr(n_elements(param_type)), $
	eigen_param:fltarr(n_elements(eigen_param[*,0]),n_elements(eigen_param[0,*])), $
	priors:dblarr(n_elements(priors)), noise_param:dblarr(n_elements(noise_param)), $
	noise_s2:dblarr(n_elements(noise_s2[*,0]), n_elements(noise_s2[0,*]))}

struc.param=param
struc.param_type=param_type
struc.eigen_param=eigen_param
struc.priors=priors
struc.noise_param=noise_param
struc.noise_s2=noise_s2



return, struc
end

; a classical analysis but we have a tag of which modes should be fitted AND we have splitting input and priors
function Read_MCMC_inputv4, file_in_res, Dnu=Dnu, C_l=C_l, mode_range=mode_range

openr, 1, file_in_res


      nl=200.
	  el=4

      n=fltarr(nl)
      param=fltarr(nl,el)
      param_type=strarr(nl)

	  a='' & i=0 & out=0d
      while out lt 3 AND eof(1) eq 0 do begin ; until we didn't reach the next # symbol...
          readf,1,a ; read data
          test=byte(a)
          if test[0] eq 33 AND test[1] ne 33 then begin
          	uu=strsplit(a,' ')
          	Dnu=double(strmid(a, uu[1]))
          endif
          if test[0] eq 33 AND test[1] eq 33 then begin
          	uu=strsplit(a,' ')
          	C_l=double(strmid(a, uu[1]))
          endif
          if test[0] eq 42 then begin
          	uu=strsplit(a,' ')
          	mode_range=double([strmid(a, uu[1], uu[2]), strmid(a, uu[2])])
		  endif
          if test[0] ne 35 AND test[0] ne 33 AND test[0] ne 42 then begin
          	uu=strsplit(a,' ')
          	N_uu=N_elements(uu)-1
          	cpt=0
          	 for j=0,N_uu-1 do begin
          	 	  val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          	 	  val0=byte(val)
          	 	  val=strtrim(val0[where(val0 ne 32)],1)
          	 	  if strtrim(val,1) eq 'p ' OR strtrim(val,1) eq 'g ' OR strtrim(val,1) eq 'co' OR $
          	 	  	strtrim(val,1) eq 'p' OR strtrim(val,1) eq 'g' OR $
          	 	  	strtrim(val,1) eq 'p  ' OR strtrim(val,1) eq 'g  ' OR strtrim(val,1) eq 'co ' then begin
          	 	  	param_type(i)=strtrim(val,1)
          	 	  	cpt=cpt+1
          	 	  endif else param(i,j)=float(val)
          	endfor
		 	param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 endif
		 if test[0] eq 35 then out=out+1
			i=i+1
      endwhile
	indexes=where(param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		param=param[indexes,*]
		param_type=param_type[indexes]
	endif else param=-1

	if indexes[0] ne -1 then param=transpose(param)
	param=param[1:*,*]  ; If you don't understand why this is here... DON'T TOUCH... We simply remove an unuseful column


	i=0d & s=0 & priors=0
	 while out lt 5 AND eof(1) eq 0 do begin ; the priors, until we reach the 5th # symbol
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
				priors=[priors, double(a)]
				s=[s,1]
			endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	priors=priors[where(s ne 0)]

	i=0d & s=0 & eigen_param=fltarr(200,10) ; last element is a tag saying if 1/0 we fix the mode
	while out lt 6 AND eof(1) eq 0 do begin ; the eigenparameters, until we reach the 6th # symbol (and last one)
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			eigen_param(i,j)=float(val)
          		endfor
		 		eigen_param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(eigen_param[*,2] ne 0) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		eigen_param=eigen_param[indexes,*]
	endif else eigen_param=-1
	if indexes[0] ne -1 then eigen_param=transpose(eigen_param)


	i=0d & s=0 & noise_param=-1d
	while out lt 7 AND eof(1) eq 0 do begin ; the priors, until we reach the 7th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		noise_tmp=fltarr(N_uu+1)
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(j)=float(val)
          		endfor
		 		if N_uu ne 0 then noise_tmp(N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))
		 		if N_uu eq 0 then noise_tmp= float(a)

		 		noise_param=[noise_param, noise_tmp]
		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_param[*] ne -1) ; remove the zero at the end of the table
	noise_param=noise_param[indexes]


	i=0d & s=0 & noise_tmp=fltarr(20, 3)-1
	while out lt 8 AND eof(1) eq 0 do begin ; the priors, until we reach the 8th # symbol OR END OF FILE
			readf,1,a ; read data
			test=byte(a)
			if test[0] ne 35 then begin
			   	uu=strsplit(a,' ')
          		N_uu=N_elements(uu)-1
          		cpt=0
          		for j=0,N_uu-1 do begin
          		 	val=strmid(a,uu(j),uu(j+1)-uu(j)-1)
          		 	val0=byte(val)
          		 	val=strtrim(val0[where(val0 ne 32)],1)
         			noise_tmp(i,j)=float(val)
          		endfor
		 		noise_tmp(i,N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)-1))

		 	endif
			if test[0] eq 35 then out=out+1
			i=i+1
	  endwhile
	indexes=where(noise_tmp[*,0] ne -1) ; remove the zero at the end of the table
	if indexes[0] ne -1 then begin
		noise_tmp=noise_tmp[indexes,*]
	endif else noise_tmp=-1
	if indexes[0] ne -1 then noise_s2=transpose(noise_tmp)

close,1

struc={param:fltarr(n_elements(param[*,0]), n_elements(param[0,*])), param_type:strarr(n_elements(param_type)), $
	eigen_param:fltarr(n_elements(eigen_param[*,0]),n_elements(eigen_param[0,*])), $
	priors:dblarr(n_elements(priors)), noise_param:dblarr(n_elements(noise_param)), $
	noise_s2:dblarr(n_elements(noise_s2[*,0]), n_elements(noise_s2[0,*]))}

struc.param=param
struc.param_type=param_type
struc.eigen_param=eigen_param
struc.priors=priors
struc.noise_param=noise_param
struc.noise_s2=noise_s2



return, struc
end
