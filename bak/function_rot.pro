Function amplitude_ratio,l,beta
; beta is in degree
        angle=!pi*beta/180.
        zz=function_rot(l,angle)
	zz=zz(*,l)^2
        norm=total(zz)
	return,zz/norm
end



Function function_rot,l,beta
        dim=2*l+1
        mat=dblarr(dim,dim)

	for i=0,l do begin
 	    for j=-i,i do begin
 	          mat(i+l,j+l)=dmm(l,i,j,beta)
 	    endfor
	endfor

	for i=-l,0 do begin
	    for j=i,-i do begin
 	          mat(i+l,j+l)=mat(-i+l,-j+l)*(-1.d0)^(i-j)
  	    endfor
	endfor

	for j=0,l do begin
	    for i=-j,j do begin
 	          mat(i+l,j+l)=dmm(l,j,i,-beta)
	    endfor
	endfor

	for j=-l,0 do begin
 	    for i=j,-j do begin
	          mat(i+l,j+l)=mat(-i+l,-j+l)*(-1.d0)^(i-j)
	    endfor
	endfor

	return,mat
end


;************************
; compute dm1,m2(beta)  for  m1+m2>=0  and m1-m2>=0
;**************************
Function dmm,l,m1,m2,beta

        co=cos(beta/2.d0)
        si=sin(beta/2.d0)

        sum=0.d0

        for s=0,l-m1 do begin
             var=0.d0
             var=combi(l+m2,l-m1-s)*combi(l-m2,s)*(-1.d0)^(l-m1-s)
             var=var*co^(2*s+m1+m2)*si^(2*l-2*s-m1-m2)
             sum=sum+var
        endfor

        sum=sum*sqrt(factorial(l+m1)*factorial(l-m1)*1.d0)
        sum=sum/sqrt(factorial(l+m2)*factorial(l-m2)*1.d0)

        return,sum
end

;************************************************
; Compute combinations n r  ,    n>r
;************************************************
Function combi,n,r
        c=factorial(n)/factorial(n-r)/factorial(r)
        return,c
end
