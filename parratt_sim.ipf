#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//=======================================================
function parratt_prepoptconst(optconstpath, e1, e2, de)
wave/t optconstpath //N rows : each row contains path to opt.const, 3 columns {engy, delta, beta}
variable e1,e2,de

make/o/c/n=((e2-e1)/de+1, dimsize(optconstpath,0)) refracindex
make/o/n=((e2-e1)/de+1) tmpb, tmpd
setscale/P x, e1, de, refracindex, tmpb, tmpd

variable m
for (m=0;m<dimsize(optconstpath,0);m+=1)
	wave engy1=$optconstpath[m][0]
	wave delta1=$optconstpath[m][1]
	wave beta1=$optconstpath[m][2]
	interpolate2/I=3/T=1/Y=tmpd engy1, delta1
	interpolate2/I=3/T=1/Y=tmpb engy1, beta1
	refracindex[][m] = cmplx(1-tmpd[p], tmpb[p])
endfor

killwaves tmpb, tmpd
end
//=======================================================
function parratt_simEQM(eqmap, structwave, refractindex, mode)
wave eqmap //dimensions must be specified: x-energy in eV; y-q_z in 1/Angstrom
wave structwave //3 columns: {thickness (A), roughnes (A), index of refractindex} // start with substrate=0 and up to top layer (do not include vacuum)
wave/c refractindex //N columns : each contain energy dependence of comoplex refr index of a compound
variable mode //
// 0: Parratt, 
// 1: no denominator
// 18: for testing - same as 1, but kzwave is not updated (used from previous calculation)
// 2: no denominator AND kz_i = kz_vac = qz/2
// 3: no denominator AND kz_i = qz/2 AND cos(alpha_i) = cos(alpha_vac)
// 4: same as 1, but use average kz
// 5: use average kz AND fresnel ~ dn*k0^2/kz_avg^2 ---- good approx, FT -like
// 6: reflection from substrate (reference layer) only
// 7: same as 5 AND kzwave averaged over energy
eqmap=0

make/o/n=(dimsize(eqmap,0), dimsize(eqmap,1)) cos_alpha_vac
make/o/c/n=(dimsize(eqmap,0), dimsize(eqmap,1), dimsize(structwave,0)) cosalpha, phshift=0, roughdamping, Fresnelcoef=0, Reflection=0
if (mode!=18)
	make/o/c/n=(dimsize(eqmap,0), dimsize(eqmap,1), dimsize(structwave,0)) kzwave
endif

make/o/c/n=(dimsize(eqmap,0), dimsize(structwave,0)) optical_const //per layer

setscale/P x, dimoffset(eqmap,0), dimdelta(eqmap,0), Fresnelcoef, Reflection, cosalpha, phshift, kzwave, cos_alpha_vac, optical_const, roughdamping
setscale/P y, dimoffset(eqmap,1), dimdelta(eqmap,1), Fresnelcoef, Reflection, cosalpha, phshift, kzwave, cos_alpha_vac, roughdamping

optical_const = refractindex(x)[structwave[q][2]] // linearly interpolate optical constants per layer

variable hbar_c = 1973.26972
cos_alpha_vac = 0.5*hbar_c*y/x
duplicate/o cos_alpha_vac, sin_alpha_vac
sin_alpha_vac = (1-cos_alpha_vac^2)^.5

variable m,k, numenergy
make/c/o/n=(dimsize(structwave,0)) temp
//layer-specific quantities
if (mode<2)
	cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
	kzwave[][][] = optical_const[p][r]*cosalpha[p][q][r]/cos_alpha_vac[p][q]*y/2 
elseif (mode==2)
	cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
	kzwave[][][] = y/2 
elseif (mode==3)
	cosalpha[][][] = cos_alpha_vac[p][q]
	kzwave[][][] = y/2 
elseif (mode==4 || mode==5 || mode==6)
	cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
	kzwave[][][] = optical_const[p][r]*cosalpha[p][q][r]/cos_alpha_vac[p][q]*y/2 
	for (k=0; k<dimsize(cosalpha,0); k+=1)
		for (m=0; m<dimsize(cosalpha,1); m+=1)
			temp = kzwave[k][m][p]
			kzwave[k][m][] = mean(temp)
		endfor
	endfor	
elseif (mode==7)
	make/o/c/n=(dimsize(kzwave,0), dimsize(kzwave,2)) temp2
	cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
	kzwave[][][] = optical_const[p][r]*cosalpha[p][q][r]/cos_alpha_vac[p][q]*y/2 
	numenergy = dimsize(kzwave,0)
	for (m=0; m<dimsize(kzwave,1); m+=1)
		temp2 = kzwave[p][m][q]
		kzwave[][m][] = mean(temp2)
	endfor	
endif
phshift[][][1,dimsize(structwave,0)-1] = exp(2*structwave[r][0]*kzwave[p][q][r]*cmplx(0,1)) // exp(2*kz*d*i)

//interface-specific quantities (m'th value = interface between layers m and m+1
roughdamping[][][0,dimsize(structwave,0)-2] = exp(-2*kzwave[p][q][r]*kzwave[p][q][r+1]*structwave[r][1]^2)
roughdamping[][][dimsize(structwave,0)-1] = exp(-2*kzwave[p][q][r]*y/2*structwave[r][1]^2) //top layer and vacuum k_z_vac==q_z/2=y/2
//roughdamping[][][0,dimsize(structwave,0)-2] = 1
//Fresnel coef. : calculation for sigma light
if (mode==5)
	Fresnelcoef[][][0,dimsize(structwave,0)-2] = 0.5*(x/hbar_c)^2/kzwave[p][q][r]^2*roughdamping[p][q][r]*(-optical_const[p][r] + optical_const[p][r+1])
	Fresnelcoef[][][dimsize(structwave,0)-1] = 0.5*(x/hbar_c)^2/kzwave[p][q][r]^2*roughdamping[p][q][r]*(1-optical_const[p][r]) //top layer and vacuum
//	Fresnelcoef[][][dimsize(structwave,0)-1] = roughdamping[p][q][r]*(cos_alpha_vac[p][q]-optical_const[p][r]*cosalpha[p][q][r])/(cos_alpha_vac[p][q] + optical_const[p][r]*cosalpha[p][q][r]) //top layer and vacuum
elseif (mode==6)
//	Fresnelcoef[][][0] = 0.5*(x/hbar_c)^2/kzwave[p][q][r]^2*roughdamping[p][q][r]*(1 - optical_const[p][r])
//	Fresnelcoef[][][1,dimsize(structwave,0)-1] = 0
//substrate not good - high freq oscillations in q-space (z_sub is high)
	Fresnelcoef[][][0,dimsize(structwave,0)-2] = 0
	Fresnelcoef[][][dimsize(structwave,0)-1] = 0.5*(x/hbar_c)^2/kzwave[p][q][r]^2*roughdamping[p][q][r]*(1-optical_const[p][r]) //top layer and vacuum
//elseif (mode==6)
//	Fresnelcoef[][][0,dimsize(structwave,0)-2] = 0
//	Fresnelcoef[][][dimsize(structwave,0)-1] = roughdamping[p][q][r]*(cos_alpha_vac[p][q]-optical_const[p][r]*cosalpha[p][q][r])/(cos_alpha_vac[p][q] + optical_const[p][r]*cosalpha[p][q][r]) //top layer and vacuum
else
	Fresnelcoef[][][0,dimsize(structwave,0)-2] = roughdamping[p][q][r]*(-optical_const[p][r]*cosalpha[p][q][r] + optical_const[p][r+1]*cosalpha[p][q][r+1])/(optical_const[p][r]*cosalpha[p][q][r] + optical_const[p][r+1]*cosalpha[p][q][r+1])
	Fresnelcoef[][][dimsize(structwave,0)-1] = roughdamping[p][q][r]*(cos_alpha_vac[p][q]-optical_const[p][r]*cosalpha[p][q][r])/(cos_alpha_vac[p][q] + optical_const[p][r]*cosalpha[p][q][r]) //top layer and vacuum
endif

// reflection between substrate and first layer
Reflection[][][0] = Fresnelcoef[p][q][0]
for (m=1; m<dimsize(structwave,0); m=m+1) //start with dimsize(structwave,0)-2 since no reflection from bottom of substrate
	if (mode==0)
		Reflection[][][m] = (Fresnelcoef[p][q][r] + Reflection[p][q][r-1]*phshift[p][q][r])/(1+Fresnelcoef[p][q][r]*Reflection[p][q][r-1]*phshift[p][q][r])
	else
		Reflection[][][m] = (Fresnelcoef[p][q][r] + Reflection[p][q][r-1]*phshift[p][q][r])
	endif
endfor

eqmap =log(cabs(Reflection[p][q][dimsize(structwave,0)-1]))*2
duplicate/o eqmap $nameofwave(eqmap)+"_phase"
wave eqmaphase=$nameofwave(eqmap)+"_phase"
eqmaphase = 0
eqmaphase = imag(r2polar(Reflection[p][q][dimsize(structwave,0)-1]))

end
//=======================================================
function parratt_simEQW(eqmap, ewave, qwave, structwave, optical_const, mode) //same as parratt_simEQW escept E and q_z defined separately
//function parratt_simEQW(eqmap, ewave, qwave, structwave, refractindex, mode) //same as parratt_simEQW escept E and q_z defined separately
wave ewave// 1D wave with energy values for the rows ofthe output, eV
wave eqmap, qwave // qwave same size as output eqmap, q_vacuum values in 1/A
wave structwave //3 columns: {thickness (A), roughnes (A), index of refractindex} // start with substrate=0 and up to top layer (do not include vacuum)
wave/c optical_const // size = num_energies x num_layers; to prepare use : optical_const = refractindex(ewave[p])[structwave[q][2]]
variable mode //

//wave/c refractindex //N columns : each contain energy dependence of comoplex refr index of a compound
// 0: Parratt, 
// 1: no denominator
// 5: use average kz AND fresnel ~ dn*k0^2/kz_avg^2 ---- good approx, FT -like
eqmap=0

make/o/n=(dimsize(eqmap,0), dimsize(eqmap,1)) cos_alpha_vac
make/o/c/n=(dimsize(eqmap,0), dimsize(eqmap,1), dimsize(structwave,0)) cosalpha, kzwave, phshift=0, roughdamping, Fresnelcoef=0, Reflection=0
//make/o/c/n=(dimsize(eqmap,0), dimsize(structwave,0)) optical_const //per layer

//optical_const = refractindex(ewave[p])[structwave[q][2]] // linearly interpolate optical constants per layer

variable hbar_c = 1973.26972
cos_alpha_vac = 0.5*hbar_c*qwave[p][q]/ewave[p]
duplicate/o cos_alpha_vac, sin_alpha_vac
sin_alpha_vac = (1-cos_alpha_vac^2)^.5

variable m,k, numenergy
make/c/o/n=(dimsize(structwave,0)) temp
//layer-specific quantities
if (mode<2)
	cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
	kzwave[][][] = optical_const[p][r]*cosalpha[p][q][r]/cos_alpha_vac[p][q]*qwave[p][q]/2 
elseif (mode==5)
	cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
	kzwave[][][] = optical_const[p][r]*cosalpha[p][q][r]/cos_alpha_vac[p][q]*qwave[p][q]/2 
	for (k=0; k<dimsize(cosalpha,0); k+=1)
		for (m=0; m<dimsize(cosalpha,1); m+=1)
			temp = kzwave[k][m][p]
			kzwave[k][m][] = mean(temp)
		endfor
	endfor	
endif
phshift[][][1,dimsize(structwave,0)-1] = exp(2*structwave[r][0]*kzwave[p][q][r]*cmplx(0,1)) // exp(2*kz*d*i)

//interface-specific quantities (m'th value = interface between layers m and m+1
roughdamping[][][0,dimsize(structwave,0)-2] = exp(-2*kzwave[p][q][r]*kzwave[p][q][r+1]*structwave[r][1]^2)
roughdamping[][][dimsize(structwave,0)-1] = exp(-2*kzwave[p][q][r]*qwave[p][q]/2*structwave[r][1]^2) //top layer and vacuum k_z_vac==q_z/2=y/2
//Fresnel coef. : calculation for sigma light
if (mode==5)
	Fresnelcoef[][][0,dimsize(structwave,0)-2] = 0.5*(ewave[p]/hbar_c)^2/kzwave[p][q][r]^2*roughdamping[p][q][r]*(-optical_const[p][r] + optical_const[p][r+1])
	Fresnelcoef[][][dimsize(structwave,0)-1] = 0.5*(ewave[p]/hbar_c)^2/kzwave[p][q][r]^2*roughdamping[p][q][r]*(1-optical_const[p][r]) //top layer and vacuum
else
	Fresnelcoef[][][0,dimsize(structwave,0)-2] = roughdamping[p][q][r]*(-optical_const[p][r]*cosalpha[p][q][r] + optical_const[p][r+1]*cosalpha[p][q][r+1])/(optical_const[p][r]*cosalpha[p][q][r] + optical_const[p][r+1]*cosalpha[p][q][r+1])
	Fresnelcoef[][][dimsize(structwave,0)-1] = roughdamping[p][q][r]*(cos_alpha_vac[p][q]-optical_const[p][r]*cosalpha[p][q][r])/(cos_alpha_vac[p][q] + optical_const[p][r]*cosalpha[p][q][r]) //top layer and vacuum
endif

// reflection between substrate and first layer
Reflection[][][0] = Fresnelcoef[p][q][0]
for (m=1; m<dimsize(structwave,0); m=m+1) //start with dimsize(structwave,0)-2 since no reflection from bottom of substrate
	if (mode==0)
		Reflection[][][m] = (Fresnelcoef[p][q][r] + Reflection[p][q][r-1]*phshift[p][q][r])/(1+Fresnelcoef[p][q][r]*Reflection[p][q][r-1]*phshift[p][q][r])
	else
		Reflection[][][m] = (Fresnelcoef[p][q][r] + Reflection[p][q][r-1]*phshift[p][q][r])
	endif
endfor

//eqmap =log(cabs(Reflection[p][q][dimsize(structwave,0)-1]))*2
eqmap =cabs(Reflection[p][q][dimsize(structwave,0)-1]*qwave[p][q]^2)^2

//duplicate/o eqmap $nameofwave(eqmap)+"_phase"
//wave eqmaphase=$nameofwave(eqmap)+"_phase"
//eqmaphase = 0
//eqmaphase = imag(r2polar(Reflection[p][q][dimsize(structwave,0)-1]))

end

//=======================================================
function Parratt_MatrixMultmode5(eqmap, structwave, refractindex) //will do same as parratt mode=5, but as a matrix multiplication
wave eqmap //dimensions must be specified: x-energy in eV; y-q_z in 1/Angstrom
wave structwave //3 columns: {thickness (A), roughnes (A), index of refractindex} // start with substrate=0 and up to top layer (do not include vacuum)
wave/c refractindex //N columns : each contain energy dependence of comoplex refr index of a compound

eqmap=0

make/o/n=(dimsize(eqmap,0), dimsize(eqmap,1)) cos_alpha_vac
make/o/c/n=(dimsize(eqmap,0), dimsize(eqmap,1), dimsize(structwave,0)) cosalpha, phshift=0
make/o/c/n=(dimsize(eqmap,0), dimsize(eqmap,1)) kzwave1, Reflection=0

make/o/n= (dimsize(structwave,0)) interface_z=0
make/o/c/n=(dimsize(eqmap,0), dimsize(structwave,0)) optical_const, optical_const_diff //per layer And per interface

setscale/P x, dimoffset(eqmap,0), dimdelta(eqmap,0), Reflection, cosalpha, phshift, kzwave1, cos_alpha_vac, optical_const, optical_const_diff
setscale/P y, dimoffset(eqmap,1), dimdelta(eqmap,1), Reflection, cosalpha, phshift, kzwave1, cos_alpha_vac

optical_const = refractindex(x)[structwave[q][2]] // linearly interpolate optical constants per layer
optical_const_diff[][0, dimsize(optical_const_diff,1)-2] = optical_const[p][q+1] - optical_const[p][q]
optical_const_diff[][ dimsize(optical_const_diff,1)-1] = 1 - optical_const[p][q]

variable m,k
for (k=dimsize(interface_z,0)-2; k>=0; k-=1)
	interface_z[k] = interface_z[k+1]+structwave[k+1][0]
endfor

variable hbar_c = 1973.26972
cos_alpha_vac = 0.5*hbar_c*y/x
duplicate/o cos_alpha_vac, sin_alpha_vac
sin_alpha_vac = (1-cos_alpha_vac^2)^.5

make/c/o/n=(dimsize(structwave,0)) temp
//average kz in medium
	cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
	for (k=0; k<dimsize(cosalpha,0); k+=1) //loop energies
		for (m=0; m<dimsize(cosalpha,1); m+=1) //loop momenta
			temp = optical_const[k][p]*cosalpha[k][m][p]/cos_alpha_vac[k][m]
			kzwave1[k][m] = mean(temp)*y/2
		endfor
	endfor	

//interface-specific quantities (m'th value = interface between layers m and m+1
phshift[][][0,dimsize(structwave,0)-2] = exp(-2*kzwave1[p][q]^2*structwave[r][1]^2 + 2*interface_z[r]*kzwave1[p][q]*cmplx(0,1))
phshift[][][dimsize(structwave,0)-1] = exp(-2*kzwave1[p][q]*y/2*structwave[r][1]^2) //top layer & vacuum

for (k=0; k<dimsize(Reflection,0); k+=1) //loop energies
	for (m=0; m<dimsize(Reflection,1); m+=1) //loop momenta
		temp = phshift[k][m][p]*optical_const_diff[k][p]
		Reflection[k][m] = 0.5*(x/hbar_c)^2/kzwave1[p][q]^2*sum(temp)
	endfor
endfor

//eqmap =log(cabs(Reflection[p][q]))*2
eqmap =log(cabs(Reflection[p][q]))*2
duplicate/o eqmap $nameofwave(eqmap)+"_phase"
wave eqmaphase=$nameofwave(eqmap)+"_phase"
eqmaphase = 0
eqmaphase = imag(r2polar(Reflection[p][q][dimsize(structwave,0)-1]))
end
//=======================================================

function Parratt_MatrixMultmode7(eqmap, structwave, refractindex) //will do same as parratt mode=7, but as a matrix multiplication
wave eqmap //dimensions must be specified: x-energy in eV; y-q_z in 1/Angstrom
wave structwave //3 columns: {thickness (A), roughnes (A), index of refractindex} // start with substrate=0 and up to top layer (do not include vacuum)
wave/c refractindex //N columns : each contain energy dependence of comoplex refr index of a compound

eqmap=0

make/o/n=(dimsize(eqmap,0), dimsize(eqmap,1)) cos_alpha_vac
make/o/c/n=(dimsize(eqmap,0), dimsize(eqmap,1), dimsize(structwave,0)) cosalpha, phshift=0
make/o/c/n=(dimsize(structwave,0), dimsize(eqmap,1)) phshift=0
make/o/c/n=(dimsize(eqmap,0), dimsize(eqmap,1)) Reflection=0
make/o/c/n=(dimsize(eqmap,1)) kzwave2

make/o/n= (dimsize(structwave,0)) interface_z=0
make/o/c/n=(dimsize(eqmap,0), dimsize(structwave,0)) optical_const, optical_const_diff, temp2 //per layer And per interface

setscale/P x, dimoffset(eqmap,0), dimdelta(eqmap,0), Reflection, cosalpha, cos_alpha_vac, optical_const, optical_const_diff, temp2
setscale/P y, dimoffset(eqmap,1), dimdelta(eqmap,1), Reflection, cosalpha, phshift, cos_alpha_vac
//setscale/P x, dimoffset(eqmap,1), dimdelta(eqmap,1), kzwave2

variable m,k, qval
for (k=dimsize(interface_z,0)-2; k>=0; k-=1)
	interface_z[k] = interface_z[k+1]+structwave[k+1][0]
endfor

variable hbar_c = 1973.26972
cos_alpha_vac = 0.5*hbar_c*y/x
duplicate/o cos_alpha_vac, sin_alpha_vac
sin_alpha_vac = (1-cos_alpha_vac^2)^.5

optical_const = refractindex(x)[structwave[q][2]] // linearly interpolate optical constants per layer
optical_const_diff[][0, dimsize(optical_const_diff,1)-2] = (optical_const[p][q+1] - optical_const[p][q])*0.5*(x/hbar_c)^2
optical_const_diff[][ dimsize(optical_const_diff,1)-1] = (1 - optical_const[p][q])*0.5*(x/hbar_c)^2

//average kz in medium & energy
cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
for (m=0; m<dimsize(eqmap,1); m+=1)
	temp2 =  optical_const[p][q]*cosalpha[p][m][q]/cos_alpha_vac[p][m]
	qval = (dimoffset(eqmap,1)+ m*dimdelta(eqmap,1))/2
	kzwave2[m] = mean(temp2)*qval
	//
endfor

//interface-specific quantities (m'th value = interface between layers m and m+1
phshift[0,dimsize(structwave,0)-2][] = exp(-2*kzwave2[q]^2*structwave[p][1]^2 + 2*interface_z[p]*kzwave2[q]*cmplx(0,1))/kzwave2[q]^2
phshift[dimsize(structwave,0)-1][] = exp(-2*kzwave2[q]*y/2*structwave[p][1]^2)/kzwave2[q]^2 //top layer & vacuum

MatrixOp/C/O/S Reflection = optical_const_diff x phshift

eqmap =log(cabs(Reflection[p][q]))*2
duplicate/o eqmap $nameofwave(eqmap)+"_phase"
wave eqmaphase=$nameofwave(eqmap)+"_phase"
eqmaphase = 0
eqmaphase = imag(r2polar(Reflection[p][q][dimsize(structwave,0)-1]))
end
//=======================================================
function/c Parratt_phasor(k_z, interface_z, roughness, engy)
variable interface_z, roughness, engy
variable/c k_z
return 0.5*(engy/1973.26972)^2/k_z^2*exp(-2*k_z^2*roughness^2 + 2*interface_z*k_z*cmplx(0,1))
end
//=======================================================
function Parratt_struct_fit(eqscans, energyvals, structwave, refractindex)
//mode 1 used

wave eqscans //dimension y must be defined: q_z in 1/Angstrom; energy values given by enegrgyvals wave
wave energyvals //1D wave with energy values
wave structwave //3 columns: {thickness (A), roughnes (A), index of refractindex} // start with substrate=0 and up to top layer (do not include vacuum)
wave/c refractindex //N columns : each contain energy dependence of comoplex refr index of a compound

//eqscans=0
make/o/n=(dimsize(eqscans,0), dimsize(eqscans,1)) cos_alpha_vac, eqscans_deriv
make/o/c/n=(dimsize(eqscans,0), dimsize(eqscans,1), dimsize(structwave,0)) cosalpha, kzwave, phshift=0, roughdamping, Fresnelcoef=0, Fresnelcoef0=0, Reflection=0

make/o/c/n=(dimsize(eqscans,0), dimsize(structwave,0)) optical_const //per layer

setscale/P y, dimoffset(eqscans,1), dimdelta(eqscans,1), Fresnelcoef, Fresnelcoef0, Reflection, cosalpha, phshift, kzwave, cos_alpha_vac, roughdamping
duplicate/O/C Reflection, dReflection_dRough, dReflection_dThick //derivatives of top layer Reflection wrt layer thicknesses and roughnesses
duplicate/O/C phshift, phshiftint //integrated 

make/o/n=(dimsize(Reflection,0),dimsize(Reflection,1), dimsize(structwave,0)) dI_dRough, dI_dThick

optical_const = refractindex(energyvals[p])[structwave[q][2]] // linearly interpolate optical constants per layer

//layer-specific quantities
variable hbar_c = 1973.26972
cos_alpha_vac = 0.5*hbar_c*y/energyvals[p]
duplicate/o cos_alpha_vac, sin_alpha_vac
sin_alpha_vac = (1-cos_alpha_vac^2)^.5
variable m,k, numenergy
make/c/o/n=(dimsize(structwave,0)) temp
cosalpha[][][] = (1-(sin_alpha_vac[p][q]/optical_const[p][r])^2)^0.5
kzwave[][][] = optical_const[p][r]*cosalpha[p][q][r]/cos_alpha_vac[p][q]*y/2 

Fresnelcoef0[][][0,dimsize(structwave,0)-2] = (-optical_const[p][r]*cosalpha[p][q][r] + optical_const[p][r+1]*cosalpha[p][q][r+1])/(optical_const[p][r]*cosalpha[p][q][r] + optical_const[p][r+1]*cosalpha[p][q][r+1])
Fresnelcoef0[][][dimsize(structwave,0)-1] = (cos_alpha_vac[p][q]-optical_const[p][r]*cosalpha[p][q][r])/(cos_alpha_vac[p][q] + optical_const[p][r]*cosalpha[p][q][r]) //top layer and vacuum
duplicate/o/c Fresnelcoef0, Fresnelcoef

end
//------------------------
Function Parratt_struct_grad() // calculates gradient of the objective function = (I_expt-I_theo)^2 wrt struct parameters
wave structwave = structwaveSL
wave/c kzwave
wave/c phshift
wave/c phshiftint
wave/c Fresnelcoef0
wave/c Fresnelcoef
wave/c Reflection
wave/c dReflection_dRough
wave/c dReflection_dThick
wave dI_dRough
wave dI_dThick
wave/c roughdamping

//phshift[][][1,dimsize(structwave,0)-1] = exp(2*structwave[r][0]*kzwave[p][q][r]*cmplx(0,1)) // exp(2*kz*d*i)
phshift[][][] = exp(2*cmplx(0,1)*structwave[r][0]*kzwave[p][q][r]) // exp(2*kz*d*i)

//interface-specific quantities (m'th value = interface between layers m and m+1
roughdamping[][][0,dimsize(structwave,0)-2] = exp(-2*kzwave[p][q][r]*kzwave[p][q][r+1]*structwave[r][1]^2)
roughdamping[][][dimsize(structwave,0)-1] = exp(-2*kzwave[p][q][r]*y/2*structwave[r][1]^2) //top layer and vacuum k_z_vac==q_z/2=y/2

Fresnelcoef[][][] = roughdamping[p][q][r]*Fresnelcoef0[p][q][r]

// reflection between substrate and first layer
Reflection[][][0] = Fresnelcoef[p][q][0]

phshiftint[][][dimsize(structwave,0)-1] = 1
variable m
for (m=dimsize(structwave,0)-2; m>=0; m=m-1) 
	phshiftint[][][m] = phshiftint[p][q][r+1]* phshift[p][q][r+1]
endfor

for (m=1; m<dimsize(structwave,0); m=m+1) 
	Reflection[][][m] = Fresnelcoef[p][q][r] + Reflection[p][q][r-1]*phshift[p][q][r]
endfor

dReflection_dRough[][][dimsize(structwave,0)-1] = -2*kzwave[p][q][r]*y/2*Fresnelcoef[p][q][r]*phshiftint[p][q][r]
dReflection_dRough[][][0,dimsize(structwave,0)-2] = -2*kzwave[p][q][r]*kzwave[p][q][r+1]*Fresnelcoef[p][q][r]*phshiftint[p][q][r]
dReflection_dThick[	][][1,dimsize(structwave,0)-1] = 2*cmplx(0,1)*kzwave[p][q][r]*phshiftint[p][q][r-1]*Reflection[p][q][r-1]

dI_dRough[][][] = 2*real(conj(Reflection[p][q][dimsize(structwave,0)-1])*dReflection_dRough[p][q][r])
dI_dThick[][][] = 2*real(conj(Reflection[p][q][dimsize(structwave,0)-1])*dReflection_dThick[p][q][r])

//eqmap =log(cabs(Reflection[p][q][dimsize(structwave,0)-1]))*2

End
//=======================================================
Function parratt_myfunc1(w,xw)
	Wave w, xw
	//xw[0] - roughness 0
	//xw[1] - roughness 4
	//xw[2] - roughness 32
	//xw[3] - roughness 66
	//xw[4] - roughness 98
	//xw[5] - roughness 105
	//xw[6] - overall scale
	//xw[7] - z 4
	//xw[8] - z 32
	//xw[9] - z 66
	//xw[10] - z 98
	//xw[11] - z 105
	wave refldata = refl51n
	wave refldataQ = refl51Q
	wave kz = kz4fit //kz4fit = kzwave2res[4](refl56Q[p])
	duplicate/o refldata, reflmodel, refmodresidue
	make/c/o/n=(dimsize(refldata,0)) resph_a, resph_b, resph_c, resph_d, resph_e
	
	variable engy = 1050
	wave refracindex = refracindex1
	variable/c nSTO = refracindex(engy)[5]
	variable/c nNTO = refracindex(engy)[6]
	variable/c nLSAT = refracindex(engy)[0]
	variable/c nSrO = refracindex(engy)[1]
	
	variable/c Dn_a = cmplx(real(1-nSTO)*xw[12], imag(1-nSTO)*xw[13])
	variable/c Dn_b = cmplx(real(1-nSTO)*xw[14], imag(1-nSTO)*xw[15])
	variable/c Dn_c = cmplx(real(nSTO-nNTO)*xw[16], imag(nSTO-nNTO)*xw[17])
	variable/c Dn_d = cmplx(real(nSTO-nLSAT)*xw[18], imag(1-nSTO)*xw[19])
	variable/c Dn_e = cmplx(real(1-nSTO)*xw[20], imag(1-nSTO)*xw[21])
	
	resph_a = Parratt_phasor(kz[p], 0, xw[0], engy)
	resph_b = Parratt_phasor(kz[p], 4+xw[7], xw[1], engy)
	resph_c = Parratt_phasor(kz[p], 32+xw[8], xw[2], engy) - Parratt_phasor(kz[p], 66+xw[9], xw[3], engy)
	resph_d = Parratt_phasor(kz[p], 98+xw[10], xw[4], engy)
	resph_e = Parratt_phasor(kz[p], 105+xw[11], xw[5], engy)
//	reflmodel = xw[8]*cabs(cmplx(dSTO,bSTO)*resph_a[p] + cmplx(dNTOSTO*xw[4],bNTOSTO*xw[5])*resph_b[p] + cmplx(dLSATSTO*xw[6],bLSATSTO*xw[7])*resph_c[p])^2
	reflmodel = xw[6]*cabs(Dn_a*resph_a[p] + Dn_b*resph_b[p] + Dn_c*resph_c[p] + Dn_d*resph_d[p] + Dn_e*resph_e[p])^2*refldataQ[p]^4
	refmodresidue = (reflmodel[p] - 2e-4*refldata[p])^2*1e9
	refmodresidue[0,6]=0
	
	return sum(refmodresidue)
End
//=======================================================
Function parratt_mode1_myfunc1(w,xw)
	Wave w, xw

	make/o/n=(6,3) structwave_fit
	structwave_fit[5][1] = xw[0]
	structwave_fit[4][1] = xw[1]
	structwave_fit[3][1] = xw[2]
	structwave_fit[2][1] = xw[3]
	structwave_fit[1][1] = xw[4]
	structwave_fit[0][1] = xw[5]
	structwave_fit[1][0] = 105+xw[11]-(98+xw[10])
	structwave_fit[2][0] = 98+xw[10] - 66-xw[9]
	structwave_fit[3][0] = 66+xw[9] - 32 - xw[8]
	structwave_fit[4][0] = 32 + xw[8] - 4 - xw[7]
	structwave_fit[5][0] = 4 + xw[7]
	structwave_fit[0,1][2] = 0
	structwave_fit[2][2] = 5
	structwave_fit[4,5][2] = 5
	structwave_fit[3][2] = 6
	
	variable engy = 1050
	wave refracindex = refracindex1
	wave/c opt_const_fit0, opt_const_fit
	//make/o/c/n=(1,6) opt_const_fit, opt_const_fit0 = refracindex(engy)[structwave_fit[q][2]]-1
	opt_const_fit[0][0] = 1+opt_const_fit0[p][q]
	opt_const_fit[0][1,5] = 1 + cmplx(xw[12+2*(q-1)]*real(opt_const_fit0[0][q]), xw[13+2*(q-1)]*imag(opt_const_fit0[0][q]))
	
	wave refldata = refl51n //maybe 2D?
	wave refldataQ = refl51Q
	duplicate/o refldata, reflmodel, refmodresidue
	make/o/n=(1,dimsize(refldata,0)) eqmap_fit, qwave_fit
	make/o/n=(1) ewave_fit
	ewave_fit[0][0] = engy
	qwave_fit = max(refldataQ[q],1e-3)
	parratt_simEQW(eqmap_fit, ewave_fit, qwave_fit, structwave_fit, opt_const_fit, 1)
	
	reflmodel = xw[6]*eqmap_fit[0][p]
	refmodresidue = (reflmodel[p] - 2e-4*refldata[p])^2*1e9
//	refmodresidue = (log(reflmodel[p]) - log(2e-4*refldata[p]))^2*1e9
	
	return sum(refmodresidue)
End
//=======================================================
