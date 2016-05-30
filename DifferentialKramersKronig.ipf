#pragma rtGlobals=1		// Use modern global access method.


Function SelfEdkk(Expf1, Expf2)
	Wave Expf1, Expf2	// The experimental f2 data, non-resonant f1 and f2, and energy wave (for all three f's)

	Variable i, j
	
	Duplicate/O Expf2, NRf1, NRf2, Energy
	Energy=x
	NRf1=0
	NRf2=0
	// Make sure all waves are the same size
	Variable wavesize = DimSize(Expf2, 0)
	If(DimSize(Expf1, 0) == wavesize)
	Else
		Print "Error: All input waves must be the same length. Make sure you've interpolated the non-resonant data onto the experimental energy wave."
		return 0
	EndIf

	// Make the output wave and linear interpolation waves
	Make /O /N = (wavesize) bi, mi
	bi = 0
	mi = 0
	
	// The output wave will be slightly offset in energy from the input to avoid numerical instabilities
//	Make /O /N = (wavesize) Expf1, Energy_out
	Make /O /N = (wavesize) Energy_out
	Expf1 = 0
	Energy_out = 0
	
	//SetScale /I x, 0, (wavesize-1), Energy_out
	//Interpolate2 /I=3  /t=1 /Y = Energy_out Energy
	Energy_out = Energy[p+0.0001]
	
	// Calculate the linear interpolation coefficients
	mi[0, wavesize-2] = ((Expf2[p+1] - NRf2[p+1]) - (Expf2[p] - NRf2[p]))/(Energy[p+1] - Energy[p])
	bi[0, wavesize-2] = (Expf2[p] - NRf2[p]) - mi[p] * Energy[p]
      
      // Perform the Diff KK transform (where I linearly interpolate between the supported points)
	//For(j =0; j < wavesize-1; j+=1)
	For(j =0; j < wavesize-1; j+=1)
			Expf1 += mi[j] *Energy[j] - mi[j] *Energy[j+1]
			Expf1 += 0.5*Energy_out[p] * mi[j] * Ln(Abs((1+(Energy[j+1]/Energy_out[p]))/(1-(Energy[j+1]/Energy_out[p])))) - 0.5*Energy_out[p] * mi[j] * Ln(Abs((1+(Energy[j]/Energy_out[p]))/(1-(Energy[j]/Energy_out[p]))))
			Expf1 += 0.5*bi[j] * Ln(Abs(Energy_out[p]^2 - Energy[j]^2)) - 0.5*bi[j] * Ln(Abs(Energy_out[p]^2 - Energy[j+1]^2))
			
      	EndFor
      	
      	Expf1 *= 2/Pi
      	Expf1 += NRf1[p]

	KillWaves bi, mi, Energy_out, NRf1, NRf2

End
//=====================================================
function mySEDKK(SE1, Yshift, NN, cont)
wave SE1 //XY -wave
variable Yshift
variable NN, cont //cont=0 - continue with zeroes; 1-constant; 2-linear descent to zero

make/O/n=(dimsize(SE1,0)) se1x, se1y
se1x=SE1[p][0]
se1y=SE1[p][1] + Yshift
sort se1x, se1y, se1x
make/O/n=(NN) se1int
setscale/I x, se1x[0], 0, se1int
Interpolate2/T=1/I=3/Y= se1int se1x, se1y
make/O/n=(NN*6) se1intfull
setscale/I x, 3*se1x[0], -3*se1x[0], se1intfull
if (cont==0)
	se1intfull[0,799] = 0
elseif (cont==1)
	se1intfull[0,2*NN-1] = se1y[0]
elseif (cont==2)
	se1intfull[0,799] = p/799*se1y[0]
endif
se1intfull[2*NN,3*NN-1] = se1int[p-2*NN]
se1intfull[3*NN,6*NN-1] = -se1intfull(-x)

duplicate/o se1intfull, se2intfull, integrsum

variable kk, x0
for (kk=0; kk<dimsize(se1intfull,0);kk+=1)
	x0 = dimoffset(se1intfull,0) + kk*DimDelta(se1intfull,0)
	integrsum = se1intfull(x)*Ln(Abs(x - x0 +0.5*DimDelta(se1intfull,0))/Abs(x - x0 - 0.5*DimDelta(se1intfull,0)))
	se2intfull[kk] = 1/pi*sum(integrsum)
endfor
x0 = se2intfull(0)
se2intfull -= x0

KillWaves se1x, se1y, se1int, integrsum
End
//=====================================================