#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//===========================================
Function myFSsym(kFS, xscaling, yscaling, X0, X1, Y0, Y1, dx, dy)
Wave kFs
Variable xscaling, yscaling, X0, X1, Y0, Y1, dx,dy// scaling and range of points to use

//duplicate/O kFS, kfsout
duplicate/O kFS, $"kFSsym"
wave kfsout = $"kFSsym"
kfsout[][0] = kfs[p][q]/xscaling + dx
kfsout[][1] = kfs[p][q]/yscaling + dy
//
variable kk, mm=0
for (kk=0; kk< dimsize(kFS,0); kk+=1)
	if (numtype(kfsout[mm][0])==2)
		DeletePoints mm, 1, kfsout
	elseif (sign(kfsout[mm][0]-X0)+sign(X1-kfsout[mm][0]) + sign(kfsout[mm][1]-Y0)+sign(Y1-kfsout[mm][1])<4)
		DeletePoints mm, 1, kfsout
	else
		mm+=1
	endif
endfor
//to first zone
kfsout = kfsout[p][q] - round(kfsout[p][q]/2)*2
//to first quadrant
kfsout = abs(kfsout[p][q])

End
//===========================================
function mySymFS2TB()
//variable bandind0 // 0- alpha, 1-beta, 2- gamma
//NVAR bandind = bandind0
make/O/n=3 TBparout

//Optimize /M = {3, 0}/X={-1.1, 0.1, 0.0} myBandFunc, TBparout
Optimize/M={0,0}/X={1.05, 0.15, 0.05} myBandFunc, TBparout

//print mySRO_TBFSVOL(200, 200, 2)

end
//===========================================
Function myBandFunc(w,x1, x2, x3)
	Wave w
	Variable x1, x2, x3 // x1=mu1d/mu2d; x2 = t3/t4; x3 = x5/0

	Wave kXY = $"kFSsym"
	NVAR bandind
	make/O/n=(dimsize(kXY,0)) tmptmp
	if (bandind==0)
		tmptmp = -x1-(1+x2)*(Cos(pi*kXY[p][0])+Cos(pi*kXY[p][1]))-0.5*sqrt(4*(1+x2)^2*(Cos(pi*kXY[p][0])+Cos(pi*kXY[p][1]))^2 - 8*(2*x2 - 2*x3^2+(x2+2*x3^2)*Cos(2*kXY[p][0]*pi) +(1+x2^2)*Cos(pi*(kXY[p][0]-kXY[p][1])) - x3^2*Cos(2*pi*(kXY[p][0]-kXY[p][1])) + (x2 + 2*x3^2)*Cos(2*pi*kXY[p][1]) + (1+x2^2)*Cos(pi*kXY[p][0]+pi*kXY[p][1]) - x3^2*Cos(2*pi*(kXY[p][0]+kXY[p][1])) ))
		tmptmp = tmptmp[p]^2 + 0.02*x2^2 + 5*x3^2
	elseif (bandind==1)
		tmptmp = -x1-(1+x2)*(Cos(pi*kXY[p][0])+Cos(pi*kXY[p][1]))+0.5*sqrt(4*(1+x2)^2*(Cos(pi*kXY[p][0])+Cos(pi*kXY[p][1]))^2 - 8*(2*x2 - 2*x3^2+(x2+2*x3^2)*Cos(2*kXY[p][0]*pi) +(1+x2^2)*Cos(pi*(kXY[p][0]-kXY[p][1])) - x3^2*Cos(2*pi*(kXY[p][0]-kXY[p][1])) + (x2 + 2*x3^2)*Cos(2*pi*kXY[p][1]) + (1+x2^2)*Cos(pi*kXY[p][0]+pi*kXY[p][1]) - x3^2*Cos(2*pi*(kXY[p][0]+kXY[p][1])) ))
		tmptmp = tmptmp[p]^2 + 0.03*x2^2 + 0.03*x3^2
	elseif (bandind==2)
		tmptmp = (-2*(Cos(Pi*kXY[p][0])+Cos(Pi*kXY[p][1])) - 4*x2*Cos(Pi*kXY[p][0])*Cos(Pi*kXY[p][1]) - x1)^2 + x3^2
	endif
	return sum(tmptmp)	// an expression...
End
//===========================================
Function mySRO_TB()

Make/O/n=(200,200, 2) tbbandout
Make/O/n=(200,200) tbbandoutIMG
SetScale/I x, -1, 1, tbbandout, tbbandoutIMG
SetScale/I y, -1, 1, tbbandout, tbbandoutIMG
SetScale/I z, -1, 1, tbbandout, tbbandoutIMG

wave tbpars = W_Extremum
NVAR bandind
Variable t1, t2, t3, t4, t5, mu1d, mu2d
t1 = 1
t2=1
mu1d = tbpars[0]
mu2d = tbpars[0]
t3 = tbpars[1]
t4 = tbpars[1]
t5= tbpars[2]

if (bandind==0)
//	tbbandout = -mu1d - t2*Cos(Pi*x) - t3*Cos(Pi*x) - t2*Cos(Pi*y) -  t3*Cos(Pi*y) - (1/Sqrt(2))*Sqrt(Abs(2*t2^2 - 4*t2*t3 + 2*t3^2 + 8*t5^2 - 4*(t2 - t3)^2*Cos(Pi*x)*Cos(Pi*y) + t2^2*Cos(2*Pi*y) - 2*t2*t3*Cos(2*Pi*y) + t3^2*Cos(2*Pi*y) - 8*t5^2*Cos(2*Pi*y) + Cos(2*Pi*x)*(t2^2 - 2*t2*t3 + t3^2 - 8*t5^2 + 8*t5^2*Cos(2*Pi*y) )) )
	tbbandout = -mu1d-(t2+t3)*(Cos(pi*x)+Cos(pi*y))-0.5*sqrt(4*(t2+t3)^2*(Cos(pi*x)+Cos(pi*y))^2 - 8*(2*t2*t3 - 2*t5^2+(t2*t3+2*t5^2)*Cos(2*x*pi) +(t2^2+t3^2)*Cos(pi*(x-y)) - t5^2*Cos(2*pi*(x-y)) + (t2*t3 + 2*t5^2)*Cos(2*pi*y) + (t2^2+t3^2)*Cos(pi*x+pi*y) - t5^2*Cos(2*pi*(x+y)) ))
elseif (bandind==1)
//	tbbandout = -mu1d-(t2+t3)*(Cos(pi*x)+Cos(pi*y))+0.5*sqrt(4*(t2+t3)^2*(Cos(pi*x)+Cos(pi*y))^2 - 8*(2*t2*t3 - 2*t5^2+(t2*t3+2*t5^2)*Cos(2*x*pi) +(t2^2+t3^2)*Cos(pi*(x-y)) - t5^2*Cos(2*pi*(x-y)) + (t2*t3 + 2*t5^2)*Cos(2*pi*y) + (t2^2+t3^2)*Cos(pi*x+pi*y) - t5^2*Cos(2*pi*(x+y)) ))
	tbbandout = -mu1d - t2*Cos(Pi*x) - t3*Cos(Pi*x) - t2*Cos(Pi*y) - t3*Cos(Pi*y) + (1/Sqrt(2))*Sqrt(Abs(2*t2^2 - 4*t2*t3 + 2*t3^2 + 8*t5^2 - 4*(t2 - t3)^2*Cos(Pi*x)*Cos(Pi*y) + t2^2*Cos(2*Pi*y) - 2*t2*t3*Cos(2*Pi*y) + t3^2*Cos(2*Pi*y) - 8*t5^2*Cos(2*Pi*y) + Cos(2*Pi*x)*(t2^2 - 2*t2*t3 + t3^2 - 8*t5^2 + 8*t5^2*Cos(2*Pi*y))))
elseif (bandind==2)
	tbbandout = -2*t1*(Cos(Pi*x)+Cos(Pi*y)) - 4*t4*Cos(Pi*x)*Cos(Pi*y) - mu2d
endif

tbbandoutIMG = exp(-(tbbandout[p][q][0]/0.02)^2)

End
//===========================================
Function mySRO_TBFSVOL(Nx, Ny, Nz)
Variable Nx,Ny,Nz
NVAR bandind
wave W_Extremum

SRO_prepwaves("tmptb", Nx,Ny,Nz)
Wave wvEa = $"tmptbEa"
Wave wvEb = $"tmptbEb"
Wave wvEg = $"tmptbEg"

if (bandind==0)
	SRO_tightbinding2D(wvEa, wvEb, wvEg, 1, 1, W_Extremum[1], 0, W_Extremum[2], W_Extremum[0], 0)
	return imag(DNosX(wvEa, 0))
elseif (bandind==1)
	SRO_tightbinding2D(wvEa, wvEb, wvEg, 1, 1, W_Extremum[1], 0, W_Extremum[2], W_Extremum[0], 0)
	return imag(DNosX(wvEb, 0))
elseif (bandind==2)
	SRO_tightbinding2D(wvEa, wvEb, wvEg, 1, 1, 0, W_Extremum[1], 0, 0, W_Extremum[0])
	return imag(DNosX(wvEg, 0))
endif

End
//===========================================
function mySROtbfitlist(list, dx,dy)
wave/T list
variable dx,dy
make/o/n=(dimsize(list,0),4) tbfitout

variable kk

for (kk=0; kk<dimsize(list,0);kk+=1)
	myFSsym($list(kk), 10, 10, -100, 100, -100, 100,dx,dy);
	mySymFS2TB();
	wave W_Extremum
	tbfitout[kk][0,2] = abs(W_Extremum[q])
	tbfitout[kk][3] = mySRO_TBFSVOL(200, 200, 2)
endfor

end
//===========================================
function mylistplot(list)
wave/T list
variable kk
for (kk=0; kk<dimsize(list,0);kk+=1)
	wave tmpout=$list[kk]
	appendtograph tmpout[][1] vs tmpout[][0]
endfor

end
//============================================
function myplottbfs(tbpar, bandind)
wave tbpar
variable bandind




end
//============================================
function myTBGdnos(t1, t4, mu2, t2n, Nx, engy)
//myTBGdnos(0.081, 0.039, 0.128, 0.005, 200, engy) - borisenko2013
//myTBGdnos(0.0694, 0.0289, 0.105, 0.00, 200, engy) - kms2007
//myTBGdnos(0.0875, 0.036, 0.130, 0.00, 200, engy) - my

wave engy
variable t1, t4, mu2, t2n, Nx //Ny=Nx, Nz=2; t2n- cos(2x) as in borisenko2013

make/o/n=(Nx, Nx, 2) tmpEg
setscale/I x, -1, 1, tmpEg
setscale/I y, -1, 1, tmpEg
setscale/I z, -1, 1, tmpEg

tmpEg = -mu2 - 2*t1*(Cos(Pi*x)+Cos(Pi*y)) - 4*t4*Cos(Pi*x)*Cos(Pi*y) - 2*t2n*(Cos(2*Pi*x)+Cos(2*Pi*y))
return 0
make/o/n=(dimsize(engy,0))/C TBGdnosOUT

TBGdnosOUT = DNosX(tmpEg, engy[p])

killwaves tmpEg
end
//============================================
function/C myTBGdnosPT(t1, t4, mu2, t2n, Nx, engy)
//myTBGdnos(0.081, 0.039, 0.128, 0.005, 200, engy) - borisenko2013
//myTBGdnos(0.0694, 0.0289, 0.105, 0.00, 200, engy) - kms2007
//myTBGdnos(0.0875, 0.036, 0.130, 0.00, 200, engy) - my

variable engy
variable t1, t4, mu2, t2n, Nx //Ny=Nx, Nz=2; t2n- cos(2x) as in borisenko2013

make/o/n=(Nx, Nx, 2) tmpEg
setscale/I x, -1, 1, tmpEg
setscale/I y, -1, 1, tmpEg
setscale/I z, -1, 1, tmpEg

tmpEg = -mu2 - 2*t1*(Cos(Pi*x)+Cos(Pi*y)) - 4*t4*Cos(Pi*x)*Cos(Pi*y) - 2*t2n*(Cos(2*Pi*x)+Cos(2*Pi*y))

return DNosX(tmpEg, engy)

killwaves tmpEg
end
//============================================