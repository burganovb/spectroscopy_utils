#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//============================================
//============================================
//============================================
//============================================
function mySum2dDim(wvin, dim, doabs)
wave wvin
variable doabs, dim //dim=0 sum rows
variable kk,mm, N, S=0
if (dim==0)
	N = DimSize(wvin,1)
	make/o/n=(N) sumout=0
	setscale/P x, dimoffset(wvin,1), dimdelta(wvin,1), sumout
	for (kk=0;kk<dimsize(wvin,0);kk+=1)
		for (mm=0;mm<dimsize(wvin,1);mm+=1)
			if (doabs==1)
				sumout[mm]+=abs(wvin[kk][mm])
			else
				sumout[mm]+=wvin[kk][mm]
			endif
		endfor
	endfor
else
	N = DimSize(wvin,0)
	make/o/n=(N) sumout=0
	setscale/P x, dimoffset(wvin,0), dimdelta(wvin,0), sumout
	for (kk=0;kk<dimsize(wvin,0);kk+=1)
		for (mm=0;mm<dimsize(wvin,1);mm+=1)
			if (doabs==1)
				sumout[kk]+=abs(wvin[kk][mm])
			else
				sumout[kk]+=wvin[kk][mm]
			endif
		endfor
	endfor

endif
end
//============================================
function mySum2d(wvin, doabs)
wave wvin
variable doabs
variable kk,mm, S=0
for (kk=0;kk<dimsize(wvin,0);kk+=1)
	for (mm=0;mm<dimsize(wvin,1);mm+=1)
		if (doabs==1)
			S+=abs(wvin[kk][mm])
		else
			S+=wvin[kk][mm]
		endif
	endfor
endfor
return S
end
//============================================
function myeffinteract(kfpts, vfvals, gapvalues, chiW, eitype, massren, U)
wave kfpts, vfvals, gapvalues, chiW 
// vfvals - Abs() values of fermi velocity at kfpts (N-by-1); kfpts - kFermi points (N-by-3, 3rd clmn for dk); 
//gapvalues - gapvalues at kfpts (N-by-1 or 2 for cmplx); chiW - suscept (M-by-M)
variable eitype // 0-sum of values; 1 - div by vf^.5; 2 - div by vf; 
variable massren // 0 - no renorm; 1 - div by renpar; 2 - div by (1+renpar). Renpar calculated as effinteraction but without gap function
variable U // for Stoner enhancement

variable N= dimsize(kfpts,0)
variable kk,mm,Qx,Qy,chiV, eisum=0, rensum = 0, gapsum=0, OneFS
variable deisum=0, dgapsum=0, drensum, dOneFS

for (kk=0;kk<N;kk+=1)
	for (mm=0;mm<N;mm+=1)
//		Qx = mod(5 + kfpts[kk][0] - kfpts[mm][0],1)
//		Qy = mod(5 + kfpts[kk][1] - kfpts[mm][1],1)
		Qx = mod(5 + kfpts[kk][0] - kfpts[mm][0],2)-1
		Qy = mod(5 + kfpts[kk][1] - kfpts[mm][1],2)-1
		chiV = chiW(Qx)(Qy)/(1-U*chiW(Qx)(Qy))
		deisum = gapvalues[kk][0]*kfpts[kk][2]*gapvalues[mm][0]*kfpts[mm][2]*chiV
		//exclude umklapp:
//		deisum *= (sign(gapvalues[kk][0])*sign(gapvalues[mm][0])+1)/2
		if (massren>0)
			drensum = kfpts[kk][2]*kfpts[mm][2]*chiV
		endif
		if (eitype ==1)
			deisum/=(vfvals[kk]*vfvals[mm])^.5
			drensum/=(vfvals[kk]*vfvals[mm])^.5
		elseif (eitype ==2)
			deisum/=(vfvals[kk]*vfvals[mm])
			drensum/=(vfvals[kk]*vfvals[mm])
		endif
		eisum += deisum
		rensum += drensum
	endfor
	if (eitype ==0)
		gapsum += gapvalues[kk][0]^2*kfpts[kk][2]
		OneFS += kfpts[kk][2]
	elseif (eitype ==1)
		gapsum += gapvalues[kk][0]^2*kfpts[kk][2]/vfvals[kk]^.5
		OneFS += kfpts[kk][2]/vfvals[kk]^.5
	elseif (eitype ==2)
		gapsum += gapvalues[kk][0]^2*kfpts[kk][2]/vfvals[kk]
		OneFS += kfpts[kk][2]/vfvals[kk]
	endif	
endfor
//print gapint1
//print gapint2
//print sumval1
//print sumval2
//print sqrt(sumval1^2+sumval2^2)
//return sumval1/gapint1
if (massren==0)
	return eisum/gapsum
elseif (massren==1)
	return eisum/gapsum/drensum*OneFS
else
	return eisum/gapsum/(drensum/OneFS + 1)
endif

end
//============================================
function mycalcgapfunc(kfpts, type)
wave kfpts
variable type //0 - sinX+i*sinY; 1 - 

make/o/n=(dimsize(kfpts,0),2) gapvalues
if (type==0)
	gapvalues = sin(Pi*kfpts[p][q])
elseif (type==801)
	gapvalues = sin(Pi*cos(atan2(kfpts[p][1],kfpts[p][0])))
elseif (type==901)
	gapvalues = sin(2*Pi*kfpts[p][q])
elseif (type==902)
	gapvalues = sin(3*Pi*kfpts[p][q])
elseif (type==1) //Wang, 2013
	gapvalues = sin(Pi*kfpts[p][q])*(-0.4375+cos(Pi*kfpts[p][1-q]))
elseif (type==2) //Raghu, 2010 for XZ
	gapvalues = sin(Pi*kfpts[p][q])*cos(Pi*kfpts[p][1-q])
elseif (type==21)
	gapvalues = sin(2*Pi*kfpts[p][q])*cos(Pi*kfpts[p][1-q])
elseif (type==22)
	gapvalues = sin(3*Pi*kfpts[p][q])*cos(Pi*kfpts[p][1-q])
elseif (type==3) //d-wave x2-y2
	gapvalues = cos(Pi*kfpts[p][0]) - cos(Pi*kfpts[p][1])
elseif (type==4) //d-wave xy
	gapvalues = sin(Pi*kfpts[p][0])*sin(Pi*kfpts[p][1])
elseif (type==5) //d-wave xy
	gapvalues = 1
endif

end
//============================================
function MyFSpointsG(N, t3,t4,mu1d,mu2d) // generates kF points for gamma band and fermi velocities
variable N, t3,t4,mu1d,mu2d

variable kk, mm, kfy, kfx, kfx0, kfx1, cval
make/o/n=(N) kfXpts=0, kfYpts=0, kfDpts=0
//make/o/n=(N,2) kfspts
variable band = 2
if (mu2d-4*t4==0)
	mu2d -= 0.00001
endif
if (mu2d-4*t4<0)
	kfYpts=p/(N-1)
else
	kfYpts=1-p/(N-1)
endif

for (kk=0; kk<N; kk+=1)
	kfy = kfYpts[kk]
	kfx0 = 0; kfx1 = 1
	for (mm=0; mm<20; mm+=1)
		kfx = 0.5*(kfx0+kfx1)
		 if (band==2)
		 	cval = -2*(Cos(Pi*kfx)+Cos(Pi*kfy)) - 4*t4*Cos(Pi*kfx)*Cos(Pi*kfy) - mu2d
		 else
		 	cval = - mu1d - Cos(Pi*kfx) - t3*Cos(Pi*kfy)
		 endif
		 if (cval<0)
		 	kfx0 = kfx
		 else
		 	kfx1 = kfx
		 endif
	endfor
	if ((kfy-kfx)*(mu2d-4*t4)<0)
		DeletePoints kk, N-kk, kfYpts, kfXpts
		break
	endif
	kfXpts[kk] = kfx
endfor
InsertPoints kk, kk, kfYpts, kfXpts
kfXpts[kk, DimSize(kfXpts,0)-1] = kfYpts[DimSize(kfXpts,0)-1 - p]
kfYpts[kk, DimSize(kfXpts,0)-1] = kfXpts[DimSize(kfXpts,0)-1 - p]

N = DimSize(kfXpts,0)-1

Make/o/n=(4*N, 3) kfpts
kfpts[0,N-1][0] = 0.5*(kfXpts[p] + kfXpts[p+1])
kfpts[0,N-1][1] = 0.5*(kfYpts[p] + kfYpts[p+1])
kfpts[0,N-1][2] = sqrt( (kfXpts[p] - kfXpts[p+1])^2 + (kfYpts[p] - kfYpts[p+1])^2 )

kfpts[N,2*N-1][0] = kfpts[p-N][q]
kfpts[N,2*N-1][1] = -kfpts[p-N][q]
kfpts[N,2*N-1][2] = kfpts[p-N][q]

kfpts[2*N,3*N-1][0] = -kfpts[p-2*N][q]
kfpts[2*N,3*N-1][1] = -kfpts[p-2*N][q]
kfpts[2*N,3*N-1][2] = kfpts[p-2*N][q]

kfpts[3*N,4*N-1][0] = -kfpts[p-3*N][q]
kfpts[3*N,4*N-1][1] = kfpts[p-3*N][q]
kfpts[3*N,4*N-1][2] = kfpts[p-3*N][q]

make/o/n=(dimsize(kfpts,0)) vfvals
vfvals = 2*pi*sqrt( (sin(pi*kfpts[p][0])*(1+2*t4*cos(pi*kfpts[p][1])))^2 + (sin(pi*kfpts[p][1])*(1+2*t4*cos(pi*kfpts[p][0])))^2 )
variable meanvf=mean(vfvals)
vfvals/=meanvf

End
//============================================
function MyFSpointsG_f(Eg, mu2d, t4, T, thresh) // generates kF points for gamma band and fermi velocities
Wave Eg
variable T, thresh, mu2d, t4

Eg = -2*(Cos(Pi*x)+Cos(Pi*y)) - 4*t4*Cos(Pi*x)*Cos(Pi*y) - mu2d
duplicate/o Eg, tmpKfsSQ
tmpKfsSQ[][][0] = x
tmpKfsSQ[][][1] = y
variable N = dimsize(Eg,0)
make/o/n=(N*N) tmpKfsX, tmpKfsY, tmpKfsD, tmpKfsSGN
tmpKfsX = tmpKfsSQ[mod(p,N)][floor(p/N)][0]
tmpKfsY = tmpKfsSQ[mod(p,N)][floor(p/N)][1]
tmpKfsD = (1+exp(Eg[mod(p,N)][floor(p/N)][0]/T))^-1*(1-(1+exp(Eg[mod(p,N)][floor(p/N)][0]/T))^-1)

Sort/R tmpKfsD, tmpKfsD, tmpKfsX, tmpKfsY
tmpKfsSGN = (sign(tmpKfsD[p]-thresh)+1)*0.5
variable M = sum(tmpKfsSGN)
print M
DeletePoints M, N*N-M, tmpKfsX, tmpKfsY, tmpKfsD

make/o/n=(M,3) kfptsF
kfptsF[][0] = tmpKfsX[p]
kfptsF[][1] = tmpKfsY[p]
kfptsF[][2] = tmpKfsD[p]

//duplicate/o tmpDKfs, tmpsgnDKfs = sign(tmpDKfs - thresh)
//variable sumDkf = mySum2d(tmpDKfs, 0)
//variable sumDkf = mySum2d(tmpDKfs, 0)
//tmpDKfs/=sumDkf

/////////////////////////////////////
make/o/n=(dimsize(kfptsF,0)) vfvals
vfvals = 2*pi*sqrt( (sin(pi*kfptsF[p][0])*(1+2*t4*cos(pi*kfptsF[p][1])))^2 + (sin(pi*kfptsF[p][1])*(1+2*t4*cos(pi*kfptsF[p][0])))^2 )
variable meanvf=mean(vfvals)
vfvals/=meanvf
//
KillWaves tmpKfsSQ, tmpKfsSGN, tmpKfsX, tmpKfsY, tmpKfsD
End
//============================================
function MyFSpointsG1(N8, mu2d, t4) // generates kF points for gamma band and fermi velocities
variable N8, mu2d, t4

variable N2 = N8
variable Nmat = 2*N2
make/O/N=(8*N2, 3) kfpts
make/O/N=(8*N2) vfvals, DK

variable Y0 = acos(min(1, (2-mu2d)/(2-4*t4)))/pi
variable Y1 = acos((-1+sqrt(1-mu2d*t4))/t4/2)/pi
kfpts[0,N2-1][1] = Y0+(Y1-Y0)/N2*(0.5 + p)
kfpts[0,N2-1][0] = acos(-(mu2d+2*cos(pi*kfpts[p][1]))/(2+4*t4*cos(pi*kfpts[p][1])))/pi
kfpts[N2,2*N2-1][1] = kfpts[2*N2-p-1][0]
kfpts[N2,2*N2-1][0] = kfpts[2*N2-p-1][1]
kfpts[Nmat,2*Nmat-1][0] = -kfpts[p-Nmat][1]
kfpts[Nmat,2*Nmat-1][1] = kfpts[p-Nmat][0]
kfpts[2*Nmat,3*Nmat-1][0] = -kfpts[p-2*Nmat][0]
kfpts[2*Nmat,3*Nmat-1][1] = -kfpts[p-2*Nmat][1]
kfpts[3*Nmat,4*Nmat-1][0] = kfpts[p-3*Nmat][1]
kfpts[3*Nmat,4*Nmat-1][1] = -kfpts[p-3*Nmat][0]

DK[0, 8*N8-2] =  sqrt((kfpts[p][0]-kfpts[p+1][0])^2+(kfpts[p][1]-kfpts[p+1][1])^2)
DK[8*N8-1] = DK[0]
variable DKs = sum(DK)
kfpts[][2] = DK[p]/DKs

vfvals = 2*pi*sqrt( (sin(pi*kfpts[p][0])*(1+2*t4*cos(pi*kfpts[p][1])))^2 + (sin(pi*kfpts[p][1])*(1+2*t4*cos(pi*kfpts[p][0])))^2 )
variable meanvf=mean(vfvals)
vfvals/=meanvf

KillWaves DK
End
//============================================
function MyFSpointsXZ(N, t3,t4,mu1d,mu2d)
variable N, t3,t4,mu1d,mu2d

variable kk, mm, kfy, kfx, kfx0, kfx1, cval
make/o/n=(N) kfXpts=0, kfYpts=0
//make/o/n=(N,2) kfspts
kfYpts=p/(N-1)

variable band = 2
for (kk=0; kk<N; kk+=1)
	kfy = kfYpts[kk]
	kfx0 = 0; kfx1 = 1
	for (mm=0; mm<20; mm+=1)
		kfx = 0.5*(kfx0+kfx1)
		 	cval = - mu1d - 2*Cos(Pi*kfx) - 2*t3*Cos(Pi*kfy)
		 if (cval<0)
		 	kfx0 = kfx
		 else
		 	kfx1 = kfx
		 endif
	endfor
	kfXpts[kk] = kfx
endfor

N = DimSize(kfXpts,0)-1

Make/o/n=(4*N, 3) kfpts
kfpts[0,N-1][0] = 0.5*(kfXpts[p] + kfXpts[p+1])
kfpts[0,N-1][1] = 0.5*(kfYpts[p] + kfYpts[p+1])
kfpts[0,N-1][2] = sqrt( (kfXpts[p] - kfXpts[p+1])^2 + (kfYpts[p] - kfYpts[p+1])^2 )

kfpts[N,2*N-1][0] = kfpts[p-N][q]
kfpts[N,2*N-1][1] = -kfpts[p-N][q]
kfpts[N,2*N-1][2] = kfpts[p-N][q]

kfpts[2*N,3*N-1][0] = -kfpts[p-2*N][q]
kfpts[2*N,3*N-1][1] = -kfpts[p-2*N][q]
kfpts[2*N,3*N-1][2] = kfpts[p-2*N][q]

kfpts[3*N,4*N-1][0] = -kfpts[p-3*N][q]
kfpts[3*N,4*N-1][1] = kfpts[p-3*N][q]
kfpts[3*N,4*N-1][2] = kfpts[p-3*N][q]
End
//============================================
function myEffInterListXY(chiBZ, mu2d, gaptype, eitype, massren, Ueff)
wave chiBZ, mu2d
variable gaptype, Ueff, eitype, massren
make/o/n=(dimsize(chiBZ,2)) EffInterOut

make/o/n=(dimsize(chiBZ,0),dimsize(chiBZ,1)) chiCurrBZ
setscale/I x, -1,1,chiCurrBZ
setscale/I y, -1,1,chiCurrBZ
variable kk
/////////
for (kk=0; kk<dimsize(chiBZ,2);kk+=1)
	MyFSpointsG(150, 0.1, 0.405, 1.08, mu2d[kk]);
	wave kfpts
	wave vfvals
	mycalcgapfunc(kfpts, gaptype);
	chiCurrBZ=chiBZ[p][q][kk]
	wave gapvalues
	EffInterOut[kk]=myeffinteract(kfpts, vfvals, gapvalues, chiCurrBZ, eitype, massren, Ueff);
endfor
return 0
///////////////////////////
wave Eg = $"tmp_Eg"
for (kk=0; kk<dimsize(chiBZ,2);kk+=1)
	MyFSpointsG_f(Eg, mu2d[kk], 0.405, 0.002, 0.01)
	wave kfptsF
	wave vfvals
	mycalcgapfunc(kfptsF, gaptype);
	chiCurrBZ=chiBZ[p][q][kk]
	wave gapvalues
	EffInterOut[kk]=myeffinteract(kfptsF, vfvals, gapvalues, chiCurrBZ, eitype, massren, Ueff);
endfor
end
//============================================
function myEffInterListXZ(chiBZ, mu1d, t3, gaptype, eitype, massren, Ueff)
wave chiBZ, mu1d, t3
variable gaptype, Ueff, eitype, massren
make/o/n=(dimsize(chiBZ,2)) EffInterOut

make/o/n=(dimsize(chiBZ,0),dimsize(chiBZ,1)) chiCurrBZ
setscale/I x, -1,1,chiCurrBZ
setscale/I y, -1,1,chiCurrBZ
variable kk
for (kk=0; kk<dimsize(chiBZ,2);kk+=1)
	MyFSpointsXZ(150, 0.1, 0.405, mu1d[kk], 1.48)
	wave kfpts
	wave vfvals
	mycalcgapfunc(kfpts, gaptype);
	chiCurrBZ=chiBZ(abs(x))(abs(y))[kk]
	wave gapvalues
	EffInterOut[kk]=myeffinteract(kfpts, vfvals, gapvalues, chiCurrBZ, eitype, massren, Ueff);
endfor
end
//============================================
function MyMBSusceptibility(SE1, SE2, band1, band2, Qx, Qy, Engy)
wave SE1, SE2, band1, band2
variable Qx, Qy, Engy



End
//============================================
function MyTSusceptibility(T, band1, band2, Qx, Qy, eps)
wave band1, band2
variable T, Qx, Qy, eps //eps>~0 small number; T is in same units as band1/band2
//Qx, Qy, and axis of band1, band2, are in pi/a units; BZ is [-1 to 1]
//band1 and  band2 must have same axis

make/o/n=(dimsize(band2,0),dimsize(band1,0)) chivals, band2vals, band2valsQ
setscale/P x, dimoffset(band2,0),dimdelta(band2,0), chivals, band2vals, band2valsQ
setscale/P y, dimoffset(band2,1),dimdelta(band2,1), chivals, band2vals, band2valsQ

band2vals = band2[p][q][0]
band2valsQ = Interp2d(band2vals, mod(x+Qx + 5,2)-1, mod(y+Qy+5,2)-1)
chivals = (Abs( (1+exp(band1[p][q][0]/T))^-1 - (1+exp(band2valsQ[p][q]/T))^-1 )+eps*T^-1*(exp(0.5*band1[p][q][0]/T)+exp(-0.5*band1[p][q][0]/T))^-2 )/(Abs(band1[p][q][0] - band2valsQ[p][q]) + eps)
//chivals = (Abs( (1+exp(band1[p][q][0]/T))^-1 - (1+exp(band2valsQ[p][q]/T))^-1 )+eps )/(Abs(band1[p][q][0] - band2valsQ[p][q]) + eps)

return mySum2d(chivals, 0)*DimDelta(band1,0)*DimDelta(band1,1)
//Interp2D (srcWaveName,  xValue, yValue )

//KillWaves band2vals, chivals
End

//============================================
function MyListGSusc(T, mu2, t4, t1, Eg, Qwedge, eps)
wave mu2, t4, t1, Eg, Qwedge
variable T, eps //1e-6

make/o/n=(dimsize(Qwedge,0), dimsize(mu2,0)) ChiWedgeList
variable kk
variable t0=ticks
for (kk=0;kk<dimsize(mu2,0);kk+=1)
	Eg = (-2*(Cos(Pi*x)+Cos(Pi*y)) - 4*t4[kk]*Cos(Pi*x)*Cos(Pi*y) - mu2[kk])
	ChiWedgeList[][kk] = MyTSusceptibility(T, Eg, Eg, Qwedge[p][0], Qwedge[p][1], eps)
	print kk
	print ticks-t0
endfor
End
//============================================
function MyWedge2Quad3D(qw, chiW, rind)
wave qw, chiW
variable rind

make/o/n=(dimsize(chiW,0),dimsize(chiW,1)) chiWtmp
chiWtmp = chiW[p][q][rind]
myWedge2Quad(qW, chiWtmp)

killwaves chiWtmp
end
//=================================================
function MyInteractMatrix(mu2d, t4, chiBZ, kpts)
variable mu2d, t4, kpts
wave chiBZ

MyFSpointsG1(kpts, mu2d, t4)
wave kfpts
wave vfvals
variable Nmat = DimSize(kfpts, 0)
make/o/n=(Nmat, Nmat) InteractMatrixOut
//InteractMatrixOut = Interp2d(chiBZ, mod(abs(kfpts[p][0]-kfpts[q][0]),1), mod(abs(kfpts[p][1]-kfpts[q][1]),1))*kfpts[p][2]*kfpts[q][2]/(vfvals[p]*vfvals[q])
InteractMatrixOut = Interp2d(chiBZ, mod(abs(kfpts[p][0]-kfpts[q][0]),1), mod(abs(kfpts[p][1]-kfpts[q][1]),1))/(vfvals[p]*vfvals[q])/kpts/8
End
//=================================================
function MyEigenVs(mu2d, t4, chiBZ, kpts8)
variable kpts8
wave chiBZ, mu2d, t4

make/o/n=(DimSize(mu2d,0),3) OutEVal //lowest 3 eigenvals
make/o/n=(DimSize(mu2d,0),8*kpts8,3) OutEVec
wave chiBZcurr
variable kk, Npts
for (kk=0; kk<DimSize(mu2d,0);kk+=1)
	chiBZcurr = chiBZ[p][q][kk]
	MyInteractMatrix(mu2d[kk], t4[kk], chiBZcurr, kpts8);
	wave InteractMatrixOut
	MatrixEigenV/SYM/EVEC/RNG={2,1,3} InteractMatrixOut
	wave M_eigenVectors
	wave W_eigenValues
	OutEVal[kk][] = W_eigenValues[q]
	OutEVec[kk][][] = M_eigenVectors[q][r]
endfor
End
//=================================================
function MyEffIntTestF(mu2d, t4, chiBZ, kpts8)
variable kpts8
wave chiBZ, mu2d, t4

make/o/n=(DimSize(mu2d,0)) OutEffInt //lowest 3 eigenvals
wave chiBZcurr
variable kk, Npts, testsum
make/o/n=(kpts8*8) tmptestfunc, tmptestfunc2, prodresult
for (kk=0; kk<DimSize(mu2d,0);kk+=1)
	chiBZcurr = chiBZ[p][q][kk]
	MyInteractMatrix(mu2d[kk], t4[kk], chiBZcurr, kpts8);
	wave InteractMatrixOut
	wave kfpts
	tmptestfunc = sin(pi*kfpts[p][0])*cos(pi*kfpts[p][1])
	tmptestfunc2 = tmptestfunc^2
	testsum = sqrt(sum(tmptestfunc2))
	tmptestfunc = tmptestfunc[p]/testsum
	MatrixMultiply InteractMatrixOut, tmptestfunc
	wave M_product
	prodresult = M_product[p]
	MatrixMultiply prodresult/T, tmptestfunc
	wave M_Product
	OutEffInt[kk] = M_Product[0][0]
endfor
KillWaves tmptestfunc, tmptestfunc2
End
//=================================================












