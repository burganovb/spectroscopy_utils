#pragma rtGlobals=3		// Use modern global access method and strict wave access.


function mydisp_QtoK(disp,EF,theta,phi,omega,orientation,EFcorr)
Wave disp //1st column - angle Q, 2nd - kinetic energy
Variable EF//in eV
Variable EFcorr //correction to EF to correct dilation effects, e.g. if want to use 17.87 instead of 16.87 => EFcorr=1
	Variable theta		// sample manipulator theta angle (degrees)
	Variable phi			// sample manipulator phi angle (degrees)
	Variable omega		// sample manipulator omega angle (degrees)
	Variable orientation	// analyzer orientation (1 = vertical angle dispersion, 2 = horizontal angle dispersion)

	Variable KE			// electron kinetic energy (eV)
	Variable Q			// analyzer measurement angle (degrees)

Make/O/N=(DimSize(disp,0),3) $NameOfWave(disp)+"_EK"
Wave kdisp=$NameOfWave(disp)+"_EK"

kdisp[][0] = real(BlueZoneAngleToMomentum(disp[p][1]+EFcorr,theta,phi,omega,disp[p][0],orientation))
kdisp[][1] = imag(BlueZoneAngleToMomentum(disp[p][1]+EFcorr,theta,phi,omega,disp[p][0],orientation))
kdisp[][2] = -EF + disp[p][1]

End

//==============================================================
Function MySplitDisp(disp, splitE)
//splits a wave with dispersion (3-column wave: kx,ky, -BE) into two disp waves containing energies only below or above splitE
Wave disp
Variable splitE

duplicate/O disp, $NameOfWave(disp)+"_0", $NameOfWave(disp)+"_1"
Wave disp0=$NameOfWave(disp)+"_0"
Wave disp1=$NameOfWave(disp)+"_1"

variable k,k0=0,k1=0
For (k=0; k<DimSize(disp,0);k+=1)
	if (disp[k][2]>splitE)
		DeletePoints k0, 1, disp0
		k1+=1
	else
		DeletePoints k1, 1, disp1
		k0+=1	
	endif
endfor

End
//==============================================================

//==============================================================
Function myGetParam(parnum, parwv)
// extracts parameter from the fitting output. the output of this function is an XY-wave: energy-parameter for MDC fits or momentum-parameter for EDC
Wave parwv
Variable parnum
Variable k,m=0
Make/O/N=(DimSize(parwv,1),2) $NameOfWave(parwv)+num2str(parnum)
Wave wvout = $NameOfWave(parwv)+num2str(parnum)
wvout[][0] = DimOffset(parwv,1) + p*DimDelta(parwv,1)
wvout[][1] = parwv[parnum][p]
for (k=0;k<DimSize(parwv,1);k+=1)
	If (numtype(wvout[m][1])>0)
		DeletePoints m, 1, wvout
	Else
		m+=1
	Endif
Endfor
End
//==============================================================
Function my10toKE(wvin, EF, KE, theta, phi, omega, edc)
Wave wvin
Variable EF, KE, theta, phi, omega// analyzer orientation (1 = vertical angle dispersion, 2 = horizontal angle dispersion)
Variable edc // 0 for MDC , 1 for EDC
//EF- measured EF (e.g. 16.864 eV); KE - value to use for Q->k conversion for electrons at EF (e.g. 18.2eV or 16.864eV)
Make/O/N=(DimSize(wvin,0),3) $NameOfWave(wvin)+"_k"
Wave wvout = $NameOfWave(wvin)+"_k"

wvout[][2]=-EF + wvin[p][edc - 0]
wvout[][0] = real(BlueZoneAngleToMomentum(wvin[p][edc - 0]+KE-EF, theta,phi,omega,wvin[p][1 - edc],1))
wvout[][1] = imag(BlueZoneAngleToMomentum(wvin[p][edc - 0]+KE-EF, theta,phi,omega,wvin[p][1 - edc],1))
//wvout[][0] = real(BlueZoneAngleToMomentum(wvin[p][edc], theta,phi,omega,wvin[p][1 - edc],1))
//wvout[][1] = imag(BlueZoneAngleToMomentum(wvin[p][edc], theta,phi,omega,wvin[p][1 - edc],1))

End
//==============================================================
Function calcSE1(data, band)
// data is xy-wave : (k, BE) pairs from ARPES
// band is same format from DFT
Wave data, band
Duplicate/O data, $NameOfWave(data)+"SE1"
Wave wSE = $NameOfWave(data)+"SE1"
wSE[][0] = data[p][1]
make/o/n=(DimSize(band,0)) bandk, bandE
bandk=band[p][0]
bandE=band[p][1]
make/o/n = (abs(band[DimSize(band,0)-1][0] - band[0][0])/0.001+1) band1
SetScale /P x, min(band[0][0],band[DimSize(band,0)-1][0]), 0.001, band1
Interpolate2/I=3/Y=band1 bandk, bandE

wSE[][1] = data[p][1] - band1(data[p][0])

End
//==============================================================
Function EKtoEK2(spec)
Wave spec
Make/O/N=(DimSize(spec,1)) Kin, K2in, mdck1
Kin = DimOffset(spec,1)+p*DimDelta(spec,1)
K2in=Kin[p]^2
SetScale/P x, DimOffset(spec,1), DimDelta(spec,1), K2in
K2in[0,x2pnt(K2in,0)]*=-1
Make/O/N=(DimSize(spec,0), DimSize(spec,1)*5) $NameOfWave(spec)+"vK2"
Wave wout = $NameOfWave(spec)+"vK2"
SetScale/P x, DimOffset(spec,0), DimDelta(spec,0), wout
SetScale/I y, K2in[0], K2in[DimSize(K2in,0)-1], wout
Make/O/N=(DimSize(spec,0), DimSize(spec,1)*5) mdck2
SetScale/I x, K2in[0], K2in[DimSize(K2in,0)-1], mdck2

Variable n, x2, x1, m
For (n=0; n<DimSize(spec,0); n+=1)
	mdck1 = spec[n][p]
	//Interpolate2/I=3/T=1/Y=mdck2 K2in, mdck1
	For (m=0;m<DimSize(wout,1);m+=1)
		x2 = DimOffset(wout,1)+m*DimDelta(wout,1)
		x1 = sqrt(abs(x2))*sign(x2)
		wout[n][m] = interp(x1, Kin, mdck1)
	EndFor
EndFor

//KillWaves mdck1, mdck2, Kin, K2in
End
//==============================================================
Function EKspecFillZeros(spec)
// the output of BlueAngleToMOmentum sometimes has streaks of zeros. This function fills those areas with neighbour averages
Wave spec

Variable m,k
For (m=1; m<DimSize(spec,0)-1;m+=1)
	For (k=1;k<DimSize(spec,1)-1;k+=1)
		if (spec[m][k]==0)
			spec[m][k] = (spec[m-1][k]+spec[m+1][k]+spec[m][k-1]+spec[m][k+1])/(sign(spec[m-1][k]-1e-6)+sign(spec[m+1][k]-1e-6)+sign(spec[m][k-1]-1e-6)+sign(spec[m][k+1]-1e-6)+4+1e-8)*2
		endif
	endfor
endfor
End
//==============================================================
Function EKspecShiftE(spec, KE, EF) //spec in Blue format
Wave spec
Variable KE,EF
SetScale/P x, DimOffset(spec,0)+KE-EF, DimDelta(spec,0),spec
End
//==============================================================
Function FitRenormBandToData(data, band,n1,n2)
//finds best renorm factor and {x,y} shift for the ARPES data and LDA band
wave band,data
variable n1,n2

make/O/n=(n2-n1+1,2) datatmp = data[p+n1][q]
make/O/n=(Dimsize(band,0)) bandX,bandY
bandX = band[p][0]
bandY = band[p][1]
make/O/N=0 temp
//NVAR x0=0,y0=0,scaling=1
Optimize/X= {0,0,1} BandToDataDist, temp

//KillWaves bandX,bandY, datatmp, temp
KillWaves temp
End
//------------------------------------
Function BandToDataDist(w, x0, y0, a)
Wave w // not used
Variable x0,y0,a
Wave bandX,bandY, datatmp
Variable k,s=0

For (k=0;k<DimSize(datatmp,0);k+=1)
	s+= (datatmp[k][1] - a*(interp(datatmp[k][0]+x0,bandX,bandY) + y0))^2
EndFor

return s
End
//=============================================================
Function DispFromCP(Dx, Dy, x0,y0, x2y,xstep,iter)
// finds a trace on an image by gradient ascent

//DM1(spec, Dx, 2, 9, 2)  -  to prepare Dx and Dy
//DM1(spec, Dy, 2, 9, 1)
//DispFromCP(Dx, Dy, -0.260, -1.0279, 0.5, 0.001, 1000)    -  for a EK wave in BLUE format

Wave Dx, Dy
Variable x2y,xstep,iter,x0,y0
Make/o/n=(iter,2) dispout
Variable k, delX, delY
dispout[0][0] = x0
dispout[0][1] = y0

for (k=1;k<iter;k+=1)
	
	delX= x2y*xstep*Dx(x0)(y0)/sqrt( x2y^2*Dx(x0)(y0)^2 + Dy(x0)(y0)^2 )
	delY= xstep*Dy(x0)(y0)/sqrt( x2y^2*Dx(x0)(y0)^2 + Dy(x0)(y0)^2 )

	delX= xstep
	delY= delX/Dx(x0)(y0)*Dy(x0)(y0)/x2y
	
	x0+=delX
	y0+=delY
	dispout[k][0] = x0
	dispout[k][1] = y0
endfor

End
//=============================================================
Function MyXYZ2spline(XYZwv, nX, nY, x0,y0,dx,dy)
Wave XYZwv
Variable nX, nY, x0,y0,dx,dy

Make/o/n=(nX,nY) image
variable k
for (k=0;k<dimsize(XYZwv,0);k+=1)
	image[(XYZwv[k][0]-x0)/dx][(XYZwv[k][1]-y0)/dy] = XYZwv[k][2]
endfor
SetScale/P x, x0, dx, image
SetScale/P y, y0, dy, image

ImageInterpolate/S={0,0.005,1,0,0.005,1} Spline image
Wave M_InterpolatedImage

//Duplicate/O M_InterpolatedImage, $nameofwave(XYZwv)+"_spln"

Make/o/n=(2*dimsize(M_InterpolatedImage,0)+1, 2*dimsize(M_InterpolatedImage,1)+1) $nameofwave(XYZwv)+"_spln"
Wave wout=$nameofwave(XYZwv)+"_spln"
setscale/I x, -1 ,1 , wout
setscale/I y, -1 ,1 , wout
wout = M_InterpolatedImage(abs(x))(abs(y))

Killwaves image, M_InterpolatedImage
End
//=============================================================
Function RescaleXY(win, a,b)
Wave win
Variable a,b
SetScale/P x, a*DimOffset(win,0), a*DimDelta(win,0), win
SetScale/P y, b*DimOffset(win,1), b*DimDelta(win,1), win
End
//=============================================================
Function MyGetLDAcut(band2d, theta, phi0, omega, KE)
Wave band2d
Variable theta, phi0, omega, KE

Variable k, Q
Make/o/n=(1201,3) $nameofwave(band2d)+"cut"
Wave wout = $nameofwave(band2d)+"cut"
for (k=0; k<1201;k+=1)
	Q = -30 + k*0.05
	wout[k][0] = real(BlueZoneAngleToMomentum(KE,theta,phi0,omega,Q,1))
	wout[k][1] = imag(BlueZoneAngleToMomentum(KE,theta,phi0,omega,Q,1))
	wout[k][2] = Interp2D (band2d,  wout[k][0], wout[k][1])
endfor

End
//=============================================================
function DispFromDeriv(spec, p0, q0, numind, wind, isX) //doesnt work
wave spec
variable p0, q0, numind, wind, isX
make/o/n=(2*wind+1) tmpwind
variable k
if (isX==1)
	make/o/n=(numind, 3) $nameofwave(spec)+"_disp"
	wave dispout = $nameofwave(spec)+"_disp"
	for (k=0;k<numind;k+=1)
		tmpwind = spec[p0+k][q0-wind+p]
		CurveFit/Q/M=0/W=0 line, tmpwind
		wave W_coef
		wave W_sigma
		dispout[k][0] = dimoffset(spec,1)+(q0-wind-W_coef[0]/W_coef[1])*dimdelta(spec,1)
		dispout[k][1] = dimoffset(spec,0)+(p0+k)*dimdelta(spec,0)
		q0 = q0-wind-W_coef[0]/W_coef[1]
	endfor
else
	make/o/n=(numind, 3) $nameofwave(spec)+"_disp"
	wave dispout = $nameofwave(spec)+"_disp"
	for (k=0;k<numind;k+=1)
		tmpwind = spec[q0+k][p0-wind+p]
		CurveFit/Q/M=0/W=0 line, tmpwind
		wave W_coef
		wave W_sigma
		dispout[k][1] = dimoffset(spec,0)+(p0-wind-W_coef[0]/W_coef[1])*dimdelta(spec,0)
		dispout[k][0] = dimoffset(spec,1)+(q0+k)*dimdelta(spec,1)
		p0 = p0-wind-W_coef[0]/W_coef[1]
	endfor
endif
killwaves tmpwind
end
//=============================================================
function DispFromPeak(spec, p0, q0, numind, wind, box, isX)
//DispFromPeak(bro8gm_iek, 233, 119, 210, 10, 1) //box = 11
wave spec
variable p0, q0, numind, wind, isX,box
make/o/n=(2*wind+1) tmpwind
variable k, q1=q0, p1=p0
if (isX==1)
	make/o/n=(numind, 4) $nameofwave(spec)+"_dispX"
	wave dispout = $nameofwave(spec)+"_dispX"
	for (k=0;k<numind;k+=1)
		tmpwind = spec[p0+k][q1-wind+p]
		FindPeak/B=(box)/P/Q tmpwind
		dispout[k][1] = dimoffset(spec,0)+(p0+k)*dimdelta(spec,0)
		dispout[k][0] = dimoffset(spec,1)+(q1-wind + V_PeakLoc)*dimdelta(spec,1)
		dispout[k][3] = V_PeakVal
		dispout[k][2] = V_flag
		if (V_flag==0)
			q1 += -wind + V_PeakLoc
		endif
	endfor
else
	make/o/n=(numind, 4) $nameofwave(spec)+"_dispY"
	wave dispout = $nameofwave(spec)+"_dispY"
	for (k=0;k<numind;k+=1)
		tmpwind = spec[p1-wind+p][q0+k]
		FindPeak/B=10/P/Q tmpwind
		dispout[k][0] = dimoffset(spec,1)+(q0+k)*dimdelta(spec,1)
		dispout[k][1] = dimoffset(spec,0)+(p1-wind + V_PeakLoc)*dimdelta(spec,0)
		dispout[k][3] = V_PeakVal
		dispout[k][2] = V_flag
		if (V_flag==0)
			p1 += -wind + V_PeakLoc
		endif
	endfor
endif
//killwaves tmpwind
end
//=============================================================
function myswapXY(waveXY)
wave waveXY
make/o/n=(dimsize(waveXY,0)) tmpmyswapXY
tmpmyswapXY = waveXY[p][0]
waveXY[][0] = waveXY[p][1]
waveXY[][1] = tmpmyswapXY[p]
killwaves tmpmyswapXY
end
//=============================================================
function bandfit2band(win)
wave win //format: 1D wave with Engy values, x-dimension==k values
//outputs XY wave suitable for calcSE1
make/n=(dimsize(win,0),2) $nameofwave(win)+"band"
wave wout=$nameofwave(win)+"band"
wout[][0] = dimoffset(win,0) + p*dimdelta(win,0)
wout[][1] = win[p]
End

//==========================