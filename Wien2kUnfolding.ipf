#pragma rtGlobals=3		// Use modern global access method and strict wave access.

///////////////////////////////////////////
// depends on HyperCube.ipf to run properly
///////////////////////////////////////////


function loadunfolded(EF, BandMin, BandMax, latticeType, N1, N2, N3, caxis, primitiveBZ, superx, supery, superz)
	variable EF, BandMin, BandMax, latticeType, N1, N2, N3, caxis, primitiveBZ //same as inputs for loadhypercube
	variable superx, supery, superz // number of u.c. repeats in supercell

	wave kx,ky,kz,engy, w // these need to be loaded (update this part)---------------------------

//making hypercubes consistent with the WienHyperCube
	If(primitiveBZ < 1)
		print "primitiveBZ was <1, I'll assume you mean primitiveBZ = 1"
		primitiveBZ = 1
	EndIf

	If(primitiveBZ == 2)
		If(caxis == 1)
			N3 = N2
		Elseif(caxis ==2)
			N3 = N1
		Else
			N2 = N1
		EndIf
	EndIf
	// Convert N1, N2, N3 to values consistent with GenerateHCkmesh's treatment (nearest, but smaller, odd number)
	N1 = 2*floor((N1+1) /2) - 1
	N2 = 2*floor((N2+1) /2) - 1
	N3 = 2*floor((N3+1) /2) - 1
	Variable BandNum = BandMax - BandMin + 1

	// Make 4D waves
	Make /O /N = ((N1-1)*superx+1, (N2-1)*supery+1, (N3-1)*superz+1, BandNum) UnfoldedHCenergy, UnfoldedHCweight
	
	// Set the scale for the hypercubes (-1 to 1 in units of pi/a)
	SetScale /I x, -1, 1, UnfoldedHCenergy, UnfoldedHCweight
	SetScale /I y, -1, 1, UnfoldedHCenergy, UnfoldedHCweight
	SetScale /I z, -1, 1, UnfoldedHCenergy, UnfoldedHCweight
	SetScale /I t, BandMin, BandMax, UnfoldedHCenergy, UnfoldedHCweight
	UnfoldedHCenergy = NaN
	UnfoldedHCweight = NaN
	
	Sort {kx,ky,kz,engy} kx,ky,kz,engy,w
	
	variable bandind
	variable dkx=dimdelta(UnfoldedHCenergy,0)
	variable dky=dimdelta(UnfoldedHCenergy,1)
	variable dkz=dimdelta(UnfoldedHCenergy,2)
	variable kxcurr,kxprev=inf
	variable kycurr,kyprev=inf
	variable kzcurr,kzprev=inf
	variable k
	for (k=0;k<dimsize(kx,0);k+=1)
		kxcurr=kx[k]
		kycurr=ky[k]
		kzcurr=kz[k]
		if (kxcurr==kxprev && kycurr==kyprev && kzcurr==kzprev)
			bandind+=1
		else
			bandind=1
		endif
		if (bandind>=BandMin && bandind<=BandMax)
		//UnfoldedHCenergy(kxcurr)(kycurr)(kzcurr)(bandind
			HC_AssignWithSymm(UnfoldedHCenergy, 2*kxcurr/dkx,2*kycurr/dky,2*abs(kzcurr/dkz), bandind-BandMin, engy[k], latticetype, caxis, primitiveBZ)
			HC_AssignWithSymm(UnfoldedHCweight, 2*kxcurr/dkx,2*kycurr/dky,2*abs(kzcurr/dkz), bandind-BandMin, w[k], latticetype, caxis, primitiveBZ)
		endif
		kxprev=kxcurr;kyprev=kycurr;kzprev=kzcurr
	endfor

	UnfoldedHCenergy -= EF
	UnfoldedHCenergy *= 13.605698
//HC_AssignWithSymm(currentcase, WBC_klist[kindex][0], WBC_klist[kindex][1], WBC_klist[kindex][2], bandindex, currval, latticetype, caxis, primitiveBZ)

end
//=============================
function getunfoldedbandchar(bandcharcube, superx, supery, superz)
wave bandcharcube //bandcharacter cube for folded hypercube
variable superx, supery, superz //# of repeats in a supercell

make/n=((dimsize(bandcharcube,0)-1)*superx+1, (dimsize(bandcharcube,1)-1)*supery+1, (dimsize(bandcharcube,2)-1)*superz+1, dimsize(bandcharcube,3)) $"unfolded_"+nameofwave(bandcharcube)
wave unfbandchar=$"unfolded_"+nameofwave(bandcharcube)
setscale/P x, dimoffset(bandcharcube,0), dimdelta(bandcharcube,0)/superx, unfbandchar
setscale/P y, dimoffset(bandcharcube,1), dimdelta(bandcharcube,1)/supery, unfbandchar
setscale/P z, dimoffset(bandcharcube,2), dimdelta(bandcharcube,2)/superz, unfbandchar

unfbandchar = bandcharcube[mod(p, dimsize(bandcharcube,0))][mod(q, dimsize(bandcharcube,1))][mod(r, dimsize(bandcharcube,2))][s]

end
//=============================