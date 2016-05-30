#pragma rtGlobals=3		// Use modern global access method and strict wave access.


////////////////////////////////////////////
// 0. LoadHyperCube  - not necessary
// 1. run WBCLoadKList with the same parameters as LoadHyperCube, open the same .energy file
// 2(a). Look into .qtl file and choose which atoms and orbitals you want.  orbital=0 is "tot"=total per atom"; atom=numatoms()+1 is interstitial; only "tot" is available for interstitial
// 2(b). make/n=(xx, 2) WBC_atomsorbitals contains pairs of (atoms, orbital) in each row
// 3. WBCGetBandCharacters
// ---------------------------------------------
// uses HyperCube.ipf : HC_AssignWithSymm
////////////////////////////////////////////

function WBCLoadKList(EF, BandMin, BandMax, latticeType, N1, N2, N3, caxis, primitiveBZ)
	Variable EF					// Fermi energy (Ry); output will be converted to eV
	Variable BandMin, BandMax 	// The indicies of the lowest and highest band the user wants to load (lowest band index is 1). Check "case.energy" for band numbers.
	Variable latticeType, N1, N2, N3, caxis, primitiveBZ // Parameters from GenerateHCkmesh()
	
	Variable fileRef
	String line
	Variable i
	Variable BandNum = BandMax - BandMin + 1
	Variable kx, ky, kz, bandindex
	Variable cubehalf = 0
	Variable G1, G2, G3

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

	// Make 4D wave
	Make /O /N = (N1, N2, N3, BandNum) $"WBC_cube0"
	Wave WienHyperCube = $"WBC_cube0"
	
	// Make klist wave compatible with HyperCube
	Make /O /N = (0, 3) WBC_klist
	
	// Set the scale for the hypercube (-1 to 1 in units of pi/a)
	SetScale /I x, -1, 1, WienHyperCube
	SetScale /I y, -1, 1, WienHyperCube
	SetScale /I z, -1, 1, WienHyperCube
	SetScale /I t, BandMin, BandMax, WienHyperCube
	
	WienHyperCube = 0
		
	// Rescale for primitiveBZ = 3
	If(primitiveBZ == 3)
		If(caxis == 1)
			SetScale /I y, -2, 2, WienHyperCube
			SetScale /I z, -2, 2, WienHyperCube
		elseIf(caxis == 2)
			SetScale /I x, -2, 2, WienHyperCube
			SetScale /I z, -2, 2, WienHyperCube
		Else
			SetScale /I x, -2, 2, WienHyperCube
			SetScale /I y, -2, 2, WienHyperCube
		EndIf
	EndIf

	Open /T="????" /R fileRef
	FReadLine fileRef, line

	variable wbc_klist_index=0
	do
		
		// If this line is a k-point;
		if(cmpstr(line[2,2],".") == 0)
		
			// set kx, ky, and kz
			If(primitiveBZ == 2)
				If(caxis == 1)
					kx = Round(2*str2num(line[1,18])/DimDelta(WienHyperCube,0))
					ky =  Round((str2num(line[20,37])+str2num(line[39,56]))/DimDelta(WienHyperCube,1))
					kz = Round((str2num(line[39,56])-str2num(line[20,37]))/DimDelta(WienHyperCube,2))
				ElseIf(caxis == 2)
					kx = Round((str2num(line[1,18])+str2num(line[39,56]))/DimDelta(WienHyperCube,0))
					ky = Round(2*str2num(line[20,37])/DimDelta(WienHyperCube,1))
					kz = Round((str2num(line[39,56])-str2num(line[1,18]))/DimDelta(WienHyperCube,2))
				Else
					kx = Round((str2num(line[1,18])+str2num(line[20,37]))/DimDelta(WienHyperCube,0))
					ky =  Round((str2num(line[20,37])-str2num(line[1,18]))/DimDelta(WienHyperCube,1))
					kz = Round(2*str2num(line[39,56])/DimDelta(WienHyperCube,2))
				EndIf
			Else
				kx = Round(2*str2num(line[1,18])/DimDelta(WienHyperCube,0))
				ky = Round(2*str2num(line[20,37])/DimDelta(WienHyperCube,1))
				kz = Round(2*str2num(line[39,56])/DimDelta(WienHyperCube,2))
			EndIf
			
			InsertPoints (wbc_klist_index+1),1, WBC_klist
			WBC_klist[wbc_klist_index][0] = kx
			WBC_klist[wbc_klist_index][1] = ky
			WBC_klist[wbc_klist_index][2] = kz
			
			wbc_klist_index += 1
			
			// read next line
			FReadLine fileRef, line
		
			// loop until you hit end of file or next k-point
			do
				FReadLine fileRef, line
			while(strlen(line) > 0 && cmpstr(line[2,2],".") != 0)
			
		else
			FReadLine fileRef, line
		endif
		
	while(strlen(line) > 0)
		
	Note WienHyperCube "This is a template for band character waves"
	Note WienHyperCube "Loaded on " + Date()+ " using:"
	Note WienHyperCube "EF = " + num2str(EF)
	Note WienHyperCube "BandMin = " + num2istr(BandMin)
	Note WienHyperCube "BandMax = " + num2istr(BandMax)
	Note WienHyperCube "latticeType = " + num2istr(latticetype)
	Note WienHyperCube "(N1,N2,N3) = (" + num2istr(N1) + ", " + num2istr(N2) + ", " + num2istr(N3) + ")"
	Note WienHyperCube "caxis = " + num2istr(caxis)
	Note WienHyperCube "primitiveBZ = " + num2istr(primitiveBZ)
	
	Close fileRef
	
	SetScale /I x, -1, 1, WienHyperCube
	SetScale /I y, -1, 1, WienHyperCube
	SetScale /I z, -1, 1, WienHyperCube

End

//------------------------------------------------------------------------

function WBCGetBandCharacters(WBC_atomsorbitals, BandMin, BandMax, latticeType, caxis, primitiveBZ)
	wave WBC_atomsorbitals
	Variable BandMin, BandMax,latticeType, caxis, primitiveBZ // Parameters from GenerateHCkmesh()
	
	wave WBC_cube0
	wave WBC_klist

//	Variable BandMin, BandMax 	// The indicies of the lowest and highest band the user wants to load (lowest band index is 1). Check "case.energy" for band numbers.
	
	Variable fileRef
	String line
	Variable kindex, atomindex, atomcount,k, caseindex
	variable kptscount = dimsize(WBC_klist,0)
	Variable BandNum = BandMax - BandMin + 1
	Variable kx, ky, kz, bandindex=-1
	variable casenum=dimsize(WBC_atomsorbitals,0)
	variable currval
	
	WBC_cube0=0
	for (caseindex=0;caseindex<casenum;caseindex+=1)
		duplicate/o WBC_cube0, $"WBC_bandchar"+num2str(caseindex)
	endfor
	
///////////////////////////////////////////////
	Open /T="????" /R fileRef

	line="bla"
	do
		FReadLine fileRef, line
		/// get number of atoms
		if(cmpstr(line[31,34],"NAT=") == 0)
			atomcount = str2num(line[37,40])
			make/o/n =(atomcount,10) WGC_kptdata

		elseif(cmpstr(line[1,4],"BAND") == 0) // If this line is a start of band data;
			if (str2num(line[7,12])<BandMin || str2num(line[7,12])>BandMax)
				continue
			endif
//			print line
			bandindex += 1
//			print bandindex
			
			for (kindex=0; kindex<kptscount; kindex+=1)
				for (atomindex=1; atomindex<=atomcount+1;atomindex+=1) //+1 for interstitial line
					FReadLine fileRef, line
					for (caseindex=0;caseindex<casenum; caseindex+=1)
						if (str2num(line[12,12]) == WBC_atomsorbitals[caseindex][0])
							if (WBC_atomsorbitals[caseindex][1]==0)
								currval = str2num(line[14, 20])
							else
								currval = str2num(line[25+(WBC_atomsorbitals[caseindex][1]-1)*8, 31+(WBC_atomsorbitals[caseindex][1]-1)*8])
							endif
							wave currentcase = $"WBC_bandchar"+num2str(caseindex)
							HC_AssignWithSymm(currentcase, WBC_klist[kindex][0], WBC_klist[kindex][1], WBC_klist[kindex][2], bandindex, currval, latticetype, caxis, primitiveBZ)
						endif	
					endfor
				endfor
			endfor
		endif
		
	while(strlen(line) > 0)

//	WBC_cube0=0
//	for (caseindex=0;caseindex<casenum;caseindex+=1)
//		wave currentcase = $"WBC_bandchar"+num2str(caseindex)
//		WBC_cube0+=currentcase
//	endfor
//	for (caseindex=0;caseindex<casenum;caseindex+=1)
//		wave currentcase = $"WBC_bandchar"+num2str(caseindex)
//		currentcase/=(WBC_cube0+1e-12)
//	endfor

////////////////////////////////////////////////////////////
	Close fileRef

End


