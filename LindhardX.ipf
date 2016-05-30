#pragma rtGlobals=1		// Use modern global access method.



////////////////////////////////////////////////////////////////////////////////////////
// Eric Monkman, July 2012
//
//
// This function computes the contribution to the Lindhard susceptibilty (assuming the matrix element == 1)
// following the formalism of Rath & Freeman PRB 11, 2109 (1975).
// The suscpetibility is calculated between two bands (which can be the same or distinct) for a given q vector.
// band# are 3D waves with x,y,z == kx,ky,kz and covering one full BZ
// They both must be the same size.
// kx, ky, kz should span 2pi/a
// Only works for systems with reciprocal lattice vectors along {x, y, z} respectively.

// band1_in and band2_in can be 3D Fermi surfaces from the Hypercube package
//
// Please note: this function is NOT symmetric in band1 <-> band2 (band1 is the occupied part, band2 the unnocupied)
////////////////////////////////////////////////////////////////////////////////////////
Function LindhardXq(band1_in, band2_in, qx, qy, qz, EF, [symmetry])

	Wave band1_in, band2_in
	Variable qx, qy, qz	// q should be in the 1st BZ (ie: between -1 and 1)
	Variable EF
	Variable symmetry	// A flag in case I attempt to include symmetries in the future
	symmetry = 0
	
	Variable nancount = 0
	Variable zerocount = 0
	Variable Omega, a, c
	Variable multiplier = 1
	
	Variable kx, ky, kz
	Variable D1, D2, D3
	Variable intT
	
	Variable i, j, k, l, count, Chi
	
	Duplicate /O band1_in band1
	Duplicate /O band2_in band2
	
	Variable numkx = DimSize(band1,0)
	Variable numky = DimSize(band1,1)
	Variable numkz = DimSize(band1,2)
	
	//If(qx == 0 && qy == 0 && qz == 0)
	//	print "Integral is ill defined for q == 0, I will ignore this q-point."
	//	return NaN
	//EndIf
	
	// Check that both bands are the same size and span the same k-window
	if(DimSize(band1,0) != DimSize(band2,0))
		Print "Band waves must both be the same size"
		return 0	
	endif
	if(DimSize(band1,1) != DimSize(band2,1))
		Print "Band waves must both be the same size"
		return 0	
	endif
	if(DimSize(band1,2) != DimSize(band2,2))
		Print "Band waves must both be the same size"
		return 0	
	endif
	if(qx > 1 || qx < -1)
		Print "q must lie in the 1st BZ"
		return 0	
	endif
	if(qy > 1 || qy < -1)
		Print "q must lie in the 1st BZ"
		return 0	
	endif
	if(qz > 1 || qz < -1)
		Print "q must lie in the 1st BZ"
		return 0	
	endif

	// Offset band2 by -q, including appropriate reciprocal lattice vectors
	Duplicate /O band2 band2_q, band2_qt
	
	LXShiftByQ(band2, band2_q, qx, qy, qz)
	band2 = band2_q
	
	// Integral is evaluated over a tetrahedral mesh.
	
	// Handle symmetry
	//If(symmetry == 1)	 // Reflection about all axes, use only one quadrant of BZ
	//	numkx = Ceil(DimSize(band1_in,0)/2)
	//	numky = Ceil(DimSize(band1_in,1)/2)
	//	numkz = Ceil(DimSize(band1_in,2)/2)
		
	//	Redimension /N = (numkx, numky, numkz) band1, band2_q
	//	SetScale /P x, (DimOffset(band1_in, 0)+(DimSize(band1_in,0)-numkx)*DimDelta(band1_in,0)), (DimDelta(band1_in,0)), band1
	//	SetScale /P y, (DimOffset(band1_in, 1)+(DimSize(band1_in,1)-numky)*DimDelta(band1_in,1)), (DimDelta(band1_in,1)), band1
	//	SetScale /P z, (DimOffset(band1_in, 2)+(DimSize(band1_in,2)-numkz)*DimDelta(band1_in,2)), (DimDelta(band1_in,2)), band1
		
	//	SetScale /P x, (DimOffset(band2, 0)+(DimSize(band2,0)-numkx)*DimDelta(band2,0)), (DimDelta(band2,0)), band2_q
	//	SetScale /P y, (DimOffset(band2, 1)+(DimSize(band2,1)-numky)*DimDelta(band2,1)), (DimDelta(band2,1)), band2_q
	//	SetScale /P z, (DimOffset(band2, 2)+(DimSize(band2,2)-numkz)*DimDelta(band2,2)), (DimDelta(band2,2)), band2_q
		
		
	//	band1 = band1_in[(DimSize(band1_in,0)-numkx+p)][(DimSize(band1_in,1)-numky+q)][(DimSize(band1_in,2)-numkz+r)]
	//	band2_q = band2[(DimSize(band1_in,0)-numkx+p)][(DimSize(band1_in,1)-numky+q)][(DimSize(band1_in,2)-numkz+r)]

	//	multiplier = 8
	//EndIf
	
	
	
	// Make an Ntetx12 wave with the coordinates of each tetrahedron's corners
	Variable ntet, NumHex
	NumHex = (numkx-1)*(numky-1)*(numkz-1) // number of hexahedra
	ntet = 5 * NumHex	// 5 tetrahedra make up each hexahedra
	
	Make /O /N = (NumHex,3) hexcoord		// Origin coordinate for each hexahedra
	Make /O /N = (ntet,12) tetcoord			// (kx,ky,kz) for each corner of all tetrahedra
	
	//Define hexahedra coordinates
	count = 0
	For(i = 0; i < numkx-1; i+=1)
		For(j = 0; j < numky-1; j+=1)
			For(k = 0; k < numkz-1; k+=1)
				hexcoord[count][0] =  i
				hexcoord[count][1] =  j
				hexcoord[count][2] =  k
				count += 1
			EndFor
		EndFor
	EndFor
	
	if(count != NumHex)
		Print "Error in defining hexcoord wave"
		return 0
	EndIf
	
	// Define tetrahedra coordinates
	For(i=0; i<NumHex; i+=1)
		LXBreakupCube(tetcoord, hexcoord, i)
	EndFor
	
	// Make a Nx1 waves to store the results of the integral for each tetrahedron
	Make /O /N = (ntet) TetInt
	TetInt = 0

	Variable occtets, unocctets
	Make /O /N=(4, 3) koc_1, koc_2, koc_3, kunoc_1, kunoc_2, kunoc_3
	Make /O /N = 4 Vi,  ET, kxT, kyT, kzT
	Make /O /N=3 k_0, k_1, k_2, k_3
	Vi = 0
	
	// Loop over each tetrahedron
	For(i = 0; i < ntet; i+=1)
		occtets = 0
		unocctets = 0
	
		// K values from the local master tetrahedron
		kxT = tetcoord[i][3*p]
		kyT = tetcoord[i][3*p+1] 
		kzT = tetcoord[i][3*p+2] 
				
		// Order corners of tetrahedra to have E4 <= E3 <= E2 <= E1 for the unshifted band
		ET[0] = band1[kxT[0]][kyT[0]][kzT[0]]
		ET[1] = band1[kxT[1]][kyT[1]][kzT[1]]
		ET[2] = band1[kxT[2]][kyT[2]][kzT[2]]
		ET[3] = band1[kxT[3]][kyT[3]][kzT[3]]
		
		//Reverse order sort, ET[3] has the lowest energy
		Sort /R ET, ET, kxT, kyT, kzT
		
		// Find the occupied sub-tetrahedra of the local master tetrahedron
		If(EF <= ET[3])
			//Whole tetrahedron is unoccupied, skip onto next one.
			koc_1 = NaN
			koc_2 = NaN
			koc_3 = NaN
			occtets = 0
		
		ElseIf(EF <= ET[2])
			// Occupied region is a single tetrahedron whose corners are:
			koc_1[0,2][0] = kxT[3] + (kxT[p] - kxT[3]) * (EF-ET[3])/(ET[p] - ET[3])
			koc_1[0,2][1] = kyT[3] + (kyT[p] - kyT[3]) * (EF-ET[3])/(ET[p] - ET[3])
			koc_1[0,2][2] = kzT[3] + (kzT[p] - kzT[3]) * (EF-ET[3])/(ET[p] - ET[3])
			koc_1[3][0] = kxT[3]
			koc_1[3][1] = kyT[3]
			koc_1[3][2] = kzT[3] 
		
			koc_2 = NaN
			koc_3 = NaN
			occtets = 1
		
		ElseIf(EF <= ET[1])
			// Occupied region is a sum of three tetrahedra whose corners are:
			koc_1[0,1][0] = kxT[3] + (kxT[p] - kxT[3]) * (EF-ET[3])/(ET[p] - ET[3])
			koc_1[0,1][1] = kyT[3] + (kyT[p] - kyT[3]) * (EF-ET[3])/(ET[p] - ET[3])
			koc_1[0,1][2] = kzT[3] + (kzT[p] - kzT[3]) * (EF-ET[3])/(ET[p] - ET[3])
			koc_1[2][0] = kxT[2]
			koc_1[2][1] = kyT[2]
			koc_1[2][2] = kzT[2]
			koc_1[3][0] = kxT[3]
			koc_1[3][1] = kyT[3]
			koc_1[3][2] = kzT[3]
			
			koc_2[0,1][0] = kxT[2] + (kxT[p] - kxT[2]) * (EF-ET[2])/(ET[p] - ET[2])
			koc_2[0,1][1] = kyT[2] + (kyT[p] - kyT[2]) * (EF-ET[2])/(ET[p] - ET[2])
			koc_2[0,1][2] = kzT[2] + (kzT[p] - kzT[2]) * (EF-ET[2])/(ET[p] - ET[2])
			koc_2[2][0] = kxT[3] + (kxT[1] - kxT[3]) * (EF-ET[3])/(ET[1] - ET[3])
			koc_2[2][1] = kyT[3] + (kyT[1] - kyT[3]) * (EF-ET[3])/(ET[1] - ET[3])
			koc_2[2][2] = kzT[3] + (kzT[1] - kzT[3]) * (EF-ET[3])/(ET[1] - ET[3])
			koc_2[3][0] = kxT[2]
			koc_2[3][1] = kyT[2]
			koc_2[3][2] = kzT[2]			
			
			koc_3[0][0] = kxT[2] + (kxT[0] - kxT[2]) * (EF-ET[2])/(ET[0] - ET[2])
			koc_3[0][1] = kyT[2] + (kyT[0] - kyT[2]) * (EF-ET[2])/(ET[0] - ET[2])
			koc_3[0][2] = kzT[2] + (kzT[0] - kzT[2]) * (EF-ET[2])/(ET[0] - ET[2])
			koc_3[1][0] = kxT[3] + (kxT[1] - kxT[3]) * (EF-ET[3])/(ET[1] - ET[3])
			koc_3[1][1] = kyT[3] + (kyT[1] - kyT[3]) * (EF-ET[3])/(ET[1] - ET[3])
			koc_3[1][2] = kzT[3] + (kzT[1] - kzT[3]) * (EF-ET[3])/(ET[1] - ET[3])	
			koc_3[2][0] = kxT[3] + (kxT[0] - kxT[3]) * (EF-ET[3])/(ET[0] - ET[3])
			koc_3[2][1] = kyT[3] + (kyT[0] - kyT[3]) * (EF-ET[3])/(ET[0] - ET[3])
			koc_3[2][2] = kzT[3] + (kzT[0] - kzT[3]) * (EF-ET[3])/(ET[0] - ET[3])
			koc_3[3][0] = kxT[2]
			koc_3[3][1] = kyT[2]
			koc_3[3][2] = kzT[2]	
			
			occtets = 3
		
		ElseIf(EF < ET[0])
			// Occupied region is a sum of three tetrahedra whose corners are:
			koc_1[0,2][0] = kxT[0] + (kxT[3-p] - kxT[0]) * (EF-ET[0])/(ET[3-p] - ET[0])
			koc_1[0,2][1] = kyT[0] + (kyT[3-p] - kyT[0]) * (EF-ET[0])/(ET[3-p] - ET[0])
			koc_1[0,2][2] = kzT[0] + (kzT[3-p] - kzT[0]) * (EF-ET[0])/(ET[3-p] - ET[0])
			koc_1[3][0] = kxT[3]
			koc_1[3][1] = kyT[3]
			koc_1[3][2] = kzT[3]
			
			koc_2[0,1][] = koc_1[p+1][q]
			koc_2[2][0] = kxT[3]
			koc_2[2][1] = kyT[3]
			koc_2[2][2] = kzT[3]
			koc_2[3][0] = kxT[1]
			koc_2[3][1] = kyT[1]
			koc_2[3][2] = kzT[1]
			
			koc_3[0][] = koc_1[1][q]
			koc_3[1][0] = kxT[3]
			koc_3[1][1] = kyT[3]
			koc_3[1][2] = kzT[3]
			koc_3[2][0] = kxT[2]
			koc_3[2][1] = kyT[2]
			koc_3[2][2] = kzT[2]
			koc_3[3][0] = kxT[1]
			koc_3[3][1] = kyT[1]
			koc_3[3][2] = kzT[1]
			
			occtets = 3
		
		Else
			//Whole tetrahedron is occupied
			koc_1[][0] = kxT[p]
			koc_1[][1] = kyT[p]
			koc_1[][2] = kzT[p] 
			
			koc_2 = NaN
			koc_3 = NaN
			occtets = 1
			
		EndIf
		
		//Within each sub-tetrahedra, find the unnocupied volume of q-shifted band.		
		For(j = 1; j <= occtets; j+=1)
			Wave k_oc = $("koc_"+num2str(j))
		
			// Calculate the energies at the corners of the tetrahedra
			For(k=0; k<4; k+=1)
				kx =  DimOffset(band2_q,0)+DimDelta(band2_q,0)*k_oc[k][0]
				ky =  DimOffset(band2_q,1)+DimDelta(band2_q,1)*k_oc[k][1]
				kz =  DimOffset(band2_q,2)+DimDelta(band2_q,2)*k_oc[k][2]
			
				ET[k] = interp3D(band2_q, kx, ky, kz)
			endfor
			kxT = k_oc[p][0]
			kyT = k_oc[p][1]
			kzT = k_oc[p][2]
			
			// Order corners of tetrahedra to have E4 <= E3 <= E2 <= E1 for the shifted band
			//Reverse order sort, ET[3] has the lowest energy
			Sort /R ET, ET, kxT, kyT, kzT
		
			// Break up the occupied tetrahedra into tetrahedra that are unoccupied for the q-shifted band
			unocctets = LXSplitupUnoccTet(EF, ET, kxT, kyT, kzT, kunoc_1, kunoc_2, kunoc_3)
		
			// For each unoccupied tetrahedra, perform the volume integral
			For(k = 1; k <= unocctets; k+=1)
				Wave kunoc = $("kunoc_"+num2str(k))
			
				// Convert k-points of corners into correct units and evaluate volume
				k_0 =  DimOffset(band1,p)+DimDelta(band1,p)*kunoc[0][p]
				k_1 =  DimOffset(band1,p)+DimDelta(band1,p)*kunoc[1][p] - k_0[p]
				k_2 =  DimOffset(band1,p)+DimDelta(band1,p)*kunoc[2][p] - k_0[p]
				k_3 =  DimOffset(band1,p)+DimDelta(band1,p)*kunoc[3][p] - k_0[p]
		
				Cross k_2, k_3
				Wave W_cross
				Omega = Abs((k_1[0] * W_cross[0]+k_1[1] * W_cross[1]+k_1[2] * W_cross[2])/6)

				// Calculate energy difference at each corner
				For(l=0; l<4; l+=1)
					kx =  DimOffset(band2_q,0)+DimDelta(band2_q,0)*kunoc[l][0]
					ky =  DimOffset(band2_q,1)+DimDelta(band2_q,1)*kunoc[l][1]
					kz =  DimOffset(band2_q,2)+DimDelta(band2_q,2)*kunoc[l][2]
			
					Vi[l] = interp3D(band2_q, kx, ky, kz)
					
					kx =  DimOffset(band1,0)+DimDelta(band1,0)*kunoc[l][0]
					ky =  DimOffset(band1,1)+DimDelta(band1,1)*kunoc[l][1]
					kz =  DimOffset(band1,2)+DimDelta(band1,2)*kunoc[l][2]
			
					Vi[l] -= interp3D(band1, kx, ky, kz)
	
				endfor
				
				Sort Vi, Vi
				D1 = 0
				D2 = 0
				D3 = 0
				
				// Go through possible special cases
				If(Vi[0]*Vi[1]*Vi[2]*Vi[3] == 0)	// If any Vi are zero, we neglect this tetrahedron
					intT = 0
					zerocount += 1
				
				ElseIf(Vi[0] == Vi[1] && Vi[1] == Vi[2] && Vi[2] == Vi[3]) // all four are equal
					intT = Omega/Vi[0]
				
				ElseIf(Vi[0] == Vi[1] && Vi[1] == Vi[2])	// three are equal
					intT = 3*Omega*( (Vi[3]^2/(Vi[0]-Vi[3])^3)*Ln(Abs(Vi[0]/Vi[3])) + (1.5*Vi[3]^2 + 0.5*Vi[0]^2 - 2*Vi[0]*Vi[3])/ (Vi[0] - Vi[3])^3 )
				
				ElseIf(Vi[1] == Vi[2] && Vi[2] == Vi[3])	// three are equal
					intT = 3*Omega*( (Vi[0]^2/(Vi[3]-Vi[0])^3)*Ln(Abs(Vi[3]/Vi[0])) + (1.5*Vi[0]^2 + 0.5*Vi[3]^2 - 2*Vi[3]*Vi[0])/ (Vi[3] - Vi[0])^3 )
				
				ElseIf(Vi[0] == Vi[1] && Vi[2] == Vi[3])	// two pairs
					intT = 3*Omega*( (2*Vi[1]*Vi[2]/(Vi[1]-Vi[2])^3)*Ln(Abs(Vi[2]/Vi[1])) + (Vi[1] + Vi[2])/(Vi[1] - Vi[2])^2  )
				
				ElseIf(Vi[0] == Vi[1])	// One pair
					intT = 3*Omega*(  Vi[2]^2/( (Vi[2]-Vi[0])^2*(Vi[2] - Vi[3]) )*Ln(Abs(Vi[2]/Vi[0])) + Vi[3]^2/( (Vi[3]-Vi[0])^2*(Vi[3] - Vi[2]) )*Ln(Abs(Vi[3]/Vi[0])) + Vi[0]/( (Vi[2]-Vi[0])*(Vi[3]-Vi[0]) )   )
				
				ElseIf(Vi[1] == Vi[2])	// One pair
					intT = 3*Omega*(  Vi[0]^2/( (Vi[0]-Vi[1])^2*(Vi[0] - Vi[3]) )*Ln(Abs(Vi[0]/Vi[1])) + Vi[3]^2/( (Vi[3]-Vi[1])^2*(Vi[3] - Vi[0]) )*Ln(Abs(Vi[3]/Vi[1])) + Vi[1]/( (Vi[0]-Vi[1])*(Vi[3]-Vi[1]) )   )
				
				ElseIf(Vi[2] == Vi[3])	// One pair
					intT = 3*Omega*(  Vi[0]^2/( (Vi[0]-Vi[2])^2*(Vi[0] - Vi[1]) )*Ln(Abs(Vi[0]/Vi[2])) + Vi[1]^2/( (Vi[1]-Vi[2])^2*(Vi[1] - Vi[0]) )*Ln(Abs(Vi[1]/Vi[2])) + Vi[2]/( (Vi[0]-Vi[2])*(Vi[1]-Vi[2]) )   )
				
				Else	// Generic formula
					D1 = (Vi[0] - Vi[3])*(Vi[0]-Vi[2])*(Vi[0]-Vi[1])
					D2 = (Vi[1] - Vi[3])*(Vi[1]-Vi[2])*(Vi[1]-Vi[0])
					D3 = (Vi[2] - Vi[3])*(Vi[2]-Vi[1])*(Vi[2]-Vi[0])
				
					intT = 3*Omega*( (Vi[0]^2/D1)*Ln(Abs(Vi[0]/Vi[3])) +  (Vi[1]^2/D2)*Ln(Abs(Vi[1]/Vi[3])) +  (Vi[2]^2/D3)*Ln(Abs(Vi[2]/Vi[3])) )
			
				EndIf
				
				// Ignore (but keep track of) NaN's
				if(intT*0 == 0)
					TetInt[i] += intT
				Else
					//print "Next cycle"
					//print kunoc
					//Print ET, kxT, kyT, kzT
					//print Vi
					//print D1, D2, D3
					//print i, j, k
					nancount +=1
				EndIf
				
			EndFor
		EndFor
	EndFor
	
	// Sum over the integral-wave to obtain the resultant chi
	Chi = 0
	For(i = 0; i < ntet; i+=1)
		Chi += TetInt[i]
	EndFor
	
	KillWaves tetcoord, hexcoord, W_cross, k_0, k_1, k_2, k_3, ET, koc_1, koc_2, koc_3, kunoc_1, kunoc_2, kunoc_3, kxT, kyT, kzT, TetInt, Vi, band1, band2, band2_q
	//If(nancount > 0 || zerocount > 0)
	//	Print "Found " + num2str(nancount) + " NaN values and " + num2str(zerocount)+ " neglected sub-tetrahedra out of " + num2str(2*ntet) + " total tetrahedra"
	//Endif
	
	return Chi*multiplier

End


// Divides a cube into five tetrahedra
Function LXBreakupCube(tetcoord, hexcoord, i)
	wave tetcoord, hexcoord
	variable i

	//Tet 1
	tetcoord[5*i][0] = hexcoord[i][0]
	tetcoord[5*i][1] = hexcoord[i][1]
	tetcoord[5*i][2] = hexcoord[i][2]
	tetcoord[5*i][3] = hexcoord[i][0] + 1
	tetcoord[5*i][4] = hexcoord[i][1]
	tetcoord[5*i][5] = hexcoord[i][2] + 1
	tetcoord[5*i][6] = hexcoord[i][0]
	tetcoord[5*i][7] = hexcoord[i][1]
	tetcoord[5*i][8] = hexcoord[i][2] + 1
	tetcoord[5*i][9] = hexcoord[i][0]
	tetcoord[5*i][10] = hexcoord[i][1] + 1
	tetcoord[5*i][11] = hexcoord[i][2] + 1

	//Tet 2
	tetcoord[5*i+1][0] = hexcoord[i][0]
	tetcoord[5*i+1][1] = hexcoord[i][1]
	tetcoord[5*i+1][2] = hexcoord[i][2]
	tetcoord[5*i+1][3] = hexcoord[i][0]+1
	tetcoord[5*i+1][4] = hexcoord[i][1]
	tetcoord[5*i+1][5] = hexcoord[i][2]
	tetcoord[5*i+1][6] = hexcoord[i][0]+1
	tetcoord[5*i+1][7] = hexcoord[i][1]+1
	tetcoord[5*i+1][8] = hexcoord[i][2]
	tetcoord[5*i+1][9] = hexcoord[i][0]+1
	tetcoord[5*i+1][10] = hexcoord[i][1]
	tetcoord[5*i+1][11] = hexcoord[i][2]+1
	
	//Tet 3
	tetcoord[5*i+2][0] = hexcoord[i][0]
	tetcoord[5*i+2][1] = hexcoord[i][1]
	tetcoord[5*i+2][2] = hexcoord[i][2]
	tetcoord[5*i+2][3] = hexcoord[i][0]
	tetcoord[5*i+2][4] = hexcoord[i][1]+1
	tetcoord[5*i+2][5] = hexcoord[i][2]
	tetcoord[5*i+2][6] = hexcoord[i][0]+1
	tetcoord[5*i+2][7] = hexcoord[i][1]+1
	tetcoord[5*i+2][8] = hexcoord[i][2]
	tetcoord[5*i+2][9] = hexcoord[i][0]
	tetcoord[5*i+2][10] = hexcoord[i][1]+1
	tetcoord[5*i+2][11] = hexcoord[i][2]+1
	
	//Tet 4
	tetcoord[5*i+3][0] = hexcoord[i][0]+1
	tetcoord[5*i+3][1] = hexcoord[i][1]+1
	tetcoord[5*i+3][2] = hexcoord[i][2]
	tetcoord[5*i+3][3] = hexcoord[i][0]+1
	tetcoord[5*i+3][4] = hexcoord[i][1]+1
	tetcoord[5*i+3][5] = hexcoord[i][2]+1
	tetcoord[5*i+3][6] = hexcoord[i][0]+1
	tetcoord[5*i+3][7] = hexcoord[i][1]
	tetcoord[5*i+3][8] = hexcoord[i][2]+1
	tetcoord[5*i+3][9] = hexcoord[i][0]
	tetcoord[5*i+3][10] = hexcoord[i][1]+1
	tetcoord[5*i+3][11] = hexcoord[i][2]+1
	
	//Tet 5
	tetcoord[5*i+4][0] = hexcoord[i][0]
	tetcoord[5*i+4][1] = hexcoord[i][1]
	tetcoord[5*i+4][2] = hexcoord[i][2]
	tetcoord[5*i+4][3] = hexcoord[i][0]+1
	tetcoord[5*i+4][4] = hexcoord[i][1]+1
	tetcoord[5*i+4][5] = hexcoord[i][2]
	tetcoord[5*i+4][6] = hexcoord[i][0]+1
	tetcoord[5*i+4][7] = hexcoord[i][1]
	tetcoord[5*i+4][8] = hexcoord[i][2]+1
	tetcoord[5*i+4][9] = hexcoord[i][0]
	tetcoord[5*i+4][10] = hexcoord[i][1]+1
	tetcoord[5*i+4][11] = hexcoord[i][2]+1

End


// Shifts a 3D band by given q vector, linearly interpolating the final result
Function LXShiftByQ(band2, band2_q, qx, qy, qz)
	Wave band2, band2_q
	Variable qx, qy, qz
	
	Duplicate /O band2 band2_qmin, band2_qmax, band2_qt
	band2_q = band2
	
	Variable qxpmin, qypmin, qzpmin, qxpmax, qypmax, qzpmax	// Translation vectors to use for the data
	qxpmin = floor(qx/DimDelta(band2, 0))
	qypmin = floor(qy/DimDelta(band2, 1))
	qzpmin = floor(qz/DimDelta(band2, 2))
	qxpmax = ceil(qx/DimDelta(band2, 0))
	qypmax = ceil(qy/DimDelta(band2, 1))
	qzpmax = ceil(qz/DimDelta(band2, 2))
	
	Variable weightx, weighty, weightz
	weightx = qx/DimDelta(band2, 0) - floor(qx/DimDelta(band2, 0))
	weighty = qy/DimDelta(band2, 1) - floor(qy/DimDelta(band2, 1))
	weightz = qz/DimDelta(band2, 2) - floor(qz/DimDelta(band2, 2))
	
	Variable numkx = DimSize(band2,0)
	Variable numky = DimSize(band2,1)
	Variable numkz = DimSize(band2,2)
	
	//Translate by qxmin
	If(qxpmin > 0)
		band2_qmin[0,(numkx-qxpmin-1)][][] = band2_qt[p+qxpmin][q][r]
		band2_qmin[(numkx-qxpmin),][][] = band2_qt[p+qxpmin - numkx][q][r]
	Else
		band2_qmin[0,(-qxpmin-1)][][] = band2_qt[p+qxpmin +numkx][q][r]
		band2_qmin[(-qxpmin),][][] = band2_qt[p+qxpmin][q][r]
	EndIf
	
	//Translate by qxmax
	If(qxpmax > 0)
		band2_qmax[0,(numkx-qxpmax-1)][][] = band2_qt[p+qxpmax][q][r]
		band2_qmax[(numkx-qxpmax),][][] = band2_qt[p+qxpmax - numkx][q][r]
	Else
		band2_qmax[0,(-qxpmax-1)][][] = band2_qt[p+qxpmax +numkx][q][r]
		band2_qmax[(-qxpmax),][][] = band2_qt[p+qxpmax][q][r]
	EndIf
		
	//Interpolate in x direction
	band2_q = band2_qmin + (band2_qmax - band2_qmin)*(weightx)
	band2_qt = band2_q
	
	//Translate by qypmin
	If(qypmin > 0)
		band2_qmin[][0,(numky-qypmin-1)][] = band2_qt[p][q+qypmin][r]
		band2_qmin[][(numky-qypmin),][] = band2_qt[p][q+qypmin - numky][r]
	Else
		band2_qmin[][0,(-qypmin-1)][] = band2_qt[p][q+qypmin +numky][r]
		band2_qmin[][(-qypmin),][] = band2_qt[p][q+qypmin][r]
	EndIf
	
	//Translate by qypmax
	If(qypmax > 0)
		band2_qmax[][0,(numky-qypmax-1)][] = band2_qt[p][q+qypmax][r]
		band2_qmax[][(numky-qypmax),][] = band2_qt[p][q+qypmax - numky][r]
	Else
		band2_qmax[][0,(-qypmax-1)][] = band2_qt[p][q+qypmax +numky][r]
		band2_qmax[][(-qypmax),][] = band2_qt[p][q+qypmax][r]
	EndIf
	
	//Interpolate in y direction
	band2_q = band2_qmin + (band2_qmax - band2_qmin)*(weighty)
	band2_qt = band2_q
	
	//Translate by qzpmin
	If(qzpmin > 0)
		band2_qmin[][][0,(numkz-qzpmin-1)] = band2_qt[p][q][r+qzpmin]
		band2_qmin[][][(numkz-qzpmin),] = band2_qt[p][q][r+qzpmin - numkz]
	Else
		band2_qmin[][][0,(-qzpmin-1)] = band2_qt[p][q][r+qzpmin +numkz]
		band2_qmin[][][(-qzpmin),] = band2_qt[p][q][r+qzpmin]
	EndIf
	
	//Translate by qzpmax
	If(qzpmax > 0)
		band2_qmax[][][0,(numkz-qzpmax-1)] = band2_qt[p][q][r+qzpmax]
		band2_qmax[][][(numkz-qzpmax),] = band2_qt[p][q][r+qzpmax - numkz]
	Else
		band2_qmax[][][0,(-qzpmax-1)] = band2_qt[p][q][r+qzpmax +numkz]
		band2_qmax[][][(-qzpmax),] = band2_qt[p][q][r+qzpmax]
	EndIf

	//Interpolate in z direction
	band2_q = band2_qmin + (band2_qmax - band2_qmin)*(weightz)

	Killwaves band2_qmin, band2_qmax, band2_qt
	
End


// Breaks up a tetrahedron into at most three unnocupied tetrahedra
Function LXSplitupUnoccTet(EF, ET, kxT, kyT, kzT, kunoc_1, kunoc_2, kunoc_3)
	Variable EF
	Wave ET, kxT, kyT, kzT, kunoc_1, kunoc_2, kunoc_3
	Variable unocctets
	
	// Find the occupied sub-tetrahedra of the local master tetrahedron
	If(EF <= ET[3])
		//Whole tetrahedron is unoccupied.
		kunoc_1[][0] = kxT[p]
		kunoc_1[][1] = kyT[p]
		kunoc_1[][2] = kzT[p] 
			
		kunoc_2 = NaN
		kunoc_3 = NaN
		unocctets = 1
		
	ElseIf(EF <= ET[2]) //ET0, ET1, and ET2 could all be equal, but must all be > ET3
		// unoccupied region is a sum of three tetrahedra whose corners are:
		kunoc_1[0][0] = kxT[0]
		kunoc_1[0][1] = kyT[0]
		kunoc_1[0][2] = kzT[0] 
		kunoc_1[1][0] = kxT[0] + (kxT[3] - kxT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[1][1] = kyT[0] + (kyT[3] - kyT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[1][2] = kzT[0] + (kzT[3] - kzT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[2][0] = kxT[3] + (kxT[2] - kxT[3]) * (EF-ET[3])/(ET[2] - ET[3])
		kunoc_1[2][1] = kyT[3] + (kyT[2] - kyT[3]) * (EF-ET[3])/(ET[2] - ET[3])
		kunoc_1[2][2] = kzT[3] + (kzT[2] - kzT[3]) * (EF-ET[3])/(ET[2] - ET[3])
		kunoc_1[3][0] = kxT[3] + (kxT[1] - kxT[3]) * (EF-ET[3])/(ET[1] - ET[3])
		kunoc_1[3][1] = kyT[3] + (kyT[1] - kyT[3]) * (EF-ET[3])/(ET[1] - ET[3])
		kunoc_1[3][2] = kzT[3] + (kzT[1] - kzT[3]) * (EF-ET[3])/(ET[1] - ET[3])
		
		kunoc_2[0][0] = kxT[0]
		kunoc_2[0][1] = kyT[0]
		kunoc_2[0][2] = kzT[0] 
		kunoc_2[1][0] = kxT[2]
		kunoc_2[1][1] = kyT[2]
		kunoc_2[1][2] = kzT[2] 
		kunoc_2[2][] = kunoc_1[3][q]
		kunoc_2[3][] = kunoc_1[2][q]

		kunoc_3[0][0] = kxT[0]
		kunoc_3[0][1] = kyT[0]
		kunoc_3[0][2] = kzT[0] 
		kunoc_3[1][0] = kxT[1]
		kunoc_3[1][1] = kyT[1]
		kunoc_3[1][2] = kzT[1] 
		kunoc_3[2][0] = kxT[2]
		kunoc_3[2][1] = kyT[2]
		kunoc_3[2][2] = kzT[2] 
		kunoc_3[3][] = kunoc_1[3][q]
			
		unocctets = 3
		
	ElseIf(EF <= ET[1])//ET2 and ET3 could be equal, and ET1 and ET0 could be equal
		// Unoccupied region is a sum of three tetrahedra whose corners are:
		kunoc_1[0][0] = kxT[0]
		kunoc_1[0][1] = kyT[0]
		kunoc_1[0][2] = kzT[0] 
		kunoc_1[1][0] = kxT[0] + (kxT[3] - kxT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[1][1] = kyT[0] + (kyT[3] - kyT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[1][2] = kzT[0] + (kzT[3] - kzT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[2][0] = kxT[3] + (kxT[1] - kxT[3]) * (EF-ET[3])/(ET[1] - ET[3])
		kunoc_1[2][1] = kyT[3] + (kyT[1] - kyT[3]) * (EF-ET[3])/(ET[1] - ET[3])
		kunoc_1[2][2] = kzT[3] + (kzT[1] - kzT[3]) * (EF-ET[3])/(ET[1] - ET[3])
		kunoc_1[3][0] = kxT[0] + (kxT[2] - kxT[0]) * (EF-ET[0])/(ET[2] - ET[0])
		kunoc_1[3][1] = kyT[0] + (kyT[2] - kyT[0]) * (EF-ET[0])/(ET[2] - ET[0])
		kunoc_1[3][2] = kzT[0] + (kzT[2] - kzT[0]) * (EF-ET[0])/(ET[2] - ET[0])

		kunoc_2[0][0] = kxT[0]
		kunoc_2[0][1] = kyT[0]
		kunoc_2[0][2] = kzT[0] 
		kunoc_2[1][0] = kxT[1]
		kunoc_2[1][1] = kyT[1]
		kunoc_2[1][2] = kzT[1] 
		kunoc_2[2][] = kunoc_1[2][q]
		kunoc_2[3][] = kunoc_1[3][q]	
			
		kunoc_3[0][0] = kxT[1]
		kunoc_3[0][1] = kyT[1]
		kunoc_3[0][2] = kzT[1] 
		kunoc_3[1][] = kunoc_1[3][q]
		kunoc_3[2][] = kunoc_1[2][q]
		kunoc_3[3][0] = kxT[1] + (kxT[2] - kxT[1]) * (EF-ET[1])/(ET[2] - ET[1])
		kunoc_3[3][1] = kyT[1] + (kyT[2] - kyT[1]) * (EF-ET[1])/(ET[2] - ET[1])
		kunoc_3[3][2] = kzT[1] + (kzT[2] - kzT[1]) * (EF-ET[1])/(ET[2] - ET[1])
			
		unocctets = 3
		
	ElseIf(EF < ET[0])//ET1, ET2, and ET3 could all be equal, but must all be < ET0
		// unccupied region is a single tetrahedron whose corners are:
		kunoc_1[0][0] = kxT[0]
		kunoc_1[0][1] = kyT[0]
		kunoc_1[0][2] = kzT[0] 
		kunoc_1[1][0] = kxT[0] + (kxT[3] - kxT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[1][1] = kyT[0] + (kyT[3] - kyT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[1][2] = kzT[0] + (kzT[3] - kzT[0]) * (EF-ET[0])/(ET[3] - ET[0])
		kunoc_1[2][0] = kxT[0] + (kxT[1] - kxT[0]) * (EF-ET[0])/(ET[1] - ET[0])
		kunoc_1[2][1] = kyT[0] + (kyT[1] - kyT[0]) * (EF-ET[0])/(ET[1] - ET[0])
		kunoc_1[2][2] = kzT[0] + (kzT[1] - kzT[0]) * (EF-ET[0])/(ET[1] - ET[0])
		kunoc_1[3][0] = kxT[0] + (kxT[2] - kxT[0]) * (EF-ET[0])/(ET[2] - ET[0])
		kunoc_1[3][1] = kyT[0] + (kyT[2] - kyT[0]) * (EF-ET[0])/(ET[2] - ET[0])
		kunoc_1[3][2] = kzT[0] + (kzT[2] - kzT[0]) * (EF-ET[0])/(ET[2] - ET[0])
			
		kunoc_2 = NaN
		kunoc_3 = NaN
		unocctets = 1
		
	Else
		//Whole tetrahedron is occupied
		kunoc_1 = NaN
		kunoc_2 = NaN
		kunoc_3 = NaN
		unocctets = 0
			
	EndIf
		
	return unocctets

End