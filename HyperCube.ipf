#pragma rtGlobals=1		// Use modern global access method.

/////////////////////////////////////////////////////////////////
// A package to load and analyse Wien2K data on 3D materials.
//
// Supports cubic, tetragonal, and orthorhombic lattices with inversion symmetry.
//
// Contents:
// GenerateHCkmesh - generates a kmesh for calculating the hypercube
// LoadHyperCube - loads data
// HyperCubeBandsAtEF - finds all of the bands that cut through a given energy

// HyperCube3DFS - extracts bands as 3D waves from the HyperCube (use for FS plotting)
// HyperCube2DFS - extracts 2D cuts at fixed kz from the HyperCube (use for FS plotting)
// HyperCubeProjectFS - extracts a 2D kz-smeared Fermi Surface from a HyperCube3DFS wave

// HCGenerateKlist - generates a klist along a curved path in k-space corresponding to an ARPES cut
// HyperCubeSlice_kz - extracts bandstructure for spaghetti plots along arbitrary k-space cuts at fixed kz.
// HyperCubePlotSpaghetti - plots the bandstructure from HyperCubeSlice_kz
// HyperCubeSlice_kzRange - HyperCubeSlice_kz for a range in kz
// HyperCubePlotSpaghetti3D - plots the bandstructure from HyperCubeSlice_kzRange
// HCsimulateARPES - Includes kz-broadening and a simple lifetime correction to simulate an ARPES spectrum
//
// Tools for reconstructed data:
// RotateHC - rotates a hypercube (ie: mainly for systems whose c-axis is not axis #3)
// HCc2x2Reconstruct - Expands a hypercube that has a c2x2 reconstruction in the xy plane --- Untested, use caution
// HClinearReconstruct - Expands a hypercube that has a reconstruction along one axis (ie: 1x1x2) --- Untested, use caution
// c2x2Reconstruct - Expands a c2x2 reconstructed 2D FS back into the unreconstructed BZ
// Batchc2x2 - Same as c2x2Reconstruct, but for batch-processing of multiple FSs
// ShadowBandsc2x2 - Breaks up contour plots of a FS into 'real' and 'shadow' bands
//
// Eric Monkman, April 2012
/////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////
// GenerateHCkmesh() and LoadHyperCube() import a 'hypercube' of wien2k data into IGOR.
// The result is a 4D data file, with kx, ky, kz, and band index as the dimensions, and energy as the 'value'
// Energy is in eV with EF = 0. The range of the hypercube in k-space is -pi/a to pi/a.
//
// Currently supported: Cubic, tetragonal, and orthorhombic lattices.
//
// To use:
// (1) Run an SCF calculation to convergence.
// (2) Use xcrysden to verify that your primitive and conventional Brilluoin zones are the same. If they aren't, one of the "primitiveBZ" options below may be necessary or your case may not be supported.
// (3) Using GenerateHCkmesh(), generate the file xcrysden.klist.
// (4) Upload xcrysden.klist to the case directory.
// (5) In w2web, under Tasks >> Bandstructure:														
//    	(a) Select "from xcrysden" in the dropdown box													
//    	(b) Click the "create case.klist_band" button														
//    	(c) Click the following buttons:																	
//         	(i)   "xlapw1 -band" for a non-spin-polarized calculation, OR									
//         	(ii)  "xlapw1 -band -up" and/or "xlapw1 -band -dn" for a spin-polarized calculation, AND			
//         	(iii) "xlapwso -band" if your calculation includes spin-orbit coupling	
// (6) In w2web, download the file(s): "case.energy", "case.energyup", "case.energydn", "case.energyso"
// (7) Use LoadHyperCube() to load the file(s)
//
//////////////////////////////////////////////////
Function GenerateHCkmesh(latticeType, N1, N2, N3, caxis, primitiveBZ)
	Variable latticeType	// 0 = cubic,1 = tetragonal (a = b != c), 2 = orthorhombic (a != b != c)
	Variable N1, N2, N3	// Max number of points along axes 1, 2, and 3 (usually x,y,z). For cubic; N2, N3 are ignored. For tetragonal; one N is ignored based on c-axis.
	Variable cAxis		// c-axis in WIEN2k k-space (1, 2, or 3; usually 3). Use caution here: if caxis != 3 in struct file, it may still be = 3 in k space. The best way to tell is by using xcrysden to generate a few k points.
	Variable primitiveBZ	// 1 if primitive BZ = conventional zone (the typical case); 2 if primitive BZ = c2x2 of conventional in xy-plane (Nx must = Ny); 3 to cover a BZ that spans 4 pi/a in xy-plane.
	
	Variable xextend, yextend
	xextend = 1
	yextend = 1
	
	Make /O /T /N = 3 GenerateHCkmesh_record
	GenerateHCkmesh_record[0] = "Generated kmesh with (lattice, N1, N2, N3, caxis, primitiveBZ) as follows: "
	GenerateHCkmesh_record[1] = num2istr(latticeType)  + ", "+ num2istr(N1) + ", "+ num2istr(N2) + ", "+ num2istr(N3)+ ", "+ num2istr(caxis) +", "+ num2istr(primitiveBZ)
	GenerateHCkmesh_record[2] = date() + ", " + time()

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
	
	If(primitiveBZ == 3)
		xextend = 2
		yextend = 2
		If(latticetype != 2)
			Print "PrimitiveBZ = 3 is only enabled for orthorhombic lattices."
			return 0
		EndIf
	EndIf
	
	Variable i, j, l, k
	Variable ip, jp
	String a, b, c, d
	Variable Nx, Ny, Nz
	
	// N1, N2, and N3 should be odd for this kmesh
	N1 = 2*floor((N1+1) /2) - 1
	N2 = 2*floor((N2+1) /2) - 1
	N3 = 2*floor((N3+1) /2) - 1
	
	
	//////////////////////////////////
	If(primitiveBZ == 1 || primitiveBZ == 3)		
	//////////////////////////////////
	switch(latticeType)
		case 0: // cubic (a = b = c)
			Nx = floor((N1+1) /2)
			Ny = Nx
			Nz = Nx
		
			Make /O /T /N=(Round(Nx*(Nx+1)*(Nx+2)/6)) HCklist_tmp
			
			Print "Number of non-equivalent k-points in BZ = ", num2istr(DimSize(HCklist_tmp,0))
			
			if(2*(Nx-1) < 10)
				d = "   " + num2istr((Nx-1)*10*2)
			elseif(2*(Nx-1) < 100)
				d = "  " + num2istr((Nx-1)*10*2)
			elseif(2*(Nx-1) < 1000)
				d = " " + num2istr((Nx-1)*10*2)
			else
				Print "Error: kmesh this fine is not supported."
				return(0)
			endif
			
			k = 0
			
			for(i = 0; i <= Nx-1; i += 1)
			
				if(i == 0)
					a = "    0"
				elseif(i < 10)
					a = "   " + num2istr(i*10)
				elseif(i < 100)
					a = "  " + num2istr(i*10)
				else
					a = " " + num2istr(i*10)
				endif
				
				for(j = 0; j <= i; j += 1)
				
					if(j == 0)
						b = "    0"
					elseif(j < 10)
						b = "   " + num2istr(j*10)
					elseif(j < 100)
						b = "  " + num2istr(j*10)
					else
						b = " " + num2istr(j*10)
					endif
									
					for(l = 0; l <= j; l += 1)
					
						if(l == 0)
							c = "    0"
						elseif(l < 10)
							c = "   " + num2istr(l*10)
						elseif(l < 100)
							c = "  " + num2istr(l*10)
						else
							c = " " + num2istr(l*10)
						endif
					
						HCklist_tmp[k] = "          " + a + b + c + d + "  1.0"
					
						k += 1
					endfor
				endfor
			endfor
			
			break
		case 1: // tetragonal (a = b != c)
			
			If(caxis == 1)
				Nx = floor((N2+1) /2)
				Ny = Nx
				Nz = floor((N1+1) /2)
			elseif(caxis ==2)
				Nx = floor((N1+1) /2)
				Ny = Nx
				Nz = floor((N2+1) /2)
			elseif(caxis ==3)
				Nx = floor((N1+1) /2)
				Ny = Nx
				Nz = floor((N3+1) /2)
			endif
		
			Make /O /T /N=(Round(Nz*Nx*(Nx+1)/2)) HCklist_tmp
			
			Print "Number of non-equivalent k-points in BZ = ", num2istr(DimSize(HCklist_tmp,0))
			
			if((Nx-1)*(Nz-1) < 10)
				d = "         " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 100)
				d = "        " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 1000)
				d = "       " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 10000)
				d = "      " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 100000)
				d = "     " + num2istr((Nx-1)*(Nz-1)*2)
			else
				Print "Error: kmesh this fine is not supported."
				return(0)
			endif
			
			k = 0
			
			for(i = 0; i <= Nx-1; i += 1)
			
				if(i == 0)
					a = "         0"
				elseif(i*(Nz-1) < 10)
					a = "         " + num2istr(i*(Nz-1))
				elseif(i*(Nz-1) < 100)
					a = "        " + num2istr(i*(Nz-1))
				elseif(i*(Nz-1) < 1000)
					a = "       " + num2istr(i*(Nz-1))
				elseif(i*(Nz-1) < 10000)
					a = "      " + num2istr(i*(Nz-1))
				elseif(i*(Nz-1) < 100000)
					a = "     " + num2istr(i*(Nz-1))
				endif
				
				for(j = 0; j <= i; j += 1)
				
					if(j == 0)
						b = "         0"
					elseif(j*(Nz-1) < 10)
						b = "         " + num2istr(j*(Nz-1))
					elseif(j*(Nz-1) < 100)
						b = "        " + num2istr(j*(Nz-1))
					elseif(j*(Nz-1) < 1000)
						b = "       " + num2istr(j*(Nz-1))
					elseif(j*(Nz-1) < 10000)
						b = "      " + num2istr(j*(Nz-1))
					elseif(j*(Nz-1) < 100000)
						b = "     " + num2istr(j*(Nz-1))
					endif
									
					for(l = 0; l <= Nz-1; l += 1)
					
						if(l == 0)
							c = "         0"
						elseif(l*(Nx-1) < 10)
							c = "         " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 100)
							c = "        " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 1000)
							c = "       " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 10000)
							c = "      " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 100000)
							c = "     " + num2istr(l*(Nx-1))
						endif
					
						if(cAxis == 1)
							HCklist_tmp[k] = "          " + c + a + b + d + "  1.0"
						elseif(cAxis == 2)
							HCklist_tmp[k] = "          " + a + c + b + d + "  1.0"
						elseif(cAxis == 3)
							HCklist_tmp[k] = "          " + a + b + c + d + "  1.0"
						endif
					
						k += 1
					endfor
				endfor
			endfor
			
			break
		case 2: // orthorhombic (a != b != c)

			If(caxis == 1)
				Nx = floor((N2+1) /2)
				Ny = floor((N3+1) /2)
				Nz = floor((N1+1) /2)
			elseif(caxis ==2)
				Nx = floor((N1+1) /2)
				Ny = floor((N3+1) /2)
				Nz = floor((N2+1) /2)
			elseif(caxis ==3)
				Nx = floor((N1+1) /2)
				Ny = floor((N2+1) /2)
				Nz = floor((N3+1) /2)
			endif
		
			Make /O /T /N=(Round(Nx*Ny*Nz)) HCklist_tmp
			
			Print "Number of non-equivalent k-points in BZ = ", num2istr(DimSize(HCklist_tmp,0))
			
			if((Nx-1)*(Ny-1)*(Nz-1) < 10)
				d = "         " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 100)
				d = "        " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 1000)
				d = "       " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 10000)
				d = "      " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 100000)
				d = "     " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 1000000)
				d = "    " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 10000000)
				d = "   " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			else
				Print "Error: kmesh this fine is not supported."
				return(0)
			endif
			
			k = 0
			
			for(i = 0; i <= Nx-1; i += 1)
			
				if(i == 0)
					a = "         0"
				elseif(i*(Ny-1)*(Nz-1)*xextend < 10)
					a = "         " + num2istr(i*(Ny-1)*(Nz-1)*xextend)
				elseif(i*(Ny-1)*(Nz-1)*xextend < 100)
					a = "        " + num2istr(i*(Ny-1)*(Nz-1)*xextend)
				elseif(i*(Ny-1)*(Nz-1)*xextend < 1000)
					a = "       " + num2istr(i*(Ny-1)*(Nz-1)*xextend)
				elseif(i*(Ny-1)*(Nz-1)*xextend < 10000)
					a = "      " + num2istr(i*(Ny-1)*(Nz-1)*xextend)
				elseif(i*(Ny-1)*(Nz-1)*xextend < 100000)
					a = "     " + num2istr(i*(Ny-1)*(Nz-1)*xextend)
				elseif(i*(Ny-1)*(Nz-1)*xextend < 1000000)
					a = "    " + num2istr(i*(Ny-1)*(Nz-1)*xextend)
				elseif(i*(Ny-1)*(Nz-1)*xextend < 10000000)
					a = "   " + num2istr(i*(Ny-1)*(Nz-1)*xextend)
				endif
				
				for(j = 0; j <= Ny-1; j += 1)
				
					if(j == 0)
						b = "         0"
					elseif(j*(Nx-1)*(Nz-1)*yextend < 10)
						b = "         " + num2istr(j*(Nx-1)*(Nz-1)*yextend)
					elseif(j*(Nx-1)*(Nz-1)*yextend < 100)
						b = "        " + num2istr(j*(Nx-1)*(Nz-1)*yextend)
					elseif(j*(Nx-1)*(Nz-1)*yextend < 1000)
						b = "       " + num2istr(j*(Nx-1)*(Nz-1)*yextend)
					elseif(j*(Nx-1)*(Nz-1)*yextend < 10000)
						b = "      " + num2istr(j*(Nx-1)*(Nz-1)*yextend)
					elseif(j*(Nx-1)*(Nz-1)*yextend < 100000)
						b = "     " + num2istr(j*(Nx-1)*(Nz-1)*yextend)
					elseif(j*(Nx-1)*(Nz-1)*yextend < 1000000)
						b = "    " + num2istr(j*(Nx-1)*(Nz-1)*yextend)
					elseif(j*(Nx-1)*(Nz-1)*yextend < 10000000)
						b = "   " + num2istr(j*(Nx-1)*(Nz-1)*yextend)
					endif
									
					for(l = 0; l <= Nz-1; l += 1)
					
						if(l == 0)
							c = "         0"
						elseif(l*(Nx-1)*(Ny-1) < 10)
							c = "         " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 100)
							c = "        " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 1000)
							c = "       " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 10000)
							c = "      " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 100000)
							c = "     " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 1000000)
							c = "    " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 10000000)
							c = "   " + num2istr(l*(Nx-1)*(Ny-1))
						endif
					
						if(cAxis == 1)
							HCklist_tmp[k] = "          " + c + a + b + d + "  1.0"
						elseif(cAxis == 2)
							HCklist_tmp[k] = "          " + a + c + b + d + "  1.0"
						elseif(cAxis == 3)
							HCklist_tmp[k] = "          " + a + b + c + d + "  1.0"
						endif
					
						k += 1
					endfor
				endfor
			endfor
			
			break
	endswitch
	
	//////////////////////////////////
	ElseIf(primitiveBZ == 2)     // primitive BZ = c2x2 of conventional BZ in xy-plane
	//////////////////////////////////
		switch(latticeType)
		case 0: // cubic (a = b = c)
			Nx = floor((N1+1) /2)
			Ny = Nx
			Nz = Nx
		
			Make /O /T /N=(Round(Nx*(Nx+1)*(Nx+2)/6)) HCklist_tmp
			
			Print "Number of non-equivalent k-points in BZ = ", num2istr(DimSize(HCklist_tmp,0))
			
			if((Nx-1) < 10)
				d = "   " + num2istr((Nx-1)*10*2)
			elseif((Nx-1) < 100)
				d = "  " + num2istr((Nx-1)*10*2)
			elseif((Nx-1) < 499)
				d = " " + num2istr((Nx-1)*10*2)
			else
				Print "Error: kmesh this fine is not supported."
				return(0)
			endif
			
			k = 0
			
			for(i = 0; i <= Nx-1; i += 1)
				for(j = 0; j <= i; j += 1)
					// Transform into the conventional Brilluoin zone coordinates
					ip = i - j
					jp = i + j
				
					if(ip == 0)
						a = "    0"
					elseif(ip < 10)
						a = "   " + num2istr(ip*10)
					elseif(ip < 100)
						a = "  " + num2istr(ip*10)
					else
						a = " " + num2istr(ip*10)
					endif
				
					if(jp == 0)
						b = "    0"
					elseif(jp < 10)
						b = "   " + num2istr(jp*10)
					elseif(jp < 100)
						b = "  " + num2istr(jp*10)
					else
						b = " " + num2istr(jp*10)
					endif
									
					for(l = 0; l <= Nz-1; l += 1)
					
						if(l == 0)
							c = "    0"
						elseif(l < 10)
							c = "   " + num2istr(l*10)
						elseif(l < 100)
							c = "  " + num2istr(l*10)
						else
							c = " " + num2istr(l*10)
						endif
					
						if(cAxis == 1)
							HCklist_tmp[k] = "          " + c + a + b + d + "  1.0"
						elseif(cAxis == 2)
							HCklist_tmp[k] = "          " + a + c + b + d + "  1.0"
						elseif(cAxis == 3)
							HCklist_tmp[k] = "          " + a + b + c + d + "  1.0"
						endif
					
						k += 1
					endfor
				endfor
			endfor
			
			break
		case 1: // tetragonal (a = b != c)
			
			If(caxis == 1)
				Nx = floor((N2+1) /2)
				Ny = Nx
				Nz = floor((N1+1) /2)
			elseif(caxis ==2)
				Nx = floor((N1+1) /2)
				Ny = Nx
				Nz = floor((N2+1) /2)
			elseif(caxis ==3)
				Nx = floor((N1+1) /2)
				Ny = Nx
				Nz = floor((N3+1) /2)
			endif
		
			Make /O /T /N=(Round(Nz*Nx*(Nx+1)/2)) HCklist_tmp
			
			Print "Number of non-equivalent k-points in BZ = ", num2istr(DimSize(HCklist_tmp,0))
			
			if((Nx-1)*(Nz-1) < 10)
				d = "         " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 100)
				d = "        " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 1000)
				d = "       " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 10000)
				d = "      " + num2istr((Nx-1)*(Nz-1)*2)
			elseif((Nx-1)*(Nz-1) < 49999)
				d = "     " + num2istr((Nx-1)*(Nz-1)*2)
			else
				Print "Error: kmesh this fine is not supported."
				return(0)
			endif
			
			k = 0
			
			for(i = 0; i <= Nx-1; i += 1)
				for(j = 0; j <= i; j += 1)
					// Transform into the conventional Brilluoin zone coordinates
					ip = i - j
					jp = i + j
					
					if(ip == 0)
						a = "         0"
					elseif(ip*(Nz-1) < 10)
						a = "         " + num2istr(ip*(Nz-1))
					elseif(ip*(Nz-1) < 100)
						a = "        " + num2istr(ip*(Nz-1))
					elseif(ip*(Nz-1) < 1000)
						a = "       " + num2istr(ip*(Nz-1))
					elseif(ip*(Nz-1) < 10000)
						a = "      " + num2istr(ip*(Nz-1))
					elseif(ip*(Nz-1) < 100000)
						a = "     " + num2istr(ip*(Nz-1))
					endif
						
					if(jp == 0)
						b = "         0"
					elseif(jp*(Nz-1) < 10)
						b = "         " + num2istr(jp*(Nz-1))
					elseif(jp*(Nz-1) < 100)
						b = "        " + num2istr(jp*(Nz-1))
					elseif(jp*(Nz-1) < 1000)
						b = "       " + num2istr(jp*(Nz-1))
					elseif(jp*(Nz-1) < 10000)
						b = "      " + num2istr(jp*(Nz-1))
					elseif(jp*(Nz-1) < 100000)
						b = "     " + num2istr(jp*(Nz-1))
					endif
									
					for(l = 0; l <= Nz-1; l += 1)
					
						if(l == 0)
							c = "         0"
						elseif(l*(Nx-1) < 10)
							c = "         " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 100)
							c = "        " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 1000)
							c = "       " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 10000)
							c = "      " + num2istr(l*(Nx-1))
						elseif(l*(Nx-1) < 100000)
							c = "     " + num2istr(l*(Nx-1))
						endif
					
						if(cAxis == 1)
							HCklist_tmp[k] = "          " + c + a + b + d + "  1.0"
						elseif(cAxis == 2)
							HCklist_tmp[k] = "          " + a + c + b + d + "  1.0"
						elseif(cAxis == 3)
							HCklist_tmp[k] = "          " + a + b + c + d + "  1.0"
						endif
					
						k += 1
					endfor
				endfor
			endfor
			
			break
		case 2: // orthorhombic (a != b != c)
			
			If(caxis == 1)
				Nx = floor((N2+1) /2)
				Ny = floor((N3+1) /2)
				Nz = floor((N1+1) /2)
			elseif(caxis ==2)
				Nx = floor((N1+1) /2)
				Ny = floor((N3+1) /2)
				Nz = floor((N2+1) /2)
			elseif(caxis ==3)
				Nx = floor((N1+1) /2)
				Ny = floor((N2+1) /2)
				Nz = floor((N3+1) /2)
			endif
		
			Make /O /T /N=(Round(Nx*Ny*Nz)) HCklist_tmp
			
			Print "Number of non-equivalent k-points in BZ = ", num2istr(DimSize(HCklist_tmp,0))
			
			if((Nx-1)*(Ny-1)*(Nz-1) < 10)
				d = "         " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 100)
				d = "        " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 1000)
				d = "       " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 10000)
				d = "      " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 100000)
				d = "     " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 1000000)
				d = "    " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			elseif((Nx-1)*(Ny-1)*(Nz-1) < 4999999)
				d = "   " + num2istr((Nx-1)*(Ny-1)*(Nz-1)*2)
			else
				Print "Error: kmesh this fine is not supported."
				return(0)
			endif
			
			k = 0
			
			for(i = 0; i <= Nx-1; i += 1)
				for(j = -i; j <= i; j += 1)
					// Transform into the conventional Brilluoin zone coordinates
					ip = i - j
					jp = i + j
				
					if(ip == 0)
						a = "         0"
					elseif(ip*(Ny-1)*(Nz-1) < 10)
						a = "         " + num2istr(ip*(Ny-1)*(Nz-1))
					elseif(ip*(Ny-1)*(Nz-1) < 100)
						a = "        " + num2istr(ip*(Ny-1)*(Nz-1))
					elseif(ip*(Ny-1)*(Nz-1) < 1000)
						a = "       " + num2istr(ip*(Ny-1)*(Nz-1))
					elseif(ip*(Ny-1)*(Nz-1) < 10000)
						a = "      " + num2istr(ip*(Ny-1)*(Nz-1))
					elseif(ip*(Ny-1)*(Nz-1) < 100000)
						a = "     " + num2istr(ip*(Ny-1)*(Nz-1))
					elseif(ip*(Ny-1)*(Nz-1) < 1000000)
						a = "    " + num2istr(ip*(Ny-1)*(Nz-1))
					elseif(ip*(Ny-1)*(Nz-1) < 10000000)
						a = "   " + num2istr(ip*(Ny-1)*(Nz-1))
					endif
				
					if(jp == 0)
						b = "         0"
					elseif(jp*(Nx-1)*(Nz-1) < 10)
						b = "         " + num2istr(jp*(Nx-1)*(Nz-1))
					elseif(jp*(Nx-1)*(Nz-1) < 100)
						b = "        " + num2istr(jp*(Nx-1)*(Nz-1))
					elseif(jp*(Nx-1)*(Nz-1) < 1000)
						b = "       " + num2istr(jp*(Nx-1)*(Nz-1))
					elseif(jp*(Nx-1)*(Nz-1) < 10000)
						b = "      " + num2istr(jp*(Nx-1)*(Nz-1))
					elseif(jp*(Nx-1)*(Nz-1) < 100000)
						b = "     " + num2istr(jp*(Nx-1)*(Nz-1))
					elseif(jp*(Nx-1)*(Nz-1) < 1000000)
						b = "    " + num2istr(jp*(Nx-1)*(Nz-1))
					elseif(jp*(Nx-1)*(Nz-1) < 10000000)
						b = "   " + num2istr(jp*(Nx-1)*(Nz-1))
					endif
									
					for(l = 0; l <= Nz-1; l += 1)
					
						if(l == 0)
							c = "         0"
						elseif(l*(Nx-1)*(Ny-1) < 10)
							c = "         " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 100)
							c = "        " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 1000)
							c = "       " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 10000)
							c = "      " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 100000)
							c = "     " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 1000000)
							c = "    " + num2istr(l*(Nx-1)*(Ny-1))
						elseif(l*(Nx-1)*(Ny-1) < 10000000)
							c = "   " + num2istr(l*(Nx-1)*(Ny-1))
						endif
					
						if(cAxis == 1)
							HCklist_tmp[k] = "          " + c + a + b + d + "  1.0"
						elseif(cAxis == 2)
							HCklist_tmp[k] = "          " + a + c + b + d + "  1.0"
						elseif(cAxis == 3)
							HCklist_tmp[k] = "          " + a + b + c + d + "  1.0"
						endif
					
						k += 1
					endfor
				endfor
			endfor
			
			break
	endswitch
	///////////////////////////////////
	EndIf
	///////////////////////////////////
	
	
	HCklist_tmp[0] = HCklist_tmp[0] + "-1.00 1.00"
	Redimension /N=(DimSize(HCklist_tmp,0)+1) HCklist_tmp
	HCklist_tmp[numpnts(HCklist_tmp) - 1] = "END"
	
	Save /J /M="\n" HCklist_tmp as "xcrysden.klist"
	
	KillWaves HCklist_tmp

End


//////////////////////////////////////////////////
// Loads a hypercube generated by Wien2k.
//////////////////////////////////////////////////
Function LoadHyperCube(EF, BandMin, BandMax, latticeType, N1, N2, N3, caxis, primitiveBZ)
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
	Make /O /N = (N1, N2, N3, BandNum) WienHyperCube
	
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

	
	// This wave for testing purposes
	Make /O /N = (N1, N2, N3) kmesh_tmp
	SetScale /P x, DimOffset(WienHyperCube,0), DimDelta(WienHyperCube,0), kmesh_tmp
	SetScale /I y, DimOffset(WienHyperCube,1), DimDelta(WienHyperCube,1), kmesh_tmp
	SetScale /I z, DimOffset(WienHyperCube,2), DimDelta(WienHyperCube,2), kmesh_tmp
	
	kmesh_tmp = 0
	
	Open /T="????" /R fileRef
	FReadLine fileRef, line

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

			//kmesh_tmp[kx + floor(DimSize(WienHyperCube, 0)/2) ][ky+floor(DimSize(WienHyperCube, 1)/2)][kz+ floor(DimSize(WienHyperCube, 2)/2)] += 1
			//HC_AssignWithSymm(kmesh_tmp, kx, ky, kz, 0, 1, latticetype, caxis, primitiveBZ)
			
			// All bands that don't have values assigned at this k-point are set to NaN
			For(i = 0; i < bandnum; i+=1)
				HC_AssignWithSymm(WienHyperCube, kx, ky, kz, i, NaN, latticetype, caxis, primitiveBZ)
			EndFor
		
			// read next line
			FReadLine fileRef, line
		
			// loop until you hit end of file or next k-point
			do
		
				// if this line is a correct band, save it's data
				if(str2num(line[0,12]) <= BandMax && str2num(line[0,12]) >= BandMin)
				
					bandindex = Round((str2num(line[0,12]) - DimOffset(WienHyperCube,3))/DimDelta(WienHyperCube,3))
		
					//G1 = floor(DimSize(WienHyperCube, 0)/2)
					//G2 =  floor(DimSize(WienHyperCube, 1)/2)
					//G3 =  floor(DimSize(WienHyperCube, 2)/2)
					//WienHyperCube[G1+kx][G2-ky][G3-kz][bandindex] = (str2num(line[13,40]))
		
					HC_AssignWithSymm(WienHyperCube, kx, ky, kz, bandindex, (str2num(line[13,40])), latticetype, caxis, primitiveBZ)
					
				endif
		
				FReadLine fileRef, line
			
			while(strlen(line) > 0 && cmpstr(line[2,2],".") != 0)
			
		else
			FReadLine fileRef, line
		endif
		
	while(strlen(line) > 0)
		
	
	// Re-scale energy units (EF = 0, energy in eV)
	WienHyperCube -= EF
	WienHyperCube *= 13.605698

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

	killwaves kmesh_tmp

End


//////////////////////////////////////////////////////
// Function used by LoadHyperCube() to symmetrize the hypercube
//////////////////////////////////////////////////////
Function HC_AssignWithSymm(WienHyperCube, kx, ky, kz, bandindex, Eval, latticetype, caxis, primitiveBZ)
	Wave WienHyperCube
	Variable kx, ky, kz, bandindex, Eval, latticetype, caxis, primitiveBZ
	
	Variable G1, G2, G3

	If(primitiveBZ == 2 && latticetype == 0)
		latticetype = 1	// Cubic symmetry operations don't make sense for primitiveBZ == 2.
	EndIf

	G1 = floor(DimSize(WienHyperCube, 0)/2)
	G2 =  floor(DimSize(WienHyperCube, 1)/2)
	G3 =  floor(DimSize(WienHyperCube, 2)/2)
	
	// Inversion
	If(primitiveBZ == 1 || primitiveBZ == 3)
	WienHyperCube[G1+kx][G2+ky][G3+kz][bandindex] = Eval
	WienHyperCube[G1-kx][G2+ky][G3+kz][bandindex] = Eval
	WienHyperCube[G1+kx][G2-ky][G3+kz][bandindex] = Eval
	WienHyperCube[G1-kx][G2-ky][G3+kz][bandindex] = Eval
	WienHyperCube[G1+kx][G2+ky][G3-kz][bandindex] = Eval
	WienHyperCube[G1-kx][G2+ky][G3-kz][bandindex] = Eval
	WienHyperCube[G1+kx][G2-ky][G3-kz][bandindex] = Eval
	WienHyperCube[G1-kx][G2-ky][G3-kz][bandindex] = Eval
	ElseIf(primitiveBZ == 2)	// Inversion about zone diagonals in kx,ky
		If(caxis == 1)
			WienHyperCube[G1+kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2-kz][G3-ky][bandindex] = Eval
			WienHyperCube[G1+kx][G2+kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1+kx][G2-ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2-kz][G3-ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2+kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3-kz][bandindex] = Eval
		Elseif(caxis == 2)
			WienHyperCube[G1+kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kz][G2+ky][G3-kx][bandindex] = Eval
			WienHyperCube[G1+kz][G2+ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1-kx][G2+ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2-ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kz][G2-ky][G3-kx][bandindex] = Eval
			WienHyperCube[G1+kz][G2-ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3-kz][bandindex] = Eval
		Else
			WienHyperCube[G1+kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2-kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1+ky][G2+kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2+ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2-kx][G3-kz][bandindex] = Eval
			WienHyperCube[G1+ky][G2+kx][G3-kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3-kz][bandindex] = Eval
		EndIf
	EndIf
	
	
	If(latticeType == 0) // cubic
		// x <-> y
		WienHyperCube[G1+ky][G2+kx][G3+kz][bandindex] = Eval
		WienHyperCube[G1-ky][G2+kx][G3+kz][bandindex] = Eval
		WienHyperCube[G1+ky][G2-kx][G3+kz][bandindex] = Eval
		WienHyperCube[G1-ky][G2-kx][G3+kz][bandindex] = Eval
		WienHyperCube[G1+ky][G2+kx][G3-kz][bandindex] = Eval
		WienHyperCube[G1-ky][G2+kx][G3-kz][bandindex] = Eval
		WienHyperCube[G1+ky][G2-kx][G3-kz][bandindex] = Eval
		WienHyperCube[G1-ky][G2-kx][G3-kz][bandindex] = Eval
		// x <-> z
		WienHyperCube[G1+kz][G2+ky][G3+kx][bandindex] = Eval
		WienHyperCube[G1-kz][G2+ky][G3+kx][bandindex] = Eval
		WienHyperCube[G1+kz][G2-ky][G3+kx][bandindex] = Eval
		WienHyperCube[G1-kz][G2-ky][G3+kx][bandindex] = Eval
		WienHyperCube[G1+kz][G2+ky][G3-kx][bandindex] = Eval
		WienHyperCube[G1-kz][G2+ky][G3-kx][bandindex] = Eval
		WienHyperCube[G1+kz][G2-ky][G3-kx][bandindex] = Eval
		WienHyperCube[G1-kz][G2-ky][G3-kx][bandindex] = Eval
		// y <-> z
		WienHyperCube[G1+kx][G2+kz][G3+ky][bandindex] = Eval
		WienHyperCube[G1-kx][G2+kz][G3+ky][bandindex] = Eval
		WienHyperCube[G1+kx][G2-kz][G3+ky][bandindex] = Eval
		WienHyperCube[G1-kx][G2-kz][G3+ky][bandindex] = Eval
		WienHyperCube[G1+kx][G2+kz][G3-ky][bandindex] = Eval
		WienHyperCube[G1-kx][G2+kz][G3-ky][bandindex] = Eval
		WienHyperCube[G1+kx][G2-kz][G3-ky][bandindex] = Eval
		WienHyperCube[G1-kx][G2-kz][G3-ky][bandindex] = Eval
		// x -> y -> z
		WienHyperCube[G1+kz][G2+kx][G3+ky][bandindex] = Eval
		WienHyperCube[G1-kz][G2+kx][G3+ky][bandindex] = Eval
		WienHyperCube[G1+kz][G2-kx][G3+ky][bandindex] = Eval
		WienHyperCube[G1-kz][G2-kx][G3+ky][bandindex] = Eval
		WienHyperCube[G1+kz][G2+kx][G3-ky][bandindex] = Eval
		WienHyperCube[G1-kz][G2+kx][G3-ky][bandindex] = Eval
		WienHyperCube[G1+kz][G2-kx][G3-ky][bandindex] = Eval
		WienHyperCube[G1-kz][G2-kx][G3-ky][bandindex] = Eval
		// x -> z -> y
		WienHyperCube[G1+ky][G2+kz][G3+kx][bandindex] = Eval
		WienHyperCube[G1-ky][G2+kz][G3+kx][bandindex] = Eval
		WienHyperCube[G1+ky][G2-kz][G3+kx][bandindex] = Eval
		WienHyperCube[G1-ky][G2-kz][G3+kx][bandindex] = Eval
		WienHyperCube[G1+ky][G2+kz][G3-kx][bandindex] = Eval
		WienHyperCube[G1-ky][G2+kz][G3-kx][bandindex] = Eval
		WienHyperCube[G1+ky][G2-kz][G3-kx][bandindex] = Eval
		WienHyperCube[G1-ky][G2-kz][G3-kx][bandindex] = Eval
	
	ElseIf(latticeType == 1) // tetragonal
		If(primitiveBZ == 1)
		If(caxis == 1)
			// y <-> z
			WienHyperCube[G1+kx][G2+kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2+kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1+kx][G2-kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2-kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1+kx][G2+kz][G3-ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2+kz][G3-ky][bandindex] = Eval
			WienHyperCube[G1+kx][G2-kz][G3-ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2-kz][G3-ky][bandindex] = Eval
		Elseif(caxis ==2)
			// x <-> z
			WienHyperCube[G1+kz][G2+ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1-kz][G2+ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1+kz][G2-ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1-kz][G2-ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1+kz][G2+ky][G3-kx][bandindex] = Eval
			WienHyperCube[G1-kz][G2+ky][G3-kx][bandindex] = Eval
			WienHyperCube[G1+kz][G2-ky][G3-kx][bandindex] = Eval
			WienHyperCube[G1-kz][G2-ky][G3-kx][bandindex] = Eval
		Else
			// x <-> y
			WienHyperCube[G1+ky][G2+kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2+kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1+ky][G2-kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2-kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1+ky][G2+kx][G3-kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2+kx][G3-kz][bandindex] = Eval
			WienHyperCube[G1+ky][G2-kx][G3-kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2-kx][G3-kz][bandindex] = Eval
		EndIf
		
		ElseIf(primitiveBZ == 2)
		If(caxis == 1)
			// y <-> z
			WienHyperCube[G1+kx][G2+kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1+kx][G2-ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2-kz][G3-ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2+kz][G3+ky][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2-kz][G3-ky][bandindex] = Eval
		Elseif(caxis==2)
			// x <-> z
			WienHyperCube[G1+kz][G2+ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1-kx][G2+ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kz][G2+ky][G3-kx][bandindex] = Eval
			WienHyperCube[G1+kz][G2-ky][G3+kx][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2-ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kz][G2-ky][G3-kx][bandindex] = Eval
		Else
			// x <-> y
			WienHyperCube[G1+ky][G2+kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2+ky][G3+kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2-kx][G3+kz][bandindex] = Eval
			WienHyperCube[G1+ky][G2+kx][G3-kz][bandindex] = Eval
			WienHyperCube[G1-kx][G2-ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1+kx][G2+ky][G3-kz][bandindex] = Eval
			WienHyperCube[G1-ky][G2-kx][G3-kz][bandindex] = Eval
		EndIf
		
		EndIf
	
	Elseif(latticeType == 2) // orthorhombic
		// No more symmetries
	EndIf

End





//////////////////////////////////
// Finds the bands that pass through the energy specified by 'EF' (EF = 0 for the Fermi level)

// Note that currently this doesn't correctly handle waves that are entirely NaN (ie: If you load more bands than
// were output by your calculation, it may tell you that these empty bands cross the Fermi level).
//////////////////////////////////
Function HyperCubeBandsAtEF(WienHyperCube, EF)
	Wave WienHyperCube
	Variable EF

	Variable l, count
	
	Make /O /N = (DimSize(WienHyperCube, 3)) BandsAtEF_List
	BandsAtEF_List = 0

	l = 0
	count = 0
	do
		Make /O /N = (DimSize(WienHyperCube, 0), DimSize(WienHyperCube, 1), DimSize(WienHyperCube, 2)) BandsAtEF_Cube

		BandsAtEF_Cube = WienHyperCube[p][q][r][l]
	
		if(WaveMin(BandsAtEF_Cube) <= EF && WaveMax(BandsAtEF_Cube) >= EF)
			BandsAtEF_List[count] = l+1
			count += 1 
			Print l+1
		endif
		
		l+=1
	while(l < DimSize(WienHyperCube, 3))

	Redimension /N = (count) BandsAtEF_List
	KillWaves BandsAtEF_Cube

End




//////////////////////////////////
// This function pulls 3D Fermi-Surface data from a HyperCube
// The FS files can be plotted with Gizmo using a contour plot at '0' (or a different energy for a rigid-shift of EF)
// Use the BandsAtEF function to find bands that cut through the Fermi level
//////////////////////////////////
Function HyperCube3DFS(WienHyperCube, Bandlist)
	Wave WienHyperCube	// The hypercube to draw data from
	Wave Bandlist			// A wave with the band indicies to draw FS's from (ie: "BandsAtEF_List")
	
	Variable i
	
	For(i = 0; i < DimSize(bandlist, 0) ; i+=1)
	
		Make /O /N = (DimSize(WienHyperCube, 0),DimSize(WienHyperCube, 1), DimSize(WienHyperCube, 2)) $("FS_Band" + num2str(bandlist[i]))
	
		Wave FS = $("FS_Band" + num2str(bandlist[i]))
		
		SetScale /P x, (DimOffset(WienHyperCube, 0)), (DimDelta(WienHyperCube, 0)), FS
		SetScale /P y, (DimOffset(WienHyperCube, 1)), (DimDelta(WienHyperCube, 1)), FS
		SetScale /P z, (DimOffset(WienHyperCube, 2)), (DimDelta(WienHyperCube, 2)), FS	
		
		
		FS = WienHyperCube[p][q][r][Bandlist[i] - 1]
		
	EndFor
End



//////////////////////////////////
// This function pulls 2D Fermi-Surface data from a HyperCube (at fixed kz)
// The FS files can be plotted with a contour plot at '0' (or a different energy for a rigid-shift of EF)
// Use the BandsAtEF function to find bands that cut through the Fermi level
//////////////////////////////////
Function HyperCube2DFS(WienHyperCube, Bandlist, kz)
	Wave WienHyperCube	// The hypercube to draw data from
	Wave Bandlist			// A wave with the band indicies to draw FS's from (ie: "BandsAtEF_List")
	Variable kz				// kz value for cut. In units of pi/c and in first BZ.
	
	Variable kzName = kz
	
	Variable i
	
	Make /O /N = (DimSize(WienHyperCube, 0), DimSize(WienHyperCube, 1), DimSize(WienHyperCube, 2)) HyperCubeFS_Temp
	SetScale /P z, DimOffset(WienHyperCube, 2), DimDelta(WienHyperCube,2), HyperCubeFS_Temp
	
	If(kz == 1)
		kz = -1
	EndIf
	
	
	For(i = 0; i < DimSize(bandlist, 0) ; i+=1)
	
		Make /O /N = (DimSize(WienHyperCube, 0),DimSize(WienHyperCube, 1)) $("FS_Band" + num2str(bandlist[i]) + "_kz" + num2str(kzName))
	
		Wave FS = $("FS_Band" + num2str(bandlist[i]) + "_kz" + num2str(kzName))
		//SetScale /I x, -1, 1, FS
		//SetScale /I y, -1, 1, FS
		SetScale /P x, (DimOffset(WienHyperCube, 0)), (DimDelta(WienHyperCube, 0)), FS
		SetScale /P y, (DimOffset(WienHyperCube, 1)), (DimDelta(WienHyperCube, 1)), FS	
		
		//If(floor((DimSize(WienHyperCube, 0)-1)/2) != (DimSize(WienHyperCube, 0)-1)/2)
		//	SetScale /I x, (-1 + 1/(2*(DimSize(WienHyperCube, 0)))), (1 - 1/(2*(DimSize(WienHyperCube, 0)))), FS
		//	SetScale /I y, (-1 + 1/(2*(DimSize(WienHyperCube, 0)))), (1 - 1/(2*(DimSize(WienHyperCube, 0)))), FS
		//EndIf
		

		HyperCubeFS_Temp = WienHyperCube[p][q][r][bandlist[i] - 1]
		
		//FS[][] =  interp3D(HyperCubeFS_Temp, p-0.0001, q-0.0001, kz)
		FS[][] =  interp3D(HyperCubeFS_Temp, p, q, kz)
		
	EndFor
	
	KillWaves HyperCubeFS_Temp
	
End



////////////////////////////////////////
// This function takes a kz projected (and smeared) Fermi-surface from a 3D Fermi Surface
////////////////////////////////////////
Function HyperCubeProjectFS(EF, Ewindow, kz_center, kz_width, upsample, FS3D_in)
	Variable EF, Ewindow				// The Fermi level (eV) and the half-width (meV) to use for integrating
	Variable kz_center, kz_width		// Both in pi/a, kz_width is the FWHM of a Lorentzian
	Variable upsample				// Up-sample the Fermi surface in k-space by this integer factor (1==no upsampling)
	Wave FS3D_in					// A 3D Fermi surface to project out

	Variable i, kz, kz_start, kz_num
	Variable Nwidth = 10				// Sets the number of FWHM to truncate the integration at

	// Up-sample the Fermi Surface
	Duplicate/O $(NameOfWave(FS3D_in)),$(NameOfWave(FS3D_in) + "_samp")
	Wave FS3D = $(NameOfWave(FS3D_in) + "_samp")
	upsample = Round(upsample)
	If(upsample != 1)
		Resample/DIM=0/UP=(upsample) FS3D
		Resample/DIM=1/UP=(upsample) FS3D
		Resample/DIM=2/UP=(upsample) FS3D
	EndIf

	// Make 2D wave to store FS
	Make /O /N = (DimSize(FS3D, 0), DimSize(FS3D, 1)) $(NameOfWave(FS3D_in) + "_2D")
	Wave FS2D = $(NameOfWave(FS3D_in) + "_2D")
	Make /O /N = (DimSize(FS3D, 0), DimSize(FS3D, 1)) FS2D_temp
	SetScale /I x, -1, 1, FS2D
	SetScale /I y, -1, 1, FS2D
	SetScale /I x, -1, 1, FS2D_temp
	SetScale /I y, -1, 1, FS2D_temp
	FS2D = 0

	// Make a wave that defines the FS window
	Make /O /N=2 HCProjectFS_twave
	HCProjectFS_twave[0] = EF - Ewindow/1000
	HCProjectFS_twave[1] = EF + Ewindow/1000

	// Convert the kz center and window into a wave that contains z-indicies for FS3D and weights
		kz_num = Ceil(Nwidth*kz_width * ((DimSize(FS3D,2)+1)/2))
		kz_start = kz_center - Nwidth*kz_width/2

		// Get kz_start into first BZ
		kz_start -=2
		do
			kz_start += 2
		while(kz_start < -1)
		kz_start +=2
		do
			kz_start -= 2
		while(kz_start > 1)

		Make /O /N=(kz_num,2) HCProjectFS_zwave
		HCProjectFS_zwave = 0
	
		//Convert kz_start to a point number and assign its weight
		HCProjectFS_zwave[0][0] = Round((kz_start - DimOffset(FS3D, 2))/DimDelta(FS3D,2))
		HCProjectFS_zwave[0][1] = ((kz_width/2)^2/((kz_start-kz_center)^2 + (kz_width/2)^2))/Pi
	

	For(i = 1; i < kz_num; i +=1)
		//Convert kz to a point number and assign its weight
		HCProjectFS_zwave[i][0] = Round(mod((i+HCProjectFS_zwave[0][0]),DimSize(FS3D,2)))
		HCProjectFS_zwave[i][1] = ((kz_num/(2*Nwidth))^2/((i-kz_num/2)^2 + (kz_num/(2*Nwidth))^2))/Pi
	EndFor
	
	
	// Loop over all kz
	For(i = 0; i < kz_num; i +=1)
		
		// Pull the 2DFS from the 3DFS at this kz
		FS2D_temp[][] = FS3D[p][q][HCProjectFS_zwave[i][0]]
		
		// Convert this slice of hypercube into a binary array, and then rescale to floating point (pulls out values within Ewindow of EF)
		ImageThreshold /O /W=HCProjectFS_twave FS2D_temp
		FS2D_temp /= 255
		Redimension /S FS2D_temp
			
		// Weight it by a lorentzian in kz and then add it to the FS waves
		FS2D_temp *= HCProjectFS_zwave[i][1]
		
		FS2D += FS2D_temp
	
	EndFor	
	Note FS2D "Generated by HyperCubeProjectFS on " + Date() + ", " + Time()
	Note FS2D "With (EF, Ewindow, kz_center, kz_width, upsample, FS3D_in) = (" + num2str(EF) + ", " + num2str(Ewindow) + ", " + num2str(kz_center) + ", " + num2str(kz_width) + ", " + num2istr(upsample) + ", " + NameOfWave(FS3D_in) + ")"

	KillWaves FS2D_temp, HCProjectFS_twave, HCProjectFS_zwave, FS3D

End



/////////////////////////////////////
// Generates a Klist along an ARPES cut for use in spaghetti plots and ARPES spectra simulation
// Output klist is scaled so that its projection onto a straight line is evenly spaced

// Uses a function from the Blue package (BlueZoneAngleToMomentum)
/////////////////////////////////////
Function HCGenerateKlist(klist_name, KE, theta, phi, omega, orientation, NumofPnts, thetaRange)
	Variable KE			// electron kinetic energy (eV)
	Variable theta		// sample manipulator theta angle (degrees)
	Variable phi			// sample manipulator phi angle (degrees)
	Variable omega		// sample manipulator omega angle (degrees)
	Variable orientation	// analyzer orientation (1 = vertical angle dispersion, 2 = horizontal angle dispersion)
	
	String klist_name		// Name of the klist to output
	Variable NumofPnts	// Number of points in the klist
	Variable thetaRange	// Range in theta that cut will span (+/- thetarange/2)

	Variable i, j, ux, uy, a, b

	Make /O /N = (NumofPnts, 2), $(klist_name)
	Wave klist =  $(klist_name)
	
	Make /O /N = (NumofPnts, 2), HCklist_tmp
	SetScale /I x, (-thetaRange/2), (thetaRange/2), HCklist_tmp
	
	// Generate points along the theta/phi cut, with constant spacing in theta
	HCklist_tmp[][0] = real(BlueZoneAngleToMomentum(KE,theta,phi,omega,x,orientation))
	HCklist_tmp[][1] = imag(BlueZoneAngleToMomentum(KE,theta,phi,omega,x,orientation))

	//Convert into a klist that is evenly spaced in terms of its projection onto a straight line
	ux = HCklist_tmp[DimSize(HCklist_tmp, 0)-1][0] - HCklist_tmp[0][0]
	uy = HCklist_tmp[DimSize(HCklist_tmp, 0)-1][1] - HCklist_tmp[0][1] 
	ux /= sqrt(ux^2 + uy^2)
	uy = sqrt(1 - ux^2)
	
	Make /O /N = (DimSize(HCklist_tmp, 0), 2) HCktmp
	HCktmp[][0] = ux*(HCklist_tmp[p][0]) + uy*(HCklist_tmp[p][1])
	HCktmp[][1] =  HCktmp[0][0] + (HCktmp[DimSize(HCktmp,0)-1][0]- HCktmp[0][0])*p/(DimSize(HCktmp, 0)-1)
	
	klist = NaN
	klist[0][] = HCklist_tmp[0][q]
	j = 1
	For(i = 1; i < DimSize(HCklist_tmp, 0); i+=1)
		If(HCktmp[i][0] >= HCktmp[j][1])
			a = (HCktmp[i][0]-HCktmp[j][1])/(HCktmp[i][0] -  HCktmp[i-1][0])
			b = (HCktmp[j][1]-HCktmp[i-1][0])/(HCktmp[i][0] -  HCktmp[i-1][0])
			klist[j][] = (a*HCklist_tmp[i-1][q] + b*HCklist_tmp[i][q])
			j+=1
			If(HCktmp[i][0] >= HCktmp[j][1] && j < DimSize(HCktmp, 0))
				i -= 1
			EndIf
		EndIf
	EndFor
	
	Note klist, "Klist Parameters Used:"
	Note klist, "KE=" + num2str(KE)
	Note klist, "theta=" + num2str(theta)
	Note klist, "phi=" + num2str(phi)
	Note klist, "omega=" + num2str(omega)
	Note klist, "orientation=" + num2str(orientation)

	KillWaves HCktmp, HCklist_tmp

End



///////////////////////////////////
// Extracts spaghetti-plot data from the Hypercube along a cut at a given kz
// Compatible with the k-list exported by BlueZone (but do not use 'Rotate zone with omega')
//
// kList should be an Nx2 wave, with the first column kx and the second column ky (see below for units)
// bandlist is a 1D wave listing the band indicies of the bands to extract
//
// Regardless of kList, kz must be in units of pi/c
//
// The output is a 2D wave, with the x-axis scaled from the first to the last k-point supplied (linearly), and
// with the same units as the supplied klist. The first column of this wave is a projected k-list.
// Note: If your k-list contains more than one line through k-space, the first column of the output wave
// will not be meaningful, and you should scale your k-axis accordingly (ie: set Spaghetti[][0] = p or similar)
///////////////////////////////////
Function HyperCubeSlice_kz(WienHyperCube, kList, bandlist, kz, a, b)
	Wave WienHyperCube
	Wave kList				// kList in units of pi/a or inv. Ang. (but see note for a and b below)
	Wave bandlist
	Variable kz				// kz IN UNITS OF pi/c.
	Variable a, b				// Lattice constants (in Angstroms). If b = 0, cubic/tetragonal symmetry is assumed.
							// If a = 0, kList is assumed to already be in units of pi/a
	Variable ux, uy, i
	Variable bandindex
	
	If(kz == 1)
		kz = -1
	EndIf
	
	
	If(trunc(kz)/2 == trunc(trunc(kz)/2))
		kz -= trunc(kz)
	Else
		kz -= trunc(kz) + 1*Sign(kz)
	EndIf
	
	
	Make /O /N = (DimSize(WienHyperCube, 0), DimSize(WienHyperCube, 1), DimSize(WienHyperCube, 2)) HyperCubeSlice_Temp
	SetScale /P x, DimOffset(WienHyperCube, 0), DimDelta(WienHyperCube, 0), HyperCubeSlice_Temp
	SetScale /P y, DimOffset(WienHyperCube, 1), DimDelta(WienHyperCube, 1), HyperCubeSlice_Temp
	SetScale /P z, DimOffset(WienHyperCube, 2), DimDelta(WienHyperCube, 2), HyperCubeSlice_Temp
	
	
	Make /O /N = (DimSize(kList, 0), DimSize(kList, 1)) kList_temp
	Make /O /N = (DimSize(klist, 0), (DimSize(bandlist,0)+1)) Spaghetti
	
	HyperCubeSlice_Temp = 0
	Spaghetti = 0
	kList_temp = 0
	
	// Make a temporary klist in units of pi/a
	If(a == 0)
		kList_temp = kList
	Elseif(b == 0)
		b = a
		kList_temp[][] = kList[p][q]  * a / Pi
	Else
		kList_temp[][0] = kList[p][0]  * a / Pi
		kList_temp[][1] = kList[p][1]  * b / Pi
	EndIf
	
	// Fold all points back into 1st BZ
	For(i = 0; i < DimSize(kList_temp, 0); i+= 1)
	
		If(trunc(kList_temp[i][0])/2 == trunc(trunc(kList_temp[i][0])/2))
			kList_temp[i][0] -= trunc(kList_temp[i][0])
		Else
			kList_temp[i][0] -= trunc(kList_temp[i][0]) + 1*Sign(kList_temp[i][0])
		EndIf
		
		If(trunc(kList_temp[i][1])/2 == trunc(trunc(kList_temp[i][1])/2))
			kList_temp[i][1] -= trunc(kList_temp[i][1])
		Else
			kList_temp[i][1] -= trunc(kList_temp[i][1]) + 1*Sign(kList_temp[i][1])
		EndIf	 
		
	EndFor
	
	
	// Fill the Spaghetti plot with interpolated data
	For(i = 0; i < DimSize(bandlist, 0); i+=1)
		bandindex = bandlist[i]
		HyperCubeSlice_Temp = WienHyperCube[p][q][r][bandindex - 1]
		
		Spaghetti[][i+1] = interp3D(HyperCubeSlice_Temp, kList_temp[p][0], kList_temp[p][1], kz)
	
	EndFor
	
	// Scale the output wave
	ux = kList[DimSize(kList, 0)-1][0] - kList[0][0]
	uy = kList[DimSize(kList, 0)-1][1] - kList[0][1] 
	ux /= sqrt(ux^2 + uy^2)
	uy = sqrt(1 - ux^2)
	
	Spaghetti[][0] = ux*(klist[p][0]) + uy*(klist[p][1])
	
	// Set an approximate k-scale for the wave
	SetScale /I x, (kList[0][0]*ux + kList[0][1]*uy), (kList[DimSize(kList, 0)-1][0]*ux + kList[DimSize(kList, 0)-1][1]*uy), Spaghetti

	//A fix for klists that start and end at the same point
	If(DimDelta(Spaghetti, 0)*0 != 0)
		SetScale /I x, 0, 1, Spaghetti
	EndIf

	KillWaves HyperCubeSlice_Temp, kList_Temp

End


/////////////////////////////////
// Plots a bandstructure (Spaghetti) wave
//
// Note: To overlay a spectrum from Blue onto this plot, you'll need to use MatrixTranspose on the spectrum
/////////////////////////////////
Function HyperCubePlotSpaghetti(Spaghetti, App, Xaxis)
	Wave Spaghetti
	Variable App // 0 == create a new graph, 1 == append to top graph
	Variable Xaxis // 0 == plot against column 1 of the spaghetti wave, 1 == plot against calculated x axis
	Variable i
	
	If(Xaxis == 0)
		If(App == 0)
			Display Spaghetti[][1] vs Spaghetti[][0]
		ElseIf(App == 1)
			AppendToGraph Spaghetti[][1] vs Spaghetti[][0]
		EndIf
	
		For(i=1; i< DimSize(Spaghetti, 1)-1; i+=1)
			AppendToGraph Spaghetti[][i+1] vs Spaghetti[][0]
		EndFor
	ElseIf(Xaxis == 1)
		If(App == 0)
			Display Spaghetti[][1]
		ElseIf(App == 1)
			AppendToGraph Spaghetti[][1]
		EndIf
	
		For(i=1; i< DimSize(Spaghetti, 1)-1; i+=1)
			AppendToGraph Spaghetti[][i+1]
		EndFor
	EndIf

End


//////////////////////////////////////
// This function generates bandstructure over a range in Kz. Output is a 3D wave, with the first 2 dimensions
// giving a 'Spaghetti' plot, and the 3rd dimension indexing kz.
//
// See 'HyperCubeSlice_kz' for notes on input.
//////////////////////////////////////
Function HyperCubeSlice_kzRange(WienHyperCube, kList, bandlist, kz_center, kz_window, kz_number, a, b)
	Wave WienHyperCube
	Wave kList
	Wave bandlist
	Variable a, b
	Variable kz_center		// Center kz value to base the 'smear' off of (in pi/a)
	Variable kz_window		// kz window (full-width) to include (in pi/a)
	Variable kz_number		// number of steps to take in kz over the entire window (in pi/a)

	Variable i, kz
	Variable ux, uy

	Make /O /N = (DimSize(klist, 0), DimSize(bandlist,0)+1, kz_number) Spaghetti_3D
	Make /O /N = 1 Spaghetti
		
	// Scale the output wave
	ux = kList[DimSize(kList, 0)-1][0] - kList[0][0]
	uy = kList[DimSize(kList, 0)-1][1] - kList[0][1] 
	ux /= sqrt(ux^2 + uy^2)
	uy = sqrt(1 - ux^2)
	
	SetScale /I x, (kList[0][0]*ux + kList[0][1]*uy), (kList[DimSize(kList, 0)-1][0]*ux + kList[DimSize(kList, 0)-1][1]*uy), Spaghetti_3D
	
	//A fix for klists that start and end at the same point
	If(DimDelta(Spaghetti_3D, 0)*0 != 0)
		SetScale /I x, 0, 1, Spaghetti_3D
	EndIf

	For(i = 0; i < kz_number; i+=1)
		kz =  (kz_center - kz_window/2) + i * kz_window/(kz_number-1)
		
		HyperCubeSlice_kz(WienHyperCube, kList, bandlist, kz, a, b)
		Spaghetti_3D[][][i] = Spaghetti[p][q]
	EndFor
		
	KillWaves Spaghetti

End



/////////////////////////////////
// Plots a 3D bandstructure (Spaghetti_3D) wave
//
// Note: To overlay a spectrum from Blue onto this plot, you'll need to use MatrixTranspose on the spectrum
/////////////////////////////////
Function HyperCubePlotSpaghetti3D(Spaghetti_3D, EF, App, Xaxis)
	Wave Spaghetti_3D
	Variable EF // Fermi level for plot (set EF = 0 unless you want to shift the chemical potential)
	Variable App // App == 0; make a new graph. App == 1; Append to top graph
	Variable Xaxis // 0 == plot against column 1 of the spaghetti wave, 1 == plot against calculated x axis
	Variable i, j
	
	If(Xaxis == 0)
		If(App == 0)
			Display Spaghetti_3D[][1][0] vs Spaghetti_3D[][0][0]
			RemoveFromGraph $(NameOfWave(Spaghetti_3D))
		EndIf
	
		For(j=0; j< DimSize(Spaghetti_3D, 2); j+=1)
			For(i=0; i< DimSize(Spaghetti_3D, 1); i+=1)
				AppendToGraph Spaghetti_3D[][i+1][j] vs Spaghetti_3D[][0][j]
			EndFor
		EndFor
	ElseIf(Xaxis == 1)
		If(App == 0)
			Display Spaghetti_3D[][1][0]
			RemoveFromGraph $(NameOfWave(Spaghetti_3D))
		EndIf
	
		For(j=0; j< DimSize(Spaghetti_3D, 2); j+=1)
			For(i=0; i< DimSize(Spaghetti_3D, 1); i+=1)
				AppendToGraph Spaghetti_3D[][i+1][j]
			EndFor
		EndFor
	
	EndIf
	
	ModifyGraph offset={0, -EF}
	TextBox/C/N=text0/A=MC NameOfWave(Spaghetti_3D)
	
End





//////////////////////////////////////////////
// This function simulates an ARPES spectrum, smearing in kz and adding a linear-in-energy lifetime.
// 

// For accurate spectra, klist should either be an evenly spaced line or generated by "HCGenerateKlist".
// (The k-axis of the output spectrum is scaled to be evenly spaced along a line connecting the first and last kpoint.
// This introduces no error for the two above options, but for an improperly scaled k-wave can introduce some
// artificial compression or expansion of the spectrum)

// The quality of the spectrum can be very sensitive to 'k_number' for very kz-dispersive features.

// Use BlueBlur to simulate experimental resolution effects

//Notes:
// See PRB 77, 165120 (2008) for a discussion of the broadening being done here
//////////////////////////////////////////////
Function HCsimulateARPES(WienHyperCube, kList, bandlist, Elow, Estep, NumEsteps, EF, kz_center, kz_width, kz_number, Gam, a, b)
	Variable Elow, Estep, NumEsteps		// Specify the lowest energy in the spectrum, the energy step size, and the total number of steps in energy (all eV)
	Variable EF							// The Fermi-level of WienHyperCube (usually 0)
	Variable kz_center					// The kz-value of the nominal direct-transition in units of pi/c
	Variable kz_width						// The FWHM of the Lorentzian broadening in kz (= 1/escape depth) in pi/c
	Variable kz_number					// The total number of points to sample in the kz direction over kz_range
	Variable Gam						// The coefficient for linear lifetime broadening in Energy of the DFT bands (width = Gam*(EF-E))
	Variable a, b							// Lattice constants (in Angstroms). If b = 0, cubic/tetragonal symmetry is assumed (b=a).
										// If a = 0, kList is assumed to already be in units of pi/a 

	Wave WienHyperCube, kList, bandlist	// See "HyperCubeSlice_kz" notes	
	
	Variable kz_range = 10*kz_width 	// The total range in kz to integrate over (Full width)


	Variable i, j, ux, uy
	Variable Normalize

	// A wave used to store a single DFT band at a time
	Make /O /N = (DimSize(klist, 0),kz_number) epsilon
	
	// Energy scale
	Make /O /N = (NumEsteps) Energy
	Energy = Elow + Estep * x
	Variable EF_pntnum = Floor((EF-Elow)/Estep)		// Determines the last point in my energy wave below EF
	
	//kz scale
	Make /O /N =(kz_number) kz
	kz = kz_center - kz_range/2 + x*(kz_range/(kz_number-1))

	// Generate bandstructure over a wide range in kz (relative to delta kz)
	//HyperCubeSlice_kzRange(WienHyperCube, kList, bandlist, kz_center, kz_range, kz_number, a, b)
	//Wave Spaghetti_3D
	
	/////////////////////////////////
	// Generate the list of kpoints (collapsed to first BZ) that will be used for the calculation (in pi/a)
	Make /O /N = (kz_number*DimSize(klist,0),3) MasterKlist
	MasterKlist = 0
	
	For(i = 0; i < kz_number; i+=1)
		If(a == 0)
			MasterKlist[i*DimSize(klist,0), (i+1)*DimSize(klist,0)-1][0] = klist[p-i*DimSize(klist,0)][0]
			MasterKlist[i*DimSize(klist,0), (i+1)*DimSize(klist,0)-1][1] = klist[p-i*DimSize(klist,0)][1]
		ElseIf(b==0)
			MasterKlist[i*DimSize(klist,0), (i+1)*DimSize(klist,0)-1][0] = klist[p-i*DimSize(klist,0)][0] * a / Pi
			MasterKlist[i*DimSize(klist,0), (i+1)*DimSize(klist,0)-1][1] = klist[p-i*DimSize(klist,0)][1] * a / Pi
		Else
			MasterKlist[i*DimSize(klist,0), (i+1)*DimSize(klist,0)-1][0] = klist[p-i*DimSize(klist,0)][0] * a / Pi
			MasterKlist[i*DimSize(klist,0), (i+1)*DimSize(klist,0)-1][1] = klist[p-i*DimSize(klist,0)][1] * b / Pi
		EndIf
		MasterKlist[i*DimSize(klist,0), (i+1)*DimSize(klist,0)-1][2] = kz[i]
	EndFor
	
	MasterKlist /= 2 // Convert into 2 Pi/a
	MasterKlist += 1/2  // Center my coord system at 'M'
	
	// Collapse into 1st BZ
	MasterKlist[][][] = Mod(MasterKlist[p][q][r],1)	
	MasterKlist += 1
	MasterKlist[][][] = Mod(MasterKlist[p][q][r],1)	
	
	MasterKlist -= 1/2	// Shift center back to Gamma
	MasterKlist *= 2 // Convert back into Pi/a
	
	Make /O /N = (DimSize(WienHyperCube, 0), DimSize(WienHyperCube, 1), DimSize(WienHyperCube, 2)) HCxtractBnd
	SetScale /P x, DimOffset(WienHyperCube, 0), DimDelta(WienHyperCube, 0), HCxtractBnd
	SetScale /P y, DimOffset(WienHyperCube, 1), DimDelta(WienHyperCube, 1), HCxtractBnd
	SetScale /P z, DimOffset(WienHyperCube, 2), DimDelta(WienHyperCube, 2), HCxtractBnd
	//////////////////////////////////
	
	
	///////////////////////////////////
	// Introduce a lifetime broadening to get a spectral function
	Make /O /N = (DimSize(kList,0),kz_number,NumEsteps) Atemp
	SetScale /P x, DimOffset(Spaghetti_3D, 0), DimDelta(Spaghetti_3D, 0), Atemp
	SetScale /P y, (kz_center - kz_range/2), (kz_range/(kz_number-1)), Atemp
	SetScale /P z, (Elow-EF), Estep, Atemp		// The Energy scale of the spectral function has EF == 0
	
	Atemp = 0
	Duplicate /O Atemp Atotal
	
	// Sum over each band
	For(i = 0; i < DimSize(bandlist, 0); i+=1)

		// Extract the DFT bandstructure for band i
		//epsilon[][] = Spaghetti_3D[p][i+1][q]
		HCxtractBnd[][][] = WienHyperCube[p][q][r][(bandlist[i])-1]
		Interp3Dpath HCxtractBnd MasterKlist
		Wave W_interpolated
		For(j = 0; j < kz_number; j+=1)
			epsilon[][j] = W_interpolated[p + j*DimSize(klist, 0)]
		EndFor
		

		// The spectral function is a lorentzian centered at the DFT dispersion, with a width that is linear in binding energy
		//Atemp[][][0,EF_pntnum] = (Abs(EF-epsilon[p][q])/(2* Pi)) / ((Energy[r]-epsilon[p][q])^2 + ((EF-epsilon[p][q])*Gam)^2)
		
		// Each energy interval is the average of the lorentzian over its energy range (more stable when lifetime << Estep)
		Atemp[][][0, EF_pntnum] = (atan((Energy[r] + Estep/2 - epsilon[p][q])/ (ABS(EF-epsilon[p][q])*Gam)) - atan((Energy[r] - Estep/2 - epsilon[p][q])/ (ABS(EF-epsilon[p][q])*Gam)))/ (Pi*Estep)
		
		// Rather than take the value of A(E, k) at each k and E point, the lorentzian is integrated over the energy and k-range spanned by each point (more stable for sharp features)
		// Not yet functional
		//Atemp[][][0, EF_pntnum] = (Energy[r]+Estep/2) * atan((Energy[r] + Estep/2 - (epsilon[p][q]+epsilon[p+1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))
		//Atemp[][][0, EF_pntnum] -= 2*epsilon[p][q]*Gam*(((epsilon[p][q]+epsilon[p+1][q])/2)*atan((Energy[r] + Estep/2 - (epsilon[p][q]+epsilon[p+1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))/(4*ABS(EF-epsilon[p][q])*Gam) + (1/(8*Log(2.718)) * Log(4*(ABS(EF-epsilon[p][q])*Gam)^2 + 4*(Energy[r] + Estep/2 - (epsilon[p][q]+epsilon[p+1][q])/2)^2)))

		//Atemp[][][0, EF_pntnum] -= (Energy[r]-Estep/2) * atan((Energy[r] - Estep/2 - (epsilon[p][q]+epsilon[p+1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))
		//Atemp[][][0, EF_pntnum] += 2*epsilon[p][q]*Gam*(((epsilon[p][q]+epsilon[p+1][q])/2)*atan((Energy[r] - Estep/2 - (epsilon[p][q]+epsilon[p+1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))/(4*ABS(EF-epsilon[p][q])*Gam) + (1/(8*Log(2.718)) * Log(4*(ABS(EF-epsilon[p][q])*Gam)^2 + 4*(Energy[r] - Estep/2 - (epsilon[p][q]+epsilon[p+1][q])/2)^2)))

		//Atemp[][][0, EF_pntnum] -= (Energy[r]+Estep/2) * atan((Energy[r] + Estep/2 - (epsilon[p][q]+epsilon[p-1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))
		//Atemp[][][0, EF_pntnum] += 2*epsilon[p][q]*Gam*(((epsilon[p][q]+epsilon[p-1][q])/2)*atan((Energy[r] + Estep/2 - (epsilon[p][q]+epsilon[p-1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))/(4*ABS(EF-epsilon[p][q])*Gam) + (1/(8*Log(2.718)) * Log(4*(ABS(EF-epsilon[p][q])*Gam)^2 + 4*(Energy[r] + Estep/2 - (epsilon[p][q]+epsilon[p-1][q])/2)^2)))

		//Atemp[][][0, EF_pntnum] = (Energy[r]-Estep/2) * atan((Energy[r] - Estep/2 - (epsilon[p][q]+epsilon[p-1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))
		//Atemp[][][0, EF_pntnum] -= 2*epsilon[p][q]*Gam*(((epsilon[p][q]+epsilon[p-1][q])/2)*atan((Energy[r] - Estep/2 - (epsilon[p][q]+epsilon[p-1][q])/2)/ (ABS(EF-epsilon[p][q])*Gam))/(4*ABS(EF-epsilon[p][q])*Gam) + (1/(8*Log(2.718)) * Log(4*(ABS(EF-epsilon[p][q])*Gam)^2 + 4*(Energy[r] - Estep/2 - (epsilon[p][q]+epsilon[p-1][q])/2)^2)))

		//Atemp[][][0, EF_pntnum] *= -2/Pi
		//Atemp[][][0, EF_pntnum] *= -2/(epsilon[p-1][q] - epsilon[p+1][q])
	
		// A zero-temperature Fermi-function
		If(EF_pntnum+1 < DimSize(Atemp, 2))
			Atemp[][][EF_pntnum+1, ] = 0
		EndIf
	
		Atotal += Atemp
	EndFor
	/////////////////////////////////////


	/////////////////////////////////////
	// Integrate in kz to get the projected 'spectral function', x axis is energy with EF == 0, y axis is k
	Make /O /N = (NumEsteps, DimSize(kList,0)) SimulatedSpectrum
	SetScale /P x, (Elow-EF), Estep, SimulatedSpectrum	
	ux = kList[DimSize(kList, 0)-1][0] - kList[0][0]
	uy = kList[DimSize(kList, 0)-1][1] - kList[0][1] 
	ux /= sqrt(ux^2 + uy^2)
	uy = sqrt(1 - ux^2)
	SetScale /I y, (kList[0][0]*ux + kList[0][1]*uy), (kList[DimSize(kList, 0)-1][0]*ux + kList[DimSize(kList, 0)-1][1]*uy), SimulatedSpectrum
	SimulatedSpectrum = 0
	
	For(i=0; i < kz_number; i+=1)
		SimulatedSpectrum[][] += ((kz_width/kz_number)*Atotal[q][i][p]/(2* Pi)) / ((kz[i] - kz_center)^2 + (kz_width/2)^2)
	EndFor
	//////////////////////////////////////


	Note SimulatedSpectrum, "Simulated ARPES spectrum from HCsimulateARPES"
	Note SimulatedSpectrum, "EF=" + num2str(EF)
	Note SimulatedSpectrum, "kz_center=" + num2str(kz_center)
	Note SimulatedSpectrum, "kz_width=" + num2str(kz_width)
	Note SimulatedSpectrum, "Gam=" + num2str(Gam)
	Note SimulatedSpectrum, "kz_number=" + num2str(kz_number)
	Note SimulatedSpectrum, note($(nameofwave(klist)))

	KillWaves Atemp, epsilon, kz, Energy, Atotal, W_interpolated, HCxtractBnd, MasterKlist

End

///////////////////////////////////
// Rotates the coordinate system of a hypercube. Mainly for systems whose c-axis is not axis #3. 
///////////////////////////////////
Function RotateHC(WienHyperCube, newX, newY, newZ)
	Wave WienHyperCube
	Variable newX, newY, newZ	// Which current axis will become x, y, and z of the rotated cube. Axes numbers are 1, 2, and 3.
								// Ie: For a wave whose c-axis is #2, you might choose newX=1, newY=3, newZ=2.

	// Check input
	If(newX==1)
		If(newY==2)
			If(newZ != 3)
				Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
				Return 0
			EndIf
		ElseIf(newY==3)
			If(newZ != 2)
				Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
				Return 0
			EndIf
		Else
			Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
			Return 0
		EndIf
	Elseif(newX==2)
		If(newY==1)
			If(newZ != 3)
				Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
				Return 0
			EndIf
		ElseIf(newY==3)
			If(newZ != 1)
				Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
				Return 0
			EndIf
		Else
			Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
			Return 0
		EndIf
	ElseIf(newX==3)
		If(newY==1)
			If(newZ != 2)
				Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
				Return 0
			EndIf
		ElseIf(newY==2)
			If(newZ != 1)
				Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
				Return 0
			EndIf
		Else
			Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
			Return 0
		EndIf
	Else
		Print "Error: newX, newY, and newZ must contain all three numbers: 1, 2, 3."
		Return 0
	Endif
	
	Duplicate /O WienHyperCube tempHC
	tempHC = 0
	
	// Since the system is assumed inversion symmetric, the axes can be simply permuted
	Redimension /N = ((DimSize(WienHyperCube,newX-1)),(DimSize(WienHyperCube,newY-1)),(DimSize(WienHyperCube,newZ-1)),(DimSize(WienHyperCube,3))) tempHC
	SetScale /P x, (DimOffset(WienHyperCube, newX-1)), (DimDelta(WienHyperCube, newX-1)), tempHC
	SetScale /P y,  (DimOffset(WienHyperCube, newY-1)), (DimDelta(WienHyperCube, newY-1)), tempHC
	SetScale /P z,  (DimOffset(WienHyperCube, newZ-1)), (DimDelta(WienHyperCube, newZ-1)), tempHC
	
	If(newX==1)
		If(newY==2)
			tempHC = WienHyperCube[p][q][r][s]
		ElseIf(newY==3)
			tempHC = WienHyperCube[p][r][q][s]
		EndIf
	Elseif(newX==2)
		If(newY==1)
			tempHC = WienHyperCube[q][p][r][s]
		ElseIf(newY==3)
			tempHC = WienHyperCube[r][p][q][s]
		EndIf
	ElseIf(newX==3)
		If(newY==1)
			tempHC = WienHyperCube[q][r][p][s]
		ElseIf(newY==2)
			tempHC = WienHyperCube[r][q][p][s]
		EndIf
	Endif
	
	Duplicate /O tempHC WienHyperCube
	//Add a note to the hypercube
	Note WienHyperCube "Rotated by RotateHC() on " + Date()+ " using:"
	Note WienHyperCube "newX = " + num2istr(newX)
	Note WienHyperCube "newY = " + num2istr(newY)
	Note WienHyperCube "newZ = " + num2istr(newZ)
	
	KillWaves tempHC

End


/////////////////////////////////////
// This function expands a hypercube into a BZ that is twice as big in the xy-plane. Use it for calculations with a
// root2-root2 reconstruction. This program assumes that kz is along axis 3, if this isn't true use "RotateHC" first.
// Caution: HC must run from -pi/a to pi/a in kx and ky for this to function correctly.
//
// This function may take a few minutes to run.
/////////////////////////////////////
Function HCc2x2Reconstruct(WienHyperCube)
	Wave WienHypercube
	
	Variable j, k, l, m
	Variable kx, ky, kz
	
	Duplicate /O WienHyperCube tempHC
	
	Redimension /N = (3*DimSize(WienHyperCube,0), 3*DimSize(WienHyperCube,1), DimSize(WienHyperCube,2),DimSize(WienHyperCube,3))  tempHC
	SetScale /I x, -3, 3, tempHC
	SetScale /I y, -3, 3, tempHC

	For(j = 0; j < 3; j+=1)
		For(k = 0; k < 3; k+=1)
			tempHC[(j*DimSize(WienHyperCube,0)), ((j+1)*DimSize(WienHyperCube,0))][(k*DimSize(WienHyperCube,1)), ((k+1)*DimSize(WienHyperCube,1))][][] = WienHyperCube[p-(j*DimSize(WienHyperCube,0))][q-(k*DimSize(WienHyperCube,1))][r][s]
		EndFor
	EndFor

	Make /O /N = (Round(2*DimSize(WienHyperCube, 0)), Round(2*DimSize(WienHyperCube, 1)), DimSize(WienHyperCube,2),DimSize(WienHyperCube,3)) $(NameofWave(WienHyperCube) + "_c2x2")
	
	Wave HC_c2x2 = $(NameofWave(WienHyperCube) + "_c2x2")
	SetScale /I x, -1, 1, HC_c2x2
	SetScale /I y, -1, 1, HC_c2x2
	SetScale /I z, -1, 1, HC_c2x2
	SetScale /P t, DimOffset(WienHyperCube, 3), DimDelta(WienHyperCube,3), HC_c2x2
	
	Make/O /N = (DimSize(tempHC, 0), DimSize(tempHC, 1)) HC_slice
	SetScale /I x, -3, 3, HC_slice
	SetScale /I y, -3, 3, HC_slice
		
	For(l = 0; l < DimSize(HC_c2x2,3); l+=1)
		For(m=0;  m < DimSize(HC_c2x2,2); m+=1)
			For(j = 0; j < DimSize(HC_c2x2,1); j+=1)
				For(k = 0; k < DimSize(HC_c2x2,0); k+=1)
					kx =  (DimOffset(HC_c2x2, 0) + k *DimDelta(HC_c2x2,0) + DimOffset(HC_c2x2, 1) + j *DimDelta(HC_c2x2,1))
					ky =  (DimOffset(HC_c2x2, 1) + j *DimDelta(HC_c2x2,1) - (DimOffset(HC_c2x2, 0) + k *DimDelta(HC_c2x2,0)))
	
					HC_slice = tempHC[p][q][m][l]
	
					HC_c2x2[k][j][m][l] = Interp2D(HC_slice, kx, ky)
				EndFor
			EndFor
		EndFor
	EndFor
		
	KillWaves tempHC, HC_slice

End


/////////////////////////////////////
// This function expands a 2D FS into a BZ that is twice as big. Use it for calculations with a
// root2-root2 reconstruction
/////////////////////////////////////
Function c2x2Reconstruct(FS)
	Wave FS
	Variable j, k

	Variable kx, ky
	
	Make /O /N = (3*DimSize(FS, 0), 3*DimSize(FS, 1)) FS_temp
	SetScale /I x, -3, 3, FS_temp
	SetScale /I y, -3, 3, FS_temp

	For(j = 0; j < 3; j+=1)
		For(k = 0; k < 3; k+=1)
			FS_temp[(j*DimSize(FS,0)), ((j+1)*DimSize(FS,0))][(k*DimSize(FS,1)), ((k+1)*DimSize(FS,1))] = FS[p-(j*DimSize(FS,0))][q-(k*DimSize(FS,1))]
		EndFor
	EndFor

	Make /O /N = (2*DimSize(FS, 0), 2*DimSize(FS, 1)) $(NameofWave(FS) + "_c2x2")
	
	Wave FS_c2x2 = $(NameofWave(FS) + "_c2x2")
	
	SetScale /I x, -1, 1, FS_c2x2
	SetScale /I y, -1, 1, FS_c2x2
	
		For(j = 0; j < DimSize(FS_c2x2,1); j+=1)
			For(k = 0; k < DimSize(FS_c2x2,0); k+=1)
				kx =  (DimOffset(FS_c2x2, 0) + k *DimDelta(FS_c2x2,0) + DimOffset(FS_c2x2, 1) + j *DimDelta(FS_c2x2,1))
				ky =  (DimOffset(FS_c2x2, 1) + j *DimDelta(FS_c2x2,1) - (DimOffset(FS_c2x2, 0) + k *DimDelta(FS_c2x2,0)))
			
				FS_c2x2[k][j] = Interp2D(FS_temp, kx, ky)
			EndFor
		EndFor
		
		KillWaves FS_temp
End

////////////////////////////
// Run c2x2 on a batch of 2D FS's generated from a bandlist (see HyperCube2DFS)
////////////////////////////
Function Batchc2x2(bandlist, kz)
	Wave bandlist
	Variable kz
	Variable i
	
	For(i = 0; i < DimSize(bandlist,0); i+=1)
		Wave FS = $("FS_Band" + num2str(bandlist[i]) + "_kz" + num2str(kz))
	
		c2x2Reconstruct(FS)
	EndFor

End


///////////////////////////////
// This function separates bands from a c2x2 reconstructed calculation into 'real' and 'shadow' bands.
// Use this function after using c2x2Reconstruct.
//
// Input is in the form of a contour-plot pair (plot FS contour on a graph, and use IGOR's extract-contour
// tool to extract the xy pair; Input for this function is the extracted xy pair).
///////////////////////////////
Function ShadowBandsc2x2(FSX, FSY)
	Wave FSX, FSY

	Variable i
	Variable pnt_1 = 0
	Variable pnt_2 = 0
	Variable bnd = 0

	//make waves to store output
	Make /O /N = (1, 2) $(NameofWave(FSX)[0,14] + "SB1"), $(NameofWave(FSX)[0,14] + "SB2")
	
	Wave FS_1 = $(NameofWave(FSX)[0,14] + "SB1")
	Wave FS_2 = $(NameofWave(FSX)[0,14] + "SB2")

	//loop over all points
	For(i = 0; i < DimSize(FSX, 0); i+=1)
	
		//if point is NaN, add NaN to both waves
		If(FSX[i] == NaN)
			Redimension /N = ((pnt_1 + 1), 2) FS_1
			Redimension /N = ((pnt_2 + 1), 2) FS_2
			FS_1[pnt_1][] = NaN
			FS_2[pnt_2][] = NaN
			pnt_1 += 1
			pnt_2 += 1
			
			bnd = 0
			
		//if point is inside c2x2 BZ, append to wave 1, otherwise to wave 2			
		ElseIf(ABS(FSX[i]) + ABS(FSY[i]) < 1)
			// Insert NaNs to break up 'lines' if last point was for other band
			If(bnd ==2)
				Redimension /N = ((pnt_1 + 1), 2) FS_1
				FS_1[pnt_1][] = NaN
				pnt_1 += 1
			EndIf
		
			Redimension /N = ((pnt_1 + 1), 2) FS_1
			FS_1[pnt_1][0] = FSX[i]
			FS_1[pnt_1][1] = FSY[i]
			pnt_1 += 1
			
			bnd = 1
		
		Else
			// Insert NaNs to break up 'lines' if last point was for other band
			If(bnd ==1)
				Redimension /N = ((pnt_2 + 1), 2) FS_2
				FS_2[pnt_2][] = NaN
				pnt_2 += 1
			EndIf
			
			Redimension /N = ((pnt_2 + 1), 2) FS_2
			FS_2[pnt_2][0] = FSX[i]
			FS_2[pnt_2][1] = FSY[i]
			pnt_2 += 1
			
			bnd = 2
			
		EndIf
	EndFor

End


//////////////////////////////////
// Expands a FS contour into multiple BZs.
//
// FS is a 2 dimensional wave, with the first column kx values and the second column ky values.
//////////////////////////////////
Function ExpandFS(FS, NX, NY)
	Wave FS
	Variable NX, NY     // The number of BZs in the x and y directions desired for output (must be integers)

	Variable i, j, cnt
	
	Make /O /N = (((NX*NY)*DimSize(FS,0))+NX*NY-1,2) $(NameOfWave(FS) + "_" + num2str(NX) + "x" + num2str(NY))
	Wave FS_new = $(NameOfWave(FS) + "_" + num2str(NX) + "x" + num2str(NY))
	FS_new = NaN

	For(i = 0; i < NX; i+=1)
		For(j = 0; j < NY; j+=1)
			FS_new[cnt*DimSize(FS,0)+cnt, (cnt+1)* DimSize(FS,0) + cnt -1][0] = FS[(p - cnt*DimSize(FS,0) - cnt)][0] + 2*i
			FS_new[cnt*DimSize(FS,0)+cnt, (cnt+1)* DimSize(FS,0) -1 + cnt][1] = FS[(p - cnt*DimSize(FS,0) - cnt)][1] + 2*j
			cnt +=1
		EndFor
	EndFor

End


///////////////////////////////
// Batch processing for ExpandFS
///////////////////////////////
Function BatchExpandFS(FSList, NX,NY, App)
	Wave /T FSList     // A text file containing the names of all waves to expand
	Variable NX, NY
	Variable App // If App == 1, append FS to top graph
	Variable i
	
	For(i = 0; i<DimSize(FSList,0); i+=1)
		ExpandFS($(FSList[i]), NX,NY)
		
		If(App == 1)
			AppendToGraph $(FSList[i] + "_" + num2str(NX) + "x" + num2str(NY))[][1] vs $(FSList[i] + "_" + num2str(NX) + "x" + num2str(NY))[][0]
		EndIf
	EndFor

End


