#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//=============================================================
function wannier2ek_hs(hisymk, Hmatr)
wave hisymk //hi symmetry path
wave Hmatr //hmatrix from wannier; cloumns: 0,1,2-hopping range; 3,4-position in matrix; 5-hopping value

variable Hsize=3
make/o/n=(dimsize(hisymk,0),Hsize) hisymE

variable kind,m, kx, ky, kz
make/c/o/n=(Hsize,Hsize) tmphamE
make/o/n=(Hsize) mxevreal

for (kind=0;kind<dimsize(hisymk,0);kind+=1)
			tmphamE=0
			kx =hisymk[kind][0]
			ky =hisymk[kind][1]
			kz =hisymk[kind][2]

			for (m=0;m<dimsize(Hmatr,0);m+=1)
				if ((Hmatr[m][3]<=Hsize)*(Hmatr[m][4]<=Hsize))
					tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])*Exp(cmplx(0,pi)*kz*Hmatr[m][2])
				endif
			endfor
			if (1)
				MatrixEigenV tmphamE
				wave/C W_eigenvalues
				mxevreal = real(W_eigenvalues[p])
				sort mxevreal,mxevreal
				hisymE[kind][] = mxevreal[q]
			else
				hisymE[kind][] = tmphamE[q][q]
			endif

//			for (m=0;m<dimsize(Hmatr,0);m+=1)
//				tmpham[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*Cos(pi*kx*Hmatr[m][0])*Cos(pi*ky*Hmatr[m][1])*Cos(pi*kz*Hmatr[m][2])*2^((Hmatr[m][0]!=0)+(Hmatr[m][1]!=0)+(Hmatr[m][2]!=0))
//			endfor
endfor

end
//============================================================
function wannier2generateQW_hs(hisymk, Hmatr, N)
//only diagonal elements are supported (no hybridization)
wave hisymk //hi symmetry path
wave Hmatr //hmatrix from wannier; cloumns: 0,1,2-hopping range; 3,4-position in matrix; 5-hopping value
variable N // # of sites in z-direction

make/o/n=(dimsize(hisymk,0),3*N) hisymEqw

variable kind,m, kx, ky, kz, nn
make/o/n=(3,3) tmphamE, tmphamT
make/o/n=3 mxevreal

for (kind=0;kind<dimsize(hisymk,0);kind+=1)
			tmphamE=0
			tmphamT=0
			kx =hisymk[kind][0]
			ky =hisymk[kind][1]
			kz =hisymk[kind][2]
			for (m=0;m<dimsize(Hmatr,0);m+=1)
				if (Hmatr[m][2]==0)
					tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*Cos(pi*kx*Hmatr[m][0])*Cos(pi*ky*Hmatr[m][1])*2^((Hmatr[m][0]!=0)+(Hmatr[m][1]!=0))
				elseif (Hmatr[m][2]==1)
					tmphamT[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*Cos(pi*kx*Hmatr[m][0])*Cos(pi*ky*Hmatr[m][1])*2^((Hmatr[m][0]!=0)+(Hmatr[m][1]!=0))
				endif	
				//Cos(pi*(nn+1*0)/(N+1))
			endfor
			for(nn=0;nn<N; nn+=1)
				hisymEqw[kind][3*nn] = tmphamE[0][0] + 2*tmphamT[0][0]*cos(pi*(nn+1)/(N+1))
				hisymEqw[kind][3*nn+1] = tmphamE[1][1] + 2*tmphamT[1][1]*cos(pi*(nn+1)/(N+1))
				hisymEqw[kind][3*nn+2] = tmphamE[2][2] + 2*tmphamT[2][2]*cos(pi*(nn+1)/(N+1))
			endfor
endfor

end
//============================================================
function wannier2generateQW(nx, ny, nz, Hmatr, N)
//only diagonal elements are supported (no hybridization)
wave Hmatr //hmatrix from wannier; cloumns: 0,1,2-hopping range; 3,4-position in matrix; 5-hopping value
variable N // # of sites in z-direction
variable nx, ny, nz //number of k points in 3 directions

make/o/n=(nx, ny, nz, N*3) TBQWcube
setscale/I x, -1, 1, TBQWcube
setscale/I y, -1, 1, TBQWcube
setscale/I z, -1, 1, TBQWcube

variable kind,m, kx, ky, kz, nn, k0,k1,k2
make/o/n=(3,3) tmphamE, tmphamT
make/o/n=3 mxevreal

for (k0=0;k0<dimsize(TBQWcube, 0);k0+=1)
	for (k1=0;k1<dimsize(TBQWcube, 1);k1+=1)
		for (k2=0;k2<dimsize(TBQWcube, 2);k2+=1)
			tmphamE=0
			tmphamT=0
			kx = Dimoffset(TBQWcube,0) + dimdelta(TBQWcube,0)*k0
			ky = Dimoffset(TBQWcube,1) + dimdelta(TBQWcube,1)*k1
			kz = Dimoffset(TBQWcube,2) + dimdelta(TBQWcube,2)*k2

			for (m=0;m<dimsize(Hmatr,0);m+=1)
				if (Hmatr[m][2]==0)
					tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*Cos(pi*kx*Hmatr[m][0])*Cos(pi*ky*Hmatr[m][1])*2^((Hmatr[m][0]!=0)+(Hmatr[m][1]!=0))
				elseif (Hmatr[m][2]==1)
					tmphamT[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*Cos(pi*kx*Hmatr[m][0])*Cos(pi*ky*Hmatr[m][1])*2^((Hmatr[m][0]!=0)+(Hmatr[m][1]!=0))
				endif	
				//Cos(pi*(nn+1*0)/(N+1))
			endfor
			for(nn=0;nn<N; nn+=1)
				TBQWcube[k0][k1][k2][3*nn] = tmphamE[0][0] + 2*tmphamT[0][0]*cos(pi*(nn+1)/(N+1))
				TBQWcube[k0][k1][k2][3*nn+1] = tmphamE[1][1] + 2*tmphamT[1][1]*cos(pi*(nn+1)/(N+1))
				TBQWcube[k0][k1][k2][3*nn+2] = tmphamE[2][2] + 2*tmphamT[2][2]*cos(pi*(nn+1)/(N+1))
			endfor
		endfor
	endfor
endfor

end
//=========================================================
//=========================================================
function wannier2cube(nx, ny, nz, Hmatr)
wave Hmatr //hmatrix from wannier; cloumns: 0,1,2-hopping range; 3,4-position in matrix; 5-hopping value
variable nx, ny, nz //number of k points in 3 directions

variable Hsize=5
variable NoHyb = 0
make/o/n=(nx, ny, nz, Hsize) TBcube
setscale/I x, -1, 1, TBcube
setscale/I y, -1, 1, TBcube
setscale/I z, -1, 1, TBcube

variable kind,m, kx, ky, kz, nn, k0,k1,k2
make/o/n=(Hsize,Hsize)/C tmphamE
make/o/n=(Hsize) mxevreal

for (k0=0;k0<dimsize(TBcube, 0);k0+=1)
	for (k1=0;k1<dimsize(TBcube, 1);k1+=1)
		for (k2=0;k2<dimsize(TBcube, 2);k2+=1)
			tmphamE=0
			kx = Dimoffset(TBcube,0) + dimdelta(TBcube,0)*k0
			ky = Dimoffset(TBcube,1) + dimdelta(TBcube,1)*k1
			kz = Dimoffset(TBcube,2) + dimdelta(TBcube,2)*k2

			for (m=0;m<dimsize(Hmatr,0);m+=1)
				if ((Hmatr[m][3]<=Hsize)*(Hmatr[m][4]<=Hsize))
					if (NoHyb==1)
						if (Hmatr[m][3]==Hmatr[m][4])
							tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])*Exp(cmplx(0,pi)*kz*Hmatr[m][2])
						endif
					else
						tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])*Exp(cmplx(0,pi)*kz*Hmatr[m][2])
					endif
				endif
			endfor
			MatrixEigenV tmphamE
			wave/C W_eigenvalues
			mxevreal = real(W_eigenvalues[p])
			sort mxevreal, mxevreal			
			TBcube[k0][k1][k2][] = mxevreal[s]

		endfor
	endfor
endfor

end
//=========================================================
function wannier2QWcube(nx, ny, nz, Hmatr, N)
// 
wave Hmatr //hmatrix from wannier; cloumns: 0,1,2-hopping range; 3,4-position in matrix; 5-hopping value
variable nx, ny, nz //number of k points in 3 directions
variable N // # qw states

variable Hsize=5
make/o/n=(nx, ny, nz, N*Hsize) TBQWcube
setscale/I x, -1, 1, TBQWcube
setscale/I y, -1, 1, TBQWcube
setscale/I z, -1, 1, TBQWcube

variable kind,m, kx, ky, kz, nn, k0,k1,k2
make/o/n=(Hsize*N,Hsize*N)/C tmphamE
make/o/n=(Hsize*N) mxevreal

for (k0=0;k0<nx;k0+=1)
	for (k1=0;k1<ny;k1+=1)
			tmphamE=0
			kx = Dimoffset(TBQWcube,0) + dimdelta(TBQWcube,0)*k0
			ky = Dimoffset(TBQWcube,1) + dimdelta(TBQWcube,1)*k1
//			kz = Dimoffset(TBQWcube,2) + dimdelta(TBQWcube,2)*k2

			for (m=0;m<dimsize(Hmatr,0);m+=1)
				if ((Hmatr[m][3]<=Hsize)*(Hmatr[m][4]<=Hsize))
					if (Hmatr[m][2]==0)
						for (nn=0; nn<N; nn+=1)
							tmphamE[nn*Hsize+Hmatr[m][3]-1][nn*Hsize+Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])
						endfor						
					elseif (Hmatr[m][2]==1)
						for (nn=0; nn<N-1; nn+=1)
							tmphamE[nn*Hsize+Hmatr[m][3]-1][(nn+1)*Hsize+Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])
						endfor
					elseif (Hmatr[m][2]==-1)
						for (nn=1; nn<N; nn+=1)
							tmphamE[nn*Hsize+Hmatr[m][3]-1][(nn-1)*Hsize+Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])
						endfor
					endif
				endif
			endfor
			MatrixEigenV tmphamE
			wave/C W_eigenvalues
			mxevreal = real(W_eigenvalues[p])
			sort mxevreal, mxevreal			
			TBQWcube[k0][k1][][] = mxevreal[s]

	endfor
endfor

end
//=========================================================
//=========================================================
function wannier2slabcube(nx, ny, nz, Hmatr)
wave Hmatr //hmatrix from wannier; cloumns: 0,1,2-hopping range; 3,4-position in matrix; 5-hopping value
variable nx, ny, nz //number of k points in 3 directions

variable Hsize=10
variable NoHyb = 0
make/o/n=(nx, ny, nz, Hsize) TBcube
setscale/I x, -1, 1, TBcube
setscale/I y, -1, 1, TBcube
setscale/I z, 0, 1, TBcube

variable kind,m, kx, ky, kz, nn, k0,k1,k2
make/o/n=(Hsize,Hsize)/C tmphamE
make/o/n=(Hsize) mxevreal

for (k0=0;k0<dimsize(TBcube, 0);k0+=1)
	for (k1=0;k1<dimsize(TBcube, 1);k1+=1)
		for (k2=0;k2<dimsize(TBcube, 2);k2+=1)
			tmphamE=0
			kx = Dimoffset(TBcube,0) + dimdelta(TBcube,0)*k0
			ky = Dimoffset(TBcube,1) + dimdelta(TBcube,1)*k1
			kz = Dimoffset(TBcube,2) + dimdelta(TBcube,2)*k2

			for (m=0;m<dimsize(Hmatr,0);m+=1)
				if ((Hmatr[m][3]<=Hsize)*(Hmatr[m][4]<=Hsize))
					if (NoHyb==1)
						if (Hmatr[m][3]==Hmatr[m][4])
//							tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])*Exp(cmplx(0,pi)*kz*Hmatr[m][2])
							tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])
						endif
					else
//						tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])*Exp(cmplx(0,pi)*kz*Hmatr[m][2])
						tmphamE[Hmatr[m][3]-1][Hmatr[m][4]-1]+=Hmatr[m][5]*exp(cmplx(0,pi)*kx*Hmatr[m][0])*Exp(cmplx(0,pi)*ky*Hmatr[m][1])
					endif
				endif
			endfor
			MatrixEigenV tmphamE
			wave/C W_eigenvalues
			mxevreal = real(W_eigenvalues[p])
			sort mxevreal, mxevreal			
			TBcube[k0][k1][k2][] = mxevreal[s]

		endfor
	endfor
endfor

end
//=========================================================




