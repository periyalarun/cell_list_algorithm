subroutine force_calc(TotAtom,Box,Rcut,r,Sig,Eps,Force,PE,NCell,ll,hoc)	! Arun 1) New parameters added and to be checked in other callers and declarations
 use general, only: dp,atom1,atom2
 implicit none 
 integer,intent(in) :: TotAtom, NCell	!Arun 2) Ncell variable declared
 real(kind=dp),intent(in) :: r(TotAtom,3)
 real(kind=dp) :: fac2,fac6,r2,df,fc(3),dr(3),Ecut,R2cut
 real(kind=dp),intent(out) :: Force(TotAtom,3),PE

 real(kind=dp) :: Box,Rcut,Eps,Sig
 integer :: ll(TotAtom), hoc(0:Ncell-1, 0:Ncell-1, 0:NCell-1) ! Arun 3) ll and hoc declared
 integer :: i, j, k, icel(3)
 R2cut=Rcut*Rcut 

 ! calculating Ecut at Rcut 
 fac2=Sig*Sig/R2cut
 fac6=fac2*fac2*fac2
 Ecut=4.d0*Eps*fac6*(fac6-1)

 ! calculting force and energy 
 do atom1=1,TotAtom-1	! Arun 4) New algorithm for cell list
	!call atom_energy(atom1,r,Force,PE)
	icel=int(r(atom1,:)/RCut)
	do i=(icel(1)-1), (icel(1)+1)
		do j=(icel(2)-1), (icel(2)+1)
			do k=(icel(3)-1), (icel(3)+1)
				atom2=hoc(i,j,k)
				do while (atom2 .ne. 0)
					write(*,*) atom2
					 if ( atom1 .ne. atom2 ) then
					 	dr=r(atom1,:)-r(atom2,:)
						dr=dr-Box*anint(dr/Box)
						r2=dot_product(dr,dr)
						r2=1/r2
						fac2=r2*Sig*Sig
						fac6=fac2*fac2*fac2
						df=48.d0*Eps*r2*fac6*(fac6-0.5d0)
						fc(2)=df*dr(2)
						fc(3)=df*dr(3)
						Force(atom1,:)=Force(atom1,:)+fc(:) 
						Force(atom2,:)=Force(atom2,:)-fc(:) 
						PE=PE+4.d0*Eps*fac6*(fac6-1)-Ecut        !shifted to zero at cutoff 
					else
				 		return
				 	endif
				 	atom2 = ll(atom2)
				enddo
			enddo
		enddo
	enddo
 enddo
 
 PE=0.d0
 Force=0.d0
 do atom1=1,TotAtom-1
   do atom2=atom1+1,TotAtom
     dr=r(atom1,:)-r(atom2,:)
     dr=dr-Box*anint(dr/Box)
     r2=dot_product(dr,dr)
     if(r2<=R2cut) then          ! r2cut  is square of rcut  
       r2=1/r2 
       fac2=r2*Sig*Sig 
       fac6=fac2*fac2*fac2 
       df=48.d0*Eps*r2*fac6*(fac6-0.5d0)
       fc(1)=df*dr(1)
       fc(2)=df*dr(2)
       fc(3)=df*dr(3)
       Force(atom1,:)=Force(atom1,:)+fc(:) 
       Force(atom2,:)=Force(atom2,:)-fc(:) 
       PE=PE+4.d0*Eps*fac6*(fac6-1)-Ecut        !shifted to zero at cutoff 
     endif 
   enddo
 enddo 

!print

! writing forces 
! open(unit=200,file='md.force',action='write') 
!   write(200,*)  
! do atom1=1,TotAtom
!   write(200,"(3F15.8)") Force(atom1,:)
! enddo

end subroutine force_calc 
