subroutine atom_energy(atom1,r(atom1,:),Force,PE)
	use general, only: dp,atom1, atom2
	implicit none
	integer :: icel(3)
	real(kind=dp),intent(in) :: r(TotAtom,3)
	real(kind=dp) :: fac2,fac6,r2,df,fc(3),dr(3)
	real(kind=dp),intent(out) :: Force(TotAtom,3),PE
	
	icel=int(r(atom1,:)/RCut)
	do i=icel(1)-1, icel(1)+1
		do j=icel(2)-1, icel(2)+1
			do k=icel(3)-1, icel(3)+1
				atom2=hoc(i,j,k)
				do while (atom2 .ne. 0)
					 if (atom1 .ne. atom2)
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
	return

end subroutine atom_energy
