subroutine new_nlist(TotAtom, r, Box, Rn, NCell, ll, hoc)
        use general, only: dp
        use conversions
        
        implicit none
        integer, intent(in):: TotAtom, Ncell
        real(kind=dp),intent(inout) :: r(TotAtom,3)
        real(kind=dp) :: Box
        real, intent(inout) :: Rn
        integer :: ll(TotAtom), hoc(0:Ncell-1, 0:Ncell-1, 0:NCell-1)
        integer :: atom, icel(3), i, j, k
        
        
        do i=0, Ncell -1
        	do j=0, NCell-1
        		do k=0, NCell-1
                		hoc(i, j, k) = 0
                	enddo
                enddo
        enddo
        do atom=1, TotAtom
        	icel = r(atom,:)/rn
        	ll(atom) = hoc(icel(1), icel(2), icel(3))
        	hoc(icel(1), icel(2), icel(3)) = atom
        enddo
        return
end subroutine new_nlist
