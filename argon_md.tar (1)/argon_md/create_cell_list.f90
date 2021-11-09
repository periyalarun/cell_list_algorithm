subroutine new_nlist(TotAtom, Box, NCell, ll, hoc)
        use general, only: dp
        use conversions
        
        implicit none
        integer, intent(in):: TotAtom, Ncell
        real(kind=dp) :: Box
        integer :: ll(TotAtom), hoc(Ncell, Ncell)
        rn = Box/int(Box/rc)
        do icel=0, ncel -1
                hoc(icel) = 0
        enddo
        do i=1, TotAtom
                icel = int(r/rn)
                ll(i) = hoc(icel)
                hoc(icel) = i
        enddo
        return
end subroutine new_nlist
