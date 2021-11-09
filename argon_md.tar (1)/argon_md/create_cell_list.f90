subroutine new_nlist(rcell, Box, TotAtom, ll, hoc)
        use general, only: dp
        use conversions
        
        implicit none
        integer, intent(in):: TotAtom
        real(kind=dp) :: Box

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
