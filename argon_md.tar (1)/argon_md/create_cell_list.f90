subroutine new_nlist(rcell, BOX, TotAtom, ll, hoc)
        rn = BOX/int(BOX/rc)
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
