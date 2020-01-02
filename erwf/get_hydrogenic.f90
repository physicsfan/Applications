      SUBROUTINE get_hydrogenic(rlist,nrlist)

        IMPLICIT NONE
        
        ! arguments
        INTEGER, DIMENSION(:), INTENT(in) :: rlist
        INTEGER, INTENT(in) :: nrlist
        
        ! local variables
        INTEGER :: i, j, k
        REAL(kind=8) :: zeff
        
        WRITE (istdo, '(/,a)') '***** Screening parameters ******' 

        DO j = 1, nrlist
           i = rlist(j)
           WRITE (istdo, '(I2,a2,F10.2)') np(i), nh(i), sigma(i)
           zeff = z - sigma(i)                
           CALL DCWF(np(i), nak(i), zeff, e(i), npx, r, p(:,i), q(:,i))

           !   Determine the effective maximum tabulation point
           k=npx
           DO
              IF (ABS(p(k,i)) .GT. 1.d-16) EXIT
              k = k - 1
           END DO
           npt(i) = k
           
        END DO
        
      END SUBROUTINE get_hydrogenic
