SUBROUTINE writewf
  USE radial  
!  *
!   Open, write a header and all subshell radial wavefunctions, and    *
!   close the  .rwf  file.                                             *    
  IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: i, j, IERR
!-----------------------------------------------
!
!
  
      OPEN (rwfnout, FILE='rwfn.out', FORM='unformatted', STATUS='replace', iostat=ierr)  

      IF (IERR == 1) THEN 
         WRITE (ISTDE, *) 'Error when opening ''rwfn.inp'''
         STOP  
      ENDIF 
      
      !   Binary file output
      WRITE (rwfnout) 'G92RWF'
      DO J = 1, NW
         WRITE (rwfnout) NP(J), NAK(J), E(J), npt(j) 
         WRITE (rwfnout) (P(I,J),I=1,npt(j)), (Q(I,J),I=1,npt(j))
         WRITE (rwfnout) (R(I),I=1,npt(j)) 
      END DO 
 
      CLOSE(rwfnout)

      RETURN  

END SUBROUTINE writewf
