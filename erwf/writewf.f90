SUBROUTINE writewf
  USE radial  
!  *
!   Open, write a header and all subshell radial wavefunctions, and    *
!   close the  .rwf  file.                                             *    
  IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: i, j, IERR, npxj 
!-----------------------------------------------
!
!

      OPEN (rwfnout, FILE='rwfn.out', FORM='unformatted', STATUS='replace', iostat=ierr)  
      OPEN (rwfndat, FILE='rwfn.dat', FORM='formatted', STATUS='replace', iostat=ierr)  

      
      IF (IERR == 1) THEN 
         WRITE (ISTDE, *) 'Error when opening ''rwfn.inp'''
         STOP  
      ENDIF 
!
!   Write the file header
!
      WRITE (rwfndat,'(a)') 'OOGRWF'
      WRITE (rwfndat, '(//,a6,2x,i5)') 'npx:', npx
      WRITE (rwfndat, '(a4,6x,d12.6)') 'h:', h

!      
!   Write out the radial wavefunctions 
!
      DO I = 1, NW 
         WRITE (rwfndat, '(//,2x,a7,3X,i2,a2,2x,5D12.4)') 'state:', np(i), nh(i)
         WRITE (rwfndat, '(2x,a7,2x,d12.4)') 'energy:', e(i)
         WRITE (rwfndat, '(//,a7,8x,a4,i2,a,a,2x,a7,i2,a,a)') 'r','p(',np(i),nh(i),')','q(',np(i),nh(i),')'
         DO j = 1, npt(i)
            WRITE (rwfndat, '(d12.4,2x,D12.4,2x,d12.4,2x,i4)') r(j), p(j,i), q(j,i), npt(i)
         END DO
      END DO
      
      !   Binary file output
      WRITE (rwfnout) 'G92RWF'
      DO J = 1, NW
         WRITE (rwfnout) NP(J), NAK(J), E(J), npt(j) 
         WRITE (rwfnout) (P(I,J),I=1,npxj), (Q(I,J),I=1,npxj)
         WRITE (rwfnout) (R(I),I=1,npxj) 
      END DO 
 
      CLOSE(rwfndat)
      CLOSE(rwfnout)

      RETURN  

END SUBROUTINE writewf
