      SUBROUTINE get_rwf(rlist, nrlist, infile)

        USE radial
        
        IMPLICIT NONE
        
        ! arguments
        INTEGER, DIMENSION(:), INTENT(in) :: rlist
        INTEGER, INTENT(in) :: nrlist
        CHARACTER(len=*), INTENT(in) :: infile

        ! local variables
        INTEGER :: i, ierr, ios, j, k, npxy
        INTEGER(kind=2) :: npy, naky
        REAL(kind=8) :: ey, dnorm, accy
        CHARACTER(len=6) :: g92rwf
        REAL(kind=8), DIMENSION(:), ALLOCATABLE :: py, qy, ry

        accy = h**6
        
        OPEN (rwfnout, FILE=infile, FORM='unformatted', STATUS='old', iostat=ierr)  

        IF (IERR == 1) THEN 
           WRITE (ISTDE, *) 'Error openning file "', infile, &
                '"' 
           CLOSE(rwfnout) 
           STOP  
        ENDIF

!   Check the file; if not as expected, try again        
        
        READ (rwfnout, IOSTAT=ios) g92rwf 
        IF (IOS/=0 .OR. G92RWF/='G92RWF') THEN 
           WRITE (ISTDE, *) 'This is not a Radial WaveFunction File;' 
           CLOSE(rwfnout) 
           STOP  
        ENDIF

!   Read orbital information from Read Orbitals File
      
        DO
           READ (rwfnout, iostat=ios) npy, naky, ey, npxy
           IF(IOS .NE. 0) EXIT
           ALLOCATE(py(npxy),qy(npxy),ry(npxy))
           READ (rwfnout) (py(k),k=1,npxy), (qy(k),k=1,npxy)
           READ (rwfnout) (ry(k),k=1,npxy)
           DO j = 1, nrlist 
              i = rlist(j)
              IF (NP(i)==NPY .AND. NAK(i)==NAKY) THEN
                 IF ((npxy .NE. npx).AND.(ABS(r(2)-ry(2)).GT.accy)) THEN
                    e(i) = ey
                    CALL INTRPQ (py, qy, npxy, ry, i, DNORM) 
                 ELSE
                    e(i) = ey
                    p(:,i) = py(:)
                    q(:,i) = qy(:)
                 END IF
              END IF
           END DO
           DEALLOCATE(py,qy,ry)
        END DO
        
        CLOSE(rwfnout)

    END SUBROUTINE get_rwf
