!> generate radial wave functions
SUBROUTINE genrwf
  USE radial

  IMPLICIT NONE
  
  INTEGER  :: i, j, loc, nradial, nelectron, nlist, nrlist
  LOGICAL :: all
  CHARACTER(len=12) :: convert, infile
  CHARACTER(len=4*nnnw) :: str
  LOGICAL, DIMENSION(:), ALLOCATABLE  :: SET
  INTEGER, DIMENSION(:), ALLOCATABLE  :: list, rlist
  
  !> Get info for run (may be put in separate routine)
  WRITE (istdo, '(a,/)') 'Loading CSF file ... Header only'
  CALL load_orbs(nfile)
  WRITE (istdo, *) 'There are/is ', NW, ' relativistic subshells;'
  WRITE (istdo, *) ''
  CALL load_isodata
  CALL radial_grid
  
  !> Allocate all needed arrays
  CALL allocate_radials
  ALLOCATE(set(nw),list(nw),rlist(nw))
  ALLOCATE(sigma(nw),zz(npx))
  
  !> Initializations
  set = .FALSE.
  P = 0.d0; q = 0.d0
  nelectron = 0; nlist = 0
  
  !> calculate screening parameters
  !...Core orbitals
  DO i = 1, ncore 
     sigma(I) = nelectron + (nkj(i)+1)/2 
     nelectron = nelectron + nkj(i) + 1 
  END DO
  
  !...Peel orbitals
  sigma(ncore+1:nw) = nelectron 
  
  
  !> calcilate point nucleus potential
  zz(:npx) = z
  
  
  ! write out complete list of subshell radial wave functions
  ! repeat until all subshells are written out
  DO
     ! print the remaining orbitals that need to be determined
     nlist=0
     DO i= 1, nw
        IF (set(i) ) CYCLE
        nlist = nlist +1
        list(nlist) = i
     END DO
     IF (nlist .GT. 0) THEN          
        PRINT '(/,a,100l2)', 'set ', (set(i), i=1,nw)
        WRITE(istdo, '(/,a)') "The remaining orbitals are:"
        WRITE(istdo, '(127(A,1x))') (el(list(j)), j=1,nlist)
        DO
           WRITE (istdo, '(/,a)') 'Choose one below' 
           WRITE (istdo, '(a)') '    1 -- GRASP92 File' 
           WRITE (istdo, '(a)') '    2 -- Screened Hydrogenic' 
           READ (istdi, *) nradial
           IF (nradial==1) THEN
              WRITE (ISTDE, *) 'Enter the file name (Null then "rwfn.out")' 
              READ (ISTDI, '(A)') INFILE
              infile = trim(adjustl(infile))
              IF (LEN_TRIM(INFILE) == 0) INFILE = 'rwfn.out'
              EXIT
           ELSE IF (nradial .NE. 2) THEN
              WRITE (istdo, *) nradial, 'is not a valid choice'
              STOP
           ELSE
              EXIT
           ENDIF
        END DO
        
        ! determin the set of orbitals selected
        WRITE(istdo,'(a)') 'Enter the list of relativistic subshells:'
        READ(istdi, '(a)') str
        CALL getrlist(nlist, list, str, nrlist, rlist)
        PRINT *, 'nrlist:', nrlist
        PRINT *, 'rlist:', rlist
        IF (nradial == 1) THEN 
           CALL get_rwf(rlist, nrlist, infile)  
        ELSE 
           CALL get_hydrogenic(rlist, nrlist)
        ENDIF
        DO i = 1, nrlist
           set(rlist(i)) = .TRUE.
        END DO
     ELSE
        EXIT
     END IF
  END DO
  
  
  WRITE (istdo,'(/,a,a)') 'All required subshell radial wavefunctions ', &
       'have been estimated:' 
  
  ! Summary
  WRITE (istdo, '(3x,a5,6x,a,11x,a4,7x,a4,7x,a4,7x,a4,/)')'state','e',&
       'r(2)','P(2)','Q(2)' 
  !
  DO I = 1, NW 
     WRITE (istdo, '(2X,i2,a2,2x,5D12.4)') np(i), nh(i), e(i), r(2), &
          p(2,i), q(2,i) 
  END DO
  
  
END SUBROUTINE genrwf
