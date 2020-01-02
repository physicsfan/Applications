SUBROUTINE load_isodata
  USE CASE

  IMPLICIT NONE
  
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !----------------------------------------------- 
  CHARACTER(LEN=14), PARAMETER :: signature = 'Atomic number:' 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  INTEGER ::  ios
  CHARACTER(LEN=14) :: str 
  !-----------------------------------------------
  !
  
  !   Check the first record of the file
  !
  READ (isofile, '(A)', iostat=ios) str 
  IF (ios/=0 .OR. str/=signature) THEN 
     WRITE (istdo, *) 'Not an Isotope Data File;' 
     RETURN  
  ENDIF
  !
  !   Atomic number
  !
  READ (isofile, *) z 

!   Nuclear geometry
!
!  READ (isofile, *) 
!  READ (isofile, *) A 
!  READ (isofile, *) 
!  READ (isofile, *) APARM 
!  READ (isofile, *) 
!  READ (isofile, *) CPARM 

!   Nuclear mass  
!  READ (isofile, *) 
!  READ (isofile, *) EMNAMU 

!
!  IF (EMNAMU /= 0.D0) THEN 
!     EMN = EMNAMU/AUMAMU 
!  ELSE 
!     EMN = 0.D0 
!  ENDIF


!
!   Nuclear spin and moments
!
!  READ (isofile, *) 
!  READ (isofile, *) SQN 
!  READ (isofile, *) 
!  READ (isofile, *) DMOMNM 
!  READ (isofile, *) 
!  READ (isofile, *) QMOMB 
    
  RETURN  
  
END SUBROUTINE load_isodata
