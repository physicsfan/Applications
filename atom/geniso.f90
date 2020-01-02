!=========================================================================
!                                                                      
!>   Generates the isotope data file for the GRASP92 suite of codes.   
!!                                                                      
!                                                                       
      PROGRAM RNUCLEUS
!==========================================================================
!   M o d u l e s 
!-----------------------------------------------
      USE case
      USE nucleus_m 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR, NENEU 
      REAL(kind=8) ::  TPARM, AMAMU, EBIND, EMNAMU, SQN, DMOMNM, QMOMB 
      LOGICAL :: YES 
!-----------------------------------------------
!
!   Startup message

       WRITE (istdout, *) 'RNUCLEUS'
       WRITE (istdout, *) 'This program defines the nucleus'
       WRITE (istdout, *) 'Output file: isodata' 
!
!   File  isodata   is FORMATTED
!
      CALL OPEN (isodatafile, FILE='isodata', FORM='fomatted', STATUS='new', IOSTATUS=ierr) 
!
      IF (IERR = 0) THEN 
         rewind(isodatafile)
      ELSE
         WRITE (ISTDE, *) 'Error when opening isodata' 
         STOP  
      ENDIF 
!
      WRITE (ISTDIN, *) 'Enter the atomic number:' 
      READ (istdin, *) Z 
      WRITE (isodatafile, 300) 'Atomic number:' 
      WRITE (isodatafile, *) Z 
!
      WRITE (ISTDIN, *) 'Enter the mass number (0 if the', &
         ' nucleus is to be modelled as a point source:' 
      READ (istdin, *) A 
      WRITE (isodatafile, 300) 'Mass number (integer) :' 
      WRITE (isodatafile, *) A 
 
      IF (A == 0.0D00) THEN 
         CPARM = 0.0D00 
         APARM = 0.0D00 
      ELSE 
         RRMS = rrms_value(int(z),int(a))
         WRITE (ISTDIN, *) 'The default root mean squared', ' radius is ', &
                 RRMS,'fm;  (', trim(rrms_source(int(z),int(a))), ')'
         TPARM = 2.30D00 
         WRITE (ISTDIN, *) '  the default nuclear skin thickness', ' is ', TPARM,'fm;'
         WRITE (ISTDIN, *) 'Revise these values?' 
         YES = GETYN() 
         IF (YES) THEN 
            WRITE (ISTDIN, *) 'Enter the root mean squared', &
               ' radius of the nucleus (in fm):' 
            READ (istdin, *) RRMS 
            WRITE (ISTDIN, *) 'Enter the skin thickness of', &
               ' the nucleus (in fm):' 
            READ (istdin, *) TPARM 
         ENDIF 
         APARM = TPARM/(4.0D00*LOG(3.0D00)) 
         CALL GETCPR (RRMS, APARM, CPARM) 
      ENDIF 
      WRITE (isodatafile, 300) 'Fermi distribution parameter a:' 
      WRITE (isodatafile, *) APARM 
      WRITE (isodatafile, 300) 'Fermi distribution parameter c:' 
      WRITE (isodatafile, *) CPARM 
!
      WRITE (ISTDIN, *) 'Enter the mass of the neutral', &
         ' atom (in amu) (0 if the nucleus is to be static):' 
      READ (istdin, *) AMAMU 
      IF (AMAMU /= 0.0D00) THEN 
        NENEU = NINT(Z) 
        EBIND = 0.D0
        EMNAMU = AMAMU - EMEAMU*DBLE(NENEU) - EMEAMU*EBIND/ALFAI**2 
      ELSE 
         EMNAMU = 0.0D00 
      ENDIF 
!
      WRITE (isodatafile, 300) 'Mass of nucleus (in amu):' 
      WRITE (isodatafile, *) EMNAMU 
!
      WRITE (ISTDIN, *) 'Enter the nuclear spin quantum', &
         ' number (I) (in units of h / 2 pi):' 
      READ (istdin, *) SQN 
      WRITE (isodatafile, 300) 'Nuclear spin (I) (in units of h / 2 pi):' 
      WRITE (isodatafile, *) SQN 
!
      WRITE (ISTDIN, *) 'Enter the nuclear dipole moment', &
         ' (in nuclear magnetons):' 
      READ (istdin, *) DMOMNM 
      WRITE (isodatafile, '(A)') 'Nuclear dipole moment (in nuclear magnetons):' 
      WRITE (isodatafile, *) DMOMNM 
!
      WRITE (ISTDIN, *) 'Enter the nuclear quadrupole', ' moment (in barns):' 
      READ (istdin, *) QMOMB 
      WRITE (isodatafile, 300) 'Nuclear quadrupole moment (in barns):' 
      WRITE (isodatafile, *) QMOMB 
!
      CLOSE(isodatafile) 
!
      END PROGRAM RNUCLEUS
!======================================================================
!                                                                      *
      REAL(dp) FUNCTION ESTRMS (APARM, CPARM) 
!                                                                      *
!>   Determines the root mean square radius for a Fermi nucleus given   *
!>  the parameters `c' (CPARM) and `a' (APARM). We use the formalism   *
!>  developed in F. A. Parpia and A. K. Mohanty ``Relativistic basis   *
!>  set calculations for atoms with Fermi nuclei'' Phys Rev A (1992)   *
!>  in press.                                                          *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
!                                                                      *
!----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp) , INTENT(IN) :: APARM, cparm
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: PI, SQTBF, ABC, PABC, CBAM, DNUMER, DDENOM 
!-----------------------------------------------
!
      PI = 4.0D00*ATAN(1.0D00) 
      SQTBF = SQRT(3.0D00/5.0D00) 
!
      ABC = APARM/CPARM 
      PABC = PI*ABC 
      CBAM = -CPARM/APARM 
      DNUMER = 1.0D00 + (10.0D00/3.0D00)*PABC**2 + (7.0D00/3.0D00)*PABC**4 - &
         120.0D00*ABC**5*SKFUN(5,CBAM) 
      DDENOM = 1.0D00 + PABC**2 - 6.0D00*ABC**3*SKFUN(3,CBAM) 
      ESTRMS = CPARM*SQTBF*SQRT(DNUMER/DDENOM) 
!
      RETURN  
      END FUNCTION ESTRMS 
!======================================================================
!                                                                      *
      SUBROUTINE GETCPR(RRMS, APARM, CPARM) 
!                                                                      *
!   Determines the parameter `c' (CPARM) for a Fermi nucleus,  given   *
!   the root mean square radius (RRMS) and the parameter `a' (APARM).  *
!   We use the formalism developed in F. A. Parpia and A. K. Mohanty   *
!   ``Relativistic basis set  calculations for atoms with  Fermi nu-   *
!   clei'' Phys Rev A (1992) in press.                                 *
!                                                                      *
!   Call(s) to: ESTRMS.                                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
!----------------------------------------------------------------------                                                                      *
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: RRMS, APARM
      REAL(DOUBLE) , INTENT(OUT) :: CPARM 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: ACCY, CPMIN, CPMAX, CPTRY, RMSTRY 
!-----------------------------------------------
!
      ACCY = 1.0D-12 
!
!   Bracket CPARM with a lower and upper limit
!
!   Lower limit
!
      CPMIN = 0.5D00*RRMS
      DO WHILE(ESTRMS(APARM,CPMIN) > RRMS) 
         CPMIN = 0.5D00*CPMIN 
      END DO 
!
!   Upper limit
!
      CPMAX = 2.0D00*RRMS
      DO WHILE(ESTRMS(APARM,CPMAX) < RRMS) 
         CPMAX = 2.0D00*CPMAX 
      END DO 
!
!   Find CPARM by the method of bisection
!
      CPTRY = 0.5D00*(CPMAX + CPMIN) 
!
      RMSTRY = ESTRMS(APARM,CPTRY) 
!
      IF (RMSTRY > RRMS) THEN 
         CPMAX = CPTRY 
      ELSE 
         CPMIN = CPTRY 
      ENDIF 
      DO WHILE((CPMAX - CPMIN)/(CPMAX + CPMIN)>ACCY .AND. ABS(RMSTRY-RRMS)/RRMS&
         >ACCY) 
         CPTRY = 0.5D00*(CPMAX + CPMIN) 
!
         RMSTRY = ESTRMS(APARM,CPTRY) 
!
         IF (RMSTRY > RRMS) THEN 
            CPMAX = CPTRY 
         ELSE 
            CPMIN = CPTRY 
         ENDIF 
!
      END DO 
!
      CPARM = CPTRY 
!
      RETURN  
      END SUBROUTINE GETCPR 
!======================================================================
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION SKFUN (K, X) 
!                                                                      *
!   Computes the function                                              *
!                                            n  nx                     *
!                             infinity   (-1)  e                       *
!                     S (x) =   Sum      -------                       *
!                               n=1          k                         *
!                                           n                          *
!                                                                      *
!   See, for instance, F. A. Parpia and A. K. Mohanty ``Relativistic   *
!   basis set calculations for atoms with Fermi nuclei'', Phys Rev A   *
!   (1992) in press.                                                   *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
!                                                                      *
!----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K 
      REAL(dp) , INTENT(IN) :: X 
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: DNUMER, EN, BASE, DELTA 
      REAL(dp), PARAMETER :: QUASIZERO = 1.D-15 
!-----------------------------------------------
 
      BASE = -EXP(X) 
      DNUMER = BASE 
      EN = 1.0D00 
      DELTA = DNUMER
      SKFUN = DELTA 
      DO WHILE(ABS(DELTA/SKFUN) > QUASIZERO) 
         DNUMER = DNUMER*BASE 
         EN = EN + 1.0D00 
         DELTA = DNUMER/EN**K 
         SKFUN = SKFUN + DELTA 
      END DO 
!
      RETURN  
      END FUNCTION SKFUN 
