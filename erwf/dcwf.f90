!=======================================================================
      SUBROUTINE DCWF (n,kappa,Z,E,NPX,R,P,Q)
!=======================================================================
!   This subroutine computes the  Dirac-Coulomb  bound-state orbital
!   radial wavefunction.
!
!   Input:
!
!      n          The (usual) principal quantum number
!      kappa      The relativistic angular quantum number
!      Z          The effective nuclear charge
!      R(1:NTP)   Radial grid
!
!   Output:
!
!      E          The Dirac-Coulomb Eigenenergy (E-mc^2)
!
!      P          r times the large component wavefunction of
!
!      Q          r times the small component wavefunction of
!
!   Call(s) to: CGAMMA
!
!   Written by Farid A Parpia, at Oxford    Last Update: 14 Oct 1992
!
!=======================================================================
!   For a given value of n > 0, there are 2*n-1 eigenfunctions:
!   n  with  kappa = -1,-2,...,-n
!   n-1 with kappa = 1,2,...,n-1
!
!   k = iabs(kappa);  gamma = sqrt[k^2-(alfa*Z)^2]
!
!   N = sqrt[n^2 - 2(n-k)(k-gamma)]
!
!   E(n,k) = c^2 / sqrt[1 + (alfa*Z)^2 / (gamma+n-k)^2 ]
!
!   x = 2*Z/N*r
!
!   P_nk(r) = sqrt[1+E(n,k)/c^2]  N(n,k)  e^-x/2  x^gamma *
!             [(N-kappa) F2  -  (n-k) F1]
!
!   Q_nk(r) = sqrt[1-E(n,k)/c^2]  N(n,k)  e^-x/2  x^gamma *
!             [(N-kappa) F2  +  (n-k) F1]
!
!   N(n,k) = 1/[N*G(2*gamma+1)]  *
!            sqrt [ [Z*G(2*gamma+1+n-k)] / [2(n-k)!(N-kappa)] ]
!
!   G(x) - GAMMA function
!
!   F2 = F(-n+k,2*gamma+1,x);  F1 = F(-n+k+1,2*gamma+1,x)
!
!   F(a,b,x) = 1 + a/b x^1/1!  + [a(a+1)]/[b(b+1)] x^2/2! + ...
!=======================================================================

        IMPLICIT NONE

        !     Arguments
        INTEGER(kind=2) :: n, kappa 
        REAL(kind=8)    :: e, z
        REAL(kind=8)    :: R(*), P(*), Q(*)
        !     Local variables
        INTEGER      :: i, k, npx, nr
        REAL(kind=8) :: a, a1, a2, an1, an2, argi, argr, b, bn, bign
        REAL(kind=8) :: c, dummy, f1, f2, fac, facn, fden, ff, fg
        REAL(kind=8) :: eps, facnr, fkappa, fk, fn, fnr, gamma, gg
        REAL(kind=8) :: nrfac, nt1, nt2, ovlfac, rgamm1, rgamm2, x, y,za
        REAL(kind=8) :: T1(n), T2(n)
      
        C = c_speed
      
! ... Check the input arguments:

      IF(n.le.0)     Stop 'DCWF: Principal quantum number < 0'
      IF(kappa.eq.0) Stop 'DCWF: Kappa quantum number = 0'
      IF(kappa.eq.n) Stop 'DCWF: Kappa quantum number = n'
      IF(kappa.gt.n) Stop 'DCWF: Kappa quantum number > n'
      IF(Z.le.0.d0)  Stop 'DCWF: Nuclear charge is too small, <= 0'
      IF(Z.gt.C)     Stop 'DCWF: Nuclear charge exceeds limit, c'

! ... Now determine all the parameters:

      fn = DBLE (n)
      fkappa = DBLE (kappa)
      k = IABS(INT(kappa));  fk = DBLE (k)
      nr = n-k;  fnr = DBLE (nr)
      Za = Z*alpha
      gamma = SQRT (fk*fk-Za*Za)
      gg = gamma + gamma + 1.d0
      BIGN = SQRT (fn*fn-2.d0*fnr*(fk-gamma))
      EPS = 1.d0 /SQRT(1.d0+(Za/(gamma+fnr))**2)

! ... EPS is the total energy divided by C*C

      E = (1.d0-EPS)*C*C         !  E => E-mc^2

! ... normalization constant N(n,k):

      NRFAC=1; Do I=1,NR; NRFAC = NRFAC*I; End do

      ARGI = 0.d0
      ARGR = gg+FNR;    CALL CGAMMA (ARGR,ARGI,RGAMM1,DUMMY)
      ARGR = gg;        CALL CGAMMA (ARGR,ARGI,RGAMM2,DUMMY)

      FAC = - SQRT (RGAMM1)/(RGAMM2*SQRT (DBLE(NRFAC))) &
            * SQRT (Z/(2.d0*BIGN*BIGN*(BIGN-FKAPPA)))

! ... Ensure that the slope of the large-component function is
! ... positive at the origin:

      IF (KAPPA .GT. 0) FAC = -FAC

      FG = FAC * SQRT(1.d0+EPS)
      FF = FAC * SQRT(1.d0-EPS)

! ...  Now set up the coefficients of the confluent hypergeometric
! ...  functions  F(-NR+1,2*GAMMA+1;RHO)  and  F(-NR,2*GAMMA+1;RHO)
! ...  in the workspace arrays  TA  and  TB, respectively

      nt1=0; if(nr.gt.0) nt1=nr-1;  nt2=nr

      FAC = 1.d0;  FACN = FAC
      A   = -fnr;  AN1  = A + 1.d0;  AN2 = A
      B   =   gg;  BN   = B
      K = 0
      
    2 K = K+1
      FDEN = 1.d0/(FACN*BN)
      if (K .LE. nt1)      T1(K) = AN1*FDEN
      if (K .LE. nt2) then
                           T2(K) = AN2*FDEN
         A = A + 1.d0;     AN1  = AN1*(A+1.d0); AN2 = AN2*A
         B = B + 1.d0;     BN   = BN*B
         FAC = FAC + 1.d0; FACN = FACN*FAC
         go to 2
      end if

! ...  Now tabulate the function over the entire grid

      FAC = (Z+Z)/BIGN
      a1 = FNR; a2 = BIGN-FKAPPA
      Do i = 1,NPX
       x = FAC*R(i); y = x
       F1 = 1.d0; F2 = 1.d0
       K = 0
    3  K = K+1
       if (K .LE. nt1) F1 = F1+T1(K)*y
       if (K .LE. nt2) then;  F2 = F2+T2(K)*y; y=y*x; go to 3;  end if

       F1 = a1*F1;  F2 = a2*F2
       OVLFAC = EXP(-0.5d0*x) * x**gamma

       P(I) = FG*OVLFAC*(F1-F2)
       Q(I) = FF*OVLFAC*(F1+F2)

      End do

      END SUBROUTINE DCWF
