! File containing the simulator subroutines (Schrodinger-Poisson solvers)
MODULE fonctions

  IMPLICIT NONE

  INTEGER, PARAMETER :: wp = KIND(1.D0)
  INTEGER, PARAMETER :: nMax = 10005, nd = 301, pRes = 301, rap = 13, pMax = (pRes-1)*rap+1
  INTEGER, PARAMETER  :: n0 = 500
  REAL(wp), PARAMETER :: L = 1.35E-7_wp
  REAL(wp), PARAMETER :: e = 1.602189E-19_wp
  REAL(wp), PARAMETER :: me = 9.10953E-31_wp
  REAL(wp), PARAMETER :: m = 0.067_wp*me
  REAL(wp), PARAMETER :: hBar = 1.054589E-34_wp
  REAL(wp), PARAMETER :: kb = 1.380662E-23_wp    
  REAL(wp), PARAMETER :: Ef = 6.709741104586385E-21_wp
  REAL(wp), PARAMETER :: T = 300._wp
  REAL(wp), PARAMETER :: epsilon = 8.854187817E-12_wp*11.44_wp
  REAL(wp), PARAMETER :: Vref = T*kb/e
  REAL(wp), PARAMETER :: V1 = 0.3_wp
  REAL(wp), PARAMETER :: densD2 = 5.E21_wp, densD1 = 1.E24_wp
  REAL(wp), PARAMETER :: pi = 3.141592653589793238462643_wp
  REAL(wp), DIMENSION(6), PARAMETER :: bar = (/ 5.E-8_wp, 6.E-8_wp, 6.5E-8_wp, 7.E-8_wp, 7.5E-8_wp, 8.5E-8_wp /)
  COMPLEX(wp), PARAMETER :: i=CMPLX(0._wp,1._wp)
  INTEGER :: j, n1, n2, n3, n4, nm, ncu, info
  REAL(wp) :: dx, dt
  REAL(wp) :: dV, dVI
  REAL(wp) :: R, w, R0, R1
  REAL(wp) :: k, kMax, dk, dkA, dkR
  REAL(wp), DIMENSION(nd-3) :: DL1
  REAL(wp), DIMENSION(nd-2) :: DI1
  REAL(wp), DIMENSION(nd) :: tx, Ve, Vfill, densD
  REAL(wp), DIMENSION(pRes) :: tEr
  REAL(wp), DIMENSION(pMax) :: tE
  REAL(wp), DIMENSION(:), ALLOCATABLE :: DL0, DI0
  COMPLEX(wp), DIMENSION(pRes) :: Lambdagr, Lambdadr
  COMPLEX(wp), DIMENSION(pMax) :: omega, qomega, qomegan, Lambdag, Lambdad
  COMPLEX(wp), DIMENSION(nMax) :: sg, sd
  COMPLEX(wp), DIMENSION(nd-1) :: DL, DU
  COMPLEX(wp), DIMENSION(nd) :: DI
  COMPLEX(wp), DIMENSION(nd+1) :: Sc
  COMPLEX(wp), DIMENSION(nd+1,nd+1) :: Ma
  COMPLEX(wp), DIMENSION(nMax,pRes) :: Gaur, Dter
  COMPLEX(wp), DIMENSION(nMax,pMax) :: Gaud, Dted

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Resolution du probleme de Schrodinger/Poisson instationnaire, coulage du potentiel : methode directe

  SUBROUTINE SOL_INSTAT(dens,d01,Vnm1,V01,psi,ep,z,nf)
    
    !ENTREE
    INTEGER :: nf
    REAL(wp), DIMENSION(nd) :: dens, d01, Vnm1, V01
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep
    COMPLEX(wp), DIMENSION(nd,pMax) :: psi

    !LOCAL
    INTEGER :: n, p
    REAL(wp), DIMENSION(nd) :: Vn, Vdemi, Vtd
    COMPLEX(wp), DIMENSION(pMax) :: psiIg, psiId

    ! Compute the boundary conditions parameters
    psiIg = sg(1)*psi(1,:)*Lambdag - psi(2,:)*(1._wp+Lambdag)
    psiId = sd(1)*psi(nd,:)*Lambdad - psi(nd-1,:)*(1._wp+Lambdad)
    Vdemi = Vnm1

    OPEN(unit=11, file='valeurs/energie')
    OPEN(unit=12, file='valeurs/ecart_dens')
    OPEN(unit=13, file='valeurs/ecart_pot')
    OPEN(unit=14, file='valeurs/puits')

    ! For all time steps
    DO n = 1, nf
       ncu = MIN(n-1,n0)
       Vtd = Ve+Vdemi
       ! Compute the half-time step resonance
       CALL VPNL_DF(Vtd,ep,z)
       k = -kMax
       ! Compute the wave function
       CALL SCHRODI(Lambdag(1),Lambdad(1),psiIg(1),psiId(1),Gaud(:,1),Dted(:,1),Vtd,psi(:,1),n)
       dens = G(k)/2._wp*ABS(psi(:,1))**2
       ! For all the frequency points
       DO p = 2, pMax-1 
          k = k + dk
          ! Compute the wave function
          CALL SCHRODI(Lambdag(p),Lambdad(p),psiIg(p),psiId(p),Gaud(:,p),Dted(:,p),Vtd,psi(:,p),n)
          ! Update the density (sum over the frequency points)
          dens = dens + G(k)*ABS(psi(:,p))**2
       END DO
       ! Compute the wave function
       CALL SCHRODI(Lambdag(pMax),Lambdad(pMax),psiIg(pMax),psiId(pMax),Gaud(:,pMax),Dted(:,pMax),Vtd,psi(:,pMax),n)
       dens = dens + G(kMax)/2._wp*ABS(psi(:,pMax))**2
       dens = dk*m*kb*T/(2._wp*pi**2*hBar**2)*dens
       ! Compute the Poisson potential at the new time step
       CALL POT(dens,Vn)
       ! Compute the half-time step Poisson potential
       Vdemi = 0.5_wp*(3._wp*Vn - Vnm1)
       Vnm1 = Vn
       ! Write in the output files
       WRITE(11,*) n*dt, REAL(z)/e
       WRITE(12,*) n*dt, 100._wp*NORM2(dens(n1:n2)-d01(n1:n2))/NORM2(d01(n1:n2))
       WRITE(13,*) n*dt, 100._wp*NORM2(Vn-V01)/NORM2(V01)
       WRITE(14,*) n*dt, SUM(dens(n1:n2))*dx
    END DO
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    CLOSE(14)

  END SUBROUTINE SOL_INSTAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Resolution du probleme de Schrodinger/Poisson instationnaire, coulage du potentiel : methode de projection
  
  SUBROUTINE SOL_INSTAT_RES(dens,d01,Vnm1,V01,psiE,theta,enm1,z,nf)
    
    !ENTREE
    INTEGER :: nf
    REAL(wp), DIMENSION(nd) :: dens, d01, Vnm1, V01
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: enm1
    COMPLEX(wp), DIMENSION(pMax) :: theta
    COMPLEX(wp), DIMENSION(nd,pRes) :: psiE
    
    !LOCAL
    INTEGER :: n, p, pr
    REAL(wp) :: Inthet
    REAL(wp), DIMENSION(nd) :: Vn, Vdemi, Vtdf
    COMPLEX(wp), DIMENSION(nd) :: en, psiEnm1
    COMPLEX(wp), DIMENSION(pRes) :: psiIg, psiId
     
    ! Compute the wave function
    psiIg = sg(1)*psiE(1,:)*Lambdagr - psiE(2,:)*(1._wp+Lambdagr)
    psiId = sd(1)*psiE(nd,:)*Lambdadr - psiE(nd-1,:)*(1._wp+Lambdadr)
    Vdemi = Vnm1
    en = enm1
    
    OPEN(unit=11, file='valeurs/energieR')
    OPEN(unit=13, file='valeurs/ecart_densR')
    OPEN(unit=14, file='valeurs/ecart_potR')
    OPEN(unit=15, file='valeurs/puitsR')

    ! For all time steps
    DO n = 1, nf
       ncu = MIN(n-1,n0)
       ! Compute the half-time step resonance
       CALL VPNL_DF(Ve+Vdemi,en,z)
       CALL phaseRes(enm1,en)
       Vtdf = Vfill+Vdemi
       k = -kMax
       psiEnm1 = psiE(:,1)
       ! Compute the non resonant part of the wave function
       CALL SCHRODI(Lambdagr(1),Lambdadr(1),psiIg(1),psiId(1),Gaur(:,1),Dter(:,1),Vtdf,psiE(:,1),n)
       dens = dkR/2._wp*G(k)*ABS(psiE(:,1))**2
       ! Compute the projection parameter
       CALL thetaCN(theta(1),en,z,psiEnm1,psiE(:,1),qomega(1),qomegan(1))
       Inthet = dk/2._wp*G(k)*ABS(theta(1))**2
       ! For all the frequency points
       DO p = 2, pMax-1
          k = k + dk
          pr = (p-1)/rap + 1
          ! Compute the non resonant part of the wave function (only for the coarse mesh)
          IF ( p == (pr-1)*rap+1 ) THEN
             psiEnm1 = psiE(:,pr)
             CALL SCHRODI(Lambdagr(pr),Lambdadr(pr),psiIg(pr),psiId(pr),Gaur(:,pr),Dter(:,pr),Vtdf,psiE(:,pr),n)
             dens = dens +  dkR*G(k)*ABS(psiE(:,pr))**2
          END IF
          ! Compute the projection parameter
          CALL thetaCN(theta(p),en,z,psiEnm1,psiE(:,pr),qomega(p),qomegan(p))
          ! Compute additional projection method parameters (interpolation of the non resonant part of the wave function)
          qomegan(p) = qomega(p)*qomegan(p)
          ! Update the resonant part of the density
          Inthet = Inthet + dk*G(k)*ABS(theta(p))**2
       END DO
       psiEnm1 = psiE(:,pRes)
       ! Compute the non resonant part of the wave function
       CALL SCHRODI(Lambdagr(pRes),Lambdadr(pRes),psiIg(pRes),psiId(pRes),Gaur(:,pRes),Dter(:,pRes),Vtdf,psiE(:,pRes),n)
       ! Update the density with non-resonant part of the wave-function
       dens = dens +  dkR/2._wp*G(kMax)*ABS(psiE(:,pRes))**2
       ! Compute the projection parameter
       CALL thetaCN(theta(pMax),en,z,psiEnm1,psiE(:,pRes),qomega(pMax),qomegan(pMax))
       Inthet = Inthet + dk/2._wp*G(kMax)*ABS(theta(pMax))**2
       dens = dens + Inthet*ABS(en)**2 
       dens = m*kb*T/(2._wp*pi**2*hBar**2)*dens
       ! Compute the Poisson potential at the new time step
       CALL POT(dens,Vn)
       ! Compute the half-time step Poisson potential
       Vdemi = 0.5_wp*(3._wp*Vn-Vnm1)
       Vnm1 = Vn
       enm1 = en
       ! Write in the output files
       WRITE(11,*) n*dt, REAL(z)/e
       WRITE(13,*) n*dt, 100._wp*NORM2(dens(n1:n2)-d01(n1:n2))/NORM2(d01(n1:n2))
       WRITE(14,*) n*dt, 100._wp*NORM2(Vn-V01)/NORM2(V01)
       WRITE(15,*) n*dt, SUM(dens(n1:n2))*dx
    END DO
    CLOSE(11)
    CLOSE(13)
    CLOSE(14)
    CLOSE(15)

  END SUBROUTINE SOL_INSTAT_RES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Scale the resonant state with a phase factor
  SUBROUTINE phaseRes(enm1,en)
    
    COMPLEX(wp), DIMENSION(nd) :: enm1, en
    COMPLEX(wp) :: zn
    !zn = SUM(en(1:nd-1)*CONJG(enm1(1:nd-1)))
    zn = SUM(en*CONJG(enm1))
    en = CONJG(zn)/ABS(zn)*en

  END SUBROUTINE phaseRes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Calcul de theta en instationnaire  

  SUBROUTINE thetaCN(th,en,z,psiEnm1,psiEn,qomega,qomegan)
    
    ! ENTREE
    COMPLEX(wp) :: th, z, qomega, qomegan
    COMPLEX(wp), DIMENSION(nd) :: en, psiEnm1, psiEn

    !LOCAL
    COMPLEX(wp) :: an, bn

    an = 0.5_wp*i*dt/hBar*z
    bn =  0.5_wp*i*e*V1*dt/hBar*qomegan*SUM( (psiEnm1(n3:n4)+psiEn(n3:n4)*qomega)*CONJG(en(n3:n4)) )*dx
    th = ( (1._wp-an)*th + bn )/(1._wp+an)

  END SUBROUTINE thetaCN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Resolution du probleme Schrodinger-Poisson stationnaire : methode directe
  SUBROUTINE SOL_STAT_DF(Vs,dens,psi,itMax,tol,cv,it)
   
    !ENTREE
    INTEGER :: itMax, cv, it
    REAL(wp) :: tol
    REAL(wp), DIMENSION(nd) :: Vs, dens
    COMPLEX(wp), DIMENSION(nd,pMax) :: psi
    
    !LOCAL
    REAL(wp) :: erreur
    REAL(wp), DIMENSION(nd) :: VsPrec
    CHARACTER(len=40) :: ch
    
    VsPrec = Vs
    erreur = tol + 1._wp
    it = 0
    cv = 0
    WRITE(ch,*) dV
    OPEN(unit=11,file='valeurs/erreur_B'//trim(adjustl(ch)))
    DO WHILE (erreur.GE.tol)
       IF (it.GE.itMax) THEN
          cv = 1
          EXIT
       END IF
       it = it + 1
       ! Calcul de la densite
       CALL CALC_DENS_DF(dens,Ve+VsPrec,psi)
       ! Calcul du nouveau potentiel
       CALL POT_GN(dens,Vs,VsPrec)
       erreur = NORM2(VsPrec-Vs)/NORM2(Vs)
       WRITE(11,*) it, erreur
       VsPrec = Vs
    END DO
    CLOSE(11)
    
  END SUBROUTINE SOL_STAT_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Resolution du probleme Schrodinger-Poisson stationnaire : maillage adaptatif
  SUBROUTINE SOL_STAT_ADA_DF(Vs,dens,itMax,tol,cv,it)
   
    !ENTREE
    INTEGER :: itMax, cv, it
    REAL(wp) :: tol
    REAL(wp), DIMENSION(nd) :: Vs, dens
    
    !LOCAL
    REAL(wp) :: erreur
    REAL(wp), DIMENSION(nd) :: VsPrec
    CHARACTER(len=40) :: ch
    
    VsPrec = Vs
    erreur = tol + 1._wp
    it = 0
    cv = 0
    WRITE(ch,*) dV
    OPEN(unit=11,file='valeurs/erreurA_B'//trim(adjustl(ch)))
    DO WHILE (erreur.GE.tol)
       IF (it.GE.itMax) THEN
          cv = 1
          EXIT
       END IF
       it = it + 1
       ! Calcul de la densite
       CALL CALC_DENS_ADA_DF(dens,Ve+VsPrec,it)
       ! Calcul du nouveau potentiel
       CALL POT_GN(dens,Vs,VsPrec)
       erreur = NORM2(VsPrec-Vs)/NORM2(Vs)
       WRITE(11,*)it, erreur
       VsPrec = Vs
    END DO
    CLOSE(11)

  END SUBROUTINE SOL_STAT_ADA_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Resolution du probleme Schrodinger-Poisson stationnaire : methode de projection
  SUBROUTINE SOL_STAT_RESO_DF(Vs,dens,psiE,theta,ep,z,itMax,tol,cv,it)
   
    !ENTREE
    INTEGER :: itMax, cv, it
    REAL(wp) :: tol
    REAL(wp), DIMENSION(nd) :: Vs, dens
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep
    COMPLEX(wp), DIMENSION(pMax) :: theta
    COMPLEX(wp), DIMENSION(nd,pRes) :: psiE
    
    !LOCAL
    REAL(wp) :: erreur
    REAL(wp), DIMENSION(nd) :: VsPrec
    CHARACTER(len=40) :: ch

    CALL VPDIR_SX(Vs,ep,z)   
    !Resolution du probleme stationnaire, dV = 0
    VsPrec = Vs
    erreur = tol + 1._wp
    it = 0
    cv = 0
    WRITE(ch,*) dV
    OPEN(unit=11,file='valeurs/erreurR_B'//trim(adjustl(ch)))
    DO WHILE (erreur.GE.tol)
       IF (it.GE.itMax) THEN
          cv = 1
          EXIT
       END IF
       it = it + 1
       ! Calcul de la densite
       CALL CALC_DENS_RESO_DF(dens,Ve+VsPrec,Vfill+VsPrec,psiE,theta,ep,z)
       ! Calcul du nouveau potentiel
       CALL POT_GN(dens,Vs,VsPrec)
       erreur = NORM2(VsPrec-Vs)/NORM2(Vs)
       WRITE(11,*)it, erreur
       VsPrec = Vs
    END DO
    CLOSE(11)
  END SUBROUTINE SOL_STAT_RESO_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Resolution du probleme Schrodinger-Poisson stationnaire : methode de projection (integration explicite de la partie resonante)
  SUBROUTINE SOL_STAT_RESO_INT_DF(Vs,dens,psiE,ep,z,itMax,tol,cv,it)
   
    !ENTREE
    INTEGER :: itMax, cv, it
    REAL(wp) :: tol
    REAL(wp), DIMENSION(nd) :: Vs, dens
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep
    COMPLEX(wp), DIMENSION(nd,pRes) :: psiE
    
    !LOCAL
    REAL(wp) :: erreur
    REAL(wp), DIMENSION(nd) :: VsPrec
    CHARACTER(len=40) :: ch

    !Resolution du probleme stationnaire, dV = 0
    CALL VPDIR_SX(Vs,ep,z)
    VsPrec = Vs
    erreur = tol + 1._wp
    it = 0
    cv = 0
    WRITE(ch,*) dV
    OPEN(unit=11,file='valeurs/erreurRI_B'//trim(adjustl(ch)))
    DO WHILE (erreur.GE.tol)
       IF (it.GE.itMax) THEN
          cv = 1
          EXIT
       END IF
       it = it + 1
       ! Calcul de la densite
       CALL CALC_DENS_RESO_INT_DF(dens,Ve+VsPrec,Vfill+VsPrec,psiE,ep,z)
       ! Calcul du nouveau potentiel
       CALL POT_GN(dens,Vs,VsPrec)
       erreur = NORM2(VsPrec-Vs)/NORM2(Vs)
       WRITE(11,*)it, erreur
       VsPrec = Vs
    END DO
    CLOSE(11)
  END SUBROUTINE SOL_STAT_RESO_INT_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcul de la densite de probabilite a partir du potentiel V, la somme en frequence est realisee de -kMax a kMax, le maillage en frequence est fixe

  SUBROUTINE CALC_DENS_DF(dens,V,psi)
    
    !ENTREE
    REAL(wp), DIMENSION(nd) :: dens, V
    COMPLEX(wp), DIMENSION(nd,pMax) :: psi
  
    ! LOCAL
    INTEGER :: p

    ! Initialisation
    k = -kMax
    CALL SCHROG_DF(V,k,psi(:,1))
    dens = G(k)/2._wp*ABS(psi(:,1))**2
    ! Integration en k
    DO p = 2, pMax-1
       k = k + dk
       ! Resolution Schrodinger
       IF (k.LT.0_wp) THEN
          CALL SCHROG_DF(V,k,psi(:,p))
       ELSE 
          CALL SCHROD_DF(V,k,psi(:,p))
       END IF
       ! Incrementation de l integration pour chaque x
       dens = dens + G(k)*ABS(psi(:,p))**2
    END DO
    CALL SCHROD_DF(V,kMax,psi(:,pMax))
    dens = dens + G(kMax)/2._wp*ABS(psi(:,pMax))**2
    ! Renormalisation
    dens = dk*m*kb*T/(2._wp*pi**2*hBar**2)*dens
    
  END SUBROUTINE CALC_DENS_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcul de la densite de probabilite a partir du potentiel V, la somme en frequence est realisee de -kMax a kMax, le maillage en frequence est adaptatif

  SUBROUTINE CALC_DENS_ADA_DF(dens,V,pl)
    
    !ENTREE
    INTEGER :: pl
    REAL(wp), DIMENSION(nd) :: dens, V

    ! LOCAL
    INTEGER :: p
    REAL(wp) :: dkpm1, dkp, Der
    COMPLEX(wp) :: psipm1
    COMPLEX(wp), DIMENSION(nd) :: psi

    ! Initialisation
    p = 1
    k = -kMax
    CALL SCHROG_DF(V,k,psi)
    dkpm1 = dkA
    dens = dkpm1/2._wp*G(k)*ABS(psi)**2
    DO WHILE ( k + dkpm1 .LE. kMax - dkA )
       p = p + 1
       k = k + dkpm1
       IF (p.GE.5000) THEN
          WRITE(*,*) k, kMax
          WRITE(*,*)'Integration inachevee, iteration ', pl
          STOP
       END IF
       ! Resolution Schrodinger
       IF (k.LT.0_wp) THEN
          psipm1 = psi(1)
          CALL SCHROG_DF(V,k,psi)
          CALL DER_G(Der,psipm1,psi(1),k,dkpm1)
       ELSE 
          psipm1 = psi(nd)
          CALL SCHROD_DF(V,k,psi)
          CALL DER_D(Der,psipm1,psi(nd),k,dkpm1)
       END IF
       IF ( Der.GT.1.3_wp*L ) THEN
          dkp = dk/1.2_wp
       ELSE
          dkp = dkA
       END IF
       dens = dens + (dkpm1+dkp)/2._wp*G(k)*ABS(psi)**2
       dkpm1 = dkp
    END DO
    k = k + dkpm1
    CALL SCHROD_DF(V,k,psi)
    dens = dens + dkpm1/2._wp*G(k)*ABS(psi)**2
    ! Renormalisation
    dens = m*kb*T/(2._wp*pi**2*hBar**2)*dens
    
  END SUBROUTINE CALC_DENS_ADA_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcul de la densite de probabilite a partir du potentiel Vs en decomposant la fonction propre en une partie non resonante psiExt et une partie resonante psiInt, passe par le calcul de l etat resonant. Permet un maillage large en frequence.

  SUBROUTINE CALC_DENS_RESO_DF(dens,Vt,Vtf,psiE,theta,ep,z)
    
    !ENTREE
    REAL(wp), DIMENSION(nd) :: dens, Vt, Vtf
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep
    COMPLEX(wp), DIMENSION(pMax) :: theta
    COMPLEX(wp), DIMENSION(nd,pRes) :: psiE
  
    ! LOCAL
    INTEGER :: p, pr
    REAL(wp) :: Inthet
    
    ! Calcul de l etat resonant
    CALL VPNL_DF(Vt,ep,z)

    ! Initialisation
    k = -kMax
    CALL SCHROG_DF(Vtf,k,psiE(:,1))
    dens = dkR/2._wp*G(k)*ABS(psiE(:,1))**2
    theta(1) = e*V1*SUM(psiE(n3:n4,1)*CONJG(ep(n3:n4)))*dx/(z-tE(1))
    Inthet = dk/2._wp*G(k)*ABS(theta(1))**2
    ! Integration en k
    DO p = 2, pMax-1 
       k = k + dk
       pr = (p-1)/rap + 1
       IF ( p == (pr-1)*rap+1 ) THEN
          IF ( k.LT.0_wp ) THEN
             CALL SCHROG_DF(Vtf,k,psiE(:,pr))
          ELSE
             CALL SCHROD_DF(Vtf,k,psiE(:,pr))
          END IF
          dens = dens + dkR*G(k)*ABS(psiE(:,pr))**2
       END IF
       theta(p) = e*V1*SUM(psiE(n3:n4,pr)*CONJG(ep(n3:n4)))*dx/(z-tE(p))
       Inthet = Inthet + dk*G(k)*ABS(theta(p))**2
    END DO
    CALL SCHROD_DF(Vtf,kMax,psiE(:,pRes))
    dens = dens + dkR/2._wp*G(kMax)*ABS(psiE(:,pRes))**2
    theta(pMax) = e*V1*SUM(psiE(n3:n4,pRes)*CONJG(ep(n3:n4)))*dx/(z-tE(pMax))
    Inthet = Inthet + dk/2._wp*G(kMax)*ABS(theta(pMax))**2
    ! Renormalisation
    dens = dens + Inthet*ABS(ep)**2
    dens = m*kb*T/(2._wp*pi**2*hBar**2)*dens
   
  END SUBROUTINE CALC_DENS_RESO_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcul de la densite de probabilite a partir du potentiel Vs en decomposant la fonction propre en une partie non resonante psiExt et une partie resonante psiInt, passe par le calcul de l etat resonant. Permet un maillage large en frequence. Ici l integrale de la partie resonnante est calculee explicitement

  SUBROUTINE CALC_DENS_RESO_INT_DF(dens,Vt,Vtf,psiE,ep,z)
    
    !ENTREE
    REAL(wp), DIMENSION(nd) :: dens, Vt, Vtf
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep
    COMPLEX(wp), DIMENSION(nd,pRes) :: psiE
  
    ! LOCAL
    INTEGER :: p
    REAL(wp) :: Er, Gamma, Inthet
    COMPLEX(wp) :: Sp

    ! Calcul de l etat resonant
    CALL VPNL_DF(Vt,ep,z)
    Er = REAL(z)
    Gamma = -AIMAG(z)

    ! Initialisation
    k = -kMax
    CALL SCHROG_DF(Vtf,k,psiE(:,1))
    dens = dkR/2._wp*G(k)*ABS(psiE(:,1))**2
    Sp = ABS(e*V1*SUM(psiE(n3:n4,1)*CONJG(ep(n3:n4)))*dx)**2
    Inthet = m*G(k)*Sp*(ATAN((tEr(2)-Er)/Gamma)-ATAN((tEr(1)-Er)/Gamma))/(hBar**2*k*Gamma)
    ! Integration en k
    DO p = 2, pRes-1
       ! Resolution Schrodinger
       k = k + dkR
       IF (k.LT.0_wp) THEN
          CALL SCHROG_DF(Vtf,k,psiE(:,p))
       ELSE
          CALL SCHROD_DF(Vtf,k,psiE(:,p))
       END IF
       dens = dens + dkR*G(k)*ABS(psiE(:,p))**2
       Sp = ABS(e*V1*SUM(psiE(n3:n4,p)*CONJG(ep(n3:n4)))*dx)**2
       Inthet=Inthet+m*G(k)*Sp*(ATAN((tEr(p+1)-Er)/Gamma)-ATAN((tEr(p)-Er)/Gamma))/(hBar**2*k*Gamma)
    END DO
    CALL SCHROD_DF(Vtf,kMax,psiE(:,pRes))
    dens = dens + dkR/2._wp*G(kMax)*ABS(psiE(:,pRes))**2
    ! Renormalisation
    dens = dens + Inthet*ABS(ep)**2
    dens = m*kb*T/(2._wp*pi**2*hBar**2)*dens
   
  END SUBROUTINE CALC_DENS_RESO_INT_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcul direct du potentiel : resolution de l equation de poisson a partir de la densite

  SUBROUTINE POT(dens,V)
    
    ! ENTREE
    REAL(wp), DIMENSION(nd) :: dens, V
    
    ! LOCAL
    REAL(wp), DIMENSION(nd-2) :: U
    
    U = R1*(densD(2:nd-1)-dens(2:nd-1))
    DI1 = 2._wp
    DL1 = -1._wp
    CALL DPTSV(nd-2,1,DI1,DL1,U,nd-2,info)
    V(1) = 0._wp
    V(2:nd-1) = U
    V(nd) = 0._wp
    
  END SUBROUTINE POT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Compute the potential corresponding to the Gummel iteration
  
  SUBROUTINE POT_GN(dens,V,W)
    
    ! ENTREE
    REAL(wp), DIMENSION(nd) :: dens, V, W
    
    ! LOCAL
    INTEGER :: cpt, cptMax
    REAL(wp) :: erreur, tol 
    REAL(wp), DIMENSION(nd-2) :: U, Q, ne
    ! Newton
    cptMax = 100
    cpt = 0
    tol = 1.E-15_wp
    erreur = tol + 1._wp
    U = W(2:nd-1)
    DO WHILE(erreur.GE.tol)
       cpt = cpt + 1
       ne = dens(2:nd-1)*EXP((U - W(2:nd-1))/Vref)
       DL1 = -1._wp
       DI1 = 2._wp + R1/Vref*ne
       Q(1) = 2._wp*U(1)-U(2) + R1*(ne(1)-densD(2))
       DO j = 2, nd-3
          Q(j) = -U(j-1)+2._wp*U(j)-U(j+1) + R1*(ne(j)-densD(j+1))
       END DO
       Q(nd-2) = -U(nd-3)+2._wp*U(nd-2) + R1*(ne(nd-2)-densD(nd-1))
       CALL DPTSV(nd-2,1,DI1,DL1,Q,nd-2,info)
       erreur = NORM2(Q)
       U = U - Q
       IF (cpt.GE.cptMax) THEN
          WRITE(*,*)'Gummel, Non convergence Newton, iteration', cpt
          STOP
       END IF
    END DO
    V(1) = 0._wp
    V(2:nd-1) = U
    V(nd) = 0._wp
  END SUBROUTINE POT_GN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Resolution de l equation de Schrodinger avec la methode des Differences Finies pour k positif

  SUBROUTINE SCHROD_DF(V,k,psi)

    ! ENTREE
    REAL(wp) :: k
    REAL(wp), DIMENSION(nd) :: V
    COMPLEX(wp), DIMENSION(nd) :: psi

    ! LOCAL
    REAL(wp) :: el
    COMPLEX(wp) :: lambda, beta

    DL = -1._wp
    DU = -1._wp
    psi = 0._wp
    el = k**2*dx**2
    lambda = 1._wp - el/2._wp + i*SQRT(el-el**2/4._wp)
    el = (k**2+2._wp*m*e/hBar**2*dV)*dx**2
    beta = 1._wp - el/2._wp + i*SQRT(el-el**2/4._wp)
    DI(1) = 1._wp/lambda
    psi(1) = 1._wp/lambda - lambda
    DI(2:nd-1) = 2._wp + dx**2*(-2._wp*m*e/hBar**2*V(2:nd-1)-k**2)
    DI(nd) = 1._wp/beta
    CALL ZGTSV(nd,1,DL,DI,DU,psi,nd,info)
    
  END SUBROUTINE SCHROD_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Resolution de l equation de Schrodinger avec la methode des Differences Finies pour k negatif
  
  SUBROUTINE SCHROG_DF(V,k,psi)

    ! ENTREE
    REAL(wp) :: k
    REAL(wp), DIMENSION(nd) :: V
    COMPLEX(wp), DIMENSION(nd) :: psi
    
    ! LOCAL
    REAL(wp) :: el
    COMPLEX(wp) :: lambda, beta
    
    DL = -1._wp
    DU = -1._wp
    psi = 0._wp
    el = k**2*dx**2
    lambda = 1._wp - el/2._wp - i*SQRT(el-el**2/4._wp)
    el = (k**2-2._wp*m*e/hBar**2*dV)*dx**2
    beta = 1._wp - el/2._wp - i*SQRT( CMPLX(el-el**2/4._wp,0._wp) )
    DI(1) = beta
    DI(2:nd-1) = 2._wp + dx**2*(2._wp*m*e/hBar**2*(-V(2:nd-1)+dV)-k**2)
    DI(nd) = lambda
    psi(nd) = lambda - 1._wp/lambda
    CALL ZGTSV(nd,1,DL,DI,DU,psi,nd,info)

  END SUBROUTINE SCHROG_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Resolution de l equation de Schrodinger instationnaire en utilisant les coefficients sg, sd permettant de definir des conditions aux limites transparentes
  
  SUBROUTINE SCHRODI(Lg,Ld,psiIg,psiId,Gau,Dte,V,psiPrec,n)
    ! ENTREE
    INTEGER :: n
    REAL(wp), DIMENSION(nd) :: V
    COMPLEX(wp) :: Lg, Ld, psiIg, psiId
    COMPLEX(wp), DIMENSION(nd) :: psiPrec
    COMPLEX(wp), DIMENSION(nMax) :: Gau, Dte
    ! LOCAL
    COMPLEX(wp), DIMENSION(nd) :: psiSuiv

    DL = -1._wp
    DU = -1._wp
    DI(1) = sg(1)
    DI(2:nd-1) = 2._wp -i*R - w*V(2:nd-1)
    DI(nd) = sd(1)
    psiSuiv(1) = psiPrec(2) + SUM(sg(ncu+1:2:-1)*Gau(n-ncu+1:n)) + psiIg
    psiSuiv(nd) = psiPrec(nd-1) + SUM(sd(ncu+1:2:-1)*Dte(n-ncu+1:n)) + psiId
    DO j = 2, nd-1
       psiSuiv(j) = psiPrec(j-1) + ( - 2._wp -i*R + w*V(j) )*psiPrec(j) + psiPrec(j+1)
    END DO
    CALL ZGTSV(nd,1,DL,DI,DU,psiSuiv,nd,info)
    psiPrec = psiSuiv
    Gau(n+1) = Gau(n+1) - psiPrec(1)
    Dte(n+1) = Dte(n+1) - psiPrec(nd)
    psiIg = psiIg*Lg
    psiId = psiId*Ld
  END SUBROUTINE SCHRODI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcul des coefficients permettant de definir les conditions aux limites transparentes
  
  SUBROUTINE SN(R,V,s,nMax)
    ! ENTREE
    INTEGER :: nMax
    REAL(wp) :: R, V
    ! SORTIE
    COMPLEX(wp), DIMENSION(nMax) :: s
    !LOCAL
    REAL(wp) :: Mu, Phi, Sigma
    COMPLEX(wp) :: Alpha, Lambda
    ! Calcul des parametres
    Sigma = 2._wp*e*m*dx**2/hBar**2*V
    Mu = (R**2+4._wp*Sigma+Sigma**2)/SQRT((R**2+Sigma**2)*(R**2+(Sigma+4._wp)**2))
    Phi = ATAN( 2._wp*R*(Sigma+2._wp)/(R**2-4._wp*Sigma-Sigma**2) )
    Alpha = i/2._wp*((R**2+Sigma**2)*(R**2+(Sigma+4._wp)**2))**(0.25_wp)*EXP(i*Phi/2._wp) 
    Lambda = EXP(i*Phi)
    ! Calcul des coefficients
    s(1) = 1._wp - i*R/2._wp + Sigma/2._wp - Alpha
    s(2) = 1._wp + i*R/2._wp + Sigma/2._wp + Alpha*Mu/Lambda
    s(3) = Alpha*(Mu**2-1._wp)/(2._wp*Lambda**2)
    DO j = 2, nMax-2
       s(j+2) = (2._wp*REAL(j)-1._wp)/(REAL(j)+1._wp)*Mu/Lambda*s(j+1) - (REAL(j)-2._wp)/(REAL(j)+1._wp)/Lambda**2*s(j)
    END DO
  END SUBROUTINE SN
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Calcul du vecteur propre avec conditions de Dirichlet sur une partie du domaine
  SUBROUTINE VPDIR_PT(Vs,ep,z)
    
    !ENTREE
    REAL(wp), DIMENSION(nd) :: Vs
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep

    ! LOCAL
    REAL(wp) :: Er
    REAL(wp), DIMENSION(nm) :: Vloc
    REAL(wp), DIMENSION(4*nm) :: wo
    REAL(wp), DIMENSION(nm,nm) :: S

    ep = 0._wp
    DL0 = -1._wp
    Vloc = -e*(Ve(n1+1:n2-1)+Vs(n1+1:n2-1))
    DI0 = 2._wp + 2._wp*m*dx**2/hBar**2*(Vloc+ABS(MINVAL(Vloc)))
    CALL DPTEQR('I',nm,DI0,DL0,S,nm,wo,info)
    ep(n1+1:n2-1) = S(:,nm)
    Er = DI0(nm)*hBar**2/(2._wp*m*dx**2) - ABS(MINVAL(Vloc))
    z = CMPLX(Er,0._wp)

  END SUBROUTINE VPDIR_PT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calcul du vecteur propre avec conditions de Dirichlet sur une partie du domaine
  SUBROUTINE VPDIR_SX(Vs,ep,z)
    
    !ENTREE
    REAL(wp), DIMENSION(nd) :: Vs
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep

    ! LOCAL
    INTEGER :: nb
    INTEGER, DIMENSION(nm) :: ifail
    INTEGER, DIMENSION(5*nm) :: iwo
    REAL(wp) :: Er, vl, vu
    REAL(wp), DIMENSION(nm) :: Vloc, ev, S
    REAL(wp), DIMENSION(5*nm) :: wo

    ep = 0._wp
    DL0 = -1._wp
    Vloc = -e*(Ve(n1+1:n2-1)+Vs(n1+1:n2-1))
    DI0 = 2._wp + 2._wp*m*dx**2/hBar**2*(Vloc+ABS(MINVAL(Vloc))) 
    CALL DSTEVX('V','I',nm,DI0,DL0,vl,vu,1,1,0._wp,nb,ev,S,nm,wo,iwo,ifail,info)
    ep(n1+1:n2-1) = S
    Er = ev(1)*hBar**2/(2._wp*m*dx**2) - ABS(MINVAL(Vloc))
    z = CMPLX(Er,0._wp)
    
  END SUBROUTINE VPDIR_SX
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calcul du vecteur propre avec conditions de Dirichlet sur une partie du domaine
  
  SUBROUTINE VPDIR_SR(Vs,ep,z)
    
    !ENTREE
    REAL(wp), DIMENSION(nd) :: Vs
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: ep

    ! LOCAL
    INTEGER :: nb
    INTEGER, DIMENSION(2) :: isz
    INTEGER, DIMENSION(10*nm) :: iwo
    REAL(wp) :: Er, vl, vu
    REAL(wp), DIMENSION(nm) :: Vloc, ev, S
    REAL(wp), DIMENSION(20*nm) :: wo

    ep = 0._wp
    DL0 = -1._wp
    Vloc = -e*(Ve(n1+1:n2-1)+Vs(n1+1:n2-1))
    DI0 = 2._wp + 2._wp*m*dx**2/hBar**2*(Vloc+ABS(MINVAL(Vloc)))
    CALL DSTEVR('V','I',nm,DI0,DL0,vl,vu,1,1,0._wp,nb,ev,S,nm,isz,wo,20*nm,iwo,10*nm,info)
    ep(n1+1:n2-1) = S
    Er = ev(1)*hBar**2/(2._wp*m*dx**2) - ABS(MINVAL(Vloc))
    z = CMPLX(Er,0._wp)
    
  END SUBROUTINE VPDIR_SR
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Recherche de l etat resonant : methode de Newton pour resoudre le probleme de vecteur propre non lineaire M(z)*u = 0

  SUBROUTINE VPNL_DF(V,u,z)

    !ENTREE
    REAL(wp), DIMENSION(nd) :: V
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: u

    !LOCAL
    INTEGER :: it, itMax
    INTEGER, DIMENSION(nd+1) :: ipiv  
    REAL(wp) :: tol, erreur
    
    itMax = 10
    it = 0
    tol = 1.E-15_wp

    u = u/SQRT(SUM(u*CONJG(u)))
    CALL CALC_FDF_DF(V,u,z)
    erreur = NORM2(ABS(Sc))
    DO WHILE (erreur.GE.tol)
       it = it + 1
       CALL ZGESV(nd+1,1,Ma,nd+1,ipiv,Sc,nd+1,info)
       u = u - Sc(1:nd)
       z = z - Sc(nd+1)
       CALL CALC_FDF_DF(V,u,z)
       erreur = NORM2(ABS(Sc))
       IF (it.GE.itMax) THEN
          WRITE(*,*) 'VPNL, non convergence Newton, it =', it
          WRITE(*,*) 'erreur =', erreur
          EXIT
       END IF
    END DO
    u = u/SQRT(SUM(u*CONJG(u))*dx)
    
  END SUBROUTINE VPNL_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Recherche de l etat resonant : Calcul de la differetielle de la fonction a annuler pour la methode de Newton
  
  SUBROUTINE CALC_FDF_DF(V,u,z)
    
    !ENTREE
    REAL(wp), DIMENSION(nd) :: V
    COMPLEX(wp) :: z
    COMPLEX(wp), DIMENSION(nd) :: u

    !LOCAL
    COMPLEX(wp) :: el
    COMPLEX(wp) :: beta, alpha, betap, alphap

    el = dx**2*m/hBar**2*z
    beta = 1._wp - el - i*SQRT(2._wp*el-el**2) 
    betap = -dx**2*m/hBar**2*( 1._wp + i*(1._wp - el)/SQRT(2._wp*el-el**2)  )
    el = dx**2*m/hBar**2*(z+e*dV)
    alpha = 1._wp - el - i*SQRT(2._wp*el-el**2) 
    alphap = -dx**2*m/hBar**2*( 1._wp + i*(1._wp - el)/SQRT(2._wp*el-el**2)  )
    Ma = 0._wp
    Ma(1,1) = beta
    Ma(1,2) = -1._wp
    Ma(1,nd+1) = betap*u(1)
    Sc(1) = beta*u(1)-u(2)
    DO j = 2, nd-1
       Ma(j,j-1) = -1._wp
       Ma(j,j) = 2._wp + 2._wp*m*dx**2/hBar**2*(-e*V(j)-z)
       Ma(j,j+1) = -1._wp
       Ma(j,nd+1) = -2._wp*m*dx**2/hBar**2*u(j)
       Sc(j) = -u(j-1)+Ma(j,j)*u(j)-u(j+1)
    END DO
    Ma(nd,nd-1) = -1._wp
    Ma(nd,nd) = alpha
    Ma(nd,nd+1) = alphap*u(nd)
    Sc(nd) = -u(nd-1)+alpha*u(nd)
    Ma(nd+1,1:nd) = CONJG(u)
    Sc(nd+1) = 0._wp

  END SUBROUTINE CALC_FDF_DF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcule la barriere de potentiel Ve
  FUNCTION CVe(x) RESULT(Pext)
    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: Pext
    Pext = 0._wp
    IF ( ((x.GE.bar(2)).AND.(x.LT.bar(3))).OR.((x.GT.bar(4)).AND.(x.LE.bar(5))) ) THEN
       Pext = -V1
    ELSE IF (x.GE.bar(6)) THEN
       Pext = dV
    END IF
    IF ( (x.GE.bar(1)).AND.(x.LT.bar(6)) ) THEN
       Pext = Pext + dV*(x-bar(1))/(bar(6)-bar(1))
    END IF
  END FUNCTION CVe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calcule le potentiel rempli Vfill
  FUNCTION CVfill(x) RESULT(Pfill)
    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: Pfill
    Pfill = 0._wp
    IF ( (x.GE.bar(2)).AND.(x.LE.bar(5)) ) THEN
       Pfill = -V1
    ELSE IF (x.GE.bar(6)) THEN
       Pfill = dV
    END IF
    IF ( (x.GE.bar(1)).AND.(x.LT.bar(6)) ) THEN
       Pfill = Pfill + dV*(x-bar(1))/(bar(6)-bar(1))
    END IF
  END FUNCTION CVfill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Statistique de fermi

  FUNCTION G(k) RESULT(fermi)
    REAL(wp), INTENT(IN) :: k
    REAL(wp) :: fermi
    fermi = LOG( 1._wp + EXP( (Ef-hBar**2*k**2/(2._wp*m))/(kb*T) ) )
  END FUNCTION G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcule la derivee du coefficient de transmition pour p positif  
  SUBROUTINE DER_D(Der,psiPrec,psiSuiv,k,dkloc)
    ! ENTREE
    COMPLEX(wp) :: psiPrec, psiSuiv
    REAL(wp) :: k, dkloc, Der
    Der = ABS( -1._wp/k + k/(k**2+2._wp*e*m/hBar**2*dV) + 2._wp*REAL(CONJG(psiSuiv)*(psiSuiv-psiPrec)/dkloc)/ABS(psiSuiv)**2 )
  END SUBROUTINE DER_D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! Calcule le coefficient de transmition et sa derivee pour p negatif
  SUBROUTINE DER_G(Der,psiPrec,psiSuiv,k,dkloc)
    ! ENTREE
    COMPLEX(wp) :: psiPrec, psiSuiv
    REAL(wp) :: k, Der, dkloc
    !LOCAL
    REAL(wp) :: Qte
    Qte = k**2 - 2._wp*e*m/hBar**2*dV
    IF ( Qte.GT.0._wp) THEN
       Der = ABS(1._wp/k - k/Qte - 2._wp*REAL(CONJG(psiSuiv)*(psiSuiv-psiPrec)/dkloc)/ABS(psiSuiv)**2) 
    ELSE
       Der = 0._wp
    END IF
  END SUBROUTINE DER_G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE fonctions
