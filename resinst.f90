! Main file. Calls the time dependent Schrodinger Poisson solvers (file fonctions.f90): the projection method and the direct resolution. Compares the efficiencies
PROGRAM resinst

  USE fonctions

  IMPLICIT NONE

  INTEGER :: n, nf, nf0, itMax, cv, it, p, pr
  REAL(wp) :: dVMax, tol, debut, fin, dt0
  REAL(wp), DIMENSION(nd) :: V0, V01, d0, d01
  REAL(wp), DIMENSION(nd) :: VR0, VR01, dR0, dR01
  REAL(wp), DIMENSION(pRes) :: wgr, wdr
  REAL(wp), DIMENSION(pMax) :: wg, wd, w0
  COMPLEX(wp) :: znm1
  COMPLEX(wp), DIMENSION(nd) :: enm1
  COMPLEX(wp), DIMENSION(pMax) :: thetanm1
  COMPLEX(wp), DIMENSION(nd,pRes) :: psiEnm1
  COMPLEX(wp), DIMENSION(nd,pMax) :: psinm1
  CHARACTER(len=40) :: ch
  
  !Open the information output file
  OPEN(unit=111, file='valeurs/resultats')
  WRITE(111,*)'** Info'
  WRITE(111,*)
  WRITE(111,*)'pRes', pRes
  WRITE(111,*)'pMax', pMax
  WRITE(111,*)

  !! Compute the numerical parameters
  
  ! Frequency mesh parameters
  kMax = SQRT(2._wp*m*(Ef+7._wp*kb*T))/hBar
  dk = 2._wp*kMax/REAL(pMax-1)
  dkR = 2._wp*kMax/REAL(pRes-1)
  dkA = 10._wp*dk
  ! Time and space discretization parameters
  dx = L/(REAL(nd)-1._wp)
  nf0 = 10000
  dt0 = 1.E-15_wp
  dt = 1.E-15_wp

  dVMax = 0.1_wp

  ! Validate the space discretization size (compare to the characteristic wave length)
  IF (dx.GE.0.5_wp/SQRT(kMax**2+2._wp*m*e/hBar**2*dVMax)) THEN
     WRITE(*,*) 'dx trop grand'
     STOP
  END IF
  
  ! Define the external potential parameters
  n1 = FLOOR(bar(2)/dx) + 2
  n2 = FLOOR(bar(5)/dx) + 1
  nm = n2 - n1 - 1
  ALLOCATE(DL0(nm-1))
  ALLOCATE(DI0(nm))
  n3 = FLOOR(bar(3)/dx) + 2
  n4 = FLOOR(bar(4)/dx) + 1
  IF ( bar(3)/dx-FLOOR(bar(3)/dx) .LT. 1.E-14_wp ) THEN
     WRITE(*,*) 'probleme bar3'
     STOP
  END IF

  ! Additional discretization parameters
  R0 = hBar**2/(2._wp*m*dx**2)
  R1 = e*dx**2/epsilon
  
  tx = (/ ( REAL(j)*dx , j = 0, nd-1 ) /)
  dV = 0._wp

  !Compute the external potential (given)
  densD = densD1
  DO j = 1, nd 
     Ve(j) = CVe(tx(j))
     Vfill(j) = CVfill(tx(j))
     IF ( (tx(j).GE.bar(1)).AND.(tx(j).LE.bar(6)) ) THEN
        densD(j) = densD2
     END IF
  END DO

  ! Compute the frequency mesh
  k = -kMax
  tE(1) = hBar**2*k**2/(2._wp*m) - e*dV
  DO p = 2, pMax
     k = k + dk
     IF (k.LT.0_wp) THEN
        tE(p) = hBar**2*k**2/(2._wp*m)-e*dV
     ELSE
        tE(p) = hBar**2*k**2/(2._wp*m)
     END IF
  END DO
  k = -kMax
  tEr(1) = hBar**2*k**2/(2._wp*m) - e*dV
  DO p = 2, pRes
     k = k + dkR
     IF (k.LT.0_wp) THEN
        tEr(p) = hBar**2*k**2/(2._wp*m)-e*dV
     ELSE
        tEr(p) = hBar**2*k**2/(2._wp*m)
     END IF
  END DO

  !Open the output files for the potential
  WRITE(ch,*) dV
  OPEN(unit=11,file='valeurs/pot_B'//trim(adjustl(ch)))
  OPEN(unit=12,file='valeurs/potR_B'//trim(adjustl(ch)))
  DO j = 1, nd
     READ(11,*) V0(j)
     READ(12,*) VR0(j)
  END DO
  CLOSE(11)
  CLOSE(12)

  !Compute the initial densities by solving the time-independent Schrodinger Poisson problem
  CALL CALC_DENS_DF(d0,Ve+V0,psinm1)
  CALL VPDIR_SX(VR0,enm1,znm1)
  CALL CALC_DENS_RESO_DF(dR0,Ve+VR0,Vfill+VR0,psiEnm1,thetanm1,enm1,znm1)

  !Compute the parameters related to the boundary conditions
  R = 4._wp*m*dx**2/(hBar*dt)
  w = 2._wp*e*m*dx**2/hBar**2

  !! Apply a potential shift at initial time
  ! Update the potential bias
  dVI = dV
  dV = 0.1_wp

  !Compute additional parameters related to the boundary conditions
  wg = 2._wp*ATAN(tE*dt/(2._wp*hBar))/dt
  Lambdag = EXP(-i*wg*dt)
  wgr = 2._wp*ATAN(tEr*dt/(2._wp*hBar))/dt
  Lambdagr = EXP(-i*wgr*dt)
  wd = 2._wp*ATAN((tE+e*(dVI-dV))*dt/(2._wp*hBar))/dt
  Lambdad = EXP(-i*wd*dt)
  wdr = 2._wp*ATAN((tEr+e*(dVI-dV))*dt/(2._wp*hBar))/dt
  Lambdadr = EXP(-i*wdr*dt)
  Gaud(1,:) = psinm1(1,:)
  Dted(1,:) = psinm1(nd,:)
  Gaur(1,:) = psiEnm1(1,:)
  Dter(1,:) = psiEnm1(nd,:)
  DO n = 2, nMax
     Gaud(n,:) = Lambdag*Gaud(n-1,:)
     Dted(n,:) = Lambdad*Dted(n-1,:)
     Gaur(n,:) = Lambdagr*Gaur(n-1,:)
     Dter(n,:) = Lambdadr*Dter(n-1,:)
  END DO

  ! Update the external potential corresponding to the new potential bias
  DO j = 1, nd 
     Ve(j) = CVe(tx(j))
     Vfill(j) = CVfill(tx(j))
  END DO

  ! Update the frequency mesh corresponding to the new potential bias
  k = -kMax
  tE(1) = hBar**2*k**2/(2._wp*m) - e*dV
  DO p = 2, pMax
     k = k + dk
     IF (k.LT.0_wp) THEN
        tE(p) = hBar**2*k**2/(2._wp*m) - e*dV
     ELSE
        tE(p) = hBar**2*k**2/(2._wp*m)
     END IF 
  END DO
  !w0 = tE/hBar
  w0 = 2._wp*ATAN(tE*dt/(2._wp*hBar))/dt
  omega = EXP(i*w0*dt)
  !omega = 1._wp
  DO p = 1, pMax
     pr = (p-1)/rap
     qomega(p) = omega(pr*rap+1)/omega(p)
  END DO
  qomegan = 1._wp
  k = -kMax
  tEr(1) = hBar**2*k**2/(2._wp*m) - e*dV
  DO p = 2, pRes
     k = k + dkR
     IF (k.LT.0_wp) THEN
        tEr(p) = hBar**2*k**2/(2._wp*m) - e*dV
     ELSE
        tEr(p) = hBar**2*k**2/(2._wp*m)
     END IF
  END DO
  
  ! Read the stationnary potential corresponding to the new potential bias
  WRITE(ch,*) dV
  OPEN(unit=11,file='valeurs/pot_B'//trim(adjustl(ch)))
  OPEN(unit=12,file='valeurs/dens_B'//trim(adjustl(ch)))
  OPEN(unit=13,file='valeurs/potR_B'//trim(adjustl(ch)))
  OPEN(unit=14,file='valeurs/densR_B'//trim(adjustl(ch)))
  DO j = 1, nd
     READ(11,*) V01(j)
     READ(12,*) d01(j)
     READ(13,*) VR01(j)
     READ(14,*) dR01(j)
  END DO
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)

  !Compute the coefficients related to the transparent boundary conditions
  CALL SN(R,0._wp,sg,nMax)
  CALL SN(R,-dV,sd,nMax)
  
  WRITE(111,*)'** Instationnaire'
  WRITE(111,*)
  
  ! Number of time steps
  nf = FLOOR(dt0/dt*REAL(nf0))

  ! Solve the time dependent Schrodinger Poisson model by the projection method and write the distance to the time-independent potential corresponding to the new potential bias
  CALL CPU_TIME(debut)
  CALL SOL_INSTAT_RES(dR0,dR01,VR0,VR01,psiEnm1,thetanm1,enm1,znm1,nf)
  CALL CPU_TIME(fin)
  WRITE(111,*)'* Fin calcul resonance iteration', nf
  WRITE(111,*)'Temps calcul', fin-debut
  WRITE(111,*)'Ecart potentiel', 100._wp*NORM2(VR01-VR0)/NORM2(VR01)
  WRITE(111,*)'Ecart densite', 100._wp*NORM2(dR01(n1:n2)-dR0(n1:n2))/NORM2(dR01(n1:n2))
  WRITE(111,*)

  ! Write the Schrodinger Poisson final time potential and density in output files 
  OPEN(unit=11,file='valeurs/potRf')
  OPEN(unit=12,file='valeurs/densRf')
  DO j = 1, nd
     WRITE(11,*) tx(j), -VR0(j)
     WRITE(12,*) tx(j), dR0(j)
  END DO
  CLOSE(11)
  CLOSE(12)

  ! Number of time steps
  nf = FLOOR(dt0/dt*REAL(nf0))  

  ! Solve the time dependent Schrodinger Poisson model with a refined mesh (reference method) - instead of the projection method - and write the distance to the time-independent potential corresponding to the new potential bias
  CALL CPU_TIME(debut)
  CALL SOL_INSTAT(d0,d01,V0,V01,psinm1,enm1,znm1,nf)
  CALL CPU_TIME(fin)
  WRITE(111,*)'* Fin calcul direct iteration', nf
  WRITE(111,*)'Temps calcul', fin-debut
  WRITE(111,*)'Ecart potentiel', 100._wp*NORM2(V01-V0)/NORM2(V01)
  WRITE(111,*)'Ecart densite', 100._wp*NORM2(d01(n1:n2)-d0(n1:n2))/NORM2(d01(n1:n2))

  ! Write the Schrodinger Poisson final time potential and density in output files
  OPEN(unit=11,file='valeurs/potf')
  OPEN(unit=12,file='valeurs/densf')
  DO j = 1, nd
     WRITE(11,*) tx(j), -V0(j)
     WRITE(12,*) tx(j), d0(j)
  END DO
  CLOSE(11)
  CLOSE(12)

  CLOSE(111)

END PROGRAM resinst
