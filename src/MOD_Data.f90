!*******************************************************************************
! Module containing subroutines used by OrbData
!
! Copyright (C) 2022 Jon May, Peter Bird, Michele Carafa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!*******************************************************************************

module DATA_subs

use SharedVars
use ShellSetSubs

IMPLICIT NONE

CONTAINS ! (only to avoid need for INTERFACEs; no use of global variables is intended)

       REAL*8 FUNCTION ATan2F (y, x)

!   Corrects for problem of two zero arguments

       IMPLICIT NONE
       REAL*8, INTENT(IN) :: y, x ! Note: Order is intentional: (motivator, restrainer).
       IF ((y /= 0.0D0).OR.(x /= 0.0D0)) THEN
            ATan2F = ATAN2(y, x)
       ELSE
            ATan2F = 0.0D0
       END IF
       END FUNCTION ATan2F


       SUBROUTINE Squeez (alphaT, density_anomaly_kgpm3, elevat, &  ! INTENT(IN)
     &                    geoth1, geoth2, geoth3, geoth4, &         ! INTENT(IN)
     &                    geoth5, geoth6, geoth7, geoth8, &         ! INTENT(IN)
     &                     gMean, &                                 ! INTENT(IN)
     &                    iUnitT,  oneKm, rhoAst, rhoBar, rhoH2O, & ! INTENT(IN)
     &                    temLim,     zM,  zStop, &                 ! INTENT(IN)
     &                     tauZZ, sigZZB)                           ! INTENT(OUT)

!   Calculates "tauZZ", the vertical integral through the plate
!      of the vertical standardized stress anomaly, which is
!      relative to a column of mantle at asthenosphere temperature
!      with a 5-km crust and a 2.7-km ocean on top, like a mid-ocean
!      rise.  The integral is from either the land surface or the
!      sea surface, down to a depth of "zStop" below the top of
!      the crust.
!      If "zStop" exceeds Moho depth "zM", then properties of the mantle
!      will be used in the lower part of the integral.
!   Also returns "sigZZB", the standardized vertical stress anomaly
!      at depth "zStop" below the solid rock surface.
!   Note: This version is different from the version found in the -Laramy-
!      program package.  First, it acts on only a single point.
!      Second, it infers sub-plate normal-stress anomalies from
!      the given topography, instead of vice-versa.

       IMPLICIT NONE
       REAL*8, INTENT(IN) :: alphaT, density_anomaly_kgpm3, elevat, &
                           & geoth1, geoth2, geoth3, geoth4,        &
                           & geoth5, geoth6, geoth7, geoth8,        &
                           &  gMean
       INTEGER, INTENT(IN) :: iUnitT
       REAL*8, INTENT(IN) :: oneKm, rhoAst, rhoBar, rhoH2O,         &
                           & temLim,    zM,  zStop
       REAL*8, INTENT(OUT) :: tauZZ, sigZZB
!   Argument arrays:
       DIMENSION alphaT(2), rhoBar(2), temLim(2)

       INTEGER, PARAMETER :: nDRef = 300
       INTEGER :: i, j, lastDR, layer1, layer2, n1, n2, nStep
       LOGICAL :: called =.FALSE.
       REAL*8 :: dense, dense1, dense2, frac, frac1, frac2, h, &
               & oldPr, oldSZZ, Pr, resid, rhoTop, sigZZ, T, z, zBase, zTop
       REAL*8 :: tempC, tempM ! <-- statement functions, to be defined below...
       REAL*8, DIMENSION(0:nDRef) :: dRef, PRef

!   Statement functions:
       tempC(h) = MIN(temLim(1), geoth1 + geoth2 * h + geoth3 * h**2 &
     &                                       + geoth4 * h**3)
       tempM(h) = MIN(temLim(2), geoth5 + geoth6 * h + geoth7 * h**2 &
     &                                       + geoth8 * h**3)

!   Create reference temperature & density profiles to depth of nDRef km:

       IF (.NOT.called) THEN
            rhoTop = rhoBar(1) * (1.0D0 - alphaT(1) * geoth1)
            dRef(1) = rhoH2O
            dRef(2) = rhoH2O
            dRef(3) = 0.7D0 * rhoH2O + 0.3D0 * rhoTop
            dRef(4) = rhoTop
            dRef(5) = rhoTop
            dRef(6) = rhoTop
            dRef(7) = rhoTop
            dRef(8) = 0.7D0 * rhoTop + 0.3D0 * rhoAst
            DO 50 j = 9, nDRef
                 dRef(j) = rhoAst
   50       CONTINUE
            PRef(0) = 0.0D0
            DO 100 i = 1, nDRef
                 PRef(i) = PRef(i - 1) + dRef(i) * gMean * oneKm
  100       CONTINUE
       END IF

!   Routine processing (in every CALL):

       IF (elevat > 0.0D0) THEN
!        Land:
            zTop = -elevat
            zBase = zStop - elevat
            dense1 = rhoBar(1) * (1.0D0 - geoth1 * alphaT(1)) + density_anomaly_kgpm3
            h = 0.0D0
            layer1 = 1
       ELSE
!         OCEAN
            zTop = 0.0D0
            zBase = zStop + (-elevat)
            dense1 = rhoH2O
            h = elevat
            layer1 = 0
       END IF
       lastDR = zBase / oneKm
       IF (zBase > (oneKm * lastDR)) lastDR = lastDR + 1
       IF (lastDR > nDRef) THEN
	     write(ErrorMsg,'(A,I10)') "IN SUBPROGRAM SQUEEZ, PARAMETER nDRef MUST BE INCREASED TO AT LEAST ",lastDR
		 call FatalError(ErrorMsg,ThID)
       END IF
       nStep = (zBase - zTop) / oneKm
       oldSZZ = 0.0D0
       oldPr = 0.0D0
       sigZZ = 0.0D0
       tauZZ = 0.0D0
       z = zTop
       DO 200 i = 1, nStep
            z = z + oneKm
            h = h + oneKm
            IF (h > 0.0D0) THEN
                 IF (h <= zM) THEN
                      T = tempC(h)
                      dense2 = rhoBar(1) * (1. - T * alphaT(1)) + density_anomaly_kgpm3
                      layer2 = 1
                 ELSE
                      T = tempM(h - zM)
                      dense2 = rhoBar(2) * (1. - T * alphaT(2)) + density_anomaly_kgpm3
                      layer2 = 2
                 END IF
            ELSE
                 dense2 = rhoH2O
                 layer2 = 0
            END IF
            IF ((layer1 == 0).AND.(layer2 == 1)) THEN
                 frac2 = h / oneKm
                 frac1 = 1.0D0 - frac2
            ELSE IF ((layer1 == 1).AND.(layer2 == 2)) THEN
                 frac2 = (h - zM) / oneKm
                 frac1 = 1.0D0 - frac2
            ELSE
                 frac1 = 0.5D0
                 frac2 = 0.5D0
            END IF
            dense = frac1 * dense1 + frac2 * dense2
            IF (z > 0.0D0) THEN
                 n1 = z / oneKm
                 n2 = n1 + 1
                 frac = (z / oneKm) - n1
                 Pr = PRef(n1) + frac * (PRef(n2) - PRef(n1))
            ELSE
                 Pr = 0.0D0
            END IF
            sigZZ = sigZZ - dense * gMean * oneKm + (pr - oldPr)
            tauZZ = tauZZ + 0.5D0 * (sigZZ + oldSZZ) * oneKm
            dense1 = dense2
            oldSZZ = sigZZ
            oldPr = Pr
            layer1 = layer2
  200  CONTINUE
       resid = zBase - z
       h = zStop
       z = zBase
       IF (zStop <= zM) THEN
            T = tempC(h)
            dense2 = rhoBar(1) * (1.0D0 - T * alphaT(1)) + density_anomaly_kgpm3
       ELSE
            T = tempM(h - zM)
            dense2 = rhoBar(2) * (1.0D0 - T * alphaT(2)) + density_anomaly_kgpm3
       END IF
       dense = 0.5D0 * (dense1 + dense2)
       IF (z > 0.0D0) THEN
            n1 = z / oneKm
            n2 = n1 + 1
            frac = (z / oneKm) - n1
            Pr = PRef(n1) + frac * (PRef(n2) - PRef(n1))
       ELSE
            Pr = 0.0D0
       END IF
       sigZZB = sigZZ - dense * gMean * resid + (Pr - oldPr)
       tauZZ = tauZZ + 0.5D0 * (sigZZB + oldSZZ) * resid
       END SUBROUTINE Squeez

       SUBROUTINE Assign (aArray,    aX1,    aDX,    aX2,    nAX,    aDY,    aY2,    nAY, & ! INTENT(IN)
     &                    alphaT, cLimit, conduc, &                                         ! INTENT(IN)
     &                    cArray,    cX1,    cDX,    cX2,    nCX,    cDY,    cY2,    nCY, & ! INTENT(IN)
     &                    delta_rho_limit, &                                                ! INTENT(IN)
     &                    eArray,    eX1,    eDX,    eX2,    nEX,    eDY,    eY2,    nEY, & ! INTENT(IN)
     &                     gMean,  hCMax,  hLMax, &                                         ! INTENT(IN)
     &                    iUnitL, iUnitT, &                                                 ! INTENT(IN)
     &                     oneKm, &                                                         ! INTENT(IN)
     &                      pLon,   pLat, &                                                 ! INTENT(IN)
     &                     qLim0, dQL_dE,  qLim1, &                                         ! INTENT(IN)
     &                    qArray,    qX1,    qDX,    qX2,    nQX,    qDY,    qY2,    nQY, & ! INTENT(IN)
     &                     radio, rhoAst, rhoBar, rhoH2O, &                                 ! INTENT(IN)
     &                    sArray,    sX1,    sDX,    sX2,    nSX,    sDY,    sY2,    nSY, & ! INTENT(IN)
     &                    TAsthK, temLim,  TSurf, &                                         ! INTENT(IN)
     &                    elevat, heatFl, &                                                 ! INTENT(INOUT)
     &                    thickC, thickM, chemical_delta_rho, cooling_curvature, &          ! INTENT(OUT)
     &                    ModNum)

!   Determines elevation (if 0.0), heat-flow (if 0.0), and
!   values of crustal and mantle-lithosphere thickness, at one point:
!   (pLon = East longitude, pLat = North latitude), both in degrees.
!  (Unlike subprogram -Assign- in ALGriddr, this version never
!   does regional averaging.  Neither does it interpolate by
!   Kriging; simple interpolation within rectangular grids is
!   used instead.)

!   This version of Assign was completely re-written for OrbData5.
!   Instead of solving for crustal thickness from isostasy, it
!   reads crustal thickness from a grid of values.
!   Instead of computing lithosphere thickness from heat flow,
!   it computes it from seafloor age or from vertical S-wave
!   travel-time anomaly under continents.  Then extra degrees of
!   freedom "chemical_Delta_rho" and "cooling_curvature" are
!   used to create a self-consistent structure.

!   Note that the result may *FAIL* to be isostatic with mid-ocean
!   ridges if that would require chemical_Delta_rho to go beyond
!   the prescribed amplitude.  Also, if the requested lithosphere
!   thickness is more than ~twice that expected based on heat-flow
!  (and a steady-state assumption), then the lithosphere thickness
!   returned will be limited so as to prevent geotherm overshoot
!   and negative temperature gradients.

!   Algorithm revision by Peter Bird, UCLA, June 2005;
!   syntax revision to Fortran 90 by Peter Bird, UCLA, December 2018.
       IMPLICIT NONE
       REAL*8, INTENT(IN) :: aArray,    aX1,    aDX,    aX2,    aDY,    aY2, &
                           & alphaT, cLimit, conduc, &
                           & cArray,    cX1,    cDX,    cX2,    cDY,    cY2, &
                           & delta_rho_limit, &
                           & eArray,    eX1,    eDX,    eX2,    eDY,    eY2, &
                           &  gMean,  hCMax,  hLMax, &
                           &  oneKm, &
                           &   pLon,   pLat, &
                           &  qLim0, dQL_dE,  qLim1, &
                           & qArray,    qX1,    qDX,    qX2,    qDY,    qY2, &
                           &  radio, rhoAst, rhoBar, rhoH2O, &
                           & sArray,    sX1,    sDX,    sX2,    sDY,    sY2, &
                           & TAsthK, temLim,  TSurf
       INTEGER, INTENT(IN) ::    nAX,    nAY,    nCX,    nCY,    nEX,    nEY, &
                            & iUnitL, iUnitT, &
                            &    nQX,    nQY,    nSX,    nSY
        REAL*8, INTENT(INOUT) :: elevat, heatFl
        REAL*8, INTENT(OUT) :: thickC, thickM, chemical_delta_rho, cooling_curvature
        !Argument arrays ALLOCATED and dimensioned in calling program:
        DIMENSION aArray(:, :), cArray(:, :), eArray(:, :), qArray(:, :), sArray(:, :)
        !Argument arrays with (crust:mantle) values:
        DIMENSION alphaT(2), conduc(2), radio(2), rhoBar(2), temLim(2)
!---------------------------------------------------------------------
        !Internal variables:
        INTEGER :: ic1, ic2, ir1, ir2
        LOGICAL :: badP, badT, outsid, needE, needQ, &
      &            warnC1, warnC2, warnM1, warnL2, wayOut
        REAL*8  :: ageMa, bot, c0_of_mantle_gradient, c1_of_mantle_gradient, &
                 & deltaQ, delta_quadratic, delta_T, delta_tS, fc, fr, &
                 & geoth1, geoth2, geoth3, geoth4, geoth5, geoth6, geoth7, geoth8, &
                 & h_Earth3, h_Earth5, h_plate, m_star, &
                 & q_gdh1, q_radioactivity, qLimit, qRed, ridge_delta_tS, sigZZB, slope, &
                 & t_geoth1, t_geoth2, t_geoth3, t_geoth5, t_geoth6, t_geoth7, t_qRed, t_TMoho, &
                 & tauZZ, TErr0r, test, TMoho, top, total_lithosphere, useLon, &
                 & z_of_maximum, z_star
        character(len=100) :: msg
        integer,intent(in) :: ModNum


	   outsid = .FALSE.
       wayOut = .FALSE.
       warnC1 = .FALSE.
       warnC2 = .FALSE.
       warnM1 = .FALSE.
       warnL2 = .FALSE.
!   Determine the elevation, if needed:
       needE = (elevat == 0.0D0)
       IF (needE) THEN
            ir1 = ((eY2 - pLat) / eDY) + 1.00001D0 ! truncated to INTEGER
            ir1 = MAX(ir1, 1)
            ir1 = MIN(ir1, nEY - 1)
            ir2 = ir1 + 1
            fr = ((eY2 - eDY * (ir1 - 1)) - pLat) / eDY

            IF ((pLon < eX1).OR.(pLon > eX2)) THEN
                 IF (((pLon + 360.0D0) >= eX1).AND. &
     &               ((pLon + 360.0D0) <= eX2)) THEN
                      useLon = pLon + 360.0D0
                 ELSE IF (((pLon - 360.0D0) >= eX1).AND. &
     &                    ((pLon - 360.0D0) <= eX2)) THEN
                      useLon = pLon - 360.0D0
                 ELSE
                      useLon = pLon
                 END IF
            ELSE
                 useLon = pLon
            END IF
            ic1 = ((useLon - eX1) / eDX) + 1.00001D0 ! truncated to INTEGER
            ic1 = MAX(ic1, 1)
            ic1 = MIN(ic1, nEX - 1)
            ic2 = ic1 + 1
            fc = (useLon - (eX1 + eDX * (ic1 - 1))) / eDX

            outsid = outsid.OR.(fr < -0.01D0).OR.(fr > 1.01D0).OR. &
     &                         (fc < -0.01D0).OR.(fc > 1.01D0)
            wayOut = wayOut.OR.(fr < -1.01D0).OR.(fr > 2.01D0).OR. &
     &                         (fc < -1.01D0).OR.(fc > 2.01D0)
            fr = MIN(1.0D0, MAX(0.0D0, fr))
            fc = MIN(1.0D0, MAX(0.0D0, fc))
            top = eArray(ir1, ic1) + fc * (eArray(ir1, ic2) - eArray(ir1, ic1))
            bot = eArray(ir2, ic1) + fc * (eArray(ir2, ic2) - eArray(ir2, ic1))
            elevat = top + fr * (bot - top)
       END IF

!   Determine age of sea-floor, which is always needed:

       ir1 = ((aY2 - pLat) / aDY) + 1.00001D0 ! truncated to INTEGER
       ir1 = MAX(ir1, 1)
       ir1 = MIN(ir1, nAY - 1)
       ir2 = ir1 + 1
       fr = ((aY2 - aDY * (ir1 - 1)) - pLat) / aDY

       IF ((pLon < aX1).OR.(pLon > aX2)) THEN
            IF (((pLon + 360.0D0) >= aX1).AND. &
     &          ((pLon + 360.0D0) <= aX2)) THEN
                 useLon = pLon + 360.0D0
            ELSE IF (((pLon - 360.0D0) >= aX1).AND. &
     &               ((pLon - 360.0D0) <= aX2)) THEN
                 useLon = pLon - 360.0D0
            ELSE
                 useLon = pLon
            END IF
       ELSE
            useLon = pLon
       END IF
       ic1 = ((useLon - aX1) / aDX) + 1.00001D0 ! truncated to INTEGER
       ic1 = MAX(ic1, 1)
       ic1 = MIN(ic1, nAX - 1)
       ic2 = ic1 + 1
       fc = (useLon - (aX1 + aDX * (ic1 - 1))) / aDX

       outsid = outsid.OR.(fr < -0.01D0).OR.(fr > 1.01D0).OR. &
     &                    (fc < -0.01D0).OR.(fc > 1.01D0)
       wayOut = wayOut.OR.(fr < -1.01D0).OR.(fr > 2.01D0).OR. &
     &                    (fc < -1.01D0).OR.(fc > 2.01D0)
       fr = MIN(1.0D0, MAX(0.0D0, fr))
       fc = MIN(1.0D0, MAX(0.0D0, fc))
       top = aArray(ir1, ic1) + fc * (aArray(ir1, ic2) - aArray(ir1, ic1))
       bot = aArray(ir2, ic1) + fc * (aArray(ir2, ic2) - aArray(ir2, ic1))
       ageMa = top + fr * (bot - top)

!   Determine heat-flow, in needed:

       needQ = (heatFl == 0.0D0)
       IF (needQ) THEN
            ir1 = ((qY2 - pLat) / qDY) + 1.00001D0 ! truncated to INTEGER
            ir1 = MAX(ir1, 1)
            ir1 = MIN(ir1, nQY - 1)
            ir2 = ir1 + 1
            fr = ((qY2 - qDY * (ir1 - 1)) - pLat) / qDY

            IF ((pLon < qX1).OR.(pLon > qX2)) THEN
                 IF (((pLon + 360.0D0) >= qX1).AND. &
     &               ((pLon + 360.0D0) <= qX2)) THEN
                      useLon = pLon + 360.0D0
                 ELSE IF (((pLon - 360.0D0) >= qX1).AND. &
     &                    ((pLon - 360.0D0) <= qX2)) THEN
                      useLon = pLon - 360.0D0
                 ELSE
                      useLon = pLon
                 END IF
            ELSE
                 useLon = pLon
            END IF
            ic1 = ((useLon - qX1) / qDX) + 1.00001
            ic1 = MAX(ic1, 1)
            ic1 = MIN(ic1, nQX - 1)
            ic2 = ic1 + 1
            fc = (useLon - (qX1 + qDX * (ic1 - 1))) / qDX

            outsid = outsid.OR.(fr < -0.01D0).OR.(fr > 1.01D0).OR. &
     &                         (fc < -0.01D0).OR.(fc > 1.01D0)
            wayOut = wayOut.OR.(fr < -1.01D0).OR.(fr > 2.01D0).OR. &
     &                         (fc < -1.01D0).OR.(fc > 2.01D0)
            fr = MIN(1.0D0, MAX(0.0D0, fr))
            fc = MIN(1.0D0, MAX(0.0D0, fc))
            top = qArray(ir1, ic1) + fc * (qArray(ir1, ic2) - qArray(ir1, ic1))
            bot = qArray(ir2, ic1) + fc * (qArray(ir2, ic2) - qArray(ir2, ic1))
            heatFl = top + fr * (bot - top)

!           Check for seafloor-age overriding heat-flow grid value:

            IF (ageMa < 200.0D0) THEN
!              Seafloor age is valid (not "unknown" or "continental")
                 IF (ageMa <= 0.0D0) THEN
                      heatFl = qLim1
                 ELSE
!                     Carol A. Stein & Seth Stein [1992]
!                     A model for the global variation in oceanic
!                     depth and heat flow with lithospheric age,
!                     Nature, v. 359, 10 September, p. 123-129.
!                     According to their preferred GDH1 model:
                      IF (ageMa <= 55.0D0) THEN
                           heatFl = 0.510D0 / SQRT(ageMa)
                      ELSE
                           heatFl = 0.048D0 + 0.096D0 * EXP(-0.0278D0 * ageMa)
                      END IF
                 END IF
            END IF

       END IF

!  Apply limits on Q:

       qLimit = qLim0 + dQL_dE * elevat
       heatFl = MAX(heatFl, qLimit)
       heatFl = MIN(heatFl, qLim1)

!  Obtain crustal thickness from grid cArray:

       ir1 = ((cY2 - pLat) / cDY) + 1.00001D0 ! truncated to INTEGER
       ir1 = MAX(ir1, 1)
       ir1 = MIN(ir1, nCY - 1)
       ir2 = ir1 + 1
       fr = ((cY2 - cDY * (ir1 - 1)) - pLat) / cDY
       IF ((pLon < cX1).OR.(pLon > cX2)) THEN
            IF (((pLon + 360.0D0) >= cX1).AND. &
     &          ((pLon + 360.0D0) <= cX2)) THEN
                 useLon = pLon + 360.0D0
            ELSE IF (((pLon - 360.0D0) >= cX1).AND. &
     &               ((pLon - 360.0D0) <= cX2)) THEN
                 useLon = pLon - 360.0D0
            ELSE
                 useLon = pLon
            END IF
       ELSE
            useLon = pLon
       END IF
       ic1 = ((useLon - cX1) / cDX) + 1.00001D0 ! truncated to INTEGER
       ic1 = MAX(ic1, 1)
       ic1 = MIN(ic1, nCX - 1)
       ic2 = ic1 + 1
       fc = (useLon - (cX1 + cDX * (ic1 - 1))) / cDX
       outsid = outsid.OR.(fr < -0.01D0).OR.(fr > 1.01D0).OR. &
     &                    (fc < -0.01D0).OR.(fc > 1.01D0)
       wayOut = wayOut.OR.(fr < -1.01D0).OR.(fr > 2.01D0).OR. &
     &                    (fc < -1.01D0).OR.(fc > 2.01D0)
       fr = MIN(1.0D0, MAX(0.0D0, fr))
       fc = MIN(1.0D0, MAX(0.0D0, fc))
       top = cArray(ir1, ic1) + fc * (cArray(ir1, ic2) - cArray(ir1, ic1))
       bot = cArray(ir2, ic1) + fc * (cArray(ir2, ic2) - cArray(ir2, ic1))
       thickC = top + fr * (bot - top)
       thickC = MAX(thickC, cLimit)
       thickC = MIN(thickC, hCMax)

!     Give warning if any necessary data could not be found:

       IF (wayOut) THEN
		  write(ErrorMsg,'(A,F8.3,A,F7.3,A)') "WARNING: (",pLon,"E,",pLat,"N) IS FAR OUTSIDE DATA GRID(S). "!Necessary data will be extrapolated outside grid edge(s)."
		  call ModelError(ErrorMsg,ThID)
       ELSE IF (outsid) THEN
		  write(ErrorMsg,'(A,F8.3,A,F7.3,A)') "WARNING: (",pLon,"E,",pLat,"N) IS SLIGHTLY OUTSIDE DATA GRID(S)."
		  call ModelError(ErrorMsg,ThID)
       END IF

!  Method of determining lithosphere thickness depends on whether
!     node has known (and valid) sea-floor age:

       IF (ageMa < 200.0D0) THEN

!           This node is in oceanic lithosphere.

!        (1)Determine age-dependent model Q:

!           Carol A. Stein & Seth Stein [1992]
!           A model for the global variation in oceanic
!           depth and heat flow with lithospheric age,
!           Nature, v. 359, 10 September, p. 123-129.
!           According to their preferred GDH1 model:
            IF (ageMa <= 55.0D0) THEN
                 q_gdh1 = 0.510D0 / SQRT(ageMa)
            ELSE
                 q_gdh1 = 0.048D0 + 0.096D0 * EXP(-0.0278D0 * ageMa)
            END IF

!        (2)Find total lithosphere thickness expected according to
!           a steady-state conduction model (like old OrbData; see
!           spreadsheet S20RTS_delta_ts_vs_age.xls for calibration).

            q_gdh1 = MAX(q_gdh1, qLimit)
            q_gdh1 = MIN(q_gdh1, qLim1)
            q_radioactivity = 0.007D0
            delta_T = TAsthK - TSurf
            h_Earth3 = (delta_T * conduc(2)) / (q_gdh1 - q_radioactivity)

!        (3)Take geometric mean of this "Earth3" or "old OrbData"
!           thickness with the constant plate thickness of Stein &
!           Stein [1992]:

            h_plate = 95000.0D0
            h_Earth5 = SQRT(h_Earth3 * h_plate)

!        (4)Mantle lithosphere excludes crust:

            thickM = h_Earth5 - thickC

!        (5)Apply lower and upper limits:

            thickM = MAX(thickM, 0.0D0)
            thickM = MIN(thickM, (hLMax - thickC))

       ELSE

!           This node is in continental (or unknown) lithosphere.

!        (1)Read vertical S-wave travel-time anomaly (in the upper
!           mantle, from Moho to 400 km depth) from array sArray:

            ir1 = ((sY2 - pLat) / sDY) + 1.00001D0 ! truncated to INTEGER
            ir1 = MAX(ir1, 1)
            ir1 = MIN(ir1, nSY - 1)
            ir2 = ir1 + 1
            fr = ((sY2 - sDY * (ir1 - 1)) - pLat) / sDY
            IF ((pLon < sX1).OR.(pLon > sX2)) THEN
                 IF (((pLon + 360.0D0) >= sX1).AND. &
     &               ((pLon + 360.0D0) <= sX2)) THEN
                      useLon = pLon + 360.0D0
                 ELSE IF (((pLon - 360.0D0) >= sX1).AND. &
     &                    ((pLon - 360.0D0) <= sX2)) THEN
                      useLon = pLon - 360.0D0
                 ELSE
                      useLon = pLon
                 END IF
            ELSE
                 useLon = pLon
            END IF
            ic1 = ((useLon - sX1) / sDX) + 1.00001D0 ! truncated to INTEGER
            ic1 = MAX(ic1, 1)
            ic1 = MIN(ic1, nSX - 1)
            ic2 = ic1 + 1
            fc = (useLon - (sX1 + sDX * (ic1 - 1))) / sDX
            outsid = outsid.OR.(fr < -0.01D0).OR.(fr > 1.01D0).OR. &
     &                         (fc < -0.01D0).OR.(fc > 1.01D0)
            wayOut = wayOut.OR.(fr < -1.01D0).OR.(fr > 2.01D0).OR. &
     &                         (fc < -1.01D0).OR.(fc > 2.01D0)
            fr = MIN(1.0D0, MAX(0.0D0, fr))
            fc = MIN(1.0D0, MAX(0.0D0, fc))
            top = sArray(ir1, ic1) + fc * (sArray(ir1, ic2) - sArray(ir1, ic1))
            bot = sArray(ir2, ic1) + fc * (sArray(ir2, ic2) - sArray(ir2, ic1))
            delta_tS = top + fr * (bot - top)

!        (2)Scale to mantle lithosphere thickness, based on
!           calibration performed in oceanic lithosphere of known age
!          (see diary.doc entry for 2005.06.02):

            slope = (1.494D0) / (85000.0D0)
!          (in units of seconds per meter of lithosphere)
            ridge_delta_tS = 1.636D0
!          (seconds; see diary.doc entry for 2005.06.02)

            thickM = (ridge_delta_ts - delta_ts) / slope

!        (3)Apply lower and upper limits:

            thickM = MAX(thickM, 0.0D0)
            thickM = MIN(thickM, (hLMax - thickC))

       END IF
!     (End of branch on oceanic vs. continental/unknown lithosphere)

!   Compute steady-state GEOTHerm, as in old OrbData:

       geoth1 = TSurf
       geoth2 = heatFl / conduc(1)
       geoth3 = -radio(1) / (2.0D0 * conduc(1))
       geoth4 = 0.0D0
       TMoho = geoth1 + geoth2 * thickC + geoth3 * thickC**2
       geoth5 = TMoho
       qRed = heatFl - thickC * radio(1)
       geoth6 = qRed / conduc(2)
       geoth7 = -radio(2) / (2.0D0 * conduc(2))
       geoth8 = 0.0D0

!   Check for Moho hotter than asthenosphere, and
!      reduce heat flow if necessary.
!      N.B. This typically occurs where spreading ridges
!      with rule-based high Q intersect continental
!      margins, where CRUST2 shows thick crust!

       test = geoth1 + geoth2 * thickC + geoth3 * thickC**2
       IF (test > TAsthK) THEN
            TErr0r = test - TAsthK
            deltaQ = -TErr0r * conduc(1) / thickC
            heatFl = heatFl + deltaQ
            qLimit = qLim0 + dQL_dE * elevat
            heatFl = MAX(heatFl, qLimit)
            heatFl = MIN(heatFl, qLim1)
            geoth2 = heatFl / conduc(1)
            TMoho = TAsthK
            geoth5 = TMoho
            qRed = heatFl - thickC * radio(1)
            geoth6 = qRed / conduc(2)
       END IF

!   Compute trial value of cooling_curvature
!     (referring to GEOTHerm without this effect):

       test = geoth5 + geoth6 * thickM + geoth7 * thickM**2
       TErr0r = test - TAsthK
       total_lithosphere = thickC + thickM
       delta_quadratic = -TErr0r / total_lithosphere**2
       cooling_curvature = -2.0D0 * delta_quadratic

!     (Note that we are postponing the addition of
!      delta_quadratic to geoth3 and geoth7 until this
!      value passes some tests!)

!   Test whether this value would cause any temperature
!      maximum within the lithosphere (forbidden!)?:

       IF (cooling_curvature > 0.0D0) THEN
            t_geoth1 = TSurf
            t_geoth2 = heatFl / conduc(1)
            t_geoth3 = delta_quadratic - radio(1) / (2.0D0 * conduc(1))
            t_TMoho = t_geoth1 + t_geoth2 * thickC + t_geoth3 * thickC**2
            t_geoth5 = t_TMoho
            t_qRed = heatFl - thickC * radio(1) - cooling_curvature * thickC * conduc(2)
            t_geoth6 = qRed / conduc(2)
            t_geoth7 = delta_quadratic - radio(2) / (2.0D0 * conduc(2))

!           Check for temperature maximum within mantle lithosphere:
            c0_of_mantle_gradient = t_geoth6
            c1_of_mantle_gradient = 2. * t_geoth7
            z_of_maximum = c0_of_mantle_gradient / (-c1_of_mantle_gradient)
            IF ((z_of_maximum < thickM) .AND. (z_of_maximum > 0.0D0)) THEN

!                Must take corrective action; reduce thickM
!                   and also change cooling_curvature,
!                   delta_quadratic, etc..  For algebra, see Peter
!                   Bird's notes of 2005.06.15.

                 z_star = (geoth7 * thickC**2 + geoth5 - &
     &                     geoth6 * thickC - TAsthK) / &
     &                    (geoth7 * thickC - 0.5D0 * geoth6)
                 m_star = z_star - thickC
                 thickM = MAX(MIN(m_star, thickM), 0.0D0)
                 total_lithosphere = thickC + thickM
                 TErr0r = (geoth5 + geoth6 * thickM + geoth7 * thickM**2) - TAsthK
                 delta_quadratic = -TErr0r / total_lithosphere**2
                 cooling_curvature = -2.0D0 * delta_quadratic
            END IF
       END IF

!   Build geotherm again with final cooling_curvature:

       geoth1 = TSurf
       geoth2 = heatFl / conduc(1)
       geoth3 = delta_quadratic - radio(1) / (2.0D0 * conduc(1))
       geoth4 = 0.0D0
       TMoho = geoth1 + geoth2 * thickC + geoth3 * thickC**2
       geoth5 = TMoho
       qRed = heatFl - thickC * radio(1) - cooling_curvature * thickC * conduc(2)
       geoth6 = qRed / conduc(2)
       geoth7 = delta_quadratic - radio(2) / (2.0D0 * conduc(2))
       geoth8 = 0.0D0

       chemical_delta_rho = 0.0D0
!     (Try this case first; adjust density_anomaly below.)

       CALL Squeez (alphaT, chemical_delta_rho, elevat, & ! INTENT(IN)
     &              geoth1, geoth2, geoth3, geoth4, &     ! INTENT(IN)
     &              geoth5, geoth6, geoth7, geoth8, &     ! INTENT(IN)
     &               gMean, iUnitT, &                     ! INTENT(IN)
     &               oneKm, rhoAst, rhoBar, rhoH2O, &     ! INTENT(IN)
     &              temLim, thickC, thickC + thickM, &    ! INTENT(IN)
     &               tauZZ, sigZZB)                       ! INTENT(OUT)

!   Trial value of chemical_density_anomaly:

       chemical_delta_rho = sigZZB / (gMean * total_lithosphere)

!   Apply limits:

       chemical_delta_rho = MIN(chemical_delta_rho,  delta_rho_limit)
       chemical_delta_rho = MAX(chemical_delta_rho, -delta_rho_limit)

!   Repeat isostasy test:

       CALL Squeez (alphaT, chemical_delta_rho, elevat, & ! INTENT(IN)
     &              geoth1, geoth2, geoth3, geoth4, &     ! INTENT(IN)
     &              geoth5, geoth6, geoth7, geoth8, &     ! INTENT(IN)
     &              gMean, iUnitT, &                      ! INTENT(IN)
     &              oneKm, rhoAst, rhoBar, rhoH2O, &      ! INTENT(IN)
     &              temLim, thickC, thickC + thickM, &    ! INTENT(IN)
     &              tauZZ, sigZZB)                        ! INTENT(OUT)

!   Test for successful calculation:

!   Isostasy is bad if basal pressure anomaly is large:
       badP = (ABS(sigZZB) > 40.0D6) ! in Pascals; so, 40 MPa.

!   Geotherm might not connect to asthenosphere adiabat:
       test = geoth5 + geoth6 * thickM + geoth7 * thickM**2
       TErr0r = test - TAsthK
       badT = (ABS(TErr0r) > 20.0D0) ! 20 degrees Kelvin

       IF (badP) THEN
            WRITE (iUnitL, 102) pLon, pLat, ABS(sigZZB)
  102       FORMAT(/'           ERROR: (', F8.3, 'E, ', F7.3, 'N) ABS(sigZZB):', F16.3,' IS NOT ISOSTATIC WITH RISE:')
	   END IF

       IF (badT) THEN
            WRITE (iUnitL, 103) pLon, pLat
  103       FORMAT(/'           ERROR: (', F8.3, 'E, ', F7.3, 'N) DOES NOT CONNECT TO ASTHENOSPHERE ADIABAT')
	   END IF
       END SUBROUTINE Assign

       SUBROUTINE Square (brief, fDip, iUnitT, &  ! input
     &                    log_strike_adjustments, &
     &                    mxBn, mxEl, mxFEl, mxNode, &
     &                    mxStar, nFl, nodeF, nodes, &
     &                    numEl, numNod, skipBC, radius, wedge, &
     &                    xNode, yNode, &         ! modify
     &                    area, detJ, &           ! output
     &                    dXS, dYS, dXSP, dYSP, edgeFS, &
     &                    edgeTS, fLen, fPFlt, fPSfer, &
     &                    fArg, nCond, nodCon, sita, &
     &                    checkN, list)           ! work

!  Check, correct, and complete the geometry of the finite element grid.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       LOGICAL, INTENT(IN) :: brief                                                           ! input
       REAL*8, INTENT(IN) :: fDip                                                             ! input
       INTEGER, INTENT(IN) :: iUnitT                                                        ! input
       LOGICAL, INTENT(IN) :: log_strike_adjustments                                          ! input
       INTEGER, INTENT(IN) :: mxBn, mxEl, mxFEl, mxNode, mxStar, nFl, nodeF, nodes, &         ! input
             & numEl, numNod                                                                  ! input
       LOGICAL, INTENT(IN) :: skipBC                                                          ! input
       REAL*8, INTENT(IN) :: radius, wedge                                                    ! input
       REAL*8, INTENT(INOUT) :: xNode, yNode                                                  ! modify
       REAL*8, INTENT(OUT) :: area, detJ, dXS, dYS, dXSP, dYSP                                ! output
       LOGICAL, INTENT(OUT) :: edgeFS, edgeTs                                                 ! output
       REAL*8, INTENT(OUT) :: fLen, fPFlt, fPSfer, fArg                                       ! output
       INTEGER, INTENT(OUT) :: nCond, nodCon                                                  ! output
       REAL*8, INTENT(OUT) :: sita                                                            ! output
       LOGICAL checkN                                                          ! work
       INTEGER list                                                            ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION fPhi, fPoint, fGauss
       COMMON / SFault / fPoint
       COMMON / FPhis /  fPhi
       COMMON / FGList / fGauss
       DIMENSION fPhi(4, 7), fPoint(7), fGauss(7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CHARACTER*21 obliqu, tag1, tag2, vertic
       INTEGER i, j, j1, j2, k, kEle, kFault, l, m, &
             & n, n1, n2, n4, nazi, nazl, nDone, nGood, nInSum, nj1, nj2, &
             & nl1, nl2, nl3, nl4, nLeft, node, node1, node4, np1, np4, number, nvpart
       REAL*8 azi, azimut, azl, cosz, &
            & daz, dazi, dazl, deld, dELon, deld1, deld2, dip1, dip2, dNLat, dx, dy, &
            & fAngle, phi, phi1, phi2, r, r2, rmax, short, sinz, &
            & t2, test, theLat, theLon, theta, theta1, theta2, toler, &
            & x, xmean, xsum, y, ymean, ysum
       LOGICAL agreed, allOK, found, switch, vert1, vert2
       DIMENSION fAngle(2), phi(2), theta(2)
       DIMENSION area(mxEl), checkN(mxNode), &
     &           detJ(7, mxEl), &
     &           dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
     &           dXSP(3, 7, mxEl), dYSP(3, 7, mxEl), &
     &           edgeFS(2, mxFEl), edgeTS(3, mxEl), fDip(2, mxFEl), &
     &           fLen(mxFEl), &
     &           fPFlt(2, 2, 2, 7, mxFEl), &
     &           fPSfer(2, 2, 3, 7, mxEl), fArg(2, mxFEl), &
     &           list(mxStar), nodCon(mxBn), &
     &           nodeF(4, mxFEl), nodes(3, mxEl), &
     &           sita(7, mxEl), xNode(mxNode), yNode(mxNode)
       DATA obliqu / '(DIP SLIP IS ALLOWED)' /
       DATA vertic / '(STRIKE-SLIP ONLY)   ' /

integer,dimension(:),allocatable :: Numbers

!  (1) Check that all nodes are connected to at least one
!      continuum (triangular) element or fault element:

       DO 110 i = 1, numNod
            checkN(i) = .FALSE.
  110  CONTINUE
       DO 130 i = 1, numEl
            DO 120 j = 1, 3
                 checkN(nodes(j, i)) = .TRUE.
  120       CONTINUE
  130  CONTINUE
       DO 136 i = 1, nFl
            DO 134 j = 1, 4
                 checkN(nodeF(j, i)) = .TRUE.
  134       CONTINUE
  136  CONTINUE
       allOK = .TRUE.
       DO 140 i = 1, numNod
            allOK = allOK.AND.checkN(i)
  140  CONTINUE
       IF (.NOT.allOK) THEN
		  write(ErrorMsg,'(A/,A)') "BAD GRID TOPOLOGY: FOLLOWING REAL NODES DO NOT ",&
		  &  "BELONG TO ANY TRIANGULAR CONTINUUM ELEMENT OR FAULT ELEMENT:"
		  allocate(Numbers(numNod))
          Numbers(1:numNod) = (/ (i, i=1,numNod) /)
		  allocate(ErrorArrayInt(count(.not.checkN)))
		  ErrorArrayInt = pack(Numbers, .not.checkN)
		  call FatalError(ErrorMsg,ThID,ErrArrInt=ErrorArrayInt)
       END IF

!  (2) Average together the coordinates of all nodes at one "point":

       DO 410 i = 1, numNod
            checkN(i) = .FALSE.
!           (Means "not yet involved in averaging")
  410  CONTINUE
       DO 490 i = 1, nFl
            DO 480 j1 = 1, 2
                 nj1 = nodeF(j1, i)
!               (Fault ends are the only places that can have problems.)
                 IF (.NOT.checkN(nj1)) THEN
                      list(1) = nj1
                      checkN(nj1) = .TRUE.
!                    Begin list of neighbors with paired node:
                      j2 = 5 - j1
                      nj2 = nodeF(j2, i)
                      list(2) = nj2
                      checkN(nj2) = .TRUE.
                      nInSum = 2
!                    Find shortest fault connected to either one:
                      dx = xNode(nj1) - xNode(nj2)
                      dy = yNode(nj1) - yNode(nj2)
                      IF (dy > 3.14159265358979D0) dy = dy - 6.28318530717959D0
                      IF (dy < -3.14159265358979D0) dy = dy + 6.28318530717959D0
                      dy = dy * SIN(xNode(nj1))
                      short = SQRT(dx**2 + dy**2)
                      DO 470 k = 1, nFl
                           nl1 = nodeF(1, k)
                           nl2 = nodeF(2, k)
                           nl3 = nodeF(3, k)
                           nl4 = nodeF(4, k)
                           IF ((nj1 == nl1).OR.(nj2 == nl1).OR. &
     &                         (nj1 == nl2).OR.(nj2 == nl2).OR. &
     &                         (nj1 == nl3).OR.(nj2 == nl3).OR. &
     &                         (nj1 == nl4).OR.(nj2 == nl4)) THEN
                                dx = xNode(nl1) - xNode(nl2)
                                dy = yNode(nl1) - yNode(nl2)
                                IF (dy > 3.14159265358979D0) dy = dy - 6.28318530717959D0
                                IF (dy < -3.14159265358979D0) dy = dy + 6.28318530717959D0
                                dy = dy * SIN(xNode(nl1))
                                test = SQRT(dx**2 + dy**2)
                                short = MIN(short, test)
                           END IF
  470                 CONTINUE
!                    Collect all corner nodes within 10% of this:
                      toler = short / 10.0D0
                      t2 = toler**2
                      DO 471 k = 1, numNod
                           IF (.NOT.checkN(k)) THEN
                                dx = xNode(nj1) - xNode(k)
                                dy = yNode(nj1) - yNode(k)
                                IF (dy > 3.14159265358979D0) dy = dy - 6.28318530717959D0
                                IF (dy < -3.14159265358979D0) dy = dy + 6.28318530717959D0
                                dy = dy * SIN(xNode(nj1))
                                r2 = dx**2 + dy**2
                                IF (r2 < t2) THEN
                                     nInSum = nInSum + 1
                                          IF (nInSum > mxStar) THEN
	                                         write(ErrorMsg,'(A)') "INCREASE VALUE OF PARAMETER mxStar."
		                                     call FatalError(ErrorMsg,ThID)
                                          END IF
                                     list(nInSum) = k
                                     checkN(k) = .TRUE.
                                END IF
                           END IF
  471                 CONTINUE
!                    (Quick EXIT if all nodes in same place)
                      agreed = .TRUE.
                      DO 472 k = 2, nInSum
                           agreed = agreed.AND. &
     &                           (xNode(list(k)) == xNode(list(1))).AND. &
     &                           (yNode(list(k)) == yNode(list(1)))
  472                 CONTINUE
                      IF (agreed) GO TO 480
                      xsum = 0.0D0
                      ysum = 0.0D0
                      DO 473 k = 1, nInSum
                           xsum = xsum + xNode(list(k))
                           ysum = ysum + yNode(list(k))
  473                 CONTINUE
                      xmean = xsum / nInSum
                      ymean = ysum / nInSum
                      rmax = 0.0D0
                      DO 474 k = 1, nInSum
                           r = SQRT((xNode(list(k)) - xmean)**2 + &
     &                              (yNode(list(k)) - ymean)**2)
                           rmax = MAX(rmax, r)
  474                 CONTINUE
                      DO 475 k = 1, nInSum
                           xNode(list(k)) = xmean
                           yNode(list(k)) = ymean
  475                 CONTINUE
                      IF (.NOT.brief) THEN
                           IF (rmax > 0.0D0) THEN
                                IF(Verbose) WRITE(iUnitT,476) nInSum, &
     &                               (list(n), n = 1, nInSum)
  476                           FORMAT(/ &
     &                           ' AVERAGING TOGETHER THE POSITIONS OF', &
     &                               ' THESE ',I6,' NODES:',(/' ',12I6))
                                IF(Verbose) WRITE(iUnitT,477) rmax
  477                           FORMAT (' MAXIMUM CORRECTION TO ', &
     &                                   'ANY POSITION IS',1P,E10.2/ &
     &                                  ' YOU ARE RESPONSIBLE FOR ', &
     &                                  ' DECIDING WHETHER THIS IS A', &
     &                                  ' SERIOUS ERR0R!')
                           END IF
                      END IF
                 END IF
  480       CONTINUE
  490  CONTINUE

!  (3) Compute derivitives of nodal
!      functions at integration points;
!      then check for negative areas:

       CALL Deriv (iUnitT, mxEl, mxNode, nodes, numEl, & ! input
     &             radius, xNode, yNode, &
     &             area, detJ, &                         ! output
     &             dXS, dYS, dXSP, dYSP, fPSfer, sita)
       allOK = .TRUE.
       DO 620 i = 1, numEl
            DO 610 m = 1, 7
                 test = area(i) * detJ(m, i)
                 IF (test <= 0.0D0) THEN
                      IF(Verbose) THEN
					    WRITE(iUnitT,605) m, i
  605                   FORMAT(/' EXCESSIVELY DISTORTED ELEMENT LEADS TO ' &
     &                       ,'NEGATIVE AREA AT POINT ',I1,' IN ELEMENT ',I5)
                        WRITE(iUnitT,606) area(i), detJ(m, i)
  606                   FORMAT('AREA = ',1P,E12.4,'   detJ: ',0P,F12.6)
					  END IF
					allOK = .FALSE.
                 END IF
  610       CONTINUE
  620  CONTINUE
       IF (.NOT.allOK) THEN
	     write(ErrorMsg,'(A)') "allOK was False in MOD_DataSubs"
		 call FatalError(ErrorMsg,ThID)
       END IF

!  (4) Compute lengths of fault elements:

       DO 750 i = 1, nFl
            n1 = nodeF(1, i)
            n2 = nodeF(2, i)
            theta1 = xNode(n1)
            theta2 = xNode(n2)
            phi1  = yNode(n1)
            phi2  = yNode(n2)
            fLen(i) = FltLen (phi1, phi2, radius, theta1, theta2)
  750  CONTINUE

!  (5) Make a list of nodes that are on the boundary and require
!      boundary conditions (nodCon); these are in counterclockwise
!      order.  Also make lists of element sides which contain these
!      nodes: edgeTS and edgeFS.

       nCond = 0
       DO 801 i = 1, numNod
            checkN(i) = .FALSE.
  801  CONTINUE
       DO 802 i = 1, nFl
            edgeFS(1, i) = .FALSE.
            edgeFS(2, i) = .FALSE.
  802  CONTINUE
       DO 810 i = 1, numEl
            DO 809 j = 1, 3
                 CALL Next (i, j, mxEl, mxFEl, nFl, nodeF, nodes, numEl, & ! input
     &                      kFault, kEle)                                  ! output
                 IF (kEle > 0) THEN
!                    (ordinary interior side)
                      edgeTS(j, i) = .FALSE.
                 ELSE IF (kFault == 0) THEN
!                    (exterior side)
                      edgeTS(j, i) = .TRUE.
                      n1 = nodes(MOD(j,  3) + 1, i)
                      n2 = nodes(MOD(j + 1, 3) + 1, i)
                      IF (.NOT.checkN(n1)) THEN
                           nCond = nCond + 1
                           checkN(n1) = .TRUE.
                      END IF
                      IF (.NOT.checkN(n2)) THEN
                           nCond = nCond + 1
                           checkN(n2) = .TRUE.
                      END IF
                 ELSE
!                    (triangular element has an exterior fault element
!                     adjacent to it)
                      edgeTS(j, i) = .FALSE.
                      n1 = nodes(MOD(j,  3) + 1, i)
                      IF (nodeF(2, kFault) == n1) THEN
                           edgeFS(2, kFault) = .TRUE.
                           DO 806 k = 3, 4
                                n = nodeF(k, kFault)
                                IF (.NOT.checkN(n)) THEN
                                     nCond = nCond + 1
                                     checkN(n) = .TRUE.
                                END IF
  806                      CONTINUE
                      ELSE
                           edgeFS(1, kFault) = .TRUE.
                           DO 808 k = 1, 2
                                n = nodeF(k, kFault)
                                IF (.NOT.checkN(n)) THEN
                                     nCond = nCond + 1
                                     checkN(n) = .TRUE.
                                END IF
  808                      CONTINUE
                      END IF
                 END IF
  809       CONTINUE
  810  CONTINUE
       IF (nCond > mxBn) THEN
	     write(ErrorMsg,'(A,I6,A)') "Increase array-size mxBn to at least ",nCond," (by adjusting formula) and recompile."
		 call FatalError(ErrorMsg,ThID)
       END IF

!   Stop work if no boundary nodes found (global grid):

       IF (nCond == 0) GO TO 899

!   Begin circuit with lowest-numbered boundary node
       DO 830 i = 1, numNod
            IF (checkN(i)) GO TO 831
  830  CONTINUE
  831  nodCon(1) = i
       nDone = 1
       nLeft = nCond
!      Beginning of indefinate loop which traces around the perimeter.
!      Each time, it progresses by one of 3 steps:
!      -1 node at a time along a triangle side, OR
!      -1 node at a time along a fault element side, or
!      -by finding another node which shares the same location.
!      Beginning of main indefinate loop:
  840       node = nodCon(ndone)

!           Important: Check that we are not revisiting a node!
!           This would mean that there are too many boundary nodes
!           to fit in the simply-connected loop, and that there
!           are excess boundary nodes somewhere, unconnected!
            IF (.NOT.checkN(node)) THEN
                nGood = nDone - 2
	            write(ErrorMsg,'(A/,A/,A,I6,A/,A,I6,A/,A/,A)') "ERROR IN GRID, reported by -Square-:" , &
         &                 "BOUNDARY IS NOT SIMPLY-CONNECTED."                                     , &
         &                 "Closed loop of ",nGood," nodes does not"                               , &
         &                 "include all ",nCond," boundary nodes."                                 , &
         &                 "Run command PerimeterTest in -OrbWin-"                                 , &
         &                 "for a map of the bad nodes."
		        call FatalError(ErrorMsg,ThID)
            END IF
            IF (nDone > 1) checkN(node) = .FALSE.

            x = xNode(node)
            y = yNode(node)
!           Look for a triangular element with an external
!           side that begins with this node:
            DO 844 i = 1, numEl
                 DO 842 j = 1, 3
                      IF (edgeTS(j, i)) THEN
                           n1 = nodes(MOD(j, 3) + 1, i)
                           IF (n1 == node) GO TO 846
                      END IF
  842            CONTINUE
  844       CONTINUE
            GO TO 850
  846       n2 = nodes(MOD(j + 1, 3) + 1, i)
!           Success by element method: n2 is next boundary node
            nDone = nDone + 1
            IF (nDone <= nCond) nodCon(nDone) = n2
            nLeft = nLeft - 1
            IF (nLeft > 0) THEN
                 GO TO 840
            ELSE
                 GO TO 870
            END IF
!           Else, look for an adjacent fault element using this node:
  850       DO 854 i = 1, nFl
                 IF (edgeFS(1, i)) THEN
                      IF (nodeF(1, i) == node) THEN
                           n2 = nodeF(2, i)
                           GO TO 856
                      END IF
                 ELSE IF (edgeFS(2, i)) THEN
                      IF (nodeF(3, i) == node) THEN
                           n2 = nodeF(4, i)
                           GO TO 856
                      END IF
                 END IF
  854       CONTINUE
            GO TO 860
  856       nDone = nDone + 1
!           Success by fault method: n2 is next boundary node:
            IF (nDone <= nCond) nodCon(nDone) = n2
            nLeft = nLeft - 1
            IF (nLeft > 0) THEN
                 GO TO 840
            ELSE
                 GO TO 870
            END IF
!           Else, look for another exterior corner node at same location:
  860       DO 865 i = 1, numNod
                 IF ((i /= node).AND.checkN(i)) THEN
                      IF ( (ABS(xNode(i) - x) < 1.D-6) .AND. &
     &                     (ABS(yNode(i) - y) < 1.D-6) ) GO TO 867
                 END IF
  865       CONTINUE
	        write(ErrorMsg,'(A/,A,I6/,A/,A/,A)') "BAD GRID TOPOLOGY: WHILE TRACING PERIMETER,", &
     &          "COULD NOT FIND ANY WAY TO CONTINUE FROM NODE ", node                  , &
     &          "EITHER THROUGH SHARED BOUNDARY ELEMENTS, OR"                          , &
     &          "THROUGH OTHER BOUNDARY NODES SHARING THE SAME"                        , &
     &          "POSITION."
		    call FatalError(ErrorMsg,ThID)
  867       nDone = nDone + 1
!           Success by location method: I is the next boundary node
            IF (nDone <= nCond) nodCon(nDone) = i
            nLeft = nLeft - 1
            IF (nLeft > 0) GO TO 840
!      End of indefinite loop which traces around perimeter.
  870  IF (.NOT.skipBC) THEN
            IF(Verbose) WRITE(iUnitT,880)
  880       FORMAT(/ /' Here follows a list, in consecutive order,'/ &
     &                ' of the nodes which define the perimeter'/ &
     &                ' of the model; these nodes require boundary', &
     &                ' conditions:'/'    BC#  Node          ', &
     &                '  Latitude Longitude')
            DO 890 i = 1, nCond
                 n = nodCon(i)
                 theLon = yNode(n) * 57.2957795130823D0
                 theLat = 90.0D0 - xNode(n) * 57.2957795130823D0
                 IF(Verbose) WRITE(iUnitT,882) i, n, theLat, theLon
  882            FORMAT(' ',2I6,10X,2F10.3)
  890       CONTINUE
            n = nodCon(1)
            IF(Verbose) WRITE(iUnitT,892) n
  892       FORMAT(' (Note: Initial node ',I6,' completes the loop,', &
     &              ' but is not listed again.)')
       END IF
  899  CONTINUE

!  (6)  Survey fault elements and issue warning if any element is of
!       mixed type (part strike-slip, and part shallow-dipping):

       DO 920 i = 1, nFl
            deld1 = fDip(1, i) - 1.57079632679490D0
            deld2 = fDip(2, i) - 1.57079632679490D0
            vert1 = ABS(deld1) <= wedge
            vert2 = ABS(deld2) <= wedge
            nvpart = 0
            IF (vert1) THEN
                 nvpart = nvpart + 1
                 tag1 = vertic
            ELSE
                 tag1 = obliqu
            END IF
            IF (vert2) THEN
                 nvpart = nvpart + 1
                 tag2 = vertic
            ELSE
                 tag2 = obliqu
            END IF
            switch = ((nvpart > 0).AND.(nvpart < 2))
            IF (switch) THEN
                 dip1 = fDip(1, i) * 57.2957795130823D0
                 IF (dip1 > 90.0D0) dip1 = dip1 - 180.0D0
                 dip2 = fDip(2, i) * 57.2957795130823D0
                 IF (dip2 > 90.0D0) dip2 = dip2 - 180.0D0
                 IF(Verbose) WRITE(iUnitT,905) i, dip1, tag1, dip2, tag2
  905            FORMAT(/ /' CAUTION:'/ &
     &                ' FAULT ELEMENT ',I6,' HAS DIPS OF '/ &
     &                ' ',F7.2,' DEGREES ',A21/ &
     &                ' ',F7.2,' DEGREES ',A21/ &
     &                ' WHICH MAKES IT MIXED-MODE.'/ &
     &             ' SUCH ELEMENTS ARE INACCURATE AND NOT RECOMMENDED.'/ &
     &           ' PREFERABLY EACH ELEMENT SHOULD BE OF A SINGLE TYPE.'/ &
     &           ' (REMEMBER, DIP NEED NOT BE CONTINUOUS FROM ONE', &
     &                ' FAULT ELEMENT TO THE NEXT.)')
            ELSE
                 nvpart = 0
                 DO 910 m = 1, 7
                      deld = deld1 * fPhi(1, m) + deld2 * fPhi(2, m)
                      IF (ABS(deld) <= wedge) nvpart = nvpart + 1
  910            CONTINUE
                 IF ((nvpart > 0).AND.(nvpart < 7)) THEN
                      IF (nvpart >= 4) THEN
                           IF(Verbose) WRITE(iUnitT,912) i, dip1, dip2
  912                      FORMAT(/ /' CAUTION:'/ &
     &                          ' FAULT ELEMENT ',I6,' HAS DIPS OF '/ &
     &                          ' ',F7.2,' DEGREES, AND'/ &
     &                          ' ',F7.2,' DEGREES'/ &
     &                     ' WHICH APPEAR TO MAKE IT STRIKE-SLIP.'/ &
     &                ' HOWEVER, THESE VALUES ARE SUCH THAT DIP-SLIP'/ &
     &               ' IS PERMITTED AT ONE OR MORE INTEGRATION POINTS.'/ &
     &             ' SUCH ELEMENTS ARE INACCURATE AND NOT RECOMMENDED.'/ &
     &           ' PREFERABLY EACH ELEMENT SHOULD BE OF A SINGLE TYPE.'/ &
     &           ' (REMEMBER, DIP NEED NOT BE CONTINUOUS FROM ONE', &
     &                               ' FAULT ELEMENT TO THE NEXT.)')
                      ELSE
                           IF(Verbose) WRITE(iUnitT,914) i, dip1, dip2
  914                      FORMAT(/ /' CAUTION:'/ &
     &                          ' FAULT ELEMENT ',I6,' HAS DIPS OF '/ &
     &                          ' ',F7.2,' DEGREES, AND'/ &
     &                          ' ',F7.2,' DEGREES'/ &
     &                     ' WHICH APPEAR TO MAKE IT FREE-SLIPPING.'/ &
     &                ' HOWEVER, THESE VALUES ARE SUCH THAT DIP-SLIP'/ &
     &              ' IS PROHIBITED AT ONE OR MORE INTEGRATION POINTS.'/ &
     &             ' SUCH ELEMENTS ARE INACCURATE AND NOT RECOMMENDED.'/ &
     &           ' PREFERABLY EACH ELEMENT SHOULD BE OF A SINGLE TYPE.'/ &
     &           ' (REMEMBER, DIP NEED NOT BE CONTINUOUS FROM ONE', &
     &                     ' FAULT ELEMENT TO THE NEXT.)')
                      END IF
                 END IF
            END IF
  920  CONTINUE

!  (7)  Calculate fault argument (in radians, measured counterclockwise
!       from +Theta = South) at each end of each fault element.

       DO 1000 i = 1, nFl
            n1 = nodeF(1, i)
            n2 = nodeF(2, i)
            theta(1) = xNode(n1)
            theta(2) = xNode(n2)
            phi(1)  = yNode(n1)
            phi(2)  = yNode(n2)
            CALL FAngls(phi, theta, & ! input
     &                  fAngle)       ! output
            DO 900 j = 1, 2
                 fArg(j, i) = fAngle(j)
  900       CONTINUE
 1000  CONTINUE

!  (8) Survey strike-slip (vertical) faults to check for conflicts in
!      argument that would lock the fault:

       IF (log_strike_adjustments .AND. Verbose) WRITE(iUnitT,1001)
 1001  FORMAT(/ /' The following tightly-connected pairs of strike-slip' &
     &          /' fault elements had their azimuths averaged at the' &
     &          /' connection point for purposes of computing the' &
     &          /' constraint on the direction of strike-slip:' &
     &        / /' Fault#1   Fault#2    Node#A    Node#B   ', &
     &            '  Latitude Longitude    Azim#1    Azim#2   Azimuth' &
     &          /' ----------------------------------------', &
     &            '--------------------------------------------------')
!      Loop on all fault elements (I):
       DO 2000 i = 1, nFl
!         Loop on 2 terminal node pairs, 1-4, 2-3 (J = 1 or 2):
          DO 1900 j = 1, 2
!            Dip must be within "wedge" of vertical for constraint:
             IF (ABS(fDip(j, i) - 1.57079632679490D0) <= wedge) THEN
                nazi = j
                n1 = j
                IF(j == 1) THEN
                   n4 = 4
                ELSE
                   n4 = 3
                END IF
                node1 = nodeF(n1, i)
                node4 = nodeF(n4, i)
!               No constraint applied where a fault ends:
                IF (node1 /= node4) THEN
!                  Endpoint pairs must be checkEd for duplication:
!                  Look for other strike-slip faults sharing this
!                  pair of nodes, at either end:
                   found = .FALSE.
                   DO 1600 l = 1, nFl
                      IF (l /= i) THEN
                         IF (ABS(fDip(1, l) - 1.57079632679490D0) <= wedge) THEN
                            IF (((node1 == nodeF(1, l)).AND. &
     &                           (node4 == nodeF(4, l))).OR. &
     &                          ((node1 == nodeF(4, l)).AND. &
     &                           (node4 == nodeF(1, l)))) THEN
                               found = .TRUE.
                               number = l
                               nazl = 1
                               GO TO 1601
                            END IF
                         END IF
                         IF (ABS(fDip(2, l) - 1.57079632679490D0) <= wedge) THEN
                            IF (((node1 == nodeF(2, l)).AND. &
     &                           (node4 == nodeF(3, l))).OR. &
     &                          ((node1 == nodeF(3, l)).AND. &
     &                           (node4 == nodeF(2, l)))) THEN
                               found = .TRUE.
                               number = l
                               nazl = 2
                               GO TO 1601
                            END IF
                         END IF
                      END IF
 1600              CONTINUE
!                  Don't worry if this pair already checkEd!
 1601              IF (found.AND.(number > i)) THEN
!                     Average arguments together (avoid cycle shifts):
                      IF(nazi == nazl) THEN
                         azl = fArg(nazl, number) + 3.14159265358979D0
                      ELSE
                         azl = fArg(nazl, number)
                      END IF
                      azi = fArg(nazi, i)
                      cosz = 0.5D0 * (COS(azi) + COS(azl))
                      sinz = 0.5D0 * (SIN(azi) + SIN(azl))
                      azimut = ATan2F(sinz, cosz)
                      fArg(nazi, i) = azimut
                      IF(nazl == nazi) THEN
                         fArg(nazl, number) = azimut - 3.14159265358979D0
                      ELSE
                         fArg(nazl, number) = azimut
                      END IF
!                    Print a warning:
                      dazi = azi * 57.2957795130823D0
                      dazl = azl * 57.2957795130823D0
                      daz = azimut * 57.2957795130823D0
                      np1 = node1
                      np4 = node4
                      dELon = 57.2957795130823D0 * yNode(node1)
                      dNLat = 90.0D0 - 57.2957795130823D0 * xNode(node1)
                      IF (log_strike_adjustments .AND. Verbose) WRITE(iUnitT,1610) &
     &                       i, number, np1, np4, &
     &                       dNLat, dELon, dazi, dazl, daz
 1610                 FORMAT(' ',I7,3X,I7,3X, &
     &                           I7,3X,I7,3X, &
     &                           2X,F8.3,1X,F9.3, &
     &                           4X,F6.1,4X,F6.1,4X,F6.1)
                   END IF
!                  ^End block which looks for constraints
                END IF
!               ^End block which checks for distinct node numbers
             END IF
!            ^End block which checks for dip of over 75 degrees
 1900     CONTINUE
!         ^End loop on 2 node pairs in fault element
 2000  CONTINUE

!  (9) Calculate nodal functions at integration points on faults:

       CALL FNodal (mxFEl, mxNode, nFl, nodeF, xNode, yNode, & ! input
     &              fPFlt)                                     ! output

       IF (Verbose) WRITE(iUnitT,9999)
 9999  FORMAT (' --------------------------------------------------', &
     &          '-----------------------------')
       END SUBROUTINE Square

       SUBROUTINE Deriv (iUnitT, mxEl, mxNode, &  ! input
     &                   nodes, numEl, &
     &                   radius, xNode, yNode, &
     &                   area, detJ, &            ! output
     &                   dXS, dYS, dXSP, dYSP, fPSfer, sita)

!   Sets up 6 vector nodal functions (fPSfer) of each spherical
!     triangle finite element, at each of its 7 integration points.
!   Calculates dXS and dYS, the Theta-derivitive and Phi-derivitive
!     of each of these 6 vector nodal functions.
!   Also computes "area", the areas of the plane triangles.
!   Also computes "detJ", the local ratio of areas on the sphere
!     to areas on the plane triangles.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
       INTEGER, INTENT(IN) :: iUnitT, mxEl, mxNode, nodes, numEl                               ! input
       REAL*8, INTENT(IN) :: radius, xNode, yNode                                              ! input
       REAL*8, INTENT(OUT) :: area, detJ, dXS, dYS, dXSP, dYSP, fPSfer, sita                   ! output
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION points
       DOUBLE PRECISION fff, skkc, skke, sncsne, snccse, csccse, cscsne
       DOUBLE PRECISION xa, xb, xc, ya, yb, yc, za, zb, zc, xyzp
       INTEGER i, j, m
       REAL*8 a, areaP, b, c, cka, cosm, cscs, cse, cssn, &
            & dd, dd1, dd2, dd3, ddpn, dpdc, dpde, &
            & pfq, phaij, phi, pnx, pny, pnz, pp, rn, rr1, rr2, rr3, &
            & sitaj, sitami, snc, sne, theta, tx, ty, &
            & x21, x31, y21, y31, z21, z31
       DIMENSION xNode(mxNode), yNode(mxNode), nodes(3, mxEl), area(mxEl)
       DIMENSION detJ(7, mxEl)
       DIMENSION dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl)
       DIMENSION dXSP(3, 7, mxEl), dYSP(3, 7, mxEl), points(3, 7)
       DIMENSION phi(3), theta(3), skkc(3), skke(3), fff(3), &
               & sita(7, mxEl), fPSfer(2, 2, 3, 7, mxEl)

       COMMON / S1S2S3 / points

      DO 900 i = 1, numEl
         DO 100 j = 1, 3
            theta(j) = xNode(nodes(j, i))
            phi(j)  = yNode(nodes(j, i))
  100    CONTINUE
         x21 = SIN(theta(2)) * COS(phi(2)) - SIN(theta(1)) * COS(phi(1))
         x31 = SIN(theta(3)) * COS(phi(3)) - SIN(theta(1)) * COS(phi(1))
         y21 = SIN(theta(2)) * SIN(phi(2)) - SIN(theta(1)) * SIN(phi(1))
         y31 = SIN(theta(3)) * SIN(phi(3)) - SIN(theta(1)) * SIN(phi(1))
         z21 = COS(theta(2)) - COS(theta(1))
         z31 = COS(theta(3)) - COS(theta(1))
         a = y21 * z31 - y31 * z21
         b = z21 * x31 - z31 * x21
         c = x21 * y31 - x31 * y21
         areaP = SQRT(a * a + b * b + c * c)
         area(i) = radius * radius * (0.5D0 * areaP)
         pnx = a / areaP
         pny = b / areaP
         pnz = c / areaP
         dd1 = SIN(theta(1)) * COS(phi(1)) * pnx
         dd2 = SIN(theta(1)) * SIN(phi(1)) * pny
         dd3 = COS(theta(1)) * pnz
         dd = dd1 + dd2 + dd3

! This part is to test if Kong's method and Bird's method give the same
! results for the derivitive:

         xa = SIN(theta(1)) * COS(phi(1))
         xb = SIN(theta(2)) * COS(phi(2))
         xc = SIN(theta(3)) * COS(phi(3))
         ya = SIN(theta(1)) * SIN(phi(1))
         yb = SIN(theta(2)) * SIN(phi(2))
         yc = SIN(theta(3)) * SIN(phi(3))
         za = COS(theta(1))
         zb = COS(theta(2))
         zc = COS(theta(3))
         cka = (yb * zc - zb * yc) * xa + (zb * xc - xb * zc) * ya + (xb * yc - yb * xc) * za

         DO 800 m = 1, 7
            snccse = 0.0D0
            sncsne = 0.0D0
            cosm = 0.0D0
            DO 200 j = 1, 3
               snccse = snccse + points(j, m) * SIN(theta(j)) * COS(phi(j))
               sncsne = sncsne + points(j, m) * SIN(theta(j)) * SIN(phi(j))
               cosm = cosm + points(j, m) * COS(theta(j))
  200       CONTINUE
            xyzp = SQRT(snccse * snccse + sncsne * sncsne + cosm * cosm)
            snccse = snccse / xyzp
            sncsne = sncsne / xyzp
            cosm = cosm / xyzp
            sitaj = ACOS(cosm)
            ty = sncsne
            tx = snccse
            phaij = ATan2F(ty, tx)
            csccse = COS(sitaj) * COS(phaij)
            cscsne = COS(sitaj) * SIN(phaij)

!   Bird's method:

            fff(1) = ((yb * zc - zb * yc) * snccse + (zb * xc - xb * zc) * sncsne &
     &             + (xb * yc - yb * xc) * cosm) / cka
            fff(2) = ((yc * za - zc * ya) * snccse + (zc * xa - xc * za) * sncsne &
     &             + (xc * ya - yc * xa) * cosm) / cka
            fff(3) = ((ya * zb - za * yb) * snccse + (za * xb - xa * zb) * sncsne &
     &             + (xa * yb - ya * xb) * cosm) / cka
            skkc(1) = ((yb * zc - zb * yc) * csccse &
     &              + (zb * xc - xb * zc) * cscsne &
     &              - (xb * yc - yb * xc) * SIN(sitaj)) / cka
            skkc(2) = ((yc * za - zc * ya) * csccse &
     &              + (zc * xa - xc * za) * cscsne &
     &              - (xc * ya - yc * xa) * SIN(sitaj)) / cka
            skkc(3) = ((ya * zb - za * yb) * csccse &
     &              + (za * xb - xa * zb) * cscsne &
     &              - (xa * yb - ya * xb) * SIN(sitaj)) / cka
            skke(1) = (-(yb * zc - zb * yc) * sncsne &
     &               + (zb * xc - xb * zc) * snccse) / cka
            skke(2) = (-(yc * za - zc * ya) * sncsne &
     &               + (zc * xa - xc * za) * snccse) / cka
            skke(3) = (-(ya * zb - za * yb) * sncsne &
     &               + (za * xb - xa * zb) * snccse) / cka

            sita(m, i) = sitaj
            rr1 = SIN(sitaj) * COS(phaij)
            rr2 = SIN(sitaj) * SIN(phaij)
            rr3 = COS(sitaj)
            rn = rr1 * pnx + rr2 * pny + rr3 * pnz
            pp = dd / rn
            dpdc = (COS(sitaj) * COS(phaij) * pnx + COS(sitaj) * SIN(phaij) * pny &
     &           - SIN(sitaj) * pnz)
            dpde = (-SIN(sitaj) * SIN(phaij) * pnx + &
     &             SIN(sitaj) * COS(phaij) * pny)
            ddpn = pp / rn
            dpdc = -ddpn * dpdc
            dpde = -ddpn * dpde
            IF(sita(m, i) <= 0.0D0.OR.sita(m, i) >= 3.14159265358979D0) THEN
               sitami = sita(m, i) * 57.2957795130823D0
	           write(ErrorMsg,'(A,I5,A,I5,A,D14.4)') "COLATITUDE OF INTEGRATION POINT ",m," OF ELEMENT ",i," IS OUT RANGE ",sitami
		       call FatalError(ErrorMsg,ThID)
            END IF
            DO 500 j = 1, 3
               dXSP(j, m, i) = dpdc * fff(j) + pp * skkc(j)
               dYSP(j, m, i) = dpde * fff(j) + pp * skke(j)
               cscs = COS(theta(j)) * COS(phi(j))
               cssn = COS(theta(j)) * SIN(phi(j))
               snc = SIN(theta(j))
               sne = SIN(phi(j))
               cse = COS(phi(j))
               fPSfer(1, 1, j, m, i) = cscs * csccse + cssn * cscsne &
     &                          + snc * SIN(sitaj)
               fPSfer(2, 1, j, m, i) = -sne * csccse + cse * cscsne
               fPSfer(1, 2, j, m, i) = -cscs * SIN(phaij) + cssn * COS(phaij)
               fPSfer(2, 2, j, m, i) = sne * SIN(phaij) + cse * COS(phaij)
               dXS(1, 1, j, m, i) = (-cscs * snccse - cssn * sncsne &
     &                           + snc * COS(sitaj)) * fff(j) &
     &                          + fPSfer(1, 1, j, m, i) * skkc(j)
               dXS(2, 1, j, m, i) = (sne * snccse - cse * sncsne) * fff(j) &
     &                          + fPSfer(2, 1, j, m, i) * skkc(j)
               dYS(1, 1, j, m, i) = (-cscs * cscsne + cssn * csccse) * fff(j) &
     &                          + fPSfer(1, 1, j, m, i) * skke(j)
               dYS(2, 1, j, m, i) = (sne * cscsne + cse * csccse) * fff(j) &
     &                          + fPSfer(2, 1, j, m, i) * skke(j)
               dXS(1, 2, j, m, i) = fPSfer(1, 2, j, m, i) * skkc(j)
               dXS(2, 2, j, m, i) = fPSfer(2, 2, j, m, i) * skkc(j)
               dYS(1, 2, j, m, i) = (-cscs * COS(phaij) - cssn * SIN(phaij)) &
     &                          * fff(j) &
     &                          + fPSfer(1, 2, j, m, i) * skke(j)
               dYS(2, 2, j, m, i) = (sne * COS(phaij) - cse * SIN(phaij)) &
     &                          * fff(j) &
     &                          + fPSfer(2, 2, j, m, i) * skke(j)
               fPSfer(1, 1, j, m, i) = fPSfer(1, 1, j, m, i) * fff(j)
               fPSfer(2, 1, j, m, i) = fPSfer(2, 1, j, m, i) * fff(j)
               fPSfer(1, 2, j, m, i) = fPSfer(1, 2, j, m, i) * fff(j)
               fPSfer(2, 2, j, m, i) = fPSfer(2, 2, j, m, i) * fff(j)
  500       CONTINUE
            pfq = fff(1) + fff(2) + fff(3) ! orphan statement, left over from some test?  (pfq does not seem to be used.)
            detJ(m, i) = rn**3 / (dd * dd)
  800    CONTINUE
  900 CONTINUE
      END SUBROUTINE Deriv

       SUBROUTINE Next (i, j, mxEl, mxFEl, nFl, nodeF, nodes, numEl, & ! input
     &                  kFault, kEle)                                  ! output

!  Determine whether there are more elements adjacent to side #j of
!  triangular continuum element #i.
!    j = 1 means the side opposite node # nodes(1, i).
!    j = 2 means the side opposite node # nodes(2, i).
!    j = 3 means the side opposite node # nodes(3, i).
!  If a fault element is adjacent, its number is kFault; otherwise,
!  kFault is set to zero.
!  If another triangular continuum element is adjacent (even across
!  fault element kFault !) then its number is kEle; otherwise, kEle = 0.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER, INTENT(IN) :: i, j, mxEl, mxFEl, nFl, nodeF, nodes, numEl                     ! input
       INTEGER, INTENT(OUT) :: kFault, kEle                                                   ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER k, khigh, klow, l, m1, m2, m3, m4, n1, n2
       LOGICAL foundF
       DIMENSION nodeF(4, mxFEl), nodes(3, mxEl)

!  Two node numbers along the side of interest, counterclockwise:
       n1 = nodes(MOD(j,  3) + 1, i)
       n2 = nodes(MOD(j + 1, 3) + 1, i)
!  Check for adjacent fault element first:
       foundF = .FALSE.
       kFault = 0
       IF (nFl > 0) THEN
            DO 10 k = 1, nFl
                 m1 = nodeF(1, k)
                 m2 = nodeF(2, k)
                 m3 = nodeF(3, k)
                 m4 = nodeF(4, k)
                 IF (((m1 == n2).AND.(m2 == n1)).OR. &
     &               ((m3 == n2).AND.(m4 == n1))) THEN
                      foundF = .TRUE.
                      kFault = k
                      GO TO 11
                 END IF
   10       CONTINUE
   11       IF (.NOT.foundF) kFault = 0
!  If there was a fault, replace 2 node numbers that we search for:
            IF (foundF) THEN
                 IF (m2 == n1) THEN
                      n1 = m3
                      n2 = m4
                 ELSE
                      n1 = m1
                      n2 = m2
                 END IF
            END IF
       END IF
!  Search for adjacent triangular continuum element:
       kEle = 0
       klow = i
       khigh = i
!  --- Begin irregular loop, to search out nearest elements first ---
  100  klow = klow - 1
            IF (klow >= 1) THEN
                 DO 110 l = 1, 3
                      m1 = nodes(MOD(l,  3) + 1, klow)
                      m2 = nodes(MOD(l + 1, 3) + 1, klow)
                      IF ((m2 == n1).AND.(m1 == n2)) THEN
                           kEle = klow
                           RETURN
                      END IF
  110            CONTINUE
            END IF
            khigh = khigh + 1
            IF (khigh <= numEl) THEN
                 DO 120 l = 1, 3
                      m1 = nodes(MOD(l,  3) + 1, khigh)
                      m2 = nodes(MOD(l + 1, 3) + 1, khigh)
                      IF ((m2 == n1).AND.(m1 == n2)) THEN
                           kEle = khigh
                           RETURN
                      END IF
  120            CONTINUE
            END IF
       IF ((klow > 1).OR.(khigh < numEl)) GO TO 100
       END SUBROUTINE Next

       REAL*8 FUNCTION FltLen (phi1, phi2, radius, theta1, theta2) ! input

!      Calculates length of great circle segment between
!      point (theta1, phi1) and point (theta2, phi2),
!      in physical length units (radians * radius).
!      Note that theta is colatitude (from North pole),
!      and phi is East longitude; both in radians.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: phi1, phi2, radius, theta1, theta2                                 ! inputs
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION ab
       ab = SIN(theta1) * SIN(theta2) * COS(phi1) * COS(phi2) + &
     &      SIN(theta1) * SIN(theta2) * SIN(phi1) * SIN(phi2) + &
     &      COS(theta1) * COS(theta2)
       ab = ACOS(ab)
       FltLen = ab * radius
       END FUNCTION FltLen


      SUBROUTINE FAngls (phi, theta, & ! input
     &                   fAngle)       ! output

!  Calculate the arguments (angles counterclockwise from +Theta)
!  at both ends of an arc of a great circle.  Results in radians.

      IMPLICIT NONE
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL*8, INTENT(IN) :: phi, theta                                                       ! input
      REAL*8, INTENT(OUT) :: fAngle                                                          ! output
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION fPoint
      COMMON / SFault / fPoint
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL*8 a1, a2, a3, b1, b2, b3, dg180, dx, dxx, dy, dyy, dz, &
           & p1, p2, phai, s, s1, s2, s3, sita, xval, xx, yy, zz
      DIMENSION fAngle(2), fPoint(7), phi(2), theta(2)

      dg180 = 3.14159265358979D0
      a1 = SIN(theta(1)) * COS(phi(1))
      a2 = SIN(theta(1)) * SIN(phi(1))
      a3 = COS(theta(1))
      b1 = SIN(theta(2)) * COS(phi(2))
      b2 = SIN(theta(2)) * SIN(phi(2))
      b3 = COS(theta(2))

      s = 0.99D0
      xx = s * a1 + (1.0D0 - s) * b1
      yy = s * a2 + (1.0D0 - s) * b2
      zz = s * a3 + (1.0D0 - s) * b3
      xval = SQRT(xx * xx + yy * yy + zz * zz)
      xx = xx / xval
      yy = yy / xval
      zz = zz / xval
      dx = xx - a1
      dy = yy - a2
      dz = zz - a3
      sita = theta(1)
      phai = phi(1)
      s1 = COS(sita) * COS(phai)
      s2 = COS(sita) * SIN(phai)
      s3 = -SIN(sita)
      p1 = -SIN(phai)
      p2 = COS(phai)
      dxx = dx * s1 + dy * s2 + dz * s3
      dyy = dx * p1 + dy * p2
      fAngle(1) = ATan2F(dyy, dxx)

      s = 0.01D0
      xx = s * a1 + (1.0D0 - s) * b1
      yy = s * a2 + (1.0D0 - s) * b2
      zz = s * a3 + (1.0D0 - s) * b3
      xval = SQRT(xx * xx + yy * yy + zz * zz)
      xx = xx / xval
      yy = yy / xval
      zz = zz / xval
      dx = b1 - xx
      dy = b2 - yy
      dz = b3 - zz
      sita = ACOS(zz)
      phai = ATan2F(yy, xx)
      IF(phai < 0.0) phai = 2.0D0 * dg180 + phai
      s1 = COS(sita) * COS(phai)
      s2 = COS(sita) * SIN(phai)
      s3 = -SIN(sita)
      p1 = -SIN(phai)
      p2 = COS(phai)
      dxx = dx * s1 + dy * s2 + dz * s3
      dyy = dx * p1 + dy * p2
      fAngle(2) = ATan2F(dyy, dxx)

      END SUBROUTINE FAngls

       SUBROUTINE FNodal (mxFEl, mxNode, nFl, nodeF, & ! input
     &                    xNode, yNode, &
     &                    fPFlt)                       ! output

!   Calculates vector nodal functions at all integration points
!   on all arc-of-great-circle fault elements.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER, INTENT(IN) :: mxFEl, mxNode, nFl, nodeF                        ! input
       REAL*8, INTENT(IN) :: xNode, yNode                                      ! input
       REAL*8, INTENT(OUT) :: fPFlt                                            ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION fPhi
       COMMON / FPhis / fPhi
       DIMENSION fPhi(4, 7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER i, j, k, l, m, n1, n2
       REAL*8 fpp, phi, theta
       DIMENSION fPFlt(2, 2, 2, 7, mxFEl), fpp(2, 2, 2, 7), &
     &           nodeF(4, mxFEl), phi(2), theta(2), &
     &           xNode(mxNode), yNode(mxNode)

       DO 900 i = 1, nFl
            n1 = nodeF(1, i)
            n2 = nodeF(2, i)
            theta(1) = xNode(n1)
            theta(2) = xNode(n2)
            phi(1) = yNode(n1)
            phi(2) = yNode(n2)
            CALL SNodal (phi, theta, & ! input
                      &  fpp)          ! output
            DO 800 m = 1, 7
                 DO 500 j = 1, 2
                      DO 400 k = 1, 2
                           DO 300 l = 1, 2
                                fPFlt(l, k, j, m, i) = fpp(l, k, j, m)
  300                      CONTINUE ! l = 1:2
  400                 CONTINUE ! k = 1:s
  500            CONTINUE ! j = 1:2
  800       CONTINUE ! m = 1:7
  900  CONTINUE ! i = 1:nFl
       END SUBROUTINE FNodal

       SUBROUTINE SNodal (phi, theta, & ! input
     &                    fpp)          ! output

!   Calculates all (vector) nodal functions at all integration points along an
!   arc-of-great-circle side of any single finite element.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: phi, theta                                                       ! input
       REAL*8, INTENT(OUT) :: fpp                                                             ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION fPhi
       COMMON / FPhis / fPhi
       DIMENSION fPhi(4, 7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER j, m
       REAL*8 cscs, csp, cssn, dd, fp, &
            & phaij, ppm, rn, sitaj, snc, snp, snsn, &
            & x1, x2, xn, xx, xyzn, y1, y2, yn, yy, z1, z2, zn, zz
       DOUBLE PRECISION pp
       DIMENSION fpp(2, 2, 2, 7), phi(2), theta(2)

       x1 = SIN(theta(1)) * COS(phi(1))
       y1 = SIN(theta(1)) * SIN(phi(1))
       z1 = COS(theta(1))
       x2 = SIN(theta(2)) * COS(phi(2))
       y2 = SIN(theta(2)) * SIN(phi(2))
       z2 = COS(theta(2))
       xn = x1 + x2
       yn = y1 + y2
       zn = z1 + z2
       xyzn = SQRT(xn * xn + yn * yn + zn * zn)
       xn = xn / xyzn
       yn = yn / xyzn
       zn = zn / xyzn
       dd = x1 * xn + y1 * yn + z1 * zn
       DO 800 m = 1, 7
          xx = fPhi(1, m) * x1 + fPhi(2, m) * x2
          yy = fPhi(1, m) * y1 + fPhi(2, m) * y2
          zz = fPhi(1, m) * z1 + fPhi(2, m) * z2
          pp = SQRT(xx * xx + yy * yy + zz * zz)
          xx = xx / pp
          yy = yy / pp
          zz = zz / pp
          sitaj = ACOS(zz)
          phaij = ATan2F(yy, xx)
          rn = xx * xn + yy * yn + zz * zn
          ppm = rn / dd
          cscs = COS(sitaj) * COS(phaij)
          cssn = COS(sitaj) * SIN(phaij)
          snsn = SIN(sitaj) * SIN(phaij)
          snc = SIN(sitaj)
          snp = SIN(phaij)
          csp = COS(phaij)
          DO 500 j = 1, 2
             fp = fPhi(j, m) * ppm
             fpp(1, 1, j, m) = ( COS(theta(j)) * COS(phi(j)) * cscs &
     &                         + COS(theta(j)) * SIN(phi(j)) * cssn &
     &                         + SIN(theta(j)) * snc) * fp
             fpp(2, 1, j, m) = (-SIN(phi(j)) * cscs + COS(phi(j)) * cssn) * fp
             fpp(1, 2, j, m) = (-COS(theta(j)) * COS(phi(j)) * snp &
     &                         + COS(theta(j)) * SIN(phi(j)) * csp) * fp
             fpp(2, 2, j, m) = ( SIN(phi(j)) * snp &
     &                         + COS(phi(j)) * csp) * fp
  500     CONTINUE
  800  CONTINUE
       END SUBROUTINE SNodal

       SUBROUTINE GetNet (iUnit7, iUnitT, &                         ! input
     &                    mxDOF, mxEl, mxFEl, mxNode, &             ! input
     &                    brief, &                                  ! output
     &                    continuum_LRi, &                          ! output
     &                    dQdTdA, elev, &                           ! output
     &                    fault_LRi, fDip, &                        ! output
     &                    nFakeN, nFl, nodeF, nodes, nRealN, &      ! output
     &                    numEl, numNod, n1000, offMax, offset, &   ! output
     &                    title1, xNode, yNode, &                   ! output
     &                    checkE, checkF, checkN)                   ! work

!   Read finite element grid from unit iUnit7 (assumed already OPENed).
!   Echoes the important values to unit iUnitT (assumed already OPENed).

!   NOTE that this version is different from GetNet in - Shells_v4.1+ -
!   because it does NOT even TRY to read the last 4 nodal data
!   (crustal thickness, mantle-lithosphere thicknesss,
!    density_anomaly, & cooling_curvature)
!   which may follow the nodal values of elevation and heat-flow.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER, INTENT(IN) :: iUnit7, iUnitT, mxDOF, mxEl, mxFEl, mxNode                       ! input
       LOGICAL, INTENT(OUT) :: brief                                                           ! output
       INTEGER, INTENT(OUT) :: continuum_LRi                                                   ! output
       REAL*8,  INTENT(OUT) :: dQdTdA, elev                                                    ! output
       INTEGER, INTENT(OUT) :: fault_LRi                                                       ! output
       REAL*8,  INTENT(OUT) :: fDip                                                            ! output
       INTEGER, INTENT(OUT) :: nFakeN, nFl, nodeF, nodes, nRealN, numEl, numNod, n1000         ! output
       REAL*8,  INTENT(OUT) :: offMax, offset                                                  ! output
       CHARACTER*100, INTENT(OUT) :: title1                                                     ! output
       REAL*8,  INTENT(OUT) :: xNode, yNode                                                    ! output
       LOGICAL checkE, checkF, checkN                                                          ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CHARACTER*80 :: longer_line, shorter_line
       LOGICAL allOK
       INTEGER i, index, j, k, l, LRi, n, nrt2
       REAL*8 dips, elevi, &
            & off, pLat, pLon, qi, vector, xi, yi
       DIMENSION checkE(mxEl), checkF(mxFEl), checkN(mxNode), &
     &           continuum_LRi(mxEl), &
     &           dQdTdA(mxNode), elev(mxNode), &
     &           fault_LRi(mxFEl), &
     &           fDip(2, mxFEl), nodeF(4, mxFEl), &
     &           nodes(3, mxEl), offset(mxFEl), &
     &           xNode(mxNode), yNode(mxNode)
       DIMENSION dips(3), vector(9)

       IF(Verbose) WRITE(iUnitT,1) iunit7
    1  FORMAT(//' Attempting to read finite-element grid (.FEG) file from unit ', I3/)
       title1 = ' '
       READ (iunit7, 2) title1
    2  FORMAT (A80)
       IF(Verbose) WRITE (iUnitT, 3) title1
    3  FORMAT(/' Title of finite-element grid ='/' ',A80)

!   Read number of nodes, plus out-dated parameters that once
!     permitted boundary nodes to be specially numbered as
!    "fake" nodes with numbers from n1000+1 ... n1000+nFakeN.
!     This option is no longer supported by my programs!
!    (Option "brief" suppresses most output.)

       READ (iunit7, * ) numNod, nRealN, nFakeN, n1000, brief
       brief = .TRUE.

       IF (numNod /= (nRealN + nFakeN)) THEN
		 write(ErrorMsg,'(A,I0,A,I0,A,I0,A)') "numNod (",numNod,") IS NOT EQUAL TO SUM OF nRealN (",nRealN,") AND nFakeN (",nFakeN,")."
		 call FatalError(ErrorMsg,ThID)
       END IF

       IF (nRealN > n1000) THEN
		 write(ErrorMsg,'(A,I0,A,I0,A)') "nRealN (",nRealN,") IS GREATER THAN n1000 (",n1000,")."
		 call FatalError(ErrorMsg,ThID)
       END IF

       IF (numNod > mxNode) THEN
		 write(ErrorMsg,'(A,I0,A)') "INCREASE ARRAY-SIZE maxNod TO BE AT LEAST THE NUMBER OF NODES (",numNod,") AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
       END IF

       nrt2 = nRealN * 2
       IF (nrt2 > mxDOF) THEN
		 write(ErrorMsg,'(A,I0,A)') "INCREASE ARRAY-SIZE maxDOF TO ",nrt2," AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
       END IF

       IF (brief) THEN
            IF(Verbose) WRITE (iUnitT, 35)
   35       FORMAT(/' (Since option ""brief"" = .TRUE., grid will not be echoed here.)')
       ELSE
            IF(Verbose) WRITE (iUnitT, 40) numNod
   40       FORMAT (/' There are',I5,' nodes in the grid')
            IF(Verbose) WRITE (iUnitT, 50)
   50       FORMAT (/ &
     &                     77X,'                mantle'/ &
     &                     77X,'   crustal lithosphere'/ &
     &               '       node E-longitude  N-latitude', &
     &               '      theta        phi elevation', &
     &               ' heat-flow thickness   thickness'/)
       END IF
       DO 90 k = 1, numNod
            checkN(k) = .FALSE.
   90  CONTINUE

       DO 100 k = 1, numNod

            CALL ReadN (iunit7, iUnitT, 5, & ! input  <=== NOTE! Only reading first 5 numbers in each nodal record (not 9).
     &                  vector)              ! output
            index = vector(1) + 0.5D0 ! truncated to INTEGER
            IF (index > nRealN) THEN
                 IF ((index <= n1000).OR. (index > (n1000 + nFakeN))) THEN
				   write(ErrorMsg,'(A,I0)') "ILLEGAL NODE NUMBER: ",index
				   call FatalError(ErrorMsg,ThID)
                 END IF
            END IF
            pLon = vector(2)
            pLat = vector(3)
            IF (ABS(pLat) > 90.01) THEN
			  write(ErrorMsg,'(A,I0)') "ABS(latitude) > 90 AT NODE ",index
			  call FatalError(ErrorMsg,ThID)
            END IF
            IF (ABS(pLat) > 89.99D0) THEN
			  write(ErrorMsg,'(A,I0,A/,A/,A)') "NODE ",index," LIES ON A POLE. ",&
			    & "THIS IS A SINGULAR POINT OF THE SPHERICAL COORDINATE SYSTEM. ", &
			    & "MOVE THIS NODE, AT LEAST SLIGHTLY."
			  call FatalError(ErrorMsg,ThID)
            END IF
            xi = (90.0D0 - pLat) * 0.0174532925199433D0
            yi = pLon * 0.0174532925199433D0
            elevi = vector(4)
            qi = vector(5)
            IF (index <= nRealN) THEN
                 i = index
            ELSE
                 i = nRealN + index - n1000
            END IF
            checkN(i) = .TRUE.
            xNode(i) = xi
            yNode(i) = yi
            elev(i) = elevi
            dQdTdA(i) = qi
            IF (qi < 0.0D0) THEN
			  write(ErrorMsg,'(A)') "NEGATIVE HEAT-FLOW IS NON-PHYSICAL."
			  call FatalError(ErrorMsg,ThID)
            END IF
            IF (.NOT.brief) THEN
                 IF(Verbose) WRITE (iUnitT, 99) INDEX, pLon, pLat, xi, yi, elevi, qi
   99            FORMAT (' ', I10, 2F12.3, 2F11.5, 2ES10.2)
            END IF

  100  CONTINUE

       allOK = .TRUE.
       DO 101 i = 1, numNod
            allOK = allOK.AND.checkN(i)
  101  CONTINUE
       IF (.NOT.allOK) THEN
		 ErrorMsg = "THE FOLLOWING NODES WERE NEVER READ:"
		    j = 1
			allocate(ErrorArray(size(checkN)-count(checkN)))
            do i = 1, numNod
              if(i <= nRealN) THEN
                index = i
              else
                index = n1000 + i - nRealN
              end if
              if(.NOT.checkN(i)) then
			    ErrorArray(j) = index
				j = j+1
			  end if
            end do
		 call FatalError(ErrorMsg,ThID,ErrArr=ErrorArray)
       END IF

!  Read triangular elements:

       READ (iunit7, * ) numEl
       IF (numEl > mxEl) THEN
		 write(ErrorMsg,'(A,I0,A)') "INCREASE PARAMETER maxEl TO BE AT LEAST EQUAL TO THE NUMBER OF ELEMENTS (",numEl,") AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
       END IF
       DO 109 k = 1, numEl
            checkE(k) = .FALSE.
  109  CONTINUE
       IF (.NOT.brief .AND. Verbose) WRITE (iUnitT, 110) numEl
  110       FORMAT(/' There are ',I6,' triangular continuum elements.'/ &
     &      ' (Node numbers for each are given at corners, counter', &
     &      'clockwise'/ / &
     &     ' element        C1        C2        C3')

       DO 200 k = 1, numEl

            !(Elements need not be input in order, but must all be present.)
            READ (iUnit7, "(A)") longer_line
            CALL Extract_LRi (longer_line, &     ! input
                            & LRi, shorter_line) ! output
            READ (shorter_line, * ) i, (nodes(j, i), j = 1, 3)
            IF ((i < 1).OR.(i > numEl)) THEN
			  write(ErrorMsg,'(A,I0)') "ILLEGAL ELEMENT NUMBER: ", i
			  call FatalError(ErrorMsg,ThID)
            END IF
            checkE(i) = .TRUE.
            continuum_LRi(i) = LRi
            IF (.NOT.brief) THEN
                 IF (LRi == 0) THEN
                      IF(Verbose) WRITE (iUnitT, 120) i, (nodes(j, i), j = 1, 3)
  120                 FORMAT (' ', I8, ':', 3I8)
                 ELSE
                      IF(Verbose) WRITE (iUnitT, "(' ', I8, ':', 3I8, ' LR', I8)") i, (nodes(j, i), j = 1, 3), LRi
                 END IF
            END IF
            DO 130 j = 1, 3
                 n = nodes(j, i)
                 IF (n > nRealN) n = nRealN + (n - n1000)
                 IF ((n <= 0).OR.(n > numNod)) THEN
				   write(ErrorMsg,'(A,I0,A)') "NODE NUMBER ", nodes(j, i)," IS ILLEGAL."
				   call FatalError(ErrorMsg,ThID)
                 END IF
                 nodes(j, i) = n
  130       CONTINUE

  200  CONTINUE

       allOK = .TRUE.
       DO 201 i = 1, numEl
            allOK = allOK.AND.checkE(i)
  201  CONTINUE
       IF (.NOT.allOK) THEN
		 write(ErrorMsg,'(A)') "THE FOLLOWING ELEMENTS WERE NEVER READ:"
		 allocate(ErrorArray(size(checkN)-count(checkN)))
		 j = 1
		 do i = 1, numEl
           if(.NOT.checkE(i)) then
		     ErrorArray(j) = i
             j = j+1
		   end if
		 end do
		 call FatalError(ErrorMsg,ThID,ErrArr=ErrorArray)
       END IF

!   Read fault elements:

       READ (iunit7, * ) nFl
       IF (nFl > mxFEl) THEN
		 write(ErrorMsg,'(A,I0,A)') "INCREASE PARAMETER maxFEl TO BE AT LEAST EQUAL TO THE NUMBER OF FAULTS (",nFl,") AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
       END IF
       offMax = 0.0D0
       DO 222 i = 1, nFl
            checkF(i) = .FALSE.
  222  CONTINUE
       IF (.NOT.brief .AND. Verbose) WRITE(iUnitT, 230) nFl
  230  FORMAT(/ /' There are ',I6,' great-circle fault elements.')
       IF ((.NOT.brief).AND.(nFl > 0) .AND. Verbose) WRITE(iUnitT, 231)
  231  FORMAT (/' (The 4 node numbers defining each element must be', &
     &      ' in a counterclockwise order:'/ &
     &      '  n1, and n2 are in left-to-right seguence on the', &
     &      ' near side,'/ &
     &      '  then n3 is opposite n2, and n4 is opposite n1.'/, &
     &      '  (Fault dips are given at n1, n2,', &
     &       ' in degrees from horizontal;'/ &
     &      '  positive dips are toward n1 and n2, respectively, '/ &
     &      '  while negative dips are toward n4 and n3.)'/ &
     &      '  Offset is the total past slip of the fault.'/ / &
     &      ' Element   n1   n2   n3   n4   dip1  dip2', &
     &      '    offset'/)
  240       FORMAT (' ', I8, ':', 4I8, 1X, 2F6.1, 1X, F9.0)

       DO 300 k = 1, nFl

            off = 0.0D0
            READ (iUnit7, "(A)") longer_line
            CALL Extract_LRi (longer_line, &     ! input
                            & LRi, shorter_line) ! output
            READ(shorter_line, * ) i, (nodeF(j, k), j = 1, 4), (dips(l), l = 1, 2), off
            IF ((i < 1).OR.(i > nFl)) THEN
			  write(ErrorMsg,'(A,I0)') "ILLEGAL FAULT NUMBER: ", i
			  call FatalError(ErrorMsg,ThID)
            END IF
            checkF(i) = .TRUE.
            fault_LRi(i) = LRi
            IF (.NOT.brief) THEN
                 IF (LRi == 0) THEN
                      IF(Verbose) WRITE (iUnitT, 240) i, (nodeF(j, i), j = 1, 4), (dips(l), l = 1, 2), off
                 ELSE
                      IF(Verbose) WRITE (iUnitT, "(' ', I8, ':', 4I8, 1X, 2F6.1, 1X, F9.0, ' LR', I8)") i, (nodeF(j, i), j = 1, 4), (dips(l), l = 1, 2), off, LRi
                 END IF
            END IF
            DO 250 j = 1, 4
                 n = nodeF(j, i)
                 IF (n > nRealN) n = nRealN + (n - n1000)
                 IF ((n <= 0).OR.(n > numNod)) THEN
				   write(ErrorMsg,'(A,I0,A,I0)') "ILLEGAL NODE NUMBER ", nodeF(j, i)," IN FAULT", i
				   call FatalError(ErrorMsg,ThID)
                 END IF
                 nodeF(j, i) = n
  250       CONTINUE
            DO 260 l = 1, 2
                 IF (ABS(dips(l)) > 90.0D0) THEN
		           write(ErrorMsg,'(A,F10.4,A/,A/,A)') "ILLEGAL DIP OF ", dips(l)," ; SHOULD BE IN RANGE OF -90. TO +90. DEGREES.",&
				   & " (NOTE: ALL DIPS ARE IN DEGREES FROM THE HORIZONAL; A + pRefIX (OR NONE) INDICATES A DIP TOWARD THE n1-n2 SIDE;", &
				   & " A - pRefIX INDICATES A DIP TOWARD THE n4-n3 SIDE.)"
		           call FatalError(ErrorMsg,ThID)
                 END IF
                 IF (dips(l) < 0.0D0) dips(l) = 180.0D0 + dips(l)
                 fDip(l, i) = dips(l) * 0.0174532925199433D0
  260       CONTINUE
            IF (off < 0.0D0) THEN
		      write(ErrorMsg,'(A,1P,E10.2,A,I0,A)') "ILLEGAL FAULT OFFSET OF ",off," FOR FAULT ELEMENT ", k,". OFFSETS MAY NOT BE NEGATIVE"
		      call FatalError(ErrorMsg,ThID)
            END IF
            offset(i) = off
            offMax = MAX(offMax, off)

  300  CONTINUE

       allOK = .TRUE.
       DO 301 i = 1, nFl
            allOK = allOK.AND.checkF(i)
  301  CONTINUE
       IF (.NOT.allOK) THEN
	     ErrorMsg = "THE FOLLOWING FAULTS WERE NEVER READ:"
		 allocate(ErrorArray(size(checkN)-count(checkN)))
		 j = 1
		 do i = 1, nFl
           if(.not.checkF(i)) then
		     ErrorArray(j) = i
		     j = j+1
		   end if
		 end do
		 call FatalError(ErrorMsg,ThID,ErrArr=ErrorArray)
       ELSE
            IF (offMax > 0.0D0) THEN
                 IF(Verbose) WRITE (iUnitT, 400) offMax
  400            FORMAT (/' Greatest fault offset read was ',1P,D10.2)
            ELSE
                 IF(Verbose) WRITE (iUnitT, 401)
  401            FORMAT (/' Since fault offsets are all zero,', &
     &               ' input parameter Byerly will have no effect.')
            END IF
       END IF
       IF (Verbose) WRITE (iUnitT, 999)
  999  FORMAT (' -------------------------------------------------------------------------------')
       END SUBROUTINE GetNet

       SUBROUTINE PutNet (iUnitO, &                               ! INTENT(IN)
     &                    brief, &                                ! INTENT(IN)
     &                    continuum_LRi, &                        ! INTENT(IN)
     &                    dQdTdA, elev, &                         ! INTENT(IN)
     &                    fault_LRi, fDip, &                      ! INTENT(IN)
     &                    mxEl, mxFEl, mxNode, n1000, &           ! INTENT(IN)
     &                    nFakeN, nFl, nodeF, nodes, &            ! INTENT(IN)
     &                    nRealN, numEl, numNod, offset, &        ! INTENT(IN)
     &                    title1, tLNode, xNode, yNode, zMNode, & ! INTENT(IN)
     &                    chemical_delta_rho_list, &              ! INTENT(IN)
     &                    cooling_curvature_list)                 ! INTENT(IN)

!   Writes finite-element grid (.FEG file) to unit "iUnitO".

       IMPLICIT NONE
       !Arguments:
       INTEGER,      INTENT(IN) :: iUnitO
       LOGICAL,      INTENT(IN) :: brief
       INTEGER,      INTENT(IN) :: continuum_LRi
       REAL*8,       INTENT(IN) :: dQdTdA, elev
       INTEGER,      INTENT(IN) :: fault_LRi
       REAL*8,       INTENT(IN) :: fDip
       INTEGER,      INTENT(IN) :: mxEl, mxFEl, mxNode, n1000, &
                                 & nFakeN, nFl, nodeF, nodes,  &
                                 & nRealN, numEl, numNod
       REAL*8,       INTENT(IN) :: offset
       CHARACTER*100, INTENT(IN) :: title1
       REAL*8,       INTENT(IN) :: tLNode, xNode, yNode, zMNode, &
                                 & chemical_delta_rho_list, &
                                 & cooling_curvature_list
       DIMENSION chemical_delta_rho_list(mxNode), &
     &           continuum_LRi(mxEl), &
     &           cooling_curvature_list(mxNode), &
     &           dQdTdA(mxNode), elev(mxNode), &
     &           fault_LRi(mxFEl), &
     &           fDip(2, mxFEl), &
     &           nodeF(4, mxFEl), nodes(3, mxEl), &
     &           offset(mxFEl), tLNode(mxNode), &
     &           xNode(mxNode), yNode(mxNode), zMNode(mxNode)

       !Internal variables:
       INTEGER :: i, iPrint, k, LRi
       INTEGER, DIMENSION(4) :: nP
       REAL*8 :: pLat, pLon
       REAL*8, DIMENSION(2) :: dips

       IF(Verbose) WRITE(iUnitVerb,"(/ /' READY TO CREATE OUTPUT .FEG FILE ON UNIT ', I3)") iUnitO
       WRITE (iUnitO, 1) title1
    1  FORMAT (A80)

       WRITE (iUnitO, 2) numNod, nRealN, nFakeN, n1000, brief
    2  FORMAT(4I8, L8, ' (numNod, nRealN, nFakeN, N1000, BRIEF)')

       DO 100 i = 1, numNod

            IF (i <= nRealN) THEN
                 iPrint = i
            ELSE
                 iPrint = n1000 + (i - nRealN)
            END IF
            pLat = 90.0D0 - xNode(i) * 57.2957795130823D0
            pLon = yNode(i) * 57.2957795130823D0
            WRITE (iUnitO, 91) i, pLon, pLat, elev(i), dQdTdA(i), zMNode(i), &
     &                         tLNode(i), chemical_delta_rho_list(i), &
     &                         cooling_curvature_list(i)
   91       FORMAT (I8, 2F11.5, 6ES10.2)

  100  CONTINUE

       WRITE (iUnitO, 110) numEl
  110  FORMAT (I10,' (numEl = NUMBER OF TRIANGULAR CONTINUUM ELEMENTS)')

       DO 200 i = 1, numEl

            DO 150 k = 1, 3
                 IF (nodes(k, i) <= nRealN) THEN
                      nP(k) = nodes(k, i)
                 ELSE
                      nP(k) = n1000 + (nodes(k, i) - nRealN)
                 END IF
  150       CONTINUE
            LRi = continuum_LRi(i)
            IF (LRi == 0) THEN
                 WRITE (iUnitO, 160) i, (nP(k), k = 1, 3)
  160            FORMAT (I8, 3I8)
            ELSE
                 WRITE (iUnitO, "(I8, 3I8, ' LR', I8)") i, (nP(k), k = 1, 3), LRi
            END IF

  200  CONTINUE

       WRITE (iUnitO, 210) nFl
  210  FORMAT (I10,' (nFl =  NUMBER OF CURVILINEAR FAULT ELEMENTS)')

       DO 300 i = 1, nFl

            DO 220 k = 1, 4
                 IF (nodeF(k, i) <= nRealN) THEN
                      nP(k) = nodeF(k, i)
                 ELSE
                      nP(k) = n1000 + (nodeF(k, i) - nRealN)
                 END IF
  220       CONTINUE
            DO 230 k = 1, 2
                 dips(k) = fDip(k, i)
                 dips(k) = dips(k) * 57.2957795130823D0
                 IF (dips(k) > 90.01D0) dips(k) = dips(k) - 180.0D0
  230       CONTINUE
            LRi = fault_LRi(i)
            IF (LRi == 0) THEN
                 WRITE (iUnitO, 250) i, (nP(k), k = 1, 4), (dips(k), k = 1, 2), offset(i)
  250            FORMAT (I8, 4I8, 2F6.1, ES10.2)
            ELSE
                 WRITE (iUnitO, "(I8, 4I8, 2F6.1, ES10.2, ' LR', I8)") i, (nP(k), k = 1, 4), (dips(k), k = 1, 2), offset(i), LRi
            END IF

  300  CONTINUE
       END SUBROUTINE PutNet

      SUBROUTINE Extract_LRi (longer_line, &     ! input
                            & LRi, shorter_line) ! output
      ! New routine added for Shells_v5.0+ to support multiple
      !"Lithospheric Rheology" (abbreviated as "LR") integer codes,
      ! in any line of the input .feg file which define an element
      !(either a triangular continuum element, or a
      !          linear fault element).
      ! CHARACTER*80, INTENT(IN) :: longer_line is the whole
      !                             element-definition line from the .feg file.
      ! INTEGER, INTENT(OUT) :: LRi is the rheologic code
      !                        (or 0, if no such code was found).
      ! CHARACTER*80, INTENT(OUT) :: shorter_line has the " LRi" portion removed (if any),
      !                              so it can be interpreted by the same code as in Shells_v4.1-.
      IMPLICIT NONE
      CHARACTER*80, INTENT(IN) :: longer_line
      INTEGER, INTENT(OUT) :: LRi
      CHARACTER*80, INTENT(OUT) :: shorter_line
      CHARACTER*80 :: string
      INTEGER :: longer_length, LR_start_byte
      longer_length = LEN_TRIM(longer_line)
      LR_start_byte = INDEX(longer_line, "LR")
      IF (LR_start_byte > 0) THEN ! the "LR" flag was found
          IF (longer_length > (LR_start_byte + 1)) THEN ! some byte(s) follow the "LR"
              string = longer_line((LR_start_byte + 2):longer_length)
              READ (string, *) LRi
              shorter_line = longer_line(1:(LR_start_byte - 1))
          ELSE ! "LR" is present, but nothing follows it; infer 0.
              LRi = 0
              shorter_line = longer_line(1:(LR_start_byte - 1))
          END IF
      ELSE ! no "LR" flag is present
          LRi = 0
          shorter_line = longer_line
      END IF
      END SUBROUTINE Extract_LRi


end module
