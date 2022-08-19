!*******************************************************************************
! Module containing subroutines used by OrbScore
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

module ScoreSubs

use ShellSetSubs
use DSphere
use SharedVars

implicit none

contains


       SUBROUTINE Adjust (data, eLon, iUnitT, mxAdj, n, nLat, &
     &                    radius, weights, &                   ! INTENT(IN)
     &                    predic, &                            ! INTENT(INOUT)
     &                    eLonP, nLatP, rate)                  ! INTENT(OUT)

!   Compares predicted (horizontal) velocities "predic"
!   (whose first component is the Southward or +Theta velocity,
!   and whose second component is the Eastward or +Phi velocity)
!   with geodetic horizontal velocities "data" (in the same format)
!   which have assigned (dimensionless) "weights" of mean value 1.0,
!   at a set of points i=1, ... ,"n"  <= "mxAdj"
!   (defined by "eLon(i)" and "nLat(i)", in degrees East & North)
!   and finds the rotation correction (to add) which minimizes the
!   RMS of the magnitude of the velocity differences.
!   It then adds this correction to "predic", and reports
!   the pole in user-friendly coordinates "eLonP" and "nLatP"
!   (in degrees East and North)
!   and gives the rotation rate "rate" in radians/second.

       USE DSphere ! in Peter Bird's file DSphere.f90 (just for radians_per_degree, degrees_per_radian)
       USE MKL95_PRECISION ! Intel's Math Kernel Library
       USE MKL95_LAPACK    ! Intel's Math Kernel Library, LAPACK (Linear Analysis Package) portion
       IMPLICIT NONE
       INTEGER,                     INTENT(IN)    :: iUnitT, mxAdj, n
       REAL*8, DIMENSION(2, mxAdj), INTENT(IN)    :: data   ! components (Theta, Phi) = (S, E)
       REAL*8, DIMENSION(mxAdj),    INTENT(IN)    :: eLon, nLat, weights
       REAL*8,                      INTENT(IN)    :: radius
       REAL*8, DIMENSION(2, mxAdj), INTENT(INOUT) :: predic ! components (Theta, Phi) = (S, E)
       REAL*8,                      INTENT(OUT)   :: eLonP, nLatP, rate

       CHARACTER*1 :: uplo
       INTEGER :: i, info, lwork
       INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
       REAL*8 :: DeltaVPhi, DeltaVTheta, dPhi, dTheta, dx, dy, dz, &
               & equat, &
               & pPhi, pTheta, px, py, pz, &
               & rx, ry, rz, &
               & w
       REAL*8, DIMENSION(3) :: DeltaVXyz, dXyz, pXyz, rMeters, uPhi, uTheta, ur
       REAL*8, DIMENSION(3, 3) :: coef
       REAL*8, DIMENSION(3, 1) :: right ! N.B. Array "right" requires 2 subscripts, although only the first column is used.
                                        !      This is an abritrary rule imposed by MKL routine dsysv, which only accepts 2-subscript arrays as inputs.
       REAL*8, DIMENSION(3)    :: solut
       REAL*8, DIMENSION(:), ALLOCATABLE :: work

       IF (n > mxAdj) THEN
	     write(ErrorMsg,'(A,I8,A,I8)') "ERR0R: PARAMETER mxAdj (NOW ",mxAdj,") MUST BE INCREASED TO AT LEAST",n
		 call FatalError(ErrorMsg,ThID)
       END IF
       coef  = 0.0D0 ! all (1:3, 1:3) components
       right = 0.0D0 ! all (1:3, 1:1) components
       DO i = 1, n ! build up coefficient matrix and right-hand-side of linear system by summing over benchmarks:
           CALL DLonLat_2_Uvec(eLon(i), nLat(i), ur)
           CALL DLocal_Theta(ur, uTheta) ! N.B. This will crash if any geodetic benchmark is exactly on the N or S pole.
           CALL DLocal_Phi  (ur, uPhi)   ! N.B. This will crash if any geodetic benchmark is exactly on the N or S pole.
           rx = radius * ur(1)
           ry = radius * ur(2)
           rz = radius * ur(3)
           dTheta = data(1, i)
           dPhi   = data(2, i)
           dXyz(1:3) = dTheta * uTheta(1:3) + dPhi * uPhi(1:3)
           dx = dXyz(1)
           dy = dXyz(2)
           dz = dXyz(3)
           pTheta = predic(1, i)
           pPhi   = predic(2, i)
           pXyz(1:3) = pTheta * uTheta(1:3) + pPhi * uPhi(1:3)
           px = pXyz(1)
           py = pXyz(2)
           pz = pXyz(3)
           w = weights(i)
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           coef(1, 1) = coef(1, 1)   + w * ( 2.0D0 * rz**2 + 2.0D0 * ry**2) ! coefficient of Omega_x in equation #1 (minimization w.r.t. Omega_x)
           coef(1, 2) = coef(1, 2)   + w * (-2.0D0 * ry * rx)               ! coefficient of Omega_y in equation #1 (minimization w.r.t. Omega_x)
           coef(1, 3) = coef(1, 3)   + w * (-2.0D0 * rz * rx)               ! coefficient of Omega_z in equation #1 (minimization w.r.t. Omega_x)
           right(1, 1) = right(1, 1) - w * (-2.0D0 * rz * (py - dy) + 2.0D0 * ry * (pz - dz))   ! right-hand side of equation #1 (minimization w.r.t. Omega_x)
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           coef(2, 1) = coef(2, 1)   + w * (-2.0D0 * rx * ry)               ! coefficient of Omega_x in equation #2 (minimization w.r.t. Omega_y)
           coef(2, 2) = coef(2, 2)   + w * ( 2.0D0 * rz**2 + 2.0D0 * rx**2) ! coefficient of Omega_y in equation #2 (minimization w.r.t. Omega_y)
           coef(2, 3) = coef(2, 3)   + w * (-2.0D0 * rz * ry)               ! coefficient of Omega_z in equation #2 (minimization w.r.t. Omega_y)
           right(2, 1) = right(2, 1) - w * ( 2.0D0 * rz * (px - dx) - 2.0D0 * rx * (pz - dz))   ! right-hand side of equation #2 (minimization w.r.t. Omega_y)
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           coef(3, 1) = coef(3, 1)   + w * (-2.0D0 * rx * rz)               ! coefficient of Omega_x in equation #3 (minimization w.r.t. Omega_z)
           coef(3, 2) = coef(3, 2)   + w * (-2.0D0 * ry * rz)               ! coefficient of Omega_y in equation #3 (minimization w.r.t. Omega_z)
           coef(3, 3) = coef(3, 3)   + w * ( 2.0D0 * ry**2 + 2.0D0 * rx**2) ! coefficient of Omega_z in equation #3 (minimization w.r.t. Omega_z)
           right(3, 1) = right(3, 1) - w * (-2.0D0 * ry * (px - dx) + 2.0D0 * rx * (py - dy))   ! right-hand side of equation #3 (minimization w.r.t. Omega_z)
          !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       END DO

       !Option to WRITE-out this linear system before attempting a solution:
       IF(Verbose) WRITE(iUnitT,*)
       IF(Verbose) WRITE(iUnitT,"(' Linear system to be solved in subprogram Adjust:')")
       DO i = 1, 3
           IF (i == 2) THEN
               IF(Verbose) WRITE(iUnitT,"(/ ' ', 3ES11.3, ' X Omega_', I1, ' = ', ES11.3)") (coef(i, j), j = 1, 3), i, right(i, 1)
           ELSE
               IF(Verbose) WRITE(iUnitT,"(/ ' ', 3ES11.3, '   Omega_', I1, '   ', ES11.3)") (coef(i, j), j = 1, 3), i, right(i, 1)
           END IF
       END DO
       IF(Verbose) WRITE(iUnitT,*)

       !Using a routine from Intel's Math Kernel Library, Linear Analysis Package (LAPACK) portion:
       !Using dsysv for this SYMMETRIC INDEFINITE linear system:
       uplo = 'U' ! problem is stated in the Upper Triangle (plus diagonal) of the coefficient matrix.
       ALLOCATE ( ipiv(3) )
      !lwork = -1 ! First time: "lwork = -1" signals a workspace-size-query to dsysv; answer to be placed in work(1). Answer was "3.0000000D0".
       lwork = 3  ! Based on test described above.
       ALLOCATE ( work(lwork) )
      !call dsysv(uplo,    n, nrhs,    a,  lda, ipiv,     b,  ldb, work, lwork, info ) ! <=== This is the sample Fortran77 CALL in the manual.
       CALL dsysv(uplo,    3,    1, coef,    3, ipiv, right,    3, work, lwork, info)  ! N.B. I'm using the Fortran77 CALL because the F95 CALL is buggy.
       IF (info /= 0) THEN
	     write(ErrorMsg,'(A,I12)') "Error when Adjust calls dsysv of MKL_LAPACK: info = ", info
		 call FatalError(ErrorMsg,ThID)
       END IF
       DEALLOCATE ( work )
       DEALLOCATE ( ipiv )
       solut(1:3) = right(1:3, 1) ! unpack the solution.  Note that only the first column has been used.

       !Characterize the Euler pole for the rotation that will be added:
       rate = SQRT(solut(1)**2 + solut(2)**2 + solut(3)**2)
       IF (rate > 0.0D0) THEN
           equat = SQRT(solut(1)**2 + solut(2)**2)
           eLonP = degrees_per_radian * ATAN2(solut(2), solut(1))
           nLatP = degrees_per_radian * ATAN2(solut(3), equat)
           !Add this rotation to all velocities in predic:
           DO i = 1, n
                CALL DLonLat_2_Uvec(eLon(i), nLat(i), ur)
                CALL DLocal_Theta(ur, uTheta) ! N.B. This will crash if any geodetic benchmark is exactly on the N or S pole.
                CALL DLocal_Phi  (ur, uPhi)   ! N.B. This will crash if any geodetic benchmark is exactly on the N or S pole.
                rMeters(1:3) = radius * ur(1:3)
                CALL DCross(solut, rMeters, DeltaVXyz)
                DeltaVTheta = DDot(DeltaVXyz, uTheta)
                DeltaVPhi   = DDot(DeltaVXyz, uPhi)
                predic(1, i) = predic(1, i) + DeltaVTheta
                predic(2, i) = predic(2, i) + DeltaVPhi
           END DO
       ELSE
           eLonP = 0.0D0
           nLatP = 0.0D0
       END IF
       END SUBROUTINE Adjust


       REAL*8 FUNCTION ATanPV (y, x)

!   Returns principal value of inverse tangent of (y, x).
!   N.B. Order of these two arguments is NOT an error;
!        this order corresponds to ("motivator", "restrainer").

       IMPLICIT NONE
       REAL*8, INTENT(IN) :: y, x
       REAL*8 :: xx, yy
       IF (x >= 0.0D0) THEN
            xx = x
            yy = y
       ELSE
            xx = -x
            yy = -y
       END IF
       IF ((yy == 0.0D0).OR.(xx == 0.0D0)) THEN
          IF (yy == 0.0D0) ATanPV = 0.0D0
          IF (xx == 0.0D0) ATanPV = 1.57079632679490D0
       ELSE
          ATanPV = ATAN2(yy, xx)
       END IF
       END FUNCTION ATanPV


       REAL*8 FUNCTION Chord (angle1, s, angle2)

!   Returns an angle obtained by interpolation between angle1
!   and angle2.  The interpolation method is NOT sensitive to
!   possible cycle shifts (of 2*n*Pi) between angle1 and angle2.

!   Unit vectors are constructed for angle1 and angle2, and a
!   linear chord is drawn between their tips.

!   S is the internal coordinate along the chord;
!   it is dimensionless, with value 0.0D0 at angle1 and 1.0D0 at
!   angle2.  (The user may input s values outside this range
!   to get results outside the (smaller) angle between angle1 and
!   angle2, if desired.)  The angle returned is that from the
!   origin to this chord point.

!   This algorithm should work equally well for angles measured
!   either clockwise or counterslockwise from any reference, as
!   long as the usage is consistent.

!   Both the input angles and the result "chord" are in radians.

       IMPLICIT NONE
       REAL*8, INTENT(IN) :: angle1, s, angle2
       REAL*8 :: c1, c2
       REAL*8, DIMENSION(2) :: uvec1, uvec2, vecs
       uvec1(1) = COS(angle1)
       uvec1(2) = SIN(angle1)
       uvec2(1) = COS(angle2)
       uvec2(2) = SIN(angle2)
       vecs(1) = (1.0D0 - s) * uvec1(1) + s * uvec2(1)
       vecs(2) = (1.0D0 - s) * uvec1(2) + s * uvec2(2)
       c1 = vecs(1)
       c2 = vecs(2)
       Chord = ATan2F(c2, c1)
       END FUNCTION Chord


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


       SUBROUTINE Deriv3 (iUnitT, mxEl, mxNode, &
     &                    nodes, numEl, &
     &                    radius, xNode, yNode, & ! INTENT(IN)
     &                    area, detJ, &           ! INTENT(OUT)
     &                    dxs, dys, dxsp, dysp, fpsfer, sita)

!   Sets up 6 vector nodal functions (fpsfer) of each spherical
!     triangle finite element, at each of its 7 integration points.
!   Calculates dxs and dys, the theta-derivitive and phi-derivitive
!     of each of these 6 vector nodel functions.
!   Also computes area, the areas of the plane triangles.
!   Also computes detJ, the local ratio of areas on the sphere
!     to areas on the plane triangles.

      USE DSphere ! in Peter Bird's file DSphere.f90
      IMPLICIT NONE
      INTEGER,                     INTENT(IN) :: iUnitT, mxEl, mxNode, numEl
      INTEGER, DIMENSION(3, mxEl), INTENT(IN) :: nodes
      REAL*8,                      INTENT(IN) :: radius
      REAL*8, DIMENSION(mxNode),   INTENT(IN) :: xNode, yNode
      REAL*8, DIMENSION(mxEl),             INTENT(OUT) :: area
      REAL*8, DIMENSION(7, mxEl),          INTENT(OUT) :: detJ
      REAL*8, DIMENSION(2, 2, 3, 7, mxEl), INTENT(OUT) :: dxs, dys
      REAL*8, DIMENSION(3, 7, mxEl),       INTENT(OUT) :: dxsp, dysp
      REAL*8, DIMENSION(2, 2, 3, 7, mxEl), INTENT(OUT) :: fpsfer
      REAL*8, DIMENSION(7, mxEl),          INTENT(OUT) :: sita

      INTEGER :: i, j, m
      REAL*8 :: a, areap, b, &
              & c, cka, cosm, csccse, cscs, cscsne, cse, cssn, &
              & dd, dd1, dd2, dd3, ddpn, dpdc, dpde, &
              & fff(3), &
              & phaij, pfq, phi(3), pnx, pny, pnz, pp, &
              & rn, rr1, rr2, rr3, &
              & sitaj, sitami, skkc(3), skke(3), snc, sne, sncsne, snccse, &
              & theta(3), tx, ty, &
              & x21, x31, xa, xb, xc, xyzp, &
              & y21, y31, ya, yb, yc, &
              & z21, z31, za, zb, zc

      REAL*8, DIMENSION(3, 7) :: points
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
         areap = SQRT(a * a + b * b + c * c)
         area(i) = radius * radius * (0.5 * areap)
         pnx = a / areap
         pny = b / areap
         pnz = c / areap
         dd1 = SIN(theta(1)) * COS(phi(1)) * pnx
         dd2 = SIN(theta(1)) * SIN(phi(1)) * pny
         dd3 = COS(theta(1)) * pnz
         dd = dd1 + dd2 + dd3

! This part is to test if Kong's method and Bird's method are same for
! calculating derivative:

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
            IF ((sita(m, i) <= 0.0D0) .OR. (sita(m, i) >= 3.14159265358979D0)) THEN
               sitami = sita(m, i) * degrees_per_radian
	           write(ErrorMsg,'(A,I5,A,I5,A,E14.4)') "COLATITUDE OF INTEGRATION POINT ",m," OF ELEMENT ",i," IS OUT RANGE ",sitami
		       call FatalError(ErrorMsg,ThID)
            END IF
            DO 500 j = 1, 3
               dxsp(j, m, i) = dpdc * fff(j) + pp * skkc(j)
               dysp(j, m, i) = dpde * fff(j) + pp * skke(j)
               cscs = COS(theta(j)) * COS(phi(j))
               cssn = COS(theta(j)) * SIN(phi(j))
               snc = SIN(theta(j))
               sne = SIN(phi(j))
               cse = COS(phi(j))
               fpsfer(1, 1, j, m, i) = cscs * csccse + cssn * cscsne + snc * SIN(sitaj)
               fpsfer(2, 1, j, m, i) = -sne * csccse + cse * cscsne
               fpsfer(1, 2, j, m, i) = -cscs * SIN(phaij) + cssn * COS(phaij)
               fpsfer(2, 2, j, m, i) = sne * SIN(phaij) + cse * COS(phaij)
               dxs(1, 1, j, m, i) = (-cscs * snccse - cssn * sncsne &
     &                           + snc * COS(sitaj)) * fff(j) &
     &                          + fpsfer(1, 1, j, m, i) * skkc(j)
               dxs(2, 1, j, m, i) = (sne * snccse - cse * sncsne) * fff(j) &
     &                          + fpsfer(2, 1, j, m, i) * skkc(j)
               dys(1, 1, j, m, i) = (-cscs * cscsne + cssn * csccse) * fff(j) &
     &                          + fpsfer(1, 1, j, m, i) * skke(j)
               dys(2, 1, j, m, i) = (sne * cscsne + cse * csccse) * fff(j) &
     &                          + fpsfer(2, 1, j, m, i) * skke(j)
               dxs(1, 2, j, m, i) = fpsfer(1, 2, j, m, i) * skkc(j)
               dxs(2, 2, j, m, i) = fpsfer(2, 2, j, m, i) * skkc(j)
               dys(1, 2, j, m, i) = (-cscs * COS(phaij) - cssn * SIN(phaij)) &
     &                          * fff(j) &
     &                          + fpsfer(1, 2, j, m, i) * skke(j)
               dys(2, 2, j, m, i) = (sne * COS(phaij) - cse * SIN(phaij)) &
     &                          * fff(j) &
     &                          + fpsfer(2, 2, j, m, i) * skke(j)
               fpsfer(1, 1, j, m, i) = fpsfer(1, 1, j, m, i) * fff(j)
               fpsfer(2, 1, j, m, i) = fpsfer(2, 1, j, m, i) * fff(j)
               fpsfer(1, 2, j, m, i) = fpsfer(1, 2, j, m, i) * fff(j)
               fpsfer(2, 2, j, m, i) = fpsfer(2, 2, j, m, i) * fff(j)
  500       CONTINUE
            pfq = fff(1) + fff(2) + fff(3)
            detJ(m, i) = rn**3 / (dd * dd)
  800    CONTINUE
  900 CONTINUE
      END SUBROUTINE Deriv3


       SUBROUTINE EDot (dXS, dYS, &
     &                  fPSfer, mxEl, &
     &                  mxNode, nodes, numEl, radius, sita, v, & ! INTENT(IN)
     &                  eRate)                                   ! INTENT(OUT)

!   Compute strain rate components EDot_xx, EDot_yy, and
!   EDot_xy (tensor form; equal to
!            (1/2) * ((dVx/dY)+(dVy/dX))
!   at the integration points of triangular continuum elements.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: dXS, dYS, fPSfer                                                 ! input
       INTEGER, INTENT(IN) :: mxEl, mxNode, nodes, numEl                                      ! input
       REAL*8, INTENT(IN) :: radius, sita                                                     ! input
       DOUBLE PRECISION, INTENT(IN) :: v                                                      ! input
       REAL*8, INTENT(OUT) :: eRate                                                           ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION points
       COMMON / S1S2S3 / points
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER i, j, m, node
       REAL*8 dy11, dy21, dy12, dy22, exx, exy, eyy, fp11, fp21, fp12, fp22, vx, vy
       DIMENSION dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
     &           eRate(3, 7, mxEl), &
     &           fPSfer(2, 2, 3, 7, mxEl), &
     &           nodes(3, mxEl), points(3, 7), &
     &           sita(7, mxEl), v(2, mxNode)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       DO 1000 m = 1, 7
            DO 900 i = 1, numEl
                 exx = 0.0D0
                 eyy = 0.0D0
                 exy = 0.0D0
                 DO 800 j = 1, 3
                      node = nodes(j, i)
                      vx = v(1, node)
                      vy = v(2, node)
                      dy11 = dYS(1, 1, j, m, i) / SIN(sita(m, i))
                      dy21 = dYS(2, 1, j, m, i) / SIN(sita(m, i))
                      dy12 = dYS(1, 2, j, m, i) / SIN(sita(m, i))
                      dy22 = dYS(2, 2, j, m, i) / SIN(sita(m, i))
                      fp11 = fPSfer(1, 1, j, m, i) / TAN(sita(m, i))
                      fp21 = fPSfer(2, 1, j, m, i) / TAN(sita(m, i))
                      fp12 = fPSfer(1, 2, j, m, i) / TAN(sita(m, i))
                      fp22 = fPSfer(2, 2, j, m, i) / TAN(sita(m, i))
                      exx = exx + vx * dXS(1, 1, j, m, i) + vy * dXS(2, 1, j, m, i)
                      eyy = eyy + vx * dy12 + vy * dy22 + vx * fp11 + vy * fp21
                      exy = exy + vx * dy11 + vy * dy21 &
     &                          + vx * dXS(1, 2, j, m, i) + vy * dXS(2, 2, j, m, i) &
     &                          - vx * fp12 - vy * fp22
  800            CONTINUE
                 eRate(1, m, i) = exx / radius
                 eRate(2, m, i) = eyy / radius
                 eRate(3, m, i) = 0.50D0 * exy / radius
  900       CONTINUE
 1000  CONTINUE
       RETURN
       END SUBROUTINE EDot


      SUBROUTINE Extract_LRi (longer_line, &     ! input
                            & LRi, shorter_line) ! output
      ! New routine added for Shells_v5.0+ to support multiple
      !"Lithospheric Rheology" (abbreviated as "LR") integer codes,
      ! in any line of the input .feg file which defines an element
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


      SUBROUTINE FAngls (phi, theta, & ! INTENT(IN)
     &                   fAngle)       ! INTENT(OUT)

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

      RETURN
      END SUBROUTINE FAngls


       REAL*8 FUNCTION FltLen (phi1, phi2, radius, theta1, theta2) ! INTENT(IN)

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
       RETURN
       END FUNCTION FltLen


      SUBROUTINE FindIt (mxEl, mxNode, nodes, numEl, xNode, yNode, &
     &                   theta, phi, &         ! INTENT(IN)
     &                   iEle, s1, s2, s3)     ! INTENT(OUT)

!   Determine which element (iEle) contains the surface point
!   (theta, phi), and also its internal triangular coordinates
!   (s1, s2, s3) in the plate triangle; s1 + s2 + s3 = 1.0D0
!   If the point is not in any element, iEle = 0 is returned.

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: mxEl, mxNode, numEl
       INTEGER, DIMENSION(3, mxEl), INTENT(IN) :: nodes
       REAL*8, DIMENSION(mxNode), INTENT(IN) :: xNode, yNode
       REAL*8 :: theta, phi

       INTEGER, INTENT(OUT) :: iEle
       REAL*8,  INTENT(OUT) :: s1, s2, s3

       INTEGER :: i, k
       REAL*8 :: dot, dotp, eps, growth, &
               & perp(3), pole(3), rA(3), rB(3), rn(3, 3), rP(3), rS(3), &
               & s1num, s1den, s2num, s2den
       DATA eps / 1.0D-4 /

!   Cartesian (x, y, z) radius vector to point on sphere:
       rS(1) = SIN(theta) * COS(phi)
       rS(2) = SIN(theta) * SIN(phi)
       rS(3) = COS(theta)
       DO 100 i = 1, numEl
            DO 10 k = 1, 3
!               Find Cartesian radius vector for each node:
                 rn(1, k) = SIN(xNode(nodes(k, i))) * COS(yNode(nodes(k, i)))
                 rn(2, k) = SIN(xNode(nodes(k, i))) * SIN(yNode(nodes(k, i)))
                 rn(3, k) = COS(xNode(nodes(k, i)))
   10       CONTINUE
!          Find normal vector of plane element:
            rA(1) = rn(1, 2) - rn(1, 1)
            rA(2) = rn(2, 2) - rn(2, 1)
            rA(3) = rn(3, 2) - rn(3, 1)
            rB(1) = rn(1, 3) - rn(1, 2)
            rB(2) = rn(2, 3) - rn(2, 2)
            rB(3) = rn(3, 3) - rn(3, 2)
            CALL DCross (rA, rB, perp)
            CALL Unit (perp)
!          Shorten/lengthen radius vector so it points to plane element:
            dot = rS(1) * perp(1) + rS(2) * perp(2) + rS(3) * perp(3)
            IF (dot <= 0.0D0) GO TO 100
            dotp = rn(1, 1) * perp(1) + rn(2, 1) * perp(2) + rn(3, 1) * perp(3)
            growth = dotp / dot
            rP(1) = rS(1) * growth
            rP(2) = rS(2) * growth
            rP(3) = rS(3) * growth
!          Compute and test s1, then s2, then s3...
            CALL DCross (rn(1:3, 2), rn(1:3, 3), pole)
            s1num = rP(1)  * pole(1) + rP(2)  * pole(2) + rP(3)  * pole(3)
            s1den = rn(1, 1) * pole(1) + rn(2, 1) * pole(2) + rn(3, 1) * pole(3)
            IF (s1den /= 0.0D0) THEN
                 s1 = s1num / s1den
            ELSE
                 GO TO 100
            END IF
            IF ((s1 >= -eps).AND.(s1 <= (1.0D0 + eps))) THEN
                 CALL DCross (rn(1:3, 3), rn(1:3, 1), pole)
                 s2num = rP(1)  * pole(1) + rP(2)  * pole(2) + rP(3)  * pole(3)
                 s2den = rn(1, 2) * pole(1) + rn(2, 2) * pole(2) + rn(3, 2) * pole(3)
                 IF (s2den /= 0.0D0) THEN
                      s2 = s2num / s2den
                 ELSE
                      GO TO 100
                 END IF
                 IF ((s2 >= -eps).AND.(s2 <= (1.0D0 + eps))) THEN
                      s3 = 1.00D0 - s1 - s2
                      IF ((s3 >= -eps).AND.(s3 <= (1.0D0 + eps))) THEN
                           iEle = i
                           RETURN
                      END IF
                 END IF
            END IF
  100  CONTINUE
       iEle = 0
       END SUBROUTINE FindIt


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
       RETURN
       END SUBROUTINE FNodal


      SUBROUTINE GetGEO (iUnitB, iUnitT, iUnitY, iUnitC, mxGeo, &
     &                   pltGeo, pltCos, pltInt, iUnitI, &                        ! INTENT(IN)
     &                   geoTag, geoPhi, geoThe, &                ! INTENT(OUT)
     &                   geoVel, geoSig, geoAzi, &
     &                   gpsFMT, numGeo)

!   Reads relative horizontal velocity vectors of benchmarks,
!   determined by geodesy.

!   The input file format expected is:
!   Peter Bird's .gps format of 2002.08.01:
!   3 lines of text/format/headers,
!   then one line per benchmark, giving
!   coordinates, E & N velocity components in mm/a,
!   E & N standard deviations of velocity components in mm/a,
!   correlation between E and N components of velocity,
!   reference frame, and benchmark identification.

!   This subprogram will return:
!      geoPhi = East longitude (radians)
!      geoThe = co-latitude (radians)
!      geoVel = velocity magnitude, m/s
!      geoAzi = azimuth of velocity, in radians counterclockwise from South.
!      geoSig = standard deviation of velocity (assumes circular uncertainty ellipse, and no correlation of errors)
!      geoTag = benchmark name

       USE DSphere ! in Peter Bird's file DSphere.f90
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: iUnitB, iUnitT, iUnitY, iUnitC, mxGeo
       INTEGER,  INTENT(IN) :: iUnitI
       LOGICAL,  INTENT(IN) :: pltGeo, pltCos, pltInt ! Create(?) new .gps-formatted output files for plotting:
                                                      !    *(Shells model error) and/or
                                                      !    *(coseismic part of Shells model velocities) and/or
                                                      !    *(interseismic part of Shells model velocities)
                                                      !      with -FiniteMap-.
       CHARACTER*20, DIMENSION(mxGeo), INTENT(OUT) :: geoTag
       REAL*8,       DIMENSION(mxGeo), INTENT(OUT) :: geoPhi, geoThe
       REAL*8,       DIMENSION(mxGeo), INTENT(OUT) :: geoVel, geoSig, geoAzi
       CHARACTER*132,                  INTENT(OUT) :: gpsFMT
       INTEGER,                        INTENT(OUT) :: numGeo

       CHARACTER*15  :: frame
       CHARACTER*20  :: tag
       CHARACTER*132 :: header
       INTEGER       :: ios
       REAL*8        :: eLonDg, nLatDg, vEMmpa, vNMmpa, vESigm, vNSigm, correl
       REAL*8        :: secPYr
	   CHARACTER(LEN=100) ::  filename,outdir
	   INTEGER :: ModNum=1,iter,ThID


       secPYr = 365.25D0 * 24.0D0 * 60.0D0 * 60.0D0

!   If first READ fails (no file connected), return numGeo = 0

       IF(Verbose) WRITE(iUnitT,1) iUnitB
    1  FORMAT(/' ATTEMPTING TO READ GEODETIC BENCHMARK VELOCITIES (.gps file) FROM UNIT', I3/)
       numGeo = 0

!      Read, but do not use, the first line of text (file mnemonic).
       READ (iUnitB, * , IOSTAT = ios)
       IF (ios /= 0) RETURN
       IF (pltGeo) WRITE(iUnitY, "('ERROR in geodetic velocity predictions')")
       IF (pltCos) WRITE(iUnitC, "('COSEISMIC part of model geodetic velocities')")
       IF (pltInt) WRITE(iUnitI, "('INTERSEISMIC part of model geodetic velocities')")

!      Read and store the FORMAT in line 2
       READ (iUnitB, "(A)", IOSTAT = ios) gpsFMT
       IF (pltGeo) WRITE(iUnitY, "(A)") TRIM(gpsFMT)
       IF (pltCos) WRITE(iUnitC, "(A)") TRIM(gpsFMT)
       IF (pltInt) WRITE(iUnitI, "(A)") TRIM(gpsFMT)

!      Read the column headers in line 3
       READ (iUnitB, "(A)", IOSTAT = ios) header
       IF (ios /= 0) RETURN
       IF (pltGeo) WRITE(iUnitY, "(A)") TRIM(header)
       IF (pltCos) WRITE(iUnitC, "(A)") TRIM(header)
       IF (pltInt) WRITE(iUnitI, "(A)") TRIM(header)

   10  READ (iUnitB, gpsFMT, END = 100) eLonDg, nLatDg, vEMmpa, vNMmpa, &
     &                                  vESigm, vNSigm, correl, frame, tag
            numGeo = numGeo + 1
            geoPhi(numGeo) = eLonDg * radians_per_degree
            geoThe(numGeo) = (90.0D0 - nLatDg) * radians_per_degree
            geoVel(numGeo) = SQRT(vEMmpa**2 + vNMmpa**2) / (1000.0D0 * secPYr)
            geoAzi(numGeo) = ATan2F(vEMmpa, -vNMmpa)
            geoSig(numGeo) = SQRT(vESigm**2 + vNSigm**2) / (1000.0D0 * secPYr)
            IF (geoSig(numGeo) <= 0.0D0) THEN
	     write(ErrorMsg,'(A,I5/,A)') "ERROR: Standard deviation of velocity is <= 0.0 for benchmark ",numGeo, &
	     & "All values must be positive."
		 call FatalError(ErrorMsg,ThID)
            END IF
            geoTag(numGeo) = TRIM(tag)
       IF (numGeo < mxGeo) GO TO 10
	     write(ErrorMsg,'(A,I6,A/,A)') "SUBPROGRAM GetGEO READ ONLY ",mxGeo," BENCHMARKS." , &
	     & "INCREASE PARAMETER maxMOR AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
  100  RETURN
       END SUBROUTINE GetGEO


      SUBROUTINE GetMOR (iUnitM, iUnitT, mxMOR, &  ! INTENT(IN)
     &                   MORTag, MORPhi, MORThe, & ! INTENT(OUT)
     &                   MORVel, MORSig, numMOR)

!   Reads total seafloor spreading rates at mid-ocean rises.

!   Current FORMAT matches file "SPREADIN.NUVEL1".
!      MORTag = ridge-site name
!      MORPhi = East longitude (radians)
!      MORThe = co-latitude (radians)
!      MORVel = velocity magnitude, m/s
!      MORSig = standard deviation of velocity

!   IF FIRST READ FAILS (NO FILE CONNECTED), RETURNS numMOR=0

       USE DSphere ! in Peter Bird's file DSphere.f90
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iUnitM, iUnitT, mxMOR
       CHARACTER*5, DIMENSION(mxMOR), INTENT(OUT) :: MORTag
       REAL*8,      DIMENSION(mxMOR), INTENT(OUT) :: MORPhi, MORThe, MORVel, MORSig
       INTEGER,                       INTENT(OUT) :: numMOR

       CHARACTER*5 :: tag
       INTEGER     :: ios
       REAL*8      :: eLon, nLat, v, sigma
       REAL*8      :: secPYr

       secPYr = 365.25D0 * 24.0D0 * 60.0D0 * 60.0D0

       IF(Verbose) WRITE(iUnitT,1) iUnitM
    1  FORMAT(/' ATTEMPTING TO READ SEAFLOOR SPREADING RATES FROM UNIT', I3/)

       numMOR = 0
       READ (iUnitM, * , IOSTAT = ios)
       IF (ios /= 0) RETURN
       READ (iUnitM, * )
       READ (iUnitM, * )
   10  READ (iUnitM, 11, END = 100) tag, nLat, eLon, v, sigma
   11       FORMAT (A5, F7.2, F8.2, F6.1, F5.1)
            numMOR = numMOR + 1
            MORTag(numMOR) = tag
            MORPhi(numMOR) = eLon * radians_per_degree
            MORThe(numMOR) = (90.0D0 - nLat) * radians_per_degree
            MORVel(numMOR) = v / (1000.0D0 * secPYr)
            MORSig(numMOR) = sigma / (secPYr * 1000.0D0)
       IF (numMOR < mxMOR) GO TO 10
	     write(ErrorMsg,'(A,I6,A/,A)') "SUBPROGRAM GetMOR READ ONLY ",mxMOR," BENCHMARKS." , &
	     & "INCREASE PARAMETER maxMOR AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
  100  RETURN
       END SUBROUTINE GetMOR


       SUBROUTINE GetNet (iUnit7, iUnitT, &           ! input
     &                    mxDOF, mxEl, mxFEl, mxNode, &
     &                    brief, continuum_LRi, cooling_curvature, & ! output
     &                    density_anomaly, &
     &                    dQdTdA, elev, fault_LRi, fDip, &
     &                    nFakeN, nFl, nodeF, nodes, nRealN, &
     &                    numEl, numNod, n1000, offMax, offset, &
     &                    title1, tLNode, xNode, yNode, zMNode, &
     &                    checkE, checkF, checkN)     ! work

!   Read finite element grid from unit iUnit7 (assumed already OPENed).
!   Echoes the important values to unit iUnitT (assumed already OPENed).

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER, INTENT(IN) :: iUnit7, iUnitT, mxDOF, mxEl, mxFEl, mxNode                       ! input
       LOGICAL, INTENT(OUT) :: brief                                                           ! output
       INTEGER, INTENT(OUT) :: continuum_LRi                                                   ! output
       REAL*8, INTENT(OUT) :: cooling_curvature, density_anomaly, dQdTdA, elev                 ! output
       INTEGER, INTENT(OUT) :: fault_LRi                                                       ! output
       REAL*8, INTENT(OUT) :: fDip                                                             ! output
       INTEGER, INTENT(OUT) :: nFakeN, nFl, nodeF, nodes, nRealN, numEl, numNod, n1000         ! output
       REAL*8, INTENT(OUT) :: offMax, offset                                                   ! output
       CHARACTER*100, INTENT(OUT) :: title1                                                     ! output
       REAL*8, INTENT(OUT) :: tLNode, xNode, yNode, zMNode                                     ! output
       LOGICAL checkE, checkF, checkN                                          ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CHARACTER*80 :: longer_line, shorter_line
       LOGICAL allOK
       INTEGER i, index, j, k, l, LRi, n, nrt2
       REAL*8 cooling_curvature_cpm2, density_anomaly_kgpm3, dips, elevi, &
            & off, pLat, pLon, qi, tli, vector, xi, yi, zmi
       DIMENSION checkE(mxEl), checkF(mxFEl), checkN(mxNode), &
     &           continuum_LRi(mxEl), &
     &           cooling_curvature(mxNode), &
     &           density_anomaly(mxNode), &
     &           dQdTdA(mxNode), elev(mxNode), &
     &           fault_LRi(mxFEl), &
     &           fDip(2, mxFEl), nodeF(4, mxFEl), &
     &           nodes(3, mxEl), offset(mxFEl), tLNode(mxNode), &
     &           xNode(mxNode), yNode(mxNode), zMNode(mxNode)
       DIMENSION dips(3), vector(9)

       title1 = ' '
       READ (iUnit7, 2) title1
    2  FORMAT (A80)
       IF(Verbose) WRITE (iUnitT, 3) title1
    3  FORMAT(/' Title of finite-element grid ='/' ',A80)

!   Read number of nodes, plus out-dated parameters that once
!     permitted boundary nodes to be specially numbered as
!    "fake" nodes with numbers from n1000+1 ... n1000+nFakeN.
!     This option is no longer supported by my programs!
!    (Option "brief" suppresses most output.)

       READ (iUnit7, * ) numNod, nRealN, nFakeN, n1000, brief

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
            CALL ReadN (iUnit7, iUnitT, 9, & ! input
     &                  vector)              ! output
            index = vector(1) + 0.5D0
            IF (index > nRealN) THEN
                 IF ((index <= n1000).OR.(index > (n1000 + nFakeN))) THEN
				   write(ErrorMsg,'(A,I0)') "ILLEGAL NODE NUMBER:", index
				   call FatalError(ErrorMsg,ThID)
                 END IF
            END IF
            pLon = vector(2)
            pLat = vector(3)
            IF (ABS(pLat) > 90.01) THEN
			  write(ErrorMsg,'(A,I0)') "ABS(latitude) > 90 AT NODE:", index
			  call FatalError(ErrorMsg,ThID)
            END IF
            IF (ABS(pLat) > 89.99D0) THEN
			  write(ErrorMsg,'(A,I0,A/,A)') "ERR0R: NODE ",index," LIES ON A POLE. ",&
			    &  "THIS IS A SINGULAR POINT OF THE SPHERICAL COORDINATE SYSTEM. MOVE THIS NODE, AT LEAST SLIGHTLY."
			  call FatalError(ErrorMsg,ThID)
            END IF
            xi = (90.0D0 - pLat) * 0.0174532925199433D0
            yi = pLon * 0.0174532925199433D0
            elevi = vector(4)
            qi = vector(5)
            zmi = vector(6)
            tli = vector(7)
            density_anomaly_kgpm3 = vector(8)
            cooling_curvature_cpm2 = vector(9)
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
            IF (zmi < 0.0D0) THEN
			  write(ErrorMsg,'(A)') "NEGATIVE CRUSTAL THICKNESS IS NON-PHYSICAL."
			  call FatalError(ErrorMsg,ThID)
            END IF
            zMNode(i) = zmi
            IF (tli < 0.0D0) THEN
			  write(ErrorMsg,'(A)') "NEGATIVE MANTLE LITHOSPHERE THICKNESS IS NON-PHYSICAL."
			  call FatalError(ErrorMsg,ThID)
            END IF
            tLNode(i) = tli
            IF (.NOT.brief) THEN
                 IF(Verbose) WRITE (iUnitT, 99) INDEX, pLon, pLat, xi, yi, elevi, &
     &                             qi, zmi, tli
   99            FORMAT (' ',I10,0P,2F12.3,2F11.5,1P,3E10.2,E12.2)
            END IF
            density_anomaly(i) = density_anomaly_kgpm3
            cooling_curvature(i) = cooling_curvature_cpm2
  100  CONTINUE
       allOK = .TRUE.
       DO 101 i = 1, numNod
            allOK = allOK.AND.checkN(i)
  101  CONTINUE
       IF (.NOT.allOK) THEN
		 ErrorMsg = "THE FOLLOWING NODES WERE NEVER READ:"
		 j = 1
         do i = 1, numNod
           if(i <= nRealN) THEN
             index = i
           else
             index = n1000 + i - nRealN
           end if
           if(.not.checkN(i)) then
		     ErrorArray(j) = index
			 j = j+1
           end if
		 end do
		 call FatalError(ErrorMsg,ThID,ErrArr=ErrorArray)
       END IF

!  Read triangular elements:

       READ (iUnit7, * ) numEl
       IF (numEl > mxEl) THEN
		 write(ErrorMsg,'(A,I0,A)') "INCREASE PARAMETER maxEl TO BE AT LEAST EQUAL TO THE NUMBER OF ELEMENTS (",numEl,") AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
       END IF
       DO 109 k = 1, numEl
            checkE(k) = .FALSE.
  109  CONTINUE
       IF (.NOT.brief) THEN
            IF(Verbose) WRITE (iUnitT, 110) numEl
  110       FORMAT(/' There are ',I6,' triangular continuum elements.'/ &
     &      ' (Node numbers for each are given at corners, counter', &
     &      'clockwise'/ / &
     &     ' element        C1        C2        C3 Lithospheric_Rheology#')
       END IF
       DO 200 k = 1, numEl
!   (Elements need not be input in order, but must all be present.)
            READ (iUnit7, "(A)") longer_line
            CALL Extract_LRi(longer_line, LRi, shorter_line)
            continuum_LRi(k) = LRi
            READ (shorter_line, * ) i, (nodes(j, i), j = 1, 3)
            IF ((i < 1).OR.(i > numEl)) THEN
			  write(ErrorMsg,'(A,I0)') "ERR0R: ILLEGAL ELEMENT NUMBER: ",i
			  call FatalError(ErrorMsg,ThID)
            END IF
            checkE(i) = .TRUE.
            IF (.NOT.brief) THEN
                 IF (LRi == 0) THEN
                     IF(Verbose) WRITE (iUnitT, 120) i, (nodes(j, i), j = 1, 3)
  120                FORMAT (' ', I6, ':', 3I10)
                 ELSE
                     IF(Verbose) WRITE (iUnitT, 121) i, (nodes(j, i), j = 1, 3), LRi
  121                FORMAT (' ', I6, ':', 3I10, ' LR', I8)
                 END IF
            END IF
            DO 130 j = 1, 3
                 n = nodes(j, i)
                 IF (n > nRealN) n = nRealN + (n - n1000)
                 IF ((n <= 0).OR.(n > numNod)) THEN
		           write(ErrorMsg,'(A,I0,A)') "NODE NUMBER ",nodes(j, i)," IS ILLEGAL."
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
		 ErrorMsg = "THE FOLLOWING ELEMENTS WERE NEVER READ:"
		 j = 1
         do i = 1, numEl
           if(.not.checkE(i)) then
		     ErrorArray(j) = i
			 j = j+1
		 end if
       end do
	   call FatalError(ErrorMsg,ThID,ErrArr=ErrorArray)
       END IF

!   Read fault elements:

       READ (iUnit7, * ) nFl
       IF (nFl > mxFEl) THEN
		 write(ErrorMsg,'(A,I0,A)') "INCREASE PARAMETER maxFEl TO BE AT LEAST EQUAL TO THE NUMBER OF FAULTS (",nFl,") AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
       END IF
       offMax = 0.0D0
       DO 222 i = 1, nFl
            checkF(i) = .FALSE.
  222  CONTINUE
       IF (.NOT.brief .AND. Verbose) WRITE(iUnitT, 230) nFl
  230  FORMAT(/ /' There are ', I6, ' great-circle fault elements.')
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
     &      '    offset Lithospheric_Rheology#'/)
  240       FORMAT (' ', I6, ':', 4I5, 1X, 2F6.1, 1X, F9.0)
       DO 300 k = 1, nFl
            off = 0.0D0
            READ(iUnit7, "(A)") longer_line
            CALL Extract_LRi(longer_line, LRi, shorter_line)
            fault_LRi(k) = LRi
            READ(shorter_line, * ) i, (nodeF(j, k), j = 1, 4), (dips(l), l = 1, 2), off
            IF ((i < 1).OR.(i > nFl)) THEN
			  write(ErrorMsg,'(A,I0,A,I0)') "ILLEGAL NODE NUMBER ", nodeF(j, i)," IN FAULT", i
			  call FatalError(ErrorMsg,ThID)
            END IF
            checkF(i) = .TRUE.
            IF (.NOT.brief) THEN
                 IF (LRi == 0) THEN
                     IF(Verbose) WRITE (iUnitT, 240) i, (nodeF(j, i), j = 1, 4), (dips(l), l = 1, 2), off
                 ELSE
                     IF(Verbose) WRITE (iUnitT, 242) i, (nodeF(j, i), j = 1, 4), (dips(l), l = 1, 2), off, LRi
  242                FORMAT (' ', I6, ':', 4I5, 1X, 2F6.1, 1X, F9.0, " LR", I8)
                 END IF
            END IF
            DO 250 j = 1, 4
                 n = nodeF(j, i)
                 IF (n > nRealN) n = nRealN + (n - n1000)
                 IF ((n <= 0).OR.(n > numNod)) THEN
			      write(ErrorMsg,'(A,I0,A,I0)') "ERR0R: ILLEGAL NODE NUMBER (", nodeF(j, i), ") IN FAULT ",i
			      call FatalError(ErrorMsg,ThID)
                 END IF
                 nodeF(j, i) = n
  250       CONTINUE
            DO 260 l = 1, 2
                 IF (ABS(dips(l)) > 90.0D0) THEN
		           write(ErrorMsg,'(A,F10.4,A/,A/,A)') "ILLEGAL DIP OF ", dips(l)," ; SHOULD BE IN RANGE OF -90. TO +90. DEGREES. ",&
				   & "(NOTE: ALL DIPS ARE IN DEGREES FROM THE HORIZONAL; A + pRefIX (OR NONE) INDICATES A DIP TOWARD THE n1-n2 SIDE; ", &
				   & "A - pRefIX INDICATES A DIP TOWARD THE n4-n3 SIDE.)"
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
		  write(ErrorMsg,'(A)') "THE FOLLOWING FAULTS WERE NEVER READ:"
		  j = 1
            do i = 1, nFl
              if(.not.checkF(i)) ErrorArray(j) = i
			  j = j+1
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
       IF (.NOT. brief .AND. Verbose) WRITE (iUnitT, 999)
  999  FORMAT (' --------------------------------------------------', &
     &          '-----------------------------')
       END SUBROUTINE GetNet


       SUBROUTINE GetPBx (iUnitM, iUnitT, names, nPBnd, nPlate, & ! input
     &                    nDPlat, pLat, pLon)                     ! output

!   Sets up arrays defining the plates in the PB2002 model of:
!      Bird [2003; G**3].

!  (The rotation vectors of the plates are contained in DATA
!      statements in the main PROGRAM.)

!   The digitized boundaries of the plates (continuous closed curves,
!      always circling counterclockwise, and redundantly describing
!      each plate boundary *twice*, from each side)
!      are read here, from an input file such as 'PB2002_plates.dig',
!      on Fortran input device iUnitM.

!   The convention for identifying the plates is a 2-character symbol.
!   See array "names" in the main PROGRAM.

!-------------------------------------------------------
       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER, INTENT(IN) :: iUnitM, iUnitT                                                 ! input
       CHARACTER*2, INTENT(IN) :: names                                                       ! input
       INTEGER, INTENT(IN) :: nPBnd, nPlate                                                   ! input
       INTEGER, INTENT(OUT) :: nDPlat                                                         ! output
       REAL*8, INTENT(OUT) :: pLat, pLon                                                      ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CHARACTER*2 symbol
       CHARACTER*3 stars
       INTEGER i, ios, ip, l, nRead
       DIMENSION names(nPlate), nDPlat(nPlate)
       DIMENSION pLat(nPlate, nPBnd), pLon(nPlate, nPBnd)
!------------------------------------------------------

       nRead = 0
  100  READ (iUnitM, 101, END = 201, IOSTAT = ios) symbol
            IF((nRead == 0).AND.(ios /= 0)) THEN
	           write(ErrorMsg,'(A)') "ERROR: File not found, or file empty."
		       call FatalError(ErrorMsg,ThID)
            END IF
  101       FORMAT (A2)
            DO 120 l = 1, nPlate
               IF(symbol == names(l)) THEN
                  ip = l
                  GO TO 140
               END IF
  120       CONTINUE
	        write(ErrorMsg,'(A,A,A,I3)') "ERR0R: BAD PLATE NAME",symbol," ON INPUT DEVICE ",iUnitM
		    call FatalError(ErrorMsg,ThID)
  140       nRead = nRead + 1
            IF (nRead > nPlate) THEN
	           write(ErrorMsg,'(A)') "Increase nPlate and recompile."
		       call FatalError(ErrorMsg,ThID)
            END IF
            i = 0
  142       READ (iUnitM, 145, END = 201) stars
  145            FORMAT (A3)
                 IF (stars == '***') THEN
                    nDPlat(ip) = i
                    GO TO 100
                 END IF
                 BACKSPACE iUnitM
                 i = i + 1
                 IF (i > nPBnd) THEN
	                write(ErrorMsg,'(A)') "Increase nPBnd and recompile."
		            call FatalError(ErrorMsg,ThID)
                 END IF
                 READ (iUnitM, * ) pLon(ip, i), pLat(ip, i)
                 pLon(ip, i) = pLon(ip, i) * 0.0174532925199433D0
                 pLat(ip, i) = pLat(ip, i) * 0.0174532925199433D0
            GO TO 142
  201  IF (nRead < nPlate) THEN
	     write(ErrorMsg,'(A,I3,A,I3)') "ERROR: Expecting ",nPlate," plates but read outlines of only ",nRead
		 call FatalError(ErrorMsg,ThID)
       END IF

       RETURN
       END SUBROUTINE GetPBx


      SUBROUTINE GetSLF(iUnitD, iUnitT, mxGSR, &              ! INTENT(IN)
     &                  fName, rLat, delVZ, fltLon, fltLat, & ! INTENT(OUT)
     &                  nFData)

       !READ in fault slip rates, for scoring purposes.

       IMPLICIT NONE
       INTEGER,                           INTENT(IN) :: iUnitD, iUnitT, mxGSR
       CHARACTER*80, DIMENSION(mxGSR),    INTENT(OUT) :: fName
       REAL*8,       DIMENSION(2, mxGSR), INTENT(OUT) :: rLat, delVZ
       REAL*8,       DIMENSION(mxGSR),    INTENT(OUT) :: fltLon, fltLat
       INTEGER,                           INTENT(OUT) :: nFData

       CHARACTER*8 :: c8r1, c8r2, c8v1, c8v2
       INTEGER     :: i, ios
       LOGICAL     :: dummy

!      Detech whether the file exists?

       nFData = 0
       IF(Verbose) WRITE(iUnitT,101) iUnitD
  101  FORMAT (/' ATTEMPTING TO READ FAULT SLIP RATES FROM UNIT', I3/)
       READ (iUnitD, * , IOSTAT = ios)
       IF (ios == 0) THEN
            REWIND(iUnitD)
       ELSE
            RETURN
       END IF

       DO 106 i = 1, mxGSR
            fName(nFData + 1) = ' '
            READ(iUnitD, "(A)", END = 111) fName(nFData + 1)

!     NOTE: The line with four slip-rate numbers (in mm/year)
!           is marked with initial letter "F".
!           This is to prevent program -Projector-
!           from interpreting such lines as (x, y) data.

!           When we read such lines here, we just assign the letter "F"
!           to LOGICAL DUMMY (and never use it).

            READ(iUnitD, * , END = 111) dummy, &
     &                                  rLat(1, nFData + 1), rLat(2, nFData + 1), &
     &                                  delVZ(1, nFData + 1), delVZ(2, nFData + 1)
            READ(iUnitD, * , END = 111) fltLon(nFData + 1), fltLat(nFData + 1)
            READ(iUnitD, * , END = 111)
            nFData = nFData + 1
  106  CONTINUE
	     write(ErrorMsg,'(A,I6,A)') "ERR0R: PARAMETER maxDat = ",mxGSR," IS NOT LARGE ENOUGH FOR THE NUMBER OF GEOLOGIC SLIP-RATE DATA."
		 call FatalError(ErrorMsg,ThID)
  111  IF (nFData == 0) THEN
            IF(Verbose) WRITE(iUnitT,112) iUnitD
  112       FORMAT (/' NO FAULT SLIP RATE DATA COULD BE READ FROM UNIT ', I2)
       ELSE
            IF(Verbose) WRITE(iUnitT,113) nFData, iUnitD
  113       FORMAT (/' ', I4, ' data on fault slip rates were read from unit ', I2, ':'/)
            IF(Verbose) WRITE(iUnitT,114)
  114       FORMAT ( &
     &     '                            ', &
     &     '       -in millimeters/year-        -in degrees-'/ &
     &     '                            ', &
     &     '    Min.    Max.    Min.    Max.      East    North'/ &
     &     ' Fault                       ', &
     &     'Dextral Dextral   Throw   Throw Longitude Latitude'/ &
     &     ' --------------------------- ', &
     &     '------- ------- ------- ------- --------- --------')
            DO 150 i = 1, nFData
                 IF (rLat(1, i) == 0.0D0) THEN
                      c8r1 = '       ?'
                 ELSE
                      WRITE (c8r1, "(F8.2)") rLat(1, i)
                 END IF
                 IF (rLat(2, i) == 0.0D0) THEN
                      c8r2 = '       ?'
                 ELSE
                      WRITE (c8r2, "(F8.2)") rLat(2, i)
                 END IF
                 IF (delVZ(1, i) == 0.0D0) THEN
                      c8v1 = '       ?'
                 ELSE
                      WRITE (c8v1, "(F8.2)") delVZ(1, i)
                 END IF
                 IF (delVZ(2, i) == 0.0D0) THEN
                      c8v2 = '       ?'
                 ELSE
                      WRITE (c8v2, "(F8.2)") delVZ(2, i)
                 END IF
                 IF(Verbose) WRITE(iUnitT,145)fName(i)(1:27), &
     &                            c8r1, c8r2, &
     &                            c8v1, c8v2, &
     &                            fltLon(i), fltLat(i)
  145            FORMAT(' ',A,4A8,F10.3,F9.3)
  150       CONTINUE
            IF(Verbose) WRITE(iUnitT,170)
  170       FORMAT (' NOTE: Throw-rates are always positive in', &
     &              ' cases of convergence (thrusting),' / &
     &              ' and negative in cases of divergence', &
     &              ' (normal faulting).')
       END IF
       END SUBROUTINE GetSLF


      SUBROUTINE GetSKS (iUnitK, iUnitT, mxSKS,  &              ! INTENT(IN)
     &                   numSKS, SKS_tag, SKS_theta, SKS_phi, & ! INTENT(OUT)
     &                   SKS_argument, SKS_delay)

!   READs data on upper-mantle seismic anisotropy, from SKS splitting.

!   Current FORMAT matches file "Fouch_upper-mantle_SKS_splitting-2004.dat":

!   If first READ fails (no file connected), returns numSKS = 0.

       USE DSphere ! in Peter Bird's file DSphere.f90
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iUnitK ! Fortran device number to read from
       INTEGER, INTENT(IN) :: iUnitT ! Fortran devide number for error messages
       INTEGER, INTENT(IN) :: mxSKS  ! maximum number of SKS data that can be stored
       INTEGER, INTENT(OUT) :: numSKS ! number of SKS data actually read
       CHARACTER*5, DIMENSION(mxSKS), INTENT(OUT) :: SKS_tag      ! 5-byte ID code for datum
       REAL*8,      DIMENSION(mxSKS), INTENT(OUT) :: SKS_theta    ! colatitude, in radians
       REAL*8,      DIMENSION(mxSKS), INTENT(OUT) :: SKS_phi      ! colatitude, in radians
       REAL*8,      DIMENSION(mxSKS), INTENT(OUT) :: SKS_argument ! phi or fast direction, in radians counterclockwise from S
       REAL*8,      DIMENSION(mxSKS), INTENT(OUT) :: SKS_delay    ! SKS splitting time, in s

       CHARACTER*5 :: tag
       INTEGER     :: azimuth, ios
       LOGICAL     :: problem
       REAL*8      :: dt, latitude, longitude

       IF(Verbose) WRITE(iUnitT,1) iUnitK
    1  FORMAT(/' ATTEMPTING TO READ FAST SKS AZIMUTHS AND DELAYS FROM UNIT', I3/)

       numSKS = 0

!   Do some test READs to see if a file is connected:

       problem = .FALSE.
       READ (iUnitK, "(A)", IOSTAT = ios) tag ! title line?
       problem = problem.OR.(ios /= 0)
       READ (iUnitK, "(A)", IOSTAT = ios) tag ! header line?
       problem = problem.OR.(ios /= 0)
       IF (problem) RETURN

!   Now proceed to read actual datum lines:
!   Main loop:

   10  READ (iUnitK, 32, END = 100) tag, latitude, longitude, azimuth, dt
   32       FORMAT (A5, 1X, 2F10.4, I8, F8.2)
            numSKS = numSKS + 1
            SKS_tag(numSKS) = tag
            SKS_phi(numSKS) = longitude * radians_per_degree
            SKS_theta(numSKS) = (90.0D0 - latitude) * radians_per_degree
            SKS_argument(numSKS) = (180.0D0 - azimuth) * radians_per_degree
            SKS_delay(numSKS) = dt
       IF (numSKS < mxSKS) GO TO 10
	     write(ErrorMsg,'(A,I6,A/,A)') "Subprogram GetSTR READ only",mxSKS," data points." , &
     &         "Increase PARAMETER maxSKS and recompile."
		 call FatalError(ErrorMsg,ThID)
  100  RETURN
       END SUBROUTINE GetSKS


      SUBROUTINE GetSTR (iUnitS, iUnitT, mxStr, &        ! INTENT(IN)
     &                   strTag, strThe, strPhi, &       ! INTENT(OUT)
     &                   strArg, strQua, strReg, numStr)

!   Reads most-compressive horizontal principal stress azimuths.

!   Current FORMAT matches file "WSM1092.LRECL116":
!      strTag = short 5-byte identification
!      strPhi = East longitude (radians)
!      strThe = co-latitude (radians)
!      strArg = azimuth of sigme_1H, in radians counterclockwise from South.
!      strQua = relative quality, converted from letter to REAL*8 weight.
!      strReg = 2-letter abbreviation for stress regime.

!   If first READ fails (no file connected), returns numStr = 0.

       USE DSphere ! in Peter Bird's file DSphere.f90
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iUnitS, iUnitT, mxStr
       CHARACTER*5, DIMENSION(mxStr), INTENT(OUT) :: strTag
       REAL*8,      DIMENSION(mxStr), INTENT(OUT) :: strThe, strPhi, strArg, strQua
       CHARACTER*2, DIMENSION(mxStr), INTENT(OUT) :: strReg
       INTEGER,                       INTENT(OUT) :: numStr

       CHARACTER*5 :: tag
       CHARACTER*2 :: regime
       CHARACTER*1 :: qualit
       REAL*8      :: eLon, nLat
       INTEGER     :: azimut, ios

       IF(Verbose) WRITE(iUnitT,1) iUnitS
    1  FORMAT(/' ATTEMPTING TO READ STRESS DIRECTIONS FROM UNIT', I3/)

       numStr = 0

!   Do one test READ to see if a file is connected?

       READ (iUnitS, 11, IOSTAT = ios) tag, nLat, eLon, azimut, qualit, regime
       IF (ios == 0) THEN
             BACKSPACE iUnitS
       ELSE
            RETURN
       END IF

!   Main loop:

   10  READ (iUnitS, 11, END = 100) tag, nLat, eLon, azimut, qualit, regime
   11       FORMAT (A5, 2F9.3, 7X, I3, 17X, A1, 5X, A2)
            IF ((azimut.LT.-180.0D0).OR.(azimut.GT.360.0D0)) THEN
	           write(ErrorMsg,'(A,I6/,A,A5,2F9.3,7X,I3,17X,A1,5X,A2)') "ERROR: Bad azimuth in stress data: ",azimut , &
		    &    "with: ", tag, nLat, eLon, azimut, qualit, regime
		       call FatalError(ErrorMsg,ThID)
            END IF
            numStr = numStr + 1
            strTag(numStr) = tag
            strPhi(numStr) = eLon * radians_per_degree
            strThe(numStr) = (90.0D0 - nLat) * radians_per_degree
            strArg(numStr) = (180.0D0 - azimut) * radians_per_degree
            IF (qualit == 'A') THEN
                 strQua(numStr) = 4.0D0
            ELSE IF (qualit == 'B') THEN
                 strQua(numStr) = 3.0D0
            ELSE IF (qualit == 'C') THEN
                 strQua(numStr) = 2.0D0
            ELSE IF (qualit == 'D') THEN
                 strQua(numStr) = 1.0D0
            ELSE
                 strQua(numStr) = 0.0D0
            END IF
            strReg(numStr) = regime
       IF (numStr < mxStr) GO TO 10
	     write(ErrorMsg,'(A,I6,A/,A)') "SUBPROGRAM GetSTR READ ONLY ",mxStr," AZIMUTHS." , &
     &              "INCREASE PARAMETER maxStr AND RECOMPILE."
		 call FatalError(ErrorMsg,ThID)
  100  RETURN
       END SUBROUTINE GetSTR


       SUBROUTINE Interp (fAtNod, mxEl, mxNode, nodes, numEl, & ! input
     &                    fAtIP)                                ! output

!   Interpolate a scalar function known at the nodes (fAtNod)
!   to values at the 7 integration points in each triangular
!   continuum element.  Note that simple linear interpolation in
!   a plane-triangle is used.  Thus, this routine is NOT suitable
!   for interpolating velocity vectors from nodes to integration points.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: fAtNod                                                           ! input
       INTEGER, INTENT(IN) :: mxEl, mxNode, nodes, numEl                                      ! input
       REAL*8, INTENT(OUT) :: fAtIP                                                           ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION points
       COMMON / S1S2S3 / points
       DIMENSION points(3, 7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER i, m
       DIMENSION fAtNod(mxNode), fAtIP(7, mxEl), nodes(3, mxEl)

       DO 100 m = 1, 7
            DO 90 i = 1, numEl
                 fAtIP(m, i) = points(1, m) * fAtNod(nodes(1, i)) + &
     &                         points(2, m) * fAtNod(nodes(2, i)) + &
     &                         points(3, m) * fAtNod(nodes(3, i))
   90       CONTINUE
  100  CONTINUE
       RETURN
       END SUBROUTINE Interp


       SUBROUTINE Likely (vMin, vMax, v, & ! INTENT(IN)
     &                    prob)            ! INTENT(OUT)

!       Compute the probability that velocity v is consistent with limits
!       vMin and vMax (both of which contain random error).

       IMPLICIT NONE
       REAL*8, INTENT(IN) :: vMin, vMax, v
       REAL*8, INTENT(OUT) :: prob

       INTEGER :: i
       LOGICAL :: boxcar, post, point
       REAL*8  :: b, badata, pg, pgf, pv, sdmodv, sdage, spread, tpim1h, v1, v2, vcent, w, &
                & d1, d2, db, dw, temp, v1u, v1l, v2u, v2l, vu, vl
       DATA badata / 0.20 / , sdmodv / 0.10 / , sdage / 0.15 / , tpim1h / 0.398942 /

       prob = 1.0D0
       IF ((vMin == 0.0D0) .AND. (vMax == 0.0D0)) RETURN
       v1 = vMin
       IF (vMin == 0.0D0) v1 = -1.0D38
       v2 = vMax
       IF (vMax == 0.0D0) v2 = + 1.0D38
       vcent = (v1 + v2) / 2.0D0
       spread = ABS((v2 - v1) / (v1 + v2))
       boxcar = (spread >= 1.3D0 * sdage).OR.(v1 == -1.0D38).OR.(v2 == 1.0D38)
       point = spread == 0.0D0
       post = (.NOT.boxcar).AND.(.NOT.point)
       vu = v * (1.0D0 + 2.0D0 * sdmodv)
       vl = v * (1.0D0 - 2.0D0 * sdmodv)
       IF(v >= 0.0D0) GO TO 4
       temp = vu
       vu = vl
       vl = temp
    4  v1u = v1 * (1.0D0 + 2.0D0 * sdage)
       v1l = v1 * (1.0D0 - 2.0D0 * sdage)
       IF(v1 >= 0.0D0) GO TO 6
       temp = v1u
       v1u = v1l
       v1l = temp
    6  v2u = v2 * (1.0D0 + 2.0D0 * sdage)
       v2l = v2 * (1.0D0 - 2.0D0 * sdage)
       IF(v2 >= 0.0D0) GO TO 8
       temp = v2u
       v2u = v2l
       v2l = temp
    8  IF(vu > v1l) GO TO 10
       prob = badata
       RETURN
   10  IF(vl < v2u) GO TO 20
       prob = badata
       RETURN
   20  IF((vl < v1u).OR.(vu > v2l)) GO TO 30
       prob = 1.0D0
       RETURN
   30  prob = 0.0D0
       b = -3.1D0
       db = 0.10D0
       d1 = ABS(v1) * sdage * 1.41421D0
       d2 = ABS(v2) * sdage * 1.41421D0
       dw = db * sdmodv * ABS(v)
       w = v - 31.0D0 * dw
       DO 40 i = 1, 61
       b = b + db
       pv = tpim1h * EXP(-0.5D0 * b**2)
       w = w + dw
       IF (point)  pgf = EXP(-MIN(99.9D0, 0.5D0 * ((w - v1) / (sdage * v1))**2))
       IF (post)   pgf = EXP(-MIN(99.9D0, (0.832D0 * (w - vcent) / (spread * vcent))**2))
       IF (boxcar) pgf = MIN(0.5D0 * (1.0D0 + ERF((w - v1) / d1)), 1.0D0, 0.5D0 * (1.0D0 + ERF((v2 - w) / d2)))
       pgf = MAX(pgf, 1.0D-50)
       pg = badata + (1. - badata) * pgf
   40  prob = prob + pv * pg * db
       END SUBROUTINE Likely


       SUBROUTINE Limits (area, detJ, iUnitT, mxEl, numEl, & ! input
     &                    okDelV, radius, refStr, sphere, tLInt, &
     &                    trHMax, zMoho, &
     &                    constr, etaMax, fMuMax, visMax)    ! output

!  Compute area, mean thickness, and other dimensional parameters
!  of the plate, then determine values of stiffness limits needed
!  to keep velocity errors down to order okDelV at shear stress
!  level refStr.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: area, detJ                                                       ! input
       INTEGER, INTENT(IN) :: iUnitT, mxEl, numEl                                             ! input
       REAL*8, INTENT(IN) :: okDelV, radius, refStr                                           ! input
       LOGICAL, INTENT(IN) :: sphere                                                          ! input
       REAL*8, INTENT(IN) :: tLInt, trHMax, zMoho                                             ! input
       REAL*8, INTENT(OUT) :: constr, etaMax, fMuMax, visMax                                  ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION weight
       COMMON / WgtVec / weight
       DIMENSION weight(7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER i, m, nfault
       REAL*8 dA, side, thick, totalA, totalV, whole
       DIMENSION area(mxEl), detJ(7, mxEl), tLInt(7, mxEl), zMoho(7, mxEl)

       totalA = 0.0D0
       totalV = 0.0D0
       DO 20 m = 1, 7
            DO 10 i = 1, numEl
                 da = area(i) * detJ(m, i) * weight(m)
                 totalA = totalA + da
                 totalV = totalV + da * (zMoho(m, i) + tLInt(m, i))
   10       CONTINUE
   20  CONTINUE
       whole = 4.0D0 * 3.14159265358979D0 * radius**2
       IF (totalA > (1.02D0 * whole)) THEN
	     write(ErrorMsg,'(A,ES12.4,A/,A,ES12.4,A/,A)') "AREA OF GRID (",totalA,") EXCEEDS" , &
     &              "AREA OF PLANET (",whole,"), WHICH MAKES NO SENSE." , &
     &              "CHECK GRID FOR ABS(LATITUDE) > 90. AND ALSO FOR OVERLAPPING ELEMENTS."
		 call FatalError(ErrorMsg,ThID)
       END IF
       thick = totalV / totalA
       IF (sphere) THEN
            side = radius
            nfault = 1
       ELSE
            side = SQRT(totalA)
            nfault = 4
       END IF
       constr = nfault * refStr * thick / okDelV
       etaMax = refStr * thick / (side * okDelV)
       etaMax = MIN(etaMax, trHMax / okDelV)
       fMuMax = nfault * refStr / okDelV
       visMax = 0.25D0 * refStr * side / okDelV
       IF(Verbose) WRITE(iUnitT,50) totalA, totalV, thick, side, constr, etaMax, &
     &                    fMuMax, visMax
   50  FORMAT (/ /' Subprogram -Limits- performed dimensional analysis'/ &
     & ' and estimated necessary stiffness limits to balance'/1P, &
     & ' the conflicting objectives of accuracy and precision:'/ / &
     & '                 area of model = ',D10.3,' length**2'/ &
     & '               volume of model = ',D10.3,' length**3'/ &
     & '             typical thickness = ',D10.3,' length'/ &
     & '                 typical width = ',D10.3,' length'/ &
     & '    constr (constraint weight) = ',D10.3,' force s/length**2'/ &
     & '  etaMax (max. basal coupling) = ',D10.3,' force s/length**3'/ &
     & ' fMuMax (max. fault stiffness) = ',D10.3,' force s/length**3'/ &
     & ' visMax (max. block viscosity) = ',D10.3,' force s/length**2')
       END SUBROUTINE Limits


       SUBROUTINE OLDMohr (aCreep, alphaT, bCreep, Biot, Byerly, & ! input
     &                  cCreep, cFric, conduc, constr, dCreep, dQdTdA, &
     &                  eCreep, elev, fDip, fFric, fMuMax, &
     &                  fPFlt, fArg, gMean, &
     &                  mxFEl, mxNode, nFl, nodeF, &
     &                  offMax, offset, &
     &                  oneKm, radio, rhoH2O, rhoBar, &
     &                  slide, tauMax, &
     &                  tLNode, tSurf, v, wedge, &
     &                  zMNode, &
     &                  zTranF, &                               ! modify
     &                  fC, fIMuDZ, fPeakS, fSlips, fTStar)     ! output

!   This subprogram contains the nonlinear rheology of the faults.
!   For each of 7 integration points along the length of each fault
!   element, it:

!   (1) Computes the slip-rate vector on the fault surface;
!   (2) Determines the shear stress on the fault surface by Mohr/
!       Coulomb/Navier theory; (This stress is proportional to depth,
!       so the calculation is actually done at unit depth and then
!       scaled.)
!   (3) Proceeds down the dip of the fault, checking temperature,
!       strain rate, and pressure to see if frictional or creep
!       shear stress is lower;
!   (4) Reports the vertical integral of "mu" (the ratio of shear
!       stress to slip rate) down the fault as fIMuDZ;
!      (Note that the integral is vertical, not on a slant, even though
!       conditions are evaluated along a slant path.)
!   (5) For dipping, oblique-slip faults only, also reports recommended
!       tactical values for the matrix fC and the vector fTStar
!       which jointly describe a linearized rheology stiffer than
!       the actual nonlinear rheology;
!   (6) zTranF is the latest estimate of the depths (in crust, in mantle lithosphere)
!       to the brittle/ductile transitions, at the fault midpoint;
!   (7) LOGICAL variable fSlips indicates whether the fault is
!       slipping at its midpoint;  (Otherwise, it is in the artificial
!       linearized regime, with stiffness fMuMax.)
!   (8) fPeakS gives the peak shear stress at the midpoint of each
!       fault, evaluated at the brittle/ductile transition;
!   (9) Faults with dip less than "slide" (radians) are limited
!       to a maximum down-dip integral shear traction of tauMax.

!   Note that pore pressures are considered in the calculation of
!   frictional strength:
!   *Normal pore pressures reduce the effective normal stress on the
!    fault surface by the amount
!                    -Biot * gMean * rhoH20 * z
!   *IF (offMax > 0.) THEN the remaining effective frictional strength
!    of the fault is multiplied by the reducing factor
!                    *(1. - Byerly * offset(i) / offMax).
!    This may also be a pore pressure effect, because Byerlee's model is
!    that gouge layers have thickness in proportion to offset, and
!    that they support non-Darcy static pore pressure gradients which
!    allow elevated pore pressures in the core of the gouge, which
!    reduce the effective friction of the fault.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: aCreep, alphaT, bCreep, Biot, Byerly, cCreep, cFric, conduc, &   ! input
          & constr, dCreep, dQdTdA, eCreep, elev, fDip, fFric, fMuMax, &                      ! input
          & fPFlt, fArg, gMean                                                                ! input
       INTEGER, INTENT(IN) :: mxFEl, mxNode, nFl, nodeF                                       ! input
       REAL*8, INTENT(IN) :: offMax, offset, oneKm, radio, rhoH2O, rhoBar, slide              ! input
       REAL*8, INTENT(IN) :: tauMax, tLNode, tSurf                                            ! input
       DOUBLE PRECISION, INTENT(IN) :: v                                                      ! input
       REAL*8, INTENT(IN) :: wedge, zMNode                                                    ! input
       REAL*8, INTENT(INOUT) :: zTranF                                                        ! modify
       REAL*8, INTENT(OUT) :: fC, fIMuDZ, fPeakS                                              ! output
       LOGICAL, INTENT(OUT) :: fSlips                                                         ! output
       REAL*8, INTENT(OUT) :: fTStar                                                          ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION fPhi
       COMMON / FPhis / fPhi
       DIMENSION fPhi(4, 7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Following PARAMETER gives number of steps in vertical integral
!   of creep shear stress on ductile parts of faults:
       INTEGER, PARAMETER :: nStep = 30
!   Higher values obviously cost more.  On the other hand, small values
!   do not merely approximate the creep law; they also introduce
!   some random error which can put a floor on convergence
!   of the whole global velocity field.

       INTEGER i, j, kIter, layer, limit, m, n1, n2, n3, n4
       REAL*8 angle, azfull, azhalf, baseZ, cGamma, close, cosr, crust, &
            & dDPNdZ, delVx, delVy, delTau, dEPdST, dip, dippy, dLEPdC, dLEPdZ, dPMax, dSFdZ, dz, &
            & efull, ehalf, elevat, ePMoho,  fric, huge, mantle, &
            & normal, oldsc, q, q0, rake, rho, rhoC, &
            & scfull, schalf, sf0, sFMoho, shearc, shearf, shearp, sheart, sinist, sinr, slip, spread, strain, sum, &
            & t, t0, t90pc, tand, tAsth, tfull, thalf, thick, thrust, tiny, tlan, tMean, tMeanC, tMoho, topZ, ts, tTrans, tu, tune, &
            & unitAx, unitAy, unitBx, unitBy, &
            & vIMuDZ, vitdz, vUpDip, whalf, wfull, z, z0, zAbs, zfull, zhalf, zman, zp, zTrans
       DOUBLE PRECISION dpt1
       LOGICAL locked, pureSS, sloped

!      DIMENSIONs of internal convenience arrays:
       DIMENSION dLEPdZ(2), dSFdZ(2), rho(2), sheart(2), tMean(2), zTrans(2)
!      DIMENSIONs of external argument arrays:
       DIMENSION aCreep(2), alphaT(2), bCreep(2), cCreep(2), conduc(2), &
     &           dCreep(2), dQdTdA(mxNode), elev(mxNode), &
     &           fC(2, 2, 7, mxFEl), fDip(2, mxFEl), &
     &           fIMuDZ(7, mxFEl), fPeakS(2, mxFEl), &
     &           fPFlt(2, 2, 2, 7, mxFEl), fSlips(mxFEl), &
     &           fArg(2, mxFEl), fTStar(2, 7, mxFEl), nodeF(4, mxFEl), &
     &           offset(mxFEl), radio(2), rhoBar(2), &
     &           tauMax(2), tLNode(mxNode), &
     &           v(2, mxNode), zMNode(mxNode), zTranF(2, mxFEl)

!   Following two numbers are "very small" and "very large", but not
!   so extreme as to cause underflow or overflow.  They may need to
!   be adjusted, depending on the computer and compiler you use:
       DATA tiny / 2.0D-38 /
       DATA huge / 1.0D+38 /

       cgamma = (1.0D0 + SIN(ATAN(cFric))) / (1.0D0 - SIN(ATAN(cFric)))
       DO 100 i = 1, nFl
            IF (offMax <= 0.) THEN
                 fric = fFric
            ELSE
                 fric = fFric * (1.0D0 - Byerly * offset(i) / offMax)
            END IF
            n1 = nodeF(1, i)
            n2 = nodeF(2, i)
            n3 = nodeF(3, i)
            n4 = nodeF(4, i)

!           Is this a purely strike-slip fault element?
            pureSS = (ABS(fDip(1, i) - 1.57079632679490D0) <= wedge).AND. &
     &               (ABS(fDip(2, i) - 1.57079632679490D0) <= wedge)

!           If so, compute estimate of relative normal stress
!          (relative to vertical stress) by using amount of divergence
!           between average of node n1 and n2 and average of node n3
!           and node n4 (in spite on constraint equation):
            IF (pureSS) THEN

!CCCC            angle = 0.5 * (fArg(1, i) + fArg(2, i))
!CCCC            Line above was replaced due to cycle-shift problem!

                 angle = Chord(fArg(1, i), 0.5D0, fArg(2, i))

                 unitBx = SIN(angle)
                 unitBy = -COS(angle)
!                unitB points outward on the n1-n2 side (away from
!                      the n3-n4 side).
                 delVx = v(1, n1) * fPFlt(1, 1, 1, 4, i) + v(2, n1) * fPFlt(2, 1, 1, 4, i) &
     &                 + v(1, n2) * fPFlt(1, 1, 2, 4, i) + v(2, n2) * fPFlt(2, 1, 2, 4, i) &
     &                 - v(1, n3) * fPFlt(1, 1, 2, 4, i) - v(2, n3) * fPFlt(2, 1, 2, 4, i) &
     &                 - v(1, n4) * fPFlt(1, 1, 1, 4, i) - v(2, n4) * fPFlt(2, 1, 1, 4, i)
                 delVy = v(1, n1) * fPFlt(1, 2, 1, 4, i) + v(2, n1) * fPFlt(2, 2, 1, 4, i) &
     &                 + v(1, n2) * fPFlt(1, 2, 2, 4, i) + v(2, n2) * fPFlt(2, 2, 2, 4, i) &
     &                 - v(1, n3) * fPFlt(1, 2, 2, 4, i) - v(2, n3) * fPFlt(2, 2, 2, 4, i) &
     &                 - v(1, n4) * fPFlt(1, 2, 1, 4, i) - v(2, n4) * fPFlt(2, 2, 1, 4, i)
!                delVx and delVy are the velocities of the n1-n2 side
!                      relative to the n3-n4 side.
                 spread = delVx * unitBx + delVy * unitBy
                 delTau = constr * spread
                 tlan = 0.5D0 * (tLNode(n1) + tLNode(n2))
                 zman = 0.5D0 * (zMNode(n1) + zMNode(n2))
                 IF ((tlan <= 0.0D0).OR.(zTranF(2, i) <= 0.0D0)) THEN
!                   Crust alone resists convergence:
                      dPMax = -2.0D0 * deltau / zTranF(1, i)
                      dDPNdZ = dpmax / zTranF(1, i)
                 ELSE
!                   Mantle lithosphere helps to resist convergence:
                      dDPNdZ = -deltau / &
     &                       (0.50D0 * zTranF(1, i)**2 + zTranF(2, i) * zman + &
     &                        0.50D0 * zTranF(2, i)**2)
                 END IF
!                dDPNdZ is the gradient of excess normal pressure (in
!                excess of vertical pressure) with depth on this fault;
!                check that it lies within frictional limits of blocks:
                 q = 0.250D0 * (dQdTdA(n1) + dQdTdA(n2) + &
     &                          dQdTdA(n3) + dQdTdA(n4))
                 tTrans = tSurf + zTranF(1, i) * q / conduc(1) - &
     &                    zTranF(1, i)**2 * radio(1) / (2. * conduc(1))
                 tMeanC = (tSurf + tTrans) / 2.0D0
                 rhoC = rhoBar(1) * (1.0D0 - alphaT(1) * tMeanC)
                 dLEPdC = gMean * (rhoC - rhoH2O * Biot)
                 thrust = dLEPdC * cGamma
                 normal = dLEPdC / cGamma
                 dDPNdZ = MAX(dDPNdZ, normal - dLEPdC)
                 dDPNdZ = MIN(dDPNdZ, thrust - dLEPdC)

            ELSE
!                Different logic will be used; this parameter is not
!                really needed.  Zero it just to be careful.
                 dDPNdZ = 0.0D0
            END IF

            DO 90 m = 1, 7

!                elevation:
                 elevat = elev(n1) * fPhi(1, m) + elev(n2) * fPhi(2, m)

!                heat flow:
                 q = dQdTdA(n1) * fPhi(1, m) + dQdTdA(n2) * fPhi(2, m)

!                crustal thickness:
                 crust = zMNode(n1) * fPhi(1, m) + zMNode(n2) * fPhi(2, m)

!                mantle lithosphere thickness:
                 mantle = tLNode(n1) * fPhi(1, m) + tLNode(n2) * fPhi(2, m)
                 mantle = MAX(mantle, 0.0D0)

!                Moho temperature:
                 tMoho = tSurf + crust * q / conduc(1) - &
     &                       crust**2 * radio(1) / (2.0D0 * conduc(1))

!                Temperature at base of plate:
                 tAsth = tMoho + mantle * (q - crust * radio(1)) / conduc(2) - &
     &                       mantle**2 * radio(2) / (2.0D0 * conduc(2))

!                mean temperatures:
                 tMean(1) = (tSurf + tMoho) / 2.0D0
                 tMean(2) = (tMoho + tAsth) / 2.0D0

!                mean densities:
                 rho(1) = rhoBar(1) * (1.0D0 - alphaT(1) * tMean(1))
                 rho(2) = rhoBar(2) * (1.0D0 - alphaT(2) * tMean(2))

!                derivitives of lithostatic effective pressure wrt depth
                 dLEPdZ(1) = gMean * (rho(1) - rhoH2O * Biot)
                 ePMoho = dLEPdZ(1) * crust
                 dLEPdZ(2) = gMean * (rho(2) - rhoH2O * Biot)

!               "angle" is the fault strike, in radians cclkws from +X.

!CCCC            angle = fArg(1, i) * fPhi(1, m) + fArg(2, i) * fPhi(2, m)
!CCCC            Line above was replaced due to cycle-shift problem!

                 angle = Chord(fArg(1, i), fPhi(2, m), fArg(2, i))

!                unitA is a unit vector along the fault, from n1 to n2.
                 unitAx = COS(angle)
                 unitAy = SIN(angle)

!                unitB is a perpendicular unit vector, pointing out
!                toward the n4-n3 side.
                 unitBx = -unitAy
                 unitBy = + unitAx

!                relative velocities are for n1-n2 side relative to
!                the n4-n3 side:
                 delVx = v(1, n1) * fPFlt(1, 1, 1, m, i) + v(2, n1) * fPFlt(2, 1, 1, m, i) &
     &                 + v(1, n2) * fPFlt(1, 1, 2, m, i) + v(2, n2) * fPFlt(2, 1, 2, m, i) &
     &                 - v(1, n3) * fPFlt(1, 1, 2, m, i) - v(2, n3) * fPFlt(2, 1, 2, m, i) &
     &                 - v(1, n4) * fPFlt(1, 1, 1, m, i) - v(2, n4) * fPFlt(2, 1, 1, m, i)
                 delVy = v(1, n1) * fPFlt(1, 2, 1, m, i) + v(2, n1) * fPFlt(2, 2, 1, m, i) &
     &                 + v(1, n2) * fPFlt(1, 2, 2, m, i) + v(2, n2) * fPFlt(2, 2, 2, m, i) &
     &                 - v(1, n3) * fPFlt(1, 2, 2, m, i) - v(2, n3) * fPFlt(2, 2, 2, m, i) &
     &                 - v(1, n4) * fPFlt(1, 2, 1, m, i) - v(2, n4) * fPFlt(2, 2, 1, m, i)

!                sinistral strike-slip rate component:
                 sinist = delVx * unitAx + delVy * unitAy

!                convergence rate component (in horizontal plane):
                 close = delVx * unitBx + delVy * unitBy

!                dip of the fault (from horizontal on the n1-n2 side):
                 dip = fDip(1, i) * fPhi(1, m) + fDip(2, i) * fPhi(2, m)
                 sloped = ABS(dip - 1.57079632679490D0) > wedge

                 IF (.NOT.sloped) THEN
!                     case of a near-vertical fault:
                      dSFdZ(1) = (dLEPdZ(1) + dDPNdZ) * fric
                      sFMoho = dSFdZ(1) * crust
                      dSFdZ(2) = (dLEPdZ(2) + dDPNdZ) * fric
                      slip = ABS(sinist)
                      locked = .FALSE.
                 ELSE
!                     case of a shallow-dipping fault:

!                     vUpDip is the up-dip velocity component, in the
!                     fault plane, of the block on the n1-n3 side:
                      vUpDip = close / COS(dip)

!                    "rake" angle is measured counterclockwise in
!                     fault plane from horizontal & parallel to "angle":
                      rake = ATan2F(vUpDip, sinist)

!                     derivitive of effective normal pressure
!                     with respect to shear traction on fault:
                      dEPdST = TAN(dip) * SIN(rake)
!                    (Notice that when sense of dip reverses, sign
!                     change caused by TAN(dip) is cancelled by sign
!                     change caused by SIN(rake).)

!                     According to theory, the equation to solve is:
!                        d(shear traction)/dz =
!                        fric*(dLEPdZ+dEPdST*d(shear_traction)/dz)
!                     This may have a physical solution (one with
!                     positive shear traction).  If not, the fault
!                     is locked.
                      locked = (fric * dEPdST) >= 1.00D0
                      IF (locked) THEN
                           dSFdZ(1) = huge
                           dSFdZ(2) = huge
                      ELSE
                           dSFdZ(1) = fric * dLEPdZ(1) / (1.00D0 - fric * dEPdST)
                           sFMoho = dSFdZ(1) * crust
                           dSFdZ(2) = fric * dLEPdZ(2) / (1.00D0 - fric * dEPdST)
                      END IF

                      slip = SQRT(sinist**2 + vUpDip**2)
                 END IF
                 slip = MAX(slip, 1.0D8 * tiny)

!                Locate plastic/creep transition(s)
!                by iterated halving of domain:

                 IF (mantle > 0.) THEN
                      limit = 2
                 ELSE
                      limit = 1
                      zTrans(2) = 0.0D0
                      sheart(2) = 0.0D0
                 END IF
                 DO 60 layer = 1, limit
                      topZ = 0.0D0
                      IF (layer == 1) THEN
                           baseZ = crust
                           sf0 = 0.0D0
                           t0 = tSurf
                           q0 = q
                           z0 = 0.0D0
                      ELSE
                           baseZ = mantle
                           sf0 = sFMoho
                           t0 = tMoho
                           q0 = q - crust * radio(1)
                           z0 = crust
                      END IF
                      DO 50 kIter = 1, 15
                           z = 0.50D0 * (topZ + baseZ)
                           zAbs = z + z0
                           shearf = z * dSFdZ(layer) + sf0
                           shearp = MIN(shearf, dCreep(layer))
                           t = t0 + q0 * z / conduc(layer) - (radio(layer) / &
     &                                          (2.0D0 * conduc(layer))) * z**2
                           IF (zAbs <= (15.0D0 * oneKm)) THEN
                                t90pc = 0.50D0 * zAbs
                           ELSE IF (zAbs < (45.0D0 * oneKm)) THEN
                                t90pc = (405.0D0 / 8.0D0) * oneKm + &
     &                                  (-7.0D0) * zAbs + &
     &                                  (13.0D0 / 40.0D0) * oneKm * (zAbs / oneKm)**2 + &
     &                                  (-1.0D0 / 300.0D0) * oneKm * (zAbs / oneKm)**3
                           ELSE
                                t90pc = 2.0D0 * zAbs
                           END IF
!                      See Turcotte et al (1980) JGR, 85, B11, 6224-6230.
                           strain = slip / t90pc
                           shearc = aCreep(layer) * (strain**eCreep) * &
     &                          EXP((bCreep(layer) + cCreep(layer) * z) / t)
                           IF (shearc < shearp) THEN
                                baseZ = z
                           ELSE
                                topZ = z
                           END IF
   50                 CONTINUE
                      zTrans(layer) = 0.50D0 * (topZ + baseZ)
                      sheart(layer) = zTrans(layer) * dSFdZ(layer) + sf0
   60            CONTINUE

!                plastic part of vertical integral(s) of traction:
!                (A) crust:
                 IF (sheart(1) <= dCreep(1)) THEN
                      vitdz = 0.50D0 * sheart(1) * zTrans(1)
                 ELSE
                      zp = zTrans(1) * dCreep(1) / sheart(1)
                      vitdz = dCreep(1) * (zTrans(1) - 0.50D0 * zp)
                 END IF
!                (B) mantle lithosphere:
                 IF ((mantle > 0.).AND.(sheart(2) > sfmoho)) THEN
                      IF (sheart(2) <= dCreep(2)) THEN
                           vitdz = vitdz + 0.50D0 * (sfmoho + sheart(2)) * zTrans(2)
                      ELSE
                           zp = zTrans(2) * (dCreep(2) - sfmoho) / &
     &                                  (sheart(2) - sfmoho)
                           zp = MAX(zp, 0.)
                           vitdz = vitdz + 0.50D0 * (sfmoho + sheart(2)) * zp + &
     &                                 dCreep(2) * (zTrans(2) - zp)
                      END IF
                 END IF

!                Add creep part(s) of integral, using parabolic rule:

                 sum = 0.0D0
                 DO 80 layer = 1, limit
                      IF (layer == 1) THEN
                           thick = crust
                           t0 = tSurf
                           q0 = q
                           zAbs = 0.0D0
                      ELSE
                           thick = mantle
                           t0 = tMoho
                           q0 = q - crust * radio(1)
                           zAbs = crust
                      END IF
                      dz = (thick - zTrans(layer)) / nStep
                      oldsc = sheart(layer)
                      oldsc = MIN(oldsc, dCreep(layer))
                      z0 = zTrans(layer)
                      DO 70 j = 1, nStep
                           zhalf = z0 + 0.50D0 * dz
                           zfull = z0 + dz
                           azhalf = zhalf + zAbs
                           azfull = zfull + zAbs
                           thalf = t0 + q0 * zhalf / conduc(layer) - &
     &                          (radio(layer) / &
     &                          (2.0D0 * conduc(layer))) * zhalf**2
                           tfull = t0 + q0 * zfull / conduc(layer) - &
     &                          (radio(layer) / &
     &                          (2.0D0 * conduc(layer))) * zfull**2
                           IF (azhalf <= (15.0D0 * oneKm)) THEN
                                whalf = 0.50D0 * azhalf
                           ELSE IF (azhalf < (45.0D0 * oneKm)) THEN
                                whalf = (405.0D0 / 8.0D0) * oneKm + &
     &                                (-7.0D0) * azhalf + &
     &                                (13.0D0 / 40.0D0) * oneKm * (azhalf / oneKm)**2 + &
     &                                (-1.0D0 / 300.0D0) * oneKm * (azhalf / oneKm)**3
                           ELSE
                                whalf = 2.0D0 * azhalf
                           END IF
                           IF (azfull <= (15.0D0 * oneKm)) THEN
                                wfull = 0.50D0 * azfull
                           ELSE IF (azfull < (45.0D0 * oneKm)) THEN
                                wfull = (405.0D0 / 8.0D0) * oneKm + &
     &                                (-7.0D0) * azfull + &
     &                                (13.0D0 / 40.0D0) * oneKm * (azfull / oneKm)**2 + &
     &                                (-1.0D0 / 300.0D0) * oneKm * (azfull / oneKm)**3
                           ELSE
                                wfull = 2.0D0 * azhalf
                           END IF
!                     See Turcotte et al (1980) JGR, 85, B11, 6224-6230.
                           ehalf = slip / whalf
                           efull = slip / wfull
                           schalf = aCreep(layer) * (ehalf**eCreep) * &
     &                          EXP((bCreep(layer) + cCreep(layer) * zhalf) &
     &                              / thalf)
                           schalf = MIN(schalf, dCreep(layer))
                           scfull = aCreep(layer) * (efull**eCreep) * &
     &                          EXP((bCreep(layer) + cCreep(layer) * zfull) &
     &                             / tfull)
                           scfull = MIN(scfull, dCreep(layer))
                           sum = sum + dz * (0.1666667D0 * oldsc + &
     &                                 0.6666667D0 * schalf + &
     &                                 0.1666666D0 * scfull)
                           z0 = zfull
                           oldsc = scfull
   70                 CONTINUE
   80            CONTINUE

                 vitdz = vitdz + sum

!                Limit shear traction on subduction zones only:

                 dippy = MIN(dip, 3.14159265358979D0 - dip)
                 IF (dippy <= slide) THEN
                      IF (elevat < 0.0D0) THEN
!                          apply oceanic subduction zone limit:
                           vitdz = MIN(vitdz, tauMax(1) * SIN(dip))
                      ELSE
!                          apply continental subduction zone limit:
                           vitdz = MIN(vitdz, tauMax(2) * SIN(dip))
                      END IF
                 END IF

                 dpt1 = (1.0D0 * vitdz) / slip
                 vIMuDZ = MIN(dpt1, 1.0D38)

                 fIMuDZ(m, i) = MIN(vimudz, fMuMax * (crust + mantle))

!                Dipping, oblique-slip integration
!                points are also characterized
!                by fC and fTStar:

                 IF (sloped) THEN
                      ts = sinist * fIMuDZ(m, i)
                      tu = vUpDip * fIMuDZ(m, i)
                      IF (locked) THEN
                           fC(1, 1, m, i) = fIMuDZ(m, i)
                           fC(1, 2, m, i) = 0.0D0
                           fC(2, 1, m, i) = 0.0D0
                           fC(2, 2, m, i) = fIMuDZ(m, i)
                      ELSE
                           sinr = SIN(rake)
                           cosr = COS(rake)
                           tand = TAN(dip)

!                      *** IMPORTANT NOTE: ***
!                          The following 7 statements are *not* the
!                          result of theory, but a tactical choice
!                          which attempts to compromise between
!                          stability of the linear system, stability
!                          of the iteration, and efficiency.
!                          They may be changed if the program does
!                          no converge satisfactorily!

                           tune = 2.0D0
                           fC(1, 1, m, i) = fIMuDZ(m, i) * &
     &                                (1.0D0 - tune * sinr * cosr**2 * tand)
                           fC(1, 2, m, i) = fIMuDZ(m, i) * &
     &                                (tune * cosr**3 * tand)
                           fC(2, 1, m, i) = fIMuDZ(m, i) * &
     &                                (-tune * sinr**2 * cosr * tand)
                           fC(2, 2, m, i) = fIMuDZ(m, i) * &
     &                                (1.0D0 + tune * sinr * cosr**2 * tand)
!                         (Often, fC(1,2) is the biggest term.
!                          In some cases, diagonals become negative.
!                          For stability, be sure that the fC
!                          matrix remains positive-definite:
                           fC(1, 1, m, i) = MAX(fC(1, 1, m, i), ABS(fC(1, 2, m, i)))
                           fC(2, 2, m, i) = MAX(fC(2, 2, m, i), ABS(fC(1, 2, m, i)))
                      END IF
                      fTStar(1, m, i) = ts - fC(1, 1, m, i) * sinist - &
     &                                 fC(1, 2, m, i) * vUpDip
                      fTStar(2, m, i) = tu - fC(2, 1, m, i) * sinist - &
     &                                 fC(2, 2, m, i) * vUpDip
                 END IF

!   Provide interesting diagnostic data at midpoints only:

                 IF (m == 4) THEN
                      fSlips(i) = (.NOT.locked).AND. &
     &                   (fIMuDZ(m, i) < (0.99D0 * fMuMax * (crust + mantle)))
                      zTranF(1, i) = zTrans(1)
                      fPeakS(1, i) = MIN(sheart(1), dCreep(1))
                      zTranF(2, i) = zTrans(2)
                      fPeakS(2, i) = MIN(sheart(2), dCreep(2))
                 END IF

   90       CONTINUE
  100  CONTINUE
       END SUBROUTINE OLDMohr


       SUBROUTINE Mohr (alphaT, conduc, constr, &                 ! input
     &                  continuum_LRi, &
     &                  dQdTdA, elev, &
     &                  fault_LRi, fDip, fMuMax, &
     &                  fPFlt, fArg, gMean, &
     &                  LRn, LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, &
     &                  LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
     &                  mxEl, mxFEl, mxNode, nFl, nodeF, &
     &                  offMax, offset, &
     &                  oneKm, radio, rhoH2O, rhoBar, &
     &                  slide, tauMax, &
     &                  tLNode, tSurf, v, wedge, &
     &                  zMNode, &
     &                  zTranF, &                               ! modify
     &                  fC, fIMuDZ, fPeakS, fSlips, fTStar)     ! output

!   This subprogram contains the nonlinear rheology of the faults.
!   For each of 7 integration points along the length of each fault
!   element, it:

!   (1) Computes the slip-rate vector on the fault surface;
!   (2) Determines the shear stress on the fault surface by Mohr/
!       Coulomb/Navier theory; (This stress is proportional to depth,
!       so the calculation is actually done at unit depth and then
!       scaled.)
!   (3) Proceeds down the dip of the fault, checking temperature,
!       strain rate, and pressure to see if frictional or creep
!       shear stress is lower;
!   (4) Reports the vertical integral of "mu" (the ratio of shear
!       stress to slip rate) down the fault as fIMuDZ;
!      (Note that the integral is vertical, not on a slant, even though
!       conditions are evaluated along a slant path.)
!   (5) For dipping, oblique-slip faults only, also reports recommended
!       tactical values for the matrix fC and the vector fTStar
!       which jointly describe a linearized rheology stiffer than
!       the actual nonlinear rheology;
!   (6) zTranF is the latest estimate of the depth
!       to the brittle/ductile transition, at the fault midpoint;
!   (7) LOGICAL variable fSlips indicates whether the fault is
!       slipping at its midpoint;  (Otherwise, it is in the artificial
!       linearized regime, with stiffness fMuMax.)
!   (8) fPeakS gives the peak shear stress at the midpoint of each
!       fault, evaluated at the brittle/ductile transition;
!   (9) Faults with dip less than "slide" (radians) are limited
!       to a maximum down-dip integral shear traction of tauMax.

!   Note that pore pressures are considered in the calculation of
!   frictional strength:
!   *Normal pore pressures reduce the effective normal stress on the
!    fault surface by the amount
!                    -Biot * gMean * rhoH20 * z
!   *IF (offMax > 0.) THEN the remaining effective frictional strength
!    of the fault is multiplied by the reducing factor
!                    *(1. - Byerly * offset(i) / offMax).
!    This may also be a pore pressure effect, because Byerlee's model is
!    that gouge layers have thickness in proportion to offset, and
!    that they support non-Darcy static pore pressure gradients which
!    allow elevated pore pressures in the core of the gouge, which
!    reduce the effective friction of the fault.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: alphaT, conduc, constr                                           ! input
       INTEGER, INTENT(IN) :: continuum_LRi                                                   ! input
       REAL*8, INTENT(IN) :: dQdTdA, elev                                                     ! input
       INTEGER, INTENT(IN) :: fault_LRi                                                       ! input
       REAL*8, INTENT(IN) :: fDip, fMuMax, fPFlt, fArg, gMean                                 ! input
       INTEGER, INTENT(IN) :: LRn                                                             ! input
       REAL*8, INTENT(IN) :: LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, &        ! input
                           & LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep ! input
       INTEGER, INTENT(IN) :: mxEl, mxFEl, mxNode, nFl, nodeF                                 ! input
       REAL*8, INTENT(IN) :: offMax, offset, oneKm, radio, rhoH2O, rhoBar, slide              ! input
       REAL*8, INTENT(IN) :: tauMax, tLNode, tSurf                                            ! input
       DOUBLE PRECISION, INTENT(IN) :: v                                                      ! input
       REAL*8, INTENT(IN) :: wedge, zMNode                                                    ! input
       REAL*8, INTENT(INOUT) :: zTranF                                                        ! modify
       REAL*8, INTENT(OUT) :: fC, fIMuDZ, fPeakS                                              ! output
       LOGICAL, INTENT(OUT) :: fSlips                                                         ! output
       REAL*8, INTENT(OUT) :: fTStar                                                          ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       DOUBLE PRECISION fPhi
       COMMON / FPhis / fPhi
       DIMENSION fPhi(4, 7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Following PARAMETER gives number of steps in vertical integral
!   of creep shear stress on ductile parts of faults:
       INTEGER, PARAMETER :: nStep = 30
!   Higher values obviously cost more.  On the other hand, small values
!   do not merely approximate the creep law; they also introduce
!   some random error which can put a floor on convergence
!   of the whole global velocity field.

       INTEGER i, j, kIter, layer, limit, m, n1, n2, n3, n4
       REAL*8 angle, azfull, azhalf, baseZ, cGamma, close, cosr, crust, &
            & dDPNdZ, delVx, delVy, delTau, dEPdST, dip, dippy, dLEPdC, dLEPdZ, dPMax, dSFdZ, dz, &
            & efull, ehalf, elevat, ePMoho,  fric, huge, mantle, &
            & normal, oldsc, q, q0, rake, rho, rhoC, &
            & scfull, schalf, sf0, sFMoho, shearc, shearf, shearp, sheart, sinist, sinr, slip, spread, strain, sum, &
            & t, t0, t90pc, tand, tAsth, tfull, thalf, thick, thrust, tiny, tlan, tMean, tMeanC, tMoho, topZ, ts, tTrans, tu, tune, &
            & unitAx, unitAy, unitBx, unitBy, &
            & vIMuDZ, vitdz, vUpDip, whalf, wfull, z, z0, zAbs, zfull, zhalf, zman, zp, zTrans
       DOUBLE PRECISION dpt1
       LOGICAL locked, pureSS, sloped

!      DIMENSIONs of external argument arrays:
       DIMENSION alphaT(2), conduc(2), &
     &           continuum_LRi(mxEl), &
     &           dQdTdA(mxNode), elev(mxNode), &
     &           fault_LRi(mxFEl), &
     &           fC(2, 2, 7, mxFEl), fDip(2, mxFEl), &
     &           fIMuDZ(7, mxFEl), fPeakS(2, mxFEl), &
     &           fPFlt(2, 2, 2, 7, mxFEl), fSlips(mxFEl), &
     &           fArg(2, mxFEl), fTStar(2, 7, mxFEl), &
     &           LR_set_fFric(0:LRn), LR_set_cFric(0:LRn), LR_set_Biot(0:LRn), LR_set_Byerly(0:LRn), &
     &           LR_set_aCreep(1:2, 0:LRn), LR_set_bCreep(1:2, 0:LRn), LR_set_cCreep(1:2, 0:LRn), LR_set_dCreep(1:2, 0:LRn), LR_set_eCreep(0:LRn), &
     &           nodeF(4, mxFEl), &
     &           offset(mxFEl), radio(2), rhoBar(2), &
     &           tauMax(2), tLNode(mxNode), &
     &           v(2, mxNode), zMNode(mxNode), zTranF(2, mxFEl)
!      DECLARATIONS and DIMENSIONs of internal convenience arrays:
       DIMENSION dLEPdZ(2), dSFdZ(2), rho(2), sheart(2), tMean(2), zTrans(2)
       INTEGER LRi
       REAL*8 t_fFric, t_cFric, t_Biot, t_Byerly, t_aCreep(2), t_bCreep(2), t_cCreep(2), t_dCreep(2), t_eCreep

!   Following two numbers are "very small" and "very large", but not
!   so extreme as to cause underflow or overflow.  They may need to
!   be adjusted, depending on the computer and compiler you use:
       DATA tiny / 2.D-38 /
       DATA huge / 1.D+38 /

       !Use default rheology to define the environment surrounding all faults:
       cGamma = (1.0D0 + SIN(ATAN(LR_set_cFric(0)))) / (1.0D0 - SIN(ATAN(LR_set_cFric(0))))
       DO 100 i = 1, nFl
            !Extract the desired rheology for this fault element:
            LRi = fault_LRi(i)
            t_fFric  = LR_set_fFric(LRi)
            t_Biot   = LR_set_Biot(LRi)
            t_Byerly = LR_set_Byerly(LRi)
            t_aCreep(1:2) = LR_set_aCreep(1:2, LRi)
            t_bCreep(1:2) = LR_set_bCreep(1:2, LRi)
            t_cCreep(1:2) = LR_set_cCreep(1:2, LRi)
            t_dCreep(1:2) = LR_set_dCreep(1:2, LRi)
            t_eCreep      = LR_set_eCreep(LRi)
            !- - - - - - - - - - - - - - - - - - - - - - --  - - -
            IF (offMax <= 0.0D0) THEN
                 fric = t_fFric
            ELSE
                 fric = t_fFric * (1.0D0 - t_Byerly * offset(i) / offMax)
            END IF
            n1 = nodeF(1, i)
            n2 = nodeF(2, i)
            n3 = nodeF(3, i)
            n4 = nodeF(4, i)

!           Is this a purely strike-slip fault element?
            pureSS = (ABS(fDip(1, i) - 1.57079632679490D0) <= wedge).AND. &
     &               (ABS(fDip(2, i) - 1.57079632679490D0) <= wedge)

!           If so, compute estimate of relative normal stress
!          (relative to vertical stress) by using amount of divergence
!           between average of node n1 and n2 and average of node n3
!           and node n4 (in spite on constraint equation):
            IF (pureSS) THEN

!CCCC            angle = 0.5 * (fArg(1, i) + fArg(2, i))
!CCCC            Line above was replaced due to cycle-shift problem!

                 angle = Chord(fArg(1, i), 0.5D0, fArg(2, i))

                 unitBx = SIN(angle)
                 unitBy = -COS(angle)
!                unitB points outward on the n1-n2 side (away from
!                      the n3-n4 side).
                 delVx = v(1, n1) * fPFlt(1, 1, 1, 4, i) + v(2, n1) * fPFlt(2, 1, 1, 4, i) &
     &                 + v(1, n2) * fPFlt(1, 1, 2, 4, i) + v(2, n2) * fPFlt(2, 1, 2, 4, i) &
     &                 - v(1, n3) * fPFlt(1, 1, 2, 4, i) - v(2, n3) * fPFlt(2, 1, 2, 4, i) &
     &                 - v(1, n4) * fPFlt(1, 1, 1, 4, i) - v(2, n4) * fPFlt(2, 1, 1, 4, i)
                 delVy = v(1, n1) * fPFlt(1, 2, 1, 4, i) + v(2, n1) * fPFlt(2, 2, 1, 4, i) &
     &                 + v(1, n2) * fPFlt(1, 2, 2, 4, i) + v(2, n2) * fPFlt(2, 2, 2, 4, i) &
     &                 - v(1, n3) * fPFlt(1, 2, 2, 4, i) - v(2, n3) * fPFlt(2, 2, 2, 4, i) &
     &                 - v(1, n4) * fPFlt(1, 2, 1, 4, i) - v(2, n4) * fPFlt(2, 2, 1, 4, i)
!                delVx and delVy are the velocities of the n1-n2 side
!                      relative to the n3-n4 side.
                 spread = delVx * unitBx + delVy * unitBy
                 delTau = constr * spread
                 tlan = 0.5D0 * (tLNode(n1) + tLNode(n2))
                 zman = 0.5D0 * (zMNode(n1) + zMNode(n2))
                 IF ((tlan <= 0.0D0).OR.(zTranF(2, i) <= 0.0D0)) THEN
!                   Crust alone resists convergence:
                      dPMax = -2.0D0 * deltau / zTranF(1, i)
                      dDPNdZ = dpmax / zTranF(1, i)
                 ELSE
!                   Mantle lithosphere helps to resist convergence:
                      dDPNdZ = -deltau / &
     &                       (0.50D0 * zTranF(1, i)**2 + zTranF(2, i) * zman + &
     &                        0.50D0 * zTranF(2, i)**2)
                 END IF
!                dDPNdZ is the gradient of excess normal pressure (in
!                excess of vertical pressure) with depth on this fault;
!                check that it lies within frictional limits of blocks:
                 q = 0.250D0 * (dQdTdA(n1) + dQdTdA(n2) + &
     &                          dQdTdA(n3) + dQdTdA(n4))
                 tTrans = tSurf + zTranF(1, i) * q / conduc(1) - &
     &                    zTranF(1, i)**2 * radio(1) / (2.0D0 * conduc(1))
                 tMeanC = (tSurf + tTrans) / 2.0D0
                 rhoC = rhoBar(1) * (1.0D0 - alphaT(1) * tMeanC)
                 dLEPdC = gMean * (rhoC - rhoH2O * t_Biot)
                 thrust = dLEPdC * cGamma
                 normal = dLEPdC / cGamma
                 dDPNdZ = MAX(dDPNdZ, normal - dLEPdC)
                 dDPNdZ = MIN(dDPNdZ, thrust - dLEPdC)

            ELSE
!                Different logic will be used; this parameter is not
!                really needed.  Zero it just to be careful.
                 dDPNdZ = 0.0D0
            END IF

            DO 90 m = 1, 7

!                elevation:
                 elevat = elev(n1) * fPhi(1, m) + elev(n2) * fPhi(2, m)

!                heat flow:
                 q = dQdTdA(n1) * fPhi(1, m) + dQdTdA(n2) * fPhi(2, m)

!                crustal thickness:
                 crust = zMNode(n1) * fPhi(1, m) + zMNode(n2) * fPhi(2, m)

!                mantle lithosphere thickness:
                 mantle = tLNode(n1) * fPhi(1, m) + tLNode(n2) * fPhi(2, m)
                 mantle = MAX(mantle, 0.0D0)

!                Moho temperature:
                 tMoho = tSurf + crust * q / conduc(1) - &
     &                       crust**2 * radio(1) / (2.0D0 * conduc(1))

!                Temperature at base of plate:
                 tAsth = tMoho + mantle * (q - crust * radio(1)) / conduc(2) - &
     &                       mantle**2 * radio(2) / (2.0D0 * conduc(2))

!                mean temperatures:
                 tMean(1) = (tSurf + tMoho) / 2.0D0
                 tMean(2) = (tMoho + tAsth) / 2.0D0

!                mean densities:
                 rho(1) = rhoBar(1) * (1.0D0 - alphaT(1) * tMean(1))
                 rho(2) = rhoBar(2) * (1.0D0 - alphaT(2) * tMean(2))

!                derivitives of lithostatic effective pressure wrt depth
                 dLEPdZ(1) = gMean * (rho(1) - rhoH2O * t_Biot)
                 ePMoho = dLEPdZ(1) * crust
                 dLEPdZ(2) = gMean * (rho(2) - rhoH2O * t_Biot)

!               "angle" is the fault strike, in radians cclkws from +X.

!CCCC            angle = fArg(1, i) * fPhi(1, m) + fArg(2, i) * fPhi(2, m)
!CCCC            Line above was replaced due to cycle-shift problem!

                 angle = Chord(fArg(1, i), fPhi(2, m), fArg(2, i))

!                unitA is a unit vector along the fault, from n1 to n2.
                 unitAx = COS(angle)
                 unitAy = SIN(angle)

!                unitB is a perpendicular unit vector, pointing out
!                toward the n4-n3 side.
                 unitBx = -unitAy
                 unitBy = + unitAx

!                relative velocities are for n1-n2 side relative to
!                the n4-n3 side:
                 delVx = v(1, n1) * fPFlt(1, 1, 1, m, i) + v(2, n1) * fPFlt(2, 1, 1, m, i) &
     &                 + v(1, n2) * fPFlt(1, 1, 2, m, i) + v(2, n2) * fPFlt(2, 1, 2, m, i) &
     &                 - v(1, n3) * fPFlt(1, 1, 2, m, i) - v(2, n3) * fPFlt(2, 1, 2, m, i) &
     &                 - v(1, n4) * fPFlt(1, 1, 1, m, i) - v(2, n4) * fPFlt(2, 1, 1, m, i)
                 delVy = v(1, n1) * fPFlt(1, 2, 1, m, i) + v(2, n1) * fPFlt(2, 2, 1, m, i) &
     &                 + v(1, n2) * fPFlt(1, 2, 2, m, i) + v(2, n2) * fPFlt(2, 2, 2, m, i) &
     &                 - v(1, n3) * fPFlt(1, 2, 2, m, i) - v(2, n3) * fPFlt(2, 2, 2, m, i) &
     &                 - v(1, n4) * fPFlt(1, 2, 1, m, i) - v(2, n4) * fPFlt(2, 2, 1, m, i)

!                sinistral strike-slip rate component:
                 sinist = delVx * unitAx + delVy * unitAy

!                convergence rate component (in horizontal plane):
                 close = delVx * unitBx + delVy * unitBy

!                dip of the fault (from horizontal on the n1-n2 side):
                 dip = fDip(1, i) * fPhi(1, m) + fDip(2, i) * fPhi(2, m)
                 sloped = ABS(dip - 1.57079632679490D0) > wedge

                 IF (.NOT.sloped) THEN
!                     case of a near-vertical fault:
                      dSFdZ(1) = (dLEPdZ(1) + dDPNdZ) * fric
                      sFMoho = dSFdZ(1) * crust
                      dSFdZ(2) = (dLEPdZ(2) + dDPNdZ) * fric
                      slip = ABS(sinist)
                      locked = .FALSE.
                 ELSE
!                     case of a shallow-dipping fault:

!                     vUpDip is the up-dip velocity component, in the
!                     fault plane, of the block on the n1-n3 side:
                      vUpDip = close / COS(dip)

!                    "rake" angle is measured counterclockwise in
!                     fault plane from horizontal & parallel to "angle":
                      rake = ATan2F(vUpDip, sinist)

!                     derivitive of effective normal pressure
!                     with respect to shear traction on fault:
                      dEPdST = TAN(dip) * SIN(rake)
!                    (Notice that when sense of dip reverses, sign
!                     change caused by TAN(dip) is cancelled by sign
!                     change caused by SIN(rake).)

!                     According to theory, the equation to solve is:
!                        d(shear traction)/dz =
!                        fric*(dLEPdZ+dEPdST*d(shear_traction)/dz)
!                     This may have a physical solution (one with
!                     positive shear traction).  If not, the fault
!                     is locked.
                      locked = ((fric * dEPdST) >= 1.00D0)
                      IF (locked) THEN
                           dSFdZ(1) = huge
                           dSFdZ(2) = huge
                      ELSE
                           dSFdZ(1) = fric * dLEPdZ(1) / (1.00D0 - fric * dEPdST)
                           sFMoho = dSFdZ(1) * crust
                           dSFdZ(2) = fric * dLEPdZ(2) / (1.00D0 - fric * dEPdST)
                      END IF

                      slip = SQRT(sinist**2 + vUpDip**2)
                 END IF
                 slip = MAX(slip, 1.0D8 * tiny)

!                Locate plastic/creep transition(s)
!                by iterated halving of domain:

                 IF (mantle > 0.0D0) THEN
                      limit = 2
                 ELSE
                      limit = 1
                      zTrans(2) = 0.0D0
                      sheart(2) = 0.0D0
                 END IF
                 DO 60 layer = 1, limit
                      topZ = 0.0D0
                      IF (layer == 1) THEN
                           baseZ = crust
                           sf0 = 0.0D0
                           t0 = tSurf
                           q0 = q
                           z0 = 0.0D0
                      ELSE
                           baseZ = mantle
                           sf0 = sFMoho
                           t0 = tMoho
                           q0 = q - crust * radio(1)
                           z0 = crust
                      END IF
                      DO 50 kIter = 1, 15
                           z = 0.50D0 * (topZ + baseZ)
                           zAbs = z + z0
                           shearf = z * dSFdZ(layer) + sf0
                           shearp = MIN(shearf, t_dCreep(layer))
                           t = t0 + q0 * z / conduc(layer) - (radio(layer) / &
     &                                          (2.0D0 * conduc(layer))) * z**2
                           IF (zAbs <= (15.0D0 * oneKm)) THEN
                                t90pc = 0.50D0 * zAbs
                           ELSE IF (zAbs < (45.0D0 * oneKm)) THEN
                                t90pc = (405.0D0 / 8.0D0) * oneKm + &
     &                                  (-7.0D0) * zAbs + &
     &                                  (13.0D0 / 40.0D0) * oneKm * (zAbs / oneKm)**2 + &
     &                                  (-1.0D0 / 300.0D0) * oneKm * (zAbs / oneKm)**3
                           ELSE
                                t90pc = 2.0D0 * zAbs
                           END IF
!                      See Turcotte et al (1980) JGR, 85, B11, 6224-6230.
                           strain = slip / t90pc
                           shearc = t_aCreep(layer) * (strain**t_eCreep) * &
     &                          EXP((t_bCreep(layer) + t_cCreep(layer) * z) / t)
                           IF (shearc < shearp) THEN
                                baseZ = z
                           ELSE
                                topZ = z
                           END IF
   50                 CONTINUE
                      zTrans(layer) = 0.50D0 * (topZ + baseZ)
                      sheart(layer) = zTrans(layer) * dSFdZ(layer) + sf0
   60            CONTINUE

!                plastic part of vertical integral(s) of traction:
!                (A) crust:
                 IF (sheart(1) <= t_dCreep(1)) THEN
                      vitdz = 0.50D0 * sheart(1) * zTrans(1)
                 ELSE
                      zp = zTrans(1) * t_dCreep(1) / sheart(1)
                      vitdz = t_dCreep(1) * (zTrans(1) - 0.50D0 * zp)
                 END IF
!                (B) mantle lithosphere:
                 IF ((mantle > 0.).AND.(sheart(2) > sfmoho)) THEN
                      IF (sheart(2) <= t_dCreep(2)) THEN
                           vitdz = vitdz + 0.50D0 * (sfmoho + sheart(2)) * zTrans(2)
                      ELSE
                           zp = zTrans(2) * (t_dCreep(2) - sfmoho) / &
     &                                  (sheart(2) - sfmoho)
                           zp = MAX(zp, 0.)
                           vitdz = vitdz + 0.50D0 * (sfmoho + sheart(2)) * zp + &
     &                                 t_dCreep(2) * (zTrans(2) - zp)
                      END IF
                 END IF

!                Add creep part(s) of integral, using parabolic rule:

                 sum = 0.0D0
                 DO 80 layer = 1, limit
                      IF (layer == 1) THEN
                           thick = crust
                           t0 = tSurf
                           q0 = q
                           zAbs = 0.0D0
                      ELSE
                           thick = mantle
                           t0 = tMoho
                           q0 = q - crust * radio(1)
                           zAbs = crust
                      END IF
                      dz = (thick - zTrans(layer)) / nStep
                      oldsc = sheart(layer)
                      oldsc = MIN(oldsc, t_dCreep(layer))
                      z0 = zTrans(layer)
                      DO 70 j = 1, nStep
                           zhalf = z0 + 0.50D0 * dz
                           zfull = z0 + dz
                           azhalf = zhalf + zAbs
                           azfull = zfull + zAbs
                           thalf = t0 + q0 * zhalf / conduc(layer) - &
     &                          (radio(layer) / &
     &                          (2.0D0 * conduc(layer))) * zhalf**2
                           tfull = t0 + q0 * zfull / conduc(layer) - &
     &                          (radio(layer) / &
     &                          (2.0D0 * conduc(layer))) * zfull**2
                           IF (azhalf <= (15.0D0 * oneKm)) THEN
                                whalf = 0.50D0 * azhalf
                           ELSE IF (azhalf < (45.0D0 * oneKm)) THEN
                                whalf = (405.0D0 / 8.0D0) * oneKm + &
     &                                (-7.0D0) * azhalf + &
     &                                (13.0D0 / 40.0D0) * oneKm * (azhalf / oneKm)**2 + &
     &                                (-1.0D0 / 300.0D0) * oneKm * (azhalf / oneKm)**3
                           ELSE
                                whalf = 2.0D0 * azhalf
                           END IF
                           IF (azfull <= (15.0D0 * oneKm)) THEN
                                wfull = 0.50D0 * azfull
                           ELSE IF (azfull < (45.0D0 * oneKm)) THEN
                                wfull = (405.0D0 / 8.0D0) * oneKm + &
     &                                (-7.0D0) * azfull + &
     &                                (13.0D0 / 40.0D0) * oneKm * (azfull / oneKm)**2 + &
     &                                (-1.0D0 / 300.0D0) * oneKm * (azfull / oneKm)**3
                           ELSE
                                wfull = 2.0D0 * azhalf
                           END IF
!                     See Turcotte et al (1980) JGR, 85, B11, 6224-6230.
                           ehalf = slip / whalf
                           efull = slip / wfull
                           schalf = t_aCreep(layer) * (ehalf**t_eCreep) * &
     &                          EXP((t_bCreep(layer) + t_cCreep(layer) * zhalf) &
     &                              / thalf)
                           schalf = MIN(schalf, t_dCreep(layer))
                           scfull = t_aCreep(layer) * (efull**t_eCreep) * &
     &                          EXP((t_bCreep(layer) + t_cCreep(layer) * zfull) &
     &                             / tfull)
                           scfull = MIN(scfull, t_dCreep(layer))
                           sum = sum + dz * (0.1666667D0 * oldsc + &
     &                                 0.6666667D0 * schalf + &
     &                                 0.1666666D0 * scfull)
                           z0 = zfull
                           oldsc = scfull
   70                 CONTINUE
   80            CONTINUE

                 vitdz = vitdz + sum

!                Limit shear traction on subduction zones only:

                 dippy = MIN(dip, 3.14159265358979D0 - dip)
                 IF (dippy <= slide) THEN
                      IF (elevat < 0.0D0) THEN
!                          apply oceanic subduction zone limit:
                           vitdz = MIN(vitdz, tauMax(1) * SIN(dip))
                      ELSE
!                          apply continental subduction zone limit:
                           vitdz = MIN(vitdz, tauMax(2) * SIN(dip))
                      END IF
                 END IF

                 dpt1 = (1.D0 * vitdz) / slip
                 vIMuDZ = MIN(dpt1, 1.D38)

                 fIMuDZ(m, i) = MIN(vimudz, fMuMax * (crust + mantle))

!                Dipping, oblique-slip integration
!                points are also characterized
!                by fC and fTStar:

                 IF (sloped) THEN
                      ts = sinist * fIMuDZ(m, i)
                      tu = vUpDip * fIMuDZ(m, i)
                      IF (locked) THEN
                           fC(1, 1, m, i) = fIMuDZ(m, i)
                           fC(1, 2, m, i) = 0.0D0
                           fC(2, 1, m, i) = 0.0D0
                           fC(2, 2, m, i) = fIMuDZ(m, i)
                      ELSE
                           sinr = SIN(rake)
                           cosr = COS(rake)
                           tand = TAN(dip)

!                      *** IMPORTANT NOTE: ***
!                          The following 7 statements are *not* the
!                          result of theory, but a tactical choice
!                          which attempts to compromise between
!                          stability of the linear system, stability
!                          of the iteration, and efficiency.
!                          They may be changed if the program does
!                          no converge satisfactorily!

                           tune = 2.0D0
                           fC(1, 1, m, i) = fIMuDZ(m, i) * &
     &                                (1.0D0 - tune * sinr * cosr**2 * tand)
                           fC(1, 2, m, i) = fIMuDZ(m, i) * &
     &                                (tune * cosr**3 * tand)
                           fC(2, 1, m, i) = fIMuDZ(m, i) * &
     &                                (-tune * sinr**2 * cosr * tand)
                           fC(2, 2, m, i) = fIMuDZ(m, i) * &
     &                                (1.0D0 + tune * sinr * cosr**2 * tand)
!                         (Often, fC(1,2) is the biggest term.
!                          In some cases, diagonals become negative.
!                          For stability, be sure that the fC
!                          matrix remains positive-definite:
                           fC(1, 1, m, i) = MAX(fC(1, 1, m, i), ABS(fC(1, 2, m, i)))
                           fC(2, 2, m, i) = MAX(fC(2, 2, m, i), ABS(fC(1, 2, m, i)))
                      END IF
                      fTStar(1, m, i) = ts - fC(1, 1, m, i) * sinist - &
     &                                 fC(1, 2, m, i) * vUpDip
                      fTStar(2, m, i) = tu - fC(2, 1, m, i) * sinist - &
     &                                 fC(2, 2, m, i) * vUpDip
                 END IF

!   Provide interesting diagnostic data at midpoints only:

                 IF (m == 4) THEN
                      fSlips(i) = (.NOT.locked).AND. &
     &                   (fIMuDZ(m, i) < (0.99D0 * fMuMax * (crust + mantle)))
                      zTranF(1, i) = zTrans(1)
                      fPeakS(1, i) = MIN(sheart(1), t_dCreep(1))
                      zTranF(2, i) = zTrans(2)
                      fPeakS(2, i) = MIN(sheart(2), t_dCreep(2))
                 END IF

   90       CONTINUE
  100  CONTINUE
       RETURN
       END SUBROUTINE Mohr


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
       RETURN
       END SUBROUTINE Next


       SUBROUTINE OldVel (iUnitT, iUnitV, mxNode, numNod, &  ! INTENT(IN)
     &                    haveNV, title1, title2, title3, v) ! INTENT(OUT)

!   READ old velocity solution from unit iUnitV, or else set LOGICAL
!     variable 'haveNV" to .FALSE.
!   WRITEs 3 title lines to iUnitT.

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iUnitT, iUnitV, mxNode, numNod
       LOGICAL,                      INTENT(OUT) :: haveNV
       CHARACTER*100,                 INTENT(OUT) :: title1, title2, title3
       REAL*8, DIMENSION(2, mxNode), INTENT(OUT) :: v

       INTEGER :: i, j

       READ (iUnitV, '(A80)', END = 100, ERR = 100) title1
       READ (iUnitV, '(A80)', END = 100, ERR = 100) title2
       READ (iUnitV, '(A80)', END = 100, ERR = 100) title3
       READ (iUnitV, * , END = 100, ERR = 100) ((v(j, i), j = 1, 2), i = 1, numNod)
       haveNV = .TRUE.
       IF(Verbose) WRITE(iUnitT,50) iUnitV, title1, title2, title3
   50  FORMAT (/ /' Velocity solution was', &
     &           ' READ from unit', I3, '; titles were:' / 3(/ ' ', A80) )
       GO TO 900
! ------------------(This section executed only if READ fails)---------
  100  IF(Verbose) WRITE(iUnitT,110) iUnitV
  110  FORMAT (/ /' NO FURTHER VELOCITY SOLUTIONS FOUND ON UNIT', I3)
       haveNV = .FALSE.
! ---------------------------------------------------------------------
  900  IF(Verbose) WRITE(iUnitT,999)
  999  FORMAT (' -------------------------------------------------------------------------------')
       END SUBROUTINE OldVel


       SUBROUTINE OnArc (s, theta1, phi1, theta2, phi2, & ! INTENT(IN)
     &                   theta, phi)                      ! INTENT(OUT)

!   Computes coordinates (theta, phi) = (colatitude, longitude) in radians
!   for a point which lies on the great circle arc from (theta1, phi1) to (theta2, phi2).
!   Input parameter s expresses the fractional distance,
!   so input s = 0.0D0 returns (theta1, phi1) and
!      input s = 1.0D0 returns (theta2, phi2).

       IMPLICIT NONE
       REAL*8, INTENT(IN) :: s, theta1, phi1, theta2, phi2
       REAL*8, INTENT(OUT) :: theta, phi

       REAL*8 :: comple, equat, size
       REAL*8 :: tvec(3), uvec(3), uvec1(3), uvec2(3)
!     These are all unit vectors (in the unit sphere)
!     in a Cartesian coordinate system with
!     X = (0E, 0N), Y = (90E, 0N), Z = North pole.

       uvec1(1) = SIN(theta1) * COS(phi1)
       uvec1(2) = SIN(theta1) * SIN(phi1)
       uvec1(3) = COS(theta1)

       uvec2(1) = SIN(theta2) * COS(phi2)
       uvec2(2) = SIN(theta2) * SIN(phi2)
       uvec2(3) = COS(theta2)

       comple = 1.0D0 - s
       tvec(1) = comple * uvec1(1) + s * uvec2(1)
       tvec(2) = comple * uvec1(2) + s * uvec2(2)
       tvec(3) = comple * uvec1(3) + s * uvec2(3)

       size = SQRT(tvec(1)**2 + tvec(2)**2 + tvec(3)**2)
       uvec(1) = tvec(1) / size
       uvec(2) = tvec(2) / size
       uvec(3) = tvec(3) / size

       equat = SQRT(uvec(1)**2 + uvec(2)**2)
       theta = ATAN2(equat, uvec(3))
       phi = ATan2F(uvec(2), uvec(1))
       END SUBROUTINE OnArc

       SUBROUTINE Prince (e11, e22, e12, &            ! input
     &                    e1, e2, u1x, u1y, u2x, u2y) ! output

!   Find principal values (e1, e2) of the symmetric 2x2 tensor e11 e12
!                                                              e12 e22
!   and also the associated eigenvectors #1=(u1x, u1y); #2=(u2x, u2y).
!   The convention is that e1 <= e2.

       IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8, INTENT(IN) :: e11, e22, e12                                                    ! input
       REAL*8, INTENT(OUT) :: e1, e2, u1x, u1y, u2x, u2y                                      ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       REAL*8 c, r, scale, test, theta

       r = SQRT(((e11 - e22) * 0.5D0)**2 + e12**2)
       c = (e11 + e22) * 0.5D0
       e1 = c - r
       e2 = c + r
       scale = MAX(ABS(e1), ABS(e2))
       test = 0.01D0 * scale
       IF ((ABS(e12) > test).OR.(ABS(e11 - e1) > test)) THEN
            theta = ATan2F(e11 - e1, -e12)
       ELSE
            theta = ATan2F(e12,  e1 - e22)
       END IF
       u1x = COS(theta)
       u1y = SIN(theta)
       u2x = u1y
       u2y = -u1x
       RETURN
       END SUBROUTINE Prince


      SUBROUTINE Projec (iEle, mxEl, mxNode, nodes, &
     &                    s1, s2, s3, &
     &                    xNode, yNode, v, &    ! INTENT(IN)
     &                    vTheta, vPhi)         ! INTENT(OUT)

!   Computes horizontal velocity (vTheta,vPhi) at point (s1, s2, s3)
!   in triangular element iEle.

       IMPLICIT NONE
       INTEGER,                       INTENT(IN) :: iEle, mxEl, mxNode
       INTEGER, DIMENSION(3, mxEl),   INTENT(IN) :: nodes
       REAL*8,                        INTENT(IN) :: s1, s2, s3
       REAL*8,  DIMENSION(mxNode),    INTENT(IN) :: xNode, yNode
       REAL*8,  DIMENSION(2, mxNode), INTENT(IN) :: v
       REAL*8,                        INTENT(OUT) :: vTheta, vPhi

       INTEGER :: k
       REAL*8 :: growth, size
       REAL*8 :: rn(3, 3), rP(3), rS(3), &
     &           uTheta(3, 3), uPhi(3, 3), uT(3), uP(3), &
     &           vCart(3, 3), vP(3), vS(3)

       DO 10 k = 1, 3
!          Find Cartesian radius vector for each node:
            rn(1, k) = SIN(xNode(nodes(k, iEle))) * COS(yNode(nodes(k, iEle)))
            rn(2, k) = SIN(xNode(nodes(k, iEle))) * SIN(yNode(nodes(k, iEle)))
            rn(3, k) = COS(xNode(nodes(k, iEle)))
!          Find Theta-pointing Cartesian unit vector at each node:
            uTheta(1, k) = COS(xNode(nodes(k, iEle))) * COS(yNode(nodes(k, iEle)))
            uTheta(2, k) = COS(xNode(nodes(k, iEle))) * SIN(yNode(nodes(k, iEle)))
            uTheta(3, k) = -SIN(xNode(nodes(k, iEle)))
!          Find Phi-pointing Cartesian unit vector at each node:
            uPhi(1, k) = -SIN(yNode(nodes(k, iEle)))
            uPhi(2, k) =  COS(yNode(nodes(k, iEle)))
            uPhi(3, k) = 0.0D0
!          Create Cartesian velocity at each node:
            vCart(1, k) = v(1, nodes(k, iEle)) * uTheta(1, k) + &
     &                    v(2, nodes(k, iEle)) * uPhi(1, k)
            vCart(2, k) = v(1, nodes(k, iEle)) * uTheta(2, k) + &
     &                    v(2, nodes(k, iEle)) * uPhi(2, k)
            vCart(3, k) = v(1, nodes(k, iEle)) * uTheta(3, k) + &
     &                    v(2, nodes(k, iEle)) * uPhi(3, k)
   10  CONTINUE
!     Interpolate velocity in plane triangle:
       vP(1) = s1 * vCart(1, 1) + s2 * vCart(1, 2) + s3 * vCart(1, 3)
       vP(2) = s1 * vCart(2, 1) + s2 * vCart(2, 2) + s3 * vCart(2, 3)
       vP(3) = s1 * vCart(3, 1) + s2 * vCart(3, 2) + s3 * vCart(3, 3)
!     Interpolate position in plane triangle:
       rP(1) = s1 * rn(1, 1) + s2 * rn(1, 2) + s3 * rn(1, 3)
       rP(2) = s1 * rn(2, 1) + s2 * rn(2, 2) + s3 * rn(2, 3)
       rP(3) = s1 * rn(3, 1) + s2 * rn(3, 2) + s3 * rn(3, 3)
       size = SQRT(rP(1)**2 + rP(2)**2 + rP(3)**2)
       growth = 1.0D0 / size
!     Project position to sphere:
       rS(1) = rP(1) * growth
       rS(2) = rP(2) * growth
       rS(3) = rP(3) * growth
!     Project velocity to velocity on sphere:
       vS(1) = vP(1) * growth
       vS(2) = vP(2) * growth
       vS(3) = vP(3) * growth
!     Find Theta and Phi unit vectors at this point:
       uP(1) = -rS(2)
       uP(2) =  rS(1)
       uP(3) = 0.0D0
       CALL Unit (uP)
       CALL DCross (uP, rS, uT)
!     Find horizontal components of velocity:
       vTheta = vS(1) * uT(1) + vS(2) * uT(2) + vS(3) * uT(3)
       vPhi   = vS(1) * uP(1) + vS(2) * uP(2) + vS(3) * uP(3)
       END SUBROUTINE Projec


       SUBROUTINE PutNet (iUnitO, &                           ! INTENT(IN)
     &                    brief, dQdTdA, elev, fDip, &
     &                    mxEl, mxFEl, mxNode, n1000, &
     &                    nFakeN, nFl, nodeF, nodes, &
     &                    nRealN, numEl, numNod, offset, &
     &                    title1, tLNode, xNode, yNode, zMNode)

!   WRITEs finite element grid to unit "iUnitO".

       USE DSphere ! in Peter Bird's file DSphere.f90
       IMPLICIT NONE
       INTEGER,                      INTENT(IN) :: iUnitO
       LOGICAL,                      INTENT(IN) :: brief
       REAL*8,  DIMENSION(mxNode),   INTENT(IN) :: dQdTdA, elev
       REAL*8,  DIMENSION(2, mxFEl), INTENT(IN) :: fDip
       INTEGER,                      INTENT(IN) :: mxEl, mxFEl, mxNode, n1000, nFakeN, nFl
       INTEGER, DIMENSION(4, mxFEl), INTENT(IN) :: nodeF
       INTEGER, DIMENSION(3, mxEl),  INTENT(IN) :: nodes
       INTEGER,                      INTENT(IN) :: nRealN, numEl, numNod
       REAL*8,  DIMENSION(mxFEl),    INTENT(IN) :: offset
       CHARACTER*100,                 INTENT(IN) :: title1
       REAL*8,  DIMENSION(mxNode),   INTENT(IN) :: tLNode, xNode, yNode, zMNode

       INTEGER :: i, iPrint, k, np(4)
       REAL*8  :: dips(2), pLat, pLon

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
            pLat = 90.0D0 - xNode(i) * degrees_per_radian
            pLon = yNode(i) * degrees_per_radian
            WRITE (iUnitO, 91) i, pLon, pLat, elev(i), dQdTdA(i), zMNode(i), tLNode(i)
   91       FORMAT (I8, 2F11.5, 4ES10.2)
  100  CONTINUE

       WRITE (iUnitO, 110) numEl
  110  FORMAT (I10,' (numEl = NUMBER OF TRIANGULAR CONTINUUM ELEMENTS)')

       DO 200 i = 1, numEl
            DO 150 k = 1, 3
                 IF (nodes(k, i) <= nRealN) THEN
                      np(k) = nodes(k, i)
                 ELSE
                      np(k) = n1000 + (nodes(k, i) - nRealN)
                 END IF
  150       CONTINUE
            WRITE (iUnitO, 160) i, (np(k), k = 1, 3)
  160       FORMAT (I8, 3I8)
  200  CONTINUE

       WRITE (iUnitO, 210) nFl
  210  FORMAT (I10,' (nFl =  NUMBER OF CURVILINEAR FAULT ELEMENTS)')

       DO 300 i = 1, nFl
            DO 220 k = 1, 4
                 IF (nodeF(k, i) <= nRealN) THEN
                      np(k) = nodeF(k, i)
                 ELSE
                      np(k) = n1000 + (nodeF(k, i) - nRealN)
                 END IF
  220       CONTINUE
            DO 230 k = 1, 2
                 dips(k) = fDip(k, i)
                 dips(k) = dips(k) * degrees_per_radian
                 IF (dips(k) > 90.01D0) dips(k) = dips(k) - 180.0D0
  230       CONTINUE
            WRITE (iUnitO, 250) i, (np(k), k = 1, 4), (dips(k), k = 1, 2), offset(i)
  250       FORMAT (I8, 4I6, 2F5.0, ES10.2)
  300  CONTINUE
       END SUBROUTINE PutNet

       REAL*8 FUNCTION ATan2F (y, x)

!   Corrects for problem of two zero arguments.
!   N.B. Order of these two arguments is NOT an error;
!        this order corresponds to ("motivator", "restrainer").

       REAL*8, INTENT(IN) :: y, x
       IF ((y /= 0.0D0).OR.(x /= 0.0D0)) THEN
            ATan2F = ATAN2(y, x)
       ELSE
            ATan2F = 0.0D0
       END IF
       END FUNCTION ATan2F


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
       INTEGER, INTENT(IN) :: iUnitT                                                         ! input
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
                    allOK = .FALSE.
	                write(ErrorMsg,'(A,I1,A,I5/,A,1P,E12.4,A,0P,F12.6)') "EXCESSIVELY DISTORTED ELEMENT LEADS TO NEGATIVE AREA AT POINT ",m," IN ELEMENT ",i ,&
		 &                     "AREA = ",area(i)," detJ: ",detJ(m, i)
		            call FatalError(ErrorMsg,ThID)
                 END IF
  610       CONTINUE
  620  CONTINUE


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
     &                 "BOUNDARY IS NOT SIMPLY-CONNECTED." , &
     &                 "Closed loop of ",nGood," nodes does not" , &
     &                 "include all ",nCond," boundary nodes." , &
     &                 "Run command PerimeterTest in -OrbWin-" , &
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
                      IF ( (ABS(xNode(i) - x) < 1.0D-6) .AND. &
     &                     (ABS(yNode(i) - y) < 1.0D-6) ) GO TO 867
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
!      End of indefinate loop which traces around perimeter.
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

       CALL FNodal (mxFEl, &                           ! input
     &              mxNode, nFl, nodeF, xNode, yNode, &
     &              fPFlt)                             ! output

       IF (.NOT. brief .AND. Verbose) WRITE(iUnitT,9999)
 9999  FORMAT (' -------------------------------------------------------------------------------')
       END SUBROUTINE Square


       SUBROUTINE ThetaPhi_2_plate (iUnitT, &
     &                              nPBnd, nBoundaryPoints, &
     &                              nPlate, pLat, pLon, &
     &                              theta, phi, & ! INTENT(IN)
     &                              iPlate)       ! INTENT(OUT)

!   Assigns an integer plate ID# to one surface point
!   with coordinates (theta = colatitude, in radians,
!                     phi = longitude, in radians).

!   NOTE: This subprogram should NOT be called to assign
!   F-E nodes to plates if those nodes might belong to
!   fault elements along plate edges!  Results would be purely
!   random, and would not consider topological connectivity
!   to an adjacent plate.  Instead, use subprogram Assign
!   from -Shells- (versions of 2006+).

       USE DSphere ! from Peter Bird's file DSphere.f90
       IMPLICIT NONE
       INTEGER, INTENT(IN)                         :: iUnitT
       INTEGER, INTENT(IN)                         :: nPBnd
       INTEGER, DIMENSION(nPlate), INTENT(IN)      :: nBoundaryPoints
       INTEGER, INTENT(IN)                         :: nPlate
       REAL*8, DIMENSION(nPlate, nPBnd), INTENT(IN):: pLat, pLon
       REAL*8, INTENT(IN)                          :: theta, phi
       INTEGER, INTENT(OUT)                        :: iPlate

       INTEGER :: iP, j, j2, nEnd, nPoint
       REAL*8 :: a1, a2, a3, aa, ab1, ab2, ab3, angle, ao, b1, b2, b3, bb, bo, dAngle, &
               & oxyz, sTheta, tangl, xo, xPoint, yo, yPoint, zo

!      PB2002 model of Bird [2003; G**3];
!      Already has plate "names" and "omega" vectors in
!      main program (DATA statements);
!      must also have digitised plate
!      outlines in arrays pLat and pLon,
!      presumably read from file "PB2002_plates.dig".
!      That is, this routine will not read any file.

            xo = COS(phi) * SIN(theta)
            yo = SIN(phi) * SIN(theta)
            zo = COS(theta)
            oxyz = xo * xo + yo * yo + zo * zo
            oxyz = SQRT(oxyz)
            xo = xo / oxyz
            yo = yo / oxyz
            zo = zo / oxyz
            nPoint = 0
            angle = 0.0D0
            iPlate = 0
            DO 500 iP = 1, nPlate
               tangl = 0.0D0
               nEnd = nBoundaryPoints(ip)
               DO 300 j = 1, nEnd
                  j2 = j + 1
                  IF(j == nEnd) THEN
                     j2 = 1
                  END IF
                  a1 = COS(pLon(ip, j)) * COS(pLat(ip, j))
                  a2 = SIN(pLon(ip, j)) * COS(pLat(ip, j))
                  a3 = SIN(pLat(ip, j))
                  b1 = COS(pLon(ip, j2)) * COS(pLat(ip, j2))
                  b2 = SIN(pLon(ip, j2)) * COS(pLat(ip, j2))
                  b3 = SIN(pLat(ip, j2))
                  ao = xo * a1 + yo * a2 + zo * a3
                  bo = xo * b1 + yo * b2 + zo * b3
                  a1 = a1 / ao
                  a2 = a2 / ao
                  a3 = a3 / ao
                  b1 = b1 / bo
                  b2 = b2 / bo
                  b3 = b3 / bo
                  a1 = a1 - xo
                  a2 = a2 - yo
                  a3 = a3 - zo
                  b1 = b1 - xo
                  b2 = b2 - yo
                  b3 = b3 - zo
                  aa = SQRT(a1 * a1 + a2 * a2 + a3 * a3)
                  bb = SQRT(b1 * b1 + b2 * b2 + b3 * b3)
                  ab1 = a2 * b3 - a3 * b2
                  ab2 = a3 * b1 - a1 * b3
                  ab3 = a1 * b2 - a2 * b1
                  sTheta = (ab1 * xo + ab2 * yo + ab3 * zo) / (aa * bb)
!                 prevent stupid abends due to imprecision:
                  sTheta = MAX(-1.0D0, MIN(1.0D0, stheta))
                  tangl = tangl + ASIN(stheta)
  300          CONTINUE
               dAngle = tangl - Pi
               IF(dAngle >= 0.0001D0) THEN
                  nPoint = nPoint + 1
                  iPlate = iP
               END IF
  500       CONTINUE
            IF(npoint >= 3) THEN
               xpoint = 90.0D0 - theta * degrees_per_radian
               ypoint = phi * degrees_per_radian
	           write(ErrorMsg,'(A,2F10.3,A)') "POINT ",xpoint, ypoint," WAS FOUND IN MORE THAN TWO PLATES; SOMETHING IS WRONG !!!!"
		       call FatalError(ErrorMsg,ThID)
            END IF
            IF (iPlate == 0) THEN
               xPoint = 90.0D0 - theta * degrees_per_radian
               yPoint = phi * degrees_per_radian
	           write(ErrorMsg,'(A,2F10.3,A)') "THE POINT ",xpoint, ypoint," DOES NOT BELONG TO ANY PLATE !!!!"
		       call FatalError(ErrorMsg,ThID)
            END IF
       END SUBROUTINE ThetaPhi_2_plate


       SUBROUTINE Tractor(iUnitQ, iUnitT, nPlates, nPoints, &
     &                    slab_q, phi_list, theta_list, whichP, &  ! INTENT(IN)
     &                    basal_shear_tractions)                   ! INTENT(OUT)

!      Requests file name of an existing torque report
!     (including traction pole vectors for each plate,
!      created by a previous experiment with -Shells-,
!      usually one that had trHMax = 0.0 and extra internal
!      velocity boundary conditions for each slabless plate).
!      Reads this file, extracts the traction pole vectors,
!      and uses them to precompute basal_shear_tractions shear tractions
!      at each point specified in (theta_list, phi_list) =
!      (colatitude, longitude) in radians.  These points must
!      already be identified as to plate affinity, in whichP.
!      For further clarification of "traction pole vectors"
!      see subprogram -TWIST- in SHELLS.

!      Derived from subprogram Tract of Shells (summer 2006),
!      but modified to free-format, and variable names changed
!      to drop previous implication that the list of points
!      to be processed is a list of node locations.

       USE DSphere ! in Peter Bird's file DSphere.f90
       IMPLICIT NONE
       INTEGER,                        INTENT(IN) :: iUnitQ, iUnitT, nPlates, nPoints
       LOGICAL, DIMENSION(nPlates),    INTENT(IN) :: slab_q
       REAL*8,  DIMENSION(nPoints),    INTENT(IN) :: theta_list, phi_list
       INTEGER, DIMENSION(nPoints),    INTENT(IN) :: whichP
       REAL*8,  DIMENSION(2, nPoints), INTENT(OUT):: basal_shear_tractions

       CHARACTER*132 :: line
       INTEGER :: i, ios, iPlate, j
       LOGICAL, DIMENSION(:), ALLOCATABLE :: tpread
       REAL*8 :: equat, lat, length, lon, t, tequat
       REAL*8, DIMENSION(3) :: tvec, uPhi, uTheta, uvec
       REAL*8, DIMENSION(:, :), ALLOCATABLE :: tpvecs

       ALLOCATE ( tpvecs(3, nPlates) )
       ALLOCATE ( tpread(nPlates) )
!      Zero whole array; advisable because some plates may not appear in the report.
       DO 30 j = 1, nPlates
            DO 20 i = 1, 3
                 tpvecs(i, j) = 0.0D0
   20       CONTINUE
            tpread(j) = .FALSE.
   30  CONTINUE

       IF(Verbose) WRITE(iUnitT,35) iUnitQ
   35  FORMAT(/' ATTEMPTING TO READ TORQUE REPORT ON UNIT', I3/)

!      Waste first 6 lines (titles & 2 blanks & header) of torque file.
       DO 40 i = 1, 6
            READ (iUnitQ, "(A)") line
   40  CONTINUE
!      Loop on plates in report (up to nPlates for whole-Earth model):
       DO 100 j = 1, nPlates
            READ(iUnitQ, * , IOSTAT = ios)
            IF (ios == -1) GO TO 101
            READ(iUnitQ, "(8X,I6)", IOSTAT = ios) iPlate
            IF (ios == -1) GO TO 101
!           Waste 23 more lines of each plate report
            DO 50 i = 1, 23
                 READ(iUnitQ, * , IOSTAT = ios)
                 IF (ios == -1) GO TO 101
   50       CONTINUE
            READ(iUnitQ, "(56X,ES10.3,2F10.2)") t, lon, lat
!           T is magnitude, in Pa, at location 90 deg. from (LON, LAT).
            tpvecs(1, iPlate) = t * COS(lat * radians_per_degree) * COS(lon * radians_per_degree)
            tpvecs(2, iPlate) = t * COS(lat * radians_per_degree) * SIN(lon * radians_per_degree)
            tpvecs(3, iPlate) = t * SIN(lat * radians_per_degree)
            tpread(iPlate) = .TRUE.
!           Waste 14 lines to get past the "=======" at the bottom of
!              each torque report:
            DO 60 i = 1, 14
                 READ(iUnitQ, * , IOSTAT = ios)
                 IF (ios == -1) GO TO 101
   60       CONTINUE
  100  CONTINUE
  101  CLOSE(iUnitQ)

       DO 200 i = 1, nPoints
            iPlate = whichP(i)
            IF (slab_q(iPlate)) THEN
!                no need for inferred basal_shear_tractions-strength traction:
                 basal_shear_tractions(1, i) = 0.0D0
                 basal_shear_tractions(2, i) = 0.0D0
            ELSE
                 IF (tpread(iPlate)) THEN

!                     Uvec is unit vector to node location:
                      uvec(1) = SIN(theta_list(i)) * COS(phi_list(i))
                      uvec(2) = SIN(theta_list(i)) * SIN(phi_list(i))
                      uvec(3) = COS(theta_list(i))

!                     Tvec is cross-product with traction pole vector:
                      tvec(1) = tpvecs(2, iPlate) * uvec(3) - &
     &                          tpvecs(3, iPlate) * uvec(2)
                      tvec(2) = tpvecs(3, iPlate) * uvec(1) - &
     &                          tpvecs(1, iPlate) * uvec(3)
                      tvec(3) = tpvecs(1, iPlate) * uvec(2) - &
     &                          tpvecs(2, iPlate) * uvec(1)
                      t = SQRT(tvec(1)**2 + tvec(2)**2 + tvec(3)**2)

!                     Unit vectors at this site (NOT a pole):
                      uPhi(1) = -uvec(2)
                      uPhi(2) = uvec(1)
                      equat = SIN(theta_list(i))
                      uPhi(1) = uPhi(1) / equat
                      uPhi(2) = uPhi(2) / equat
                      uPhi(3) = 0.0D0

                      tequat = uvec(3)
                      uTheta(3) = -equat
                      uTheta(1) = tequat * uvec(1) / equat
                      uTheta(2) = tequat * uvec(2) / equat
                      length = SQRT(uTheta(1)**2 + uTheta(2)**2 + uTheta(3)**2)
                      IF (length > 0.0D0) THEN
                          uTheta(1) = uTheta(1) / length
                          uTheta(2) = uTheta(2) / length
                          uTheta(3) = uTheta(3) / length
                      END IF

!                     Horizontal components of shear traction:
                      basal_shear_tractions(1, i) = tvec(1) * uTheta(1) + tvec(2) * uTheta(2) + &
     &                                              tvec(3) * uTheta(3)
                      basal_shear_tractions(2, i) = tvec(1) * uPhi(1) + tvec(2) * uPhi(2) + &
     &                                              tvec(3) * uPhi(3)
                 ELSE
                      basal_shear_tractions(1, i) = 0.0D0
                      basal_shear_tractions(2, i) = 0.0D0
                 END IF
            END IF
  200  CONTINUE
       DEALLOCATE ( tpread )
       DEALLOCATE ( tpvecs )
       END SUBROUTINE Tractor


       SUBROUTINE Unit (anyVec)

!   Converts any 3-component vector to a unit vector.

       REAL*8, DIMENSION(3), INTENT(INOUT) :: anyVec
       REAL*8 :: r2, size
       r2 = anyVec(1) * anyVec(1) + anyVec(2) * anyVec(2) + anyVec(3) * anyVec(3)
       IF (r2 > 0.0D0) THEN
            size = 1.0D0 / SQRT(r2)
            anyVec(1) = anyVec(1) * size
            anyVec(2) = anyVec(2) * size
            anyVec(3) = anyVec(3) * size
       ELSE
            anyVec(1) = 1.0D0
            anyVec(2) = 0.0D0
            anyVec(3) = 0.0D0
       END IF
       END SUBROUTINE Unit


end module
