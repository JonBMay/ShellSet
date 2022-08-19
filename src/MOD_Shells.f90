!*******************************************************************************
! Module containing subroutines used by Shells
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

module ShellsSubs

use ShellSetSubs
use SharedVars

! MKL version:
USE MKL95_PRECISION
USE MKL95_LAPACK
! Intel's Math Kernel Library (MKL), LAPACK portion; these MODULEs need INTERFACEs.
! These INTERFACEs are provided in file "lapack.f90", which must be available to be
!     compiled jointly with the rest of the project.
! I am using the following two routines from MKL:
!     dgbsv:  Simple driver to solve a REAL*8 linear system with banded coefficient matrix
!             in proprietary MKL "band storage scheme for LU factorization."
!     dgesv:  Simple driver to solve a small REAL*8 linear system.
! The advantage of this library is that I can compile for either 32-bit or
! 64-bit Windows, and with (or without) parallel processing,
! with little or no change to the source code.
! ***REMEMBER*** in the Microsoft Visual Studio GUI provided with
!                Intel Parallel Studio XE 2013, you need to set
!                Project / Properties / Fortran / Libraries /
!                   Use Math Kernel Library {to some choice other than "No"}!
!======================================================================================

implicit none

integer,parameter :: iUnitT = 6 ! Unit number for optional verbose.txt file

contains

SUBROUTINE AddFSt (constr, fC, fDip, fIMuDZ, fLen, fPFlt, fArg, & ! input
&                    mxFEl, &
&                    nFl, nodeF, &
&                    wedge, &
&                    k)                                             ! modify

!    Add fault stiffness to linear system.

!    A two-step process is used:
!    -A stiffness matrix for the fault element is formed, using
!     generic node numbering, 1-4.  Each entry in this matrix is
!     a 2 x 2 submatrix because node velocities have two components.
!    -The element stiffness matrix terms are added to the global
!     stiffness matrix.  (This step involves complex indirect
!     addressing, and is very difficult to optimize.)

!    The constant "constr" is the weight used in enforcing
!    strike-slip constraint equations.  It has the same units as
!    fIMuDZ and has a value comparable to
!    fMuMax*(thickness of the plate).

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: mxFEl, nFl, nodeF                                 ! input
REAL*8, INTENT(IN) :: constr, fC, fDip, fIMuDZ, fLen, fPFlt, fArg, wedge ! input
DOUBLE PRECISION, INTENT(INOUT) :: k                                     ! modify
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
DOUBLE PRECISION fPhi, fPoint, fGauss
COMMON / SFault / fPoint
COMMON / FPhis /  fPhi
COMMON / FGList / fGauss
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, i4, iFE, ifx, ify, iq, irx, iry, &
	 & j, j4, jcx, jcy, jfx, jfy, jj, &
	 & m, nodeI, nodeJ
REAL*8 ads11, ads12, ads21, ads22, angle, b11, b12, b21, b22, &
	& cosA, dip, dS, fP, &
	& oCosD, oSinD, oSin2D, sinA, &
	& tm1111, tm1112, tm1211, tm1212, tm2121, tm2122, tm2221, tm2222
DOUBLE PRECISION elK

DIMENSION fPhi(4, 7), fPoint(7), fGauss(7), fPFlt(2, 2, 2, 7, mxFEl)
DIMENSION fC(2, 2, 7, mxFEl), fDip(2, mxFEl), fLen(mxFEl), &
&           fIMuDZ(7, mxFEl), fArg(2, mxFEl), &
&           nodeF(4, mxFEl)
DIMENSION elK(2, 2, 4, 4), fP(2, 2, 4), k(nKRows, nRank)

!   Note: Convention is that row numbers identify the force balance
!         equation, while column numbers identify the degree
!         of freedom inFluencing that force balance.

DO 500 iFE = 1, nFl

!         Zero and then build up the element stiffness matrix:

  DO 10 i = 1, 4
	 DO 9 j = 1, 4
		elK(1, 1, i, j) = 0.0D0
		elK(1, 2, i, j) = 0.0D0
		elK(2, 1, i, j) = 0.0D0
		elK(2, 2, i, j) = 0.0D0
9        CONTINUE
10     CONTINUE

  DO 60 m = 1, 7
!CCCC         angle = fArg(1, iFE) * fPhi(1, m) + fArg(2, iFE) * fPhi(2, m)
!CCCC         Line above was replaced because of cycle-shift problem!

	  angle = Chord(fArg(1, ife), fPhi(2, m), fArg(2, ife))

	  sinA = SIN(angle)
	  cosA = COS(angle)
	  dS = fLen(ife) * fGauss(m)
	  dip = fPhi(1, m) * fDip(1, ife) + fPhi(2, m) * fDip(2, ife)
	  DO 15 j = 1, 2
		 jj = 4
		 IF (j == 2) jj = 3
		 fP(1, 1, j) = fPFlt(1, 1, j, m, iFE)
		 fP(2, 1, j) = fPFlt(2, 1, j, m, iFE)
		 fP(1, 2, j) = fPFlt(1, 2, j, m, iFE)
		 fP(2, 2, j) = fPFlt(2, 2, j, m, iFE)
		 fP(1, 1, jj) = -fP(1, 1, j)
		 fP(2, 1, jj) = -fP(2, 1, j)
		 fP(1, 2, jj) = -fP(1, 2, j)
		 fP(2, 2, jj) = -fP(2, 2, j)
15         CONTINUE
	  IF (ABS(dip - 1.57079632679490D0) <= wedge) THEN

!                Vertical strike-slip fault

		 DO 30 i = 1, 4
			DO 20 j = 1, 4
			   tm1111 = (fP(1, 1, j) * cosA + fP(1, 2, j) * sinA) &
&                        * (fP(1, 1, i) * cosA + fP(1, 2, i) * sinA) &
&                        * ds * fIMuDZ(m, ife)
			   tm1112 = (fP(1, 1, j) * cosA + fP(1, 2, j) * sinA) &
&                        * (fP(2, 1, i) * cosA + fP(2, 2, i) * sinA) &
&                        * ds * fIMuDZ(m, ife)
			   tm1211 = (fP(2, 1, j) * cosA + fP(2, 2, j) * sinA) &
&                        * (fP(1, 1, i) * cosA + fP(1, 2, i) * sinA) &
&                        * ds * fIMuDZ(m, ife)
			   tm1212 = (fP(2, 1, j) * cosA + fP(2, 2, j) * sinA) &
&                        * (fP(2, 1, i) * cosA + fP(2, 2, i) * sinA) &
&                        * ds * fIMuDZ(m, ife)
			   tm2121 = (-fP(1, 1, j) * sinA + fP(1, 2, j) * cosA) &
&                        * (-fP(1, 1, i) * sinA + fP(1, 2, i) * cosA) &
&                        * ds * constr
			   tm2122 = (-fP(1, 1, j) * sinA + fP(1, 2, j) * cosA) &
&                        * (-fP(2, 1, i) * sinA + fP(2, 2, i) * cosA) &
&                        * ds * constr
			   tm2221 = (-fP(2, 1, j) * sinA + fP(2, 2, j) * cosA) &
&                        * (-fP(1, 1, i) * sinA + fP(1, 2, i) * cosA) &
&                        * ds * constr
			   tm2222 = (-fP(2, 1, j) * sinA + fP(2, 2, j) * cosA) &
&                        * (-fP(2, 1, i) * sinA + fP(2, 2, i) * cosA) &
&                        * ds * constr
			   elK(1, 1, i, j) = elK(1, 1, i, j) + tm1111 &
&                                          + tm2121
			   elK(1, 2, i, j) = elK(1, 2, i, j) + tm1211 &
&                                          + tm2221
			   elK(2, 1, i, j) = elK(2, 1, i, j) + tm1112 &
&                                          + tm2122
			   elK(2, 2, i, j) = elK(2, 2, i, j) + tm1212 &
&                                          + tm2222
20               CONTINUE
30            CONTINUE
	  ELSE

!                Dipping oblique-slip fault

		 oSinD = 1.0D0 / SIN(dip)
		 oCosD = 1.0D0 / COS(dip)
		 oSin2D = 1.0D0 / SIN(2.0D0 * dip)
		 DO 50 i = 1, 4
			DO 40 j = 1, 4
			   ads11 = (fC(1, 1, m, ife) * ( fP(1, 1, j) * cosA + &
&                                        fP(1, 2, j) * sinA) &
&                 + fC(1, 2, m, ife) * oCosD * (-fP(1, 1, j) * sinA + &
&                                        fP(1, 2, j) * cosA)) * ds
			   ads12 = (fC(1, 1, m, ife) * ( fP(2, 1, j) * cosA + &
&                                        fP(2, 2, j) * sinA) &
&                 + fC(1, 2, m, ife) * oCosD * (-fP(2, 1, j) * sinA + &
&                                        fP(2, 2, j) * cosA)) * ds
			   ads21 = (fC(2, 1, m, ife) * ( fP(1, 1, j) * cosA + &
&                                        fP(1, 2, j) * sinA) &
&                 + fC(2, 2, m, ife) * oCosD * (-fP(1, 1, j) * sinA + &
&                                        fP(1, 2, j) * cosA)) * ds
			   ads22 = (fC(2, 1, m, ife) * ( fP(2, 1, j) * cosA + &
&                                        fP(2, 2, j) * sinA) &
&                 + fC(2, 2, m, ife) * oCosD * (-fP(2, 1, j) * sinA + &
&                                        fP(2, 2, j) * cosA)) * ds
			   b11 = fP(1, 1, i) * cosA + fP(1, 2, i) * sinA
			   b12 = fP(2, 1, i) * cosA + fP(2, 2, i) * sinA
			   b21 = -fP(1, 1, i) * sinA + fP(1, 2, i) * cosA
			   b22 = -fP(2, 1, i) * sinA + fP(2, 2, i) * cosA
			   elK(1, 1, i, j) = elK(1, 1, i, j) + ads11 * b11 * oSinD &
&                                        + 2.0D0 * ads21 * b21 * oSin2D
			   elK(1, 2, i, j) = elK(1, 2, i, j) + ads12 * b11 * oSinD &
&                                        + 2.0D0 * ads22 * b21 * oSin2D
			   elK(2, 1, i, j) = elK(2, 1, i, j) + ads11 * b12 * oSinD &
&                                        + 2.0D0 * ads21 * b22 * oSin2D
			   elK(2, 2, i, j) = elK(2, 2, i, j) + ads12 * b12 * oSinD &
&                                        + 2.0D0 * ads22 * b22 * oSin2D
40               CONTINUE
50            CONTINUE
	  END IF
60     CONTINUE

!         Apply element matrix to augment global stiffness matrix K:

	DO 400 i4 = 1, 4
	   nodeI = nodeF(i4, iFE)
	   iry = 2 * nodeI
	   irx = iry - 1
	   ify = 2 * i4
	   ifx = ify - 1
	   DO 300 j4 = 1, 4
		  nodeJ = nodeF(j4, iFE)
		  jcy = 2 * nodeJ
		  jcx = jcy - 1
		  jfy = 2 * j4
		  jfx = jfy - 1
		 !matrix element(irx, jcx):
		  iq = iDiagonal + irx - jcx
		  k(iq, jcx) = k(iq, jcx) + elK(1, 1, i4, j4)
		 !matrix element(irx, jcy):
		  iq = iDiagonal + irx - jcy
		  k(iq, jcy) = k(iq, jcy) + elK(1, 2, i4, j4)
		 !matrix element(iry, jcx):
		  iq = iDiagonal + iry - jcx
		  k(iq, jcx) = k(iq, jcx) + elK(2, 1, i4, j4)
		 !matrix element(iry, jcy):
		  iq = iDiagonal + iry - jcy
		  k(iq, jcy) = k(iq, jcy) + elK(2, 2, i4, j4)
300          CONTINUE
400       CONTINUE
500  CONTINUE
RETURN
END SUBROUTINE AddFSt

SUBROUTINE Assign (iUnitT, &           ! input
&                    nPBnd, nDPlat, nFl, nodeF, nodes, &
&                    nPlate, numEl, numNod, &
&                    pLat, pLon, &
&                    xNode, yNode, &
&                    whichP, &           ! output
&                    checkN)             ! work

!   Assigns an integer plate# to each node of the grid.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iUnitT, nPBnd, nFl, nPlate, numEl, numNod ! input
INTEGER, INTENT(IN) :: nDPlat                                    ! input
INTEGER, INTENT(IN) :: nodeF                                     ! input
INTEGER, INTENT(IN) :: nodes                                     ! input
REAL*8, INTENT(IN) :: pLat, pLon                                 ! input
REAL*8, INTENT(IN) :: xNode, yNode                               ! input
INTEGER, INTENT(OUT) :: whichP                                   ! output
LOGICAL checkN, inside                                           ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER :: i, iP, iPlate, j, j1, j2, k, n1, n2, n3, nEnd, nPoint, oldIP, outline_count
LOGICAL :: gotOut
REAL*8 aa, a1, a2, a3, ab1, ab2, ab3, ac, angle, ao, &
	& b1, b2, b3, bb, bc, bo, dangle, equat, &
	& g1, g2, g3, gc, length, length2, size, sTheta, tangl,  &
	& x, x1, x2, xInPl, xO, xOff, xPB, xPoint, xVel, &
	& y, y1, y2, yInPl, yO, yOff, yPB, yPoint, yVel, &
	& z1, z2, zO, zOff, zPB
REAL*8, DIMENSION(3) :: alongv, crossv, uvec, uvec1
DIMENSION nDPlat(nPlate), nodeF(4, nFl), nodes(3, numEl), &
	   & pLat(nPlate, nPBnd), pLon(nPlate, nPBnd), &
	   & xNode(numNod), yNode(numNod), whichP(numNod), checkN(numNod)
REAL*8, DIMENSION(:, :), ALLOCATABLE :: plate_outline_uvecs

!      PB2002SCEC model of Bird [2003; G**3] + [2017.01 microplates];
!      Already has plate "names" and "omega" vectors in
!      main program (DATA statements);
!      must also have digitised plate
!      outlines in arrays pLat and pLon,
!      presumably already read from file "PB2002SCEC_plates.dig".
!     (That is, this routine will not read any file.)

ALLOCATE ( plate_outline_uvecs(3, nPBnd) ) ! using global/MAIN value of anticipated maximum length for any plate circuit

!      Check which nodes are on faults:
DO 10 i = 1, numNod
	checkN(i) = .FALSE.
10  CONTINUE
DO 30 i = 1, nFl
	DO 20 k = 1, 4
		 checkN(nodeF(k, i)) = .TRUE.
20       CONTINUE
30  CONTINUE

!      For nodes on faults, attempt to offset test position
!      which is used to determine plate membership
!      (but not position used in V = OMEGA x R ).
oldIP = 1 ! to guard against possible undefined integer, in case of failure on the first point
DO 999 i = 1, numNod
	xVel = xNode(i)
	yVel = yNode(i)
	IF (checkN(i)) THEN

!                Node is on fault; seek offset position
!                for determination of plate affiliation...

		 gotOut = .FALSE.

!                1st strategy:
!                Is there a continuum element including this node
!                which has some other node NOT on a fault?
!                If so, use that other node's position.

		 DO 100 j = 1, numEl
			  n1 = nodes(1, j)
			  n2 = nodes(2, j)
			  n3 = nodes(3, j)
			  IF ((n1 == i).OR.(n2 == i).OR.(n3 == i)) THEN
				   IF ((n1 /= i).AND.(.NOT.checkN(n1)))THEN
						gotOut = .TRUE.
						xInPl = xNode(n1)
						yInPl = yNode(n1)
						GO TO 101
				   END IF
				   IF ((n2 /= i).AND.(.NOT.checkN(n2)))THEN
						gotOut = .TRUE.
						xInPl = xNode(n2)
						yInPl = yNode(n2)
						GO TO 101
				   END IF
				   IF ((n3 /= i).AND.(.NOT.checkN(n3)))THEN
						gotOut = .TRUE.
						xInPl = xNode(n3)
						yInPl = yNode(n3)
						GO TO 101
				   END IF
			  END IF
100            CONTINUE

!                If there is still a problem, try
!                2nd strategy:
!                If any continuum element includes this node
!               (even though its other nodes are all on faults),
!                we can use the midpoint of the continuum element...

101            IF (.NOT.gotOut) THEN
			  DO 200 j = 1, numEl
				   n1 = nodes(1, j)
				   n2 = nodes(2, j)
				   n3 = nodes(3, j)
				   IF ((n1 == i).OR.(n2 == i).OR. &
&                         (n3 == i)) THEN
						gotOut = .TRUE.
						a1 = SIN(xNode(n1)) * COS(yNode(n1))
						b1 = SIN(xNode(n1)) * SIN(yNode(n1))
						g1 = COS(xNode(n1))
						a2 = SIN(xNode(n2)) * COS(yNode(n2))
						b2 = SIN(xNode(n2)) * SIN(yNode(n2))
						g2 = COS(xNode(n2))
						a3 = SIN(xNode(n3)) * COS(yNode(n3))
						b3 = SIN(xNode(n3)) * SIN(yNode(n3))
						g3 = COS(xNode(n3))
						ac = (a1 + a2 + a3) / 3.0D0
						bc = (b1 + b2 + b3) / 3.0D0
						gc = (g1 + g2 + g3) / 3.0D0
						size = SQRT(ac**2 + bc**2 + gc**2)
						ac = ac / size
						bc = bc / size
						gc = gc / size
						equat = SQRT(ac**2 + bc**2)
						xInPl = ATAN2(equat, gc)
						yInPl = ATAN2(bc, ac)
						GO TO 201
				   END IF
200                 CONTINUE
		 END IF

!                If there is still a problem, then this fault
!                node does not belong to any continuum element.
!                It must be on the outer perimeter of the model.
!                Try a small offset toward the outside...

201            IF (.NOT.gotOut) THEN
!                     Find where node #i is on the fault...
			  DO 220 j = 1, nFl
				   DO 210 k = 1, 4
						IF (nodeF(k, j) == i) THEN
!                                    N.B. k & j are what we are seeking.
							 GO TO 221
						END IF
210                      CONTINUE
220                 CONTINUE
221                 IF (k <= 2) THEN
!                          Node is on N1-N2 side of fault.
				   n1 = nodeF(1, j)
				   n2 = nodeF(2, j)
			  ELSE
!                          Node is on N3-N4 side of fault.
				   n1 = nodeF(3, j)
				   n2 = nodeF(4, j)
			  END IF

			  x1 = COS(yNode(n1)) * SIN(xNode(n1))
			  y1 = SIN(yNode(n1)) * SIN(xNode(n1))
			  z1 = COS(xNode(n1))
			  x2 = COS(yNode(n2)) * SIN(xNode(n2))
			  y2 = SIN(yNode(n2)) * SIN(xNode(n2))
			  z2 = COS(xNode(n2))
			  alongV(1) = x2 - x1
			  alongV(2) = y2 - y1
			  alongV(3) = z2 - z1
			  xOff = x1 + 0.50D0 * alongV(1)
			  yOff = y1 + 0.50D0 * alongV(2)
			  zOff = z1 + 0.50D0 * alongV(3)
			  crossV(1) = alongV(2) * zOff - alongV(3) * yOff
			  crossV(2) = alongV(3) * xOff - alongV(1) * zOff
			  crossV(3) = alongV(1) * yOff - alongV(2) * xOff
!                    "crossV: has same length as alongV,
!                     and points out of fault (to right,
!                     when looking from n1 toward n2).
			  xOff = xOff + 0.250D0 * crossV(1)
			  yOff = yOff + 0.250D0 * crossV(2)
			  zOff = zOff + 0.250D0 * crossV(3)
			  equat = SQRT(xOff**2 + yOff**2)
			  xInPl = ATAN2(equat, zOff)
			  yInPl = ATAN2(yOff, xOff)
		 END IF
	ELSE

!                Node is not on any fault;
!                no offset of position is needed:

		 xInPl = xVel
		 yInPl = yVel
	END IF
	!Convert test position to a spherical uvec:
	xO = COS(yInPl) * SIN(xInPl)
	yO = SIN(yInPl) * SIN(xInPl)
	zO = COS(xInPl)
	!including a normalization step that (in theory) should not be necessary:
	length2 = (xO * xO) + (yO * yO) + (zO * zO)
	length = SQRT(length2)
	xO = xO / length
	yO = yO / length
	zO = zO / length
	uvec(1) = xO; uvec(2) = yO; uvec(3) = zO
	!Initialize search(es) for plates enclosing this test point:
	nPoint = 0 ! number of plates enclosing this test point
	iPlate = 0 ! INTEGER index of (last) plate enclosing this test point
	DO 500 iP = 1, nPlate
	   outline_count = nDPlat(iP)
	   DO j = 1, outline_count
		   !Convert plate-boundary positions to spherical uvecs, all around this one plate:
		   !{N.B. Arrays pLat and pLon have already been divided by 57.296... to convert them to radians!}
		   x = 1.57079632679490D0 - pLat(iP, j)  ! colatitude, measured from N pole, in radians
		   y = pLon(iP, j)                       ! East longitude, in radians
		   xPB = COS(y) * SIN(x)
		   yPB = SIN(y) * SIN(x)
		   zPB = COS(x)
		   !including a normalization step that (in theory) should not be necessary:
		   length2 = (xPB * xPB) + (yPB * yPB) + (zPB * zPB)
		   length = SQRT(length2)
		   xPB = xPB / length
		   yPB = yPB / length
		   zPB = zPB / length
		   uvec1(1) = xPB; uvec1(2) = yPB; uvec1(3) = zPB
		   !and store in array (needed by function Within):
		   plate_outline_uvecs(1:3, j) = uvec1(1:3)
	   END DO
	   !=================================================================
	   inside = Within(uvec, outline_count, plate_outline_uvecs)
	   !=================================================================
	   IF(inside) THEN
		  nPoint = nPoint + 1
		  iPlate = iP
	   END IF
500       CONTINUE
	IF (iPlate == 0) THEN
	   xPoint = 90.0D0 - (xInPl * 57.2957795130823D0)
	   yPoint = yInPl * 57.2957795130823D0
	   WRITE(iUnitT, 600) xPoint, yPoint
	   IF(Verbose) WRITE(iUnitVerb, 600) xPoint, yPoint, oldIP
600          FORMAT(/' THE POINT (',F13.5,'N, ',F13.5'E) DOES NOT BELONG TO ANY PLATE !!!!' &
&               /' Arbitrarily Assign-ing to last previous plate (#',I4,').')
	   iPlate = oldIP
	   IF(Verbose) WRITE(iUnitVerb, "(' Continuing to assign all nodes to plates...')")
	END IF
	IF(nPoint > 4) THEN
	   xPoint = 90.0D0 - (xInPl * 57.2957795130823D0)
	   yPoint = yInPl * 57.2957795130823D0
       write(ErrorMsg,'(A,F13.5,A,F13.5,A)') "THE POINT (",xPoint,"N, ",yPoint,"E) WAS FOUND IN MORE THAN FOUR PLATES; SOMETHING IS WRONG!"
       call FatalError(ErrorMsg,ThID)
	END IF

	whichP(i) = iPlate
	oldIP = iPlate

999  CONTINUE
END SUBROUTINE Assign

REAL*8 FUNCTION ATan2F (y, x)
!   Corrects for problem of abend due to ATAN2(0.0D0, 0.0D0):
IMPLICIT NONE
!      - - - - - - - - - - -
REAL*8, INTENT(IN) :: x, y   ! input
!      - - - - - - - - - - -
IF ((y /= 0.0D0).OR.(x /= 0.0D0)) THEN
	ATan2F = ATAN2(y, x)
ELSE
	ATan2F = 0.0D0
END IF
RETURN
END FUNCTION ATan2F

SUBROUTINE Balanc (alphaT, area, conduc, constr, &         ! input
&                    density_anomaly, detJ, dQdTdA, dXS, &
&                    dXSP, dYS, dYSP, edgeTS, elev, eta, &
&                    fArg, fC, fDip, &
&                    fIMuDZ, fLen, fPFlt, fPSfer, fTStar, &
&                    gMean, iCond, iUnitF, &
&                    iUnitT, log_force_balance, &
&                    mxBn, mxDOF, mxEl, mxFEl, mxNode, &
&                    nCond, nFl, nodCon, nodeF, nodes, &
&                    numEl, numNod, oneKm, oVB, radio, radius, &
&                    rhoAst, rhoBar, rhoH2O, &
&                    sigZZI, sita, &
&                    tauMat, tauZZI, tauZZN, temLim, &
&                    title1, title2, title3, &
&                    tLNode, tSurf, v, wedge, xNode, yNode, zMNode, &
&                    sigHB, &                               ! modify
&                    comp, &                                ! output
&                    fBase, outVec)                         ! work

!=============
! Test #1:
!=============
!   Checks the balance of forces on each node by computing the
!   convolutions of all nodal functions with the standardised
!   traction anomalies over the surface of artificial fixed-V
!   boundaries.  For any degree of freedom that has been removed by
!   the imposition of a velocity boundary condition, this gives the
!   consistent nodal force required to impose that velocity.
!   For boundary nodes which are free (except for lithostatic
!   pressure), gives the horizontal force due to this pressure.
!   If Shells is working properly, then the convolutions
!   for internal nodes should be "small" (compared to the area
!   integral of the dot product of a nodal function (of order 1)
!   with typical traction anomalies over the surface
!   which projects vertically through the plate from a typical
!   element side on the surface).  For example, in SI units, if the
!   layer thickness is 1D5 m, the typical element side is
!   6D5 m, and the typical stress anomaly is 5.D7 Pa, then
!   an apparent force of 3.D18 N would be "large" (100% error),
!   but an apparent force OF 3.D15 N would be "small" (equilibrium
!   within 0.1%).

!=============
! Test #2:
!=============
!   Checks the global balance of forces by computing the net
!   torque on the model passed through all boundaries
!   (basal pressure and drag, faults, side boundary conditions)
!   which should, of course, be zero (except for numeric error).

!   Note: For this purpose, the array sigHB is modified from
!   the nonlinear flow-law value (of -THOnB-) to the linearized
!   value, based on the ETA of the last iteration.  If the
!   solution is well converged, this should be a small change.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alphaT, area, conduc, constr, density_anomaly, detJ, &      ! input
				   & dQdTdA, dXS, dXSP, dYS, dYSP                                ! input
LOGICAL, INTENT(IN) :: edgeTS                                                     ! input
REAL*8, INTENT(IN) :: elev, eta, &                                                ! input
				   & fArg, fC, fDip, fIMuDZ, fLen, fPFlt, fPSfer, fTStar, gMean  ! input
INTEGER, INTENT(IN) :: iCond, iUnitF, iUnitT, mxBn, mxDOF, mxEl, mxFEl, mxNode, & ! input
					& nCond, nFl, nodCon, nodeF, nodes, numEl, numNod            ! input
LOGICAL, INTENT(IN) :: log_force_balance                                          ! input
REAL*8, INTENT(IN) :: oneKm, oVB, radio, radius, rhoAst, rhoBar, rhoH2O, &        ! input
				   & sigZZI, sita, &                                             ! input
				   & tauMat, tauZZI, tauZZN, temLim                              ! input
CHARACTER*100, INTENT(IN) :: title1, title2, title3                                ! input
REAL*8, INTENT(IN) :: tLNode, tSurf                                               ! input
DOUBLE PRECISION, INTENT(IN) :: v                                                 ! input
REAL*8, INTENT(IN) :: wedge, xNode, yNode, zMNode                                 ! input
REAL*8, INTENT(INOUT) :: sigHB                                                    ! modify
REAL*8, INTENT(OUT) :: comp                                                       ! output
DOUBLE PRECISION fBase                                             ! work
REAL*8 outVec                                                      ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DOUBLE PRECISION points, weight
DOUBLE PRECISION fPhi, fPoint, fGauss
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CHARACTER*2 cx, cy, cz
CHARACTER*15 large, size, small
INTEGER i, ic, ix, iy, j, jj, k, m, n1, n2, n3, n4, nEntry, node
LOGICAL doFB1, doFB2, doFB3, doFB4, hasBC, sloped
DOUBLE PRECISION tQxx, tQyy, tQzz
REAL*8 angle, dAOR, close, &
	& dA, delVx, delVy, dip, dS, dUxx, dUxy, dUyy, &
	& fx, fy, fz, phai, rQxx, rQyy, rQzz, &
	& sinist, sinS, sumvdf, tanS, taidb, taidz, tauxx, tauxy, tauyy, &
	& tbidz, tFxx, tFyy, tFzz, tgidb, tgidz, theta, txidb, tyidb, &
	& unitAx, unitAy, unitBx, unitBy, vBCa, vUpDip

COMMON / SFault / fPoint
COMMON / FPhis /  fPhi
COMMON / FGList / fGauss
COMMON / S1S2S3 / points
COMMON / WgtVec / weight

DIMENSION points(3, 7), weight(7)
DIMENSION fPhi(4, 7), fPoint(7), fGauss(7)
DIMENSION area(mxEl), alphaT(2), &
&           comp(6, mxDOF), conduc(2), density_anomaly(mxNode), &
&           detJ(7, mxEl), dQdTdA(mxNode), &
&           dXS(2, 2, 3, 7, mxEl), dXSP(3, 7, mxEl), &
&           dYS(2, 2, 3, 7, mxEl), dYSP(3, 7, mxEl), &
&           edgeTS(3, mxEl), elev(mxNode), eta(7, mxEl), &
&           fArg(2, mxFEl), fBase(mxDOF), fC(2, 2, 7, mxFEl), &
&           fDip(2, mxFEl), fIMuDZ(7, mxFEl), &
&           fLen(mxFEl), fPFlt(2, 2, 2, 7, mxFEl), &
&           fPSfer(2, 2, 3, 7, mxEl), fTStar(2, 7, mxFEl), iCond(mxBn), &
&           nodCon(mxBn), nodeF(4, mxFEl), nodes(3, mxEl), &
&           outVec(2, 7, mxEl), oVB(2, 7, mxEl), radio(2), rhoBar(2), &
&           sigHB(2, 7, mxEl), sigZZI(7, mxEl), sita(7, mxEl), &
&           tauMat(3, 7, mxEl), tauZZI(7, mxEl), tauZZN(mxNode), &
&           temLim(2), tLNode(mxNode), &
&           v(2, mxNode), &
&           xNode(mxNode), yNode(mxNode), zMNode(mxNode)
DATA cx / '+S' / , cy / '+E' /
DATA large / 'MAY BE LARGE   ' / , small / 'should be small' /

!   Write text explaining purpose of table:

IF (log_force_balance .AND. Verbose) WRITE (iUnitVerb, 10)
10  FORMAT (' ===================================', &
&          '===================================')
IF (log_force_balance .AND. Verbose) WRITE (iUnitVerb, 1)
1  FORMAT (/ &
&/' Check the balance of forces on each node by computing the' &
&/' convolutions of all nodal functions with the standardised' &
&/' traction anomalies over the surface of artificial fixed-V' &
&/' boundaries.  For degrees of freedom that have been removed by' &
&/' the imposition of velocity boundary conditions, this gives the' &
&/' consistent nodal forces required to impose those velocities.' &
&/' For boundary nodes which are free (except for lithostatic' &
&/' pressure), gives the horizontal force due to this pressure.')
IF (log_force_balance .AND. Verbose) WRITE (iUnitVerb, 2)
2  FORMAT (/ &
&/' If the program is working correctly, then the convolutions' &
&/' for internal nodes should be "small" (compared to the area' &
&/' integral of the dot product of a nodal function (order 1)' &
&/' with typical traction anomalies over the surface' &
&/' which projects vertically through the plate from a typical' &
&/' element side on the surface).  For example, in SI units, if' &
&/' the layer thickness is 1D5 m, the typical element side is' &
&/' 6D5 m, and the typical stress anomaly is 5.D7 Pa, then' &
&/' an apparent force of 3.D18 N would be "large" (100% error),' &
&/' but an apparent force of 3.D15 would be "small" (equilibrium' &
&/' within 0.1%).')
IF (log_force_balance .AND. Verbose) WRITE (iUnitVerb, 3)
3  FORMAT (/' Explanation of the table:'/ / &
&' Each row corresponds to one degree of freedom, so'/ &
&'   row 1 gives apparent force on node 1 in the South direction,'/ &
&'   row 2 gives apparent force on node 1 in the East direction,'/ &
&'   row 3 gives apparent force on node 2 in the South direction,'/ &
&'     et cetera.'/ / &
&' The *Total* column gives the convolution of nodal functions'/ &
&'    with all traction anomalies on all surfaces of the model'/ &
&'   (external top, bottom, and sides, plus internal faults).'/ &
&' The *Fault_P*  column gives the convolution of nodal functions'/ &
&'    with pressure-anomalies on internal fault surfaces.'/ &
&' The *Fault_S*  column gives the convolution of nodal functions'/ &
&'    with deviatoric tractions on internal fault surfaces.'/ &
&' The *Base_P*  column gives the convolution of nodal functions'/ &
&'    with pressure-anomaly*grad(depth of base of lithosphere).'/ &
&' The *Base_S*  column gives the convolution of nodal functions'/ &
&'    with shear tractions on the base of the lithosphere.'/ &
&' The *Bounds*  column = *Total* - *Fault* - *Base_P* - *Base_S*,' &
&/'    and gives the convolution of nodal functions with'/ &
&'    traction anomalies on artificial fixed-V nodes.'/ &
&' The *Comment* column indicates whether *Bounds* is expected'/ &
&'    to be small or not.'/ /)
IF (log_force_balance .AND. Verbose) WRITE (iUnitVerb, 4)
4  FORMAT ('       Node Component   *Total* *Fault_P*', &
&         ' *Fault_S*  *Base_P*  *Base_S*  *Bounds* *Comment*')
5  FORMAT (' ',I10,8X,A2,6ES10.2,1X,A15)

nEntry = 2 * numNod
DO 110 i = 1, nEntry
	comp(1, i) = 0.0D0
	comp(2, i) = 0.0D0
	comp(3, i) = 0.0D0
	comp(4, i) = 0.0D0
	comp(5, i) = 0.0D0
!          (no need to zero comp(6, i), as it is not an accumulated sum)
110  CONTINUE

!  *Total* column is convolution of nodal functions with
!          total inhomogeneous term(s) in horizontal
!          stress-equilibrium equation for the plate shell(s):

DO 150 m = 1, 7
	DO 140 i = 1, numEl
		 tauxx = tauMat(1, m, i) + tauZZI(m, i)
		 tauyy = tauMat(2, m, i) + tauZZI(m, i)
		 tauxy = tauMat(3, m, i)
		 sinS = 1.00D0 / SIN(sita(m, i))
		 tanS = 1.00D0 / TAN(sita(m, i))
		 dAOR = area(i) * weight(m) * detJ(m, i) / radius
		 DO 130 j = 1, 3
			  node = nodes(j, i)
			  ix = 2 * node - 1
			  iy = ix + 1
			  dUxx = tauxx * dXS(1, 1, j, m, i)
			  dUxy = tauxy * (dXS(1, 2, j, m, i) &
&                          + dYS(1, 1, j, m, i) * sinS &
&                          - fPSfer(1, 2, j, m, i) * tanS)
			  dUyy = tauyy * (dYS(1, 2, j, m, i) * sinS &
&                          + fPSfer(1, 1, j, m, i) * tanS)
			  comp(1, ix) = comp(1, ix) + dAOR * (dUxx + dUxy + dUyy)
			  dUxx = tauxx * dXS(2, 1, j, m, i)
			  dUxy = tauxy * (dXS(2, 2, j, m, i) &
&                          + dYS(2, 1, j, m, i) * sinS &
&                          - fPSfer(2, 2, j, m, i) * tanS)
			  dUyy = tauyy * (dYS(2, 2, j, m, i) * sinS &
&                          + fPSfer(2, 1, j, m, i) * tanS)
			  comp(1, iy) = comp(1, iy) + dAOR * (dUxx + dUxy + dUyy)
130            CONTINUE
140       CONTINUE
150  CONTINUE

!   *Fault_P* column is convolution of nodal functions with
!           lithostatic pressures on fault surfaces:
!           Anomaly in lithostatic vertical compressive stress
!          (relative to standard weak-ridge structure)
!           integrated down the dip of the fault:

doFB1 = .FALSE.
doFB2 = .FALSE.
doFB3 = .FALSE.
doFB4 = .TRUE.
CALL Fixed (alphaT, area, conduc, &  ! input
&             density_anomaly, detJ, &
&             doFB1, doFB2, doFB3, doFB4, &
&             dQdTdA, dXS, dYS, &
&             dXSP, dYSP, edgeTS, elev, fDip, fLen, fPFlt, &
&             fPSfer, fArg, gMean, &
&             iCond, iUnitT, &
&             mxBn, mxDOF, mxEl, mxFEl, mxNode, &
&             nCond, nFl, nodCon, nodeF, nodes, numEl, &
&             oneKm, radio, radius, &
&             rhoAst, rhoBar, rhoH2O, sigZZI, &
&             sita, tauZZI, tauZZN, temLim, tLNode, tSurf, wedge, &
&             xNode, yNode, zMNode, &
&             fBase)                   ! output
DO 210 i = 1, nEntry
	comp(2, i) = fBase(i)
210  CONTINUE

!     *Fault_S* column:

!      tractions caused by departures of stress from the local
!      lithostatic pressure, due to strength (deviatoric stress),
!      integrated down the dip of faults:

DO 290 i = 1, nFl
	n1 = nodeF(1, i)
	n2 = nodeF(2, i)
	n3 = nodeF(3, i)
	n4 = nodeF(4, i)
	DO 280 m = 1, 7

!               "angle" is the fault strike, in radians cclkws from +X.
!CCCC            angle = fArg(1, i) * fPhi(1, m) + fArg(2, i) * fPhi(2, m)
!CCCC            Line above was replaced due to cycle-shift problem!

		 angle = Chord(fArg(1, i), fPhi(2, m), fArg(2, i))
!               "angle" is the argument of the forward ray from n1-->n2,
!                measured counterclockwise from +theta = +X = South.

!               "unitA" is a unit vector along the fault, from N1-->N2.
		 unitAx = COS(angle)
		 unitAy = SIN(angle)

!               "unitB" is a perpendicular unit vector, pointing in
!               (toward the n4-n3 side).
		 unitBx = -unitAy
		 unitBy = +unitAx

!                Relative velocities are for n1-n2 side relative to
!                the n4-n3 side:
		 delVx = v(1, n1) * fPFlt(1, 1, 1, m, i) + v(2, n1) * fPFlt(2, 1, 1, m, i) &
&                 + v(1, n2) * fPFlt(1, 1, 2, m, i) + v(2, n2) * fPFlt(2, 1, 2, m, i) &
&                 - v(1, n3) * fPFlt(1, 1, 2, m, i) - v(2, n3) * fPFlt(2, 1, 2, m, i) &
&                 - v(1, n4) * fPFlt(1, 1, 1, m, i) - v(2, n4) * fPFlt(2, 1, 1, m, i)
		 delVy = v(1, n1) * fPFlt(1, 2, 1, m, i) + v(2, n1) * fPFlt(2, 2, 1, m, i) &
&                 + v(1, n2) * fPFlt(1, 2, 2, m, i) + v(2, n2) * fPFlt(2, 2, 2, m, i) &
&                 - v(1, n3) * fPFlt(1, 2, 2, m, i) - v(2, n3) * fPFlt(2, 2, 2, m, i) &
&                 - v(1, n4) * fPFlt(1, 2, 1, m, i) - v(2, n4) * fPFlt(2, 2, 1, m, i)

!                Sinistral strike-slip rate component:
		 sinist = delVx * unitAx + delVy * unitAy

!                Convergence rate component (in horizontal plane):
		 close = delVx * unitBx + delVy * unitBy

!                Dip of the fault (from horizontal on the n1-n2 side):
		 dip = fDip(1, i) * fPhi(1, m) + fDip(2, i) * fPhi(2, m)
		 sloped = ABS(dip - 1.57079632679490D0) > wedge

!                Find traction on n3-n4 side in (alpha, gamma) coordinates
		 IF (sloped) THEN
			  vUpDip = close / COS(dip)
!                     positive for thrusting with dip < Pi/2
			  taidz = fC(1, 1, m, i) * sinist + fC(1, 2, m, i) * vUpDip + &
&                        fTStar(1, m, i)
			  tbidz = fC(2, 1, m, i) * sinist + fC(2, 2, m, i) * vUpDip + &
&                        fTStar(2, m, i)
			  tgidz = tbidz / COS(dip)
			  taidb = taidz / SIN(dip)
			  tgidb = tgidz / SIN(dip)
		 ELSE
			  taidb = fIMuDZ(m, i) * sinist
			  tgidb = constr * close
		 END IF

!                Reverse for tractions on n1-n2 side:
		 taidb = -taidb
		 tgidb = -tgidb
!                Now, positive tgidb is associated with divergence
!                and positive taidb is associated with dextral slip.

!                Express traction on n1-n2 side in (x, y) coordinates:
		 txidb = taidb * COS(angle) - tgidb * SIN(angle)
		 tyidb = tgidb * COS(angle) + taidb * SIN(angle)

		 dS = fLen(i) * fGauss(m)

		 DO 270 j = 1, 2
			  node = nodeF(j, i)
			  ix = 2 * node - 1
			  iy = ix + 1
			  comp(3, ix) = comp(3, ix) + ds * &
&                   (txidb * fPFlt(1, 1, j, m, i) + tyidb * fPFlt(1, 2, j, m, i))
			  comp(3, iy) = comp(3, iy) + ds * &
&                   (txidb * fPFlt(2, 1, j, m, i) + tyidb * fPFlt(2, 2, j, m, i))
			  jj = 5 - j
			  node = nodeF(jj, i)
			  ix = 2 * node - 1
			  iy = ix + 1
			  comp(3, ix) = comp(3, ix) - ds * &
&                   (txidb * fPFlt(1, 1, j, m, i) + tyidb * fPFlt(1, 2, j, m, i))
			  comp(3, iy) = comp(3, iy) - ds * &
&                   (txidb * fPFlt(2, 1, j, m, i) + tyidb * fPFlt(2, 2, j, m, i))
270            CONTINUE
280       CONTINUE
290  CONTINUE

!  *Base_P* column is convolution of nodal functions with
!      basal pressure-anomaly * grad(bottom depth):

doFB1 = .FALSE.
doFB2 = .TRUE.
doFB3 = .FALSE.
doFB4 = .FALSE.
CALL Fixed (alphaT, area, conduc, &   ! input
&             density_anomaly, detJ, &
&             doFB1, doFB2, doFB3, doFB4, &
&             dQdTdA, dXS, dYS, &
&             dXSP, dYSP, edgeTS, elev, fDip, fLen, fPFlt, &
&             fPSfer, fArg, gMean, &
&             iCond, iUnitT, &
&             mxBn, mxDOF, mxEl, mxFEl, mxNode, &
&             nCond, nFl, nodCon, nodeF, nodes, numEl, &
&             oneKm, radio, radius, &
&             rhoAst, rhoBar, rhoH2O, sigZZI, &
&             sita, tauZZI, tauZZN, temLim, tLNode, tSurf, wedge, &
&             xNode, yNode, zMNode, &
&             fBase)                   ! output
DO 310 i = 1, nEntry
	comp(4, i) = fBase(i)
310  CONTINUE

!  *Base_S* column is convolution of nodal functions with
!      shear tractions on base of lithosphere
!     (based on linearized form, not flow-law form;
!      should be similar if solution is converged):

CALL Flow (fPSfer, mxEl, mxNode, nodes, numEl, v, & ! input
&            outVec)                                  ! output
DO 350 m = 1, 7
	DO 340 i = 1, numEl
		 dA = area(i) * weight(m) * detJ(m, i)
		 sigHB(1, m, i) = eta(m, i) * (oVB(1, m, i) - outVec(1, m, i))
		 sigHB(2, m, i) = eta(m, i) * (oVB(2, m, i) - outVec(2, m, i))
		 DO 330 j = 1, 3
			  node = nodes(j, i)
			  ix = 2 * node - 1
			  iy = ix + 1
			  comp(5, ix) = comp(5, ix) + dA * &
&                    (sigHB(1, m, i) * fPSfer(1, 1, j, m, i) + &
&                     sigHB(2, m, i) * fPSfer(1, 2, j, m, i))
			  comp(5, iy) = comp(5, iy) + dA * &
&                    (sigHB(1, m, i) * fPSfer(2, 1, j, m, i) + &
&                     sigHB(2, m, i) * fPSfer(2, 2, j, m, i))
330            CONTINUE
340       CONTINUE
350  CONTINUE

!  *Bounds* column is inferred from preceding five columns,
!         and should represent the either:
!     (1) consistent nodal forces due to velocity boundary conditions;
!  OR (2) consistent nodal forces due to pressure anomalies on any
!        "free" lateral boundaries [not applicable in global models].

DO 400 i = 1, nEntry
	comp(6, i) = comp(1, i) - comp(2, i) - comp(3, i) - comp(4, i) - comp(5, i)
400  CONTINUE

!   Write out matrix, with annotations:
!  (Also, sum forces times velocities over boundary nodes.)

sumvdf = 0.0D0
DO 1000 i = 1, nEntry

	node = (i + 1) / 2

	IF (MOD(i, 2) == 1) THEN
		 cz = cx ! character constants
	ELSE
		 cz = cy
	END IF

!   Are we expecting large *Bounds* forces, or not?

	ic = 0
	hasBC = .FALSE.
	vBCa = 0.0D0
	DO 910 k = 1, nCond
		IF (nodCon(k) == node) THEN
			hasBC = .TRUE.
			GO TO 911
		END IF
910       CONTINUE
911       IF (hasBC) THEN
		 size = large
	ELSE
		 size = small
	END IF

	IF (log_force_balance .AND. Verbose) WRITE (iUnitVerb, 5) &
&          node, cz, (comp(j, i), j = 1, 6), size
	IF (size == large) THEN
		 IF (MOD(i, 2) == 1) THEN
			  sumvdf = sumvdf + v(1, node) * comp(6, i)
		 ELSE
			  sumvdf = sumvdf + v(2, node) * comp(6, i)
		 END IF
	END IF

1000  CONTINUE

!   Write *Bounds* forces to file in same format as velocities:

WRITE (iUnitF, 1010) title1
WRITE (iUnitF, 1010) title2
WRITE (iUnitF, 1010) title3
1010  FORMAT (A80)
WRITE (iUnitF, 1020) (comp(6, i), i = 1, nEntry)
1020  FORMAT (1P,5D16.8)

! Calculate the sum of torques, which should be zero:

tQxx = 0.0D0
tQyy = 0.0D0
tQzz = 0.0D0
tFxx = 0.0D0
tFyy = 0.0D0
tFzz = 0.0D0
DO 2000 i = 1, numNod
  ix = 2 * i - 1
  iy = ix + 1
  theta = xNode(i)
  phai = yNode(i)
  fx = COS(theta) * COS(phai) * comp(1, ix) &
&       - SIN(phai) * comp(1, iy)
  fy = COS(theta) * SIN(phai) * comp(1, ix) &
&       + COS(phai) * comp(1, iy)
  fz = -SIN(theta) * comp(1, ix)
  tQxx = tQxx + radius * (SIN(theta) * SIN(phai) * fz - &
&                      COS(theta) * fy)
  tQyy = tQyy + radius * (COS(theta) * fx - &
&                      SIN(theta) * COS(phai) * fz)
  tQzz = tQzz + radius * (SIN(theta) * COS(phai) * fy - &
&                      SIN(theta) * SIN(phai) * fx)
  tFxx = tFxx + radius * ABS(SIN(theta) * SIN(phai) * fz - &
&                         COS(theta) * fy)
  tFyy = tFyy + radius * ABS(COS(theta) * fx - &
&                         SIN(theta) * COS(phai) * fz)
  tFzz = tFzz + radius * ABS(SIN(theta) * COS(phai) * fy - &
&                         SIN(theta) * SIN(phai) * fx)
2000  CONTINUE
IF (tFxx > 0.0D0) THEN
	rQxx = ABS(tQxx / tFxx)
ELSE
	rQxx = 0.0D0
END IF
IF (tFyy > 0.0D0) THEN
	rQyy = ABS(tQyy / tFyy)
ELSE
	rQyy = 0.0D0
END IF
IF (tFzz > 0.0D0) THEN
	rQzz = ABS(tQzz / tFzz)
ELSE
	rQzz = 0.0D0
END IF
IF(Verbose) WRITE(iUnitVerb, 2001) tQxx, tQyy, tQzz, &
&                     tFxx, tFyy, tFzz, &
&                     rQxx, rQyy, rQzz
2001  FORMAT(/' Net torque from all standardized surface' &
&            ,' traction anomalies:' &
&        /'    X=0N,0E  Y=0N,90E     Z=90N  <- axes' &
&                    ,' (through center of planet)' &
&        /' ',1P,3D10.2,'  <- sum of torque' &
&        /' ',   3D10.2,'  <- sum of ABS(torque)' &
&        /' ',   3D10.2,'  <- quotient (fractional' &
&                                    ,' error)')
IF(Verbose) WRITE (iUnitVerb, 10)
RETURN
END SUBROUTINE Balanc

SUBROUTINE BuildF (area, detJ, dXS, dYS, eta, & ! input
&                    fBase, fDip, fLen, &
&                    fPFlt, fPSfer, fArg, fTStar, &
&                    mxDOF, mxEl, mxFEl, &
&                    nDOF, nFl, nodeF, nodes, &
&                    numEl, oVB, pulled, radius, &
&                    sita, tOfset, trHMax, &
&                    wedge, &
&                    force)                       ! output

!   Compute forcing vector "force": Includes fixed terms from fBase
!   (mostly gravitational spreading), plus variable terms:
!   *From triangular continuum elements:
!      'pre-stress' or intercept-stress on linearized flow-laws,
!       and basal shear stress forces.
!   *From dipping, oblique-slip fault elements:
!      'initial traction' used in linearization of rheology.

!   In both cases, a small element vector is formed first,
!   with local node numbers, and then transferred to the global
!   forcing vector; this simplifies addressing.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: area, detJ, dXS, dYS, eta                                         ! input
DOUBLE PRECISION, INTENT(IN) :: fBase                                                   ! input
REAL*8, INTENT(IN) :: fDip, fLen, fPFlt, fPSfer, fArg, fTStar                           ! input
INTEGER, INTENT(IN) :: mxDOF, mxEl, mxFEl, nDOF, nFl, nodeF, nodes, numEl               ! input
REAL*8, INTENT(IN) :: oVB                                                               ! input
LOGICAL, INTENT(IN) :: pulled                                                           ! input
REAL*8, INTENT(IN) :: radius, sita, tOfset, trHMax, wedge                               ! input
DOUBLE PRECISION, INTENT(OUT) :: force                                                  ! output
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
DOUBLE PRECISION :: fPhi, fGauss
DOUBLE PRECISION :: weight
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
COMMON / FPhis /  fPhi
COMMON / FGList / fGauss
COMMON / WgtVec / weight
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
INTEGER i, j, jj, ju, jv, m
REAL*8 ads11, ads12, ads21, ads22, angle, cosA, dA, dip, dS, &
	& fP, oSinD, oSin2D, &
	& sHx, sHy, sinA, tOxx, tOxy, tOyy
DOUBLE PRECISION :: dUxx, dUxy, dUyy, elE, elF, sinS, tanS
DIMENSION fPhi(4, 7), fGauss(7)
DIMENSION weight(7)
DIMENSION area(mxEl), detJ(7, mxEl), &
&           dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
&           ele(6), elf(8), eta(7, mxEl), &
&           fBase(mxDOF), fDip(2, mxFEl), fLen(mxFEl), fP(2, 2, 4), &
&           force(mxDOF, 1), fPFlt(2, 2, 2, 7, mxFEl), fArg(2, mxFEl), &
&           fPSfer(2, 2, 3, 7, mxEl), fTStar(2, 7, mxFEl), &
&           nodeF(4, mxFEl), nodes(3, mxEl), oVB(2, 7, mxEl), &
&           pulled(7, mxEl), &
&           sita(7, mxEl), tOfset(3, 7, mxEl)

!   Begin with constant terms (the same in each iteration):

DO 10 i = 1, nDOF
	force(i, 1) = fBase(i) ! N.B. The second (dummy) column subscript on "force" is required by MKL conventions.
10  CONTINUE

!   Contributions of triangular continuum elements:

DO 1000 i = 1, numEl
	DO 20 j = 1, 6
		 elE(j) = 0.0D0
20       CONTINUE

!        Effects of pre-stress:

	DO 100 m = 1, 7
		 dA = area(i) * weight(m) * detJ(m, i) / radius
		 sinS = 1.00D0 / SIN(sita(m, i))
		 tanS = 1.00D0 / TAN(sita(m, i))
		 tOxx = tOfset(1, m, i)
		 tOyy = tOfset(2, m, i)
		 tOxy = tOfset(3, m, i)
		 DO 90 j = 1, 3
			  ju = 2 * j - 1
			  jv = 2 * j
			  dUxx = tOxx * dXS(1, 1, j, m, i)
			  dUxy = tOxy * (dXS(1, 2, j, m, i) + dYS(1, 1, j, m, i) * sinS &
&                          - fPSfer(1, 2, j, m, i) * tanS)
			  dUyy = tOyy * (dYS(1, 2, j, m, i) * sinS &
&                          + fPSfer(1, 1, j, m, i) * tanS)
			  elE(ju) = elE(ju) - dA * (dUxx + dUxy + dUyy)
			  dUxx = tOxx * dXS(2, 1, j, m, i)
			  dUxy = tOxy * (dXS(2, 2, j, m, i) + dYS(2, 1, j, m, i) * sinS &
&                          - fPSfer(2, 2, j, m, i) * tanS)
			  dUyy = tOyy * (dYS(2, 2, j, m, i) * sinS &
&                          + fPSfer(2, 1, j, m, i) * tanS)
			  elE(jv) = elE(jv) - dA * (dUxx + dUxy + dUyy)
90            CONTINUE
100       CONTINUE

!        Basal shear stresses (if any), in case where grid doesn't move:

	IF (trHMax > 0.0D0) THEN
		 DO 200 m = 1, 7
			  IF (pulled(m, i)) THEN
				   dA = area(i) * weight(m) * detJ(m, i)
				   sHx = oVB(1, m, i) * eta(m, i)
				   sHy = oVB(2, m, i) * eta(m, i)
				   DO 190 j = 1, 3
						ju = 2 * j - 1
						jv = 2 * j
						elE(ju) = elE(ju) &
&                                 + dA * (sHx * fPSfer(1, 1, j, m, i) &
&                                       + sHy * fPSfer(1, 2, j, m, i))
						elE(jv) = elE(jv) &
&                                 + dA * (sHx * fPSfer(2, 1, j, m, i) &
&                                       + sHy * fPSfer(2, 2, j, m, i))
190                      CONTINUE
			  END IF
200            CONTINUE
	END IF

!      Move entries of continuum-element force vector into global vector

	DO 900 j = 1, 3
		 jv = 2 * nodes(j, i)
		 ju = jv - 1
		 force(ju, 1) = force(ju, 1) + elE(2 * j - 1)
		 force(jv, 1) = force(jv, 1) + elE(2 * j)
900       CONTINUE
1000  CONTINUE

!   Contributions from dipping, oblique-slip fault elements:

DO 2000 i = 1, nFl
	DO 1020 j = 1, 8
		 elF(j) = 0.0D0
1020       CONTINUE

!        Effects of artificial 'initial traction' (fTStar):

	DO 1100 m = 1, 7
		 dip = fPhi(1, m) * fDip(1, i) + fPhi(2, m) * fDip(2, i)
		 IF (ABS(dip - 1.57079632679490D0) > wedge) THEN
			  oSinD = 1.0D0 / SIN(dip)
			  oSin2D = 1.0D0 / SIN(2.0D0 * dip)

!CCCC                 angle = fArg(1, i) * fPhi(1, m) + fArg(2, i) * fPhi(2, m)
!CCCC                 Line above was replaced due to cycle-shift problem

			  angle = Chord(fArg(1, i), fPhi(2, m), fArg(2, i))

			  sinA = SIN(angle)
			  cosA = COS(angle)
			  dS = fLen(i) * fGauss(m)
			  DO 1030 j = 1, 2
				 jj = 4
				 IF (j == 2) jj = 3
				 fP(1, 1, j) = fPFlt(1, 1, j, m, i)
				 fP(2, 1, j) = fPFlt(2, 1, j, m, i)
				 fP(1, 2, j) = fPFlt(1, 2, j, m, i)
				 fP(2, 2, j) = fPFlt(2, 2, j, m, i)
				 fP(1, 1, jj) = -fP(1, 1, j)
				 fP(2, 1, jj) = -fP(2, 1, j)
				 fP(1, 2, jj) = -fP(1, 2, j)
				 fP(2, 2, jj) = -fP(2, 2, j)
1030                 CONTINUE
			  DO 1090 j = 1, 4
				   ju = 2 * j - 1
				   jv = 2 * j
				   ads11 = fP(1, 1, j) * cosA + fP(1, 2, j) * sinA
				   ads12 = fP(1, 1, j) * sinA - fP(1, 2, j) * cosA
				   ads21 = fP(2, 1, j) * cosA + fP(2, 2, j) * sinA
				   ads22 = fP(2, 1, j) * sinA - fP(2, 2, j) * cosA
				   elF(ju) = elF(ju) - dS * (ads11 * fTStar(1, m, i) * oSinD &
&                            - 2.0D0 * ads12 * oSin2D * fTStar(2, m, i))
				   elF(jv) = elF(jv) - dS * (ads21 * fTStar(1, m, i) * oSinD &
&                            - 2.0D0 * ads22 * oSin2D * fTStar(2, m, i))
1090                 CONTINUE
		 END IF
1100       CONTINUE

!       Move entries of fault-element force vector into global vector:

	DO 1900 j = 1, 4
			  jv = 2 * nodeF(j, i)
			  ju = jv - 1
			  force(ju, 1) = force(ju, 1) + elF(2 * j - 1)
			  force(jv, 1) = force(jv, 1) + elF(2 * j)
1900       CONTINUE
2000  CONTINUE
RETURN
END SUBROUTINE BuildF

SUBROUTINE BuildK (alpha, area, detJ, dXS, dYS, & ! input
&                    eta, fPSfer, &
&                    mxEl, &
&                    nodes, numEl, pulled, &
&                    radius, sita, trHMax, &
&                    stiff)                         ! output

!   Computes stiffness matrix "stiff" (alias "K" or "k" in some other subprograms)
!   which represents stiffness of triangular continuum elements,
!   from tensor "alpha" and derivitives of nodal functions
!   of the element grid.

!   Also adds diagonal stiffening associated with shear coupling to
!   the mantle beneath, if any.

!   Note that the stiffness associated with fault elements is not
!   included here (for historical reasons).  See subprogram -AddFSt-.

!    A two-step process is used:
!    -A stiffness matrix for each element is formed, using
!     generic node numbering, 1-3.  Each entry in this matrix is
!     a 2 x 2 submatrix, because node velocities have two components.
!    -The element stiffness matrix terms are added to the global
!     stiffness matrix.  (This step involves complex indirect
!     addressing, and is very difficult to optimize).

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alpha, area, detJ, dXS, dYS, eta, FPSfer                          ! input
INTEGER, INTENT(IN) :: mxEl, nodes, numEl                                               ! input
LOGICAL, INTENT(IN) :: pulled                                                           ! input
REAL*8, INTENT(IN) :: radius, sita, trHMax                                              ! input
DOUBLE PRECISION, INTENT(OUT) :: stiff                                                  ! output
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
DOUBLE PRECISION weight
COMMON / WgtVec / weight
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
INTEGER i, i3, iq, irx, iry, &
	 & j3, jcx, jcy, m, nodeI, nodeJ
REAL*8 etadA
DOUBLE PRECISION dA, elK, fl11j, fl12i, fl12j, fl22i, fl22j, sinS, sum, tanS
DIMENSION alpha(3, 3, 7, mxEl), area(mxEl), detJ(7, mxEl), &
&           dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
&           eta(7, mxEl), fPSfer(2, 2, 3, 7, mxEl), &
&           nodes(3, mxEl), pulled(7, mxEl), &
&           sita(7, mxEl)
DIMENSION elK(2, 2, 3, 3), weight(7), stiff(nKRows, nRank)

!   Begin by zeroing the matrix; all other logic will add to it.

stiff = 0.0D0 ! whole matrix (per power of Fortran 90)

!   Major loop is on triangular continuum elements:

DO 500 i = 1, numEl

!         Zero and then build up the element stiffness matrix:

	DO 10 i3 = 1, 3
		 DO 9 j3 = 1, 3
			  elK(1, 1, i3, j3) = 0.0D0
			  elK(1, 2, i3, j3) = 0.0D0
			  elK(2, 1, i3, j3) = 0.0D0
			  elK(2, 2, i3, j3) = 0.0D0
9            CONTINUE
10       CONTINUE

!           Incorporate stiffness tensors "alpha":

	DO 90 i3 = 1, 3
		 DO 80 j3 = 1, 3

!                     upper left terms: X-coefficients in X-balance:
			  sum = 0.D0
			  DO 40 m = 1, 7
				 dA = weight(m) * detJ(m, i)
				 sinS = 1.0D0 / SIN(sita(m, i))
				 tanS = 1.0D0 / TAN(sita(m, i))
				 fl11j = alpha(1, 1, m, i) * dXS(1, 1, j3, m, i) &
&                         + 0.5D0 * alpha(1, 3, m, i) * (dYS(1, 1, j3, m, i) * sinS &
&                         + dXS(1, 2, j3, m, i) - fPSfer(1, 2, j3, m, i) * tanS) &
&                         + alpha(1, 2, m, i) * (dYS(1, 2, j3, m, i) * sinS &
&                         + fPSfer(1, 1, j3, m, i) * tanS)
				 fl12j = alpha(3, 1, m, i) * dXS(1, 1, j3, m, i) &
&                         + 0.5D0 * alpha(3, 3, m, i) * (dYS(1, 1, j3, m, i) * sinS &
&                         + dXS(1, 2, j3, m, i) - fPSfer(1, 2, j3, m, i) * tanS) &
&                         + alpha(3, 2, m, i) * (dYS(1, 2, j3, m, i) * sinS &
&                         + fPSfer(1, 1, j3, m, i) * tanS)
				 fl12i = dXS(1, 2, i3, m, i) + dYS(1, 1, i3, m, i) * sinS &
&                         - fPSfer(1, 2, i3, m, i) * tanS
				 fl22j = alpha(2, 1, m, i) * dXS(1, 1, j3, m, i) &
&                         + 0.5D0 * alpha(2, 3, m, i) * (dYS(1, 1, j3, m, i) * sinS &
&                         + dXS(1, 2, j3, m, i) - fPSfer(1, 2, j3, m, i) * tanS) &
&                         + alpha(2, 2, m, i) * (dYS(1, 2, j3, m, i) * sinS &
&                         + fPSfer(1, 1, j3, m, i) * tanS)
				 fl22i = dYS(1, 2, i3, m, i) * sinS &
&                         + fPSfer(1, 1, i3, m, i) * tanS
				 sum = sum + dA * &
&                             (fl11j * dXS(1, 1, i3, m, i) + fl12j * fl12i &
&                             + fl22j * fl22i)
40                 CONTINUE
			  elK(1, 1, i3, j3) = elK(1, 1, i3, j3) + &
&                               sum * area(i) / (radius * radius)

!                     lower right terms: Y-coefficients in Y-balance:
			  sum = 0.D0
			  DO 50 m = 1, 7
				 dA = weight(m) * detJ(m, i)
				 sinS = 1.0D0 / SIN(sita(m, i))
				 tanS = 1.0D0 / TAN(sita(m, i))
				 fl11j = alpha(1, 1, m, i) * dXS(2, 1, j3, m, i) &
&                         + 0.5D0 * alpha(1, 3, m, i) * (dYS(2, 1, j3, m, i) * sinS &
&                         + dXS(2, 2, j3, m, i) - fPSfer(2, 2, j3, m, i) * tanS) &
&                         + alpha(1, 2, m, i) * (dYS(2, 2, j3, m, i) * sinS &
&                         + fPSfer(2, 1, j3, m, i) * tanS)
				 fl12j = alpha(3, 1, m, i) * dXS(2, 1, j3, m, i) &
&                         + 0.5D0 * alpha(3, 3, m, i) * (dYS(2, 1, j3, m, i) * sinS &
&                         + dXS(2, 2, j3, m, i) - fPSfer(2, 2, j3, m, i) * tanS) &
&                         + alpha(3, 2, m, i) * (dYS(2, 2, j3, m, i) * sinS &
&                         + fPSfer(2, 1, j3, m, i) * tanS)
				 fl12i = dXS(2, 2, i3, m, i) + dYS(2, 1, i3, m, i) * sinS &
&                         - fPSfer(2, 2, i3, m, i) * tanS
				 fl22j = alpha(2, 1, m, i) * dXS(2, 1, j3, m, i) &
&                         + 0.5D0 * alpha(2, 3, m, i) * (dYS(2, 1, j3, m, i) * sinS &
&                         + dXS(2, 2, j3, m, i) - fPSfer(2, 2, j3, m, i) * tanS) &
&                         + alpha(2, 2, m, i) * (dYS(2, 2, j3, m, i) * sinS &
&                         + fPSfer(2, 1, j3, m, i) * tanS)
				 fl22i = dYS(2, 2, i3, m, i) * sinS &
&                         + fPSfer(2, 1, i3, m, i) * tanS
				 sum = sum + dA * (fl11j * dXS(2, 1, i3, m, i) &
&                              + fl12j * fl12i + fl22j * fl22i)
50                 CONTINUE
			  elK(2, 2, i3, j3) = elK(2, 2, i3, j3) + &
&                               sum * area(i) / (radius * radius)

!                     upper right terms: Y-coefficients in X-balance:
			  sum = 0.0D0
			  DO 60 m = 1, 7
				 dA = weight(m) * detJ(m, i)
				 tanS = 1.0D0 / TAN(sita(m, i))
				 sinS = 1.0D0 / SIN(sita(m, i))
				 fl11j = alpha(1, 1, m, i) * dXS(2, 1, j3, m, i) &
&                         + 0.5D0 * alpha(1, 3, m, i) * (dYS(2, 1, j3, m, i) * sinS &
&                         + dXS(2, 2, j3, m, i) - fPSfer(2, 2, j3, m, i) * tanS) &
&                         + alpha(1, 2, m, i) * (dYS(2, 2, j3, m, i) * sinS &
&                         + fPSfer(2, 1, j3, m, i) * tanS)
				 fl12j = alpha(3, 1, m, i) * dXS(2, 1, j3, m, i) &
&                         + 0.5D0 * alpha(3, 3, m, i) * (dYS(2, 1, j3, m, i) * sinS &
&                         + dXS(2, 2, j3, m, i) - fPSfer(2, 2, j3, m, i) * tanS) &
&                         + alpha(3, 2, m, i) * (dYS(2, 2, j3, m, i) * sinS &
&                         + fPSfer(2, 1, j3, m, i) * tanS)
				 fl12i = dXS(1, 2, i3, m, i) + dYS(1, 1, i3, m, i) * sinS &
&                         - fPSfer(1, 2, i3, m, i) * tanS
				 fl22j = alpha(2, 1, m, i) * dXS(2, 1, j3, m, i) &
&                         + 0.5D0 * alpha(2, 3, m, i) * (dYS(2, 1, j3, m, i) * sinS &
&                         + dXS(2, 2, j3, m, i) - fPSfer(2, 2, j3, m, i) * tanS) &
&                         + alpha(2, 2, m, i) * (dYS(2, 2, j3, m, i) * sinS &
&                         + fPSfer(2, 1, j3, m, i) * tanS)
				 fl22i = dYS(1, 2, i3, m, i) * sinS &
&                         + fPSfer(1, 1, i3, m, i) * tanS
				 sum = sum + dA * (fl11j * dXS(1, 1, i3, m, i) &
&                              + fl12j * fl12i + fl22j * fl22i)
60                 CONTINUE
			  elK(1, 2, i3, j3) = elK(1, 2, i3, j3) + &
&                               sum * area(i) / (radius * radius)

!                     lower left terms: X-coefficients in Y-balance:
			  sum = 0.0D0
			  DO 70 m = 1, 7
				 dA = weight(m) * detJ(m, i)
				 sinS = 1.0D0 / SIN(sita(m, i))
				 tanS = 1.0D0 / TAN(sita(m, i))
				 fl11j = alpha(1, 1, m, i) * dXS(1, 1, j3, m, i) &
&                         + 0.5D0 * alpha(1, 3, m, i) * (dYS(1, 1, j3, m, i) * sinS &
&                         + dXS(1, 2, j3, m, i) - fPSfer(1, 2, j3, m, i) * tanS) &
&                         + alpha(1, 2, m, i) * (dYS(1, 2, j3, m, i) * sinS &
&                         + fPSfer(1, 1, j3, m, i) * tanS)
				 fl12j = alpha(3, 1, m, i) * dXS(1, 1, j3, m, i) &
&                         + 0.5D0 * alpha(3, 3, m, i) * (dYS(1, 1, j3, m, i) * sinS &
&                         + dXS(1, 2, j3, m, i) - fPSfer(1, 2, j3, m, i) * tanS) &
&                         + alpha(3, 2, m, i) * (dYS(1, 2, j3, m, i) * sinS &
&                         + fPSfer(1, 1, j3, m, i) * tanS)
				 fl12i = dXS(2, 2, i3, m, i) + dYS(2, 1, i3, m, i) * sinS &
&                         - fPSfer(2, 2, i3, m, i) * tanS
				 fl22j = alpha(2, 1, m, i) * dXS(1, 1, j3, m, i) &
&                         + 0.5D0 * alpha(2, 3, m, i) * (dYS(1, 1, j3, m, i) * sinS &
&                         + dXS(1, 2, j3, m, i) - fPSfer(1, 2, j3, m, i) * tanS) &
&                         + alpha(2, 2, m, i) * (dYS(1, 2, j3, m, i) * sinS &
&                         + fPSfer(1, 1, j3, m, i) * tanS)
				 fl22i = dYS(2, 2, i3, m, i) * sinS &
&                         + fPSfer(2, 1, i3, m, i) * tanS
				 sum = sum + dA * &
&                             (fl11j * dXS(2, 1, i3, m, i) + fl12j * fl12i &
&                             + fl22j * fl22i)
70                 CONTINUE
			  elK(2, 1, i3, j3) = elK(2, 1, i3, j3) + &
&                               sum * area(i) / (radius * radius)

80            CONTINUE
90       CONTINUE

!          Add any diagonal stiffness associated with viscous basal drag

	IF (trHMax > 0.D0) THEN
		 DO 200 m = 1, 7
			  IF (pulled(m, i)) THEN
				   etadA = eta(m, i) * weight(m) * area(i) * detJ(m, i)
				   DO 190 i3 = 1, 3
						DO 180 j3 = 1, 3
							 elK(1, 1, i3, j3) = elK(1, 1, i3, j3) + &
&                                      etadA * (fPSfer(1, 1, i3, m, i) * &
&                                               fPSfer(1, 1, j3, m, i) + &
&                                               fPSfer(1, 2, i3, m, i) * &
&                                               fPSfer(1, 2, j3, m, i))
							 elK(1, 2, i3, j3) = elK(1, 2, i3, j3) + &
&                                      etadA * (fPSfer(1, 1, i3, m, i) * &
&                                               fPSfer(2, 1, j3, m, i) + &
&                                               fPSfer(1, 2, i3, m, i) * &
&                                               fPSfer(2, 2, j3, m, i))
							 elK(2, 1, i3, j3) = elK(2, 1, i3, j3) + &
&                                      etadA * (fPSfer(2, 1, i3, m, i) * &
&                                               fPSfer(1, 1, j3, m, i) + &
&                                               fPSfer(2, 2, i3, m, i) * &
&                                               fPSfer(1, 2, j3, m, i))
							 elK(2, 2, i3, j3) = elK(2, 2, i3, j3) + &
&                                      etadA * (fPSfer(2, 1, i3, m, i) * &
&                                               fPSfer(2, 1, j3, m, i) + &
&                                               fPSfer(2, 2, i3, m, i) * &
&                                               fPSfer(2, 2, j3, m, i))
180                             CONTINUE
190                      CONTINUE
			  END IF
200            CONTINUE
	END IF

!           Apply element matrix to augment global stiffness matrix:

	DO 400 i3 = 1, 3
	   nodeI = nodes(i3, i)
	   iry = 2 * nodeI
	   irx = iry - 1
	   DO 300 j3 = 1, 3
		   nodeJ = nodes(j3, i)
		   jcy = 2 * nodeJ
		   jcx = jcy - 1
		  !matrix element(irx, jcx):
		   iq = iDiagonal + irx - jcx
		   stiff(iq, jcx) = stiff(iq, jcx) + elK(1, 1, i3, j3)
		  !matrix element(irx, jcy):
		   iq = iDiagonal + irx - jcy
		   stiff(iq, jcy) = stiff(iq, jcy) + elK(1, 2, i3, j3)
		  !matrix element(iry, jcx):
		   iq = iDiagonal + iry - jcx
		   stiff(iq, jcx) = stiff(iq, jcx) + elK(2, 1, i3, j3)
		  !matrix element(iry, jcy):
		   iq = iDiagonal + iry - jcy
		   stiff(iq, jcy) = stiff(iq, jcy) + elK(2, 2, i3, j3)
300          CONTINUE
400       CONTINUE
500  CONTINUE
RETURN
END SUBROUTINE BuildK

REAL*8 FUNCTION Chord (angle1, s, angle2)

!   Returns an angle obtained by interpolation between "angle1"
!   and "angle2".  The interpolation method is NOT sensitive to any
!   possible cycle shifts (of 2*n*Pi) between angle1 and angle2.

!   Unit vectors are constructed for angle1 and angle2, and a
!   linear chord is drawn between their tips.

!   DOUBLE PRECISION s is the internal coordinate along the chord;
!   it is dimensionless, with value 0.0D0 at angle1 and 1.0D0 at
!   angle2.  (The user may input "s" values outside this range
!   to get results outside the (smaller) angle between angle1 and
!   angle2, if desired.)  The angle returned is that from the
!   origin to this chord point.

!   This algorithm should work equally well for angles measured
!   either clockwise or counterclockwise from any reference, as
!   long as the usage is consistent.

!   Both the inputs angle1, angle2 and the result Chord are in units of radians.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: angle1, angle2                                                    ! input
DOUBLE PRECISION, INTENT(IN) :: s                                                       ! input
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
REAL*8 uvec1, uvec2, vecs
DIMENSION uvec1(2), uvec2(2), vecs(2)
uvec1(1) = COS(angle1)
uvec1(2) = SIN(angle1)
uvec2(1) = COS(angle2)
uvec2(2) = SIN(angle2)
vecs(1) = (1.0D0 - s) * uvec1(1) + s * uvec2(1)
vecs(2) = (1.0D0 - s) * uvec1(2) + s * uvec2(2)
Chord = ATan2F(vecs(2), vecs(1))
RETURN
END FUNCTION Chord

SUBROUTINE Convec (iConve, iPAfri, iPVRef, iUnitM, iUnitT, & ! input
&                    mxNode, &
&                    names, &
&                    nPlate, numNod, &
&                    omega, radius, vTimes, &
&                    whichP, xNode, yNode, &
&                    vM)                                       ! output

!   Computes lower-mantle flow velocity below asthenosphere;
!   or, if iConve == 5, computes velocity of subducting plate(s).
!   Note that no code is provided here for case of iConve == 6;
!   that is handled elsewhere.

!   Computation strategy varies by model. For many, data files
!      must be read from unit iUnitM.

!   For all models except #5, the factor vTimes is applied.

!   Velocities are initially computed in the Africa-fixed
!   reference frame (for historical reasons); then they are
!   transformed to appear in the reference frame of plate
!   #iPVRef; this is done by a common transformation at the end of
!   this routine.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iConve, iPAfri, iPVRef, iUnitM, iUnitT, mxNode                   ! input
CHARACTER*2, INTENT(IN) :: names                                                        ! input
INTEGER, INTENT(IN) :: nPlate, numNod                                                   ! input
REAL*8, INTENT(IN) :: omega, radius, vTimes                                             ! input
INTEGER, INTENT(IN) :: whichP                                                           ! input
REAL*8, INTENT(IN) :: xNode, yNode                                                      ! input
DOUBLE PRECISION, INTENT(OUT) :: vM                                                     ! output
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
REAL*8 CosDeg, SinDeg, deg
!      - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
CHARACTER*2 c2
INTEGER i, ios, iPlate, ir, irBot, irTop, iSouth, &
	 & j, jc, jcLeft, jcRigh, jEast, jVec, k, nCross, numVec
REAL*8 baum88, eLon, eLon1, eLon2, fE, frac, fS, &
	& hoc792, hx, hy, hz, lat1, lat2, lon1, lon2, &
	& nLat, nLat1, nLat2, omegax, omegay, omegaz, phi, phix, phiy, phiz, &
	& r2, r2min, test, theta, thetax, thetay, thetaz, tx, ty, tz, &
	& vBot, vPhi, vTheta, vTop, vx, vy, vz, xn, yn, zn
CHARACTER*27 endSeg
DIMENSION xNode(mxNode), yNode(mxNode), vm(2, mxNode)
DIMENSION hoc792(2, -8:8, 1:36)
DIMENSION baum88(5, 1000)
DIMENSION names(nPlate), omega(3, nPlate)
DIMENSION whichP(mxNode)
!     - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -

!  Data specific to the iConve == 5 case.
!    (Fortran 90 code added 2000.04.14):

!  ID code of subducting plate (or, in case of several subducting
!     plates, ID code of the largest plate (applied in any areas NOT
!     outlined by .dig outlines):
CHARACTER(2) :: underPlate, otherPlate
INTEGER :: iUnderPlate
!  Number (may be 0) of additional plates, each to be represented
!     with a digitised outline:
INTEGER ::n_others
!  ID codes of other subducting plates:
INTEGER, DIMENSION(:), ALLOCATABLE :: iOtherPlate
!  Counts of points in each of the "other" digitized outlines:
INTEGER, DIMENSION(:), ALLOCATABLE :: other_counts
!  Largest value found in other_counts:
INTEGER :: max_count
!  Storage for digitized outlines, in (lon, lat) format, in
!     decimal degrees (+ = N, E; - = S, W).
REAL*8, DIMENSION(:, :, :), ALLOCATABLE :: other_shapes

!     - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -

!   Statement functions:
CosDeg(deg) = COS(deg * 0.0174532925199433D0)
SinDeg(deg) = SIN(deg * 0.0174532925199433D0)

IF (iConve == 0) THEN

	DO 99 i = 1, numNod
		 vM(1, i) = 0.0D0
		 vM(2, i) = 0.0D0
!                Note: This is in Africa-fixed reference frame;
!                see below at end of routine for transformation.
99       CONTINUE

ELSE IF (iConve == 1) THEN

!           Hager and O'Connell (1979) viscosity model II

!           Read from file "HOC79II.DIG"
!           Vectors are every 10 degrees in latitude and longitude
!           Columns march East from 10E to 360E.
!           Within each column, travel is S from 80N to 80S.
!           Units of input data are degrees East and North.
!           2nd end of line segment shows where the grid point
!           will be displaced to after 50 m.y. of flow.

	IF(Verbose) WRITE(iUnitVerb, 100) iUnitM
100       FORMAT (' Attempting to read plate OUTLINES from unit ',I3/)
	DO 140 jEast = 1, 36
		 DO 130 iSouth = -8, 8
			  READ (iUnitM, * , END = 101, ERR = 101) eLon1, nLat1
			  GO TO 103
!                     -------------------- ERR0R HANDLER ----------
101  write(ErrorMsg,'(A/,A,I3/,A,I2,A,I2,A/,A)') "ERR0R IN -Convec-:", &
&                        "WHILE READING MANTLE VELOCITIES FROM UNIT ",iUnitM, &
&                        "TO FILL IN COLUMN ",jEast,", ROW ",iSouth ,"ENCOUNTERED A RECORD WHICH DOES NOT", &
&                        "HOLD TWO RECOGNIZABLE NUMBERS."
     call FatalError(ErrorMsg,ThID)
!                     ---------------------------------------------
103                 jc = (eLon1 / 10.0D0) + 0.50D0
			  IF (nLat1 >= 0.) THEN
				   ir = (nLat1 / 10.0D0) + 0.50D0
				   ir = -ir
			  ELSE
				   ir = (-nLat1 / 10.0D0) + 0.50D0
			  END IF
			  IF ((jc /= jEast).OR.(ir /= iSouth)) THEN
                 write(ErrorMsg,'(A/,A,I3/,A,I2,A,I2/,A,I2,A,I2/,A,F7.2,A,F6.2)') "ERR0R: WHILE READING LOWER-MANTLE", &
&                        "FLOW VECTORS FROM UNIT ",iUnitM, &
&                        "AND LOOKING FOR ROW ",iSouth,", COLUMN ",jEast, &
&                        "ENCOUNTERED     ROW ",ir,", COLUMN ",jc, &
&                        "(LONGITUDE ",eLon1,", LATITUDE ",nLat1
                 call FatalError(ErrorMsg,ThID)
			  END IF
			  READ (iUnitM, * , ERR = 101, END = 101) eLon2, nLat2
			  READ (iUnitM, '(A)') endSeg
			  tx = CosDeg(nLat1) * CosDeg(eLon1)
			  ty = CosDeg(nLat1) * SinDeg(eLon1)
			  tz = SinDeg(nLat1)
			  hx = CosDeg(nLat2) * CosDeg(eLon2)
			  hy = CosDeg(nLat2) * SinDeg(eLon2)
			  hz = SinDeg(nLat2)
			  vx = (hx - tx) * radius / (50.D6 * 3.15576D7)
			  vy = (hy - ty) * radius / (50.D6 * 3.15576D7)
			  vz = (hz - tz) * radius / (50.D6 * 3.15576D7)
			  thetax = SinDeg(nLat1) * CosDeg(eLon1)
			  thetay = SinDeg(nLat1) * SinDeg(eLon1)
			  thetaz = -CosDeg(nLat1)
			  vTheta = vx * thetax + vy * thetay + vz * thetaz
			  phix = -SinDeg(eLon1)
			  phiy = CosDeg(eLon1)
			  phiz = 0.0D0
			  vPhi = vx * phix + vy * phiy + vz * phiz
			  hoc792(1, ir, jc) = vTheta
			  hoc792(2, ir, jc) = vPhi
130            CONTINUE
140       CONTINUE
	DO 190 i = 1, numNod
		 nLat1 = 90.0D0 - xNode(i) * 57.2957795130823D0
		 nLat1 = MIN(nLat1, + 80.0D0)
		 nLat1 = MAX(nLat1, -80.0D0)
		 eLon1 = yNode(i) * 57.2957795130823D0
		 IF (eLon1 < 0.) eLon1 = eLon1 + 360.0D0
		 IF (eLon1 < 0.) eLon1 = eLon1 + 360.0D0
		 IF (eLon1 > 360.) eLon1 = eLon1 - 360.0D0
		 IF (nLat1 >= 0.) THEN
			  irTop = (nLat1 / 10.0D0) + 1.0D0
			  irTop = -irTop
		 ELSE
			  irTop = (-nLat1 / 10.0D0)
		 END IF
		 IF (irTop < 8) THEN
			  irBot = irTop + 1
			  fS = (-irTop * 10.0D0 - nLat1) / 10.0D0
		 ELSE
			  irBot = irTop
			  fS = 0.
		 END IF
		 jcRigh = eLon1 / 10.0D0 + 1.0D0
		 jcRigh = MIN(jcRigh, 36)
		 IF (jcRigh > 1) THEN
			  jcLeft = jcRigh - 1
			  fE = (eLon1 - 10.0D0 * jcLeft) / 10.0D0
		 ELSE
			  jcLeft = 36
			  fE = eLon1 / 10.0D0
		 END IF
		 vTop = hoc792(1, irTop, jcLeft) + &
&                 (hoc792(1, irTop, jcRigh) - hoc792(1, irTop, jcLeft)) * fe
		 vBot = hoc792(1, irBot, jcLeft) + &
&                 (hoc792(1, irBot, jcRigh) - hoc792(1, irBot, jcLeft)) * fe
		 vM(1, i) = vTop + (vBot - vTop) * fS
		 vTop = hoc792(2, irTop, jcLeft) + &
&                 (hoc792(2, irTop, jcRigh) - hoc792(2, irTop, jcLeft)) * fe
		 vBot = hoc792(2, irBot, jcLeft) + &
&                 (hoc792(2, irBot, jcRigh) - hoc792(2, irBot, jcLeft)) * fe
		 vM(2, i) = vTop + (vBot - vTop) * fS
		 vM(1, i) = vM(1, i) * vTimes
		 vM(2, i) = vM(2, i) * vTimes
190       CONTINUE

ELSE IF (iConve == 2) THEN

!           Baumgardner (1988) Figure 7, parts A-F

!           Read from file "BAUM887.DIG"
!           Vectors are in random order, about 729 in all.
!           Units of input data are degrees East and North.
!           2nd end of line segment shows where the grid point
!           will be displaced to after 11 m.y. of flow.
!          (Time would be 110 m.y., but he says to scale V up
!           *10 because Earth's Rayleigh number is higher that
!           that of the model.)

	IF(Verbose) WRITE(iUnitVerb, 200) iUnitM
200       FORMAT (' Attempting to read Baumgardner [1988] mantle', &
&              ' flow from unit ',I3/)
	numVec = 0
	DO 220 jVec = 1, 1000
		 READ (iUnitM, * , END = 221, ERR = 201) eLon1, nLat1
		 GO TO 203
!                -------------------- ERR0R HANDLER ----------
201  write(ErrorMsg,'(A/,A,I3,A,I2/,A)') "ERR0R IN -Convec-:", &
&                   " WHILE READING MANTLE VELOCITIES FROM UNIT ",iUnitM," TO FILL IN VECTOR ",jVec, &
&                   " ENCOUNTERED A RECORD WHICH DOES NOT HOLD TWO RECOGNIZABLE NUMBERS."
     call FatalError(ErrorMsg,ThID)
!                ---------------------------------------------
203      READ (iUnitM, * , ERR = 201, END = 221) eLon2, nLat2
		 READ (iUnitM, '(A)') endSeg
		 tx = CosDeg(nLat1) * CosDeg(eLon1)
		 ty = CosDeg(nLat1) * SinDeg(eLon1)
		 tz = SinDeg(nLat1)
		 hx = CosDeg(nLat2) * CosDeg(eLon2)
		 hy = CosDeg(nLat2) * SinDeg(eLon2)
		 hz = SinDeg(nLat2)
		 vx = (hx - tx) * radius / (11.D6 * 3.15576D7)
		 vy = (hy - ty) * radius / (11.D6 * 3.15576D7)
		 vz = (hz - tz) * radius / (11.D6 * 3.15576D7)
		 thetax = SinDeg(nLat1) * CosDeg(eLon1)
		 thetay = SinDeg(nLat1) * SinDeg(eLon1)
		 thetaz = -CosDeg(nLat1)
		 vTheta = vx * thetax + vy * thetay + vz * thetaz
		 phix = -SinDeg(eLon1)
		 phiy = CosDeg(eLon1)
		 phiz = 0.0D0
		 vPhi = vx * phix + vy * phiy + vz * phiz
		 baum88(1, jVec) = vTheta
		 baum88(2, jVec) = vPhi
		 baum88(3, jVec) = tx
		 baum88(4, jVec) = ty
		 baum88(5, jVec) = tz
		 numVec = numVec + 1
220       CONTINUE
221       DO 290 i = 1, numNod
		 tx = SIN(xNode(i)) * COS(yNode(i))
		 ty = SIN(xNode(i)) * SIN(yNode(i))
		 tz = COS(xNode(i))
		 r2min = 999.0D0
		 DO 280 j = 1, numVec
			  r2 = (tx - baum88(3, j))**2 + &
&                     (ty - baum88(4, j))**2 + &
&                     (tz - baum88(5, j))**2
			  IF (r2 < r2min) THEN
				   r2min = r2
				   vM(1, i) = baum88(1, j)
				   vM(2, i) = baum88(2, j)
				   vM(1, i) = vM(1, i) * vTimes
				   vM(2, i) = vM(2, i) * vTimes
			  END IF
280            CONTINUE
290       CONTINUE

ELSE IF ((iConve == 3).OR.(iConve == 4)) THEN

!          PB2002 model of Bird [2003; G**3];
!          Already has plate "names" and "omega" vectors in
!          main program (DATA statements);
!          also, plate-ID's for each node have already
!          been computed (by CALL Assign) and stored in whichP.

   DO 390 i = 1, numNod
		iPlate = whichP(i)
!               Convert to AFrica-fixed, and radians/second:
		omegax = (omega(1, iPlate) - omega(1, iPAfri)) * 3.168809D-14
		omegay = (omega(2, iPlate) - omega(2, iPAfri)) * 3.168809D-14
		omegaz = (omega(3, iPlate) - omega(3, iPAfri)) * 3.168809D-14
!               Convert to length/second:
		omegax = omegax * radius
		omegay = omegay * radius
		omegaz = omegaz * radius
!               Velocity = OMEGA x position:
		theta = xNode(i)
		phi = yNode(i)
		xn = SIN(theta) * COS(phi)
		yn = SIN(theta) * SIN(phi)
		zn = COS(theta)
		vx = omegay * zn - omegaz * yn
		vy = omegaz * xn - omegax * zn
		vz = omegax * yn - omegay * xn
!               Create unit +Theta and +Phi vectors in Cartesian:
		thetax = COS(theta) * COS(phi)
		thetay = COS(theta) * SIN(phi)
		thetaz = -SIN(theta)
		phix = -SIN(phi)
		phiy = COS(phi)
		phiz = 0.0D0
!               Find argument from dot products:
		vTheta = vx * thetax + vy * thetay + vz * thetaz
		vPhi = vx * phix + vy * phiy + vz * phiz
		vM(1, i) = vTheta * vTimes
		vM(2, i) = vPhi * vTimes
390      CONTINUE

ELSE IF (iConve == 5) THEN

!        Code added for Japan models, 2000.04; written generally so as
!        work in any case of subduction under one margin of the model,
!        where the subduction shear zone is one model boundary, which
!        is NOT represented with fault elements.  This code determines
!        which plate is subducting underneath each node, and returns
!        its velocity.  The decision about whether the subducting
!        plate is touching the model (LOGICAL pulled(m, i)) is made elsewhere.

	IF(Verbose) WRITE(iUnitVerb, 501) iUnitM
501       FORMAT ( &
&/' Attempting to read subducting plate identification code' &
&/' (and plate outlines if there is more than one subducting' &
&/' plate) from unit ',I3/)

!           Read file once for first line, and just count lengths
!           and number of digitised outlines, without saving them:

	READ(iUnitM, "(A2)", IOSTAT = ios) underPlate
	IF (ios /= 0) THEN
       write(ErrorMsg,'(A/,A/,A/,A/,A/,A/,A)') "File not found, or file empty.", &
     & "This file MUST be supplied when iConve = 5.", &
     & "First line must begin with plate code (e.g., PA).", &
     & "Second line should be ""*** END OF SEGMENT ***"".", &
     & "Then, if more than one plate is subducting, follow this", &
     & "with other plate-ID codes (e.g., PH), each followed by", &
     & "a closed outline in (lon, lat) format with decimal degrees."
       call FatalError(ErrorMsg,ThID)
	END IF
	iUnderPlate = 0
	DO 510 i = 1, nPlate
		 IF (underPlate == names(i)) THEN
			  iUnderPlate = i
			  GO TO 511
		 END IF
510       CONTINUE
511       IF (iUnderPlate == 0) THEN
            write(ErrorMsg,'(A,A2,A,I3)') "ERR0R: Illegal plate code ",underPlate," was read from unit ",iUnitM
            call FatalError(ErrorMsg,ThID)
	END IF
	n_others = 0
	max_count = 0
	ALLOCATE ( other_counts(nPlate - 1) )
	ALLOCATE ( iOtherPlate(nPlate - 1) )
!          (clear the next line, which should be *** END...)
	READ (iUnitM, * , END = 549)
	DO 540 i = 1, nPlate - 1
!                Scan file for outlines, and count lengths:
		 READ (iUnitM, "(A)", END = 549) otherPlate
		 n_others = n_others + 1
		 iOtherPlate(n_others) = 0
		 DO 520 j = 1, nPlate
			  IF (otherPlate == names(j)) THEN
				   iOtherPlate(n_others) = j
				   GO TO 521
			  END IF
520            CONTINUE
521            IF (iOtherPlate(n_others) == 0) THEN
                  write(ErrorMsg,'(A,A2,A,I3)') "ERR0R: Illegal plate code ",underPlate," was read from unit ",iUnitM
                  call FatalError(ErrorMsg,ThID)
		 END IF
		 other_counts(n_others) = 0
530            READ (iUnitM, "(A)", END = 549) c2
		 IF ((c2 == " +").OR.(c2 == " -")) THEN
			  other_counts(n_others) = other_counts(n_others) + 1
			  max_count = MAX(max_count, other_counts(n_others))
			  GO TO 530
		 END IF
540       CONTINUE
549       REWIND (iUnitM)

!           Read file again, and store digitised outlines:

	IF (n_others > 0) THEN
		 ALLOCATE ( other_shapes(n_others, 2, max_count) )
		 READ (iUnitM, *)
		 READ (iUnitM, *)
		 DO 560 i = 1, n_others
			  READ (iUnitM, *)
			  DO 555 j = 1, other_counts(i)
				   READ (iUnitM, *) other_shapes(i, 1, j), &
&                                      other_shapes(i, 2, j)
555                 CONTINUE
			  READ (iUnitM, *)
560            CONTINUE
	END IF
	CLOSE (iUnitM)

!           Now, apply the plate information for each node:

	DO 590 i = 1, numNod
		 theta = xNode(i)
		 phi = yNode(i)
		 nLat = 90.0D0 - theta * 57.2957795130823D0
		 eLon = phi * 57.2957795130823D0
		 IF (eLon < -180.0D0) eLon = eLon + 360.0D0
		 IF (eLon >  + 180.0D0) eLon = eLon - 360.0D0

!             Decide iPlate for this node,
!             by counting crossings of a line extending to South...

		 iPlate = iUnderPlate
		 IF (n_others > 0) THEN
			  DO 585 j = 1, n_others
				   nCross = 0
				   DO 580 k = 2, other_counts(j)
						lon1 = other_shapes(j, 1, k - 1)
						lon2 = other_shapes(j, 1, k)
						lat1 = other_shapes(j, 2, k - 1)
						lat2 = other_shapes(j, 2, k)
						IF (lon2 /= lon1) THEN
							 frac = (eLon - lon1) / (lon2 - lon1)
							 IF ((frac >= 0.0D0).AND. &
&                                   (frac < 1.0D0)) THEN
								  test = lat1 + frac * (lat2 - lat1)
								  IF (nLat > test) THEN
									   nCross = nCross + 1
								  END IF
							 END IF
						END IF
580                      CONTINUE
				   IF (MOD(nCross, 2) == 1) THEN
!                               odd number of crossings: inside
						iPlate = iOtherPlate(j)
						GO TO 586
				   END IF
585                 CONTINUE
586                 CONTINUE
		 END IF
!             Convert OMEGA(iPlate) to AFrica-fixed, and radians/second:
		 omegax = (omega(1, iPlate) - omega(1, iPAfri)) * 3.168809D-14
		 omegay = (omega(2, iPlate) - omega(2, iPAfri)) * 3.168809D-14
		 omegaz = (omega(3, iPlate) - omega(3, iPAfri)) * 3.168809D-14
!            Convert to length/second:
		 omegax = omegax * radius
		 omegay = omegay * radius
		 omegaz = omegaz * radius
!           Velocity = OMEGA x position:
		 xn = SIN(theta) * COS(phi)
		 yn = SIN(theta) * SIN(phi)
		 zn = COS(theta)
		 vx = omegay * zn - omegaz * yn
		 vy = omegaz * xn - omegax * zn
		 vz = omegax * yn - omegay * xn
!           Create unit +Theta and +Phi vectors in Cartesian:
		 thetax = COS(theta) * COS(phi)
		 thetay = COS(theta) * SIN(phi)
		 thetaz = -SIN(theta)
		 phix = -SIN(phi)
		 phiy = COS(phi)
		 phiz = 0.0D0
!           Find argument from dot products:
		 vTheta = vx * thetax + vy * thetay + vz * thetaz
		 vPhi = vx * phix + vy * phiy + vz * phiz
		 vM(1, i) = vTheta
		 vM(2, i) = vPhi
590       CONTINUE

ELSE
  write(ErrorMsg,'(A,I6)') "ILLEGAL INTEGER CODE FOR LOWER-MANTLE CONVECTION PATTERN (iConve): ",iConve
  call FatalError(ErrorMsg,ThID)
END IF

!      End of selection based on iConve;
!      Now apply velocity reference frame transformation from
!      AFrica-fixed to plate #iPVRef fixed:

!      Rotation of plate iPVRef wrt AFrica, in radians/second:
omegax = (omega(1, iPVRef) - omega(1, iPAfri)) * 3.168809D-14
omegay = (omega(2, iPVRef) - omega(2, iPAfri)) * 3.168809D-14
omegaz = (omega(3, iPVRef) - omega(3, iPAfri)) * 3.168809D-14
!      Convert to length/second:
omegax = omegax * radius
omegay = omegay * radius
omegaz = omegaz * radius

DO 2000 i = 1, numNod
!         Velocity of iPVRef wrt AFrica = OMEGA x position:
  theta = xNode(i)
  phi = yNode(i)
  xn = SIN(theta) * COS(phi)
  yn = SIN(theta) * SIN(phi)
  zn = COS(theta)
  vx = omegay * zn - omegaz * yn
  vy = omegaz * xn - omegax * zn
  vz = omegax * yn - omegay * xn
!         Create unit +Theta and +Phi vectors in Cartesian:
  thetax = COS(theta) * COS(phi)
  thetay = COS(theta) * SIN(phi)
  thetaz = -SIN(theta)
  phix = -SIN(phi)
  phiy = COS(phi)
  phiz = 0.0D0
!         Find argument from dot products:
  vTheta = vx * thetax + vy * thetay + vz * thetaz
  vPhi = vx * phix + vy * phiy + vz * phiz

!         Transform the velocity previously found in the
!         AFrica-fixed reference frame to one in the
!         iPVRef-fixed reference frame:
  vM(1, i) = vM(1, i) - vTheta
  vM(2, i) = vM(2, i) - vPhi

2000  CONTINUE

RETURN
END SUBROUTINE Convec

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
       write(ErrorMsg,'(A,I5/,A,I5,A,D14.4)') "COLATITUDE OF INTEGRATION POINT ",m, &
&                "OF ELEMENT ",i," IS OUT RANGE ",sitami
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
RETURN
END SUBROUTINE Deriv

SUBROUTINE Diamnd (aCreep, alphaT, bCreep, & ! input
&                    Biot, cCreep, dCreep, &
&                    eCreep, &
&                    e1, e2, fric, g, &
&                    geoth1, &
&                    geoth2, &
&                    geoth3, &
&                    geoth4, &
&                    pl0, pw0, &
&                    rhoBar, rhoH2O, sigHBi, &
&                    thick, temLim, &
&                    visMax, zOfTop, &
&                    pT1dE1, pT1dE2, &         ! output
&                    pT2dE1, pT2dE2, &
&                    pT1, pT2, zTran)

!      For one homogeneous layer (crust, *or* mantle lithosphere),
!      computes the vertical integral, through the layer, of
!      horizontal principal stresses (relative to the vertical stress);
!      reports these as pT1 (more negative) and pT2 (more positive).

!      Also reports zTran, the depth into the layer of the brittle/
!      ductile transition (greatest depth of earthquakes).

!      Finally, recommends layer partial derivitives
!            pT1dE1, pT1dE2, pT2dE1, pT2dE2
!      to be used in constructing "alpha" and tOfset (in -Viscos-),
!      according to strategy in pages 3973-3977 of Bird (1989).
!      In computing these, as in computing pT1 and pT2, the viscosity
!      limit visMax is applied to the average behavior of the whole
!      frictional layer, and again to the average behavior of the
!      whole creeping layer; it is not applied locally at each depth.

!      Necessary conditions when calling -Diamnd-:
!        -> horizontal principal strain-rates e1 and e2 not both zero;
!        -> e2 >= e1;
!        -> layer thickness "thick" is positive.

!      Note special kludge: if friction "fric" is >2.0D0, then this is
!      taken to be a signal that NO frictional layer is desired,
!      and that the whole layer should be power-law (or plastic, or
!      viscous-- whichever gives the least shear stress).

!      New version, May 5, 1998, by Peter Bird; intended to improve
!      the convergence behavior of all F-E programs which use it.
IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Arguments (*** all are scalars, even though
!      these same names may be arrays in other programs! ***):
REAL*8, INTENT(IN) :: aCreep, alphaT, bCreep, Biot, cCreep, dCreep, &                   ! input
&      eCreep, e1, e2, fric, g, &                                                         ! input
&      geoth1, geoth2, geoth3, geoth4, &                                                  ! input
&      pl0, pw0, &                                                                        ! input
&      rhoBar, rhoH2O, sigHBi, &                                                          ! input
&      thick, temLim, visMax, zOfTop                                                      ! input
REAL*8, INTENT(OUT) :: pT1, pT2, pT1dE1, pT1dE2, pT2dE1, pT2dE2, zTran                  ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Internal variables:
INTEGER n, nVStep
DOUBLE PRECISION secInv
REAL*8 angat2, angat3, angle, argume, &
&        delNeg, delPos, dSFdEV, &
&        dS1dE1, dS1dE2, dS2dE1, dS2dE2, &
&        dT1dE1, dT1dE2, dT2dE1, dT2dE2, dz, &
&        e1at1, e1at2, e1at3, e1at4, &
&        e2at1, e2at2, e2at3, e2at4, &
&        eSCrit, ez, &
&        frac, &
&        gamma, great, &
&        pH2O, &
&        r, rhoUse, &
&        sigma1, sigma2, s1Eff, s2Eff, s1rel, s2rel, &
&        sc0, sch, sc1, sf0, sfh, sf1, sTFric, sz, szEff, &
&        tau1, tau2, tecn, tecs, tect, tMean, tsfn, tsfs, tsft, &
&        t, t0, th, t1, &
&        vis, visDCr, visInf, visInt, visMin, visSHB, &
&        z, z0, zh, z1
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!      CHARACTERIZE THE STRAIN-RATE TENSOR:

ez = -(e1 + e2)
!     (Formula for vertical strain-rate ez comes from the
!      incompressibility of all permanent, anelastic strain types.)
secInv = -(e1 * e2 + e1 * ez + e2 * ez)
!     (One possible form for the second invariant of the matrix.)
!      Note that the double-precision is just to prevent underflows
!      from squaring small strain rates, not for precision.
visInf = 0.5D0 * aCreep * (2.0D0 * SQRT(secInv))**(eCreep - 1.0D0)
!      visInf is the viscosity for dislocation creep, lacking only
!      the exponential term; therefore, as a mathematical abstraction,
!      we can say that it is the viscosity at infinite temperature.

!      CHARACTERIZE THE CONTINUUM FRICTION:

sTFric = SIN(ATAN(fric))
gamma = (1.0D0 + sTFric) / (1.0D0 - sTFric)
!      Note: For thrusting, effective-sigma1h is effective-sigma1z
!            times gamma.  For normal faulting, effective-sigma2h
!            is effective-sigmaz/gamma.  For small "fric", gamma
!            is approximately equal to 1.+2.*fric

!      FIND THE BRITTLE/DUCTILE TRANSITION (zTran, measured from
!      the top of the layer):

!      In the thrusting quadrant (e1<0, e2<0) and in the normal-
!      faulting quadrant (e1>0, e2>0) the brittle/ductile transtion
!      is clear: it the greatest depth of frictional behavior
!      (possibly including earthquakes) on any fault, which is also
!      the greatest depth of frictional behavior on the most active
!      fault set.

!      However, in the strike-slip quadrant (e1<0, e2>0) the
!      transition is less clear.  I do not know of any empirical
!      field study which has determined how the transition depth
!      depends on (e1+e2) within the transtensional and transpressional
!      wedges of the strain-rate field.  Therefore, we have to choose
!      some simple rule.  The rule that the transition is at the
!      greatest depth of frictional behavior on any fault would
!      create two discontinuities (at the e1=0 line, where normal
!      faulting appears/dissapears; and at the e2=0 line, where
!      strike-slip faulting appears/dissapears).  Furthermore, the
!      transition depth near to these lines (on the deeper side) would
!      be defined by the less-active fault set, which asymptotically
!      becomes totally inactive as the line is approached!  If we
!      chose the alternate rule of taking the deepest frictional
!      behavior on the most active fault set, we would still have
!      two discontinuities, although at different places, both within
!      the strike-slip quadrant.  My F-E programs cannot converge well
!      when there is any discontinuity; therefore, I have chosen an
!      arbitrary rule which smooths the transition depth across each
!      of the transpressional and transtensional wedges, giving the
!      correct (unambiguous) depths on the lines e1=0, e1=-e2, and
!      e2=0.  In order to do this, I apply SIN(2*theta) smoothing to
!      both the frictional parameter dSFdEV and also to the creep
!      parameter eSCrit, and then compute the transition depth from
!      the combination of values.  (I do this instead of smoothing
!      the depth itself because I have no formula for the transition
!      depth on any of these three lines, and would have to locate
!      it by additional numerical searches.)

!      eSCrit is the shear strain rate (tensor type, =
!        0.5*(larger principal rate - smaller principal rate)
!        of the shear system which defines the transition
!        from the creep side (from below);
!      dSFdEV is the partial derivitive of the maximum shear
!        stress (on any plane) in the frictional domain
!        with respect to effective vertical stress
!       (vertical stress plus Biot times water pressure).

IF (e1 >= 0.0D0) THEN
	IF (e2 >= 0.0D0) THEN
!                Normal-normal; faster E2 dominates.
		 eSCrit = 0.50D0 * (e2 - ez)
		 dSFdEV = 0.50D0 * (1.0D0 - (1.0D0 / gamma))
	ELSE
!               (e1 >=0, e2 < 0)
!                e2 < e1? Should not happen!
       write(ErrorMsg,'(A,1P,D10.2,A,D10.2,A)') "ERROR in DIAMND: e1:",e1," > e2:",e2,")"
       call FatalError(ErrorMsg,ThID)
	END IF
ELSE
!          Note: (E1 < 0)
	IF (e2 >= 0.0D0) THEN
!                Note: (e1 < 0, E2 >= 0)
		 IF (ez >= 0.0D0) THEN
!                     Transpression (T/S).
!                     Enforce smooth transition in dSFdEV
!                     as the pure strike-slip line is approached.
!                    (This smoothing cannot be with visMax because
!                     zTran is not yet known; instead, use a smooth
!                     function of angle from origin of the
!                     strain-rate plane, varying over 45 degrees
!                     from the pure-strike-slip line e1=-e2
!                     to the pure-thrust line e2=0.)
			  tsft = 0.50D0 * (gamma - 1.0D0)
			  tsfs = sTFric
!                     Note: One might expect tsfs=fric, but check on
!                     a Mohr-circle diagram, remembering that the
!                     pure strike-slip condition is eZ==0 -->
!                     szzEff = 0.5 * (s1Eff + s2Eff).
!                     Also remember that the "SF" in dSFdEV is not the
!                     shear stress on the fault, but the maximum shear
!                     stress, because this is what creep will attack and
!                     lower first, at the brittle/ductile transition.
			  angle = ATan2F(e2, e1)
			  dSFdEV = tsfs + (tsft - tsfs) * SIN(2.0D0 * (angle - 2.3561945D0))

			  r = SQRT(e1**2 + e2**2)
			  tect = 1.0D0
			  tecs = 0.7071067D0
			  eSCrit = r * (tecs + (tect - tecs) * SIN(2.0D0 * (angle - 2.3561945D0)))
		 ELSE
!                     Note: (e1 < 0, e2 >= 0, eZ < 0)
!                     Transtension (N/S).
!                     Enforce smooth transition in dSFdEV
!                     as the pure strike-slip line is approached.
!                    (This smoothing cannot be with visMax because
!                     zTran is not yet known; instead, use a smooth
!                     function of angle from origin of the
!                     strain-rate plane, varying over 45 degrees
!                     from the pure-strike-slip line e1=-e2 to the
!                     pure-normal faulting line e1=0.)
			  tsfn = 0.5D0 * (1.0D0 - (1.0D0 / gamma))
			  tsfs = sTFric
!                     Note: One might expect tsfs=fric, but check on
!                     a Mohr-circle diagram, remembering that the
!                     pure strike-slip condition is ez==0 -->
!                     szzEff = 0.5 * (s1Eff + s2Eff).
!                     Also remember that the "SF" in dSFdEV is not the
!                     shear stress on the fault, but the maximum shear
!                     stress, because this is what creep will attack and
!                     lower first, at the brittle/ductile transition.
			  angle = ATan2F(e2, e1)
			  dSFdEV = tsfs + (tsfn - tsfs) * SIN(2.0D0 * (2.3561945D0 - angle))

			  r = SQRT(e1**2 + e2**2)
			  tecn = 1.0D0
			  tecs = 0.7071067D0
			  eSCrit = r * (tecs + (tecn - tecs) * SIN(2.0D0 * (2.3561945D0 - angle)))
		 END IF
	ELSE
!                Note: (e1 < 0, e2 < 0)
!                Thrust-thrust; faster (more negative) e1 dominates.
		 eSCrit = 0.5D0 * (ez - e1)
		 dSFdEV = 0.5D0 * (gamma - 1.0D0)
	END IF
END IF

!      Use eSCrit and dSFdEV to locate zTran (brittle/ductile trans.):

IF (fric > 2.0D0) THEN
!           Special kludge; no frictional layer is wanted
!          (for models with a purely power-law or linear-viscous
!           rheology, you specify an unrealistically high friction.
!           This makes the transition occur at the surface, and
!           below the surface, the friction value is irrelevant.)
	zTran = 0.
ELSE
!           Normal case; compute friction and creep at top and bottom:

	z0 = 0.0D0
	sf0 = dSFdEV * (pl0 - Biot * pw0)
	t0 = MIN(temLim, geoth1)
	argume = (bCreep + cCreep * zOfTop) / t0
!           Avoid overflow in EXP() by limiting the argument:
	argume = MAX(MIN(argume, 87.0D0), -87.0D0)
	sc0 = 2.0D0 * (visInf * eSCrit) * EXP(argume)
	sc0 = MIN(sc0, dCreep)

	z1 = thick
	tMean = geoth1 + &
&        0.50D0 * geoth2 * z1 + &
&      0.3330D0 * geoth3 * z1**2 + &
&       0.250D0 * geoth4 * z1**3
	rhoUse = rhoBar * (1.0D0 - alphaT * tMean)
	sf1 = sf0 + dSFdEV * (rhoUse - Biot * rhoH2O) * g * thick
	t1 = MIN(temLim, geoth1 + geoth2 * z1 + geoth3 * z1**2 + geoth4 * z1**3)
	argume = (bCreep + cCreep * (zOfTop + z1)) / t1
	argume = MAX(MIN(argume, 87.0D0), -87.0D0)
	sc1 = 2.0D0 * (visInf * eSCrit) * EXP(argume)
	sc1 = MIN(sc1, dCreep)
	sc1 = MAX(sc1, sigHBi)

!           Check if whole layer is frictional:
	IF (sc1 >= sf1) THEN
		 zTran = thick

!           Check if none of layer is frictional:
	ELSE IF (sc0 <= sf0) THEN
		 zTran = 0.0D0

	ELSE
!                Transition is within layer, between z0 and z1.
!                Use a binary-division search to bracket within
!                the nearest 1/128 of the layer (usually, within
!                0.5 km); then, finish with linear interpolation.
!                Note ASSUMPTION: T increases montonically with z!!!
!                Also note that linearity may fail if the
!                power-law/dCreep-limit transition falls into the
!                remaining interval; however, the error will be small.
		 DO 100 n = 1, 7
			  zh = 0.50D0 * (z0 + z1)
			  tMean = 0.50D0 * (t0 + t1)
			  rhoUse = rhoBar * (1.0D0 - alphaT * tMean)
			  sfh = sf0 + dSFdEV * (rhoUse - Biot * rhoH2O) * g * (zh - z0)
			  th = MIN(temLim, geoth1 + geoth2 * zh + geoth3 * zh**2 + &
&                              geoth4 * zh**3)
			  argume = (bCreep + cCreep * (zOfTop + zh)) / th
			  argume = MAX(MIN(argume, 87.0D0), -87.0D0)
			  sch = 2.0D0 * (visInf * eSCrit) * EXP(argume)
			  sch = MIN(sch, dCreep)
			  sch = MAX(sch, sigHBi)
			  IF (sch > sfh) THEN
!                          Transition is between zh and z1.
				   z0 = zh
				   sf0 = sfh
				   t0 = th
				   sc0 = sch
			  ELSE
!                          Transition is between z0 and zh.
				   z1 = zh
				   sf1 = sfh
				   t1 = th
				   sc1 = sch
			  END IF
100            CONTINUE
		 delNeg = sf0 - sc0
		 delPos = sf1 - sc1
		 frac = -delNeg / (delPos - delNeg)
		 IF ((frac < -0.01D0).OR.(frac > 1.01D0)) THEN
			  IF(Verbose) WRITE(iUnitVerb, "(' WARNING: Failure to bracket zTran', &
&                          ' within -Diamnd-.')")
		 END IF
		 frac = MIN(1.0D0, MAX(0.0D0, frac))
		 zTran = z0 + frac * (z1 - z0)
	END IF
END IF

!      SUM TAU (AND DERIVITIVES) OVER FRICTIONAL AND CREEP LAYERS:

!      Initialize sums over (up to) two layers:
!      -brittle layer at <= zTran from the top;
!      -creeping layer at > zTran from the top.
pT1 = 0.0D0
pT2 = 0.0D0
pT1dE1 = 0.0D0
pT1dE2 = 0.0D0
pT2dE1 = 0.0D0
pT2dE2 = 0.0D0

!      COMPUTE AND ADD STRENGTH OF FRICTIONAL PART OF LAYER:

IF (zTran > 0.0D0) THEN
!          Compute the effective vertical stress at the midpoint
!          of the frictional layer:
	tMean = geoth1 + &
&        0.5D0 * geoth2 * (zTran / 2.0D0) + &
&      0.333D0 * geoth3 * (zTran / 2.0D0)**2 + &
&       0.25D0 * geoth4 * (zTran / 2.0D0)**3
	rhoUse = rhoBar * (1.0D0 - alphaT * tMean)
	sz = -pl0 - rhoUse * g * zTran / 2.0D0
	pH2O = pw0 + rhoH2O * g * zTran / 2.0D0
	szEff = sz + Biot * pH2O

!          Compute effective horizontal principal stresses,
!          and their derivitives with respect to e1 and e2,
!          at the midpoint of the frictional layer, according
!          to the methods in Bird (1989), pages 3973-3977
!         (except, correcting the typos in the caption for
!          Figure 4):

!          Define the corner points of the diamond in the
!          ordered principal strain-rate plane:
	e1at1 = ((1.0D0 / gamma) - 1.0D0) * szEff / (6.0D0 * visMax)
	e2at1 = e1at1
	e1at2 = (1.0D0 - (1.0D0 / gamma)) * szEff / (6.0D0 * visMax)
	e2at2 = ((2.0D0 / gamma) - 2.0D0) * szEff / (6.0D0 * visMax)
	e1at3 = (2.0D0 * gamma - 2.0D0) * szEff / (6.0D0 * visMax)
	e2at3 = (1.0D0 - gamma) * szEff / (6.0D0 * visMax)
	e1at4 = (gamma - 1.0D0) * szEff / (6.0D0 * visMax)
	e2at4 = e1at4
	angat2 = ATan2F((e2 - e2at2), (e1 - e1at2))
	angat3 = ATan2F((e2 - e2at3), (e1 - e1at3))

!          Select proper segment of diagram and assign effective
!          principal stresses.
!          Also, begin definition of strategic stiffnesses
!          dS1dE1, dS1dE2, dS2dE1, and dS2dE2, by computing
!          stiffness required to give warning of local cliffs.
!          Afterward, basic minimum stiffness required to avoid
!          singularity of stiffness matrix will be imposed with
!          a formula common to all regions.
	IF (e1 > e1at1) THEN
!                Region N/N: two conjugate sets of normal faults
		 s1Eff = szEff / gamma
		 s2Eff = s1Eff

		 dS1dE1 = (0.50D0 * ((1.0D0 / gamma) - 1.0D0) * szEff) / e1
		 dS1dE2 = 0.0D0
		 dS2dE1 = 0.0D0
		 dS2dE2 = (0.50D0 * ((1.0D0 / gamma) - 1.0D0) * szEff) / e2
	ELSE IF ((e1 >= e1at2).AND.(angat2 > ATan2F((e2at1 - e2at2), (e1at1 - e1at2)))) THEN
!                Region N: single conjugate set of normal faults
		 s2Eff = szEff / gamma
		 frac = (e1 - e1at1) / (e1at2 - e1at1)
!                fraction increases in -e1 direction, from point 1 -> 2
		 s1Eff = szEff * ((1 / gamma) + frac * (1.0D0 - (1.0D0 / gamma)))

		 dS1dE1 = 4.0D0 * visMax
		 dS1dE2 = 0.0D0
		 dS2dE1 = 0.0D0
		 dS2dE2 = 0.0D0
	ELSE IF ((angat2 <= 1.9635D0).AND.(angat2 >= 1.5707D0)) THEN
!                Region N/S: transtension, dominantly normal.
		 s1Eff = szEff
		 s2Eff = szEff / gamma

		 dS1dE1 = (0.5D0 * ((1.0D0 - 1.0D0 / gamma)) * szEff) / e1
		 dS1dE2 = 0.0D0
		 dS2dE1 = 0.0D0
		 dS2dE2 = 0.0D0
	ELSE IF ((angat2 <= 2.3562D0).AND.(angat2 >= 1.9635D0)) THEN
!                Region S/N: transtension, dominantly strike-slip.
		 s1Eff = szEff
		 s2Eff = szEff / gamma

!               "great" is the value of dS1dE1 in region S:
		 great = 6.0D0 * visMax * (gamma - 1.0D0) / (gamma - (1.0D0 / gamma))
!               "frac" is also defined exactly as in S, so here it
!                   will be negative:
		 frac = ((e1 + e2) - (e1at2 + e2at2)) / &
&                  ((e1at3 + e2at3) - (e1at2 + e2at2))
!                Reduce all derivitives according to distance:
		 great = great * (-0.50D0) / (frac - 0.50D0)
!                Pattern of derivitives is the same as in S:
		 dS1dE1 = great
		 dS1dE2 = dS1dE1
		 dS2dE1 = dS1dE1 / gamma
		 dS2dE2 = dS2dE1
	ELSE IF ((angat3 <= 2.3562D0).AND.(angat3 >= ATan2F((e2at2 - e2at3), (e1at2 - e1at3)))) THEN
!                Region S: single set of conjugate strike-slip faults
		 frac = ((e1 + e2) - (e1at2 + e2at2)) / &
&                  ((e1at3 + e2at3) - (e1at2 + e2at2))
!               "frac" increases across band from the S/N (point 2) side
!                toward the S/T (point 3) side; contours of "frac" are
!                parallel to the band sides, not normal to the diamond.
		 s1Eff = szEff * (1.0D0 + frac * (gamma - 1.0D0))
		 s2Eff = szEff * ((1.0D0 / gamma) + frac * (1.0D0 - (1.0D0 / gamma)))
!                Notes: The equation of this line is s2Eff=s1Eff/gamma.
!                I used algebra to check (1998.04.21) that the
!                pure strike-slip stress (s1Eff,s2Eff)=
!                szzEff*(1.+sTFric,1.-sTFric) correctly falls on
!                this line, at the correct point (e1= -e2).

		 dS1dE1 = 6.0D0 * visMax * (gamma - 1.0D0) / (gamma - (1.0D0 / gamma))
		 dS1dE2 = dS1dE1
		 dS2dE1 = dS1dE1 / gamma
		 dS2dE2 = dS2dE1
	ELSE IF ((angat3 <= 2.7489D0).AND.(angat3 >= 2.3562D0)) THEN
!                Region S/T: transpression; strike-slip dominant.
		 s1Eff = szEff * gamma
		 s2Eff = szEff

!               "great" is the value of dS1dE1 in region S:
		 great = 6.0D0 * visMax * (gamma - 1.0D0) / (gamma - (1.0D0 / gamma))
!               "frac" is also defined exactly as in S, so here it
!                will be greater than one:
		 frac = ((e1 + e2) - (e1at2 + e2at2)) / &
&                  ((e1at3 + e2at3) - (e1at2 + e2at2))
!                Reduce all derivitives according to distance:
		 great = great * (0.50D0) / (frac - 0.50D0)
!                Pattern of derivitives is the same as in S:
		 dS1dE1 = great
		 dS1dE2 = dS1dE1
		 dS2dE1 = dS1dE1 / gamma
		 dS2dE2 = dS2dE1
	ELSE IF ((e2 >= e2at3).AND.(angat3 >= 2.7489D0)) THEN
!                Region T/S: transpression; thrusting dominant.
		 s1Eff = szEff * gamma
		 s2Eff = szEff

		 dS1dE1 = 0.0D0
		 dS1dE2 = 0.0D0
		 dS2dE1 = 0.0D0
		 dS2dE2 = (0.50D0 * (1.0D0 - gamma) * szEff) / e2
	ELSE IF ((e2 >= e2at4).AND.(angat3 <= ATan2F((e2at4 - e2at3), (e1at4 - e1at3)))) THEN
!                Region T: single conjugate thrust fault set.
		 s1Eff = szEff * gamma
		 frac = (e2 - e2at3) / (e2at4 - e2at3)
!               "frac" increases in the -e2 direction across the band.
		 s2Eff = szEff * (1.0D0 + frac * (gamma - 1.0D0))

		 dS1dE1 = 0.0D0
		 dS1dE2 = 0.0D0
		 dS2dE1 = 0.0D0
		 dS2dE2 = 4.0D0 * visMax
	ELSE IF (e2 <= e2at4) THEN
!                Region T/T: Two set of conjugate thrust faults.
		 s1Eff = szEff * gamma
		 s2Eff = s1Eff

		 dS1dE1 = (0.50D0 * (gamma - 1.0D0) * szEff) / e1
		 dS1dE2 = 0.0D0
		 dS2dE1 = 0.0D0
		 dS2dE2 = (0.50D0 * (gamma - 1.0D0) * szEff) / e2
	ELSE
!                Region V: linear viscosity
!                Note that equations are now for sigma1,2 and no
!                longer for s1Eff and s2Eff.  However, we can
!                easily compute both:
		 sigma1 = sz + visMax * (4.0D0 * e1 + 2.0D0 * e2)
		 sigma2 = sz + visMax * (2.0D0 * e1 + 4.0D0 * e2)
		 s1Eff = sigma1 + Biot * pH2O
		 s2Eff = sigma2 + Biot * pH2O

		 dS1dE1 = 0.0D0
		 dS1dE2 = 0.0D0
		 dS2dE1 = 0.0D0
		 dS2dE2 = 0.0D0
	END IF

!          Regardless of region, be sure that stiffnesses do
!          not fall below those which represent a minimum
!          effective viscosity-- one based on the weakest of
!          the active fault sets.  This is to guaruntee that
!          the linear system will not have any zero eigenvalues,
!          even if a creeping layer does not exist.
	visMin = visMax
	IF ((e1 < 0.0D0).AND.(e2 > 0.0D0)) THEN
!                strike-slip faults are active
		 visMin = MIN(visMin, 0.50D0 * (s2Eff - s1Eff) / (e2 - e1))
	END IF
	IF ((e1 < 0.0D0).AND.(ez > 0.0D0)) THEN
!                thrust faults are active
		 visMin = MIN(visMin, 0.50D0 * (szEff - s1Eff) / (ez - e1))
	END IF
	IF ((e2 > 0.0D0).AND.(ez < 0.0D0)) THEN
!                normal faults are active
		 visMin = MIN(visMin, 0.50D0 * (s2Eff - szEff) / (e2 - ez))
	END IF
	dS1dE1 = dS1dE1 + 4.0D0 * visMin
	dS1dE2 = dS1dE2 + 2.0D0 * visMin
	dS2dE1 = dS2dE1 + 2.0D0 * visMin
	dS2dE2 = dS2dE2 + 4.0D0 * visMin

!          Convert effective principal stresses at the midpoint
!          of the frictional layer into total principal stresses:
	sigma1 = s1Eff - Biot * pH2O
	sigma2 = s2Eff - Biot * pH2O
!         (Note that correcting S1 and S2 by a constant does not
!          affect the values of any of the 4 derivitives dS1dE1, ..., dS2dE2.)

!          Convert total principal stresses at the midpoint of
!          the frictional layer into relative principal stresses
!         (relative to the total vertical stress, that is):
	s1rel = sigma1 - sz
	s2rel = sigma2 - sz
!         (Note that correcting S1 and S2 by a constant does not
!          affect the values of any of the 4 derivitives dS1dE1, ..., dS2dE2.)

!          Convert values at midpoint of frictional layer to
!          integrals over the frictional layer:
	tau1 = s1rel * zTran
	tau2 = s2rel * zTran
	dT1dE1 = dS1dE1 * zTran
	dT1dE2 = dS1dE2 * zTran
	dT2dE1 = dS2dE1 * zTran
	dT2dE2 = dS2dE2 * zTran

!          Add integrals over frictional layer to layer totals:
	pT1 = pT1 + tau1
	pT2 = pT2 + tau2
	pT1dE1 = pT1dE1 + dT1dE1
	pT1dE2 = pT1dE2 + dT1dE2
	pT2dE1 = pT2dE1 + dT2dE1
	pT2dE2 = pT2dE2 + dT2dE2
END IF
!     (IF the frictional layer thickness zTran > 0)

!      COMPUTE AND ADD STRENGTH OF CREEPING PART OF LAYER:

IF (zTran < thick) THEN

!           Precompute the maximum viscosity limit imposed by the
!           requirement that creep shear stress never exceeds
!           dCreep on any plane:
	visDCr = dCreep / (MAX(e1, e2, ez) - MIN(e1, e2, ez))

!           Precompute the lower viscosity limit imposed by the
!           requirement that creep shear stress does not
!           fall below sigHBI:
	visSHB = sigHBi / (MAX(e1, e2, ez) - MIN(e1, e2, ez))

!           Compute the vertical integral of viscosity,
!           observing the local limit visDCr, and terminating
!           the integral if creep shear stress falls below
!           sigHBI (because then we are in a horizontally-
!           sheared boundary layer which does not contribute
!           anything to plate strength):

	nVStep = 50
	dz = (thick - zTran) / nVStep

	visInt = 0.0D0
	DO 200 n = 0, nVStep
		 z = zTran + n * dz
!                Note that z is measured from top of layer
!               (upper surface of hard crust, or Moho) and
!                may not be absolute depth.
		 t = geoth1 + geoth2 * z + geoth3 * z**2 + geoth4 * z**3
		 t = MIN(t, temLim)
		 argume = (bCreep + cCreep * (zOfTop + z)) / t
!                Prevent over/underflow in EXP() by limiting the argument:
		 argume = MAX(MIN(argume, 87.0D0), -87.0D0)
		 vis = visInf * EXP(argume)
		 vis = MIN(vis, visDCr)
		 IF ((n == 0).OR.(n == nVStep)) THEN
			  frac = 0.50D0
		 ELSE
			  frac = 1.0D0
		 END IF
		 IF (vis < visSHB) GO TO 201
		 visInt = visInt + frac * vis * dz
200       CONTINUE
201       CONTINUE

!           Limit the mean viscosity of the creeping layer to
!           be no more than visMax:
	visInt = MIN(visInt, visMax * (thick - zTran))

	tau1 = 4.0D0 * visInt * e1 + 2.0D0 * visInt * e2
	tau2 = 2.0D0 * visInt * e1 + 4.0D0 * visInt * e2
!           Note that these principal values of tau (the two
!           horizontal principal values, contributed by the
!           creeping layer only) are relative to tauzz, which
!           is the vertical integral of the vertical stress
!           anomaly through the creeping layer.

	dT1dE1 = 4.0D0 * visInt
	dT1dE2 = 2.0D0 * visInt
	dT2dE1 = 2.0D0 * visInt
	dT2dE2 = 4.0D0 * visInt

!          Add integrals over creeping layer to layer totals:
	pT1 = pT1 + tau1
	pT2 = pT2 + tau2
	pT1dE1 = pT1dE1 + dT1dE1
	pT1dE2 = pT1dE2 + dT1dE2
	pT2dE1 = pT2dE1 + dT2dE1
	pT2dE2 = pT2dE2 + dT2dE2
END IF
!     (IF the creeping layer thickness (thick - zTran) > 0)

RETURN
END SUBROUTINE Diamnd

SUBROUTINE Downer (brief, fDip, iUnitT, mxBn, mxFEl, mxNode, & ! input
&                    nFl, nodeF, numNod, slide, &
&                    xNode, yNode, &
&                    nCond, nodCon, &                            ! output
&                    checkN)                                     ! work

!   Surveys faults for dips less than "slide" (in radians), and
!   lists the footwall nodes as needing boundary conditions.
!  (This routine is only called for whole-Earth models.)

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
LOGICAL, INTENT(IN) :: brief                                                           ! input
REAL*8, INTENT(IN) :: fDip                                                             ! input
INTEGER, INTENT(IN) :: iUnitT, mxBn, mxFEl, mxNode, nFl, nodeF, numNod                 ! input
REAL*8, INTENT(IN) :: slide, xNode, yNode                                              ! input
INTEGER, INTENT(OUT) :: nCond, nodCon                                                  ! output
LOGICAL checkN                                                          ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, n
REAL*8 theLat, theLon, tops
DIMENSION checkN(mxNode), fDip(2, mxFEl), nodCon(mxBn), &
&           nodeF(4, mxFEl), xNode(mxNode), yNode(mxNode)

tops = 3.14159265358979D0 - slide

DO 10 i = 1, numNod
	checkN(i) = .FALSE.
10  CONTINUE

DO 100 i = 1, nFl
	IF (nodeF(1, i) /= nodeF(4, i)) THEN
		 IF (fDip(1, i) <= slide) THEN
			  checkN(nodeF(4, i)) = .TRUE.
		 ELSE IF (fDip(1, i) >= tops) THEN
			  checkN(nodeF(1, i)) = .TRUE.
		 END IF
	END IF
	IF (nodeF(2, i) /= nodeF(3, i)) THEN
		 IF (fDip(2, i) <= slide) THEN
			  checkN(nodeF(3, i)) = .TRUE.
		 ELSE IF (fDip(2, i) >= tops) THEN
			  checkN(nodeF(2, i)) = .TRUE.
		 END IF
	END IF
100  CONTINUE

nCond = 0
DO 200 i = 1, numNod
	IF (checkN(i)) THEN
		 IF (nCond < mxBn) THEN
			  nCond = nCond + 1
			  nodCon(nCond) = i
		 ELSE
           write(ErrorMsg,'(A)') "Increase the constant in the formula for mxBn, and recompile."
           call FatalError(ErrorMsg,ThID)
		 END IF
	END IF
200  CONTINUE

IF (.NOT.brief .AND. Verbose) THEN
	WRITE(iUnitVerb, 880)
880       FORMAT(/ /' Here follows a list, in consecutive order,'/ &
&                ' of the nodes in the footwalls of '/ &
&                ' SUBduction zones; these nodes require boundary', &
&                ' conditions:'/'    BC#  Node  Latitude', &
&                ' Longitude')
	DO 890 i = 1, nCond
		 n = nodCon(i)
		 theLat = 90.0D0 - 57.2957795130823D0 * xNode(n)
		 theLon = 57.2957795130823D0 * yNode(n)
		 WRITE(iUnitVerb, 882) i, n, theLat, theLon
882            FORMAT(' ',2I6,2F10.2)
890       CONTINUE
END IF
RETURN
END SUBROUTINE Downer

SUBROUTINE EdgeVs (fDip, iPVRef, iUnitD, iUnitT, mxBn, mxNode, & ! input
&                    mxFEl, names, nCond, nFl, nodCon, nodeF, nPlate, &
&                    omega, radius, slide, sphere, xNode, yNode, &
&                    iCond, vBCArg, vBCMag, &                      ! modify
&                    iEdge, r2Edge, xEdge, yEdge)                  ! work

!      Supplies velocities (vBCMag) and arguments (vBCArg)
!     (measured counterclockwise in radians from +X = +Theta = +South)
!      for all the boundary nodes listed in nodCon
!      which have iCond = 3 or 4.

!      It does this by consulting a data table of
!      rotation-rate poles with respect to the PAcific plate, from:

!      Bird, P., An updated digital model of plate boundaries,
!      Geochemistry Geophysics Geosystems, 4(3), 1027,
!      doi:10.1029/2001GC000252, 2003.

!      The output reference frame is one where plate #iPVRef is fixed.

!      To figure out which plate applies, a dataset of digitized
!      plate boundary lines (*.dig) is read on device iUnitD.
!     (As of 2005.07, the file to use is "PB2002_boundaries.dig".)

!      Thest line segments must be labeled with 2-letter plate
!      "names" and an intervening boundary symbol, e.g.
!          Column 1 ---->PA\NA [or] PA)NA
!      so that the first plate named is on the left as one
!      proceeds along the digitised boundary.
!      Boundary symbols "/" [or "("] versus "\" [or ")"] must be used
!      for subduction boundaries to show polarity.  Ridges and/or
!      transforms are marked by "-".

!      If no digitized boundary is found near the node (within
!      0.1 radians) or if the closest digitised boundary has a plate
!      name which is not found in this table,
!      then iCond is modified to make the node free.

!CC    IF (sphere) THEN ! a global model
!CC         look-up position is that of a fault midpoint.
!CC         boundary v is that of the subducting plate.
!CC    ELSE (a local model)
!CC         look-up position is the node position.
!CC         Boundary v is that of the plate on the right, looking along
!CC            digitised line (i.e., If we assume that the edge
!CC            of the model is surrounded with fault elements,
!CC            then it should be digitised counterclockwise.
!CC           (If model is not surrounded by faults, then using
!CC            boundary conditons of type 5 would be more appropriate.)
!CC    END IF
!-------------------------------------------------------
IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: fDip                                                             ! input
INTEGER, INTENT(IN) :: iPVRef, iUnitD, iUnitT, mxBn, mxNode, mxFEl                     ! input
CHARACTER*2, INTENT(IN) :: names                                                       ! input
INTEGER, INTENT(IN) :: nCond, nFl, nodCon, nodeF, nPlate                               ! input
REAL*8, INTENT(IN) :: omega, radius, slide                                             ! input
LOGICAL, INTENT(IN) :: sphere                                                          ! input
REAL*8, INTENT(IN) :: xNode, yNode                                                     ! input
INTEGER, INTENT(INOUT) :: iCond                                                        ! modify
REAL*8, INTENT(INOUT) :: vBCArg, vBCMag                                                ! modify
INTEGER iEdge                                                           ! work
REAL*8 r2Edge, xEdge, yEdge                                             ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CHARACTER*1 left1, left2, right1, right2, which
CHARACTER*2 first, second
CHARACTER*3 stars
INTEGER i, ios, iPlate, k, n, n1, n2, n3, n4, nRead
LOGICAL readIt
REAL*8 a1, a2, am, b1, b2, bm, equpar, g1, g2, gm, &
	& omegax, omegay, omegaz, &
	& phi, phix, phiy, phiz, pLat, pLon, r2, size, &
	& theta, thetax, thetay, thetaz, vPhi, vTheta, vx, vy, vz, &
	& x, x1, x2, xn, xp, y, y1, y2, yn, yp, z, zn, zp
DIMENSION fDip(2, mxFEl), iCond  (mxBn), iEdge  (mxBn), &
&           nodCon (mxBn), nodeF(4, mxFEl), r2Edge (mxBn), &
&           vBCArg (mxBn), vBCMag (mxBn), &
&           xEdge  (mxBn), yEdge  (mxBn), &
&           xNode (mxNode), yNode (mxNode)
DIMENSION names (nPlate), omega (3, nPlate)
DATA left1 / '\' / , left2 / ')' / , right1 / '/' / , right2 / '(' /
!----------------------------------------------------------------
readIt = .FALSE.
DO 90 i = 1, nCond
	n = nodCon(i)
	IF ((iCond(i) == 3).OR.(iCond(i) == 4)) THEN
		 readIt = .TRUE.
		 iEdge(i) = 0
		 r2Edge(i) = 9.99D29
		 IF (sphere) THEN
			  DO 50 k = 1, nFl
				   n1 = nodeF(1, k)
				   n2 = nodeF(2, k)
				   n3 = nodeF(3, k)
				   n4 = nodeF(4, k)
				   IF (((n == n1).AND.(fDip(1, k) >= (3.14159265358979D0 - slide))) .OR. &
			  &        ((n == n2).AND.(fDip(2, k) >= (3.14159265358979D0 - slide))) .OR. &
			  &        ((n == n3).AND.(fDip(2, k) <= slide)) .OR. &
			  &        ((n == n4).AND.(fDip(1, k) <= slide))) THEN
						x1 = xNode(n1)
						x2 = xNode(n2)
						y1 = yNode(n1)
						y2 = yNode(n2)
						a1 = COS(y1) * SIN(x1)
						a2 = COS(y2) * SIN(x2)
						b1 = SIN(y1) * SIN(x1)
						b2 = SIN(y2) * SIN(x2)
						g1 = COS(x1)
						g2 = COS(x2)
						am = 0.50D0 * (a1 + a2)
						bm = 0.50D0 * (b1 + b2)
						gm = 0.50D0 * (g1 + g2)
						size = SQRT(am * am + bm * bm + gm * gm)
						am = am / size
						bm = bm / size
						gm = gm / size
						equpar = SQRT(am * am + bm * bm)
						theta = ATan2F(equpar, gm)
						phi = ATan2F(bm, am)
						xEdge(i) = theta
						yEdge(i) = phi
						GO TO 55
				   END IF
50                 CONTINUE
                   write(ErrorMsg,'(A/,A,I6,A/,A)') "ERROR IN SUBPROGRAM -EdgeVs-:", &
&                       "SUBPROGRAM -Downer- PLACED NODE ",n," ON LIST OF SUBDUCTING-SLAB NODES,", &
&                       "BUT NO SUCH FAULT ELEMENT FOUND."
                   call FatalError(ErrorMsg,ThID)
55                 CONTINUE
		 ELSE
			  xEdge(i) = xNode(nodCon(i))
			  yEdge(i) = yNode(nodCon(i))
		 END IF
	END IF
90  CONTINUE
IF (.NOT.readIt) RETURN

IF(Verbose) WRITE(iUnitVerb, 2) iUnitD
2  FORMAT (/' Attempting to read plate BOUNDARIES from unit ',I3/)
nRead = 0
100  READ (iUnitD, 101, END = 201, IOSTAT = ios) first, which, second
101       FORMAT (A2,A1,A2)
	IF ((nRead == 0).AND.(ios /= 0)) THEN
      write(ErrorMsg,'(A)') "ERROR: File not found, or file empty."
      call FatalError(ErrorMsg,ThID)
	END IF
	nRead = nRead + 1
	IF (sphere) THEN
!                Use plate which is subducting, if any:
		 IF ((which == left1).OR.(which == left2)) THEN
			  iPlate = 0
			  DO 110 k = 1, nPlate
				   IF (first == names(k)) iPlate = k
110                 CONTINUE
		 ELSE IF ((which == right1).OR. &
&                    (which == right2)) THEN
			  iPlate = 0
			  DO 120 k = 1, nPlate
				   IF (second == names(k)) iPlate = k
120                 CONTINUE
		 ELSE
			  iPlate = 0
		 END IF
	ELSE

!                Local model; use plate on right
		 iPlate = 0
		 DO 130 k = 1, nPlate
			  IF (second == names(k)) iPlate = k
130            CONTINUE
	END IF
140       READ (iUnitD, 141, END = 201) stars ! CHARACTER*3
141            FORMAT (A3)
		 IF (stars == '***') GO TO 100
		 BACKSPACE iUnitD
		 READ (iUnitD, *) pLon, pLat
		 IF (iPlate > 0) THEN
			  pLon = pLon * 0.0174532925199433D0
			  pLat = pLat * 0.0174532925199433D0
			  x = COS(pLon) * COS(pLat)
			  y = SIN(pLon) * COS(pLat)
			  z = SIN(pLat)
			  DO 150 i = 1, nCond
				   IF ((iCond(i) == 3).OR.(iCond(i) == 4)) THEN
						theta = xEdge(i)
						phi = yEdge(i)
						xp = COS(phi) * SIN(theta)
						yp = SIN(phi) * SIN(theta)
						zp = COS(theta)
						r2 = (x - xp)**2 + (y - yp)**2 + (z - zp)**2
						IF (r2 <= r2Edge(i)) THEN
							 r2Edge(i) = r2
							 iEdge(i) = iPlate
						END IF
				   END IF
150                 CONTINUE
		 END IF
	GO TO 140

201  DO 300 i = 1, nCond
	IF ((iCond(i) == 3).OR.(iCond(i) == 4)) THEN
		 iPlate = iEdge(i)
		 IF ((iPlate == iPVRef).AND.(iCond(i) == 3)) THEN
!                  problem: subduction direction is undefined!
			  iCond(i) = 4
			  vBCMag(i) = 0.0D0
			  vBCArg(i) = 0.0D0
!                  N.B. A plate never moves in its own reference frame!
			  IF(Verbose) WRITE (iUnitVerb, 205) n, pLon, pLat
205                 FORMAT (/' SUBDUCTION DIRECTION IS UNDEFINED', &
&                         ' FOR BOUNDARY NODE ',I6/ &
&                        '    AT ',F8.3,' EAST, ',F7.3,' NORTH: ', &
&                             'NODE WILL BE GIVEN A TYPE-4 BC.')
		 ELSE IF ((iPlate == 0).OR.(r2Edge(i) > 0.01D0)) THEN
!                  plate not identified, or trench very far away:
			  iCond(i) = 0
			  n = nodCon(i)
			  theta = xNode(n)
			  phi = yNode(n)
			  pLon = 57.2957795130823D0 * phi
			  pLat = 90.0D0 - 57.2957795130823D0 * theta
			  IF(Verbose) WRITE (iUnitVerb, 209) n, pLon, pLat
209                 FORMAT (/' NO RECOGNIZABLE DIGITISED SUBDUCTION', &
&                         ' ZONE PASSES BOUNDARY NODE ',I6/ &
&                        '    AT ',F8.3,' EAST, ',F7.3,' NORTH: ', &
&                             'NODE WILL BE FREE.')
		 ELSE
!                  normal cases:
		   IF (iPlate == iPVRef) THEN
			  vBCMag(i) = 0.0D0
			  vBCArg(i) = 0.0D0
		   ELSE
!                     Convert to iPVRef-fixed, and radians/second:
			  omegax = (omega(1, iPlate) - &
&                        omega(1, iPVRef)) * 3.168809D-14
			  omegay = (omega(2, iPlate) - &
&                        omega(2, iPVRef)) * 3.168809D-14
			  omegaz = (omega(3, iPlate) - &
&                        omega(3, iPVRef)) * 3.168809D-14
!                     Convert to length/second:
			  omegax = omegax * radius
			  omegay = omegay * radius
			  omegaz = omegaz * radius
			  n = nodCon(i)
			  theta = xNode(n)
			  phi = yNode(n)
			  xn = COS(phi) * SIN(theta)
			  yn = SIN(phi) * SIN(theta)
			  zn = COS(theta)
!                     Velocity = OMEGA x position:
			  vx = omegay * zn - omegaz * yn
			  vy = omegaz * xn - omegax * zn
			  vz = omegax * yn - omegay * xn
			  vBCMag(i) = SQRT(vx**2 + vy**2 + vz**2)
!                     Create unit +Theta and +Phi vectors in Cartesian:
			  thetax = COS(theta) * COS(phi)
			  thetay = COS(theta) * SIN(phi)
			  thetaz = -SIN(theta)
			  phix = -SIN(phi)
			  phiy = COS(phi)
			  phiz = 0.0D0
!                     Find argument from dot products:
			  vTheta = vx * thetax + vy * thetay + vz * thetaz
			  vPhi = vx * phix + vy * phiy + vz * phiz
			  vBCArg(i) = ATan2F(vPhi, vTheta)
		   END IF
		 END IF
	END IF
300  CONTINUE
RETURN
END SUBROUTINE EdgeVs

SUBROUTINE EDot (dXS, dYS, & ! input
&                  fPSfer, mxEl, &
&                  mxNode, nodes, numEl, radius, sita, v, &
&                  eRate)      ! output

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

SUBROUTINE ElUvec (n1, n2, n3, numNod, xNode, yNode, & ! input
&                    phiM, thetaM, uvecM)                ! output

!   Computes all necessary unit vectors at the integration points
!   within one triangular continuum element:
!      uvecM = position vector
!      thetaM= +theta unit vector at that site
!      phiM  = +phi unit vector at that site

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: n1, n2, n3, numNod                                              ! input
REAL*8, INTENT(IN) :: xNode, yNode                                                     ! input
REAL*8, INTENT(OUT) :: phiM, thetaM, uvecM                                             ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER j, m, n
REAL*8 equat, length, tequat, uvecN, x, y
DIMENSION xNode(numNod), yNode(numNod)
DIMENSION phiM(3, 7), thetaM(3, 7), uvecM(3, 7)
DIMENSION uvecN(3, 3)

!      Named COMMON blocks hold the fixed values of the positions,
!      weights, and nodal function values at the integration points
!      in the elements (triangular elements in BLOCK  DATA  BD1,
!      and fault elements in BLOCK  DATA  BD2).
!      Entries corresponding to BD1:
DOUBLE PRECISION points
COMMON / S1S2S3 / points
DIMENSION points(3, 7)

!      Find uvecN = uvec's of corner nodes:

DO 10 j = 1, 3
	IF (j == 1) THEN
		 n = n1
	ELSE IF (j == 2) THEN
		 n = n2
	ELSE
		 n = n3
	END IF
	x = xNode(n)
	y = yNode(n)
	uvecN(1, j) = SIN(x) * COS(y)
	uvecN(2, j) = SIN(x) * SIN(y)
	uvecN(3, j) = COS(x)
10  CONTINUE

!      Create each of 7 integration points:

DO 100 m = 1, 7
!           Rough linear interpolation:
	uvecM(1, m) = points(1, m) * uvecN(1, 1) + points(2, m) * uvecN(1, 2) + &
&                    points(3, m) * uvecN(1, 3)
	uvecM(2, m) = points(1, m) * uvecN(2, 1) + points(2, m) * uvecN(2, 2) + &
&                    points(3, m) * uvecN(2, 3)
	uvecM(3, m) = points(1, m) * uvecN(3, 1) + points(2, m) * uvecN(3, 2) + &
&                    points(3, m) * uvecN(3, 3)
!           Normalization:
	length = SQRT(uvecM(1, m)**2 + uvecM(2, m)**2 + uvecM(3, m)**2)
	uvecM(1, m) = uvecM(1, m) / length
	uvecM(2, m) = uvecM(2, m) / length
	uvecM(3, m) = uvecM(3, m) / length

!           Unit vectors at this site (NOT a pole):

	phiM(1, m) = -uvecM(2, m)
	phiM(2, m) = uvecM(1, m)
	equat = SQRT(phiM(1, m)**2 + phiM(2, m)**2)
	phiM(1, m) = phiM(1, m) / equat
	phiM(2, m) = phiM(2, m) / equat
	phiM(3, m) = 0.00D0
	tequat = uvecM(3, m)
	thetaM(3, m) = -equat
	thetaM(1, m) = tequat * uvecM(1, m) / equat
	thetaM(2, m) = tequat * uvecM(2, m) / equat
	length = SQRT(thetaM(1, m)**2 + thetaM(2, m)**2 + &
&                    thetaM(3, m)**2)
	thetaM(1, m) = thetaM(1, m) / length
	thetaM(2, m) = thetaM(2, m) / length
	thetaM(3, m) = thetaM(3, m) / length
100  CONTINUE
RETURN
END SUBROUTINE ElUvec

SUBROUTINE Euler (namTag, node, &  ! input
&                   iPVRef, names, nPlate, omega, &
&                   iUnitT, radius, &
&                   mxNode, xNode, yNode, &
&                   vAz, vMag)       ! output

!  Given a 2-letter plate identifier (namTag),
!  finds this plate in the table names(nPlate) and
!  looks up its Euler rotation vector in omega(3, nPlate).
!  Then converts the rotation to the reference frame of
!  plate # iPVRef (in the same tables).
!  Finally, computes the velocity vector
! (using the theta and phi of node #node in the lists
!  xNode(mxNode) and yNode(mxNode))
!  and expresses it as:
!  vAz = azimuth clockwise from North, in degrees;
!  vMag = magnitude in SI (m/s).

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CHARACTER(2),                      INTENT(IN) :: namTag
INTEGER,                           INTENT(IN) :: node
INTEGER,                           INTENT(IN) :: iPVRef
INTEGER,                           INTENT(IN) :: nPlate
CHARACTER(2), DIMENSION(nPlate),   INTENT(IN) :: names
REAL*8,       DIMENSION(3, nPlate),INTENT(IN) :: omega
INTEGER,                           INTENT(IN) :: iUnitT
REAL*8,                            INTENT(IN) :: radius
INTEGER,                           INTENT(IN) :: mxNode
REAL*8,       DIMENSION(mxNode),   INTENT(IN) :: xNode
REAL*8,       DIMENSION(mxNode),   INTENT(IN) :: yNode
REAL*8,                            INTENT(OUT):: vAz, vMag
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER :: i, jPlate
REAL*8 :: vTheta, vPhi
REAL*8, DIMENSION(3) :: phi, rotate, rvec, theta, vMPS

jplate = 0
DO 10 i = 1, nPlate
   IF (namTag == names(i)) THEN
		jPlate = i
		GO TO 11
   END IF
10 CONTINUE
11 IF (jPlate == 0) THEN
  write(ErrorMsg,'(A,A2,A/,A)') "ERROR: plate name (",namTag,") not found in lists.", &
&           "Please correct the boundary conditions file."
  call FatalError(ErrorMsg,ThID)
END IF
!     relative Euler vector in radians/second:
DO 20 i = 1, 3
   rotate(i) = (omega(i, jplate) - omega(i, iPVRef)) * 3.168809D-14
20 CONTINUE
!     radius vector to location:
rvec(1) = radius * SIN(xNode(node)) * COS(yNode(node))
rvec(2) = radius * SIN(xNode(node)) * SIN(yNode(node))
rvec(3) = radius * COS(xNode(node))
!     cross product gives velocity in m/s in Cartesian:
vMPS(1) = rotate(2) * rvec(3) - rotate(3) * rvec(2)
vMPS(2) = rotate(3) * rvec(1) - rotate(1) * rvec(3)
vMPS(3) = rotate(1) * rvec(2) - rotate(2) * rvec(1)
vMag = SQRT(vMPS(1)**2 + vMPS(2)**2 + vMPS(3)**2)
!     create unit +Theta and +Phi vectors in Cartesian:
theta(1) = COS(xNode(node)) * COS(yNode(node))
theta(2) = COS(xNode(node)) * SIN(yNode(node))
theta(3) = -SIN(xNode(node))
phi(1) = -SIN(yNode(node))
phi(2) = COS(yNode(node))
phi(3) = 0.0D0
!     find azimuth from dot products:
vTheta = vMPS(1) * theta(1) + vMPS(2) * theta(2) + vMPS(3) * theta(3)
vPhi = vMPS(1) * phi(1) + vMPS(2) * phi(2) + vMPS(3) * phi(3)
vAz = 57.2957795130823D0 * ATan2F(vPhi, -vTheta)
END SUBROUTINE Euler

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

RETURN
END SUBROUTINE FAngls

SUBROUTINE FEM (alpha, area, constr, detJ, &  ! input
&                 dXS, dYS, eta, &
&                 everyP, fBase, fC, fDip, &
&                 fIMuDZ, fLen, fPFlt, fPSfer, fArg, &
&                 fTStar, iCond, iUnitS, iUnitT, &
&                 mxBn, mxDOF, mxEl, mxFEl, mxNode, &
&                 nCond, nDOF, nFl, nLB, nodCon, nodeF, &
&                 nodes, nUB, numEl, numNod, &
&                 oVB, pulled, radius, sita, &
&                 title1, title2, title3, tOfset, trHMax, &
&                 vBCArg, vBCMag, wedge, &
&                 lastpm, &
&                 eRate, v, &                   ! modify
&                 dV, scoreA, scoreB, tauMat, & ! output
&                 f, k, ipiv)                   ! work

!   Computes horizontal velocity of nodes in a thin-shell lithosphere,
!      based on applied forces and boundary conditions.
!   Uses the current strain rate (eRate must be input) as a basis
!      for linearizing the equations by the secant method.

!   Also returns two scores:  scoreA = max_dV, scoreB = RMS_dV/RMS_V

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alpha, area, constr, detJ, dXS, dYS, eta                         ! input
LOGICAL, INTENT(IN) :: everyP                                                          ! input
DOUBLE PRECISION, INTENT(IN) :: fBase                                                  ! input
REAL*8, INTENT(IN) :: fC, fDip, fIMuDZ, fLen, fPFlt, fPSfer, fArg, fTStar              ! input
INTEGER, INTENT(IN) :: iCond, iUnitS, iUnitT, mxBn, mxDOF, mxEl, mxFEl, mxNode         ! input
INTEGER, INTENT(IN) :: nCond, nDOF, nFl, nLB                                           ! input
INTEGER, INTENT(IN) :: nodCon, nodeF, nodes, nUB, numEl, numNod                        ! input
REAL*8, INTENT(IN) :: oVB                                                              ! input
LOGICAL, INTENT(IN) :: pulled                                                          ! input
REAL*8, INTENT(IN) :: radius, sita                                                     ! input
CHARACTER*100, INTENT(IN) :: title1, title2, title3                                     ! input
REAL*8, INTENT(IN) :: tOfset, trHMax, vBCArg, vBCMag, wedge                            ! input
INTEGER, INTENT(IN) :: lastpm                                                          ! input
REAL*8, INTENT(INOUT) :: eRate                                                         ! modify
DOUBLE PRECISION, INTENT(INOUT) :: v                                                   ! modify
REAL*8, INTENT(OUT) :: dV, scoreA, scoreB, tauMat                                      ! output
DOUBLE PRECISION f, k                                                   ! work
INTEGER ipiv                                                            ! work
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
INTEGER i, j
REAL*8 bDenom, bDen, bDenoN, dVSize

DIMENSION alpha(3, 3, 7, mxEl), area(mxEl), detJ(7, mxEl), &
&           dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
&           dV(2, mxNode), &
&           eRate(3, 7, mxEl), eta(7, mxEl), &
&           f(mxDOF, 1), fBase(mxDOF), fC(2, 2, 7, mxFEl), &
&           fDip(2, mxFEl), fIMuDZ(7, mxFEl), fLen(mxFEl), &
&           fPFlt(2, 2, 2, 7, mxFEl), &
&           fPSfer(2, 2, 3, 7, mxEl), &
&           fArg(2, mxFEl), fTStar(2, 7, mxFEl), &
&           iCond(mxBn), ipiv(nRank), &
&           nodCon(mxBn), nodeF(4, mxFEl), nodes(3, mxEl), &
&           oVB(2, 7, mxEl), pulled(7, mxEl), &
&           sita(7, mxEl), tauMat(3, 7, mxEl), &
&           tOfset(3, 7, mxEl), v(2, mxNode), vBCArg(mxBn), &
&           vBCMag(mxBn), k(nKRows, nRank)

IF(lastpm /= 999) THEN
  write(ErrorMsg,'(A)') "WRONG NUMBER OF ARGUMENTS IN CALL TO -FEM-!"
  call FatalError(ErrorMsg,ThID)
END IF
CALL BuildF (area, detJ, dXS, dYS, eta, & ! input
&              fBase, fDip, fLen, fPFlt, &
&              fPSfer, fArg, fTStar, &
&              mxDOF, mxEl, mxFEl, &
&              nDOF, nFl, nodeF, nodes, &
&              numEl, oVB, pulled, &
&              radius, sita, tOfset, trHMax, &
&              wedge, &
&              f)                           ! output
CALL BuildK (alpha, area, detJ, dXS, dYS, & ! input
&              eta, fPSfer, &
&              mxEl, &
&              nodes, numEl, &
&              pulled, radius, sita, trHMax, &
&              k)                             ! output
CALL AddFSt (constr, fC, fDip, fIMuDZ, fLen, fPFlt, fArg, & ! input
&              mxFEl, &
&              nFl, nodeF, &
&              wedge, &
&              k)                                             ! modify
CALL VBCs (iCond, mxBn, mxDOF, &           ! input
&            nCond, nDOF, nLB, nodCon, nUB, &
&            vBCArg, vBCMag, &
&            f, k)                           ! modify

CALL Solver (iUnitT, &         ! input
		  & k, f, &           ! modify (coefficient matrix and forcing vector)
		  & ipiv)             ! work

!  After this CALL, new solution is in "f", and old one in "v".

!  Compare, and compute difference dV and two scores:
bDenom = 0.0D0
bDenoN = 0.0D0
scoreA = 0.0D0
scoreB = 0.0D0
DO 90 i = 1, numNod
	bDenom = bDenom + SQRT(f(2 * i - 1, 1)**2 + f(2 * i, 1)**2)
	bDenoN = bDenoN + SQRT(v(1, i)**2 + v(2, i)**2)
	dV(1, i) = v(1, i) - f(2 * i - 1, 1)
	dV(2, i) = v(2, i) - f(2 * i    , 1)
	dVSize = SQRT(dV(1, i)**2 + dV(2, i)**2)
	scoreA = MAX(scoreA, dVSize)
	scoreB = scoreB + dVSize
90  CONTINUE
bDen = MAX(bDenom, bDenoN)
IF (bDen > 0.0D0) THEN
	scoreB = scoreB / bDen
ELSE
	scoreB = 1.0D0
END IF

!   Transfer new solution to "v", where it will be "old" during next call:

DO 100 i = 1, numNod
	v(1, i) = f(2 * i - 1, 1)
	v(2, i) = f(2 * i    , 1)
100  CONTINUE
IF (everyP) THEN
	WRITE (iUnitS, 10) title1
	WRITE (iUnitS, 10) title2
	WRITE (iUnitS, 10) title3
10       FORMAT (A80)
	WRITE (iUnitS, 20) ((v(j, i), j = 1, 2), i = 1, numNod)
20       FORMAT (1P,4D20.12)
END IF

!   Compute strain rate and stress (the latter according to the
!   current tentative linearization):

CALL EDot (dXS, dYS, &  ! input
&            fPSfer, mxEl, &
&            mxNode, nodes, numEl, radius, sita, v, &
&            eRate)       ! output
CALL TauDef (alpha, eRate, mxEl, numEl, tOfset, & ! input
&              tauMat)                              ! output

RETURN
END SUBROUTINE FEM

SUBROUTINE FillIn (alphaT, basal, conduc, &                  ! input
&                    continuum_LRi, &
&                    cooling_curvature, &
&                    density_anomaly, &
&                    dQdTdA, elev, &
&                    fPSfer, gMean, gradie, &
&                    iConve, iPAfri, iPVRef, iUnitM, iUnitT, &
&                    LRn, LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_eCreep, &
&                    mxEl, mxNode, &
&                    names, nodes, &
&                    nPlate, numEl, numNod, omega, oneKm, &
&                    radio, radius, rhoAst, rhoBar, rhoH2O, &
&                    tAdiab, temLim, tLNode, trHMax, tSurf, &
&                    vTimes, whichP, xNode, yNode, zBAsth, &
&                    zMNode, &
&                    contin, curviness, delta_rho, geothC, geothM, glue, & ! output
&                    oVB, pulled, sigZZI, &
&                    tauZZI, tauZZN, tLInt, vM, zMoho, &
&                    atNode)                                   ! work

!   Precompute and interpolate all "convenience arrays":

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alphaT                                                           ! input
DOUBLE PRECISION, INTENT(IN) :: basal                                                  ! input
REAL*8, INTENT(IN) :: conduc                                                           ! input
INTEGER, INTENT(IN) :: continuum_LRi                                                   ! input
REAL*8, INTENT(IN) :: cooling_curvature, &                                             ! input
				   & density_anomaly, dQdTdA, elev, &                                 ! input
				   & fPSfer, gMean, gradie                                            ! input
INTEGER, INTENT(IN) :: iConve, iPAfri, iPVRef, iUnitM, iUnitT, mxEl, mxNode            ! input
INTEGER, INTENT(IN) :: LRn                                                             ! input
REAL*8, INTENT(IN) :: LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_eCreep       ! input
CHARACTER*2, INTENT(IN) :: names                                                       ! input
INTEGER, INTENT(IN) :: nodes, nPlate, numEl, numNod                                    ! input
REAL*8, INTENT(IN) :: omega, oneKm, radio, radius, rhoAst, rhoBar, rhoH2O, &           ! input
				   & tAdiab, temLim, tLNode, trHMax, tSurf, vTimes                    ! input
INTEGER, INTENT(IN) :: whichP                                                          ! input
REAL*8, INTENT(IN) :: xNode, yNode                                                     ! input
REAL*8, INTENT(IN) :: zBAsth, zMNode                                                   ! input
LOGICAL, INTENT(OUT) :: contin                                                         ! output
REAL*8, INTENT(OUT) :: curviness, delta_rho, geothC, geothM, glue, oVB                 ! output
LOGICAL, INTENT(OUT) :: pulled                                                         ! output
REAL*8, INTENT(OUT) :: sigZZI, tauZZI, tauZZN, tLInt                                   ! output
DOUBLE PRECISION, INTENT(OUT) :: vM                                                    ! output
REAL*8, INTENT(OUT) :: zMoho                                                           ! output
REAL*8 atNode                                                           ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DOUBLE PRECISION points
COMMON / S1S2S3 / points
DIMENSION points(3, 7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, iconv2, m
REAL*8 baseT, delta_quadratic, difMag, dTdZC, dTdZM, &
	& geoth1, geoth2, geoth3, geoth4, geoth5, geoth6, geoth7, geoth8, &
	& huge, q, shrMag, tAsthK, test, terr0r, vtime2, z
DIMENSION alphaT(2), atNode(mxNode), &
&           basal(2, mxNode),  &
&           conduc(2), contin(7, mxEl), &
&           continuum_LRi(mxEl), &
&           cooling_curvature(mxNode), &
&           curviness(7, mxEl), &
&           delta_rho(7, mxEl), &
&           density_anomaly(mxNode), &
&           dQdTdA(mxNode), &
&           elev(mxNode), &
&           fPSfer(2, 2, 3, 7, mxEl), &
&           geothC(4, 7, mxEl), geothM(4, 7, mxEl), &
&           glue(7, mxEl), nodes(3, mxEl), &
&           LR_set_aCreep(1:2, 0:LRn), LR_set_bCreep(1:2, 0:LRn), LR_set_cCreep(1:2, 0:LRn), LR_set_eCreep(0:LRn), &
&           oVB(2, 7, mxEl), &
&           pulled(7, mxEl), radio(2), rhoBar(2), &
&           sigZZI(7, mxEl), tauZZI(7, mxEl), tauZZN(mxNode), &
&           temLim(2), tLNode(mxNode), &
&           tLInt(7, mxEl), whichP(mxNode), &
&           vm(2, mxNode), xNode(mxNode), yNode(mxNode), &
&           zMNode(mxNode), zMoho(7, mxEl)
DIMENSION names(nPlate)
DIMENSION omega(3, nPlate)

DATA huge / 1.0D+30 /

!   Lower-mantle flow at nodes (for computing basal drag)
!   (Notes: If iConve=4, mantle velocity is the same under
!           continents and oceans; however, it is only used
!           for drag computation under continents.
!           If iConve=6, a virtual mantle velocity is created,
!           which differs from velocity model PB2002 (iConve==3~4)
!           by a large differential (or shear) velocity difMag
!           in the direction given by "basal" vectors.)

IF (iConve /= 6) THEN
	CALL Convec (iConve, iPAfri, iPVRef, iUnitM, iUnitT, & ! input
&                   mxNode, &
&                   names, &
&                   nPlate, numNod, &
&                   omega, radius, vTimes, &
&                   whichP, xNode, yNode, &
&                   vM)                                       ! output
ELSE
!           The new case of iConve == 6...
	iconv2 = 3
	vtime2 = 1.0D0
!           Note that use of these parameters will give the PB2002
!           surface velocity model in vM...
	CALL Convec (iConv2, iPAfri, iPVRef, iUnitM, iUnitT, & ! input
&                   mxNode, &
&                   names, &
&                   nPlate, numNod, &
&                   omega, radius, vtime2, &
&                   whichP, xNode, yNode, &
&                   vM)                                       ! output
!           Differential or shear velocity of 100 mm/a:
	difMag = 0.1D0 / 3.15576D7
	DO 610 i = 1, numNod
		  shrMag = SQRT(basal(1, i)**2 + basal(2, i)**2)
		  IF (shrMag > 0.0D0) THEN
			   vM(1, i) = vM(1, i) + (difMag / shrmag) * basal(1, i)
			   vM(2, i) = vM(2, i) + (difMag / shrmag) * basal(2, i)
		  END IF
610       CONTINUE
END IF

!   Same field expressed as values at integration points:

CALL Flow (fPSfer, mxEl, mxNode, nodes, numEl, vm, & ! input
&            oVB)                                      ! output

!   Decide which points are "continental"
!     (a distinction that matters only if iConve=4),
!     using zMoho as temporary storage for interpolated elevation,
!     and tLInt as temporary storage for interpolated heatflow:

CALL Interp (elev, mxEl, mxNode, nodes, numEl, & ! input
&              zMoho)                              ! output
CALL Interp (dQdTdA, mxEl, mxNode, nodes, numEl, & ! input
&              tLInt)                                ! output
DO 2 m = 1, 7
	DO 1 i = 1, numEl
		 contin(m, i) = (zMoho(m, i) > -2500.0D0).AND. &
&                          (tLInt(m, i) < 0.1500D0)
!                Note: Heat-flow limit excludes Iceland, Afar.
1       CONTINUE
2  CONTINUE

!   Thickness of layers:

CALL Interp (zMNode, mxEl, mxNode, nodes, numEl, & ! input
&              zMoho)                                ! output
CALL Interp (tLNode, mxEl, mxNode, nodes, numEl, & ! input
&              tLInt)                                ! output
DO 4 m = 1, 7
	DO 3 i = 1, numEl
		 tLInt(m, i) = MAX(tLInt(m, i), 0.0D0)
3       CONTINUE
4  CONTINUE

!   Density anomaly of chemical origin (applies to whole lithosphere):

CALL Interp (density_anomaly, mxEl, mxNode, nodes, numEl, & ! input
&              delta_rho)                                     ! output

!   Geotherm:

!      -------------- The following method is easy but WRONG!-----------
!CCC   CALL Interp (cooling_curvature, mxEl, mxNode, nodes, numEl, ! input
!CCC +              curviness)                                     ! output
!      ----------The nonlinearities are too great for this approach,----
!                especially when one node of the element is on a
!                spreading ridge.
!                The correct way is to set curviness(m, i) to make the
!                geotherm of each integration point arrive at
!                temperature tAsthK = tAdiab + gradie * 100.D3
!                at depth (in lithosphere) of
!                (zMoho(M,I)+tLInt(M,I)).
!      -----------------------------------------------------------------

tAsthK = tAdiab + gradie * 100.0D3

geoth1 = tSurf
geoth3 = -0.5D0 * radio(1) / conduc(1)
geoth4 = 0.0D0
geoth7 = -0.5D0 * radio(2) / conduc(2)
geoth8 = 0.0D0
DO 90 m = 1, 7
	DO 80 i = 1, numEl

!            N.B. On first pass, omit curviness:

		 geothC(1, m, i) = geoth1
		 q = dQdTdA(nodes(1, i)) * points(1, m) + &
&               dQdTdA(nodes(2, i)) * points(2, m) + &
&               dQdTdA(nodes(3, i)) * points(3, m)
		 geothC(2, m, i) = q / conduc(1)
		 geothC(3, m, i) = geoth3
		 geothC(4, m, i) = geoth4
		 z = zMoho(m, i)
		 geothM(1, m, i) = geothC(1, m, i) + &
&                             geothC(2, m, i) * z + &
&                             geothC(3, m, i) * z**2 + &
&                             geothC(4, m, i) * z**3
		 dTdZC =   geothC(2, m, i) + &
&                2.0D0 * geothC(3, m, i) * z + &
&                3.0D0 * geothC(4, m, i) * z**2
		 dTdZM = dTdZC * conduc(1) / conduc(2)
		 geothM(2, m, i) = dTdZM
		 geothM(3, m, i) = geoth7
		 geothM(4, m, i) = geoth8

!            Now, correct geotherm to hit tAsthK:

		 IF (tLInt(m, i) > 0.0D0) THEN
			  test = geothM(1, m, i) + &
&                       geothM(2, m, i) * tLInt(m, i) + &
&                       geothM(3, m, i) * tLInt(m, i)**2 + &
&                       geothM(4, m, i) * tLInt(m, i)**3
		 ELSE
			  test = geothC(1, m, i) + &
&                       geothC(2, m, i) * zMoho(m, i) + &
&                       geothC(3, m, i) * zMoho(m, i)**2 + &
&                       geothC(4, m, i) * zMoho(m, i)**3
		 END IF
		 terr0r = test - tAsthK
		 delta_quadratic = -terr0r / (zMoho(m, i) + tLInt(m, i))**2
		 curviness(m, i) = -2.0D0 * delta_quadratic
		 geothC(3, m, i) = geoth3 + delta_quadratic
		 geothM(3, m, i) = geoth7 + delta_quadratic
		 geothM(1, m, i) = geothC(1, m, i) + &
&                             geothC(2, m, i) * zMoho(m, i) + &
&                             geothC(3, m, i) * zMoho(m, i)**2 + &
&                             geothC(4, m, i) * zMoho(m, i)**3
		 dTdZC =   geothC(2, m, i) + &
&                2.0D0 * geothC(3, m, i) * zMoho(m, i) + &
&                3.0D0 * geothC(4, m, i) * zMoho(m, i)**2
		 dTdZM = dTdZC * conduc(1) / conduc(2)
		 geothM(2, m, i) = dTdZM
80       CONTINUE
90  CONTINUE

!   Vertical integrals of vertical stress anomaly
!   (relative to a standard pressure curve, in -SQUEEZ-):

DO 100 i = 1, numNod
	geoth2 = dQdTdA(i) / conduc(1)
	geoth3 = -0.50D0 * radio(1) / conduc(1) - 0.50D0 * cooling_curvature(i)
	geoth5 = geoth1 + &
&               geoth2 * zMNode(i) + &
&               geoth3 * zMNode(i)**2 + &
&               geoth4 * zMNode(i)**3
	dTdZC =   geoth2 + &
&           2.0D0 * geoth3 * zMNode(i) + &
&           3.0D0 * geoth4 * zMNode(i)**2
	geoth6 = dTdZC * conduc(1) / conduc(2)
	geoth7 = -0.50D0 * radio(2) / conduc(2) - 0.50D0 * cooling_curvature(i)
	CALL Squeez (alphaT, density_anomaly(i), elev(i), &  ! input
&                   geoth1, geoth2, geoth3, geoth4, &
&                   geoth5, geoth6, geoth7, geoth8, &
&                   gMean, &
&                   iUnitT, oneKm, rhoAst, rhoBar, rhoH2O, &
&                   temLim, zMNode(i), zMNode(i) + tLNode(i), &
&                   tauZZN(i), atNode(i))                   ! output
100  CONTINUE
CALL Interp (atNode, mxEl, mxNode, nodes, numEl, & ! input
&              sigZZI)                               ! output
CALL Interp (tauZZN, mxEl, mxNode, nodes, numEl, & ! input
&              tauZZI)                               ! output

!  Compute strength of shearing layer in asthenosphere:

CALL OneBar (continuum_LRi, &                                                   ! input
&              geothC, geothM, gradie, &                                          ! input
&              LRn, LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_eCreep, & ! input
&              mxEl, numEl, oneKm, tAdiab, &                                      ! input
&              zBAsth, zMoho, &                                                   ! input
&              glue)                                                              ! output


!  However, asthenosphere strength is not considered when
!  iConve = 5; in this case, "glue" is not used; it is set
!  to very large values so that shear tractions will be
!  based on parameters etaMax or trHMax, instead of "glue":

IF (iConve == 5) THEN
	DO 180 m = 1, 7
		 DO 170 i = 1, numEl
			  glue(m, i) = huge
170            CONTINUE
180       CONTINUE
END IF

!   Determine which points have horizontal shear tractions:

DO 200 m = 1, 7
	DO 190 i = 1, numEl
		 IF (iConve <= 3) THEN
			  pulled(m, i) = (trHMax > 0.0D0)
		 ELSE IF (iConve == 4) THEN
			  pulled(m, i) = (trHMax > 0.0D0).AND.contin(m, i)
		 ELSE IF (iConve == 5) THEN
!                     Forearc is defined where base of plate is at less
!                     than 1000 C = 1273 K.
			  baseT = geothM(1, m, i) + &
&                        geothM(2, m, i) * tLInt(m, i) + &
&                        geothM(3, m, i) * tLInt(m, i)**2 + &
&                        geothM(4, m, i) * tLInt(m, i)**3
			  pulled(m, i) = (baseT < 1273.0D0)
		 ELSE IF (iConve == 6) THEN
			  pulled(m, i) = (trHMax > 0.0D0)
!                    (However, even when "pulled" = T for all
!                     integration points, some will still have
!                     zero traction because nodal shear-traction
!                     vectors "basal" are zero for all nodes around
!                     the element.  This happens within plates that
!                     have slab_q = T, where inferred basal traction
!                     is not needed or wanted.)
		 END IF
190       CONTINUE
200  CONTINUE

RETURN
END SUBROUTINE FillIn

SUBROUTINE FindPV (iPAfri, iUnitT, nDPlat, nPBnd, nPlate, & ! input
&                    omega,  pLat,   pLon,   radius, &
&                    xInPl,  xVel,   yInPl,  yVel, &
&                    vPhi,  vTheta)                           ! output

!  Finds out in which plate (xInPl, yInPl) is located,
!  and calculates the velocity of the point from the
!  PB2002 model of Bird [2003; G**3].

!  Requires that "names" and "omega" be pre-filled with names and
!  rotation vectors of all the plates.

!  Requires that -GetPBx- has already been called to fill in the
!  arrays with digitized plate outlines.

!  Returns vPhi (Southward velocity) and vTheta (Eastward velocity)
!  in a reference frame where the AFrica plate is fixed.


IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iPAfri, iUnitT, nDPlat, nPBnd, nPlate                           ! input
REAL*8, INTENT(IN) :: omega, pLat, pLon, radius, xInPl, xVel, yInPl, yVel              ! input
REAL*8, INTENT(OUT) :: vPhi, vTheta                                                    ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, iPlate, j, j2, nEnd, nPoint
REAL*8 a1, a2, a3, aa, ab1, ab2, ab3, ao, angle, b1, b2, b3, bb, bo, &
	& dangle, omegax, omegay, omegaz, oxyz, phi, phix, phiy, phiz, &
	& stheta, tangl, theta, thetax, thetay, thetaz, &
	& vx, vy, vz, xn, xo, xPoint, yn, yo, yPoint, zn, zo
DIMENSION pLat(nPlate, nPBnd), pLon(nPlate, nPBnd)
DIMENSION nDPlat(nPlate), omega(3, nPlate)

xo = COS(yInPl) * SIN(xInPl)
yo = SIN(yInPl) * SIN(xInPl)
zo = COS(xInPl)
oxyz = xo * xo + yo * yo + zo * zo
oxyz = SQRT(oxyz)
xo = xo / oxyz
yo = yo / oxyz
zo = zo / oxyz
nPoint = 0
angle = 0.0D0
iPlate = 0
DO 500 i = 1, nPlate
  tangl = 0.0D0
  nEnd = nDPlat(i)
  DO 300 j = 1, nEnd
	 j2 = j + 1
	 IF(j == nend) THEN
		j2 = 1
	 END IF
	 a1 = COS(pLon(i, j)) * COS(pLat(i, j))
	 a2 = SIN(pLon(i, j)) * COS(pLat(i, j))
	 a3 = SIN(pLat(i, j))
	 b1 = COS(pLon(i, j2)) * COS(pLat(i, j2))
	 b2 = SIN(pLon(i, j2)) * COS(pLat(i, j2))
	 b3 = SIN(pLat(i, j2))
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
	 stheta = (ab1 * xo + ab2 * yo + ab3 * zo) / (aa * bb)
!            prevent stupid abends due to imprecision:
	 stheta = MAX(-1.0D0, MIN(1.0D0, stheta))
	 tangl = tangl + ASIN(stheta)
300     CONTINUE
  dangle = tangl - 3.14159265358979D0
  IF(dangle >= 0.0001D0) THEN
	 nPoint = nPoint + 1
	 iPlate = i
  END IF
500  CONTINUE
IF(nPoint >= 3) THEN
  xPoint = 90.0D0 - xInPl * 57.2957795130823D0
  yPoint = yInPl * 57.2957795130823D0
  WRITE(iUnitT, 505) xPoint, yPoint
505     FORMAT(' POINT ',2F10.3,' WAS FOUND IN MORE THAN TWO PLATES;' &
&          ,' SOMETHING IS WRONG')
END IF
IF(iPlate > 0) THEN
!       Convert to AFrica-fixed, and radians/second:
  omegax = (omega(1, iPlate) - omega(1, iPAfri)) * 3.168809D-14
  omegay = (omega(2, iPlate) - omega(2, iPAfri)) * 3.168809D-14
  omegaz = (omega(3, iPlate) - omega(3, iPAfri)) * 3.168809D-14
!       Convert to length/second:
  omegax = omegax * radius
  omegay = omegay * radius
  omegaz = omegaz * radius
!      Velocity = OMEGA x position:
  theta = xvel
  phi = yvel
  xn = SIN(theta) * COS(phi)
  yn = SIN(theta) * SIN(phi)
  zn = COS(theta)
  vx = omegay * zn - omegaz * yn
  vy = omegaz * xn - omegax * zn
  vz = omegax * yn - omegay * xn
!      Create unit +Theta and +Phi vectors in Cartesian:
  thetax = COS(theta) * COS(phi)
  thetay = COS(theta) * SIN(phi)
  thetaz = -SIN(theta)
  phix = -SIN(phi)
  phiy = COS(phi)
  phiz = 0.0D0
!      Find argument from dot products:
  vTheta = vx * thetax + vy * thetay + vz * thetaz
  vPhi = vx * phix + vy * phiy + vz * phiz
ELSE
  xPoint = 90.0D0 - xInPl * 57.2957795130823D0
  yPoint = yInPl * 57.2957795130823D0
  write(ErrorMsg,'(A,2F13.5,A)') "THE POINT ",xPoint, yPoint," DOES NOT BELONG ANY PLATE !!! Therefore plate velocity is undefined."
  call FatalError(ErrorMsg,ThID)
END IF
RETURN
END SUBROUTINE FindPV

SUBROUTINE Fixed (alphaT, area, conduc, & ! input
&                   density_anomaly, detJ, &
&                   doFB1, doFB2, doFB3, doFB4, &
&                   dQdTdA, dXS, dYS, &
&                   dXSP, dYSP, edgeTS, elev, fDip, fLen, &
&                   fPFlt, fPSfer, fArg, gMean, &
&                   iCond, iUnitT, &
&                   mxBn, mxDOF, mxEl, mxFEl, mxNode, &
&                   nCond, nFl, nodCon, nodeF, nodes, numEl, &
&                   oneKm, radio, radius, rhoAst, &
&                   rhoBar, rhoH2O, sigZZI, sita, &
&                   tauZZI, tauZZN, temLim, tLNode, tSurf, wedge, &
&                   xNode, yNode, zMNode, &
&                   fBase)                  ! output

!   Precompute the fixed part(s) of the forcing vector of the linear
!     systems of equations.
!   LOGICAL switches doFB1, ..., doFB4 allow individual terms to be computed
!     separately (for example, by the multiple CALLs from subprogram -Balanc-).
!   Set all 4 switches to .TRUE. to get the full fBase.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alphaT, area, conduc, density_anomaly, detJ                      ! input
LOGICAL, INTENT(IN) :: doFB1, doFB2, doFB3, doFB4                                      ! input
REAL*8, INTENT(IN) :: dQdTdA, dXS, dYS, dXSP, dYSP                                     ! input
LOGICAL, INTENT(IN) :: edgeTS                                                          ! input
REAL*8, INTENT(IN) :: elev, fDip, fLen, fPFlt, fPSfer, fArg, gMean                     ! input
INTEGER, INTENT(IN) :: iCond, iUnitT, mxBn, mxDOF, mxEl, mxFEl, mxNode, &              ! input
	 & nCond, nFl, nodCon, nodeF, nodes, numEl                                        ! input
REAL*8, INTENT(IN) :: oneKm, radio, radius, rhoAst, rhoBar, rhoH2O, sigZZI, sita, &    ! input
  & tauZZI, tauZZN, temLim, tLNode, tSurf, wedge, xNode, yNode, zMNode                ! input
DOUBLE PRECISION, INTENT(OUT) :: fBase                                                 ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DOUBLE PRECISION fPoint, fPhi, fGauss
DOUBLE PRECISION points, weight
COMMON / SFault / fPoint
COMMON / FPhis /  fPhi
COMMON / FGList / fGauss
COMMON / S1S2S3 / points
COMMON / WgtVec / weight
DIMENSION fPhi(4, 7), fPoint(7), fGauss(7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8 PhiVal, s1, s2, s3, f1, f2, f3 ! statement function and its arguments
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, PARAMETER :: nStep = 100 ! Number of steps to use in vertical integrations:
INTEGER i, j, k, kEle, kEle12, kEle34, krowx, krowy, m, n, n1, n2, n3, n4, nd, node
LOGICAL atSea, ridge
REAL*8 angle, coss, dArea, delta_rho, dip, ds, dx, dy, dz, &
	& elevat, eLong, fAngle, fpp, ft1, ft2, &
	& geoth1, geoth2, geoth3, geoth4, geoth5, geoth6, geoth7, geoth8, &
	& phi, q, s, sideAz, sigzzb, sinA, sinn, slopex, slopey, sMid, &
	& tauzz, theta, tL, toSide, tzz, &
	& x, x0, x1, x2, xout, xta, y, y0, y1, y2, yout, yta, z, zA, zM, zta, zta1, zta2
DOUBLE PRECISION fp1, fp2
DIMENSION alphaT(2), conduc(2), &
&           radio(2),  rhoBar(2), temLim(2)
DIMENSION phi(2), points(3, 7), theta(2), weight(7)
DIMENSION area(mxEl), density_anomaly(mxNode), &
&           detJ(7, mxEl), dQdTdA(mxNode), &
&           dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
&           dXSP(3, 7, mxEl), dYSP(3, 7, mxEl), edgeTS(3, mxEl), &
&           elev(mxNode), fAngle(2), fBase(mxDOF), fDip(2, mxFEl), &
&           fLen(mxFEl), fPFlt(2, 2, 2, 7, mxFEl), fpp(2, 2, 2, 7), &
&           fPSfer(2, 2, 3, 7, mxEl), fArg(2, mxFEl), &
&           iCond(mxBn), nodCon(mxBn), &
&           nodeF(4, mxFEl), nodes(3, mxEl), &
&           sigZZI(7, mxEl), sita(7, mxEl), &
&           tauZZI(7, mxEl), tauZZN(mxNode), tLNode(mxNode), &
&           xNode(mxNode), yNode(mxNode), zMNode(mxNode)

!      Statement function:
PhiVal (s1, s2, s3, f1, f2, f3) = s1 * f1 + s2 * f2 + s3 * f3

!      Initialize accumulator (to be incremented by loops below):
DO 10 i = 1, mxDOF
	fBase(i) = 0.0D0
10  CONTINUE

IF (doFB1.OR.doFB2) THEN
!           IF (doFB1): vertically-integrated topographic stress (tauZZ), and/or
!           IF (doFB2): horizontal components of basal traction anomalies
!               on areas of the triangular continuum elements:
	DO 100 m = 1, 7
		 DO 90 i = 1, numEl
			  dArea = area(i) * detJ(m, i) * weight(m) / radius
			  slopex = 0.0D0
			  slopey = 0.0D0
			  sinA = SIN(sita(m, i))
			  DO 20 j = 1, 3
				   nd = nodes(j, i)
				   zA = zMNode(nd) + tLNode(nd) - elev(nd)
				   slopex = slopex + zA * dXSP(j, m, i)
				   slopey = slopey + zA * dYSP(j, m, i)
!                          Note: These are not dimensionless; divide by "radius"
!                                to get dimensionless slopes.
20                 CONTINUE
			  DO 80 j = 1, 3
				   node = nodes(j, i)
				   krowx = 2 * node - 1
				   krowy = krowx + 1
				   ft1 = -tauZZI(m, i) * (dXS(1, 1, j, m, i) + &
  &                                  dYS(1, 2, j, m, i) / SIN(sita(m, i)) + &
  &                               fPSfer(1, 1, j, m, i) / TAN(sita(m, i)))
				   ft2 = -tauZZI(m, i) * (dXS(2, 1, j, m, i) + &
  &                                  dYS(2, 2, j, m, i) / SIN(sita(m, i)) + &
  &                               fPSfer(2, 1, j, m, i) / TAN(sita(m, i)))
				   fp1 = -sigZZI(m, i) * ( fPSfer(1, 1, j, m, i) * slopex &
  &                                  + fPSfer(1, 2, j, m, i) * slopey)
				   fp2 = -sigZZI(m, i) * ( fPSfer(2, 1, j, m, i) * slopex &
  &                                  + fPSfer(2, 2, j, m, i) * slopey)
				   IF (doFB1) THEN
						fBase(krowx) = fBase(krowx) + dArea * ft1
						fBase(krowy) = fBase(krowy) + dArea * ft2
				   END IF
				   IF (doFB2) THEN
						fBase(krowx) = fBase(krowx) + dArea * fp1
						fBase(krowy) = fBase(krowy) + dArea * fp2
				   END IF
80                 CONTINUE ! j = 1:3
90            CONTINUE ! i = 1:numEl
100       CONTINUE ! m = 1:7
END IF ! (doFB1.OR.doFB2)

IF (doFB3) THEN

!           Effect of anomalous normal traction on exterior
!           side boundaries of triangular continuum elements:
!           since the standard for normal traction anomalies is the
!           pressure under a spreading ridge, boundaries with a ridge
!          (iCond = -1) boundary condition have none.
!           However, boundaries with a "free" (iCond = 0) boundary
!           condition have normal traction equal to vertical,
!           as if the adjacent material had the same crustal thickness
!           and geotherm, but no strength.
!          (Note: These forces will often be overwritten by velocity
!           boundary conditions, but are provided just is case this is not so.)

	DO 200 i = 1, numEl
		 DO 190 j = 1, 3
			  IF (edgeTS(j, i)) THEN
				   n1 = nodes(MOD(j,  3) + 1, i)
				   n2 = nodes(MOD(j + 1, 3) + 1, i)
				   ridge = .FALSE.
				   DO 110 n = 1, nCond
						IF (nodCon(n) == n1) ridge = ridge.OR.(iCond(n) == -1)
						IF (nodCon(n) == n2) ridge = ridge.OR.(iCond(n) == -1)
110                      CONTINUE
				   IF (.NOT.ridge) THEN
						theta(1) = xNode(n1)
						theta(2) = xNode(n2)
						phi(1)  = yNode(n1)
						phi(2)  = yNode(n2)
						eLong = FltLen (phi(1), phi(2), radius, theta(1), theta(2))
						CALL SNodal (phi, theta, & ! input
&                                       fpp)          ! output
						CALL FAngls (phi, theta, & ! input
&                                       fAngle)       ! output
						DO 180 m = 1, 7
							 s = fPoint(m)
							 tzz = tauZZN(n1) * fPhi(1, m) + tauZZN(n2) * fPhi(2, m)
							 ds = fGauss(m) * eLong
							 coss = COS(fAngle(1)) * fPhi(1, m) &
&                                    + COS(fAngle(2)) * fPhi(2, m)
							 sinn = SIN(fAngle(1)) * fPhi(1, m) &
&                                    + SIN(fAngle(2)) * fPhi(2, m)
							 angle = ATan2F(sinn, coss)
							 xout = COS(angle - 1.57079632679490D0)
							 yout = SIN(angle - 1.57079632679490D0)
							 krowx = 2 * n1 - 1
							 krowy = krowx + 1
							 fBase(krowx) = fBase(krowx) + &
&                                              ds * (fpp(1, 1, 1, m) * xout + &
&                                              fpp(1, 2, 1, m) * yout) * tzz
							 fBase(krowy) = fBase(krowy) + &
&                                         ds * (fpp(2, 1, 1, m) * xout + &
&                                         fpp(2, 2, 1, m) * yout) * tzz
							 krowx = 2 * n2 - 1
							 krowy = krowx + 1
							 fBase(krowx) = fBase(krowx) + &
&                                              ds * (fpp(1, 1, 2, m) * xout + &
&                                              fpp(1, 2, 2, m) * yout) * tzz
							 fBase(krowy) = fBase(krowy) + &
&                                              ds * (fpp(2, 1, 2, m) * xout + &
&                                              fpp(2, 2, 2, m) * yout) * tzz
180                           CONTINUE ! m = 1:7
				   END IF ! .NOT.ridge
			  END IF ! edgeTS(j, i)
190            CONTINUE ! j = 1:3
200       CONTINUE ! i = 1:numEl
END IF ! doFB3

IF (doFB4) THEN

!           Effect of vertical-stress (sigZZ) component of normal traction on
!           fault planes is obtained by integrating down dip of each
!           fault at each of the seven integration points along its length:

	DO 300 i = 1, nFl
		 n1 = nodeF(1, i)
		 n2 = nodeF(2, i)
		 n3 = nodeF(3, i)
		 n4 = nodeF(4, i)
		 x1 = xNode(n1)
		 x2 = xNode(n2)
		 y1 = yNode(n1)
		 y2 = yNode(n2)

!                Find neighboring triangular elements, if any:

		 kEle12 = 0
		 DO 210 j = 1, numEl
			  IF (((nodes(1, j) == n2).AND. &
&                     (nodes(2, j) == n1)     ) .OR. &
&                    ((nodes(3, j) == n2).AND. &
&                     (nodes(1, j) == n1)     ) .OR. &
&                    ((nodes(2, j) == n2).AND. &
&                     (nodes(3, j) == n1)     )     ) THEN
				   kEle12 = j
				   GO TO 211
			  END IF
210            CONTINUE
211            kEle34 = 0
		 DO 220 j = 1, numEl
			  IF (((nodes(1, j) == n3).AND. &
&                     (nodes(3, j) == n4)     ) .OR. &
&                    ((nodes(3, j) == n3).AND. &
&                     (nodes(2, j) == n4)     ) .OR. &
&                    ((nodes(2, j) == n3).AND. &
&                     (nodes(1, j) == n4)     )     ) THEN
				   kEle34 = j
				   GO TO 221
			  END IF
220            CONTINUE
221            DO 290 m = 1, 7
			  s = fPoint(m)
			  x0 = x1 * fPhi(1, m) + x2 * fPhi(2, m)
			  y0 = y1 * fPhi(1, m) + y2 * fPhi(2, m)

!CCCC                 angle = fArg(1, i) * fPhi(1, m) + fArg(2, i) * fPhi(2, m)
!CCCC                 Line above was replaced due to cycle-shift problem!

			  angle = Chord(fArg(1, i), fPhi(2, m), fArg(2, i))

			  xout = COS(angle - 1.57079632679490D0)
			  yout = SIN(angle - 1.57079632679490D0)
			  dip = fDip(1, i) * fPhi(1, m) + fDip(2, i) * fPhi(2, m)
			  IF (ABS(dip - 1.57079632679490D0) < wedge) THEN
!                          Case of vertical dip (within "wedge" radians):
				   tzz = tauZZN(n1) * fPhi(1, m) + tauZZN(n2) * fPhi(2, m)
			  ELSE
!                          Case of shallow dip:
				   IF (dip > 1.57079632679490D0) THEN
!                               dip is toward n3,n4 side:
						kEle = kEle34
				   ELSE
!                               dip is toward n1,n2 side:
						kEle = kEle12
				   END IF
				   IF (kEle == 0) THEN
!                               no neighboring element (at grid edge):
						tzz = tauZZN(n1) * fPhi(1, m) + &
&                                tauZZN(n2) * fPhi(2, m)
				   ELSE
!                               Integrate on a slant below neighbor elements.
!                               (1) Find intersection of fault with
!                                   asthenosphere (plate base):
						zta1 = zMNode(n1) + tLNode(n1)
						zta2 = zMNode(n2) + tLNode(n2)
						zta = zta1 * fPhi(1, m) + zta2 * fPhi(2, m)
						toSide = zta / TAN(dip)
!                               Note: toSide is <0 for DIP > Pi/2.
						sideAz = angle - 1.57079632679490D0
						xta = x0 + toSide * COS(sideAz) / radius
						yta = y0 + toSide * SIN(sideAz) / (radius * SIN(x0))
!                               (2) Subdivide slant path into steps
						dx = (xta - x0) / nStep
						dy = (yta - y0) / nStep
						dz = zta / nStep
!                               (3) Actual integration on slant path:
						s1 = 0.3333D0
						s2 = 0.3333D0
						s3 = 0.3334D0
						tzz = 0.0D0
						DO 250 k = 1, nStep
							 sMid = k - 0.5D0
							 x = x0 + smid * dx
							 y = y0 + smid * dy
							 z = sMid * dz
							 CALL Lookup (iUnitT, mxEl, mxFEl, & ! input
&                                            mxNode, nFl, nodeF, &
&                                            nodes, numEl, &
&                                            x, xNode, y, yNode, &
&                                            kEle, s1, s2, s3, &    ! modify
&                                            atSea)                 ! output
							 IF (atSea) THEN
								  tzz = tauZZN(n1) * fPhi(1, m) + &
&                                          tauZZN(n2) * fPhi(2, m)
								  GO TO 251
							 ELSE
								  elevat = PhiVal(s1, s2, s3, &
&                                                    elev(nodes(1, kEle)), &
&                                                    elev(nodes(2, kEle)), &
&                                                    elev(nodes(3, kEle)))
								  delta_rho = PhiVal(s1, s2, s3, &
&                                                       density_anomaly(nodes(1, kEle)), &
&                                                       density_anomaly(nodes(2, kEle)), &
&                                                       density_anomaly(nodes(3, kEle)))
								  q = PhiVal(s1, s2, s3, &
&                                               dQdTdA(nodes(1, kEle)), &
&                                               dQdTdA(nodes(2, kEle)), &
&                                               dQdTdA(nodes(3, kEle)))
								  zM = PhiVal(s1, s2, s3, &
&                                                zMNode(nodes(1, kEle)), &
&                                                zMNode(nodes(2, kEle)), &
&                                                zMNode(nodes(3, kEle)))
								  tL = PhiVal(s1, s2, s3, &
&                                                tLNode(nodes(1, kEle)), &
&                                                tLNode(nodes(2, kEle)), &
&                                                tLNode(nodes(3, kEle)))
!                                         (4) Terminate integral if it
!                                             emerges into asthenosphere
!                                             anywhere along slant path:
								  IF (z > (zM + tL)) GO TO 251
								  geoth1 = tSurf
								  geoth2 = q / conduc(1)
								  geoth3 = -0.5D0 * radio(1) / conduc(1)
								  geoth4 = 0.0D0
								  geoth5 = geoth1 + geoth2 * zm + &
&                                             geoth3 * zm**2
								  geoth6 = (q - zm * radio(1)) / conduc(2)
								  geoth7 = -0.5D0 * radio(2) / conduc(2)
								  geoth8 = 0.0D0
								  CALL Squeez (alphaT, &      ! input
&                                                 delta_rho, &
&                                                 elevat, &
&                                                 geoth1, geoth2, &
&                                                 geoth3, geoth4, &
&                                                 geoth5, geoth6, &
&                                                 geoth7, geoth8, &
&                                                 gMean, iUnitT, &
&                                                 oneKm, rhoAst, &
&                                                 rhoBar, rhoH2O, &
&                                                 temLim, zm, z, &
&                                                 tauzz, sigzzb) ! output
								  tzz = tzz + sigzzb * dz
							 END IF ! atSea, or NOT
250                           CONTINUE ! k = 1:nStep
251                           CONTINUE ! (possible exit from loop above)
				   END IF ! kEle == 0, or NOT
			  END IF ! vertical, or dipping fault
!                     - - - - - - - - - - - -  - - - -
			  dS = fGauss(m) * fLen(i) * tzz ! NOTE inclusion of tzz in "dS" here.
!                     - - - - - - - - - - - -  - - - -
			  krowx = 2 * n1 - 1
			  krowy = krowx + 1
			  fBase(krowx) = fBase(krowx) - &
&                               dS * (fPFlt(1, 1, 1, m, i) * xout + &
&                               fPFlt(1, 2, 1, m, i) * yout)
			  fBase(krowy) = fBase(krowy) - &
&                               dS * (fPFlt(2, 1, 1, m, i) * xout + &
&                               fPFlt(2, 2, 1, m, i) * yout)
			  krowx = 2 * n2 - 1
			  krowy = krowx + 1
			  fBase(krowx) = fBase(krowx) - &
&                               dS * (fPFlt(1, 1, 2, m, i) * xout + &
&                               fPFlt(1, 2, 2, m, i) * yout)
			  fBase(krowy) = fBase(krowy) - &
&                               dS * (fPFlt(2, 1, 2, m, i) * xout + &
&                               fPFlt(2, 2, 2, m, i) * yout)
			  krowx = 2 * n3 - 1
			  krowy = krowx + 1
			  fBase(krowx) = fBase(krowx) + &
&                               dS * (fPFlt(1, 1, 2, m, i) * xout + &
&                               fPFlt(1, 2, 2, m, i) * yout)
			  fBase(krowy) = fBase(krowy) + &
&                               dS * (fPFlt(2, 1, 2, m, i) * xout + &
&                               fPFlt(2, 2, 2, m, i) * yout)
			  krowx = 2 * n4 - 1
			  krowy = krowx + 1
			  fBase(krowx) = fBase(krowx) + &
&                               dS * (fPFlt(1, 1, 1, m, i) * xout + &
&                               fPFlt(1, 2, 1, m, i) * yout)
			  fBase(krowy) = fBase(krowy) + &
&                               dS * (fPFlt(2, 1, 1, m, i) * xout + &
&                               fPFlt(2, 2, 1, m, i) * yout)
290            CONTINUE ! ending loop on integration points m = 1:7
300       CONTINUE ! ending loop on fault elements i = 1:nFl
END IF ! doFB4
RETURN
END SUBROUTINE Fixed

SUBROUTINE Flow (fPSfer, mxEl, mxNode, nodes, numEl, v, & ! input
&                  outVec)                                  ! output

!   Calculates velocity vectors at integration points, from nodal values

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: FPSfer                                                           ! input
INTEGER, INTENT(IN) :: mxEl, mxNode, nodes, numEl                                      ! input
DOUBLE PRECISION, INTENT(IN) :: v                                                      ! input
REAL*8, INTENT(OUT) :: outVec                                                          ! ouput
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, j, m, nji
DIMENSION fPSfer(2, 2, 3, 7, mxEl), nodes(3, mxEl), outVec(2, 7, mxEl), &
&           v(2, mxNode)
DO 50 m = 1, 7
	DO 40 i = 1, numEl
		 outVec(1, m, i) = 0.0D0
		 outVec(2, m, i) = 0.0D0
40       CONTINUE
50  CONTINUE
DO 100 j = 1, 3
	DO 90 m = 1, 7
		 DO 80 i = 1, numEl
			  nji = nodes(j, i)
			  outVec(1, m, i) = outVec(1, m, i) &
&                                + v(1, nji) * fPSfer(1, 1, j, m, i) &
&                                + v(2, nji) * fPSfer(2, 1, j, m, i)
			  outVec(2, m, i) = outVec(2, m, i) &
&                                + v(1, nji) * fPSfer(1, 2, j, m, i) &
&                                + v(2, nji) * fPSfer(2, 2, j, m, i)
80            CONTINUE
90       CONTINUE
100  CONTINUE
RETURN
END SUBROUTINE Flow

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
RETURN
END FUNCTION FltLen

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
IF(Verbose) WRITE (iUnitVerb, 3) title1
3  FORMAT(/' Title of finite-element grid ='/' ',A80)

!   Read number of nodes, plus out-dated parameters that once
!     permitted boundary nodes to be specially numbered as
!    "fake" nodes with numbers from n1000+1 ... n1000+nFakeN.
!     This option is no longer supported by my programs!
!    (Option "brief" suppresses most output.)

READ (iUnit7, * ) numNod, nRealN, nFakeN, n1000, brief

IF (numNod /= (nRealN + nFakeN)) THEN
  write(ErrorMsg,'(A,I0,A,I0,A,I0)') "numNod (",numNod,") IS NOT EQUAL TO SUM OF nRealN (",nRealN,") AND nFakeN (",nFakeN,")."
  call FatalError(ErrorMsg,ThID)
END IF

IF (nRealN > n1000) THEN
  write(ErrorMsg,'(A,I0,A,I0)') "nRealN (",nRealN,") IS GREATER THAN n1000 (",n1000,")."
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
	IF(Verbose) WRITE (iUnitVerb, 35)
35       FORMAT(/' (Since option ""brief"" = .TRUE., grid will not be echoed here.)')
ELSE
	IF(Verbose) WRITE (iUnitVerb, 40) numNod
40       FORMAT (/' There are',I5,' nodes in the grid')
	IF(Verbose) WRITE (iUnitVerb, 50)
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
		   write(ErrorMsg,'(A,I0)') "ILLEGAL NODE NUMBER: ",index
		   call FatalError(ErrorMsg,ThID)
		 END IF
	END IF
	pLon = vector(2)
	pLat = vector(3)
	IF (ABS(pLat) > 90.01) THEN
	  write(ErrorMsg,'(A,I0)') "ABS(latitude) > 90 AT NODE: ",index
	  call FatalError(ErrorMsg,ThID)
	END IF
	IF (ABS(pLat) > 89.99D0) THEN
	  write(ErrorMsg,'(A,I0,A/,A)') "NODE ",index," LIES ON A POLE. THIS IS A SINGULAR POINT OF THE ",&
	    &  "SPHERICAL COORDINATE SYSTEM. MOVE THIS NODE, AT LEAST SLIGHTLY."
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
		 IF(Verbose) WRITE (iUnitVerb, 99) INDEX, pLon, pLat, xi, yi, elevi, &
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
  allocate(ErrorArray(count(.not.checkN)))
  j = 1
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

READ (iUnit7, * ) numEl
IF (numEl > mxEl) THEN
  write(ErrorMsg,'(A,I0,A)') "INCREASE PARAMETER maxEl TO BE AT LEAST EQUALTO THE NUMBER OF ELEMENTS (",numEl,") AND RECOMPILE."
  call FatalError(ErrorMsg,ThID)
END IF
DO 109 k = 1, numEl
	checkE(k) = .FALSE.
109  CONTINUE
IF (.NOT.brief) THEN
	IF(Verbose) WRITE (iUnitVerb, 110) numEl
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
	  write(ErrorMsg,'(A,I0)') "ILLEGAL ELEMENT NUMBER: ",i
	  call FatalError(ErrorMsg,ThID)
	END IF
	checkE(i) = .TRUE.
	IF (.NOT.brief) THEN
		 IF (LRi == 0) THEN
			 IF(Verbose) WRITE (iUnitVerb, 120) i, (nodes(j, i), j = 1, 3)
120                FORMAT (' ', I6, ':', 3I10)
		 ELSE
			 IF(Verbose) WRITE (iUnitVerb, 121) i, (nodes(j, i), j = 1, 3), LRi
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
  write(ErrorMsg,'(A)') "THE FOLLOWING ELEMENTS WERE NEVER READ:"
  allocate(ErrorArray(size(checkN)-count(checkN)))
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
IF (.NOT.brief) WRITE(iUnitT, 230) nFl
230  FORMAT(/ /' There are ', I6, ' great-circle fault elements.')
IF ((.NOT.brief).AND.(nFl > 0)) WRITE(iUnitT, 231)
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
	  write(ErrorMsg,'(A,I0)') "ILLEGAL FAULT NUMBER: ", i
	  call FatalError(ErrorMsg,ThID)
	END IF
	checkF(i) = .TRUE.
	IF (.NOT.brief) THEN
		 IF (LRi == 0) THEN
			 IF(Verbose) WRITE (iUnitVerb, 240) i, (nodeF(j, i), j = 1, 4), (dips(l), l = 1, 2), off
		 ELSE
			 IF(Verbose) WRITE (iUnitVerb, 242) i, (nodeF(j, i), j = 1, 4), (dips(l), l = 1, 2), off, LRi
242                FORMAT (' ', I6, ':', 4I5, 1X, 2F6.1, 1X, F9.0, " LR", I8)
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
		 IF(Verbose) WRITE (iUnitVerb, 400) offMax
400            FORMAT (/' Greatest fault offset read was ',1P,D10.2)
	ELSE
		 IF(Verbose) WRITE (iUnitVerb, 401)
401            FORMAT (/' Since fault offsets are all zero,', &
&               ' input parameter Byerly will have no effect.')
	END IF
END IF
IF (Verbose) WRITE (iUnitVerb, 999)
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
INTEGER, INTENT(IN) :: iUnitM, iUnitT                                                  ! input
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
120 CONTINUE
	write(ErrorMsg,'(A,I3)') "ERR0R: BAD PLATE NAME ON INPUT DEVICE ",iUnitM
	call FatalError(ErrorMsg,ThID)

140 nRead = nRead + 1
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
201 IF(nRead < nPlate) THEN
	    write(ErrorMsg,'(A,I3,A,I3)') "ERROR: Expecting ",nPlate," plates but read outlines of only ",nRead
		call FatalError(ErrorMsg,ThID)
    END IF

RETURN
END SUBROUTINE GetPBx

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

SUBROUTINE KSize (brief, iUnitP, iUnitT, mxEl, mxFEl, mxNode, & ! input
&                   nFl, nodeF, nodes, numEl, numNod, &
&                   nDOF, nLB, nUB, &                             ! output (+ more in un-named COMMON)
&                   jCol1, jCol2)                                 ! work

!   Characterize the size and shape of the banded linear system,
!   and compute values for the INTEGER variables in un-named COMMON.
!   Determine the lower and upper half-bandwidths of the stiffness
!   matrix by proceeding through the same loops as will be used to
!   create it.
!   The calculation is done in terms of node numbers first, and then
!   the results are (almost) doubled to account for two degrees of
!   freedom per node.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
LOGICAL, INTENT(IN) :: brief                                                           ! input
INTEGER, INTENT(IN) :: iUnitP, iUnitT, mxEl, mxFEl, mxNode, &                          ! input
					& nFl, nodeF, nodes, numEl, numNod                                ! input
INTEGER, INTENT(OUT) :: nDOF, nLB, nUB                                                 ! output (+ more in un-named COMMON)
INTEGER jCol1, jCol2                                                    ! work
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
CHARACTER*1 blank, star, asc1, asc2, ascr
INTEGER i, j, k, nc, nr
LOGICAL worst1, worst2, worstr
REAL*8 rKSize, rKGB ! to avoid any risk of integer-overflow
DIMENSION jCol1(mxNode), jCol2(mxNode), &
&           nodeF(4, mxFEl), nodes(3, mxEl)
DATA blank / ' ' / , star / '*' /

!  Initialize bandwidth to 1 node:
DO 10 i = 1, numNod
	jCol1(i) = i
	jCol2(i) = i
10  CONTINUE
!  Band widening by triangular continuum elements:
DO 50 i = 1, numEl
	DO 40 j = 1, 3
	   nr = nodes(j, i)
	   DO 30 k = 1, 3
		  nc = nodes(k, i)
		  jCol1(nr) = MIN(jCol1(nr), nc)
		  jCol2(nr) = MAX(jCol2(nr), nc)
30          CONTINUE
40       CONTINUE
50  CONTINUE
!  Band widening by linear fault elements:
DO 80 i = 1, nFl
	DO 70 j = 1, 4
		 nr = nodeF(j, i)
		 DO 60 k = 1, 4
			  nc = nodeF(k, i)
			  jCol1(nr) = MIN(jCol1(nr), nc)
			  jCol2(nr) = MAX(jCol2(nr), nc)
60            CONTINUE
70       CONTINUE
80  CONTINUE

nLB = 0
nUB = 0
DO 190 i = 1, numNod
	nLB = MAX(nLB, i - jCol1(i))
	nUB = MAX(nUB, jCol2(i) - i)
190  CONTINUE
IF (.NOT.brief) THEN
	WRITE(iUnitT, 200)
200       FORMAT(/ /' Table of most distant connections between', &
&      ' nodes'/ &
&      ' (* marks the cases which determine the bandwidth)'/ / &
&      ' Lowest-connection   Node   Highest-connection')
	DO 220 i = 1, numNod
		 worst1 = (i - jCol1(i)) == nLB
		 worst2 = (jCol2(i) - i) == nUB
		 worstr = worst1.OR.worst2
		 asc1 = blank
		 asc2 = blank
		 ascr = blank
		 IF (worst1) asc1 = star
		 IF (worst2) asc2 = star
		 IF (worstr) ascr = star
		 IF(Verbose) WRITE (iUnitVerb, 210) jCol1(i), asc1, i, ascr, jCol2(i), asc2
210            FORMAT(' ',I12,A1,I11,A1,I11,A1)
220       CONTINUE
END IF

!  Correct numbers for presence of two degrees of freedom per node:

nDOF = 2 * numNod
nRank = nDOF ! stored in un-named COMMON
nLB = 2 * nLB + 1
nUB = 2 * nUB + 1
nCodiagonals = MAX(nLB, nUB) ! stored in un-named COMMON
nKRows = 3 * nCodiagonals + 1 ! stored in un-named COMMON
iDiagonal = 2 * nCodiagonals + 1 ! stored in un-named COMMON

! The manual page for "?gbsv" of MKL/LAPACK is here: https://software.intel.com/en-us/node/468882

! Be CAREFUL because some pages in the MKL/LAPACK manual give INCORRECT descriptions (and
!    illustrations, and examples!) of the band-storage scheme, describing only kl+1+ku rows in ab,
!    and the diagonal row located at row kl+1.

! Correct documentation can be found here: http://www.netlib.no/netlib/lapack/double/dgbsv.f
! The CORRECT storage scheme has 2*kl+1+ku rows, with the diagonal row located at row 2*kl+1.
! I have checked (2016.07.08) that solutions under this correct scheme are essentially
! identical to old solutions, using an IMSL solver, in Shells_v3.9.

!  Compute and report size of stiffness matrix (in MKL-style banded storage):
rKSize = (1.0D0 * nKRows) * (1.0D0 * nRank) ! switching to REAL*8 to avoid risk of integer-overflow
rKGB = 8.0D0 * rKSize / 1073741824.0D0
IF(Verbose) WRITE (iUnitVerb, *)
IF(Verbose) WRITE (iUnitVerb, 500) nRank, nCodiagonals, nKRows, rKSize, rKGB
WRITE (iUnitP, *)
WRITE (iUnitP, 500) nRank, nCodiagonals, nKRows, rKSize, rKGB
500    FORMAT (' Size of banded linear system:'/ &
	 & ' nRank = ',I8,', nCodiagonals = ',I8,', nKRows = ',I8,', so'/ &
	 & ' number of REAL*8 values = ',F12.0 / &
	 & ' Stiffness-matrix storage requires ',F10.3,' GB of memory.')

RETURN
END SUBROUTINE KSize

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
    write(ErrorMsg,'(A,1P,D12.4,A/,A,D12.4,A/,A)') "AREA OF GRID (",totalA,") EXCEEDS " , &
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
IF(Verbose) WRITE(iUnitVerb, 50) totalA, totalV, thick, side, constr, etaMax, &
&               fMuMax, visMax
IF(Verbose) WRITE (iUnitVerb, 50) totalA, totalV, thick, side, constr, etaMax, &
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
RETURN
END SUBROUTINE Limits

SUBROUTINE Lookup (iUnitT, mxEl, mxFEl, mxNode, & ! input
&                    nFl, nodeF, nodes, numEl, &
&                    x, xNode, y, yNode, &
&                    iE, s1, s2, s3, &              ! modify
&                    atSea)                         ! output

!   Finds element and internal coordinates in grid matching location
!   of a particular point (x, y) and reports them as iE and (s1, s2, s3).

!   Note that x is colatitude (from North pole) and y is
!   East longitude.  Both are in radians.

!   A returned value of atSea indicates that point fell off edge
!   of the grid.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iUnitT, mxEl, mxFEl, mxNode, nFl, nodeF, nodes, numEl           ! input
REAL*8, INTENT(IN) :: x, xNode, y, yNode                                               ! input
INTEGER, INTENT(INOUT) :: iE                                                           ! modify
REAL*8, INTENT(INOUT) :: s1, s2, s3                                                    ! modify
LOGICAL, INTENT(OUT) :: atSea                                                          ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8 PhiVal, f1, f2, f3
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, PARAMETER :: nToTry = 50
INTEGER i1, i2, i3, iEHist, j, k, kEle, kFault, limit, n, nRefin, nTried
LOGICAL trubbl
REAL*8 cf11, cf12, cf13, cf21, cf22, cf23, &
	& coef11, coef12, coef13, coef21, coef22, coef23, det, detI, &
	& delx, dely, ds1, ds2, ds3, dStep, err, &
	& m11, m12, m13, m21, m22, m23, &
	& sHist, step11, step12, step21, step22, step31, step32, &
	& x1, x2, x3, xt, y1, y2, y3, yt
DIMENSION nodeF(4, mxFEl), nodes(3, mxEl), &
&           xNode(mxNode), yNode(mxNode)
DIMENSION iEHist(nToTry), sHist(3, nToTry)

!   Statement function:
PhiVal(s1, s2, s3, f1, f2, f3) = s1 * f1 + s2 * f2 + s3 * f3


nTried = 0

!      Loop as many times as needed:

100  nTried = nTried + 1
	iEHist(nTried) = iE
	IF (nTried >= (nToTry - 10)) THEN
		 trubbl = (iEHist(nTried) == iEHist(nTried - 2))
	ELSE
		 trubbl = .FALSE.
	END IF
	IF (trubbl) THEN
		 atSea = .TRUE.
		 RETURN
	END IF
	i1 = nodes(1, ie)
	i2 = nodes(2, ie)
	i3 = nodes(3, ie)
	x1 = xNode(i1)
	x2 = xNode(i2)
	x3 = xNode(i3)
	y1 = yNode(i1)
	y2 = yNode(i2)
	y3 = yNode(i3)
	s3 = 1.0D0 - s1 - s2
	limit = 3
	nRefin = 0

!           Loop to refine estimate of internal coordinates:

150       nRefin = nRefin + 1
		 xt = PhiVal(s1, s2, s3, x1, x2, x3)
		 yt = PhiVal(s1, s2, s3, y1, y2, y3)

!                COEF:=MAT((DXDS1,DXDS2,DXDS3),
!                          (DYDS1,DYDS2,DYDS3),(1,1,1));

		 coef11 = x1
		 coef12 = x2
		 coef13 = x3
		 coef21 = y1
		 coef22 = y2
		 coef23 = y3
		 m11 = coef22 - coef23
		 m12 = coef21 - coef23
		 m13 = coef21 - coef22
		 m21 = coef12 - coef13
		 m22 = coef11 - coef13
		 m23 = coef11 - coef12
		 cf11 = + m11
		 cf12 = -m12
		 cf13 = + m13
		 cf21 = -m21
		 cf22 = + m22
		 cf23 = -m23
		 det = coef11 * cf11 + coef12 * cf12 + coef13 * cf13
		 IF (det == 0.00D0) THEN
			  atSea = .TRUE.
			  IF(Verbose) WRITE (iUnitVerb, 151)
151                 FORMAT (' LOOKUP IS atSea.')
			  RETURN
		 END IF
		 detI = 1.0D0 / det
		 step11 = cf11
		 step12 = cf21
		 step21 = cf12
		 step22 = cf22
		 step31 = cf13
		 step32 = cf23
		 delx = x - xt
		 dely = y - yt
		 ds1 = (step11 * delx + step12 * dely) * detI
		 ds2 = (step21 * delx + step22 * dely) * detI
		 ds3 = (step31 * delx + step32 * dely) * detI
		 err = (ds1 + ds2 + ds3) / 3.0D0
		 ds1 = ds1 - err
		 ds2 = ds2 - err
		 ds3 = ds3 - err
		 dStep = MAX(ABS(ds1), ABS(ds2), ABS(ds3))
		 IF (dStep >  0.100D0) THEN
			  limit = limit + 1
			  ds1 = ds1 * 0.10D0 / dStep
			  ds2 = ds2 * 0.10D0 / dStep
			  ds3 = ds3 * 0.10D0 / dStep
		 END IF
		 s1 = s1 + ds1
		 s2 = s2 + ds2
		 s3 = s3 + ds3

!  Loop-back (with some conditions):

	IF (((nRefin < limit).AND.(limit <= (nToTry - 10))).AND. &
&          ((s1 >= -0.10D0).AND.(s1 <= 1.10D0)).AND. &
&          ((s2 >= -0.10D0).AND.(s2 <= 1.10D0)).AND. &
&          ((s3 >= -0.10D0).AND.(s3 <= 1.10D0))) GO TO 150

!  Point is now as well-located as possible "in" the current element;
!  however, the internal coordinates may not all be positive, so
!  point may be outside, and we may need to shift to a new element.

	sHist(1, nTried) = s1
	sHist(2, nTried) = s2
	sHist(3, nTried) = s3
	IF (trubbl.OR.(nTried >= nToTry)) THEN
		 WRITE(iUnitT, 201) x, y
201            FORMAT(' REQUEST FOR VALUE AT LOCATION', &
&             ' (',1P,E10.2,',',E10.2,') CAUSES ', &
&             'INFINITE LOOP IN -Lookup-.'/ &
&            ' HISTORY OF SEARCH: element     s1          s2', &
&            '          s3')
		 DO 203 n = 1, nTried - 1
			  WRITE(iUnitT, 202) iEHist(n), (sHist(k, n), k = 1, 3)
202                 FORMAT(22X,I3,2X,3F12.4)
203            CONTINUE
		 WRITE(iUnitT, 204) iEHist(nTried - 1), &
&                             (nodes(j, iEHist(nTried - 1)), j = 1, 3), &
&                             (xNode(nodes(j, iEHist(nTried - 1))), j = 1, 3), &
&                             (yNode(nodes(j, iEHist(nTried - 1))), j = 1, 3)
		 WRITE(iUnitT, 204) iEHist(nTried), &
&                             (nodes(j, iEHist(nTried)), j = 1, 3), &
&                             (xNode(nodes(j, iEHist(nTried))), j = 1, 3), &
&                             (yNode(nodes(j, iEHist(nTried))), j = 1, 3)
204            FORMAT(' Element',I3,' Nodes:',I3,2I10/ &
&                  9X,'X:',1P,3E10.2/9X,'Y:',3E10.2)
		 RETURN
	END IF
	IF (s1 > -0.03D0) THEN
		 IF (s2 > -0.03D0) THEN
			  IF (s3 > -0.03D0) THEN
!                          Point has been successfully found!
				   atSea = .FALSE.
				   RETURN
			  ELSE
				   CALL Next (iE, 3, mxEl, mxFEl, nFl, & ! input
&                                nodeF, nodes, numEl, &
&                                kFault, kEle)              ! output
			  END IF
		 ELSE
			  CALL Next (iE, 2, mxEl, mxFEl, nFl, & ! input
&                           nodeF, nodes, numEl, &
&                           kFault, kEle)              ! output
		 END IF
	ELSE
		 CALL Next (iE, 1, mxEl, mxFEl, nFl, & ! input
&                      nodeF, nodes, numEl, &
&                      kFault, kEle)              ! output
	END IF
	IF (kEle > 0) THEN
		 iE = kEle
		 s1 = 0.3333D0
		 s2 = 0.3333D0
		 s3 = 0.3334D0
		 GO TO 100
	ELSE
		 atSea = .TRUE.
		 RETURN
	END IF

!  Note: Indentation reflects indefinite loop on trial element iE.

END SUBROUTINE Lookup

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

SUBROUTINE OldVel (iUnitT, iUnitV, mxNode, numNod, & ! input
&                    v)                                ! output

!   Read old velocity solution "v" from unit iUnitV, or else fills this array
!     with zeros.  Comments are sent to unit iUnitT.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iUnitT, iUnitV, mxNode, numNod                                  ! input
DOUBLE PRECISION, INTENT(OUT) :: v                                                     ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CHARACTER*100 title1, title2, title3
INTEGER i, j
DIMENSION v(2, mxNode)

title1 = ' '
READ (iUnitV, '(A80)', END = 100, ERR = 100) title1
title2 = ' '
READ (iUnitV, '(A80)', END = 100, ERR = 100) title2
title3 = ' '
READ (iUnitV, '(A80)', END = 100, ERR = 100) title3
READ (iUnitV, * , END = 100, ERR = 100) ((v(j, i), j = 1, 2), i = 1, numNod)
IF(Verbose) WRITE (iUnitVerb, 50) iUnitV, title1, title2, title3
50  FORMAT (/ /' Old velocity solution (initial estimate) was', &
&           ' read from unit',I3,'; titles were:'/3(/' ',A80))
GO TO 900
! ------------------(This section executed only if READ fails)---------
100  IF(Verbose) WRITE (iUnitVerb, 110) iUnitV
110  FORMAT (/ /' UNABLE TO READ OLD VELOCITY SOLUTION FROM UNIT', &
&                I3/ /' VELOCITIES WILL BE INITIALIZED TO ZERO.')
DO 150 i = 1, numNod
	v(1, i) = 0.0D0
	v(2, i) = 0.0D0
150  CONTINUE
! ---------------------------------------------------------------------
900  IF(Verbose) WRITE (iUnitVerb, 999)
999  FORMAT (' --------------------------------------------------', &
&          '-----------------------------')
RETURN
END SUBROUTINE OldVel

SUBROUTINE OneBar (continuum_LRi, &                                                   ! input
&                    geothC, geothM, gradie, &                                          ! input
&                    LRn, LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_eCreep, & ! input
&                    mxEl, numEl, oneKm, tAdiab, &                                      ! input
&                    zBAsth, zMoho, &                                                   ! input
&                    glue)                                                              ! output

!   Calculates "glue" (shear stress required to create one unit of relative
!   horizontal velocity across the lithosphere+asthenosphere mantle layer, down to depth zBAsth).

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: continuum_LRi                                                   ! input
REAL*8, INTENT(IN) :: geothC, geothM, gradie                                           ! input
INTEGER, INTENT(IN) :: LRn                                                             ! input
REAL*8, INTENT(IN) :: LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_eCreep       ! input
INTEGER, INTENT(IN) :: mxEl, numEl                                                     ! input
REAL*8, INTENT(IN) :: oneKm, tAdiab, zBAsth, zMoho                                     ! input
REAL*8, INTENT(OUT) :: glue                                                            ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, layer, level, limit, m
REAL*8 ailog, arg, bi, dz, ecini, gt, t, ta, tg, v, z
!      External argument arrays:
DIMENSION continuum_LRi(mxEl), geothC(4, 7, mxEl), geothM(4, 7, mxEl), &
&           glue(7, mxEl), &
&           LR_set_aCreep(1:2, 0:LRn), LR_set_bCreep(1:2, 0:LRn), LR_set_cCreep(1:2, 0:LRn), LR_set_eCreep(0:LRn), &
&           zMoho(7, mxEl)
!      Internal variables:
INTEGER LRi
REAL*8 t_aCreep(2), t_bCreep(2), t_cCreep(2), t_eCreep
!      Internal arrays:
DIMENSION ailog(2), gt(4)

dz = oneKm
limit = zBAsth / dz + 0.5D0
DO 100 i = 1, numEl
   !retrieve desired rheology for this continuum element:
	LRi = continuum_LRi(i)
	t_aCreep(1:2) = LR_set_aCreep(1:2, LRi)
	t_bCreep(1:2) = LR_set_bCreep(1:2, LRi)
	t_cCreep(1:2) = LR_set_cCreep(1:2, LRi)
	t_eCreep      = LR_set_eCreep(LRi)
   !statements that were formerly outside the loops:
	ecini = -1.0D0 / t_eCreep
	ailog(1) = log(t_aCreep(1)) * ecini
	ailog(2) = log(t_aCreep(2)) * ecini
	DO 90 m = 1, 7
		!Integrate difference in horizontal velocity over depth:
		 v = 0.0D0
		 DO 20 level = 1, limit
			  z = (level - 0.5D0) * dz
			  IF (z < zMoho(m, i)) THEN
				   layer = 1
				   gt(1) = geothC(1, m, i)
				   gt(2) = geothC(2, m, i)
				   gt(3) = geothC(3, m, i)
				   gt(4) = geothC(4, m, i)
			  ELSE
				   layer = 2
				   gt(1) = geothM(1, m, i)
				   gt(2) = geothM(2, m, i)
!                         Note: Quadratic and cubic terms could
!                         cause lithospheric geotherm to have
!                         multiple (nonphysical) intersections
!                         with the adiabat!
				   gt(3) = 0.0D0
				   gt(4) = 0.0D0
			  END IF
			  tg = gt(1) &
&                   + gt(2) * z &
&                   + gt(3) * z * z &
&                   + gt(4) * z * z * z
			  ta = tAdiab + z * gradie
			  t = MIN(tg, ta)
			  t = MAX(t, 200.0D0)
			  bi = (t_bCreep(layer) + t_cCreep(layer) * z) * ecini
			  arg = MAX(ailog(layer) + bi / t, -87.0D0)
			  v = v + dz * EXP(arg)
20            CONTINUE
		 glue(m, i) = 1.0D0 / (v**t_eCreep)
90       CONTINUE
100  CONTINUE
RETURN
END SUBROUTINE OneBar


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

SUBROUTINE PrintK (f, iUnitT, iter, k, & ! input
&                    mxDOF, nDOF, nLB, nUB)

!   Prints out the "k" (or "K", or "stiff") stiffness matrix
!   and the "f" (or "F") forcing vector, for debugging purposes.
!   Typically, it must be printed in blocks, and then
!   these blocks must be pasted together.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DOUBLE PRECISION, INTENT(IN) :: f                                                      ! input
INTEGER, INTENT(IN) :: iUnitT, iter                                                    ! input
DOUBLE PRECISION, INTENT(IN) :: k                                                      ! input
INTEGER, INTENT(IN) :: mxDOF, nDOF, nLB, nUB                                           ! input
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
CHARACTER*4 cnode4
CHARACTER*7 cnode7
CHARACTER*9 text
INTEGER, PARAMETER :: ncol = 14
INTEGER i, i1, i2, iq, irb, j, j1, j2, jcb, l, m, mode, nblock, node
LOGICAL doit
DIMENSION f(mxDOF, 1)
DIMENSION text(ncol), k(nKRows, nRank)

!  Note: 1 CCC + I4 + 'X:' + 14D9.2 = 133 columns.
1  FORMAT(       '1',' Block Row',I2,', Block Column',I2/)
10  FORMAT(       ' ',4X,2X,     14A9)
!  11  FORMAT(/ / / /' ',I4,'X:',1P,14D9.2)
!  12  FORMAT(       ' ',I4,'Y:',1P,14D9.2)
21  FORMAT(/ / / /' ',I4,'X:',   14A9)
22  FORMAT(       ' ',I4,'Y:',   14A9)

WRITE(6, 23) iter
23  FORMAT(' ITERATION ',I5)
nblock = (nDOF + 2) / ncol
IF ((nDOF + 2) > ncol * nblock) nblock = nblock + 1
DO 100 irb = 1, nblock
  DO 90 jcb = 1, nblock
	 i2 = ncol * irb
	 i1 = i2 - ncol + 1
	 j2 = ncol * jcb
	 j1 = j2 - ncol + 1
	 doit = (i1 <= nDOF)    .AND. &
&      (   (j2 > nDOF)      .OR. &
&      ((j2 >= (i1 - nLB)).AND.(j1 <= (i2 + nUB)))    )
	 IF (.NOT. doit) GO TO 90

!   Write header for each block (page):

	 IF(Verbose) WRITE (iUnitVerb, 1) irb, jcb

!   Prepare and WRITE headers over the columns

	 DO 60 j = j1, j2
		m = j - j1 + 1
		IF (j <= nDOF) THEN
		   mode = (j + 1) / 2
		   WRITE (cnode7, '(I7)') mode
		   IF (MOD(j, 2) == 1) THEN
			  text(m) = cnode7//'X:'
		   ELSE
			  text(m) = cnode7//'Y:'
		   END IF
		ELSE
		   text(m) = '         '
		END IF
60        CONTINUE
	 IF(Verbose) WRITE (iUnitVerb, 10) (text(l), l = 1, ncol)
	 DO 80 i = i1, i2

!   Prepare text of a line within the system of equations

		node = (i + 1) / 2
		IF (i <= nDOF) THEN
		   DO 70 j = j1, j2
			  m = j - j1 + 1
			  IF (j <= nDOF) THEN
				!matrix element (i, j)
				 iq = iDiagonal + i - j
				 IF ((j >= (i - nLB)).AND.(j <= (i + nUB))) THEN
					WRITE(text(m), '(1P,D9.2)') k(iq, j)
				 ELSE
					text(m) = ' ------- '
				 END IF
			  ELSE IF (j == (nDOF + 1)) THEN
				 WRITE (cnode4, '(I4)') node
				 IF (MOD(i, 2) == 1) THEN
					text(m) = ' *'//cnode4//'X ='
				 ELSE
					text(m) = ' *'//cnode4//'Y ='
				 END IF
			  ELSE IF (j == (nDOF + 2)) THEN
				 WRITE(text(m), '(1P,D9.2)') f(i, 1)
			  ELSE
				 text(m) = '         '
			  END IF
70              CONTINUE

!   Actually print the line:

		   IF (MOD(i, 2) == 1) THEN
			  IF(Verbose) WRITE (iUnitVerb, 21) node, (text(l), l = 1, ncol)
		   ELSE
			  IF(Verbose) WRITE (iUnitVerb, 22) node, (text(l), l = 1, ncol)
		   END IF
		END IF
80        CONTINUE
90     CONTINUE
100  CONTINUE
IF(Verbose) WRITE (iUnitVerb, 101)
101  FORMAT('1----------------------------------------------------', &
&        '---------------------------')
RETURN
END SUBROUTINE PrintK

SUBROUTINE Pure (alphaT, area, &                       ! input
&                  basal, &
&                  conduc, constr, continuum_LRi, &
&                  delta_rho, detJ, dQdTdA, dXS, dYS, &
&                  elev, etaMax, everyP, &
&                  fault_LRi, fBase, fDip, fLen, fMuMax, &
&                  fPFlt, fPSfer, fArg, geothC, geothM, glue, &
&                  gMean, iCond, iConve, iUnitI, iUnitS, iUnitT, &
&                  LRn, LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, &
&                  LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
&                  maxItr, mxBn, mxDOF, mxEl, mxFEl, &
&                  mxNode, nCond, nDOF, nFl, nLB, nodCon, &
&                  nodeF, nodes, nUB, numEl, numNod, offMax, &
&                  offset, okToQt, oneKm, oVB, pulled, &
&                  radio, radius, rhoBar, rhoH2O, sita, slide, &
&                  tauMax, temLim, title1, &
&                  title2, title3, tLInt, tLNode, trHMax, &
&                  tSurf, vBCArg, vBCMag, visMax, wedge, &
&                  zMNode, zMoho, lastPm, &
&                  v, &                                  ! modify
&                  eRate, eta, fIMuDZ, fPeakS, fSlips, & ! output
&                  sigHB, tauMat, zTranC, zTranF, &
&                  alpha, dv, dVLast, f, fC, fTStar, &   ! work
&                  outVec, k, ipiv, tOfset, &
&                  ThID,SHELLSconv)

!   Create and solve thin-plate, weak-form version of equilibrium to determine
!   horizontal velocity components (using iteration to handle nonlinearities).

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alphaT, area                                                      ! input
DOUBLE PRECISION, INTENT(IN) :: basal                                                   ! input
REAL*8, INTENT(IN) :: conduc, constr                                                    ! input
INTEGER, INTENT(IN) :: continuum_LRi                                                    ! input
REAL*8, INTENT(IN) :: delta_rho, detJ, dQdTdA, dXS, dYS, elev, etaMax                   ! input
LOGICAL, INTENT(IN) :: everyP                                                           ! input
DOUBLE PRECISION, INTENT(IN) :: fBase                                                   ! input
INTEGER, INTENT(IN) :: fault_LRi                                                        ! input
REAL*8, INTENT(IN) :: fDip, fLen, fMuMax, fPFlt, fPSfer, fArg, geothC, &                ! input
  & geothM, glue, gMean                                                                ! input
INTEGER, INTENT(IN) :: iCond, iConve, iUnitI, iUnitS, iUnitT                            ! input
INTEGER, INTENT(IN) :: LRn                                                              ! input
REAL*8, INTENT(IN) :: LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, &         ! input
				   & LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep ! input
INTEGER, INTENT(IN) :: maxItr, mxBn, mxDOF, mxEl, mxFEl, &                              ! input
	 &                mxNode, nCond, nDOF, nFl, nLB, &                                 ! input
	 &                nodCon, nodeF, nodes, nUB, numEl, numNod                         ! input
REAL*8, INTENT(IN) :: offMax, offset, okToQt, oneKm, oVB                                ! input
LOGICAL, INTENT(IN) :: pulled                                                           ! input
REAL*8, INTENT(IN) :: radio, radius, rhoBar, rhoH2O, sita, slide                        ! input
REAL*8, INTENT(IN) :: tauMax, temLim                                                    ! input
CHARACTER*100, INTENT(IN) :: title1, title2, title3                                      ! input
REAL*8, INTENT(IN) :: tLInt, tLNode, trHMax, tSurf, vBCArg, vBCMag, visMax, wedge, &    ! input
				   & zMNode, zMoho                                                     ! input
INTEGER, INTENT(IN) :: lastPm                                                           ! input
DOUBLE PRECISION, INTENT(INOUT) :: v                                                    ! modify
REAL*8, INTENT(OUT) :: eRate, eta, fIMuDZ, fPeakS                                       ! output
LOGICAL, INTENT(OUT) :: fSlips                                                          ! output
REAL*8, INTENT(OUT) :: sigHB, tauMat, zTranC, zTranF                                    ! output
REAL*8 alpha, dv, dVLast                                                 ! work
DOUBLE PRECISION f                                                       ! work
REAL*8 fC, fTStar, outVec                                                ! work
DOUBLE PRECISION k                                                       ! work
INTEGER ipiv                                                             ! work
REAL*8 tOfset                                                            ! work
INTEGER,INTENT(IN) :: ThID
LOGICAL,INTENT(INOUT) :: SHELLSconv
CHARACTER(LEN=100) :: iterFile
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
INTEGER i, ios, iprint, iter, m, memory
REAL*8 dVCorr, scoreA, scoreB, scoreC, scoreD, size1, size2, sumD, sumN, vRMS
LOGICAL valid
DIMENSION alpha(3, 3, 7, mxEl), area(mxEl), &
&           basal(2, mxNode), &
&           continuum_LRi(mxEl), &
&           delta_rho(7, mxEl), detJ(7, mxEl), &
&           dQdTdA(mxNode), &
&           dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
&           dv(2, mxNode), dVLast(2, mxNode), &
&           elev(mxNode), eRate(3, 7, mxEl), eta(7, mxEl), &
&           fault_LRi(mxFEl), &
&           f(mxDOF, 1), fBase(mxDOF), fC(2, 2, 7, mxFEl), &
&           fDip(2, mxFEl), fIMuDZ(7, mxFEl), &
&           fLen(mxFEl), fPeakS(2, mxFEl), &
&           fPFlt(2, 2, 2, 7, mxFEl), &
&           fPSfer(2, 2, 3, 7, mxEl), fSlips(mxFEl), &
&           fArg(2, mxFEl), fTStar(2, 7, mxFEl), &
&           geothC(4, 7, mxEl), geothM(4, 7, mxEl), glue(7, mxEl), &
&           iCond(mxBn), ipiv(nRank), k(nKRows, nRank), &
&           LR_set_fFric(0:LRn), LR_set_cFric(0:LRn), LR_set_Biot(0:LRn), LR_set_Byerly(0:LRn), &
&           LR_set_aCreep(1:2, 0:LRn), LR_set_bCreep(1:2, 0:LRn), LR_set_cCreep(1:2, 0:LRn), LR_set_dCreep(1:2, 0:LRn), LR_set_eCreep(0:LRn), &
&           nodCon(mxBn), nodeF(4, mxFEl), nodes(3, mxEl), offset(mxFEl), &
&           outVec(2, 7, mxEl), oVB(2, 7, mxEl), pulled(7, mxEl), &
&           sigHB(2, 7, mxEl), sita(7, mxEl), tauMat(3, 7, mxEl), &
&           tauMax(2), tOfset(3, 7, mxEl), tLInt(7, mxEl), &
&           tLNode(mxNode), vBCArg(mxBn), &
&           vBCMag(mxBn), v(2, mxNode), &
&           zMNode(mxNode), zMoho(7, mxEl), &
&           zTranC(2, 7, mxEl), zTranF(2, mxFEl)
DIMENSION alphaT(2), conduc(2), &
&           radio(2),  rhoBar(2), temLim(2)

IF (lastPm /= 999) THEN
  write(ErrorMsg,'(A)') "WRONG NUMBER OF ARGUMENTS IN CALL TO -Pure-!"
  call FatalError(ErrorMsg,ThID)
END IF

!   Initialize strain rate and vertical integrals of relative stress
!   for the triangular continuum elements:

CALL EDot (dXS, dYS, & ! input
&            fPSfer, mxEl, &
&            mxNode, nodes, numEl, radius, sita, v, &
&            eRate)      ! output
DO 20 m = 1, 7
	DO 10 i = 1, numEl
		 sigHB(1, m, i) = 0.0D0
		 sigHB(2, m, i) = 0.0D0
		 tauMat(1, m, i) = 0.0D0
		 tauMat(2, m, i) = 0.0D0
		 tauMat(3, m, i) = 0.0D0
10       CONTINUE
20  CONTINUE

CALL Viscos (alphaT, &                              ! input
&              continuum_LRi, &
&              delta_rho, &
&              eRate, gMean, geothC, geothM, &
&              LRn, LR_set_cFric, LR_set_Biot, &
&              LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
&              mxEl, numEl, rhoBar, rhoH2O, &
&              sigHB, tauMat, temLim, tLInt, &
&              visMax, zMoho, &
&              alpha, scoreC, scoreD, tOfset, zTranC) ! output

CALL TauDef (alpha, eRate, mxEl, numEl, tOfset, & ! input
&              tauMat)                              ! output

!   Initialize slip rate and vertical integrals of relative stress
!   for the linear fault elements

DO 30 i = 1, nFl
	zTranF(1, i) = (zMNode(nodeF(1, i)) + &
&                      zMNode(nodeF(2, i))) / 6.0D0
	zTranF(2, i) = (tLNode(nodeF(1, i)) + &
&                      tLNode(nodeF(2, i))) / 6.0D0
30  CONTINUE
CALL Mohr (alphaT, conduc, constr, &                     ! input
&            continuum_LRi, &
&            dQdTdA, elev, &
&            fault_LRi, fDip, fMuMax, &
&            fPFlt, fArg, gMean, &
&            LRn, LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, &
&            LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
&            mxEl, mxFEl, mxNode, nFl, nodeF, &
&            offMax, offset, &
&            oneKm, radio, rhoH2O, rhoBar, &
&            slide, tauMax, &
&            tLNode, tSurf, v, wedge, &
&            zMNode, &
&            zTranF, &                             ! modify
&            fC, fIMuDZ, fPeakS, fSlips, fTStar)   ! output

!   Create "iteration permit" file

WRITE(iterFile,"('iteration_',I0,'_permit.txt')") ThID
OPEN (UNIT = iUnitI, FILE = trim(iterFile))
WRITE (iUnitI, 98)
98  FORMAT('If you delete this file, Shells will' &
&       /'stop at the end of the next iteration' &
&       /'and report the current unconverged solution.')
CLOSE (UNIT = iUnitI)

!   Major Iteration Loop of the Entire Program !!!!!

IF(Verbose) WRITE(iUnitVerb, 99)
WRITE(iUnitT, 99)

99  FORMAT (/ /' Iteration history:'/ &
&'              ', &
&'                                                 Relative'/ &
&'              ', &
&'                            Corre-     Maximum       mean'/ &
&'              ', &
&'                   Relative lation  vertically vertically'/ &
&'              ', &
&'           Maximum     mean   with  integrated integrated'/ &
&' Iter-        ', &
&'   RMS    velocity velocity   last      stress     stress'/ &
&' ation      ve', &
&'locity      change   change change       error      error'/)

DO 1000 iter = 1, maxItr
	memory = iter
	CALL THOnB (basal, continuum_LRi, etaMax, & ! input
&                  fPSfer, glue, iConve, &
&                  LRn, LR_set_eCreep, &
&                  mxEl, mxNode, nodes, numEl, &
&                  oVB, pulled, trHMax, v, &
&                  eta, sigHB, &                   ! output
&                  outVec)                         ! work
	CALL Viscos (alphaT, &                              ! input
&                   continuum_LRi, &
&                   delta_rho, &
&                   eRate, gMean, geothC, geothM, &
&                   LRn, LR_set_cFric, LR_set_Biot, &
&                   LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
&                   mxEl, numEl, rhoBar, rhoH2O, &
&                   sigHB, tauMat, temLim, tLInt, &
&                   visMax, zMoho, &
&                   alpha, scoreC, scoreD, tOfset, zTranC) ! output
	CALL Mohr (alphaT, conduc, constr, &                     ! input
&                 continuum_LRi, &
&                 dQdTdA, elev, &
&                 fault_LRi, fDip, fMuMax, &
&                 fPFlt, fArg, gMean, &
&                 LRn, LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, &
&                 LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
&                 mxEl, mxFEl, mxNode, nFl, nodeF, &
&                 offMax, offset, &
&                 oneKm, radio, rhoH2O, rhoBar, &
&                 slide, tauMax, &
&                 tLNode, tSurf, v, wedge, &
&                 zMNode, &
&                 zTranF, &                             ! modify
&                 fC, fIMuDZ, fPeakS, fSlips, fTStar)   ! output
	IF (iter > 1) THEN
		 iprint = iter - 1
		 IF (iprint == 1) THEN
			  IF(Verbose) WRITE(iUnitVerb, 101) iprint, vRMS, &
&                          scoreA, scoreB, scoreC, scoreD
			  WRITE(iUnitT, 101) iprint, vRMS, &
&                          scoreA, scoreB, scoreC, scoreD
101                 FORMAT(' ',I5,1P,E14.6,E12.4,0P,F9.6,'   ----', &
&                              1P,E12.4,0P,F11.6)
		 ELSE
			  IF(Verbose) WRITE(iUnitVerb, 102) iprint, vRMS, &
&                          scoreA, scoreB, dVCorr, scoreC, scoreD
			  WRITE(iUnitT, 102) iprint, vRMS, &
&                          scoreA, scoreB, dVCorr, scoreC, scoreD
102                 FORMAT(' ',I5,1P,E14.6,E12.4,0P,F9.6,F7.2, &
&                              1P,E12.4,0P,F11.6)
		 END IF
	END IF
	DO 110 i = 1, numNod
		 dVLast(1, i) = dv(1, i)
		 dVLast(2, i) = dv(2, i)
110       CONTINUE
	CALL FEM (alpha, area, constr, detJ, &  ! input
&                dXS, dYS, eta, &
&                everyP, fBase, fC, fDip, &
&                fIMuDZ, fLen, fPFlt, fPSfer, fArg, &
&                fTStar, iCond, iUnitS, iUnitT, &
&                mxBn, mxDOF, mxEl, mxFEl, mxNode, &
&                nCond, nDOF, nFl, nLB, nodCon, nodeF, &
&                nodes, nUB, numEl, numNod, &
&                oVB, pulled, radius, sita, &
&                title1, title2, title3, tOfset, trHMax, &
&                vBCArg, vBCMag, wedge, &
&                999, &
&                eRate, v, &                   ! modify
&                dv, scoreA, scoreB, tauMat, & ! output
&                f, k, ipiv)                   ! work
	vRMS = 0.0D0
	DO 105 i = 1, numNod
		 vRMS = vRMS + v(1, i)**2 + v(2, i)**2
105       CONTINUE
	vRMS = SQRT(vRMS / (1.0D0 * numNod))
	IF (iter >= 2) THEN
		 sumN = 0.0D0
		 sumD = 0.0D0
		 DO 107 i = 1, numNod
			  size1 = SQRT(dv(1, i)**2 + dv(2, i)**2)
			  size2 = SQRT(dVLast(1, i)**2 + dVLast(2, i)**2)
			  sumN = sumN + dv(1, i) * dVLast(1, i) + &
&                              dv(2, i) * dVLast(2, i)
			  sumD = sumD + size1 * size2
107            CONTINUE
		 IF (sumD > 0.0D0) THEN
			  dVCorr = sumN / sumD
		 ELSE
			  dVCorr = 0.0D0
		 END IF
	END IF
	IF (scoreB <= okToQt) THEN
		 IF(Verbose) WRITE(iUnitVerb, 109) iter, vRMS, scoreA, scoreB, dVCorr
		 WRITE(iUnitT, 109) iter, vRMS, scoreA, scoreB, dVCorr
109            FORMAT(' ',I5,1P,E14.6,E12.4,0P,F9.6,F7.2)
		 IF(Verbose) WRITE(iUnitVerb, 998)
		 WRITE(iUnitT, 998)
998            FORMAT (' CONVERGED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
!                Open file again, just in order to delete it:
		 OPEN (UNIT = iUnitI, FILE = trim(iterFile), &
&                 STATUS = 'OLD', IOSTAT = ios)
		 CLOSE (UNIT = iUnitI, STATUS = 'DELETE')
		 RETURN
	END IF

!           Check whether iteration permit still exists:

	OPEN (UNIT = iUnitI, FILE = trim(iterFile), STATUS = 'OLD', IOSTAT = ios)
	valid = (ios == 0)
	IF (valid) CLOSE (UNIT = iUnitI)
	IF (.NOT.valid) GO TO 1001

1000  CONTINUE
1001  IF(Verbose) WRITE(iUnitVerb, 109) memory, vRMS, scoreA, scoreB, dVCorr
WRITE(iUnitT, 109) memory, vRMS, scoreA, scoreB, dVCorr
IF (valid) THEN
    SHELLSconv = .False.
	IF(Verbose) WRITE(iUnitVerb, 1002)
	WRITE(iUnitT, 1002)
1002       FORMAT(' ITERATION LIMIT REACHED BEFORE CONVERGENCE.')
!           Open file again just in order to delete it:
	OPEN (UNIT = iUnitI, FILE = trim(iterFile), STATUS = 'OLD', &
&            IOSTAT = ios)
	CLOSE (UNIT = iUnitI, STATUS = 'DELETE')
ELSE
	IF(Verbose) WRITE(iUnitVerb, 1003)
	WRITE(iUnitT, 1003)
1003       FORMAT(' ITERATION WAS STOPPED BY OPERATOR.')
END IF
RETURN
END SUBROUTINE Pure

SUBROUTINE ReadBC(brief, fDip, iPVRef, iUnitB, iUnitD, iUnitT, &   ! input
&                   mxBn, mxFEl, mxNode, names, nFl, &
&                   nodeF, nPlate, nRealN, numNod, n1000, omega, &
&                   radius, slide, sphere, trHMax, xNode, &
&                   yNode, &
&                   nCond, &                                         ! modify
&                   iCond, nodCon, savTag, title2, vBCArg, vBCMag, & ! output
&                   iEdge, r2Edge, xEdge, yEdge)                     ! work

!   Read in velocity boundary conditions from unit iUnitB,
!      with comments sent to device iUnitT.
!   One option is to have the velocity boundary conditions set by
!      subprogram -EdgeVs-, which would read unit iUnitD.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
LOGICAL, INTENT(IN) :: brief                                                           ! input
REAL*8, INTENT(IN) :: fDip                                                             ! input
INTEGER, INTENT(IN) :: iPVRef, iUnitB, iUnitD, iUnitT, mxBn, mxFEl, mxNode             ! input
CHARACTER*2, INTENT(IN) :: names                                                       ! input
INTEGER, INTENT(IN) :: nFl, nodeF, nPlate, nRealN, numNod, n1000                       ! input
REAL*8, INTENT(IN) :: omega, radius                                                    ! input
REAL*8, INTENT(IN) :: slide                                                            ! input
LOGICAL, INTENT(IN) :: sphere                                                          ! input
REAL*8, INTENT(IN) :: trHMax, xNode, yNode                                             ! input
INTEGER, INTENT(INOUT) :: nCond                                                        ! modify
INTEGER, INTENT(OUT) :: iCond, nodCon                                                  ! output
CHARACTER*2, INTENT(OUT) :: savTag                                                     ! output
CHARACTER*100, INTENT(OUT) :: title2                                                    ! output
REAL*8, INTENT(OUT) :: vBCArg, vBCMag                                                  ! output
INTEGER iEdge                                                           ! work
REAL*8 r2Edge, xEdge, yEdge                                             ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CHARACTER*2 namTag
INTEGER i, ios, n, nFixed, node, nodExp, nRead, number
REAL*8 phi, pLat, pLon, theta, vAz, vMag
LOGICAL allOK, readIt
DIMENSION fDip(2, mxFEl), iCond(mxBn), iEdge(mxBn), nodCon(mxBn), &
&           nodeF(4, mxFEl), r2Edge(mxBn), savTag(mxBn), &
&           vBCArg(mxBn), vBCMag(mxBn), &
&           xEdge (mxBn), yEdge (mxBn), &
&           xNode(mxNode), yNode(mxNode)
DIMENSION names(nPlate), omega(3, nPlate)

IF(Verbose) WRITE(iUnitVerb, 10) iUnitB
IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 10) iUnitB
10  FORMAT(/ /' Attempting to read boundary conditions from unit', &
&      I3/)

title2 = ' '

READ (iUnitB, 12, IOSTAT = ios) title2

IF (ios /= 0) THEN
   write(ErrorMsg,'(A)') "ERROR in -ReadBC-: File not found, or empty, too short, or defective."
   call FatalError(ErrorMsg,ThID)
END IF
12  FORMAT (A80)
IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 15) title2
15  FORMAT (/' Title for set of boundary conditions ='/' ',A80)

allOK = .TRUE.
nFixed = 0
readIt = .FALSE.

!   During first pass, don't print table entries (incomplete)

!   Begin indefinate loop (at least nCond entries required, but up to
!   numNod entries might appear!

i = 0
30       i = i + 1
	IF (i <= nCond) THEN
		 nodExp = nodCon(i)
	ELSE
		 nodExp = 0
	END IF

	READ (iUnitB, * , IOSTAT = ios, END = 100) number, node, iCond(i)
	IF (ios == 24) GO TO 100
!          (jumping out of loop due to EOF condition)
	IF (ios /= 0) THEN
	     write(ErrorMsg,'(A)') "ERROR in -ReadBC-: File not found, or empty, too short, or defective."
		 call FatalError(ErrorMsg,ThID)
	END IF

	IF (number /= i) THEN
		 IF(Verbose) WRITE (iUnitVerb, 40) number, i
40            FORMAT (' ILLEGAL ORDERING OF BOUNDARY CONDITIONS:'/ &
&             ' READ CONDITION #',I6,' WHEN EXPECTING #',I6,'.'/ &
&            ' SUGGESTION: EDIT LOG FILE TABLE TO MAKE B.C. FILE.')
		 allOK = .FALSE.
	END IF

	IF (node > nRealN) node = nRealN + (node - n1000)
	IF ((node <= 0).OR.(node > numNod)) THEN
		 IF (node > nRealN) node = n1000 + (node - nRealN)
		 WRITE(iUnitT, 45) node
45            FORMAT(' ILLEGAL NODE NUMBER IN BOUNDARY', &
&                ' CONDITIONS:',I6)
		 allOK = .FALSE.
	END IF

	IF ((nodexp > 0).AND.(node /= nodexp)) THEN
		 IF (node > nRealN) node = n1000 + (node - nRealN)
		 IF (nodexp > nRealN) nodexp = n1000 + (nodexp - nRealN)
		 WRITE(iUnitT, 47) node, nodexp
47            FORMAT(/' BOUNDARY CONDITIONS INPUT IN WRONG ORDER;'/ &
&              ' (SEE LIST PREVIOUSLY WRITTEN IN OUTPUT ABOVE)', &
&              /'    ',I6,' WAS READ WHEN EXPECTING ',I6)
		 allOK = .FALSE.
	END IF
	IF (nodexp == 0) nodCon(i) = node

	IF ((iCond(i) == 0).OR.(iCond(i) == -1)) THEN
!                No action needed for free nodes (of either type)
	ELSE IF (iCond(i) == 1) THEN
		 BACKSPACE iUnitB
		 READ (iUnitB, * ) number, node, iCond(i), vMag, vAz
		 nFixed = nFixed + 1
		 vBCMag(i) = vMag
		 vBCArg(i) = (180.0D0 - vAz) * 0.0174532925199433D0
	ELSE IF (iCond(i) == 2) THEN
		 BACKSPACE iUnitB
		 READ (iUnitB, * ) number, node, iCond(i), vMag, vAz
		 nFixed = nFixed + 2
		 vBCMag(i) = vMag
		 vBCArg(i) = (180.0D0 - vAz) * 0.0174532925199433D0
	ELSE IF (iCond(i) == 3) THEN
		 readIt = .TRUE.
		 nFixed = nFixed + 1
	ELSE IF (iCond(i) == 4) THEN
		 readIt = .TRUE.
		 nFixed = nFixed + 2
	ELSE IF (iCond(i) == 5) THEN
		 BACKSPACE iUnitB
		 READ (iUnitB, * ) number, node, iCond(i), namTag
		 savTag(i) = namTag
		 CALL Euler (namTag, node, & ! input
&                       iPVRef, names, nPlate, omega, &
&                       iUnitT, radius, &
&                       mxNode, xNode, yNode, &
&                       vAz, vMag)      ! output
		 vBCMag(i) = vMag
		 vBCArg(i) = (180.0D0 - vAz) * 0.0174532925199433D0
		 nFixed = nFixed + 2
	ELSE
		 WRITE(iUnitT, 95) iCond(i)
95            FORMAT(' ILLEGAL TYPE OF BOUNDARY', &
&                ' CONDITION:',I6)
		 allOK = .FALSE.
	END IF

!      end of indefinate loop:
IF (i < numNod) GO TO 30
100  CONTINUE
nRead = i - 1
IF (nRead < nCond) THEN
	     write(ErrorMsg,'(A)') "ERROR in -ReadBC-: File not found, or empty, too short, or defective."
		 call FatalError(ErrorMsg,ThID)
ELSE IF (nRead > nCond) THEN
	nCond = MIN(nRead, numNod)
END IF

!   Do we need to complete table (by filling-in iCond=3/4 nodes)?

IF (readIt) THEN
	CALL EdgeVs (fDip, iPVRef, iUnitD, iUnitT, mxBn, mxNode, & ! input
&                   mxFEl, names, nCond, nFl, nodCon, nodeF, nPlate, &
&                   omega, radius, slide, sphere, xNode, yNode, &
&                   iCond, vBCArg, vBCMag, &                      ! modify
&                   iEdge, r2Edge, xEdge, yEdge)                  ! work
END IF

!   Now, it's OK to print the table:

IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 101) nCond
101  FORMAT(/' There are ',I6,' nodes with boundary conditions.'/ &
&     ' When describing the kind of boundary condition,', &
&     ' the code is:'/ &
&     '             -1 = no velocity constraint (ridge adjacent).'/ &
&     '              0 = no velocity constraint (weak margin).'/ &
&     '              1 = fix velocity in specified direction;'/ &
&     '                  perpendicular component remains free.'/ &
&     '              2 = fix velocity in specified direction;'/ &
&     '                  perpendicular component set to zero.'/ &
&     '              3 = fix PB2002 component at PB2002'/ &
&     '                  velocity value.'/ &
&     '              4 = fix both components at PB2002'/ &
&     '                  velocity value.'/ &
&     '              5 = fix velocity to that of named plate;'/ &
&     '                  azimuth is also fixed.'// &
&'    BC#  Node (E.lon) (N.lat)  Code    Velocity     Azimuth (deg' &
&,'rees clockwise from North)')
!     (' ', I6,   I6, F8.2,   F8.2,     I6,   1P,E12.3,    0P,F12.1)
DO 200 i = 1, nCond
	n = nodCon(i)
	IF (n <= nRealN) THEN
		 node = n
	ELSE
		 node = n1000 + n - nRealN
	END IF
	theta = xNode(n)
	phi = yNode(n)
	pLon = 57.2957795130823D0 * phi
	pLat = 90.0D0 - theta * 57.2957795130823D0
	IF (iCond(i) == -1) THEN
		 IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 159) i, node, &
&                pLon, pLat, iCond(i)
159            FORMAT(' ',2I6,2F8.2,I6,'        FREE','        FREE' &
&                 ,' (RIDGE ADJACENT)')
	ELSE IF (iCond(i) == 0) THEN
		 IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 160) i, node, &
&                pLon, pLat, iCond(i)
160            FORMAT(' ',2I6,2F8.2,I6,'        FREE','        FREE' &
&                 ,' (WEAK MATERIAL ADJACENT)')
	ELSE IF (iCond(i) == 1) THEN
		 vAz = 180.0D0 - vBCArg(i) * 57.2957795130823D0
		 IF (vaz < 0.0D0) vaz = vaz + 360.0D0
		 IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 161) i, node, &
&                pLon, pLat, iCond(i), vBCMag(i), vaz
161            FORMAT(' ',2I6,2F8.2,I6,1P,E12.3,0P,F12.1,' (PERPEN' &
&                  ,'DICULAR COMPONENT FREE)')
	ELSE IF (iCond(i) == 2) THEN
		 vAz = 180.0D0 - vBCArg(i) * 57.2957795130823D0
		 IF (vaz < 0.0D0) vaz = vaz + 360.0D0
		 IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 162) i, node, &
&                pLon, pLat, iCond(i), vBCMag(i), vaz
162            FORMAT(' ',2I6,2F8.2,I6,1P,E12.3,0P,F12.1,' (NO ' &
&                  ,'PERPENDICULAR COMPONENT)')
	ELSE IF (iCond(i) == 3) THEN
		 vAz = 180.0D0 - vBCArg(i) * 57.2957795130823D0
		 IF (vaz < 0.0D0) vaz = vaz + 360.0D0
		 IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 161) i, node, &
&                pLon, pLat, iCond(i), vBCMag(i), vaz
	ELSE IF (iCond(i) == 4) THEN
		 vAz = 180.0D0 - vBCArg(i) * 57.2957795130823D0
		 IF (vaz < 0.0D0) vaz = vaz + 360.0D0
		 IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 162) i, node, &
&                pLon, pLat, iCond(i), vBCMag(i), vaz
	ELSE IF (iCond(i) == 5) THEN
		 vAz = 180.0D0 - vBCArg(i) * 57.2957795130823D0
		 IF (vaz < 0.0D0) vaz = vaz + 360.0D0
		 IF (.NOT.brief .AND. Verbose) WRITE (iUnitVerb, 165) i, node, &
&                pLon, pLat, iCond(i), savTag(i), vBCMag(i), vaz
165            FORMAT(' ',2I6,2F8.2,I3,':',A2,1P,E12.3,0P,F12.1, &
&                  ' (NO PERPENDICULAR COMPONENT)')
	END IF
200  CONTINUE

IF ((nFixed < 3).AND.(.NOT.sphere).AND.(trHMax <= 0.0D0)) THEN
	allOK = .FALSE.
	write(ErrorMsg,'(A/,A/,A/,A)') "INSUFFICIENT BOUNDARY CONDITIONS. EVERY PROBLEM REQUIRES THAT AT LEAST 3 DEGREES" , &
&               "OF FREEDOM BE CONSTRAINED, TO PREVENT NONUNIQUENESS OF THE SOLUTION WITH" , &
&               "RESPECT TO TRANSLATION AND/OR ROTATION. YOU HAVE CONSTRAINED ONLY',I2,' DEGREES OF" , &
&               "FREEDOM; ADD MORE CONSTRAINED NODES."
	call FatalError(ErrorMsg,ThID)
END IF

IF (Verbose) WRITE (iUnitVerb, 999)
999  FORMAT (' --------------------------------------------------', &
&          '-----------------------------')
RETURN
END SUBROUTINE ReadBC


SUBROUTINE Result (alphaT, area, comp, detJ, elev, eRate, everyP, & ! input
&                    fault_LRi, &
&                    fDip, fIMuDZ, fPFlt, &
&                    fPeakS, fPSfer, fSlips, fArg, &
&                    geothC, iUnitQ, iUnitS, iUnitT, &
&                    log_node_velocities, &
&                    log_element_dynamics, &
&                    log_fault_dynamics, &
&                    LRn, LR_set_fFric, &
&                    mxDOF, mxEl, mxFEl, mxNode, names, nFl, &
&                    nodeF, nodes, nPlate, nRealN, numEl, numNod, &
&                    n1000, oneKm, &
&                    radius, rhoAst, rhoBar, rhoH2O, sigHB, &
&                    tauMat, tauMax, tauZZI, title1, title2, &
&                    title3, tLInt, tLNode, &
&                    v, wedge, whichP, xNode, yNode, &
&                    zMNode, zMoho, zTranC, zTranF, &
&                    torqBS, torqCL, torqFS, torqLP, torqMD, &        ! output
&                    torqSS, torqVB)

!   Output the solution:
!      -Node velocities to unit iUnitS,
!      -Descriptive tables to unit iUnitT:
!        * nodal velocities table
!        * element properties table
!        * fault properties table
!        * single-plate torque-balance report
!   And, compute the various components of the (balanced)
!      vector torques on each plate.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alphaT, area, comp, detJ, elev, eRate                              ! input
LOGICAL, INTENT(IN) :: everyP                                                            ! input
INTEGER, INTENT(IN) :: fault_LRi                                                         ! input
REAL*8, INTENT(IN) :: fDip, fIMuDZ, fPFlt, fPeakS, fPSfer                                ! input
LOGICAL, INTENT(IN) :: fSlips                                                            ! input
REAL*8, INTENT(IN) :: fArg, geothC                                                       ! input
INTEGER, INTENT(IN) :: iUnitQ, iUnitS, iUnitT                                            ! input
LOGICAL, INTENT(IN) :: log_node_velocities, log_element_dynamics, log_fault_dynamics     ! input
INTEGER, INTENT(IN) :: LRn                                                               ! input
REAL*8, INTENT(IN) :: LR_set_fFric                                                       ! input
INTEGER, INTENT(IN) :: mxDOF, mxEl, mxFEl, mxNode                                        ! input
CHARACTER*2, INTENT(IN) :: names                                                         ! input
INTEGER, INTENT(IN) :: nFl, nodeF, nodes, nPlate, nRealN, numEl, numNod, n1000           ! input
REAL*8, INTENT(IN) :: oneKm, radius, rhoAst, rhoBar, rhoH2O, sigHB, tauMat, tauMax, tauZZI ! input
CHARACTER*100, INTENT(IN) :: title1, title2, title3                                       ! input
REAL*8, INTENT(IN) :: tLInt, tLNode                                                      ! input
DOUBLE PRECISION, INTENT(IN) :: v                                                        ! input
REAL*8, INTENT(IN) :: wedge                                                              ! input
INTEGER, INTENT(IN) :: whichP                                                            ! input
REAL*8, INTENT(IN) :: xNode, yNode, zMNode, zMoho, zTranC, zTranF                        ! input
DOUBLE PRECISION, INTENT(OUT) :: torqBS, torqCL, torqFS, torqLP, torqMD, torqSS, torqVB  ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DOUBLE PRECISION points, weight
COMMON / S1S2S3 / points
COMMON / WgtVec / weight
DIMENSION points(3, 7), weight(7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, ip, iPlate, ix, iy, j, j1, j2, j3, j4, jb, jm, k, l, m, n, nInSum
REAL*8 angle, azim, azimHS, azimut, &
	& close, crossx, crossy, dip, dot_E, dot_N, dot_S, du, dv, &
	& e1, e2, equat, equat2, exx, exy, eyy, ezz, &
	& factor, fbsx, fbsy, ffsx, ffsy, flpx, flpy, fmdx, fmdy, &
	& force_bs, force_fs, force_lp, force_md, force_ss, force_vb, &
	& fpMax, fvbx, fvby, &
	& height, hors, latitu, longit, &
	& phi, pLat, pLon, plunge, pPhi, relV, rhoC, rLt, rvec, &
	& shear, sighx, sighy, sinist, sNet, stheta, &
	& t1, t2, tbsLat, tbsLon, tclLat, tclLon, &
	& tequat, tfsLat, tfsLon, theta, tlav, tlpLat, tlpLon, &
	& tmdLat, tmdLon, tmdMag, tMid, tssLat, tssLon, &
	& tTheta, tvbLat, tvbLon, tvbMag, tvec, twiLat, twiLon, twiMag, &
	& twistV, txx, txy, tyy, tzz, &
	& u1x, u1y, u2x, u2y, unitx, unity, uPhi, uTheta, uvec, vE, vMag, vN, vUpDip, vZ, &
	& zmav, zTranM, zTranS
REAL*8 t_fFric
DOUBLE PRECISION tbsMag, tclMag, tfsMag, tlpMag, tssMag
DOUBLE PRECISION ddot, dequat, dvec, length, sumNod
DIMENSION alphaT(2), rhoBar(2), tauMax(2)
DIMENSION area(mxEl), &
&           comp(6, mxDOF), &
&           detJ(7, mxEl), &
&           elev(mxNode),  eRate(3, 7, mxEl), &
&           fault_LRi(mxFEl), &
&           fDip(2, mxFEl), fIMuDZ(7, mxFEl), fPeakS(2, mxFEl), &
&           fPFlt(2, 2, 2, 7, mxFEl), &
&           fPSfer(2, 2, 3, 7, mxEl), &
&           fSlips(mxFEl), fArg(2, mxFEl), &
&           geothC(4, 7, mxEl), &
&           LR_set_fFric(0:LRn), &
&           names(nPlate), nodeF(4, mxFEl), nodes(3, mxEl), &
&           sigHB(2, 7, mxEl), sumNod(3), &
&           tauMat(3, 7, mxEl), tauZZI(7, mxEl), &
&           tLInt(7, mxEl), tLNode(mxNode), &
&           torqBS(3, nPlate), torqCL(3, nPlate), torqFS(3, nPlate), &
&           torqLP(3, nPlate), torqMD(3, nPlate), torqSS(3, nPlate), &
&           torqVB(3, nPlate), &
&           v(2, mxNode), &
&           whichP(mxNode), &
&           xNode(mxNode), yNode(mxNode), &
&           zMNode(mxNode), zMoho(7, mxEl), &
&           zTranC(2, 7, mxEl), zTranF(2, mxFEl)
!   Cartesian vectors:
DIMENSION force_bs(3), force_fs(3), force_lp(3), force_md(3), &
&           force_ss(3), force_vb(3), &
&           dvec(3), rvec(3), tvec(3), twistV(3), uphi(3), utheta(3), &
&           uvec(3)

IF (.NOT.everyP) THEN
	WRITE (iUnitS, 10) title1
	WRITE (iUnitS, 10) title2
	WRITE (iUnitS, 10) title3
10       FORMAT (A80)
	WRITE (iUnitS, 20) ((v(k, i), k = 1, 2), i = 1, numNod)
20       FORMAT (1P,4D20.12)
END IF
!------------------------End of report on unit iUnitS---------------
!------------------------Begin writing to unit iUnitT---------------

!  Velocities at nodes:

IF (log_node_velocities .AND. Verbose) WRITE (iUnitVerb, 30)
30  FORMAT(/ /' Velocities of the nodes:'/ &
&           '                                                    ', &
&           '     Azimuth    East   North'/ &
&           '                                                    ', &
&           '    (degrees    long.    lat.'/ &
&           '  Node      East-component     North-component Magni', &
&           'tude from North)'/)
sumNod(1) = 0.0D0
sumNod(2) = 0.0D0
sumNod(3) = 0.0D0
DO 100 i = 1, numNod
	ip = i
	IF (i > nRealN) ip = n1000 + (i - nRealN)
	theta = xNode(i)
	phi = yNode(i)
	sumNod(1) = sumNod(1) + SIN(theta) * COS(phi)
	sumNod(2) = sumNod(2) + SIN(theta) * SIN(phi)
	sumNod(3) = sumNod(3) + COS(theta)
	pLon = phi * 57.2957795130823D0
	pLat = 90.0D0 - theta * 57.2957795130823D0
	vE = v(2, i)
	vN = -v(1, i)
	azimut = ATan2F(ve, vn) * 57.2957795130823D0
	IF (azimut < 0.0D0) azimut = azimut + 360.0D0
	vMag = SQRT(v(1, i)**2 + v(2, i)**2)
	IF (log_node_velocities .AND. Verbose) WRITE (iUnitVerb, 40) &
&          ip, vE, vN, vMag, azimut, pLon, pLat
40       FORMAT(' ',I5,1P,2D20.12,E10.2,0P,3F8.2)
100  CONTINUE

!  Triangular continuum element properties at their centers:

IF (log_element_dynamics .AND. Verbose) WRITE (iUnitVerb, 110)
110  FORMAT (/ /' Continuum element properties (at center points):'/ &
&     /'                     E1=most   E2=most Isostatic  Vertic', &
&'al  Vertical  Vertical  Brittle/  Brittle/    Basal    Basal' &
&      /' Element   Azimuth  compress.   extens.   uplift  integr', &
&'al  integral  integral   ductile   ductile    shear    shear', &
&   '    East   North' &
&      /'  number     of E1      rate      rate      rate  of(Sz+', &
&'P0) of(S1+P0) of(S2+P0) in crust in mantle   stress  azimuth', &
&   '    long.   lat.'/)
120  FORMAT (' ',I7,F10.2,1P,8E10.2,E9.1,0P,F9.2,2F8.2)
121  FORMAT (' ',I7,F10.2,1P,7E10.2,'  --------',E9.1,0P,F9.2,2F8.2)
122  FORMAT (' ',I7,F10.2,1P,6E10.2,'  --------',E10.2,E9.1,0P,F9.2, &
&             2F8.2)
m = 1
DO 200 i = 1, numEl
	exx = eRate(1, m, i)
	eyy = eRate(2, m, i)
	exy = eRate(3, m, i)
	CALL Prince (exx, eyy, exy, &            ! input
&                   e1, e2, u1x, u1y, u2x, u2y) ! output
	azim = 180.0D0 - ATan2F(u1y, u1x) * 57.2957795130823D0
	IF (azim < 0.0D0) azim = azim + 360.0D0
	ezz = -(exx + eyy)
	tMid = geothC(1, m, i) + geothC(2, m, i) * zMoho(m, i) / 2.0D0 + &
&             geothC(3, m, i) * (zMoho(m, i) / 2.)**2
	rhoC = rhoBar(1) * (1.0D0 - alphaT(1) * tMid)

!         Interpolate height, position to element center:
	height = 0.0D0
	DO 140 l = 1, 3
		 tvec(l) = 0.0D0
140       CONTINUE
	DO 150 k = 1, 3
		 n = nodes(k, i)
		 height = height + elev(n) / 3.0D0
		 theta = xNode(n)
		 phi = yNode(n)
		 equat = SIN(theta)
		 uvec(1) = equat * COS(phi)
		 uvec(2) = equat * SIN(phi)
		 uvec(3) = COS(theta)
		 DO 149 l = 1, 3
			  tvec(l) = tvec(l) + uvec(l)
149            CONTINUE
150       CONTINUE
	equat2 = tvec(1)**2 + tvec(2)**2
	IF (equat2 == 0.0D0) THEN
		 pLon = 0.0D0
		 IF (tvec(3) > 0.0D0) THEN
			  pLat = 90.0D0
		 ELSE
			  pLat = -90.0D0
		 END IF
	ELSE
		 equat = SQRT(equat2)
		 pLat = 57.2957795130823D0 * ATan2F(tvec(3), equat)
		 pLon = 57.2957795130823D0 * ATan2F(tvec(2), tvec(1))
	END IF

	IF (height > 0.0D0) THEN
		 factor = (rhoAst - rhoc) / rhoAst
	ELSE
		 factor = (rhoAst - rhoc) / (rhoAst - rhoH2O)
	END IF
	vZ = ezz * zMoho(m, i) * factor
	txx = tauMat(1, m, i) + tauZZI(m, i)
	tyy = tauMat(2, m, i) + tauZZI(m, i)
	txy = tauMat(3, m, i)
	tzz = tauZZI(m, i)
	CALL Prince (txx, tyy, txy, &            ! input
&                   t1, t2, u1x, u1y, u2x, u2y) ! output
	zTranS = zTranC(1, m, i)
	sighx = sigHB(1, m, i)
	sighy = sigHB(2, m, i)
	sTheta = 180.0D0 - 57.2957795130823D0 * ATan2F(sighy, sighx)
	IF (sTheta >= 360.0D0) sTheta = sTheta - 360.0D0
	IF (sTheta < 0.0D0)   sTheta = sTheta + 360.0D0
	shear = SQRT(sighx**2 + sighy**2)
	IF ((tLInt(m, i) > 0.0D0).AND. &
&          (zTranC(2, m, i) > (0.1D0 * oneKm))) THEN
		 zTranM = zMoho(m, i) + zTranC(2, m, i)
		 IF ((zTrans / zMoho(m, i)) > 0.97D0) THEN
			  IF (log_element_dynamics .AND. Verbose) &
			      &  WRITE (iUnitVerb, 122) i, azim, e1, e2, vz, &
&                             tzz, t1, t2,       zTranm, shear, stheta, &
&                                   pLon, pLat
		 ELSE
			  IF (log_element_dynamics .AND. Verbose) &
			      &  WRITE (iUnitVerb, 120) i, azim, e1, e2, vz, &
&                             tzz, t1, t2, zTrans, zTranm, shear, stheta, &
&                                   pLon, pLat
		 END IF
	ELSE
		 IF (log_element_dynamics .AND. Verbose) &
		     &  WRITE (iUnitVerb, 121) i, azim, e1, e2, vz, &
&                             tzz, t1, t2, zTrans,       shear, stheta, &
&                                   pLon, pLat
	END IF
200  CONTINUE
IF (log_element_dynamics .AND. Verbose) WRITE (iUnitVerb, 210)
210  FORMAT ( &
& /' The figures above include vertical integrals of', &
&        ' normal stresses through the plate.  Compressive' &
& /' stresses are negative.  For convenience, normal stresses are', &
&        ' first corrected using a standard pressure curve' &
& /' P0(z), based on the structure of mid-ocean spreading', &
&        ' rises (see subprogram -SQUEEZ-).')

!   Fault element properties, also at midpoints:

IF (log_fault_dynamics .AND. Verbose) WRITE (iUnitVerb, 300)
300  FORMAT (/ / /' Fault element properties (at mid-points):'/ &
&  '                                       ', &
& '                                            ', &
&           '  Down-dip          Brittle/   Mantle         '/ &
&  '  Fault   Nodes#1,4    Horiz.  Azimuth', &
& ' Plunge    Total     Right    Perpen. Relative', &
&           ' integral     Peak  ductile brittle/ Is this '/ &
&  ' element  (N1 moves     slip        of', &
& '     of     slip   lateral shortning  vertical', &
&           ' of shear    shear    depth  ductile   fault '/ &
&  ' number   rel.to N4)    rate      slip', &
& '   slip     rate      rate      rate      rate', &
&           ' traction traction in crust    depth plastic?'/)
310  FORMAT (' ',I6,1X,I5,',',I5,1P,E9.2,0P,F10.2,F7.2, &
&             1P,E9.2,3E10.2,4E9.2,L3,I6)
311  FORMAT (' ',I6,1X,I5,',',I5,1P,E9.2,0P,F10.2,F7.2, &
&             1P,E9.2,3E10.2,3E9.2,' --------',L3,I6)
312  FORMAT (' ',I6,1X,I5,',',I5,1P,E9.2,0P,F10.2,F7.2, &
&             1P,E9.2,3E10.2,2E9.2,' --------',E9.2,L3,I6)
m = 4
DO 400 i = 1, nFl
	dip = 0.5D0 * (fDip(1, i) + fDip(2, i))
	j1 = nodeF(1, i)
	j2 = nodeF(2, i)
	j3 = nodeF(3, i)
	j4 = nodeF(4, i)
	jm = j1
	IF (jm > nRealN) jm = n1000 + (jm - nRealN)
	jb = j4
	IF (jb > nRealN) jb = n1000 + (jb - nRealN)
	du = v(1, j1) * fPFlt(1, 1, 1, 4, i) + v(2, j1) * fPFlt(2, 1, 1, 4, i) &
&         + v(1, j2) * fPFlt(1, 1, 2, 4, i) + v(2, j2) * fPFlt(2, 1, 2, 4, i) &
&         - v(1, j3) * fPFlt(1, 1, 2, 4, i) - v(2, j3) * fPFlt(2, 1, 2, 4, i) &
&         - v(1, j4) * fPFlt(1, 1, 1, 4, i) - v(2, j4) * fPFlt(2, 1, 1, 4, i)
	dv = v(1, j1) * fPFlt(1, 2, 1, 4, i) + v(2, j1) * fPFlt(2, 2, 1, 4, i) &
&         + v(1, j2) * fPFlt(1, 2, 2, 4, i) + v(2, j2) * fPFlt(2, 2, 2, 4, i) &
&         - v(1, j3) * fPFlt(1, 2, 2, 4, i) - v(2, j3) * fPFlt(2, 2, 2, 4, i) &
&         - v(1, j4) * fPFlt(1, 2, 1, 4, i) - v(2, j4) * fPFlt(2, 2, 1, 4, i)
	azimHS = 3.14159265358979D0 - ATan2F(dv, du)
	hors = SQRT(du**2 + dv**2)

!          "angle" is the fault strike, in radians cclkws from +sita (= +Theta; = South).

!CCCC       angle = 0.5 * (fArg(1, i) + fArg(2, i))
!CCCC       Line above was replaced due to cycle-shift problem

	angle = Chord(fArg(1, i), 0.5D0, fArg(2, i))
	unitx = COS(angle)
	unity = SIN(angle)
	crossx = -unity
	crossy = + unitx
	sinist = du * unitx + dv * unity
	close = du * crossx + dv * crossy
	IF (ABS(dip - 1.57079632679490D0) < wedge) THEN
		 vUpDip = 0.0D0
		 relV = 0.0D0
		 sNet = hors
		 plunge = 0.0D0
	ELSE
		 vUpDip = close / COS(dip)
		 relV = vUpDip * SIN(dip)
		 sNet = SQRT(hors**2 + relv**2)
		 plunge = -ASIN(relv / snet)
	END IF
	rLt = -sinist
	IF (ABS(dip - 1.57079632679490D0) < wedge) THEN
	   shear = fIMuDZ(4, i) * ABS(rlt)
	ELSE
	   shear = fIMuDZ(4, i) * snet / SIN(dip)
	END IF
	azimhs = azimhs * 57.2957795130823D0
	IF (azimhs >= 360.0D0) azimhs = azimhs - 360.0D0
	IF (azimhs <= -360.0D0) azimhs = azimhs + 360.0D0
	plunge = plunge * 57.2957795130823D0
	tlav = 0.5D0 * (tLNode(j1) + tLNode(j2))
	zmav = 0.5D0 * (zMNode(j1) + zMNode(j2))
	IF ((tlav > 0.0D0).AND. &
&          (zTranF(2, i) > (0.1D0 * oneKm))) THEN
		 fpMax = MAX(fPeakS(1, i), fPeakS(2, i))
		 zTranM = zmav + zTranF(2, i)
		 IF ((zTranF(1, i) / zmav) > 0.97D0) THEN
			  IF (log_fault_dynamics .AND. Verbose) WRITE (iUnitVerb, 312) i, jm, jb, hors, azimHS, plunge, &
&                     sNet, rLt, close, relV, shear, fpMax, &
&                                 zTranM, fSlips(i), i
		 ELSE
			  IF (log_fault_dynamics .AND. Verbose) WRITE (iUnitVerb, 310) i, jm, jb, hors, azimHS, plunge, &
&                     sNet, rLt, close, relv, shear, fpmax, &
&                     zTranF(1, i), zTranM, fSlips(i), i
		 END IF
	ELSE
		 IF (log_fault_dynamics .AND. Verbose) WRITE (iUnitVerb, 311) i, jm, jb, hors, azimhs, plunge, sNet, &
&                rLt, CLOSE, relV, shear, fPeakS(1, i), &
&                zTranF(1, i), fSlips(i), i
	END IF
400  CONTINUE
IF (log_fault_dynamics .AND. Verbose) WRITE (iUnitVerb, 401)
401  FORMAT(' ===========================================', &
&         '===========================================')

!----------------Begin writing to units iUnitT & iUnitQ---------------

! Single-plate torque-balance reports:

!      Zero out all torque components, prior to accumulating them:

DO 502 i = 1, 3
	DO 501 j = 1, nPlate
		 torqBS(i, j) = 0.0D0
		 torqCL(i, j) = 0.0D0
		 torqFS(i, j) = 0.0D0
		 torqLP(i, j) = 0.0D0
		 torqMD(i, j) = 0.0D0
		 torqSS(i, j) = 0.0D0
		 torqVB(i, j) = 0.0D0
501       CONTINUE
502  CONTINUE

!      Build torque components from info in "comp":

DO 510 i = 1, numNod

!           Subscript accounting:

	iPlate = whichP(i)
	ix = 2 * i - 1
	iy = 2 * i

!           Consistent nodal forces in (theta,phi) coordinates:

!           Basal Strength (2 components, and sum):
	fmdx = comp(5, ix)
	fmdy = comp(5, iy)
	fvbx = comp(6, ix)
	fvby = comp(6, iy)
	fbsx = comp(5, ix) + comp(6, ix)
	fbsy = comp(5, iy) + comp(6, iy)

!           Fault Strength
	ffsx = comp(3, ix)
	ffsy = comp(3, iy)

!           Lithostatic Pressure (sum of fault and basal):
	flpx = comp(2, ix) + comp(4, ix)
	flpy = comp(2, iy) + comp(4, iy)

!          (N.B. Sum of these consistent nodal forces
!                should be equal to comp(1).)

!           Uvec of the node:

	tTheta = xNode(i)
	pPhi = yNode(i)
	equat = SIN(ttheta)
	uvec(1) = equat * COS(pphi)
	uvec(2) = equat * SIN(pphi)
	uvec(3) = COS(ttheta)

!           Unit vectors at this site (NOT a pole):
	uPhi(1) = -uvec(2)
	uPhi(2) = uvec(1)
	uPhi(1) = uphi(1) / equat
	uPhi(2) = uphi(2) / equat
	uPhi(3) = 0.0D0
	tequat = uvec(3)
	uTheta(3) = -equat
	uTheta(1) = tequat * uvec(1) / equat
	uTheta(2) = tequat * uvec(2) / equat
	length = SQRT(utheta(1)**2 + utheta(2)**2 + utheta(3)**2)
	uTheta(1) = utheta(1) / length
	uTheta(2) = utheta(2) / length
	uTheta(3) = utheta(3) / length

!           Consistent nodal forces in (x,y,z):

	force_md(1) = fmdx * utheta(1) + fmdy * uphi(1)
	force_md(2) = fmdx * utheta(2) + fmdy * uphi(2)
	force_md(3) = fmdx * utheta(3) + fmdy * uphi(3)

	force_vb(1) = fvbx * utheta(1) + fvby * uphi(1)
	force_vb(2) = fvbx * utheta(2) + fvby * uphi(2)
	force_vb(3) = fvbx * utheta(3) + fvby * uphi(3)

	force_bs(1) = fbsx * utheta(1) + fbsy * uphi(1)
	force_bs(2) = fbsx * utheta(2) + fbsy * uphi(2)
	force_bs(3) = fbsx * utheta(3) + fbsy * uphi(3)

	force_fs(1) = ffsx * utheta(1) + ffsy * uphi(1)
	force_fs(2) = ffsx * utheta(2) + ffsy * uphi(2)
	force_fs(3) = ffsx * utheta(3) + ffsy * uphi(3)

	force_lp(1) = flpx * utheta(1) + flpy * uphi(1)
	force_lp(2) = flpx * utheta(2) + flpy * uphi(2)
	force_lp(3) = flpx * utheta(3) + flpy * uphi(3)

!           Nodal forces x moment arms:

	rvec(1) = radius * uvec(1)
	rvec(2) = radius * uvec(2)
	rvec(3) = radius * uvec(3)

	torqMD(1, iPlate) = torqMD(1, iPlate) + &
&                       rvec(2) * force_md(3) - rvec(3) * force_md(2)
	torqMD(2, iPlate) = torqMD(2, iPlate) + &
&                       rvec(3) * force_md(1) - rvec(1) * force_md(3)
	torqMD(3, iPlate) = torqMD(3, iPlate) + &
&                       rvec(1) * force_md(2) - rvec(2) * force_md(1)

	torqVB(1, iPlate) = torqVB(1, iPlate) + &
&                       rvec(2) * force_vb(3) - rvec(3) * force_vb(2)
	torqVB(2, iPlate) = torqVB(2, iPlate) + &
&                       rvec(3) * force_vb(1) - rvec(1) * force_vb(3)
	torqVB(3, iPlate) = torqVB(3, iPlate) + &
&                       rvec(1) * force_vb(2) - rvec(2) * force_vb(1)

	torqBS(1, iPlate) = torqBS(1, iPlate) + &
&                       rvec(2) * force_bs(3) - rvec(3) * force_bs(2)
	torqBS(2, iPlate) = torqBS(2, iPlate) + &
&                       rvec(3) * force_bs(1) - rvec(1) * force_bs(3)
	torqBS(3, iPlate) = torqBS(3, iPlate) + &
&                       rvec(1) * force_bs(2) - rvec(2) * force_bs(1)


	torqFS(1, iPlate) = torqFS(1, iPlate) + &
&                       rvec(2) * force_fs(3) - rvec(3) * force_fs(2)
	torqFS(2, iPlate) = torqFS(2, iPlate) + &
&                       rvec(3) * force_fs(1) - rvec(1) * force_fs(3)
	torqFS(3, iPlate) = torqFS(3, iPlate) + &
&                       rvec(1) * force_fs(2) - rvec(2) * force_fs(1)

	torqLP(1, iPlate) = torqLP(1, iPlate) + &
&                       rvec(2) * force_lp(3) - rvec(3) * force_lp(2)
	torqLP(2, iPlate) = torqLP(2, iPlate) + &
&                       rvec(3) * force_lp(1) - rvec(1) * force_lp(3)
	torqLP(3, iPlate) = torqLP(3, iPlate) + &
&                       rvec(1) * force_lp(2) - rvec(2) * force_lp(1)

510  CONTINUE

WRITE (iUnitQ, "(' ',A)") TRIM(title1)
WRITE (iUnitQ, "(' ',A)") TRIM(title2)
WRITE (iUnitQ, "(' ',A)") TRIM(title3)
WRITE (iUnitQ, * )
!Use default fFric for labelling purposes, to characterize this run:
t_fFric = LR_set_fFric(0)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
IF(Verbose) WRITE (iUnitVerb, 511) nPlate, t_fFric, tauMax(1), tauMax(2)
WRITE (iUnitQ, 511) nPlate, t_fFric, tauMax(1), tauMax(2)
511  FORMAT(/' Torque-balance reports for each of ', &
&           I3,' plates (fFric ',F5.3, &
&        ', tauMax ',ES8.1,'\',ES8.1,'):')

DO 600 n = 1, nPlate

!           Find rough center point for this plate #n,
!           defined by mean uvec of all nodes associated with it.
!          (If no nodes are associated with plate, skip to next.)

	nInSum = 0
	sumNod(1) = 0.0D0
	sumNod(2) = 0.0D0
	sumNod(3) = 0.0D0
	DO 515 i = 1, numNod
		 IF (whichP(i) == n) THEN
			  tTheta = xNode(i)
			  pPhi = yNode(i)
			  equat = SIN(tTheta)
			  uvec(1) = equat * COS(pPhi)
			  uvec(2) = equat * SIN(pPhi)
			  uvec(3) = COS(tTheta)
			  sumNod(1) = sumNod(1) + uvec(1)
			  sumNod(2) = sumNod(2) + uvec(2)
			  sumNod(3) = sumNod(3) + uvec(3)
			  nInSum = nInSum + 1
		 END IF
515       CONTINUE
	IF (nInSum == 0) GO TO 600

	IF(Verbose) WRITE (iUnitVerb, 516) n, names(n), nInSum
	WRITE(iUnitQ, 516) n, names(n), nInSum
516       FORMAT(/' plate #',I6,' =',A2,': ',I6,' nodes.')


!           ************************************************************
!                                *CRITICAL LOGIC*

!           N.B. torqBS = torqMD + torqVB (by definition),
!                but not computed here because each was computed
!                separately above.

	torqSS(1, n) = -torqBS(1, n) - torqLP(1, n)
	torqSS(2, n) = -torqBS(2, n) - torqLP(2, n)
	torqSS(3, n) = -torqBS(3, n) - torqLP(3, n)
!           so that these 3 will always add to zero (by definition).

	torqCL(1, n) = torqSS(1, n) - torqFS(1, n)
	torqCL(2, n) = torqSS(2, n) - torqFS(2, n)
	torqCL(3, n) = torqSS(3, n) - torqFS(3, n)
!           inferring continuum-link torque (if any) from residual.

	CALL Twist(area, detJ, fPSfer, & ! input
&                 iUnitT, n, nodes, nPlate, numEl, numNod, &
&                 radius, torqBS, whichP, xNode, yNode, &
&                 twistV)               ! output

!               (everything below is just reformatting and reporting)
!           ************************************************************

!           reformat Basal-Strength torque for this plate:

	tbsMag = SQRT(torqBS(1, n)**2 + torqBS(2, n)**2 + torqBS(3, n)**2)
	IF (tbsMag == 0.0D0) THEN
		tbsLon = 0.0D0
		tbsLat = 0.0D0
	ELSE
		dequat = SQRT(torqBS(1, n)**2 + torqBS(2, n)**2)
		IF (dequat == 0.0D0) THEN
			 IF (torqBS(3, n) > 0.0D0) THEN
				  tbsLat = 90.0D0
			 ELSE
				  tbsLat = -90.0D0
			 END IF
			 tbsLon = 0.0D0
		ELSE
			 tbsLat = 57.2957795130823D0 * ATAN2(torqBS(3, n), dequat)
			 tbsLon = 57.2957795130823D0 * ATAN2(torqBS(2, n), torqBS(1, n))
		END IF
	END IF

!           reformat Continuum-Link torque for this plate:

	tclMag = SQRT(torqCL(1, n)**2 + torqCL(2, n)**2 + torqCL(3, n)**2)
	IF (tclMag == 0.0D0) THEN
		tclLon = 0.0D0
		tclLat = 0.0D0
	ELSE
		dequat = SQRT(torqCL(1, n)**2 + torqCL(2, n)**2)
		IF (dequat == 0.0D0) THEN
			 IF (torqCL(3, n) > 0.0D0) THEN
				  tclLat = 90.0D0
			 ELSE
				  tclLat = -90.0D0
			 END IF
			 tclLon = 0.0D0
		ELSE
			 tclLat = 57.2957795130823D0 * ATAN2(torqCL(3, n), dequat)
			 tclLon = 57.2957795130823D0 * ATAN2(torqCL(2, n), torqCL(1, n))
		END IF
	END IF

!           reformat Fault-Strength torque for this plate:

	tfsMag = SQRT(torqFS(1, n)**2 + torqFS(2, n)**2 + torqFS(3, n)**2)
	IF (tfsMag == 0.0D0) THEN
		tfsLon = 0.0D0
		tfsLat = 0.0D0
	ELSE
		dequat = SQRT(torqFS(1, n)**2 + torqFS(2, n)**2)
		IF (dequat == 0.0D0) THEN
			 IF (torqFS(3, n) > 0.0D0) THEN
				  tfsLat = 90.0D0
			 ELSE
				  tfsLat = -90.0D0
			 END IF
			 tfsLon = 0.0D0
		ELSE
			 tfsLat = 57.2957795130823D0 * ATAN2(torqFS(3, n), dequat)
			 tfsLon = 57.2957795130823D0 * ATAN2(torqFS(2, n), torqFS(1, n))
		END IF
	END IF

!           reformat Lithostatic Pressure torque for this plate:

	tlpMag = SQRT(torqLP(1, n)**2 + torqLP(2, n)**2 + torqLP(3, n)**2)
	IF (tlpMag == 0.0D0) THEN
		tlpLon = 0.0D0
		tlpLat = 0.0D0
	ELSE
		dequat = SQRT(torqLP(1, n)**2 + torqLP(2, n)**2)
		IF (dequat == 0.0D0) THEN
			 IF (torqLP(3, n) > 0.0D0) THEN
				  tlpLat = 90.0D0
			 ELSE
				  tlpLat = -90.0D0
			 END IF
			 tlpLon = 0.0D0
		ELSE
			 tlpLat = 57.2957795130823D0 * ATAN2(torqLP(3, n), dequat)
			 tlpLon = 57.2957795130823D0 * ATAN2(torqLP(2, n), torqLP(1, n))
		END IF
	END IF

!           reformat Mantle-Drag torque for this plate:

	tmdMag = SQRT(torqMD(1, n)**2 + torqMD(2, n)**2 + torqMD(3, n)**2)
	IF (tmdMag == 0.0D0) THEN
		tmdLon = 0.0D0
		tmdLat = 0.0D0
	ELSE
		dequat = SQRT(torqMD(1, n)**2 + torqMD(2, n)**2)
		IF (dequat == 0.0D0) THEN
			 IF (torqMD(3, n) > 0.0D0) THEN
				  tmdLat = 90.0D0
			 ELSE
				  tmdLat = -90.0D0
			 END IF
			 tmdLon = 0.0D0
		ELSE
			 tmdLat = 57.2957795130823D0 * ATAN2(torqMD(3, n), dequat)
			 tmdLon = 57.2957795130823D0 * ATAN2(torqMD(2, n), torqMD(1, n))
		END IF
	END IF

!           reformat Side-Strength torque for this plate:

	tssMag = SQRT(torqSS(1, n)**2 + torqSS(2, n)**2 + torqSS(3, n)**2)
	IF (tssMag == 0.0D0) THEN
		tssLon = 0.0D0
		tssLat = 0.0D0
	ELSE
		dequat = SQRT(torqSS(1, n)**2 + torqSS(2, n)**2)
		IF (dequat == 0.0D0) THEN
			 IF (torqSS(3, n) > 0.0D0) THEN
				  tssLat = 90.0D0
			 ELSE
				  tssLat = -90.0D0
			 END IF
			 tssLon = 0.0D0
		ELSE
			 tssLat = 57.2957795130823D0 * ATAN2(torqSS(3, n), dequat)
			 tssLon = 57.2957795130823D0 * ATAN2(torqSS(2, n), torqSS(1, n))
		END IF
	END IF

!           reformat Velocity-Boundary-Condition torque for this plate:

	tvbMag = SQRT(torqVB(1, n)**2 + torqVB(2, n)**2 + torqVB(3, n)**2)
	IF (tvbMag == 0.0D0) THEN
		tvbLon = 0.0D0
		tvbLat = 0.0D0
	ELSE
		dequat = SQRT(torqVB(1, n)**2 + torqVB(2, n)**2)
		IF (dequat == 0.0D0) THEN
			 IF (torqVB(3, n) > 0.0D0) THEN
				  tvbLat = 90.0D0
			 ELSE
				  tvbLat = -90.0D0
			 END IF
			 tvbLon = 0.0D0
		ELSE
			 tvbLat = 57.2957795130823D0 * ATAN2(torqVB(3, n), dequat)
			 tvbLon = 57.2957795130823D0 * ATAN2(torqVB(2, n), torqVB(1, n))
		END IF
	END IF

!           reformat traction pole vector for this plate:

	twiMag = SQRT(twistV(1)**2 + twistV(2)**2 + twistV(3)**2)
	IF (twiMag == 0.0D0) THEN
		twiLon = 0.0D0
		twiLat = 0.0D0
	ELSE
		dequat = SQRT(twistV(1)**2 + twistV(2)**2)
		IF (dequat == 0.0D0) THEN
			 IF (twistV(3) > 0.0D0) THEN
				  twiLat = 90.0D0
			 ELSE
				  twiLat = -90.0D0
			 END IF
			 twiLon = 0.0D0
		ELSE
			 twiLat = 57.2957795130823D0 * ATAN2(twistV(3), dequat)
			 twiLon = 57.2957795130823D0 * ATAN2(twistV(2), twistV(1))
		END IF
	END IF

	IF(Verbose) WRITE (iUnitVerb, 520)
	WRITE(iUnitQ, 520)
520       FORMAT(/' Torques on plate bottoms:   X=0N,0E  Y=0N,90E' &
&             ,'     Z=90N Magnitude Longitude  Latitude' &
&             /' ------------------------- --------- ---------' &
&             ,' --------- --------- --------- ---------')
	IF(Verbose) WRITE (iUnitVerb, 521)(torqMD(i, n), i = 1, 3), tmdmag, tmdlon, tmdlat
	WRITE (iUnitQ, 521)(torqMD(i, n), i = 1, 3), tmdmag, tmdlon, tmdlat
521       FORMAT(' Mantle-Drag:             ',3ES10.2,ES10.3,2F10.2)
	IF(Verbose) WRITE (iUnitVerb, 522)(torqVB(i, n), i = 1, 3), tvbmag, tvblon, tvblat
	WRITE (iUnitQ, 522)(torqVB(i, n), i = 1, 3), tvbmag, tvblon, tvblat
522       FORMAT(' Velocity-Boundary-C.''s   ',3ES10.2,ES10.3,2F10.2)
	IF(Verbose) WRITE (iUnitVerb, 523)
	WRITE(iUnitQ, 523)
523       FORMAT(' ---------------------------------------------' &
&            ,'----------------------------------------')
	IF(Verbose) WRITE (iUnitVerb, 524)(torqBS(i, n), i = 1, 3), tbsmag, tbslon, tbslat
	WRITE (iUnitQ, 524)(torqBS(i, n), i = 1, 3), tbsmag, tbslon, tbslat
524       FORMAT(' Basal-Strength:          ',3ES10.2,ES10.3,2F10.2)

	IF(Verbose) WRITE (iUnitVerb, 530)
	WRITE(iUnitQ, 530)
530       FORMAT(/' Torques on plate sides:     X=0N,0E  Y=0N,90E' &
&             ,'     Z=90N Magnitude Longitude  Latitude' &
&             /' ------------------------- --------- ---------' &
&             ,' --------- --------- --------- ---------')
	IF(Verbose) WRITE (iUnitVerb, 531)(torqFS(i, n), i = 1, 3), tfsmag, tfslon, tfslat
	WRITE (iUnitQ, 531)(torqFS(i, n), i = 1, 3), tfsmag, tfslon, tfslat
531       FORMAT(' Fault-Strength:          ',3ES10.2,ES10.3,2F10.2)
	IF(Verbose) WRITE (iUnitVerb, 532)(torqCL(i, n), i = 1, 3), tclmag, tcllon, tcllat
	WRITE (iUnitQ, 532)(torqCL(i, n), i = 1, 3), tclmag, tcllon, tcllat
532       FORMAT(' Continuum-Link [PLUG]:   ',3ES10.2,ES10.3,2F10.2)
	IF(Verbose) WRITE (iUnitVerb, 533)
	WRITE(iUnitQ, 533)
533       FORMAT(' ---------------------------------------------' &
&            ,'----------------------------------------')
	IF(Verbose) WRITE (iUnitVerb, 534)(torqSS(i, n), i = 1, 3), tssmag, tsslon, tsslat
	WRITE (iUnitQ, 534)(torqSS(i, n), i = 1, 3), tssmag, tsslon, tsslat
534       FORMAT(' Side-Strength:           ',3ES10.2,ES10.3,2F10.2)

	IF(Verbose) WRITE (iUnitVerb, 540)
	WRITE(iUnitQ, 540)
540       FORMAT(/' Kind of torque:             X=0N,0E  Y=0N,90E' &
&             ,'     Z=90N Magnitude Longitude  Latitude' &
&             /' ------------------------- --------- ---------' &
&             ,' --------- --------- --------- ---------')
	IF(Verbose) WRITE (iUnitVerb, 541)(torqLP(i, n), i = 1, 3), tlpmag, tlpLon, tlpLat
	WRITE (iUnitQ, 541)(torqLP(i, n), i = 1, 3), tlpmag, tlpLon, tlpLat
541       FORMAT(' Lithostatic-Pressure:    ',3ES10.2,ES10.3,2F10.2)
	IF(Verbose) WRITE (iUnitVerb, 542)(torqSS(i, n), i = 1, 3), tssmag, tsslon, tsslat
	WRITE (iUnitQ, 542)(torqSS(i, n), i = 1, 3), tssmag, tsslon, tsslat
542       FORMAT(' Side-Strength:           ',3ES10.2,ES10.3,2F10.2)
	IF(Verbose) WRITE (iUnitVerb, 543)(torqBS(i, n), i = 1, 3), tbsmag, tbslon, tbslat
	WRITE (iUnitQ, 543)(torqBS(i, n), i = 1, 3), tbsmag, tbslon, tbslat
543       FORMAT(' Basal-Strength:          ',3ES10.2,ES10.3,2F10.2)

	IF(Verbose) WRITE (iUnitVerb, 550)
	WRITE (iUnitQ, 550)
550       FORMAT(/' Traction pole vector:       X=0N,0E  Y=0N,90E' &
&             ,'     Z=90N Magnitude Longitude  Latitude' &
&             /' ------------------------- --------- ---------' &
&             ,' --------- --------- --------- ---------')
	IF(Verbose) WRITE (iUnitVerb, 551)(twistV(i), i = 1, 3), twimag, twilon, twilat
	WRITE (iUnitQ, 551)(twistV(i), i = 1, 3), twimag, twilon, twilat
551       FORMAT(' Basal-strength:          ',3ES10.2,ES10.3,2F10.2)

!        Find viewpoint orthogonal to all 3 (BS,LP,SS) torque vectors:
	dvec(1) = torqSS(2, n) * torqBS(3, n) - torqSS(3, n) * torqBS(2, n)
	dvec(2) = torqSS(3, n) * torqBS(1, n) - torqSS(1, n) * torqBS(3, n)
	dvec(3) = torqSS(1, n) * torqBS(2, n) - torqSS(2, n) * torqBS(1, n)

!           Check that viewpoint is on same side of planet as plate:

	ddot = dvec(1) * sumNod(1) + dvec(2) * sumNod(2) + dvec(3) * sumNod(3)
	IF (ddot < 0.0D0) THEN
		 dvec(1) = -dvec(1)
		 dvec(2) = -dvec(2)
		 dvec(3) = -dvec(3)
	END IF
	length = SQRT(dvec(1)**2 + dvec(2)**2 + dvec(3)**2)
	IF (length == 0.0D0) THEN
		longit = 0.0D0
		latitu = 0.0D0
		uvec(1) = 1.0D0
		uvec(2) = 0.0D0
		uvec(3) = 0.0D0
	ELSE
		uvec(1) = dvec(1) / length
		uvec(2) = dvec(2) / length
		uvec(3) = dvec(3) / length
		dequat = SQRT(uvec(1)**2 + uvec(2)**2)
		IF (dequat == 0.0D0) THEN
			 IF (uvec(3) > 0.0D0) THEN
				  latitu = 90.0D0
			 ELSE
				  latitu = -90.0D0
			 END IF
			 longit = 0.0D0
		ELSE
			 latitu = 57.2957795130823D0 * ATAN2(uvec(3), dequat)
			 longit = 57.2957795130823D0 * ATAN2(uvec(2), uvec(1))
		END IF
	END IF
	IF(Verbose) WRITE (iUnitVerb, 560) longit, latitu
	WRITE (iUnitQ, 560) longit, latitu
560       FORMAT(/' Suggested viewpoint for orthographic projection' &
&            ,' of torques on this plate is: (',F7.2,'E,',F6.2,'N)' &
&             /' from which direction all 3 torque vectors will be' &
&             ,' in the plane of the figure.')

	IF(Verbose) WRITE (iUnitVerb, 570)
	WRITE (iUnitQ, 570)
570       FORMAT(/' Equivalent horizontal forces at this point:' &
&             /' Kind of force:              X=0N,0E  Y=0N,90E' &
&             ,'     Z=90N Magnitude   Azimuth' &
&             /' ------------------------- --------- ---------' &
&             ,' --------- --------- ---------')

	equat = SQRT(uvec(1)**2 + uvec(2)**2)
	uPhi(1) = -uvec(2)
	uPhi(2) = uvec(1)
	uPhi(1) = uPhi(1) / equat
	uPhi(2) = uPhi(2) / equat
	uPhi(3) = 0.0D0
	tequat = uvec(3)
	uTheta(3) = -equat
	uTheta(1) = tequat * uvec(1) / equat
	uTheta(2) = tequat * uvec(2) / equat
	length = SQRT(utheta(1)**2 + utheta(2)**2 + utheta(3)**2)
	uTheta(1) = uTheta(1) / length
	uTheta(2) = uTheta(2) / length
	uTheta(3) = uTheta(3) / length

!           Lithostatic pressure force:

	force_lp(1) = (torqLP(2, n) * uvec(3) - torqLP(3, n) * uvec(2)) / radius
	force_lp(2) = (torqLP(3, n) * uvec(1) - torqLP(1, n) * uvec(3)) / radius
	force_lp(3) = (torqLP(1, n) * uvec(2) - torqLP(2, n) * uvec(1)) / radius
	length = SQRT(torqLP(1, n)**2 + torqLP(2, n)**2 + torqLP(3, n)**2) / &
&             radius
	dot_s = force_lp(1) * utheta(1) + force_lp(2) * utheta(2) + &
&            force_lp(3) * utheta(3)
	dot_e = force_lp(1) * uphi(1) + force_lp(2) * uphi(2) + &
&            force_lp(3) * uphi(3)
	dot_n = -dot_s
	azimut = 57.2957795130823D0 * ATAN2(dot_e, dot_n)
	IF (azimut < 0.0D0) azimut = azimut + 360.0D0
	IF(Verbose) WRITE (iUnitVerb, 571)(force_lp(i), i = 1, 3), length, azimut
	WRITE (iUnitQ, 571)(force_lp(i), i = 1, 3), length, azimut
571       FORMAT(' Lithostatic pressure:    ',3ES10.2,ES10.3,F10.1)

!           Side-strength force:

	force_ss(1) = (torqSS(2, n) * uvec(3) - torqSS(3, n) * uvec(2)) / radius
	force_ss(2) = (torqSS(3, n) * uvec(1) - torqSS(1, n) * uvec(3)) / radius
	force_ss(3) = (torqSS(1, n) * uvec(2) - torqSS(2, n) * uvec(1)) / radius
	length = SQRT(torqSS(1, n)**2 + torqSS(2, n)**2 + torqSS(3, n)**2) / &
&             radius
	dot_S = force_ss(1) * utheta(1) + force_ss(2) * utheta(2) + &
&            force_ss(3) * utheta(3)
	dot_E = force_ss(1) * uphi(1) + force_ss(2) * uphi(2) + &
&            force_ss(3) * uphi(3)
	dot_N = -dot_s
	azimut = 57.2957795130823D0 * ATAN2(dot_E, dot_N)
	IF (azimut < 0.0D0) azimut = azimut + 360.0D0
	IF(Verbose) WRITE (iUnitVerb, 572)(force_ss(i), i = 1, 3), length, azimut
	WRITE (iUnitQ, 572)(force_ss(i), i = 1, 3), length, azimut
572       FORMAT(' Side-strength:           ',3ES10.2,ES10.3,F10.1)

!           Basal-strength force:

	force_bs(1) = (torqBS(2, n) * uvec(3) - torqBS(3, n) * uvec(2)) / radius
	force_bs(2) = (torqBS(3, n) * uvec(1) - torqBS(1, n) * uvec(3)) / radius
	force_bs(3) = (torqBS(1, n) * uvec(2) - torqBS(2, n) * uvec(1)) / radius
	length = SQRT(torqBS(1, n)**2 + torqBS(2, n)**2 + torqBS(3, n)**2) / &
&             radius
	dot_S = force_bs(1) * utheta(1) + force_bs(2) * utheta(2) + &
&            force_bs(3) * utheta(3)
	dot_E = force_bs(1) * uphi(1) + force_bs(2) * uphi(2) + &
&            force_bs(3) * uphi(3)
	dot_N = -dot_s
	azimut = 57.2957795130823D0 * ATAN2(dot_E, dot_N)
	IF (azimut < 0.0D0) azimut = azimut + 360.0D0
	IF(Verbose) WRITE (iUnitVerb, 573)(force_bs(i), i = 1, 3), length, azimut
	WRITE (iUnitQ, 573)(force_bs(i), i = 1, 3), length, azimut
573       FORMAT(' Basal-strength:          ',3ES10.2,ES10.3,F10.1)

	dequat = SQRT(sumNod(1)**2 + sumNod(2)**2)
	latitu = 57.2957795130823D0 * ATAN2(sumNod(3), dequat)
	longit = 57.2957795130823D0 * ATAN2(sumNod(2), sumNod(1))
	IF(Verbose) WRITE (iUnitVerb, 580) longit, latitu
	WRITE (iUnitQ, 580) longit, latitu
580       FORMAT(/' and this cluster of forces should be connected' &
&             ,' by a leader line' &
&             /' to the plate center at approximately:' &
&             ,' (',F7.2,'E, ',F6.2,'N).')

	IF(Verbose) WRITE (iUnitVerb, 401)
	WRITE (iUnitQ, 401)
600  CONTINUE
CLOSE (iUnitQ)
RETURN
END SUBROUTINE Result

SUBROUTINE Rotor (mxDOF, nDOF, nLB, node, nUB, theta, & ! input
&                   force, stiff)                         ! modify

!   Operate on two adjacent row equations of the linear system
!   (coefficient matrix "stiff" and right-side vector "force")
!   which represent the balance of forces on one node in the
!   x and y directions, respectively.
!   Rotate these equations to a new coordinate system (alpha, beta)
!   where alpha is theta radians counterclockwise from x, and
!   beta is theta radians counterclockwise from y.

!   Note: This transformation has ***no effect*** on the definitions
!   of the unknown velocities, which remain in the (x, y) system.

!   The rows operated on are #(2*node-1) and #(2*node).
!   After rotation, the alpha-equation will replace the x-equation,
!   and the beta-equation will replace the y-equation.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: mxDOF, nDOF, nLB, node, nUB                                     ! input
REAL*8, INTENT(IN) :: theta                                                            ! input
DOUBLE PRECISION, INTENT(INOUT) :: force, stiff                                        ! modify
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
INTEGER iarow, ibrow, iq, ixrow, iyrow, j1, j2, jColum
DOUBLE PRECISION cosT, dTheta, sinT, xtemp, ytemp
DIMENSION force(mxDOF, 1), stiff(nKRows, nRank)

dTheta = theta
cosT = COS(dTheta)
sinT = SIN(dTheta)
ixrow = 2 * node - 1
iyrow = 2 * node
iarow = ixrow
ibrow = iyrow
xtemp = force(ixrow, 1)
ytemp = force(iyrow, 1)
force(iarow, 1) = cosT * xtemp + sinT * ytemp
force(ibrow, 1) = cosT * ytemp - sinT * xtemp
j1 = MAX(iyrow - nLB, 1)
j2 = MIN(ixrow + nUB, nDOF)
DO 10 jColum = j1, j2
   !matrix element(ixrow, jColum):
	iq = iDiagonal + ixrow - jColum
	xtemp = stiff(iq, jColum)
   !matrix element(iyrow, jColum):
	iq = iDiagonal + iyrow - jColum
	ytemp = stiff(iq, jColum)
   !matrix element(ixrow, jColum):
	iq = iDiagonal + ixrow - jColum
	stiff(iq, jColum) = cosT * xtemp + sinT * ytemp
   !matrix element(iyrow, jColum):
	iq = iDiagonal + iyrow - jColum
	stiff(iq, jColum) = cosT * ytemp - sinT * xtemp
10  CONTINUE
RETURN
END SUBROUTINE Rotor

SUBROUTINE Sander (fDip, iCond, iUnitT, & ! input
&                    log_strike_adjustments, &
&                    mxBn, mxFEl, mxNode, nCond, nFl, &
&                    nodCon, nodeF, vBCArg, vBCMag, &
&                    wedge, xNode, yNode, &
&                    fArg)                  ! modify

!     "Rounds the angular corners" of any model edges which are
!      multi-element strike-slip fault systems, by averaging the
!      arguments at matched ends of the adjacent s-s fault elements.

!      This is only done where boundary conditions for the external
!      nodes are identical, creating one rigid plate outside the
!      model domain.

!      This correction is necessary to prevent two artifacts:
!     -Extremely large equal-but-opposite boundary force
!      vectors plotted at the same location (for the two
!      external nodes that are co-located).
!     -Artificial resistance to strike-slip, since the
!      resistance added by mismatched arguments is proportional
!      to fMuMax, but not dependent on fault or plate rheology!!!

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: fDip                                                             ! input
INTEGER, INTENT(IN) :: iCond, iUnitT                                                   ! input
LOGICAL, INTENT(IN) :: log_strike_adjustments                                          ! input
INTEGER, INTENT(IN) :: mxBn, mxFEl, mxNode, nCond, nFl, nodCon, nodeF                  ! input
REAL*8, INTENT(IN) :: vBCArg, vBCMag, wedge, xNode, yNode                              ! input
REAL*8, INTENT(INOUT) :: fArg                                                          ! modify
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER b1, b2, end1, end2, f1, f2, i, n1, n2, on1, on2
REAL*8 argMid, deg1, deg2, degM, dELon, dNLat
LOGICAL didOne
DIMENSION iCond(mxBn), nodCon(mxBn), nodeF(4, mxFEl)
DIMENSION fArg(2, mxFEl), fDip(2, mxFEl), vBCArg(mxBn), vBCMag(mxBn), &
&           xNode(mxNode), yNode(mxNode)

IF (log_strike_adjustments) WRITE(iUnitT, 1)
1  FORMAT(/ /' The following pairs of model-bounding strike-slip' &
&          /' fault elements had their azimuths averaged at the' &
&          /' connection point for purposes of computing the' &
&          /' constraint on the directino of strike-slip:' &
&        / /' Fault#1   Fault#2    Node#A    Node#B   ', &
&            '  Latitude Longitude    Azim#1    Azim#2   Azimuth' &
&          /' ----------------------------------------', &
&            '--------------------------------------------------')
didOne = .FALSE.
!      loop on all boundary nodes (referring backwards for neighbors)
b1 = nCond
DO 1000 b2 = 1, nCond
	n1 = nodCon(b1)
	n2 = nodCon(b2)

!           consider only if 2 consecutive boundary nodes are colocated
	IF ((xNode(n1) == xNode(n2)).AND. &
&          (yNode(n1) == yNode(n2))) THEN

!                consider only if both boundary nodes are type-2, 4, 5:
		 IF (((iCond(b1) == 2).AND.(iCond(b2) == 2)).OR. &
&               ((iCond(b1) == 4).AND.(iCond(b2) == 4)).OR. &
&               ((iCond(b1) == 5).AND.(iCond(b2) == 5))) THEN

!                consider only if both type-2/4/5 BC's are same velocity
			  IF ((vBCArg(b1) == vBCArg(b2)).AND. &
&                    (vBCMag(b1) == vBCMag(b2))) THEN

! <<<<<<<<<<<<<<<<<<<<<<<< shift <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!      Find fault element F1 containing node N1 (and remember the end
!      end1 and opposing inner node on1):
f1 = 0
DO 110 i = 1, nFl
	IF (n1 == nodeF(1, i)) THEN
		 f1 = i
		 end1 = 1
		 on1 = nodeF(4, i)
		 GO TO 111
	END IF
	IF (n1 == nodeF(4, i)) THEN
		 f1 = i
		 end1 = 1
		 on1 = nodeF(1, i)
		 GO TO 111
	END IF
	IF (n1 == nodeF(2, i)) THEN
		 f1 = i
		 end1 = 2
		 on1 = nodeF(3, i)
		 GO TO 111
	END IF
	IF (n1 == nodeF(3, i)) THEN
		 f1 = i
		 end1 = 2
		 on1 = nodeF(2, i)
		 GO TO 111
	END IF
110  CONTINUE
!      Find fault element f2 containing node n2 (and remember the end
!      end2 and opposing inner node on2):
111  f2 = 0
DO 120 i = 1, nFl
	IF (n2 == nodeF(1, i)) THEN
		 f2 = i
		 end2 = 1
		 on2 = nodeF(4, i)
		 GO TO 121
	END IF
	IF (n2 == nodeF(4, i)) THEN
		 f2 = i
		 end2 = 1
		 on2 = nodeF(1, i)
		 GO TO 121
	END IF
	IF (n2 == nodeF(2, i)) THEN
		 f2 = i
		 end2 = 2
		 on2 = nodeF(3, i)
		 GO TO 121
	END IF
	IF (n2 == nodeF(3, i)) THEN
		 f2 = i
		 end2 = 2
		 on2 = nodeF(2, i)
		 GO TO 121
	END IF
120  CONTINUE
!      Consider only if 2 distinct faults were found:
121  IF ((f1 > 0).AND.(f2 > 0).AND.(f1 /= f2)) THEN

!           Consider only if both faults are vertical:
	IF ((ABS(fDip(end1, f1) - 1.57079632679490D0) <= wedge).AND. &
&          (ABS(fDip(end2, f2) - 1.57079632679490D0) <= wedge)) THEN

!                Consider only if opposite/inner nodes are the same
!               (no internal fault creates a triple-junction):
		 IF (on1 == on2) THEN

			  argMid = Chord(fArg(end1, f1), 0.50D0, fArg(end2, f2))
			  deg1 = 180.0D0 - 57.2957795130823D0 * fArg(end1, f1)
			  deg2 = 180.0D0 - 57.2957795130823D0 * fArg(end2, f2)
			  degM = 180.0D0 - 57.2957795130823D0 * argMid
!                    =================== modify! =====================
			  fArg(end1, f1) = argMid
			  fArg(end2, f2) = argMid
			  didOne = .TRUE.
!                    =================== modify! =====================
!                     Write a line for the output table:
			  dELon = 57.2957795130823D0 * yNode(n1)
			  dNLat = 90.0D0 - 57.2957795130823D0 * xNode(n1)
			  IF (log_strike_adjustments .AND. Verbose) WRITE (iUnitVerb, 900) &
&                                   f1, f2, n1, n2, &
&                                   dNLat, dELon, deg1, deg2, degm
900                 FORMAT(' ',I7,3X,I7,3X, &
&                           I7,3X,I7,3X, &
&                           2X,F8.3,1X,F9.3, &
&                           4X,F6.1,4X,F6.1,4X,F6.1)

		 END IF
!                ^end of test that opposite/inner nodes are the same

	END IF
!           ^end of test that both faults are vertical

END IF
!      ^end of test that 2 distinct faults were found

! >>>>>>>>>>>>>>>>>>>>>>>> shift >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

			  END IF
!                     ^end of test that 2 type-2 VBC's are same vector

		 END IF
!                ^end of test that both nodes are type-2

	END IF
!           ^end of test that two consecutive boundary nodes are colocat

!           prepare to loop: current lead node becomes new following nod
	b1 = b2
1000  CONTINUE
!      ^end loop on all boundary nodes (referring backwards to neighbor

IF (.NOT.didone) THEN
	IF(Verbose) WRITE (iUnitVerb, 1001)
1001       FORMAT(' (No fault pairs were found which needed ' &
&              ,'this correction.)')
END IF
END SUBROUTINE Sander

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
RETURN
END SUBROUTINE SNodal

SUBROUTINE Solver (iUnitT, &  ! input
				& ab, b, &   ! modify (coefficient matrix and forcing vector)
				& ipiv)      ! work

!   Sets up for CALL to the MKL library routine which actually
!   solves the linear system

!                      |ab| |x| = |b|.

!   The left-hand coefficient matrix |ab| ("a;banded") is destroyed.
!   The answer vector |x| is written over the forcing vector |b|.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iUnitT                                           ! input
DOUBLE PRECISION, INTENT(INOUT) :: ab, b                                ! modify
INTEGER :: ipiv                                                         ! work
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
INTEGER info, kl, ku, ldab, ldb, n, nrhs
DIMENSION ab(nKRows, nRank), b(nRank, 1), ipiv(nRank)

!   ----- Name conversions -----------------------------------
!      Coefficient matrix: "stiff" or "k" in -Shells-; here called "a".
!      Right-hand forcing vector: "f" in -Shells-; here called "b".
!      Note that last argument below is the solution vector;
!      here it is overwritten onto "b" to save storage.

n = nRank
kl = nCodiagonals
ku = nCodiagonals
nrhs = 1
ldab = nKRows
ldb = nRank

!-----------------------------------------------------------
CALL dgbsv ( n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info )
!-----------------------------------------------------------

! The manual page for "?gbsv" of MKL/LAPACK is here: https://software.intel.com/en-us/node/468882

! Be CAREFUL because some pages in the MKL/LAPACK manual give INCORRECT descriptions (and
!    illustrations, and examples!) of the band-storage scheme, describing only kl+1+ku rows in ab,
!    and the diagonal row located at row kl+1.

! Correct documentation can be found here: http://www.netlib.no/netlib/lapack/double/dgbsv.f
! The CORRECT storage scheme has 2*kl+1+ku rows, with the diagonal row located at row 2*kl+1.
! I have checked (2016.07.08) that solutions under this correct scheme are essentially
! identical to old solutions, using an IMSL solver, in Shells_v3.9.

IF (info /= 0) THEN
  write(ErrorMsg,'(A,I12)') "ERROR: dgbsv (of the MKL library) reports results info = ",info
  call FatalError(ErrorMsg,ThID)
END IF

RETURN
END SUBROUTINE Solver

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
INTEGER, INTENT(IN) :: iUnitT                                                          ! input
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
						WRITE(iUnitT, 476) nInSum, &
&                               (list(n), n = 1, nInSum)
476                           FORMAT(/ &
&                           ' AVERAGING TOGETHER THE POSITIONS OF', &
&                               ' THESE ',I6,' NODES:',(/' ',12I6))
						IF(Verbose) WRITE (iUnitVerb, 477) rmax
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
			  WRITE(iUnitT, 605) m, i
605                 FORMAT(/' EXCESSIVELY DISTORTED ELEMENT LEADS TO ' &
&                     ,'NEGATIVE AREA AT POINT ',I1,' IN ELEMENT ', &
&                       I5)
			  WRITE(iUnitT, 606) area(i), detJ(m, i)
606                 FORMAT('AREA = ',1P,E12.4,'   detJ: ',0P,F12.6)
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
       	  write(ErrorMsg,'(A/,A,I6/,A/,A)') "BAD GRID TOPOLOGY: WHILE TRACING PERIMETER," , &
&          "COULD NOT FIND ANY WAY TO CONTINUE FROM NODE ",node, &
&          "EITHER THROUGH SHARED BOUNDARY ELEMENTS, OR" , &
&          "THROUGH OTHER BOUNDARY NODES SHARING THE SAME POSITION."
	      call FatalError(ErrorMsg,ThID)
867       nDone = nDone + 1
!           Success by location method: I is the next boundary node
	IF (nDone <= nCond) nodCon(nDone) = i
	nLeft = nLeft - 1
	IF (nLeft > 0) GO TO 840
!      End of indefinate loop which traces around perimeter.
870  IF (.NOT.skipBC .AND. Verbose) THEN
	   WRITE(iUnitVerb, 880)
880    FORMAT(/ /' Here follows a list, in consecutive order,'/ &
&                ' of the nodes which define the perimeter'/ &
&                ' of the model; these nodes require boundary', &
&                ' conditions:'/'    BC#  Node          ', &
&                '  Latitude Longitude')
	   DO 890 i = 1, nCond
		 n = nodCon(i)
		 theLon = yNode(n) * 57.2957795130823D0
		 theLat = 90.0D0 - xNode(n) * 57.2957795130823D0
		 WRITE(iUnitVerb, 882) i, n, theLat, theLon
882      FORMAT(' ',2I6,10X,2F10.3)
890    CONTINUE
	 n = nodCon(1)
	 IF(Verbose) WRITE (iUnitVerb, 892) n
892  FORMAT(' (Note: Initial node ',I6,' completes the loop,', &
&           ' but is not listed again.)')
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
		 IF(Verbose) WRITE (iUnitVerb, 905) i, dip1, tag1, dip2, tag2
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
				   IF(Verbose) WRITE (iUnitVerb, 912) i, dip1, dip2
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
				   IF(Verbose) WRITE (iUnitVerb, 914) i, dip1, dip2
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

IF (log_strike_adjustments) WRITE(iUnitT, 1001)
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
1601      IF (found.AND.(number > i)) THEN
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
			  IF (log_strike_adjustments .AND. Verbose) WRITE (iUnitVerb, 1610) &
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

IF (Verbose) WRITE (iUnitVerb, 9999)
9999  FORMAT (' --------------------------------------------------', &
&          '-----------------------------')
RETURN
END SUBROUTINE Square

SUBROUTINE Squeez (alphaT, density_anomaly_kgpm3, elevat, & ! input
&                    geoth1, geoth2, geoth3, geoth4, &
&                    geoth5, geoth6, geoth7, geoth8, &
&                    gMean, &
&                    iUnitT, oneKm, rhoAst, rhoBar, rhoH2O, &
&                    temLim, zM, zStop, &
&                    tauzz, sigzzb)                           ! output

!   Calculates tauzz, the vertical integral through the plate
!      of the vertical stress anomaly, which is defined
!      relative to a column of mantle with asthenosphere density
!      with a 5 km crust and a 2.7 km ocean on top, like a mid-ocean
!      spreading rise of high spreading velocity.
!      The integral is from either the land surface or the
!      sea surface, down to a depth of zStop below the top of
!      the crust.
!      If zStop exceeds Moho depth zM, then properties of the mantle
!      will be used in the lower part of the integral.
!   Also returns sigzzb, the standardized vertical stress anomaly
!      at depth zStop below the solid rock surface.
!   Note: This version is different from the version found in the -Laramy-
!      program package.  First, it acts on only a single point.
!      Second, it infers sub-plate normal-stress anomalies from
!      the given topography, instead of from model structure.
!      Finally, it was modified (in 2005, for Earth5) to accept
!      the additional input parameter density_anomaly_kgpm3,
!      which is a density anomaly of chemical origin (applying to
!      both crust and mantle lithosphere) in addition to the
!      crust/mantle density difference, and density variations
!      of thermal origin.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alphaT, density_anomaly_kgpm3, elevat, geoth1, geoth2, geoth3, & ! input
  & geoth4, geoth5, geoth6, geoth7, geoth8, gMean                                     ! input
INTEGER, INTENT(IN) :: iUnitT                                                          ! input
REAL*8, INTENT(IN) :: oneKm, rhoAst, rhoBar, rhoH2O, temLim, zM, zStop                 ! input
REAL*8, INTENT(OUT) :: tauzz, sigzzb                                                   ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8 TempC, TempM, h ! statement functions
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, PARAMETER :: ndRef = 300
INTEGER i, j, lastDR, layer1, layer2, n1, n2, nStep
REAL*8 dense, dense1, dense2, dRef, frac, frac1, frac2, &
	& oldpr, oldszz, pr, pRef, resid, rhotop, sigzz, t, z, zBase, zTop
LOGICAL :: called = .FALSE.
!   Internal arrays:
DIMENSION dRef(ndRef), pRef(0:ndRef)
!   Argument arrays:
DIMENSION alphaT(2), rhoBar(2), temLim(2)


!   Statement functions:
TempC(h) = MIN(temLim(1), geoth1 + geoth2 * h + geoth3 * h**2 + geoth4 * h**3)
TempM(h) = MIN(temLim(2), geoth5 + geoth6 * h + geoth7 * h**2 + geoth8 * h**3)

!   Create reference temperature & density profiles to depth of ndRef kilometers:

IF (.NOT.called) THEN
	rhotop = rhoBar(1) * (1.0D0 - alphaT(1) * geoth1)
	dRef(1) = rhoH2O
	dRef(2) = rhoH2O
	dRef(3) = 0.70D0 * rhoH2O + 0.30D0 * rhotop
	dRef(4) = rhotop
	dRef(5) = rhotop
	dRef(6) = rhotop
	dRef(7) = rhotop
	dRef(8) = 0.70D0 * rhotop + 0.30D0 * rhoAst
	DO 50 j = 9, ndRef
		 dRef(j) = rhoAst
50       CONTINUE
	pRef(0) = 0.0D0
	DO 100 i = 1, ndRef
		 pRef(i) = pRef(i - 1) + dRef(i) * gMean * oneKm
100       CONTINUE
END IF

!   Routine processing (in every CALL):

IF (elevat > 0.0D0) THEN
!        Land:
	zTop = -elevat
	zBase = zStop - elevat
	dense1 = rhoBar(1) * (1.0D0 - geoth1 * alphaT(1)) + &
&             density_anomaly_kgpm3
	h = 0.0D0
	layer1 = 1
ELSE
!         Ocean:
	zTop = 0.0D0
	zBase = zStop + (-elevat)
	dense1 = rhoH2O
	h = elevat
	layer1 = 0
END IF
lastDR = zBase / oneKm
IF (zBase > (oneKm * lastDR)) lastDR = lastDR + 1
IF (lastDR > ndRef) THEN
	write(ErrorMsg,'(A,I10)') "IN SUBPROGRAM SQUEEZ, PARAMETER ndRef MUST BE INCREASED TO AT LEAST ",lastDR
	call FatalError(ErrorMsg,ThID)
END IF
nStep = (zBase - zTop) / oneKm
oldszz = 0.0D0
oldpr = 0.0D0
sigzz = 0.0D0
tauzz = 0.0D0
z = zTop
DO 200 i = 1, nStep
	z = z + oneKm
	h = h + oneKm
	IF (h > 0.0D0) THEN
		 IF (h <= zM) THEN
			  t = TempC(h)
			  dense2 = rhoBar(1) * (1.0D0 - t * alphaT(1)) + &
&                       density_anomaly_kgpm3
			  layer2 = 1
		 ELSE
			  t = TempM(h - zM)
			  dense2 = rhoBar(2) * (1.0D0 - t * alphaT(2)) + &
&                       density_anomaly_kgpm3
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
		 frac1 = 0.50D0
		 frac2 = 0.50D0
	END IF
	dense = frac1 * dense1 + frac2 * dense2
	IF (z > 0.0D0) THEN
		 n1 = z / oneKm
		 n2 = n1 + 1
		 frac = z / oneKm - n1
		 pr = pRef(n1) + frac * (pRef(n2) - pRef(n1))
	ELSE
		 pr = 0.0D0
	END IF
	sigzz = sigzz - dense * gMean * oneKm + (pr - oldpr)
	tauzz = tauzz + 0.50D0 * (sigzz + oldszz) * oneKm
	dense1 = dense2
	oldszz = sigzz
	oldpr = pr
	layer1 = layer2
200  CONTINUE
resid = zBase - z
h = zStop
z = zBase
IF (zStop <= zM) THEN
	t = TempC(h)
	dense2 = rhoBar(1) * (1.0D0 - t * alphaT(1)) + &
&             density_anomaly_kgpm3
ELSE
	t = TempM(h - zM)
	dense2 = rhoBar(2) * (1.0D0 - t * alphaT(2)) + &
&             density_anomaly_kgpm3
END IF
dense = 0.50D0 * (dense1 + dense2)
IF (z > 0.0D0) THEN
	n1 = z / oneKm
	n2 = n1 + 1
	frac = z / oneKm - n1
	pr = pRef(n1) + frac * (pRef(n2) - pRef(n1))
ELSE
	pr = 0.0D0
END IF
sigzzb = sigzz - dense * gMean * resid + (pr - oldpr)
tauzz = tauzz + 0.50D0 * (sigzzb + oldszz) * resid
RETURN

END SUBROUTINE Squeez

SUBROUTINE TauDef (alpha, eRate, mxEl, numEl, tOfset, & ! input
&                    tauMat)                              ! output

!   Computes vertical integrals of relative horizontal
!   stress anomalies (relative to vertical stress): tauMat.

!   The components are:
!   tauMat(1) = vertical integral of (Sxx-Szz)
!   tauMat(2) = vertical integral of (Syy-Szz)
!   tauMat(3) = vertical integral of Sxy.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alpha, eRate                                                     ! input
INTEGER, INTENT(IN) :: mxEl, numEl                                                     ! input
REAL*8, INTENT(IN) :: tOfset                                                           ! input
REAL*8, INTENT(OUT) :: tauMat                                                          ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, m
REAL*8 exx, exy, eyy
DIMENSION alpha(3, 3, 7, mxEl), eRate(3, 7, mxEl), &
&           tauMat(3, 7, mxEl), tOfset(3, 7, mxEl)

DO 1000 m = 1, 7
	DO 900 i = 1, numEl
		 exx = eRate(1, m, i)
		 eyy = eRate(2, m, i)
		 exy = eRate(3, m, i)
		 tauMat(1, m, i) = tOfset(1, m, i) + exx * alpha(1, 1, m, i) + &
&                     eyy * alpha(1, 2, m, i) + exy * alpha(1, 3, m, i)
		 tauMat(2, m, i) = tOfset(2, m, i) + exx * alpha(2, 1, m, i) + &
&                     eyy * alpha(2, 2, m, i) + exy * alpha(2, 3, m, i)
		 tauMat(3, m, i) = tOfset(3, m, i) + exx * alpha(3, 1, m, i) + &
&                     eyy * alpha(3, 2, m, i) + exy * alpha(3, 3, m, i)
900       CONTINUE
1000  CONTINUE
RETURN
END SUBROUTINE TauDef

SUBROUTINE THOnB (basal, continuum_LRi, &                ! input
&                   etaMax, fPSfer, glue, &
&                   iConve, &
&                   LRn, LR_set_eCreep, &
&                   mxEl, mxNode, nodes, numEl, &
&                   oVB, pulled, trHMax, v, &
&                   eta, sigHB, &                          ! output
&                   outVec)                                ! work

!   Calculates shear stresses on base of plate (sigHB), using
!   the vector velocity of the layer below (oVB), and also reports
!   the linearized coupling coefficient for next iteration (eta).

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DOUBLE PRECISION, INTENT(IN) :: basal                                                  ! input
INTEGER, INTENT(IN) :: continuum_LRi                                                   ! input
REAL*8, INTENT(IN) :: etaMax, fPSfer, glue                                             ! input
INTEGER, INTENT(IN) :: iConve                                                          ! input
INTEGER, INTENT(IN) :: LRn                                                             ! input
REAL*8, INTENT(IN) :: LR_set_eCreep                                                    ! input
INTEGER, INTENT(IN) :: mxEl, mxNode, nodes, numEl                                      ! input
REAL*8, INTENT(IN) :: oVB                                                              ! input
LOGICAL, INTENT(IN) :: pulled                                                          ! input
REAL*8, INTENT(IN) :: trHMax                                                           ! input
DOUBLE PRECISION, INTENT(IN) :: v                                                      ! input
REAL*8, INTENT(OUT) :: eta, sigHB                                                      ! output
REAL*8 outVec                                                           ! work
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, LRi, m
REAL*8 shear, shear1, shear2, shear3, tiny, vcx, vcy, vMag, vmx, vmy, vrx, vry
REAL*8 t_eCreep
DIMENSION basal(2, mxNode), &
&           continuum_LRi(mxEl), &
&           eta(7, mxEl), fPSfer(2, 2, 3, 7, mxEl), &
&           glue(7, mxEl), &
&           LR_set_eCreep(0:LRn), &
&           nodes(3, mxEl), outVec(2, 7, mxEl), &
&           oVB(2, 7, mxEl), pulled(7, mxEl), &
&           sigHB(2, 7, mxEl), &
&           v(2, mxNode)

!   Small number to prevent division by zero:
DATA tiny / 2.0D-38 /

IF (iConve /= 6) THEN
!           older code, for defined lower-mantle velocity field.
!           First, interpolate surface flow to integration points:
	CALL Flow (fPSfer, mxEl, mxNode, nodes, numEl, v, & ! input
&                 outVec)                                  ! output
	DO 1000 i = 1, numEl
		 !Extract desired rheology for this continuum element:
		 LRi = continuum_LRi(i)
		 t_eCreep = LR_set_eCreep(LRi)
		 !- - - - - - - - - - - - - - - - -
		 DO 900 m = 1, 7
			  IF (pulled(m, i)) THEN
				   vcx = outVec(1, m, i)
				   vcy = outVec(2, m, i)
				   vmx = oVB(1, m, i)
				   vmy = oVB(2, m, i)
				   vrx = vmx - vcx
				   vry = vmy - vcy
				   vMag = SQRT(vrx**2 + vry**2)
				   IF (vMag > 0.0D0) THEN
						shear1 = glue(m, i) * vmag**t_eCreep
				   ELSE
						shear1 = 0.0D0
				   END IF
				   shear2 = trHMax
				   shear3 = etaMax * vMag
				   shear = MIN(shear1, shear2, shear3)
				   eta(m, i) = shear / MAX(tiny, vMag)
				   eta(m, i) = MIN(eta(m, i), etaMax)
				   sigHB(1, m, i) = eta(m, i) * vrx
				   sigHB(2, m, i) = eta(m, i) * vry
			  ELSE
				   eta(m, i) = 0.0D0
				   sigHB(1, m, i) = 0.0D0
				   sigHB(2, m, i) = 0.0D0
			  END IF
900            CONTINUE
1000       CONTINUE
ELSE
!           New code for iConve == 6: use nodal values of shear traction
!           vectors contained in BASAL, and interpolate:
	CALL Flow (fPSfer, mxEl, mxNode, nodes, numEl, basal, & ! input
&                 sigHB)                                       ! output
!           Next, interpolate surface velocity (as above) to compare
!           to values in oVB, for computation of ETA:
	CALL Flow (fPSfer, mxEl, mxNode, nodes, numEl, v, & ! input
&                 outVec)                                  ! output
	DO 2000 i = 1, numEl
		 DO 1900 m = 1, 7
			  vcx = outVec(1, m, i)
			  vcy = outVec(2, m, i)
			  vmx = oVB(1, m, i)
			  vmy = oVB(2, m, i)
			  vrx = vmx - vcx
			  vry = vmy - vcy
			  vMag = SQRT(vrx**2 + vry**2)
			  IF (vMag > 0.0D0) THEN
				   shear1 = SQRT(sigHB(1, m, i)**2 + sigHB(2, m, i)**2)
			  ELSE
				   shear1 = 0.0D0
			  END IF
			  shear2 = trHMax
			  shear3 = etaMax * vmag
			  shear = MIN(shear1, shear2, shear3)
			  eta(m, i) = shear / MAX(tiny, vmag)
			  eta(m, i) = MIN(eta(m, i), etaMax)
1900            CONTINUE
2000       CONTINUE
END IF
RETURN
END SUBROUTINE THOnB

SUBROUTINE Tract (iUnitR, iUnitT, nPlate, numNod, & ! input
&                   slab_q, whichP, xNode, yNode, &
&                   basal)                            ! output

!      Requests file name of an existing torque report
!     (including traction pole vectors for each plate,
!      created by a previous experiment with -Shells-,
!      usually one that had trHMax = 0. and extra internal
!      velocity boundary conditions for each slabless plate).
!      Reads this file, extracts the traction pole vectors,
!      and uses them to precompute basal shear tractions
!      on each node.
!      For further clarification of "traction pole vectors"
!      see subprogram -Twist- below.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iUnitR, iUnitT, nPlate, numNod                                  ! input
LOGICAL, INTENT(IN) :: slab_q                                                          ! input
INTEGER, INTENT(IN) :: whichP                                                          ! input
REAL*8, INTENT(IN) :: xNode, yNode                                                     ! input
DOUBLE PRECISION, INTENT(OUT) :: basal                                                 ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CHARACTER*132 line, trqFil
INTEGER i, ios, iPlate, j
REAL*8 equat, lat, length, lon, t, tequat, &
	& tvec, uphi, utheta, uvec

DIMENSION tvec(3), uphi(3), utheta(3), uvec(3)
DIMENSION slab_q(nPlate), whichP(numNod), &
&           xNode(numNod), yNode(numNod)
DIMENSION basal(2, numNod)
LOGICAL, DIMENSION(:), ALLOCATABLE :: tpread
REAL*8, DIMENSION(:, :), ALLOCATABLE :: tpvecs


ALLOCATE (tpvecs(3, nPlate))
ALLOCATE (tpread(nPlate))
!      Zero whole array; advisable because some plates may not
!         appear in report.
DO 30 j = 1, nPlate
	DO 20 i = 1, 3
		 tpvecs(i, j) = 0.0D0
20       CONTINUE
	tpread(j) = .FALSE.
30  CONTINUE
!      Waste first 6 lines (titles & 2 blanks & header) of torque file:
DO 40 i = 1, 6
	READ (iUnitR, "(A)") line
40  CONTINUE
!      Loop on plates in report (up to nPlate for whole-Earth model):
DO 100 j = 1, nPlate
	READ(iUnitR, * , IOSTAT = ios)
	IF (ios == -1) GO TO 101
	READ(iUnitR, "(8X,I6)", IOSTAT = ios) iPlate
	IF (ios == -1) GO TO 101
!           Waste 23 more lines of each plate report
	DO 50 i = 1, 23
		 READ(iUnitR, * , IOSTAT = ios)
		 IF (ios == -1) GO TO 101
50       CONTINUE
	READ(iUnitR, "(56X,ES10.3,2F10.2)") t, lon, lat
!           T is magnitude, in Pa, at location 90 deg. from (LON, LAT).
	tpvecs(1, iPlate) = t * COS(lat / 57.2957795130823D0) * COS(lon / 57.2957795130823D0)
	tpvecs(2, iPlate) = t * COS(lat / 57.2957795130823D0) * SIN(lon / 57.2957795130823D0)
	tpvecs(3, iPlate) = t * SIN(lat / 57.2957795130823D0)
	tpread(iPlate) = .TRUE.
!           Waste 14 lines to get past the "=======" at the bottom of
!              each torque report:
	DO 60 i = 1, 14
		 READ(iUnitR, * , IOSTAT = ios)
		 IF (ios == -1) GO TO 101
60       CONTINUE
100  CONTINUE
101  CLOSE(iUnitR)

DO 200 i = 1, numNod
	iPlate = whichP(i)
	IF (slab_q(iPlate)) THEN
!                no need for inferred basal-strength traction:
		 basal(1, i) = 0.0D0
		 basal(2, i) = 0.0D0
	ELSE
		 IF (tpread(iPlate)) THEN

!                     Uvec is unit vector to node location:
			  uvec(1) = SIN(xNode(i)) * COS(yNode(i))
			  uvec(2) = SIN(xNode(i)) * SIN(yNode(i))
			  uvec(3) = COS(xNode(i))

!                     Tvec is cross-product with traction pole vector:
			  tvec(1) = tpvecs(2, iPlate) * uvec(3) - &
&                        tpvecs(3, iPlate) * uvec(2)
			  tvec(2) = tpvecs(3, iPlate) * uvec(1) - &
&                        tpvecs(1, iPlate) * uvec(3)
			  tvec(3) = tpvecs(1, iPlate) * uvec(2) - &
&                        tpvecs(2, iPlate) * uvec(1)
			  t = SQRT(tvec(1)**2 + tvec(2)**2 + tvec(3)**2)

!                     Unit vectors at this site (NOT a pole):
			  uphi(1) = -uvec(2)
			  uphi(2) = uvec(1)
			  equat = SIN(xNode(i))
			  uphi(1) = uphi(1) / equat
			  uphi(2) = uphi(2) / equat
			  uphi(3) = 0.0D0

			  tequat = uvec(3)
			  utheta(3) = -equat
			  utheta(1) = tequat * uvec(1) / equat
			  utheta(2) = tequat * uvec(2) / equat
			  length = SQRT(utheta(1)**2 + utheta(2)**2 + &
&                            utheta(3)**2)
			  utheta(1) = utheta(1) / length
			  utheta(2) = utheta(2) / length
			  utheta(3) = utheta(3) / length

!                     Horizontal components of shear traction:
			  basal(1, i) = tvec(1) * utheta(1) + tvec(2) * utheta(2) + &
&                           tvec(3) * utheta(3)
			  basal(2, i) = tvec(1) * uphi(1) + tvec(2) * uphi(2) + &
&                           tvec(3) * uphi(3)
		 ELSE
			  basal(1, i) = 0.0D0
			  basal(2, i) = 0.0D0
		 END IF
	END IF
200  CONTINUE

DEALLOCATE (tpread)
DEALLOCATE (tpvecs)

RETURN
END SUBROUTINE Tract

SUBROUTINE Twist (area, detJ, fPSfer, & ! input
&                   iUnitT, n, nodes, nPlate, numEl, numNod, &
&                   radius, torqBS, whichP, xNode, yNode, &
&                   twistV)               ! output

!   Computes the twist pole vector twistV(3) that will apply basal-
!   strength torque torqBS(1:3, n) to plate #n, if used in a iConve==6
!   basal boundary condition in the next run of -Shells-.

!   The area, shape, and position of plate #n are represented by
!   information in "nodes" and "whichP".

!   A twist pole vector has units of shear traction (Pa, in the SI system),
!   and can be used to compute basal shear traction according to:

!       basal_shear_traction = twistV x uvec {vector cross product},

!   where uvec is a dimensionless unit vector giving position.
!   Thus the magnitude (length) of twistV represents the largest
!   basal traction, applying to points 90 degrees from the pole.
!   The 3 components of twistV are Cartesian (x, y, z) measured
!   from the center of the planet, as in any uvec.

!   The solution is achieved by setting up a 3 x 3 linear system:

!   For further clarification, read subprogram -Tract- above.

!      torqBS(1,N)   c11 c12 c13  twistV(1)
!      torqBS(2,N) = c21 c22 c23  twistV(2)
!      torqBS(3,N)   c31 c32 c33  twistV(3)

!   which is then inverted to get twistV(1:3).

!   Each column of the [c] matrix is computed by integrating a
!   hypothetical case: for example, column 1, transposed as
!   (c11, c21, c31) gives the (x, y, z) components of the
!   basal torque on plate #n that would be produced if
!   twistV(1:3) = (1., 0., 0.).

!   Linear system of equations is solved by dgesv of Intel's MKL = Math Kernel Library.

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: area, detJ, fPSfer                                               ! input
INTEGER, INTENT(IN) :: iUnitT, n, nodes, nPlate, numEl, numNod                         ! input
REAL*8, INTENT(IN) :: radius                                                           ! input
DOUBLE PRECISION, INTENT(IN) :: torqBS                                                 ! input
INTEGER, INTENT(IN) :: whichP                                                          ! input
REAL*8, INTENT(IN) :: xNode, yNode                                                     ! input
REAL*8, INTENT(OUT) :: twistV                                                          ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Named COMMON blocks hold the fixed values of the positions,
!      weights, and nodal function values at the integration points
!      in the elements (triangular elements in BLOCK  DATA  BD1,
!      and fault elements in BLOCK  DATA  BD2).
!      Entries corresponding to BD1:
DOUBLE PRECISION points, weight
COMMON / S1S2S3 / points
COMMON / WgtVec / weight
DIMENSION points(3, 7), weight(7)
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, info, ip1, ip2, ip3, j, k, l, m, n1, n2, n3, node
INTEGER, DIMENSION(3, 3) :: ipiv ! needed by dgesv of MKL library.
LOGICAL singular
REAL*8 :: dArea, equat, forcex, forcey, length, pphi, shx, shy, tequat, ttheta
REAL*8, DIMENSION(3, 3) :: c
REAL*8, DIMENSION(3) :: fxyz, rvec, tqbs, uphi, utheta, uvec
REAL*8, DIMENSION(3, 7) :: phiM, thetaM, uvecM

DIMENSION area(numEl), &
	   & detJ(7, numEl), &
	   & fPSfer(2, 2, 3, 7, numEl), &
	   & nodes(3, numEl), &
	   & torqBS(3, nPlate), &
	   & twistV(3), &
	   & whichP(numNod), &
	   & xNode(numNod), yNode(numNod)

!      Compute the 3 hypothetical cases:
DO 100 j = 1, 3
	DO 5 i = 1, 3
		c(i, j) = 0.0D0
		twistV(i) = 0.0D0
5       CONTINUE
	twistV(j) = 1.0D0

!           Integrate over elements beLonging ENTIRELY to plate #n:
	DO 90 l = 1, numEl
		 n1 = nodes(1, l)
		 n2 = nodes(2, l)
		 n3 = nodes(3, l)
		 ip1 = whichP(n1)
		 ip2 = whichP(n2)
		 ip3 = whichP(n3)
		 IF ((ip1 == n).AND.(ip2 == n).AND.(ip3 == n)) THEN
			  CALL ElUvec(n1, n2, n3, numNod, xNode, yNode, & ! input
&                            phiM, thetaM, uvecM)                ! output
!                     Numerical integration over area, with 7 Gauss
!                     integration points:
			  DO 80 m = 1, 7
				   dArea = area(l) * detJ(m, l) * weight(m)
!                          Basal shear tractions for this case,
!                          ...in 3-D (x, y, z):
				   fxyz(1) = twistV(2) * uvecM(3, m) - &
&                             twistV(3) * uvecM(2, m)
				   fxyz(2) = twistV(3) * uvecM(1, m) - &
&                             twistV(1) * uvecM(3, m)
				   fxyz(3) = twistV(1) * uvecM(2, m) - &
&                             twistV(2) * uvecM(1, m)
!                          ...in 2-D (X = +theta = S; Y = +phi = E):
				   shx = fxyz(1) * thetaM(1, m) + fxyz(2) * thetaM(2, m) + &
&                         fxyz(3) * thetaM(3, m)
				   shy = fxyz(1) * phiM(1, m) + fxyz(2) * phiM(2, m) + &
&                         fxyz(3) * phiM(3, m)
!                          Three nodal functions:
				   DO 70 k = 1, 3
						node = nodes(k, l)

!                               Contribution to consistent nodal forces:

						forcex = dArea * (shx * fPSfer(1, 1, k, m, l) &
&                                          + shy * fPSfer(1, 2, k, m, l))
						forcey = dArea * (shx * fPSfer(2, 1, k, m, l) &
&                                          + shy * fPSfer(2, 2, k, m, l))

!                               Uvec of this node:

						ttheta = xNode(node)
						pphi = yNode(node)
						equat = SIN(ttheta)
						uvec(1) = equat * COS(pphi)
						uvec(2) = equat * SIN(pphi)
						uvec(3) = COS(ttheta)

!                               Unit vectors at this site (NOT a pole):

						uphi(1) = -uvec(2)
						uphi(2) = uvec(1)
						uphi(1) = uphi(1) / equat
						uphi(2) = uphi(2) / equat
						uphi(3) = 0.0D0
						tequat = uvec(3)
						utheta(3) = -equat
						utheta(1) = tequat * uvec(1) / equat
						utheta(2) = tequat * uvec(2) / equat
						length = SQRT(utheta(1)**2 + utheta(2)**2 + &
&                                      utheta(3)**2)
						utheta(1) = utheta(1) / length
						utheta(2) = utheta(2) / length
						utheta(3) = utheta(3) / length

!                               Consistent nodal force in (x,y,z):

						fxyz(1) = forcex * utheta(1) + forcey * uphi(1)
						fxyz(2) = forcex * utheta(2) + forcey * uphi(2)
						fxyz(3) = forcex * utheta(3) + forcey * uphi(3)

!                               Nodal forces x moment arms:

						rvec(1) = radius * uvec(1)
						rvec(2) = radius * uvec(2)
						rvec(3) = radius * uvec(3)

!                               Sum up the torque for this hypothetical:

						c(1, j) = c(1, j) + &
&                                 rvec(2) * fxyz(3) - rvec(3) * fxyz(2)
						c(2, j) = c(2, j) + &
&                                 rvec(3) * fxyz(1) - rvec(1) * fxyz(3)
						c(3, j) = c(3, j) + &
&                                 rvec(1) * fxyz(2) - rvec(2) * fxyz(1)

70                      CONTINUE
80                 CONTINUE
		 END IF
90       CONTINUE
100  CONTINUE

!      Now the [c] matrix is finished.

!      Check to be sure it is not singular!
singular = (c(1, 1) <= 0.0D0).OR.(c(2, 2) <= 0.0D0).OR.(c(2, 2) <= 0.0D0)

!      If (singular) then just report a dummy torque pole at (0 E, 0 N):
IF (singular) THEN
   IF(Verbose) WRITE(iUnitVerb, *)
   IF(Verbose) WRITE(iUnitVerb, "(' CAUTION: SUBROUTINE Twist failed to define the twist pole vector of plate #',I3)") n
   IF(Verbose) WRITE (iUnitVerb, *)
   IF(Verbose) WRITE (iUnitVerb, "('CAUTION: SUBROUTINE Twist failed to define the twist pole vector of plate #',I3)") n
   twistV = (/ 1.0D0, 0.0D0, 0.0D0 /) ! with magnitude of 1 Pa, in the SI system.
   RETURN
END IF

!      Otherwise, set up the inverse problem and solve it:

tqbs(1) = torqBS(1, n)
tqbs(2) = torqBS(2, n)
tqbs(3) = torqBS(3, n)

!MKL is invoked to solve this little 3x3 linear system:
CALL dgesv( 3, 1, c, 3, ipiv, tqbs, 3, info )
!On-line reference web page:
!  https://software.intel.com/en-us/node/468876#5C18C896-C835-402E-AE63-BA7C98789A75
!    as of 2016.07.06.
twistV(1:3) = tqbs(1:3) ! move solution to desired vector for output

RETURN
END SUBROUTINE Twist

SUBROUTINE VBCs (iCond, mxBn, mxDOF, & ! input
&                  nCond, nDOF, nLB, nodCon, nUB, &
&                  vBCArg, vBCMag, &
&                  f, k)                         ! modify

!   Impose velocity boundary conditions.
!   Replace the equilibrium equation(s) for any fixed-velocity node
!   with trivial equation(s) saying that the velocity
!   is equal to that desired.  In the case of iCond(i)=1 or 3, only
!   one component is to be specified; this is done by rotating the
!   equilibrium equations to new directions (while keeping the
!   velocity variables unchanged) and replacing only the redundant
!   equation, then rotating back.  In any case, the weight used for
!   such constraint equations is equal to the largest diagonal element
!   already in the "k" matrix (to preserve its condition number).

IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, INTENT(IN) :: iCond, mxBn, mxDOF, nCond, nDOF, nLB, nodCon, nUB               ! input
REAL*8, INTENT(IN) :: vBCArg, vBCMag                                                   ! input
DOUBLE PRECISION, INTENT(INOUT) :: f, k                                                ! modify
!----------------------------------------------------------------------------
! un-named COMMON, to be placed in all programs that access the linear system:
INTEGER nRank, nCodiagonals, nKRows, iDiagonal
COMMON  nRank, nCodiagonals, nKRows, iDiagonal
!These numbers describe the shape of the banded linear system, per MKL usage.
!Values are computed by one early CALL to KSize.  Then:
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: stiff; ALLOCATE(stiff(nKRows, nRank)
!Note that element (i, j) = (row, column) of the NON-banded full linear system
!   is actually stored at: stiff(iDiagonal + i - j, j).
!----------------------------------------------------------------------------
INTEGER i, iq, ircon, iRowx, iRowy, j1, j2, jColum, node
REAL*8 vbcx, vbcy
DOUBLE PRECISION topOne
DIMENSION iCond(mxBn), f(mxDOF, 1), nodCon(mxBn), &
&           vBCArg(mxBn), vBCMag(mxBn), k(nKRows, nRank)

topOne = 0.D0
DO 10 i = 1, nDOF
   !matrix element(i, i):
	iq = iDiagonal ! == iDiagonal + i - i
	topOne = MAX(topOne, k(iq, i))
10  CONTINUE

DO 100 i = 1, nCond
	node = nodCon(i)

!    Nodes are constrained by modifying the linear system:

	IF ((iCond(i) == 1).OR.(iCond(i) == 3)) THEN
!                Impose component in the direction vBCArg,
!                but leave the perpendicular component free:
		 CALL Rotor (mxDOF, nDOF, nLB, node, & ! input
&                       nUB, vBCArg(i), &
&                       f, k)                     ! modify
		 ircon = 2 * node - 1
		 f(ircon, 1) = vBCMag(i) * topOne
		 j1 = MAX(1, ircon - nLB)
		 j2 = MIN(nDOF, ircon + nUB)
		 DO 20 jColum = j1, j2
			 !matrix element(ircon, jColum):
			  iq = iDiagonal + ircon - jColum
			  k(iq, jColum) = 0.0D0
20            CONTINUE
		!matrix element(ircon, ircon  ):
		 iq = iDiagonal    ! == iDiagonal + ircon - ircon
		 k(iq, ircon) = topOne * COS(vBCArg(i))
		!matrix element(ircon, ircon + 1):
		 iq = iDiagonal - 1 ! == iDiagonal + ircon - (ircon + 1)
		 k(iq, ircon+1) = topOne * SIN(vBCArg(i))
	ELSE IF ((iCond(i) == 2).OR.(iCond(i) == 4).OR. &
&               (iCond(i) == 5)) THEN
!                Impose both components of velocity:
		 vbcx = vBCMag(i) * COS(vBCArg(i))
		 vbcy = vBCMag(i) * SIN(vBCArg(i))
		 iRowx = 2 * node - 1
		 iRowy = 2 * node
		 f(iRowx, 1) = vbcx * topOne
		 f(iRowy, 1) = vbcy * topOne
		 j1 = MAX(1, iRowx - nLB)
		 j2 = MIN(nDOF, iRowx + nUB)
		 DO 50 jColum = j1, j2
			 !matrix element(iRowx, jColum):
			  iq = iDiagonal + iRowx - jColum
			  k(iq, jColum) = 0.0D0
50            CONTINUE
		!matrix element(iRowx, iRowx):
		 iq = iDiagonal ! == iDiagonal + iRowx - iRowx
		 k(iq, iRowx) = topOne
		 j1 = MAX(1, iRowy - nLB)
		 j2 = MIN(nDOF, iRowy + nUB)
		 DO 60 jColum = j1, j2
			 !matrix element(iRowy, jColum):
			  iq = iDiagonal + iRowy - jColum
			  k(iq, jColum) = 0.0D0
60            CONTINUE
		!matrix element(iRowy, iRowy):
		 iq = iDiagonal ! == iDiagonal + iRowy - iRowy
		 k(iq, iRowy) = topOne
	END IF
100  CONTINUE

RETURN
END SUBROUTINE VBCs

SUBROUTINE Viscos (alphaT, &                              ! input
&                    continuum_LRi, &
&                    delta_rho, eRate, &
&                    g, geothC, geothM, &
&                    LRn, LR_set_cFric, LR_set_Biot, &
&                    LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
&                    mxEl, numEl, rhoBar, rhoH2O, &
&                    sigHB, tauMat, temLim, tLInt, &
&                    visMax, zMoho, &
&                    alpha, scoreC, scoreD, tOfset, zTranC) ! output

!   Computes tactical partial-derivitive tensor alpha(1:3, 1:3, 1:7, 1:numEl)
!     (partial derivitives of vertically-integrated stresses
!      tau_ij [where normal components are relative to vertical stress]
!      with respect to strain-rates e_kl)
!      in 3 x 3 component form, from 2 x 2 principal-axis form
!      provided by -Diamnd-, at each integration point of each element.
!   Also records intercept values (tOfset(3,7,numEl)) for next iteration
!      Calculation of tauMat = tOfset + alpha * e will give model
!      relative stress integrals (relative to vertical stress integral).
!   zTranC(1:2, 1:7, 1:numEl) is the depth into the (1:crust, 2:mantle) where
!      the brittle/ductile transition occurs, for each integration point
!      of each element.  Note: "C" in the name stands for "Continuum"
!     (as opposed to fault), not for "Crust".
!   scoreC and scoreD are measures of mismatch between current
!      linearized and actual nonlinear rheologies:
!      scoreC is the maximum (absolute value) err0r in tau [N/m];
!      scoreD is the mean-err0r/mean-value [dimensionless; <=1?].

!      New version, May 5, 1998, by Peter Bird; intended to improve
!      the convergence behavior of all F-E programs which use it.
!      For an elementary (not comprehensive) test of -Viscos-,
!      see test program -Isotropy.for, 1998.4.18, which shows that
!      it preserves linear-viscous behavior in all 3 branches
!      of its code (when linear-viscous behavior is reported by -Diamnd-).
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
IMPLICIT NONE
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL*8, INTENT(IN) :: alphaT                                                              ! input
INTEGER, INTENT(IN) :: continuum_LRi                                                      ! input
REAL*8, INTENT(IN) :: delta_rho, eRate, g, geothC, geothM                                 ! input
INTEGER, INTENT(IN) :: LRn                                                                ! input
REAL*8, INTENT(IN) :: LR_set_cFric, LR_set_Biot, &                                        ! input
				   & LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep ! input
INTEGER, INTENT(IN) :: mxEl, numEl                                                        ! input
REAL*8, INTENT(IN) :: rhoBar, rhoH2O, sigHB, tauMat, temLim, tLInt, visMax, zMoho         ! input
REAL*8, INTENT(OUT) :: alpha, scoreC, scoreD, tOfset, zTranC                              ! output
!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER i, m
REAL*8 center, delp2, denom, denom0, denom1, diver, &
&        dandex, dandey, dandes, &
&        de1dex, de1dey, de1des, &
&        de2dex, de2dey, de2des, &
&        dtsde1, dtsde2, &
&        dtsdt1, dtsdt2, dtsdan, &
&        dtxde1, dtxde2, &
&        dtxdt1, dtxdt2, dtxdan, &
&        dtyde1, dtyde2, &
&        dtydt1, dtydt2, dtydan, &
&        dT1dE1, dT1dE2, dT2dE1, dT2dE2, &
&        dxx, dxy, dyy, &
&        exx, exy, eyy, e1, e2, pl0, pw0, &
&        pT1dE1, pT1dE2, pT2dE1, pT2dE2, &
&        pt1, pt2, ptxx, ptxy, ptyy, &
&        r, rho_use, rhoUse, &
&        shear, shear2, sigHBi, &
&        theta, thickC, thickM, tMean, txx, txy, tyy, &
&        zOfTop, zTran
DIMENSION alpha(3, 3, 7, mxEl), alphaT(2), &
&           continuum_LRi(mxEl), &
&           delta_rho(7, mxEl), &
&           eRate(3, 7, mxEl), &
&           geothC(4, 7, mxEl), geothM(4, 7, mxEl), &
&           LR_set_cFric(0:LRn), LR_set_Biot(0:LRn), &
&           LR_set_aCreep(1:2, 0:LRn), LR_set_bCreep(1:2, 0:LRn), LR_set_cCreep(1:2, 0:LRn), LR_set_dCreep(1:2, 0:LRn), LR_set_eCreep(0:LRn), &
&           rhoBar(2), sigHB(2, 7, mxEl), &
&           tauMat(3, 7, mxEl), temLim(2), &
&           tLInt(7, mxEl), tOfset(3, 7, mxEl), &
&           zMoho(7, mxEl), zTranC(2, 7, mxEl)
!      Internal variables:
INTEGER LRi
REAL*8 t_fric, t_Biot, t_aCreep(2), t_bCreep(2), t_cCreep(2), t_dCreep(2), t_eCreep
DIMENSION zTran(2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Initialize sums to be used in computation of scores:
scoreC = 0.0D0
scoreD = 0.0D0
denom0 = 0.0D0
denom1 = 0.0D0

DO 1000 i = 1, numEl
	!Extract desired rheology for this continuum element:
	LRi = continuum_LRi(i)
	t_fric = LR_set_cFric(LRi)
	t_Biot = LR_set_Biot(LRi)
	t_aCreep(1:2) = LR_set_aCreep(1:2, LRi)
	t_bCreep(1:2) = LR_set_bCreep(1:2, LRi)
	t_cCreep(1:2) = LR_set_cCreep(1:2, LRi)
	t_dCreep(1:2) = LR_set_dCreep(1:2, LRi)
	t_eCreep      = LR_set_eCreep(LRi)
	!Process all 7 integration points:
	DO 900 m = 1, 7

!              ----------- rheology (& zTranC) section ------------

!              Extract data for this integration point, as scalars:
		 sigHBi = SQRT(sigHB(1, m, i)**2 + &
&                         sigHB(2, m, i)**2)
		 thickC = zMoho(m, i)
		 thickM = tLInt(m, i)
		 exx = eRate(1, m, i)
		 eyy = eRate(2, m, i)
		 exy = eRate(3, m, i)

!              Guard against special case of zero strain-rate:
		 IF ((exx == 0.0D0).AND.(exy == 0.0D0).AND.(eyy == 0.0D0)) THEN
			  txx = 0.0D0
			  txy = 0.0D0
			  tyy = 0.0D0
!                   1st subscript of alpha is (1:t_xx, 2:t_yy, 3:t_xy)
!                   2nd subscript of alpha is (1:e_xx, 2:e_yy, 3:e_xy)
			  alpha(1, 1, m, i) = 4.0D0 * visMax * (thickC + thickM)
			  alpha(1, 2, m, i) = 2.0D0 * visMax * (thickC + thickM)
			  alpha(1, 3, m, i) = 0.0D0
			  alpha(2, 1, m, i) = 2.0D0 * visMax * (thickC + thickM)
			  alpha(2, 2, m, i) = 4.0D0 * visMax * (thickC + thickM)
			  alpha(2, 3, m, i) = 0.0D0
			  alpha(3, 1, m, i) = 0.0D0
			  alpha(3, 2, m, i) = 0.0D0
			  alpha(3, 3, m, i) = 2.0D0 * visMax * (thickC + thickM)
			  tOfset(1, m, i) = 0.0D0
			  tOfset(2, m, i) = 0.0D0
			  tOfset(3, m, i) = 0.0D0
			  zTranC(1, m, i) = 0.0D0
!                     Note: "C" is for Continuum (as opposed to fault), not for Crust!
!                           1st subscript is: (1:crust; 2:mantle).
			  zTranC(2, m, i) = 0.0D0
		 ELSE
!                  (Strain-rate tensor is not zero.)
!                   Find principal strain-rates (e1 <= e2)
!                   in the horizontal plane:
			  diver = exx + eyy
			  r = SQRT(exy**2 + (0.5D0 * (exx - eyy))**2)
			  e1 = 0.5D0 * diver - r
			  e2 = 0.5D0 * diver + r
			  theta = ATan2F(2.0D0 * exy, exx - eyy)
!                     See (29) of Bird (1989);
!                     theta is like angular coordinate of Mohr's circles
!                        of strain-rate and also of stress;
!                     theta = 0 when e_xx > e_yy and e_xy = 0;
!                     theta = small, + when e_xy > 0, and e_xx > e_yy;
!                     theta = Pi when e_xy = 0, and e_yy > e_xx.

!                   Prepare to sum tau (and derivitives) over layers:
			  txx = 0.0D0
			  txy = 0.0D0
			  tyy = 0.0D0
			  dT1dE1 = 0.0D0
			  dT1dE2 = 0.0D0
			  dT2dE1 = 0.0D0
			  dT2dE2 = 0.0D0

			  IF (thickC > 0) THEN
				   zOfTop = 0.0D0
				   pl0 = 0.0D0
				   pw0 = 0.0D0
				   rho_use = rhoBar(1) + delta_rho(m, i)
				   CALL Diamnd (t_aCreep(1), alphaT(1), & ! input
&                                  t_bCreep(1), t_Biot, &
&                                  t_cCreep(1), t_dCreep(1), &
&                                  t_eCreep, &
&                                  e1, e2, t_fric, g, &
&                                  geothC(1, m, i), &
&                                  geothC(2, m, i), &
&                                  geothC(3, m, i), &
&                                  geothC(4, m, i), &
&                                  pl0, pw0, &
&                                  rho_use, rhoH2O, sigHBi, &
&                                  thickC, temLim(1), &
&                                  visMax, zOfTop, &
&                                  pT1dE1, pT1dE2, &      ! output
&                                  pT2dE1, pT2dE2, &
&                                  pt1, pt2, zTran(1))
				   center = 0.5D0 * (pt1 + pt2)
				   shear = 0.5D0 * (pt2 - pt1)
				   ptxx = center + shear * COS(theta)
				   ptyy = center - shear * COS(theta)
				   ptxy = shear * SIN(theta)
!                          Add contribution of crust to total:
				   txx = txx + ptxx
				   txy = txy + ptxy
				   tyy = tyy + ptyy
				   dT1dE1 = dT1dE1 + pT1dE1
				   dT1dE2 = dT1dE2 + pT1dE2
				   dT2dE1 = dT2dE1 + pT2dE1
				   dT2dE2 = dT2dE2 + pT2dE2
				   zTranC(1, m, i) = zTran(1)
			  ELSE
				   zTranC(1, m, i) = 0.0D0
			  END IF

			  IF (thickM > 0) THEN
				   zOfTop = thickC
				   pw0 = rhoH2O * g * thickC
				   tMean = geothC(1, m, i) + &
&                       0.5D0 * geothC(2, m, i) * thickC + &
&                     0.333D0 * geothC(3, m, i) * thickC**2 + &
&                      0.25D0 * geothC(4, m, i) * thickC**3
				   rhoUse = rhoBar(1) * (1.0D0 - alphaT(1) * tMean)
				   pl0 = rhoUse * g * thickC
				   rho_use = rhoBar(2) + delta_rho(m, i)
				   CALL Diamnd (t_aCreep(2), alphaT(2), & ! input
&                                  t_bCreep(2), t_Biot, &
&                                  t_cCreep(2), t_dCreep(2), &
&                                  t_eCreep, &
&                                  e1, e2, t_fric, g, &
&                                  geothM(1, m, i), &
&                                  geothM(2, m, i), &
&                                  geothM(3, m, i), &
&                                  geothM(4, m, i), &
&                                  pl0, pw0, &
&                                  rho_use, rhoH2O, sigHBi, &
&                                  thickM, temLim(2), &
&                                  visMax, zOfTop, &
&                                  pT1dE1, pT1dE2, &       ! output
&                                  pT2dE1, pT2dE2, &
&                                  pt1, pt2, zTran(2))
				   center = 0.5D0 * (pt1 + pt2)
				   shear = 0.5D0 * (pt2 - pt1)
				   ptxx = center + shear * COS(theta)
				   ptyy = center - shear * COS(theta)
				   ptxy = shear * SIN(theta)
				   txx = txx + ptxx
				   txy = txy + ptxy
				   tyy = tyy + ptyy
				   dT1dE1 = dT1dE1 + pT1dE1
				   dT1dE2 = dT1dE2 + pT1dE2
				   dT2dE1 = dT2dE1 + pT2dE1
				   dT2dE2 = dT2dE2 + pT2dE2
				   zTranC(2, m, i) = zTran(2)
			  ELSE
				   zTranC(2, m, i) = 0.0D0
			  END IF

!              ---------- alpha and tOfset section -------------
!                      (cases of non-zero strain-rate)

			  IF (r <= 0.0D0) THEN
!                       Pathological case: e_xy = 0, and (e_xx == e_yy) /= 0.
!                       See notes from derivations of 18 April 1998;
!                       based on (28) of Bird(1989), but not using
!                       (29) because r = 0 and alpha is undefined.
!                       1st subscript of alpha is (1:t_xx, 2:t_yy, 3:t_xy)
!                       2nd subscript of alpha is (1:e_xx, 2:e_yy, 3:e_xy)
				  alpha(1, 1, m, i) = dT2dE2
				  alpha(1, 2, m, i) = dT1dE2
				  alpha(1, 3, m, i) = 0.0D0
				  alpha(2, 1, m, i) = dT1dE2
				  alpha(2, 2, m, i) = dT2dE2
				  alpha(2, 3, m, i) = 0.0D0
				  alpha(3, 1, m, i) = 0.0D0
				  alpha(3, 2, m, i) = 0.0D0
				  alpha(3, 3, m, i) = 0.5D0 * (dT1dE1 - dT2dE1 - dT1dE2 + dT2dE2)
			  ELSE
!                       Typical case, r > 0: see p. 3976 in Bird (1989).
				   de1dex = 0.5D0 - ((exx - eyy) / (4.0D0 * r))
				   de1dey = 0.5D0 + ((exx - eyy) / (4.0D0 * r))
				   de1des = -exy / r
				   de2dex = de1dey
				   de2dey = de1dex
				   de2des = -de1des
				   dandex = -SIN(theta) / (2.0D0 * r)
!                          Note: Formula above is equivalent to (29) of
!                          Bird (1989), but less likely to be singular.
				   dandey = -dandex
				   dandes = COS(theta) / r
!                          Note: Formula above is equivalent to (29) of
!                          Bird (1989), but less likely to be singular.
				   dtxdt1 = 0.5D0 * (1.0D0 - COS(theta))
				   dtxdt2 = 0.5D0 * (1.0D0 + COS(theta))
				   dtxdan = -txy
				   dtydt1 = dtxdt2
				   dtydt2 = dtxdt1
				   dtydan = txy
				   dtsdt1 = -0.5D0 * SIN(theta)
				   dtsdt2 = -dtsdt1
				   shear = SQRT(txy**2 + (0.5D0 * (txx - tyy))**2)
				   dtsdan = shear * COS(theta)
!                          1st subscript of ALPHA is (1:TXX,2:TYY,3:TXY)
!                          2nd subscript of ALPHA is (1:EXX,2:EYY,3:EXY)
				   dtxde1 = dtxdt1 * dT1dE1 + dtxdt2 * dT2dE1
				   dtxde2 = dtxdt1 * dT1dE2 + dtxdt2 * dT2dE2
				   alpha(1, 1, m, i) = &
&                        dtxde1 * de1dex + dtxde2 * de2dex + dtxdan * dandex
				   alpha(1, 2, m, i) = &
&                        dtxde1 * de1dey + dtxde2 * de2dey + dtxdan * dandey
				   alpha(1, 3, m, i) = &
&                        dtxde1 * de1des + dtxde2 * de2des + dtxdan * dandes
				   dtyde1 = dtydt1 * dT1dE1 + dtydt2 * dT2dE1
				   dtyde2 = dtydt1 * dT1dE2 + dtydt2 * dT2dE2
				   alpha(2, 1, m, i) = &
&                        dtyde1 * de1dex + dtyde2 * de2dex + dtydan * dandex
				   alpha(2, 2, m, i) = &
&                        dtyde1 * de1dey + dtyde2 * de2dey + dtydan * dandey
				   alpha(2, 3, m, i) = &
&                        dtyde1 * de1des + dtyde2 * de2des + dtydan * dandes
				   dtsde1 = dtsdt1 * dT1dE1 + dtsdt2 * dT2dE1
				   dtsde2 = dtsdt1 * dT1dE2 + dtsdt2 * dT2dE2
				   alpha(3, 1, m, i) = &
&                        dtsde1 * de1dex + dtsde2 * de2dex + dtsdan * dandex
				   alpha(3, 2, m, i) = &
&                        dtsde1 * de1dey + dtsde2 * de2dey + dtsdan * dandey
				   alpha(3, 3, m, i) = &
&                        dtsde1 * de1des + dtsde2 * de2des + dtsdan * dandes
			  END IF

!                     ----------- tOfset section ------------------
!                          (case of non-zero strain rate)
			  tOfset(1, m, i) = txx - alpha(1, 1, m, i) * exx &
&                                      - alpha(1, 2, m, i) * eyy &
&                                      - alpha(1, 3, m, i) * exy
			  tOfset(2, m, i) = tyy - alpha(2, 1, m, i) * exx &
&                                      - alpha(2, 2, m, i) * eyy &
&                                      - alpha(2, 3, m, i) * exy
			  tOfset(3, m, i) = txy - alpha(3, 1, m, i) * exx &
&                                      - alpha(3, 2, m, i) * eyy &
&                                      - alpha(3, 3, m, i) * exy
		 END IF
!C
!                 ---------- score section -----------------

!              Build tentative denominator for score, based
!              on old values of tauMat (tau relative to vertical).
		 delp2 = (0.5D0 * (tauMat(1, m, i) + tauMat(2, m, i)))**2
		 shear2 = tauMat(3, m, i)**2 + &
&            (0.5D0 * (tauMat(1, m, i) - tauMat(2, m, i)))**2
		 denom0 = denom0 + SQRT(MAX(delp2, shear2))

!              Build alternative denominator for score, based
!              on new values of t_xx, t_xy, t_yy (tau relative to vertical).
		 delp2 = (0.5D0 * (txx + tyy))**2
		 shear2 = txy**2 + (0.5D0 * (txx - tyy))**2
		 denom1 = denom1 + SQRT(MAX(delp2, shear2))

!              Evaluate difference between old and new tau:
		 dxx = tauMat(1, m, i) - txx
		 dyy = tauMat(2, m, i) - tyy
		 dxy = tauMat(3, m, i) - txy
		 delp2 = (0.5D0 * (dxx + dyy))**2
		 shear2 = (0.5D0 * (dxx - dyy))**2 + dxy**2
		 scoreC = MAX(scoreC, SQRT(delp2), SQRT(shear2))
		 scoreD = scoreD + SQRT(MAX(delp2, shear2))

900       CONTINUE
1000  CONTINUE

!      In computing scoreD, use larger of (old, new) denominators:
denom = MAX(denom0, denom1)
IF (denom > 0.0D0) THEN
	scoreD = scoreD / denom
ELSE
	scoreD = 0.0D0
END IF

!      NOTE: scoreC is already computed in loop above.

END SUBROUTINE Viscos

LOGICAL FUNCTION Within(uvec, outline_count, plate_outline_uvecs)
! Determines whether uvec is inside the circuit of plate_outline_uvecs(1:outline_count),
! where the convention is that plate_outline_uvecs(1) == plate_outline_uvecs(outline_count).
USE DSphere ! Fortran MODULE DSphere is in file DSphere.f90, provided by Peter Bird of UCLA.
IMPLICIT NONE
REAL*8, DIMENSION(3), INTENT(IN) :: uvec
INTEGER, INTENT(IN) :: outline_count
REAL*8, DIMENSION(3, outline_count), INTENT(IN) :: plate_outline_uvecs
INTEGER :: i
REAL*8, DIMENSION(3) :: tuvec_0, tuvec_1
REAL*8 :: angle_0, angle_1, angle_sum, d_angle
angle_sum = 0.0D0
tuvec_0(1:3) = plate_outline_uvecs(1:3, 1)
angle_0 = DRelative_Compass(from_uvec = uvec, to_uvec = tuvec_0)
!result is azimuth, clockwise from N, in radians
DO i = 2, outline_count
	tuvec_1(1:3) = plate_outline_uvecs(1:3, i)
	angle_1 = DRelative_Compass(from_uvec = uvec, to_uvec = tuvec_1)
	!If uvec is inside, then typically angle_1 < angle_0 (except for cycle shifts)
	d_angle = -(angle_1 - angle_0) ! reversing sign, so d_angle will typically be positive if uvec is inside.
	d_angle = ATAN2(SIN(d_angle), COS(d_angle)) ! getting rid any cycle shifts!
	angle_sum = angle_sum + d_angle
	!prepare for next plate-boundary step:
	tuvec_0 = tuvec_1
	angle_0 = angle_1
END DO
!If uvec is inside, then angle_sum should be somewhere close to 2*Pi.
Within = (angle_sum > 3.0D0) .AND. (angle_sum < 9.0D0)
!but angle_sum will be either ~0.0 or around -2.0*Pi, if point is outside.
END FUNCTION Within


end module
