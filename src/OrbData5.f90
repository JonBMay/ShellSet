!*******************************************************************************
! Module containing OrbData
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

module orbdata

contains


subroutine OrbData5(ModNum,VarNames,VarValues,ListDIR)

use ShellSetSubs
use DATA_subs
use SharedVars

implicit none

integer,intent(in) :: ModNum                            ! ETOPO file copy
character(len=10),dimension(:),intent(in) :: VarNames ! Variable update subroutine
real*8,dimension(:),intent(in) :: VarValues           ! Variable update subroutine
character(len=100),intent(in) :: ListDIR              ! ETOPO file copy
character(len=100) :: filename,dir
integer :: mxNode,mxDOF,mxBN,mxEl,mxFEl,mxStar
integer :: nFakeN,nFl,nodeF,nodes,nRealN,numEl,nEX,nEY
integer :: numNod,n1000,nCond,nodCon,list
integer :: iRow,jCol,nQX,nQY,nAX,nAY
integer :: nCX,nCY,nSX,nSY,iNode
real*8 :: TAsthK,dQdTdA,elev,fDip
real*8 :: offset,xNode,yNode,area,detJ,dXs,dYs,dXSP,dYSP,fLen,fpflt,fpsfer,fArg,sita
real*8 :: eX1,eX2,eDX,eY1,eY2,eDY
real*8 :: qX1,qX2,qDX,qY1,qY2,qDY,qLimit,aX1,aX2,aDX,aY1,aY2,aDY
real*8 :: cX1,cX2,cDX,cY1,cY2,cDY,sX1,sX2,sDX,sY1,sY2,sDY,pLon,pLat,elevat,heatFl
real*8 :: thickC,thickM,chemical_delta_rho,chemical_delta_rho_list
real*8 :: cooling_curvature,cooling_curvature_list,zMNode,tLNode

!                 PARAMETER (array-size) statements

!   Set the following PARAMETERs at least as large as your problem:

! maxNod = maximum number of nodes
!  [N.B. The use of "fake" nodes as in program -Faults- has not been
!        supported in the -Shells- family; here all nodes are "real".]
       INTEGER,PARAMETER :: maxNod = 20000

! maxBN = maximum number of boundary nodes
!  [N.B. The use of "fake" nodes as in program -Faults- has not been
!        supported in the -Shells- family; here all nodes are "real".]
       INTEGER,PARAMETER :: maxBN = 6000

! maxEl = maximum number of continuum elements (triangles).
       INTEGER,PARAMETER :: maxEl = 30000

! maxFEl = maximum number of fault elements (line segments);
       INTEGER,PARAMETER :: maxFEl = 6000

! maxAtP = maximum number of nodes which may overlap at a fault-
!          intersection point.
       INTEGER,PARAMETER :: maxAtP = 10

!---------------------------------------------------------------------
!                         TYPE statements

       INTEGER :: continuum_LRi, fault_LRi ! both DIMENSIONed below ...

       LOGICAL :: brief, log_strike_adjustments, needE, needQ, skipBC
       LOGICAL :: checkE, checkF, checkN, edgeTS, edgeFS

!   The following to agree with BLOCK  DATA  BD1:
       DOUBLE PRECISION :: points, weight
!   The following to agree with BLOCK  DATA  BD2:
       DOUBLE PRECISION :: fPhi, fPoint, fGauss

!---------------------------------------------------------------------
!                        DIMENSION statements:

!  DIMENSIONs using dynamic memory allocation:
       REAL*8, DIMENSION(:, :), ALLOCATABLE :: eArray, qArray, aArray, &
     &                                         cArray, sArray

!  DIMENSIONS using PARAMETER maxNod:
       DIMENSION checkN(maxNod), chemical_delta_rho_list(maxNod), &
     &           cooling_curvature_list(maxNod), dQdTdA(maxNod), &
     &           elev  (maxNod), &
     &           tLNode(maxNod), &
     &           xNode (maxNod), yNode (maxNod), zMNode(maxNod)

!  DIMENSIONS using PARAMETER maxBN:
       DIMENSION nodCon(maxBN)

!  DIMENSIONS USING PARAMETER maxEl:
       DIMENSION area              (maxEl), checkE            (maxEl), &
     &           continuum_LRi     (maxEl), &
     &           detJ           (7, maxEl), &
     &           dXS   (2, 2, 3, 7, maxEl), dYS   (2, 2, 3, 7, maxEl), &
     &           dXSP        (3, 7, maxEl), dYSP        (3, 7, maxEl), &
     &           edgeTS         (3, maxEl), &
     &           fpsfer(2, 2, 3, 7, maxEl), &
     &           nodes          (3, maxEl), &
     &           sita           (7, maxEl)

!  DIMENSIONS using PARAMETER maxFEl:
       DIMENSION checkF           (maxFEl), edgeFS (2, maxFEl), &
     &           fault_LRi        (maxFEl), &
     &           fDip          (2, maxFEl), &
     &           fLen             (maxFEl), &
     &           fpflt(2, 2, 2, 7, maxFEl), &
     &           fArg          (2, maxFEl), nodeF  (4, maxFEl), &
     &           offset           (maxFEl)

!  DIMENSIONS using PARAMETER maxAtP:
       DIMENSION list(maxAtP)

!  Fixed DIMENSIONs:
!   Following statement to agree with BLOCK DATA BD1:
       DIMENSION points(3, 7), weight(7)
!   Following statement to agree with BLOCK DATA BD2:
       DIMENSION fPhi(4, 7), fPoint(7), fGauss(7)

!---------------------------------------------------------------------
!                       COMMON statements

!   Following statements to agree with BLOCK DATA BD1:
       COMMON / s1s2s3 / points
       COMMON / wgtVec / weight
!   Following statements to agree with BLOCK DATA BD2:
       COMMON / sFault / fPoint
       COMMON / fPhis /  fPhi
       COMMON / fGList / fGauss
!-------------------------------------------------------------------
!                       DATA statements
!   "qLim0" is the lower limit on heat-flow for points
!      with an elevation of zero.
       REAL*8,PARAMETER :: qLim0 = 0.000D0 ! units of (watts per square meter)

!  [N.B. Former limit, in program OrbData, was:}
!      DATA qLim0 /0.045D0/ ! units of (watts per square meter)

!   "dQL_dE" is the derivitive d(qLim0)/d(elevation),
!    which adjusts the minimum heat-flow for elevation.
       REAL*8,PARAMETER :: dQL_dE = 0.00D-06 ! units of (watts per square meter)/(meter)

!  [N.B. Former limit, in program OrbData, was:}
!      DATA dQL_dE /1.43D-06/ ! units of (watts per square meter)/(meter)

!   "qLim1 is the upper limit on heat-flow for all points.
       REAL*8,PARAMETER :: qLim1 = 0.300D0 ! units of (watts per square meter)

!   "cLimit" is a lower-limit on crustal thickness:
       REAL*8,PARAMETER :: cLimit = 6570.0D0 ! meters
!    Changed from 5000. to 6570. on 2005.05.24, to agree with CRUST2.grd.

!   "hCMax" is an upper-limit on crustal thickness:
       REAL*8,PARAMETER :: hCMax = 75000.0D0
!    Already agreed with CRUST2.grd; no need to change!

!   "hLMax" is an upper-limit on total lithosphere thickness:
       REAL*8,PARAMETER :: hLMax = 400000.0D0
!    This limit prevents "lithosphere" from extending into the
!    transition zone, where new mineral phases appear,
!    and the input physical parameters would therefore not be valid.

!   "delta_rho_limit" is the maximum permitted size of chemical
!    density anomalies throughout the lithosphere:
       ! REAL*8,PARAMETER :: delta_rho_limit = 50.0D0 ! units of (kilogram per cubic meter)
       REAL*8,PARAMETER :: delta_rho_limit = 100.0D0 ! units of (kilogram per cubic meter)

!  Note that all of the above are in SI units (W/m**2, m, m, m,
!      kg/m**3, etc.)

!   "iUnitL" = Fortran device number for -Assign- log file:
       INTEGER,PARAMETER :: iUnitL = 13

!---------------------------------------------------------------------

!                   BEGINNING OF EXECUTABLE CODE

!   *** Kludge alert ***
!   Conversion of PARAMETERs (constants) to variables should logically
!   have no effect, but in fact helps to suppress some spurious
!   messages from the IBM VS-FORTRAN compiler:
       mxNode = maxNod
       mxDOF = 2 * mxNode
       mxBN  = maxBN
       mxEl  = maxEl
       mxFEl = maxFEl
       mxStar = maxAtP

       wedge = ABS(90.0D0 - ABS(dipMax)) * 0.0174532925199433D0
!---------------------------------------------------------------------


!  INTRODUCTION / HEADER SECTION:

       IF(Verbose) WRITE (iUnitVerb, 1)
    1  FORMAT ( &
     &  ' -------------------------------------------------------------' &
     & /' Program OrbData5  (version of 14 February 2019)' &
     &//' by Peter Bird, Dept. of Earth, Planetary, & Space Sciences,' &
     & /' University of California, Los Angeles, CA 90095-1567.' &
     &//' Reads a finite element grid produced by OrbWin (or OrbWeave)' &
     & /'   (that was run in Shells-mode)' &
     & /'    and a parameter file (formatted for Shells) with' &
     & /'    thermal conductivities, densities, et cetera.' &
     & /'    Then, fills-in or computes nodal data:' &
     & /'   -elevation (+ above sea level, - below);' &
     & /'   -heat-flow;' &
     & /'   -crustal thickness;' &
     & /'   -mantle lithosphere thickness (NOT including crust);' &
     & /'   -chemical density anomaly (of whole lithosphere);' &
     & /'   -geotherm curvature due to cooling (whole lithosphere).' &
     & /' This version preserves optional per-element LR#s, if present.')

       IF(Verbose) WRITE (iUnitVerb, 2)
   2   FORMAT( &
     &//' All units are SI: m, W/m**2 (NOT mW/m**2), m, m, kg/m**3,' &
     & /'    and C/m**2.' &
     & /' Note that the change in the quadratic coefficient of the' &
     & /'    geotherm caused by cooling_curvature is:' &
     & /'    delta_geoth3 = delta_geoth7 = -0.5 * cooling_curvature.' &
     &//' In the output data, the base of the mantle lithosphere' &
     & /'    will be an isothermal surface.  The isotherm value' &
     & /'    is chosen to lie on the asthenosphere adiabat' &
     & /'   (evaluated at an arbitrary depth of 100 km).')

       IF(Verbose) WRITE (iUnitVerb, 3) qLim0, dQL_dE, qLim1
    3  FORMAT ( &
     & /' If any elevation is non-zero, this elevation will be left' &
     & /'    unchanged, to preserve effects of hand-editing.' &
     & /'    If any elevation is zero, a new value will be' &
     & /'    interpolated from a gridded dataset, such as ETOPO5.grd.' &
     & /' If any nodal heat-flow is non-zero, this heat-flow will be left' &
     & /'    unchanged, to preserve effects of hand-editing.' &
     & /'    If any heat-flow is zero, a new value will be' &
     & /'    interpolated from a gridded dataset.' &
     & /' Ocean-floor heat-flow will be based on age of seafloor,' &
     & /'    where known (< 200 Ma), from a gridded-age dataset,' &
     & /'    such as age_1p5.grd derived from Mueller et al. [1997a].' &
     &//' Heat-flow (from any source) is subject to a minimum value of' &
     & /'   (',F5.3,' + ',ES10.3,' * elevation)' &
     & /'    and a maximum value of ',F5.3,'.' &
     & /' ------------------------------------------------------------')

!    Initiate log file for -Assign- output:

       IF(Verbose) WRITE (iUnitVerb, 91) iUnitL
   91  FORMAT (/ /' Attempting to create detailed log file on unit', &
     &         I3/)
       WRITE (iUnitL, "('Log file of a run of OrbData5:')")

!    Echo the limits that are compiled-in-place, for a complete record:

       IF(Verbose) WRITE (iUnitVerb, 5) qLim0, dQL_dE, qLim1, cLimit, hCMax, hLMax, &
     &                  delta_rho_limit
       WRITE (iUnitL, 5) qLim0, dQL_dE, qLim1, cLimit, hCMax, hLMax, &
     &                  delta_rho_limit
    5  FORMAT(/' The following limits apply in this run:' &
     &/'    Lower limit on heat-flow = ', F5.3, '+', ES10.3, ' * elevation' &
     &/'    Upper limit on heat-flow = ', F5.3 &
     &/'    Lower limit on crustal thickness = ', F7.0 &
     &/'    Upper limit on crustal thickness = ', F7.0 &
     &/'    Upper limit on total lithosphere thickness = ', F7.0 &
     &/'    Upper limit of chemical density variation = ', F7.0)

       WRITE (iUnitL, 92)
   92  FORMAT (/'Table of results, possibly interrupted by' &
     &         /'warning messages from subprogram -Assign-:')


!   Lithosphere/asthenosphere temperature, in Kelvin, is determined
!      from the asthenosphere adiabat in the parameter file,
!      evaluated at (rather arbitrarily) 100 km depth:

       TAsthK = TAdiab + gradie * 100.0D3 ! where 100 km is expressed in meters

!   Read finite-element grid on unit 2:

       CALL GetNet (     2, iUnitVerb, &                          ! INTENT(IN)
     &               mxDOF,  mxEl,  mxFEl, mxNode, &           ! INTENT(IN)
     &               brief, &                                  ! INTENT(OUT)
     &               continuum_LRi, &                          ! INTENT(OUT)
     &               dQdTdA,   elev, &                         ! INTENT(OUT)
     &               fault_LRi, fDip, &                        ! INTENT(OUT)
     &               nFakeN,    nFl,  nodeF,  nodes, nRealN, & ! INTENT(OUT)
     &               numEl, numNod,  n1000, offMax, offset, &  ! INTENT(OUT)
     &               title1,  xNode,  yNode, &                 ! INTENT(OUT)
     &               checkE, checkF, checkN)                   ! working arrays
!     [N.B. This version IGNORES any crustal thicknesses,
!           mantle lithosphere thicknesses, chemical density anomalies,
!           or cooling curvatures that might already be in the .FEG file.
!           New values will be computed each time OrbData5 is run.]

       IF(Verbose) WRITE (iUnitVerb, 40)
   40  FORMAT (/' Successfully read F-E grid, now verifying topology...')
       log_strike_adjustments = .FALSE.
       skipBC = .TRUE.
       CALL Square (brief, fDip, iUnitVerb, &                   ! INTENT(IN)
     &              log_strike_adjustments, &                ! INTENT(IN)
     &               mxBN,   mxEl,  mxFEl, mxNode, &         ! INTENT(IN)
     &             mxStar,    nFl,  nodeF,  nodes, &         ! INTENT(IN)
     &              numEl, numNod, skipBC, radius,  wedge, & ! INTENT(IN)
     &              xNode,  yNode,  &                        ! INTENT(INOUT)
     &               area,   detJ, &                         ! INTENT(OUT)
     &                dXS,    dYS,   dXSP,   dYSP, edgeFS, & ! INTENT(OUT)
     &             edgeTS,   fLen,  fpflt, fpsfer, &         ! INTENT(OUT)
     &               fArg,  nCond, nodCon,   sita, &         ! INTENT(OUT)
     &             checkN,   list)                           ! working arrays
       IF(Verbose) WRITE (iUnitVerb, 50)
   50  FORMAT (/' Topology of finite element grid has been verified.'/)

!   Read in elevation array on unit 3, if needed

       needE = .FALSE.
       DO 60 i = 1, numNod
            IF (elev(i) == 0.0D0) needE = .TRUE.
   60  CONTINUE
       IF (needE) THEN
	     write(dir,"(A,'/ThID_',I0,'_Data_input')") trim(ListDIR),ThID
		 write(filename,"('fort_',I0,'.3')") ModNum
	     call execute_command_line('cp INPUT/ETOPO20.grd '//trim(dir)//'/'//trim(filename))

            IF(Verbose) WRITE(iUnitVerb, 61)
   61       FORMAT(/ /' Attempting to read gridded elevations:'/)
            READ (3, * ) eX1, eDX, eX2
            READ (3, * ) eY1, eDY, eY2
            nEX = (eX2 - eX1) / eDX + 1.5D0 ! truncating to INTEGER
            nEY = (eY2 - eY1) / eDY + 1.5D0
            ALLOCATE ( eArray(nEY, nEX) )
            READ (3, * ) ((eArray(iRow, jCol), jCol = 1, nEX), iRow = 1, nEY)
       ELSE
            IF(Verbose) WRITE(iUnitVerb, 69)
   69       FORMAT (/' All nodes have non-zero elevation.' &
     &              /' No elevation grid is needed.')
       END IF

!   Read in heat-flow array on unit 4, if needed:

       needQ = .FALSE.
       DO 70 i = 1, numNod
            IF (dQdTdA(i) == 0.0D0) needQ = .TRUE.
   70  CONTINUE

       IF (needQ) THEN
            IF(Verbose) WRITE(iUnitVerb, 71)
   71       FORMAT(/ /' Attempting to read gridded heat-flow:'/)
            READ (4, * ) qX1, qDX, qX2
            READ (4, * ) qY1, qDY, qY2
            nQX = (qX2 - qX1) / qDX + 1.5D0 ! truncating to INTEGER
            nQY = (qY2 - qY1) / qDY + 1.5
            ALLOCATE ( qArray(nQY, nQX) )
            READ (4, * ) ((qArray(iRow, jCol), jCol = 1, nQX), iRow = 1, nQY)
       ELSE ! all heat-flow values at nodes are already present in .FEG file.
            IF(Verbose) WRITE(iUnitVerb, 79)
   79       FORMAT (/' All nodes have non-zero heat-flow.' &
     &              /' No heat-flow grid is needed.')
       END IF

!   Impose heat-flow limits at nodes, to prevent unreasonably
!    thick and stiff lithosphere anywhere:

       DO 80 i = 1, numNod
            qLimit = qLim0 + dQL_dE * elev(i)
            IF (dQdTdA(i) /= 0.0D0) dQdTdA(i) = MAX(dQdTdA(i), qLimit)
            dQdTdA(i) = MIN(dQdTdA(i), qLim1)
   80  CONTINUE

!      Read dataset of gridded seafloor ages, on unit 7:
!          (Note: age > 200 Ma means "unknown" or "continental".)
       IF(Verbose) WRITE(iUnitVerb, 73)
   73  FORMAT(/ /' Attempting to read gridded ages of seafloor:'/)
       READ (7, * ) aX1, aDX, aX2
       READ (7, * ) aY1, aDY, aY2
       nAX = (aX2 - aX1) / aDX + 1.5D0 ! truncating to INTEGER
       nAY = (aY2 - aY1) / aDY + 1.5D0
       ALLOCATE ( aArray(nAY, nAX) )
       READ (7, * ) ((aArray(iRow, jCol), jCol = 1, nAX), iRow = 1, nAY)

!   Read in array of crustal thickness, such as CRUST2.grd, on unit 11:

       IF(Verbose) WRITE(iUnitVerb, 81)
   81  FORMAT(/ /' Attempting to read gridded crustal thickness:'/)
       READ (11, * ) cX1, cDX, cX2
       READ (11, * ) cY1, cDY, cY2
       nCX = (cX2 - cX1) / cDX + 1.5D0 ! truncating to INTEGER
       nCY = (cY2 - cY1) / cDY + 1.5D0
       ALLOCATE ( cArray(nCY, nCX) )
       READ (11, * ) ((cArray(iRow, jCol), jCol = 1, nCX), iRow = 1, nCY)

!   Convert CRUST2.grd from km to m units:
       DO iRow = 1, nCY
           DO jCol = 1, nCX
                cArray(iRow, jCol) = cArray(iRow, jCol) * 1000.0D0
           END DO
       END DO

!   Read in array of vertical S-wave travel-time anomalies
!     (from the Moho to 400 km depth),
!      such as delta_ts.grd, on unit 12:

       IF(Verbose) WRITE(iUnitVerb, 82)
   82  FORMAT(/ /' Attempting to read gridded S-wave travel-time' &
     &          /'    anomalies (in the upper mantle):'/)
       READ (12, * ) sX1, sDX, sX2
       READ (12, * ) sY1, sDY, sY2
       nSX = (sX2 - sX1) / sDX + 1.5D0 ! truncating to INTEGER
       nSY = (sY2 - sY1) / sDY + 1.5D0
       ALLOCATE ( sArray(nSY, nSX) )
       READ (12, * ) ((sArray(iRow, jCol), jCol = 1, nSX), iRow = 1, nSY)

!   PROCESS ALL NODES EQUALLY:

       IF(Verbose) WRITE(iUnitVerb, 600)
  600  FORMAT (/' Computing layer thicknesses at all nodes...' &
     &         /'        0 nodes completed...')
       WRITE (iUnitL, 601)
  601  FORMAT ('  NODE LONGITUDE  LATITUDE      ELEV    dQdTdA', &
     &         '    zMNode    tLNode', &
     &         ' chemical_Delta_rho cooling_curvature')
       DO 680 iNode = 1, numNod
            pLon = yNode(iNode) * 57.2957795130823D0
            pLat = 90.0D0 - xNode(iNode) * 57.2957795130823D0
            elevat = elev(iNode)
            heatFl = dQdTdA(iNode)

            CALL Assign (aArray,    aX1,    aDX,    aX2,    nAX,    aDY,    aY2,    nAY, & ! INTENT(IN)
     &                   alphaT, cLimit, conduc, &                                         ! INTENT(IN)
     &                   cArray,    cX1,    cDX,    cX2,    nCX,    cDY,    cY2,    nCY, & ! INTENT(IN)
     &                   delta_rho_limit, &                                                ! INTENT(IN)
     &                   eArray,    eX1,    eDX,    eX2,    nEX,    eDY,    eY2,    nEY, & ! INTENT(IN)
     &                    gMean,  hCMax,  hLMax, &                                         ! INTENT(IN)
     &                   iUnitL, iUnitVerb, &                                                 ! INTENT(IN)
     &                    oneKm, &                                                         ! INTENT(IN)
     &                     pLon,   pLat, &                                                 ! INTENT(IN)
     &                    qLim0, dQL_dE,  qLim1, &                                         ! INTENT(IN)
     &                   qArray,    qX1,    qDX,    qX2,    nQX,    qDY,    qY2,    nQY, & ! INTENT(IN)
     &                    radio, rhoAst, rhoBar, rhoH2O, &                                 ! INTENT(IN)
     &                   sArray,    sX1,    sDX,    sX2,    nSX,    sDY,    sY2,    nSY, & ! INTENT(IN)
     &                   TAsthK, temLim,  TSurf, &                                         ! INTENT(IN)
     &                   elevat, heatFl, &                                                 ! INTENT(INOUT)
     &                   thickC, thickM, chemical_delta_rho, cooling_curvature, &          ! INTENT(OUT)
	 &         ModNum)

            elev(iNode) = elevat
            dQdTdA(iNode) = heatFl
            zMNode(iNode) = MAX(thickC, cLimit)
            tLNode(iNode) = MAX(thickM, 0.0D0)
            chemical_delta_rho_list(iNode) = chemical_delta_rho
            cooling_curvature_list(iNode) = cooling_curvature
            WRITE (iUnitL, 678) iNode, pLon, pLat, &
     &                          elev(iNode), dQdTdA(iNode), zMNode(iNode), &
     &                          tLNode(iNode), &
     &                          chemical_delta_rho_list(iNode), &
     &                          cooling_curvature_list(iNode)
  678       FORMAT (' ', I8, 2F10.4, 4ES10.2, ES19.2, ES18.2)
            IF(Verbose) WRITE(iUnitVerb, 679) iNode
  679       FORMAT ('+', I8, ' nodes completed...')
  680  CONTINUE
       IF(Verbose) WRITE(iUnitVerb, 700)
  700  FORMAT (/ /' CHECK THE LOG FILE CAREFULLY FOR PROBLEMS!')

!   OUTPUT THE MODIFIED .FEG FILE:

       call PutNet (   14, &                                   ! INTENT(IN)
     &              brief, &                                   ! INTENT(IN)
     &              continuum_LRi, &                           ! INTENT(IN)
     &              dQdTdA,   elev, &                          ! INTENT(IN)
     &              fault_LRi, fDip, &                         ! INTENT(IN)
     &              mxEl,  mxFEl, mxNode,  n1000, &            ! INTENT(IN)
     &              nFakeN,    nFl,  nodeF,  nodes, &          ! INTENT(IN)
     &              nRealN,  numEl, numNod, offset, &          ! INTENT(IN)
     &              title1, tLNode,  xNode,  yNode,  zMNode, & ! INTENT(IN)
     &              chemical_delta_rho_list, &                 ! INTENT(IN)
     &              cooling_curvature_list)                    ! INTENT(IN)

       if(Verbose) write (iUnitVerb, *)
       if(Verbose) write (iUnitVerb, "(' Job completed.')")


end subroutine OrbData5


end module
