!*******************************************************************************
! Module containing Shells
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

module SHELLS

contains


subroutine Shells_v5p0(VarNames,VarValues,rpeat,SHELLSconv,OData)

       use ShellsSubs
	   use ShellSetSubs
	   use SharedVars

       implicit none
!      All variables names must be explicitly declared.
!      This catches many typos where a variable name is locally misspelled.


INTEGER,INTENT(IN) :: rpeat
LOGICAL,INTENT(INOUT) :: SHELLSconv

CHARACTER(LEN=10),DIMENSION(:),INTENT(INOUT) :: VarNames
REAL*8,DIMENSION(:),INTENT(INOUT) :: VarValues
LOGICAL :: OData

!                      Array-size statement(s):
!---------------------------------------------------------------------

!                         TYPE statements
!   for scalar variables and fixed arrays (not ALLOCATABLE arrays):

       CHARACTER*5  :: zone
       CHARACTER*8  :: date
       CHARACTER*10 :: clock_time
       CHARACTER*100 :: logFil, longer_line, shorter_line

	   INTEGER :: ios, LRi, LRn, &
                & mxBn, mxDOF, mxEl, mxFEl, mxNode, mxStar, &
                & n, n1000, nCond, nDOF, nFakeN, nFl, nLB, nPBnd, &
                & nRealN, nUB, numEl, numNod
       INTEGER, DIMENSION(8) :: dateTimeNumber

!   The following switches control the size of the log file;
!      set them .TRUE. for maximum detail, or .FALSE. for brevity:
       LOGICAL :: log_strike_adjustments = .FALSE.
       LOGICAL :: log_force_balance      = .TRUE.
       LOGICAL :: log_node_velocities    = .FALSE.
       LOGICAL :: log_element_dynamics   = .TRUE.
       LOGICAL :: log_fault_dynamics     = .TRUE.
       LOGICAL :: brief, doFB1, doFB2, doFB3, doFB4, skipBC, slab_q, sphere

!   Following statement must agree with BLOCK  DATA  BD1:
       DOUBLE PRECISION :: points, weight
!   Following statement must agree with BLOCK  DATA  BD2:
       DOUBLE PRECISION :: fPhi, fPoint, fGauss
!   Following 3-vectors accumulate components of net torque:
       DOUBLE PRECISION :: torqBS, torqCL, torqFS, torqLP, torqMD, torqSS, torqVB
!---------------------------------------------------------------------

!                        DIMENSION statements:

!   DIMENSIONs that will be ALLOCATEd based on variable mxNode:
       INTEGER, DIMENSION(:), ALLOCATABLE :: jCol1, jCol2, whichP
       LOGICAL, DIMENSION(:), ALLOCATABLE :: checkN
       REAL*8,    DIMENSION(:), ALLOCATABLE :: atNode, dQdTdA, elev, tauZZN, &
     &                                         tLNode, xNode, yNode, zMNode
       REAL*8,    DIMENSION(:), ALLOCATABLE :: density_anomaly, &
     &                                         cooling_curvature
       REAL*8,    DIMENSION(:, :), ALLOCATABLE:: dv, dVLast
       DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: basal, v, vm

!  DIMENSIONs that will be ALLOCATEd based on variable mxDOF = nDOF = nRank:
       REAL*8, DIMENSION(:, :), ALLOCATABLE :: comp
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fBase
       DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: force
       !Note that "force" will have size of ALLOCATE ( force(nRank, 1) )
       !due to requirement of MKL software that the forcing vector
       !presented at solution-time should have 2 subscripts.
       INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv

!  DIMENSIONs that will be ALLOCATEd based on variable mxBn
!   (which on 2005.08.10 was made equal to numNod):
       CHARACTER(2), DIMENSION(:), ALLOCATABLE :: savTag
       INTEGER, DIMENSION(:), ALLOCATABLE :: iCond, iEdge, nodCon
       REAL*8,    DIMENSION(:), ALLOCATABLE :: r2Edge, vBCArg, vBCMag, &
     &                                         xEdge, yEdge

!  DIMENSIONs that will be ALLOCATEd based on variable mxEl:
       INTEGER, DIMENSION(:), ALLOCATABLE :: continuum_LRi
       INTEGER, DIMENSION(:, :), ALLOCATABLE :: nodes
       LOGICAL, DIMENSION(:), ALLOCATABLE :: checkE
       LOGICAL, DIMENSION(:, :), ALLOCATABLE :: contin, edgeTS, pulled
       REAL*8,    DIMENSION(:), ALLOCATABLE :: area
       REAL*8,    DIMENSION(:, :), ALLOCATABLE :: detJ, eta, glue, sigZZI, &
     &                                            sita, tauZZI, tLInt, zMoho
       REAL*8,    DIMENSION(:, :), ALLOCATABLE :: curviness, delta_rho
       REAL*8,    DIMENSION(:, :, :), ALLOCATABLE :: dXSP, dYSP, eRate, &
     &                                               geothC, geothM, &
     &                                               oVB, outVec, sigHB, &
     &                                               tauMat, tOfset, zTranC
       REAL*8,    DIMENSION(:, :, :, :), ALLOCATABLE :: alpha
       REAL*8,    DIMENSION(:, :, :, :, :), ALLOCATABLE :: dXS, dYS, fPSfer

!  DIMENSIONs that will be ALLOCATEd based on variable mxFEl:
       INTEGER, DIMENSION(:), ALLOCATABLE :: fault_LRi
       INTEGER, DIMENSION(:, :), ALLOCATABLE :: nodeF
       LOGICAL, DIMENSION(:), ALLOCATABLE :: checkF, fSlips
       LOGICAL, DIMENSION(:, :), ALLOCATABLE :: edgeFS
       REAL*8, DIMENSION(:), ALLOCATABLE :: fLen, offset
       REAL*8, DIMENSION(:, :), ALLOCATABLE :: fArg, fDip, fIMuDZ, fPeakS, zTranF
       REAL*8, DIMENSION(:, :, :), ALLOCATABLE :: fTStar
       REAL*8, DIMENSION(:, :, :, :), ALLOCATABLE :: fC
       REAL*8, DIMENSION(:, :, :, :, :), ALLOCATABLE :: fPFlt

!   DIMENSIONs that will be ALLOCATEd based on variable mxStar:
       INTEGER, DIMENSION(:), ALLOCATABLE :: list

!   DIMENSIONs that will be ALLOCATEd based on variable LRn, using range (0:LRn):
       LOGICAL, DIMENSION(:),   ALLOCATABLE :: LR_is_defined, LR_is_used
       REAL*8,  DIMENSION(:),   ALLOCATABLE :: LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, LR_set_eCreep
       REAL*8,  DIMENSION(:,:), ALLOCATABLE :: LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep ! These have initial subscript (1:2) for crust:mantle.

!  DIMENSIONs that will be ALLOCATEd based on variables nPlate and nPBnd
       INTEGER, DIMENSION(:), ALLOCATABLE :: nDPlat
       REAL*8, DIMENSION(:, :), ALLOCATABLE :: pLat, pLon

!  DIMENSIONs that will be ALLOCATEd based on results from subprogram KSize:
       DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: stiff

!  DIMENSIONs fixed by a PARAMETER (not adjustable at run-time):
       DIMENSION slab_q(nPlate)
!   Following vectors collect sums of torque components:
       DIMENSION torqBS(3, nPlate), torqCL(3, nPlate), torqFS(3, nPlate), &
     &           torqLP(3, nPlate), torqMD(3, nPlate), torqSS(3, nPlate), &
     &           torqVB(3, nPlate)

!  DIMENSIONs of fixed size:
!   Following statement must agree with BLOCK  DATA  BD1:
       DIMENSION points(3, 7), weight(7)
!   Following statement must agree with BLOCK  DATA  BD2:
       DIMENSION fPhi(4, 7), fPoint(7), fGauss(7)

!---------------------------------------------------------------------
!                           COMMON statements

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

!      Named COMMON blocks hold the fixed values of the positions,
!      weights, and nodal function values at the integration points
!      in the elements (triangular elements in BLOCK  DATA  BD1,
!      and fault elements in BLOCK  DATA  BD2).
!   Entries corresponding to BD1:
       COMMON / S1S2S3 / points
       COMMON / WgtVec / weight
!   Entries corresponding to BD2:
       COMMON / SFault / fPoint
       COMMON / FPhis /  fPhi
       COMMON / FGList / fGauss

!-------------------------------------------------------------------------
!   The following are the FORTRAN input and output device numbers:
!   (avoiding 5,6 which are STDIN and STDOUT on UNIX and AIX systems!)
!             ===== REQUIRED INPUT GROUP ========
!
!   iUnitG = device number associated with the grid input file:
       integer,parameter :: iUnitG = 1
!
!   iUnitB = device number associated with the boundary-condition
!            input file:
       integer,parameter :: iUnitB = 2
!
!   iUnitP = device number associated with the parameter input file:
       integer,parameter :: iUnitP = 3
!
!   iUnitC = device number associated with digitized outlines of all plates:
       integer,parameter :: iUnitC = 8
!
!   iUnitD = device number associated with the digitised plate-pair
!            boundaries, if required by subprogram EdgeVs;
!           (for example, "PB2002_boundaries.dig").
!     (Caution: This is *not* the same file as the digitised plate
!               outline file that may be read in on devices iUnitC and=or iUnitM;
!               the two have different formats, and may both be
!               needed if plate velocity boundary
!               conditions are imposed.)
       integer,parameter :: iUnitD = 9
!
!             ===== OPTIONAL INPUT GROUP ========
!      (but, may be mandatory for certain parameter values; see iConve)
!
!   iUnitV = device number associated with the approximate velocity
!            solution (optionally used to initialize):
       integer,parameter :: iUnitV = 11
!
!   iUnitM = device number associated with mantle flow vector file,
!            or perhaps plate outlines (e.g., PB2002_boundaries.dig):
       integer,parameter :: iUnitM = 12
!
!   iUnitR = device number associated with OLD torque- and force-balance
!            report for each plate.  Only used for input if iConve == 6.
!            See also iUnitQ below.
       integer,parameter :: iUnitR = 13
!
!   iUnitLR = device number associated with non-default Lithospheric Rheologies:
       integer,parameter :: iUnitLR = 14
!
!               ======= OUTPUT FILE GROUP =============
!
!   iUnitVerb = device number associated with the logFile (ASCII text output, including large tables):
       integer,parameter :: iUnitLog = 21
!
!   iUnitS = device number associated with velocity output (solution):
       integer,parameter :: iUnitS = 22
!
!   iUnitF = device number associated with force output (reactions).
       integer,parameter :: iUnitF = 23
!
!   iUnitQ = device number associated with NEW torque- and force-balance
!            report for each plate.  Always used.  See also iUnitR above.
       integer,parameter :: iUnitQ = 24
!
!   iUnitI = device number associated with the temporary file
!           "iteration permit.txt" which is used as a flag to
!            let the user interrupt an long job without crashing it:
       integer,parameter :: iUnitI = 25

!---------------------------------------------------------------------

!                   Beginning of Executable Code (!)

       slide = subDip * 0.0174532925199433D0
!      Mark most plates as LACKING extensive attached slabs...
       DO 10 i = 1, nPlate
            slab_q(i) = .FALSE.
   10  CONTINUE
!      ...except for these particular cases:
       slab_q( 8) = .TRUE. !  8 = AU = Australia
       slab_q(14) = .TRUE. ! 14 = CL = Caroline
       slab_q(15) = .TRUE. ! 15 = CO = Cocos
       slab_q(21) = .TRUE. ! 21 = IN = India
       slab_q(22) = .TRUE. ! 22 = jf = Juan de Fuca
       slab_q(34) = .TRUE. ! 34 = NZ = Nazca
       slab_q(37) = .TRUE. ! 37 = PA = Pacific
       slab_q(39) = .TRUE. ! 39 = PS = philippine Sea
       slab_q(40) = .TRUE. ! 40 = RI = Rivera
       slab_q(46) = .TRUE. ! 46 = SS = Solomon Sea

	   WRITE (iUnitLog, 501)
       IF(Verbose) WRITE(iUnitVerb, 501)
  501  FORMAT ( &
     &' =============================================================='/ &
     &' I               Output from program Shells,                    I'/ &
     &' I   a spherical-Earth, thin-shell program for computing time-  I'/ &
     &' I     averaged (non-elastic) neotectonic deformation of plates I'/ &
     &' I     with realistic frictional/dislocation-creep rheology.    I'/ &
     &' I   Distinct thicknesses and thermal and mechanical            I'/ &
     &' I     properties are read for the crust and mantle layers      I'/ &
     &' I     of the lithosphere.                                      I'/ &
     &' I   Faults may be included, with specified dip and friction.   I'/ &
     &' I   Also, different elements *may* have different rheologies.  I'/ &
     &' I   (*This is the primary new feature in version 5.0+.)        I'/ &
     &' I   The velocity below the base of the model may be fixed,     I'/ &
     &' I     (to represent subduction and other convection),          I'/ &
     &' I     or shear traction on the base of the lithosphere may     I'/ &
     &' I     be set to zero, or basal shear traction may be adjusted  I'/ &
     &' I     (in a series of runs) to obtain desired plate velocities.I'/ &
     &' I                           by                                 I'/ &
     &' I              Peter Bird & Xianghong Kong                     I'/ &
     &' I    Department of Earth, Planetary, and Space Sciences        I'/ &
     &' I                University of California                      I'/ &
     &' I                Los Angeles, CA 90095-1567                    I'/ &
     &" I      Peter Bird's version 5.0* of 29 January 2018            I"/ &
     &' ================================================================')

      WRITE (iUnitLog, "('----------------------------------------------', &
     &                '-------------')")
      CALL Date_And_Time (date, clock_time, zone, dateTimeNumber)
      WRITE (iUnitLog, "(' Run began on ',I4,'.',I2,'.',I2,' at ',I2,':', &
     &                 I2,':',I2)") &
     &   dateTimeNumber(1), dateTimeNumber(2), dateTimeNumber(3), &
     &   dateTimeNumber(5), dateTimeNumber(6), dateTimeNumber(7)
      WRITE (iUnitLog, "('----------------------------------------------', &
     &                '-------------')")

       wedge = ABS(90.0D0 - ABS(dipMax)) * 0.0174532925199433D0

!      Preview .feg file to determine array sizes:
       IF(Verbose) WRITE(iUnitVerb, 101) iUnitG
  101  FORMAT (/' Attempting to read finite element grid from unit',I3/)
       READ (iUnitG, * , IOSTAT = ios)

       IF (ios /= 0) THEN
            ErrorMsg = "ERROR: File not found, or file is empty or file is too short"
			call FatalError(ErrorMsg,ThID)
       END IF
       READ (iUnitG, * , IOSTAT = ios) numNod
       IF (ios /= 0) THEN
            ErrorMsg = "ERROR: File not found, or file is empty or file is too short"
            call FatalError(ErrorMsg,ThID)
       END IF
       mxNode = numNod
       mxDOF = 2 * mxNode ! = nDOF; = nRank (later synonyms)
       mxBn = numNod
!      Which permits any/all nodes to have boundary conditions!
!     (This is unphysical, but useful for computing reaction forces.)

       DO 102 i = 1, numNod
            READ (iUnitG, * , IOSTAT = ios)
            IF (ios /= 0) THEN
                 ErrorMsg = "ERROR: File not found, or file is empty or file is too short"
				 call FatalError(ErrorMsg,ThID)
            END IF
  102  CONTINUE
       READ (iUnitG, * , IOSTAT = ios) numEl
       IF (ios /= 0) THEN
            ErrorMsg = "ERROR: File not found, or file is empty or file is too short"
            call FatalError(ErrorMsg,ThID)
       END IF
       mxEl = numEl
       !Initialize survey to find LRn = MAX(continuum_LRi(1:mxEl), fault_LRi(1:MXFel)
       LRn = 0 ! until incremented below...
       DO 103 i = 1, numEl
            READ (iUnitG, "(A)", IOSTAT = ios) longer_line
            IF (ios /= 0) THEN
                 ErrorMsg = "ERROR: File not found, or file is empty or file is too short"
				 call FatalError(ErrorMsg,ThID)
            END IF
            CALL Extract_LRi(longer_line, LRi, shorter_line)
            LRn = MAX(LRn, LRi)
  103  CONTINUE
       nFl = 0
       READ (iUnitG, * , IOSTAT = ios) n
       IF (ios == 0) nFl = n
       nFl = MAX(nFl, 0)
       IF (nFl == 0) log_fault_dynamics = .FALSE.
       mxFEl = nFl
       mxStar = 20
       nPBnd = 1250
       DO 105 i = 1, nFl
            READ (iUnitG, "(A)", IOSTAT = ios) longer_line
            IF (ios /= 0) THEN
                 ErrorMsg = "ERROR: File not found, or file is empty or file is too short"
				 call FatalError(ErrorMsg,ThID)
            END IF
            CALL Extract_LRi(longer_line, LRi, shorter_line)
            LRn = MAX(LRn, LRi)
  105  CONTINUE
       REWIND (UNIT = iUnitG) ! to prepare for CALL GetNet, below...

!      ALLOCATE adjustable arrays (except those whose sizes
!                                  are based on results from CALL KSize):

!   DIMENSIONs using size variable mxNode:
       ALLOCATE (atNode(mxNode), &
     &           basal(2, mxNode), &
     &           checkN(mxNode), &
     &           cooling_curvature(mxNode), &
     &           density_anomaly(mxNode), &
     &           dQdTdA(mxNode), &
     &           dv(2, mxNode), dVLast(2, mxNode), &
     &           elev(mxNode), jCol1(mxNode), jCol2(mxNode), &
     &           tauZZN(mxNode), tLNode(mxNode), &
     &           v(2, mxNode), vM(2, mxNode), &
     &           whichP(mxNode), &
     &           xNode(mxNode), yNode(mxNode), zMNode(mxNode))

!  DIMENSIONs using size variable mxDOF = nDOF = nRank:
       ALLOCATE (comp(6, mxDOF), fBase(mxDOF), ipiv(mxDOF))
       ALLOCATE (force(mxDOF, 1))
!      Note that MKL software requires the forcing vector
!      presented at solution-time to have 2 subscripts.

!  DIMENSIONs using size variable mxBn:
       ALLOCATE (iCond(mxBn),  iEdge(mxBn), &
     &           nodCon(mxBn), r2Edge(mxBn), &
     &           savTag(mxBn), &
     &           vBCArg(mxBn), vBCMag(mxBn), &
     &           xEdge(mxBn),  yEdge(mxBn))

!  DIMENSIONs using size variable mxEl:
       ALLOCATE (alpha(3, 3, 7, mxEl), area(mxEl), &
     &           checkE(mxEl), contin(7, mxEl), continuum_LRi(mxEl), &
     &           curviness(7, mxEl), &
     &           delta_rho(7, mxEl), detJ(7, mxEl), &
     &           dXS(2, 2, 3, 7, mxEl), dYS(2, 2, 3, 7, mxEl), &
     &           dXSP(3, 7, mxEl), dYSP(3, 7, mxEl), edgeTS(3, mxEl), &
     &           eRate(3, 7, mxEl), eta(7, mxEl), &
     &           fPSfer(2, 2, 3, 7, mxEl), &
     &           geothC(4, 7, mxEl), geothM(4, 7, mxEl), &
     &           glue(7, mxEl), nodes(3, mxEl), &
     &           oVB(2, 7, mxEl), &
     &           outVec(2, 7, mxEl), pulled(7, mxEl), &
     &           sigHB(2, 7, mxEl), sigZZI(7, mxEl), sita(7, mxEl), &
     &           tauMat(3, 7, mxEl), tauZZI(7, mxEl), tLInt(7, mxEl), &
     &           tOfset(3, 7, mxEl), zMoho(7, mxEl), &
     &           zTranC(2, 7, mxEl))

!  DIMENSIONs using size variable mxFEl:
       ALLOCATE (checkF(mxFEl), edgeFS(2, mxFEl), &
     &           fault_LRi(mxFEl), fC(2, 2, 7, mxFEl), fDip(2, mxFEl), &
     &           fIMuDZ(7, mxFEl), fLen(mxFEl), fPeakS(2, mxFEl), &
     &           fPFlt(2, 2, 2, 7, mxFEl), fSlips(mxFEl), &
     &           fArg(2, mxFEl), fTStar(2, 7, mxFEl), nodeF(4, mxFEl), &
     &           offset(mxFEl), zTranF(2, mxFEl))

!  DIMENSIONs using size variable mxStar:
       ALLOCATE (list(mxStar))

!  DIMENSIONs using size variables nPlate and nPBnd:
       ALLOCATE (nDPlat(nPlate), &
     &           pLat(nPlate, nPBnd), pLon(nPlate, nPBnd))

!  DIMENSIONs using size variable LRn:
       ALLOCATE ( LR_is_defined(0:LRn) )
       ALLOCATE ( LR_is_used(0:LRn) )
       LR_is_defined = .FALSE. ! whole array, until information is read, below...
       LR_is_used    = .FALSE. ! whole array, until information is read, below...
       ALLOCATE ( LR_set_fFric(0:LRn) )
       ALLOCATE ( LR_set_cFric(0:LRn) )
       ALLOCATE ( LR_set_Biot(0:LRn) )
       ALLOCATE ( LR_set_Byerly(0:LRn) )
       ALLOCATE ( LR_set_aCreep(1:2, 0:LRn) )
       ALLOCATE ( LR_set_bCreep(1:2, 0:LRn) )
       ALLOCATE ( LR_set_cCreep(1:2, 0:LRn) )
       ALLOCATE ( LR_set_dCreep(1:2, 0:LRn) )
       ALLOCATE ( LR_set_eCreep(0:LRn) )
       !Just for ease in debugging, initialize all (currently) undefined array values as zero:
       LR_set_fFric  = 0.0D0
       LR_set_cFric  = 0.0D0
       LR_set_Biot   = 0.0D0
       LR_set_Byerly = 0.0D0
       LR_set_aCreep = 0.0D0
       LR_set_bCreep = 0.0D0
       LR_set_cCreep = 0.0D0
       LR_set_dCreep = 0.0D0
       LR_set_eCreep = 0.0D0

!   Input finite element grid and nodal data (up to 6 fields):
       CALL GetNet (iUnitG, iUnitLog, &            ! input
     &              mxDOF, mxEl, mxFEl, mxNode, &
     &              brief, continuum_LRi, cooling_curvature, &  ! output
     &              density_anomaly, &
     &              dQdTdA, elev, fault_LRi, fDip, &
     &              nFakeN, nFl, nodeF, nodes, nRealN, &
     &              numEl, numNod, n1000, offMax, offset, &
     &              title1, tLNode, xNode, yNode, zMNode, &
     &              checkE, checkF, checkN)      ! work
       IF(Verbose) WRITE(iUnitVerb, "(' Finite element grid file has been read.')")

      !Remember the default ("d_") Lithospheric Rheology as LR0, or LR_set_XXXX(0):
       LR_set_fFric(0)       = fFric
       LR_set_cFric(0)       = cFric
       LR_set_Biot(0)        = Biot
       LR_set_Byerly(0)      = Byerly
       LR_set_aCreep(1:2, 0) = aCreep(1:2)
       LR_set_bCreep(1:2, 0) = bCreep(1:2)
       LR_set_cCreep(1:2, 0) = cCreep(1:2)
       LR_set_dCreep(1:2, 0) = dCreep(1:2)
       LR_set_eCreep(0)      = eCreep
       LR_is_defined(0) = .TRUE.

!    Obtain extra input file with Lithospheric Rheologies from the user:
       IF (LRn > 0) THEN
           IF(Verbose) WRITE(iUnitVerb, *)
           IF(Verbose) WRITE(iUnitVerb, "(' Lithospheric Rheology indeces from 0 to ', I8, ' are used in this .feg file.')") LRn
           IF(Verbose) WRITE(iUnitVerb, 113) iUnitLR
  113      FORMAT (/' Attempting to read table of needed Lithospheric Rheologies from unit', I3/)
           READ (iUnitLR, * , IOSTAT = ios) ! READ (and discard) column-header line at top
           IF (ios /= 0) THEN
               ErrorMsg = "ERROR: File not found, or file is empty or file is too short"
			   call FatalError(ErrorMsg,ThID)
           END IF
           collect_LRs: DO
               READ (iUnitLR, *, IOSTAT = ios) i
               IF (ios /= 0) EXIT collect_LRs ! at EOF, probably
               IF ((i < 1).OR.(i > LRn)) THEN
				   WRITE(ErrorMsg,'(A,I0,A,I0,A)') "ERROR: LR# ", i," is outside the legal range of (1:", LRn," To make it legal, some element in the .feg file must use this (or higher) LR#"
				   call FatalError(ErrorMsg,ThID)
               END IF
               BACKSPACE(iUnitLR)
               READ (iUnitLR, *, IOSTAT = ios) i, LR_set_fFric(i), LR_set_cFric(i), LR_set_Biot(i), LR_set_Byerly(i), &
                                            &     LR_set_aCreep(1:2, i), LR_set_bCreep(1:2, i), LR_set_cCreep(1:2, i), LR_set_dCreep(1:2, i), &
                                            &     LR_set_eCreep(i)
               IF (ios == 0) THEN
                   LR_is_defined(i) = .TRUE.
               ELSE
				   WRITE(ErrorMsg,'(A,I0)') "ERROR: while trying to read 13 REAL*8 values that make up LR# ",i
				   call FatalError(ErrorMsg,ThID)
               END IF
           END DO collect_LRs
           CLOSE (iUnitLR)
           !Now, "stress-test" the continuum elements to be sure that each has a defined rheology:
           DO j = 1, numEl
               i = continuum_LRi(j)
               IF (.NOT.LR_is_defined(i)) THEN
				   WRITE(ErrorMsg,'(A,I0,A,I0,A)') "ERROR: Continuum element ", j," uses LR# ", i," which has NOT been defined!"
				   call FatalError(ErrorMsg,ThID)
               ELSE
                   LR_is_used(i) = .TRUE.
               END IF
           END DO
           !Now, "stress-test" the fault elements to be sure that each has a defined rheology:
           IF (nFl > 0) THEN
               DO j = 1, nFl
                   i = fault_LRi(j)
                   IF (.NOT.LR_is_defined(i)) THEN
                       WRITE(ErrorMsg,'(A,I0,A,I0,A)') "ERROR: Fault element ", j," uses LR# ", i," which has NOT been defined!"
					   call FatalError(ErrorMsg,ThID)
                   ELSE
                       LR_is_used(i) = .TRUE.
                   END IF
               END DO
           END IF
           !Write a report to the log-file, to provide a record of the LRs used:
           WRITE (iUnitLog, *)
           WRITE (iUnitLog, "('===========================================================================================================================')")
           WRITE (iUnitLog, "('Table of alternative Lithospheric Rheologies defined and used:')")
           WRITE (iUnitLog, "('     LR# fFric cFric  Biot Byerly aCreep(1) aCreep(2) bCreep(1) bCreep(2) cCreep(1) cCreep(2) dCreep(1) dCreep(2)    eCreep')")
           DO i = 0, LRn
               IF (LR_is_defined(i).AND.LR_is_used(i)) THEN
                   WRITE (iUnitLog, "(I8, F6.3, F6.3, F6.3, F7.3, ES10.2, ES10.2, F10.0, F10.0, F10.4, F10.4, ES10.2, ES10.2, F10.5)") &
                         & i, LR_set_fFric(i), LR_set_cFric(i), LR_set_Biot(i), LR_set_Byerly(i), &
                         & LR_set_aCreep(1:2, i), LR_set_bCreep(1:2, i), LR_set_cCreep(1:2, i), LR_set_dCreep(1:2, i), &
                         & LR_set_eCreep(i)
               END IF
           END DO
           WRITE (iUnitLog, "('===========================================================================================================================')")
           WRITE (iUnitLog, *)
       END IF ! LRn > 0

!   Check grid topology and compute geometric properties:

       IF(Verbose) WRITE(iUnitVerb, "(/' Analyzing grid topology for defects...')")
       CALL Square (brief, fDip, iUnitLog, &                 ! input
     &              log_strike_adjustments, &
     &              mxBn, mxEl, mxFEl, mxNode, &
     &              mxStar, nFl, nodeF, nodes, &
     &              numEl, numNod, skipBC, radius, wedge, &
     &              xNode, yNode, &                        ! modify
     &              area, detJ, &                          ! output
     &              dXS, dYS, dXSP, dYSP, edgeFS, &
     &              edgeTS, fLen, fPFlt, fPSfer, &
     &              fArg, nCond, nodCon, sita, &
     &              checkN, list)                          ! work
       IF(Verbose) WRITE(iUnitVerb, "(' Grid topology has been verified.')")

!   Read plate outlines, for -Assign-ing each node to a plate:
       IF(Verbose) WRITE(iUnitVerb, 2) iUnitC
    2  FORMAT (/' Attempting to read OUTLINES of plates from unit', I3/)
       CALL GetPBx (iUnitC, iUnitLog, names, nPBnd, nPlate, & ! input
     &              nDPlat, pLat, pLon)                     ! output

!   Assign each node of grid to a plate:
       IF(Verbose) WRITE(iUnitVerb, "(/' Assigning each node to a plate...')")
       CALL Assign (iUnitLog, &                                 ! input
     &              nPBnd, nDPlat, nFl, nodeF, nodes, &
     &              nPlate, numEl, numNod, &
     &              pLat, pLon, &
     &              xNode, yNode, &
     &              whichP, &                                 ! output
     &              checkN)                                   ! work
       IF(Verbose) WRITE(iUnitVerb, "(' Nodes have all been assigned.')")

IF(rpeat>1) iConve = 6  ! after first iteration use generated qEarth file from previous
       IF (iConve == 6) THEN
!           Note that -Tract- will request name of torque report file,
!           read it in, and compute values for Basal:
            CALL Tract(iUnitR, iUnitLog, nPlate, numNod, & ! input
     &                 slab_q, whichP, xNode, yNode, &
     &                 basal)                            ! output
       ELSE
            DO 3 i = 1, numNod
                 basal(1, i) = 0.0D0
                 basal(2, i) = 0.0D0
   3        CONTINUE
       END IF

!   Determine if grid covers whole sphere; if so, boundary
!      conditions will be required for footwalls of thrusts
!     (because they represent truncated subducting slabs):

       IF (nCond == 0) THEN
            sphere = .TRUE.
            skipBC = .FALSE.
            CALL Downer (skipBC, fDip, iUnitLog, mxBn, mxFEl, mxNode, & ! input
     &                   nFl, nodeF, numNod, slide, &
     &                   xNode, yNode, &
     &                   nCond, nodCon, &                             ! output
     &                   checkN)                                      ! work
       ELSE
            sphere = .FALSE.
       END IF

!   Attempt to read old velocity solution for initialization;
!     if this fails, set velocities to zero:

       IF(Verbose) WRITE(iUnitVerb, 4) iUnitV
    4  FORMAT (/ /' Attempting to read old velocity solution', &
     &            ' from unit ',I3 &
     &           /' (If none is available, give a', &
     &            ' non-existent filename, like X.)' &
     &         /)
       CALL OldVel (iUnitLog, iUnitV, mxNode, numNod, & ! input
     &              v)                                ! output

!   Read boundary conditions, in order determined by Square, or by Downer
!   (and perhaps using corrected node positions determined by Square).

       skipBC = .FALSE.
       CALL ReadBC (skipBC, fDip, iPVRef, iUnitB, iUnitD, iUnitLog, &  ! input
     &              mxBn, mxFEl, mxNode, names, nFl, &
     &              nodeF, nPlate, nRealN, numNod, n1000, omega, &
     &              radius, slide, sphere, trHMax, xNode, &
     &              yNode, &
     &              nCond, &                                         ! modify
     &              iCond, nodCon, savTag, title2, vBCArg, vBCMag, & ! output
     &              iEdge, r2Edge, xEdge, yEdge)                     ! work
       IF(Verbose) WRITE(iUnitVerb, "(' Boundary-conditions file has been read...')")

!   If necessary, average arguments of model-bounding strike-slip
!   fault elements:

       IF (.NOT.sphere) THEN
            CALL Sander (fDip, iCond, iUnitLog, &            ! input
     &                   log_strike_adjustments, &
     &                   mxBn, mxFEl, mxNode, nCond, nFl, &
     &                   nodCon, nodeF, vBCArg, vBCMag, &
     &                   wedge, xNode, yNode, &
     &                   fArg)                             ! modify
       END IF

!   Determine bandwidth of linear systems and compute storage needed:
       CALL KSize (brief, iUnitP, iUnitLog, mxEl, mxFEl, mxNode, & ! input
     &             nFl, nodeF, nodes, numEl, numNod, &
     &             nDOF, nLB, nUB, &                             ! output (+ more in un-named COMMON)
     &             jCol1, jCol2)                                 ! work

!   It is finally possible to allocate the stiffness matrix:
       ALLOCATE ( stiff(nKRows, nRank) ) ! <========= NOTE: This stiffness-matrix can easily amount to many GB of memory, requiring Win64 compilation!!!
                    ! A precise pre-estimate was provided from KSize to the screen and the log-file.

!   Interpolate and initialize all "convenience arrays":

       IF(Verbose) WRITE(iUnitVerb, "(/' Constant arrays are being computed.')")
       CALL FillIn (alphaT, basal, conduc, &                  ! input
     &              continuum_LRi, &
     &              cooling_curvature, &
     &              density_anomaly, &
     &              dQdTdA, elev, &
     &              fPSfer, gMean, gradie, &
     &              iConve, iPAfri, iPVRef, iUnitM, iUnitLog, &
     &              LRn, LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_eCreep, &
     &              mxEl, mxNode, &
     &              names, nodes, &
     &              nPlate, numEl, numNod, omega, oneKm, &
     &              radio, radius, rhoAst, rhoBar, rhoH2O, &
     &              tAdiab, temLim, tLNode, trHMax, tSurf, &
     &              vTimes, whichP, xNode, yNode, zBAsth, &
     &              zMNode, &
     &              contin, curviness, delta_rho, geothC, geothM, glue, & ! output
     &              oVB, pulled, sigZZI, &
     &              tauZZI, tauZZN, tLInt, vM, zMoho, &
     &              atNode)                                   ! work
       IF(Verbose) WRITE(iUnitVerb, "(' Constant arrays have been computed.')")

!   Compute tactical values of limits on viscosity, and weights for
!   imposition of constraints in linear systems:

       CALL Limits (area, detJ, iUnitLog, mxEl, numEl, &      ! input
     &              okDelV, radius, refStr, sphere, tLInt, &
     &              trHMax, zMoho, &
     &              constr, etaMax, fMuMax, visMax)         ! output

!   Precompute the fixed parts of the forcing vector of the linear
!     systems of equations:

       doFB1 = .TRUE.
       doFB2 = .TRUE.
       doFB3 = .TRUE.
       doFB4 = .TRUE.
       CALL Fixed (alphaT, area, conduc, &   ! input
     &             density_anomaly, detJ, &
     &             doFB1, doFB2, doFB3, doFB4, &
     &             dQdTdA, dXS, dYS, &
     &             dXSP, dYSP, edgeTS, elev, fDip, fLen, fPFlt, &
     &             fPSfer, fArg, gMean, &
     &             iCond, iUnitLog, &
     &             mxBn, mxDOF, mxEl, mxFEl, mxNode, &
     &             nCond, nFl, nodCon, nodeF, nodes, numEl, &
     &             oneKm, radio, radius, &
     &             rhoAst, rhoBar, rhoH2O, sigZZI, &
     &             sita, tauZZI, tauZZN, temLim, tLNode, tSurf, wedge, &
     &             xNode, yNode, zMNode, &
     &             fBase)                    ! output

!  -Create and solve a thin-plate version of equilibrium to determine the
!      horizontal velocity components (using iteration to handle
!      nonlinearities):

       IF(Verbose) WRITE(iUnitVerb, "(/' Beginning the iterative solution for velocity.')")
       CALL Pure (alphaT, area, &                               ! input
     &            basal, &
     &            conduc, constr, continuum_LRi, &
     &            delta_rho, detJ, dQdTdA, dXS, dYS, &
     &            elev, etaMax, everyP, &
     &            fault_LRi, fBase, fDip, fLen, fMuMax, &
     &            fPFlt, fPSfer, fArg, geothC, geothM, glue, &
     &            gMean, iCond, iConve, iUnitI, iUnitS, iUnitLog, &
     &            LRn, LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, &
     &            LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep, LR_set_eCreep, &
     &            maxItr, mxBn, mxDOF, mxEl, mxFEl, &
     &            mxNode, nCond, nDOF, nFl, nLB, nodCon, &
     &            nodeF, nodes, nUB, numEl, numNod, offMax, &
     &            offset, okToQt, oneKm, oVB, pulled, radio, &
     &            radius, rhoBar, rhoH2O, sita, slide, &
     &            tauMax, temLim, title1, &
     &            title2, title3, tLInt, tLNode, trHMax, &
     &            tSurf, vBCArg, vBCMag, visMax, &
     &            wedge, zMNode, zMoho, 999, &
     &            v, &                                          ! modify
     &            eRate, eta, fIMuDZ, fPeakS, fSlips, &         ! output
     &            sigHB, tauMat, zTranC, zTranF, &
     &            alpha, dv, dVLast, force, fC, fTStar, &       ! work
     &            outVec, stiff, ipiv, tOfset, ThID, SHELLSconv)
!   Test and display the equilibrium found:

       CALL Balanc (alphaT, area, conduc, constr, &        ! input
     &              density_anomaly, detJ, dQdTdA, dXS, &
     &              dXSP, dYS, dYSP, edgeTS, elev, eta, &
     &              fArg, fC, fDip, &
     &              fIMuDZ, fLen, fPFlt, fPSfer, fTStar, &
     &              gMean, iCond, iUnitF, &
     &              iUnitLog, log_force_balance, &
     &              mxBn, mxDOF, mxEl, mxFEl, mxNode, &
     &              nCond, nFl, nodCon, nodeF, nodes, &
     &              numEl, numNod, oneKm, oVB, radio, radius, &
     &              rhoAst, rhoBar, rhoH2O, &
     &              sigZZI, sita, &
     &              tauMat, tauZZI, tauZZN, temLim, &
     &              title1, title2, title3, tLNode, &
     &              tSurf, v, wedge, xNode, yNode, &
     &              zMNode, &
     &              sigHB, &                              ! modify
     &              comp, &                               ! output
     &              fBase, outVec)                        ! work

!   Output the solution:

       CALL Result (alphaT, area, comp, detJ, elev, eRate, everyP, & ! input
     &              fault_LRi, &
     &              fDip, fIMuDZ, fPFlt, fPeakS, fPSfer, fSlips, &
     &              fArg, geothC, iUnitQ, iUnitS, iUnitLog, &
     &              log_node_velocities, &
     &              log_element_dynamics, &
     &              log_fault_dynamics, &
     &              LRn, LR_set_fFric, &
     &              mxDOF, mxEl, mxFEl, mxNode, names, &
     &              nFl, nodeF, nodes, nPlate, nRealN, numEl, numNod, &
     &              n1000, oneKm, &
     &              radius, rhoAst, rhoBar, rhoH2O, &
     &              sigHB, tauMat, tauMax, &
     &              tauZZI, title1, title2, title3, tLInt, tLNode, &
     &              v, wedge, whichP, xNode, yNode, &
     &              zMNode, zMoho, zTranC, zTranF, &
     &              torqBS, torqCL, torqFS, torqLP, torqMD, torqSS, & ! output
     &              torqVB)

      WRITE (iUnitLog, "('----------------------------------------------', &
     &                '-------------')")
      CALL Date_And_Time (date, clock_time, zone, dateTimeNumber)
      WRITE (iUnitLog, "(' Run ended on ',I4,'.',I2,'.',I2,' at ',I2,':', &
     &                 I2,':',I2)") &
     &   dateTimeNumber(1), dateTimeNumber(2), dateTimeNumber(3), &
     &   dateTimeNumber(5), dateTimeNumber(6), dateTimeNumber(7)

       WRITE (iUnitLog, "('----------------------------------------------', &
     &                '-------------')")

       IF(Verbose) WRITE(iUnitVerb, "(' See the logFile for detailed output:')")
       IF(Verbose) WRITE(iUnitVerb, "(' ',A)") TRIM(logFil)
       CLOSE (UNIT = iUnitLog)

       DEALLOCATE (stiff)

END SUBROUTINE


end module
