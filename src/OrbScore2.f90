!*******************************************************************************
! Module containing OrbScore2
!
! Copyright (C) 2023 Jon May, Peter Bird, Michele Carafa
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

module OrbScore

contains


SUBROUTINE OrbScore2(VarNames,VarValues,misfits)
! Compares output from -Shells- with data from:
!     -geodetic networks,
!     -stress directions,
!     -fault slip rates,
!     -seafloor spreading rates,
!     -seismicity,
!     -upper-mantle seismic anisotropy,
! and reports both comparison tables, and summary scalar misfits.

!=========== Part of the -Shells- package of programs ===========

!   Given a finite-element grid file in the format produced by -OrbWin-
!   (or -OrbWeave-) and renumbered by -OrbNumber-, with nodal data
!   added by -OrbData5- (or -OrbData-), and node-velocity output from -Shells-,
!   computes a variety of misfit measures for the realism of the results.

!   NOTES: *Does not contain subprograms Viscos or Diamnd, hence independent
!           of changes made in May 1998, and equally compatible
!           with -Old_SHELLS- or with improved -Shells-.
!           (However, subprogram Mohr IS used to find brittle/ductile
!            transition depths in fault elements.)
!          *Reads any optional Lithospheric Rheology index numbers
!           (LR#s) in different fault and/or continuum elements
!           (a feature of Shells_v5.0+).  Uses these in CALL Mohr.

!                               by
!                          Peter Bird
!          Department of Earth, Planetary, and Space Sciences,
!    University of California, Los Angeles, California 90095-1567
! (C) Copyright 1994, 1998, 1999, 2000, 2001, 2002, 2005, 2006,
!               2015, 2018, 2019, 2021
!                by Peter Bird and the Regents of
!                 the University of California.
!        (For version date, search for "Version of" below)

!   This program was developed with support from the University of
!     California, the United States Geologic Survey, the National
!     Science Foundation, and the National Aeronautics and Space
!     Administration.
!   It is freeware, and may be copied and used without charge.
!   It may not be modified in a way which hides its origin
!     or removes this message or the copyright message.
!   It may not be resold for more than the cost of reproduction
!      and transmission.
!
!---------------------------------------------------------------------
!  *** ADDITIONAL SOFTWARE REQUIRED TO SOLVE A LINEAR SYSTEM: ***
!
! Using one routine (dsysv) from the LAPACK (Linear Analysis Package)
!   portion of the Intel Math Kernel Library (MKL):
    USE MKL95_LAPACK ! in file LAPACK.f90, available from Peter Bird (I just renamed LAPACK.inc).
!
!  {Note that program lines associated with each of these solvers are just
!   commented-out, not removed, when the solver package is changed.}
!---------------------------------------------------------------------
!  *** ADDITIONAL Fortran 90 code needed to compile:
    USE DSphere      ! provided by Peter Bird as file DSphere.f90
    USE DDislocation ! provided by Peter Bird as file DDislocation.f90
    USE DIcosahedron ! provided by Peter Bird as file DIcosahedron.f90
    USE DWeighting   ! provided by Peter Bird as file DWeighting.f90
!---------------------------------------------------------------------

use ShellSetSubs
use ScoreSubs
use SharedVars

!---------------------------------------------------------------------

!            DIFFERENCES between OrbScore and OrbScore2:
!
! (1) OrbScore was in fixed-format, due to its FORTRAN 77 heritage.  But,
!     OrbScore2 has been (nominally, minimally) converted to free-format
!     Fortran 90 code, which is clearer and easier to read and modify.
! (2) OrbScore2 uses MODULE DWeighting to weight geodesy, stress,
!          & anisotropy by geographic area (but not spreading-rates
!          or fault-slip0rates or seismicity):
!     N.B. MODULE DWeighting has the PRIVATE attribute to shield all
!          of its internal variables and arrays from accidental
!          direct access.  This is critical because MODULE DWeighting
!          creates and stores a *different* F-E grid (a uniform subdivision
!          of the icosahedron) for purposes of computing data weights.
! (3) IF (pltGeo) a new output file provides errors in predicted velocities
!          of geodetic benchmarks, for plotting (e.g., with -FiniteMap-).
!     IF (pltCos) a new output file provides coseismic part of model velocities
!          of geodetic benchmarks, for plotting (e.g., with -FiniteMap-).
!     IF (pltInt) a new output file provides interseismic part of model velocities
!          of geodetic benchmarks, for plotting (e.g., with -FiniteMap-).
! (4) OrbScore2 includes scoring for match between upper-mantle seismic
!          anisotropy (specifically, fast-direction azimuths from SKS splitting)
!          and basal shear traction vectors computed from a torque report,
!          in runs with iConve = 6.  This scoring is area- and delay-weighted.
!          N.B. This code, added 2006.11 for the Earth5 project, still needs
!          to be generalized to other iConve options for basal shear traction!
! (5) OrbScore2 reads and uses any (optional) Lithospheric Rheology index
!          codes following the definitions of continuum-triangle and/or
!          linear-fault elements (a Shells_v5.0+ feature); OrbScore does not.

!---------------------------------------------------------------------
       IMPLICIT NONE ! All variable names must be declared & typed.
!---------------------------------------------------------------------

character(LEN=10),dimension(:),intent(in) :: VarNames 
real*8,dimension(:),intent(in) :: VarValues           
real*8,dimension(7),intent(out) :: misfits
real*8 :: GeoMisfit,SprdMisfit,StrsMisfit,SlpMisfit,SeisMisfit,AniMisfit
logical :: WScores


!                 PARAMETER (ARRAY-SIZE) STATEMENTS

!   Set the following parameters at least as large as your problem:

! nPBnd = Maximim number of boundary points (per plate) in PB2002 model of Bird [2003]:
       INTEGER, PARAMETER :: nPBnd = 1250

! maxNod = Maximun number of nodes (includes both "real" & "fake"):
       INTEGER, PARAMETER :: maxNod = 17000

! maxEl = Maximun number of continuum elements (triangles):
       INTEGER, PARAMETER :: maxEl = 27000

! maxFEl = Maximun number of fault elements (lines);
       INTEGER, PARAMETER :: maxFEl = 3000

! maxBN  = Maximum number of boundary nodes:
       INTEGER, PARAMETER :: maxBN = 2000

! maxAtP = Maximum number of nodes which may overlap at a fault-
!          intersection point:
       INTEGER, PARAMETER :: maxAtP = 20

! maxGeo = Maximum number of geodetic benchmarks at which geodetic
!          horizontal velocity vectors are provided (for scoring):
       INTEGER, PARAMETER :: maxGeo = 20000

! maxStr = Maximum number of points at which most-compressive
!          horizontal principal stress direction is provided
!          (for scoring):
       INTEGER, PARAMETER :: maxStr = 10000

! maxSKS = maximum number of points at which fast-polarization
!         ("phi") azimuth and splitting time of SKS waves are
!           provided (for scoring).
       INTEGER, PARAMETER :: maxSKS = 10000

! maxDat = Maximum number of sites with data on fault slip rate:
       INTEGER, PARAMETER :: maxDat = 1000

! maxMOR = Maximum number of points at which seafloor-spreading
!          rates are provided (for scoring):
       INTEGER, PARAMETER :: maxMOR = 800

! maxEqs = Maximum number of earthquakes in catalog used for scoring:
       INTEGER, PARAMETER :: maxEqs = 30000

! maxAdj = Maximum number of (longitude/latitude/velocity-datum/
!          predicted-velocity) points passed to -Adjust- for application
!          of a rigid-body rotation-rate that minimizes the weighted
!          RMS velocity difference:
       INTEGER, PARAMETER :: maxAdj = 188000

!  "piO180" is Pi/180 (conversion from degrees to radians):
       REAL*8, PARAMETER :: piO180 = 0.0174532925199433D0
!---------------------------------------------------------------------

!         TYPE & DIMENSION STATEMENTS associated with supporting Lithospheric Rheology index #s:

       CHARACTER*80 :: longer_line, shorter_line
       INTEGER :: LRi, LRn ! <== LRn = highest LRi# found in the .FEG file; used in allocating arrays just below:
       LOGICAL, DIMENSION(:), ALLOCATABLE :: LR_is_defined, LR_is_used ! (0:LRn)
       REAL*8, DIMENSION(:), ALLOCATABLE :: LR_set_fFric, LR_set_cFric, LR_set_Biot, LR_set_Byerly, & ! (0:LRn)
                                          & LR_set_eCreep
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: LR_set_aCreep, LR_set_bCreep, LR_set_cCreep, LR_set_dCreep ! (1:2, 0:LRn)

!---------------------------------------------------------------------

!                        DIMENSION STATMENTS

!  DIMENSIONs using PARAMETER maxNod:
       LOGICAL :: checkN
       REAL*8  :: atnode, cooling_curvature, density_anomaly, dQdTdA, eDotNC, eDotNM, elev, &
                & tLNode, uvecN, xNode, yNode, zMNode
       DIMENSION atnode (maxNod), checkN (maxNod), &
     &           cooling_curvature(maxNod), density_anomaly(maxNod), &
     &           dQdTdA (maxNod), &
     &           eDotNC (maxNod), eDotNM (maxNod), &
     &           elev   (maxNod), &
     &           tLNode (maxNod), &
     &           uvecN (3, maxNod),   v (2, maxNod), &
     &           xNode (maxNod), yNode (maxNod), zMNode (maxNod)

!  DIMENSIONs using PARAMETER maxEl:
       INTEGER :: continuum_LRi(maxEl), nodes
       LOGICAL :: checkE
       REAL*8  :: area, CartR, detJ, dxs, dys, dxsp, dysp, &
                & eDotEC, eDotEM, eRate, fpsfer, sita, tLInt, xIP, yIP, zMoho
       DIMENSION area (maxEl), &
     &           CartR (3, 7, maxEl), checkE (maxEl), detJ (7, maxEl), &
     &           dxs (2, 2, 3, 7, maxEl), dys (2, 2, 3, 7, maxEl), &
     &           dxsp (3, 7, maxEl), dysp (3, 7, maxEl), edgeTs (3, maxEl), &
     &           eDotEC (maxEl), eDotEM (maxEl), &
     &           eRate (3, 7, maxEl), &
     &           fpsfer (2, 2, 3, 7, maxEl), &
     &           nodes (3, maxEl), &
     &           sita (7, maxEl), &
     &           tLInt (7, maxEl), &
     &           xIP (7, maxEl), yIP (7, maxEl), &
     &           zMoho (7, maxEl)

!  DIMENSIONs using PARAMETER maxFEl:
       INTEGER :: fault_LRi(maxFEl), nodeF
       LOGICAL :: checkF
       REAL*8  :: fc, fDip, fimudz, fLen, fpeaks, &
                & fPFlt, fArg, ftstar, offset, zTranF
       DIMENSION checkF (maxFEl), edgeFS (2, maxFEl), &
     &           fc (2, 2, 7, maxFEl), fDip (2, maxFEl), &
     &           fimudz (7, maxFEl), fLen (maxFEl), fpeaks (2, maxFEl), &
     &           fPFlt (2, 2, 2, 7, maxFEl), fSlips (maxFEl), &
     &           fArg (2, maxFEl), ftstar (2, 7, maxFEl), nodeF (4, maxFEl), &
     &           offset (maxFEl), zTranF (2, maxFEl)

!   DIMENSIONs using PARAMETER maxBN:
       INTEGER :: nodCon
       DIMENSION nodCon(maxBN)

!   DIMENSIONs using PARAMETER maxAtP:
       INTEGER :: list
       DIMENSION list (maxAtP)

!   DIMENSIONs using PARAMETER nPlate:
       INTEGER, DIMENSION(nPlate) :: nBoundaryPoints
       REAL*8, DIMENSION(nPlate, nPBnd) :: pLat, pLon

!   DIMENSIONs using PARAMETER maxGeo:
       CHARACTER*20 :: geoTag
       REAL*8 :: geoPhi, geoThe, geoVel, geoSig, geoAzi, &
               & geoUTh, geoUPh
       DIMENSION geoTag(maxGeo), geoPhi(maxGeo), geoThe(maxGeo), &
     &           geoVel(maxGeo), geoSig(maxGeo), geoAzi(maxGeo), &
     &           geoUTh(maxGeo), geoUPh(maxGeo)

!   DIMENSIONs using PARAMETER maxStr:
       CHARACTER*5 :: strTag(maxStr)
       REAL*8 :: strPhi, strThe, strQua, strArg
       DIMENSION strPhi(maxStr), strThe(maxStr), &
     &           strQua(maxStr), strArg(maxStr), strReg(maxStr)

!   DIMENSIONs using PARAMETER maxSKS:
       CHARACTER*5, DIMENSION(maxSKS) :: SKS_tag      ! 5-byte ID code for datum
       REAL*8,      DIMENSION(maxSKS) :: SKS_theta    ! colatitude, in radians
       REAL*8,      DIMENSION(maxSKS) :: SKS_phi      ! colatitude, in radians
       REAL*8,      DIMENSION(maxSKS) :: SKS_argument ! phi or fast direction, in radians counterclockwise from S
       REAL*8,      DIMENSION(maxSKS) :: SKS_delay    ! SKS splitting time, in s
       INTEGER,     DIMENSION(maxSKS) :: SKS_iPlate   ! integer identifying PB2002 plate which contains this datum
       REAL*8,   DIMENSION(2, maxSKS) :: basal_shear_tractions ! (theta = S, phi = E) components, in Pa

!   DIMENSIONs using PARAMETER maxDat:
       REAL*8 :: delVZ, fltLat, fltLon, rLat
       DIMENSION delVZ(2, maxDat), fName(maxDat), &
     &           fltLat(maxDat), fltLon(maxDat), &
     &           rLat(2, maxDat)

!   DIMENSIONs using PARAMETER maxMOR:
       CHARACTER*5 :: MORTag(maxMOR)
       REAL*8 :: MORPhi(maxMOR), MORThe(maxMOR), &
     &           MORVel(maxMOR), MORSig(maxMOR)

!   DIMENSIONs using PARAMETER maxEqs:
       REAL*8 :: eqELon, eqNLat, eqMag
       DIMENSION eqELon(maxEqs), eqNLat(maxEqs), eqMag(maxEqs)

!   DIMENSIONs using PARAMETER maxAdj:
       REAL*8 :: data, predic, longterm_at_GPS, eLon ! N.B. nLat is already declared REAL*8
       DIMENSION data(2, maxAdj), predic(2, maxAdj), longterm_at_GPS(2, maxAdj), &
     &           eLon (maxAdj), nLat (maxAdj)

!   DIMENSIONs of fixed size:
       REAL*8 :: CartVs, cUvec, equVec, rA, rB, rS, tempV, &
               & uP, uT, vf, voff, vt
       DIMENSION CartVs(3, 3), cUvec(3), &
     &           equVec(3), &
     &           rA(3), rB(3), rS(3), tempV(3), &
     &           uP(3), uT(3), vf(3), voff(3), vt(3)
!   Following statement to agree with BLOCK DATA BD1:
       DIMENSION points(3, 7), weight(7)
!   Following statement to agree with BLOCK DATA BD2:
       DIMENSION fPhi(4, 7), fPoint(7), fGauss(7)

!   ALLOCATABLE array(s) added for OrbScore2, to work with MODULE DWeighting:

       REAL*8, DIMENSION(:), ALLOCATABLE :: weights
!---------------------------------------------------------------------

!                           COMMON STATEMENTS

! Note: Un-named COMMON passes INTEGER variables used in the
!       INTEGER FUNCTION "IndexK", to avoid passing these same
!       through long sequences of subprograms.
       INTEGER :: lda, md
       COMMON lda, md

!      Named COMMON BLOCKs BD1 & BD2 hold the positions,
!      weights, and nodal function values at the integration points
!      in the elements (triangular elements in BLOCK DATA BD1, and
!      fault elements in BLOCK DATA BD2).

!   Entries corresponding to BD1:
       COMMON / S1S2S3 / points
       COMMON / WgtVec / weight

!   Entries corresponding to BD2:
       COMMON / SFault / fPoint
       COMMON / fPhis /  fPhi
       COMMON / FGList / fGauss

!--------------------------------------------------------------------

!                        DATA STATEMENTS

!  "floats" is a switch controlling whether subprogram Adjust
!   will be used to minimize the kinetic energy (RMS velocity)
!   of velocity-vector lists in two contexts:
!   (1) In computing RMS velocities of nodes;
!   (2) In computing mismatch between geodetic and model velocities.
!   For normal use, set floats =.TRUE.
!   If you are *SURE* that your geodetic data are expressed in the
!   same velocity reference frame as is used to express velocities
!   in your finite element model, select floats = .FALSE. for a
!   more discriminating (less forgiving) test of model quality.
       LOGICAL, PARAMETER :: floats = .TRUE.

!  "oezOpi" is 180/Pi  (conversion from radians to degrees):
       REAL*8, PARAMETER :: oezOpi = 57.2957795130823D0

!  "secPYr" is the number of seconds in a year:
       REAL*8, PARAMETER :: secPYr = 3.155760D7

!  "unitL" is the number of millimeters per unit length in the
!   system of units in use on Fortran devices iUnitG (finite
!   element grid), and device iUnitP (parameters),
!   and iUnitV (velocity solutions).  For example, if these are in
!   SI units (meters), then unitL = 1000.0D0
       REAL*8, PARAMETER :: unitL = 1000.0D0

!  "unitT" is the number of years per unit of time in the
!   system of units in use on Fortran devices iUnitG (finite
!   element grid), and devide iUnitP (parameters),
!   and iUnitV (velocity solutions).  For example, if these are in
!   SI units (seconds), then unitT = 3.1688087D-08
       REAL*8, PARAMETER :: unitT  = 3.1688087D-08



!   THE FOLLOWING ARE THE FORTRAN INPUT AND OUTPUT DEVICE NUMBERS:


!  "iUnitG"= device number associated with the .FEG grid input file:
       INTEGER, PARAMETER :: iUnitG = 1

!  "iUnitP"= device number associated with the .IN parameter input file.
       INTEGER, PARAMETER :: iUnitP = 2

!  "iUnitV"= device number associated with the v_____.OUT velocity
!            solution (at the nodes) file:
       INTEGER, PARAMETER :: iUnitV = 3

!  "iUnitB"= device number associated with the geodetic file of
!            benchmark velocities:
       INTEGER, PARAMETER :: iUnitB = 11

!  "iUnitS"= device number associated with the file of maximum
!            horizontal principal stress directions:
       INTEGER, PARAMETER :: iUnitS = 12

!  "iUnitD"= device number associated with the geological slip rate file
!           (in .DIG format, with two title lines per point):
       INTEGER, PARAMETER :: iUnitD = 13

!  "iUnitM"= device number associated with the file of mid-ocean
!            spreading rates:
       INTEGER, PARAMETER :: iUnitM = 14

!  "iUnitE"= device number associated with the file of earthquakes
!            to be used for seismicity scoring:
       INTEGER, PARAMETER :: iUnitE = 15

!  "iUnitO"= device number to be used for outputting the
!            .FEG file with scalar strain-rates at nodes
!            based on the earthquake catalog from iUnitE:
       INTEGER, PARAMETER :: iUnitO = 21

!  "iUnitZ"= device number to be used for outputting the
!            .FEG file with scalar strain-rates at nodes
!            based on the model (after healing faults and smoothing):
       INTEGER, PARAMETER :: iUnitZ = 22

!  "iUnitX"= Device number to be used for outputting the
!            smoothed (interseismic) velocities of each node
!            in v_____.out format, for use with FiniteMap.
!           (This special file is only created IF(EXPLOD),
!            which is true only if the "benchmarks" were
!            created by Explode3 and if they match one-for-
!            one with nodes of the .feg file.)
       INTEGER, PARAMETER :: iUnitX = 23

!  "iUnitY"= Device number to be used for outputting the
!            error in geodetic velocity predictions, in same
!            .gps format as the input dataset.
!            Note that this occurs only IF (pltGeo) .AND.
!           (some geodetic data are provided for scoring).
       INTEGER, PARAMETER :: iUnitY = 31
       LOGICAL, PARAMETER :: pltGeo = .TRUE.

!  "iUnitC"= Device number to be used for outputting the model
!            coseismic part of geodetic velocity predictions,
!            in same .gps format as the input dataset.
!            Note that this occurs only IF (pltCos) .AND.
!           (some geodetic data are provided for scoring).
       INTEGER, PARAMETER :: iUnitC = 32
       LOGICAL, PARAMETER :: pltCos = .TRUE.

!  "iUnitI"= Device number to be used for outputting the model
!            interseismic part of geodetic velocity predictions,
!            in same .gps format as the input dataset.
!            Note that this occurs only IF (pltInt) .AND.
!           (some geodetic data are provided for scoring).
       INTEGER, PARAMETER :: iUnitI = 37
       LOGICAL, PARAMETER :: pltInt = .TRUE.

!  "iUnitK"= Device number to be used for reading file with
!            upper-mantle anisotropy data, in the form of fast-
!            polarization ("phi") azimuths from SKS splitting,
!            together with SKS splitting times in s.
       INTEGER, PARAMETER :: iUnitK = 33

!  "iUnitQ"= Device number to be used for reading file with
!            torque report (obtained by running Shells version
!            of summer 2006 or later).  This is used in scoring
!            the anisotropy data read through iUnitK (above).
       INTEGER, PARAMETER :: iUnitQ = 34

!  "iUnitF"= Device number to be used for reading file with
!            outlines of the plates (e.g., PB2002_plates.dig).
!            This is used in scoring the anisotropy data read
!            through iUnitK (above).
       INTEGER, PARAMETER :: iUnitF = 35

!   iUnitLR = device number associated with non-default Lithospheric Rheologies:
       INTEGER, PARAMETER :: iUnitLR = 36

!---------------------------------------------------------------------

!                    STATEMENT FUNCTION(S):

!      Interpolation within one fault element:
       REAL*8 :: fhival, s, x1, x2
       fhival(s, x1, x2) = x1 + s * (x2 - x1)

!---------------------------------------------------------------------

!            NON-Array (scalar) specifications:
!
    CHARACTER*132 :: gpsFMT
    CHARACTER*100 :: enctit, fName, line
    CHARACTER*20  :: tag
    CHARACTER*15  :: frame
    CHARACTER*10  :: flag10
    CHARACTER*8   :: c8r1, c8r2, c8rerr, c8v1, c8v2, c8verr
    CHARACTER*7   :: c7rcha, c7vcha
    CHARACTER*4   :: outcom
    CHARACTER*3   :: c3
    CHARACTER*2   :: reg, strreg

    INTEGER :: iDepth, iEle, iHour, iMohr, ios, iPlate, iSecon, iTenth, &
             & j1, j2, jYear, k, kDay, kp1, &
             & m, mIP, minute, month, mxDOF, mxNode, mxEl, mxFEl, mxBn, mxStar, &
             &     mxGeo, mxStr, mxSKS, mxGSR, mxMOR, mxAdj, &
             & n, n1, n2, n3, n4, n_AB_data, nB, nBadRe, nBland, nCond, nFakeN, nFData, nFl, nGSDat, nInSum, nJ1, nJ2, nL1, nL2, nL3, nL4,  &
             &     node, node1, node2, nRealN, nSKS_bad, nSKS_used, numAdj, numDen, numEl, numEqs, numGeo, numMOR, numNod, numSKS, numStr, n1000, &
             & subdivision

    LOGICAL :: brief, e1Part, e2Part, eZPart, &
             & haveNV, log_strike_adjustments, showis, skipBC, sphere, vertic
    LOGICAL, DIMENSION(nPlate) :: slab_q
!   Note: The following arrays could be compressed with LOGICAL*1:
    LOGICAL :: edgeTs, edgeFS, explod, fSlips

!   N.B. DOUBLE PRECISION specifications predate the current REAL*8
!              specifications, although the meaning is the same.
    DOUBLE PRECISION :: s_dp, shrink, v
!   Following statement to agree with BLOCK DATA BD1:
    DOUBLE PRECISION :: points, weight
!   Following statement to agree with BLOCK DATA BD2:
    DOUBLE PRECISION :: fPhi, fPoint, fGauss

    REAL*8  :: a, aBadRe, aCos, angle, anisotropy, arc, arg, aRegime_denominator, aSum1, aSum1D, aSum2, aSum2D, aveC, aveM, azim, azimuth, &
             & b, bad, badMma, badSig, bAzim, bigB, &
             & chance, closes, correl, cosA, cosB, croSiz, crossP, &
             & deltaT, delVX, delVY, dExtra, dip, dot, dVAP, dVAT, dVBP, dVBT, dVP, dVT, dX, dY,  &
             & e1, e11h, e12h, e1h, e2, e22h, e2h, e2me1, e2meZ, e3, e3azim, eLarge, eLonDe, eLonP, &
             &     enB, eqM0, eqMB, eqMS, eqMW, equPar, err, exx, exy, eyy, eZ, eZme1, &
             & factor, fLong, frac, fSave, &
             & geodes, gVMma, &
             & mvmma, &
             & nLat, nLatDe, nLatP, &
             & opens, &
             & pAzim, perBad, phi, pure_bad, probRL, probVZ, &
             & r, r2, r2flt, r2min, r2n1, r2n2, rate, regime, relVZ, rErr, rho, rLt, rMin, rmsv, &
             & s1, s2, s3, seismi, short, side, sigma, sigma2, simple_bad, SIMu, sinist, slpErr, smallA, smallB, smallC, spread, sprMma, st, stress, &
             & sum, sum0s, sum1, sum1d, sum1mm, sum1s, sum2, sum2d, sum2mm, sum2s, sumB,sumB2, sumN, sumNmm, sumNs, sumSid, &
             & t, t1, t2, test, theLat, theLon, theta, thick, time1, time2, toler, tSec, &
             & u1x, u1y, u2x, u2y, unitBX, unitBY, unitX, unitY, &
             & varC, varM, vEMmpa, vESigm, vM, vMMma,  vNMmpa, vNSigm, vP, vPhi, vS, vTheta, vtl, vX1, vX2, vX3, vX4, vY1, vY2, vY3, vY4, &
             & x3, x4, xDatum, xMean, xSum, &
             & y1, y2, y3, y4, yDatum, yMean, ySum, &
             & zErr

!---------------------------------------------------------------------

!                   BEGINNING OF EXECUTABLE CODE

!   *** KLUDGE ALERT *************************************************
!   Conversion of parameters (constants) to variables should logically
!   have no effect, but in fact helped to suppress some spurious
!   messages from the old IBM vS-FORTRAN compiler:
       mxNode = maxNod
       mxDOF  = 2 * mxNode
       mxEl   = maxEl
       mxFEl  = maxFEl
       mxBn   = maxBN
       mxStar = maxAtP
       mxGeo  = maxGeo
       mxStr  = maxStr
       mxSKS  = maxSKS
       mxGSR  = maxDat
       mxMOR  = maxMOR
       mxAdj  = maxAdj
!   ******************************************************************

!      Mark most pLates as LACKING extensive attached slabs...
       DO i = 1, nPlate
            slab_q(i) = .FALSE.
       END DO
!      ...except for these particular cases:
       slab_q( 8) = .TRUE. !  8 = AU = Australia
       slab_q(14) = .TRUE. ! 14 = CL = Caroline
       slab_q(15) = .TRUE. ! 15 = CO = Cocos
       slab_q(21) = .TRUE. ! 21 = IN = India
       slab_q(22) = .TRUE. ! 22 = JF = Juan de Fuca
       slab_q(34) = .TRUE. ! 34 = NZ = Nazca
       slab_q(37) = .TRUE. ! 37 = PA = Pacific
       slab_q(39) = .TRUE. ! 39 = PS = Philippine Sea
       slab_q(40) = .TRUE. ! 40 = RI = Rivera
       slab_q(46) = .TRUE. ! 46 = SS = Solomon Sea
!      N.B. This array is used when scoring upper-mantle seismic-
!           anisotropy data, because in the Earth5 set of models,
!           plates with slab_q(iPlate) = .TRUE. {those with
!           extensive attached slabs} are assumed to have
!           minimal and/or unresolvable basal shear tractions,
!           and are excluded from the scoring process.
!           That is, if total basal-strength traction for the
!           Pacific plate is toward the Northwest, this might
!           be primarily due to the attached slabs, and it does
!           not preclude weak, distributed basal shear tractions
!           that might point in a completely different direction!


!   Write header on output file:

       IF(Verbose) WRITE(iUnitVerb, 1)
    1  FORMAT ( / &
     &' =============================================================='/ &
     &' I                    PROGRAM -OrbScore2-                     I'/ &
     &' I     A spherical-Earth, thin-plate program for scoring      I'/ &
     &' I      time-averaged (anelastic) deformation solutions       I'/ &
     &' I       computed by finite element program -Shells-.         I'/ &
     &' I                                                            I'/ &
     &' I                           by                               I'/ &
     &' I                       Peter Bird                           I'/ &
     &' I    Department of Earth, Planetary, and Space Sciences      I'/ &
     &' I                University of California                    I'/ &
     &' I               Los Angeles, CA 90095-1567                   I'/ &
     &' I               Version of 14 February 2021                  I'/ &
     &' ==============================================================')

       wedge = ABS(90.0D0 - ABS(dipMax)) * radians_per_degree ! converted to radians away from vertical

       slide = subDip * radians_per_degree ! converted to radians

! ---------------------------------------------------------------------
!      Preview .FEG file to determine array sizes:

       IF(Verbose) WRITE(iUnitVerb, 101) iUnitG
  101  FORMAT (/' Attempting to read finite element grid from unit', I3/)
       READ (iUnitG, * , IOSTAT = ios)
       IF (ios /= 0) THEN
	     write(ErrorMsg,'(A)') "File not found, or file is empty, or file is too short."
		 call FatalError(ErrorMsg,ThID)
       END IF
       READ (iUnitG, * , IOSTAT = ios) numNod
       IF (ios /= 0) THEN
	     write(ErrorMsg,'(A)') "File not found, or file is empty, or file is too short."
		 call FatalError(ErrorMsg,ThID)
       END IF
       IF (numNod > mxNode) THEN
	     write(ErrorMsg,'(A,I0,A)') "Increase PARAMETER maxNode to at least ", numNod," and recompile OrbScore2."
		 call FatalError(ErrorMsg,ThID)
       END IF
       DO 102 i = 1, numNod
            READ (iUnitG, * , IOSTAT = ios)
            IF (ios /= 0) THEN
			  write(ErrorMsg,'(A)') "File not found, or file is empty, or file is too short."
			  call FatalError(ErrorMsg,ThID)
            END IF
  102  CONTINUE
       READ (iUnitG, * , IOSTAT = ios) numEl
       IF (ios /= 0) THEN
		  write(ErrorMsg,'(A)') "File not found, or file is empty, or file is too short."
		  call FatalError(ErrorMsg,ThID)
       END IF
       IF (numeL > mxEl) THEN
	     write(ErrorMsg,'(A,I0,A)') "OrbScore ERROR: Increase PARAMETER maxEl to at least ", numEl," and recompile OrbScore2."
		 call FatalError(ErrorMsg,ThID)
       END IF
       !Initialize survey to find LRn = MAX(continuum_LRi(1:mxEl), fault_LRi(1:MXFel)
       LRn = 0 ! until incremented below...
       DO 103 i = 1, numEl
            READ (iUnitG, "(A)", IOSTAT = ios) longer_line
            IF (ios /= 0) THEN
			  write(ErrorMsg,'(A)') "File not found, or file is empty, or file is too short."
			  call FatalError(ErrorMsg,ThID)
            END IF
            CALL Extract_LRi (longer_line, &     ! input
                            & LRi, shorter_line) ! output
            LRn = MAX(LRn, LRi)
  103  CONTINUE
       nFl = 0
       READ (iUnitG, * , IOSTAT = ios) n
       IF (ios == 0) nFl = n
       nFl = MAX(nFl, 0)
       IF (nFl > mxFEl) THEN
	     write(ErrorMsg,'(A,I0,A)') "OrbScore ERROR: Increase PARAMETER maxFEl to at least ",nFl," and recompile OrbScore2."
		 call FatalError(ErrorMsg,ThID)
       END IF
       DO 105 i = 1, nFl
            READ (iUnitG, "(A)", IOSTAT = ios) longer_line
            IF (ios /= 0) THEN
			  write(ErrorMsg,'(A)') "OrbScore ERROR: FEG File not found, file is empty, or file is too short."
			  call FatalError(ErrorMsg,ThID)
            END IF
            CALL Extract_LRi (longer_line, &     ! input
                            & LRi, shorter_line) ! output
            LRn = MAX(LRn, LRi)
  105  CONTINUE
       REWIND (UNIT = iUnitG) ! to prepare for CALL GetNet, below...

!  DIMENSIONs using newly-determined size variable LRn:
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
! ---------------------------------------------------------------------

!   Input finite element grid and data values at node points:
       CALL GetNet (iUnitG, iUnitVerb, &            ! input
     &              mxDOF, mxEl, mxFEl, mxNode, &
     &              brief, continuum_LRi, cooling_curvature, &  ! output
     &              density_anomaly, &
     &              dQdTdA, elev, fault_LRi, fDip, &
     &              nFakeN, nFl, nodeF, nodes, nRealN, &
     &              numEl, numNod, n1000, offMax, offset, &
     &              title1, tLNode, xNode, yNode, zMNode, &
     &              checkE, checkF, checkN)      ! work

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
	         write(ErrorMsg,'(A)') "OrbScore ERROR: LR File not found, file is empty, or file is too short."
		     call FatalError(ErrorMsg,ThID)
           END IF
           collect_LRs: DO
               READ (iUnitLR, *, IOSTAT = ios) i
               IF (ios /= 0) EXIT collect_LRs ! at EOF, probably
               IF ((i < 1).OR.(i > LRn)) THEN
				 write(ErrorMsg,'(A,I0,A,I0,A/,A)') "OrbScore ERROR: LR# ",i," is outside the legal range of (1:",LRn,").", &
						&  "To make it legal, some element in the .feg file must use this (or higher) LR#."
				 call FatalError(ErrorMsg,ThID)
               END IF
               BACKSPACE(iUnitLR)
               READ (iUnitLR, *, IOSTAT = ios) i, LR_set_fFric(i), LR_set_cFric(i), LR_set_Biot(i), LR_set_Byerly(i), &
                                            &     LR_set_aCreep(1:2, i), LR_set_bCreep(1:2, i), LR_set_cCreep(1:2, i), LR_set_dCreep(1:2, i), &
                                            &     LR_set_eCreep(i)
               IF (ios == 0) THEN
                   LR_is_defined(i) = .TRUE.
               ELSE
				 write(ErrorMsg,'(A,I0)') "ERROR while trying to read 13 REAL*8 values that make up LR# ",i
				 call FatalError(ErrorMsg,ThID)
               END IF
           END DO collect_LRs
           CLOSE (iUnitLR)
           !Now, "stress-test" the continuum elements to be sure that each has a defined rheology:
           DO j = 1, numEl
               i = continuum_LRi(j)
               IF (.NOT.LR_is_defined(i)) THEN
				 write(ErrorMsg,'(A,I0,A,I0,A)') "OrbScore ERROR: Continuum element ", j," uses LR# ", i," which has NOT been defined!"
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
					 write(ErrorMsg,'(A,I0,A,I0,A)') "OrbScore ERROR: Fault element ", j," uses LR# ", i," which has NOT been defined!"
					 call FatalError(ErrorMsg,ThID)
                   ELSE
                       LR_is_used(i) = .TRUE.
                   END IF
               END DO
           END IF
           !Write a report to the log-file, to provide a record of the LRs used:
           IF(Verbose) WRITE(iUnitVerb, *)
           IF(Verbose) WRITE(iUnitVerb, "(' ===========================================================================================================================')")
           IF(Verbose) WRITE(iUnitVerb, "(' Table of alternative Lithospheric Rheologies defined and used:')")
           IF(Verbose) WRITE(iUnitVerb, "('      LR# fFric cFric  Biot Byerly aCreep(1) aCreep(2) bCreep(1) bCreep(2) cCreep(1) cCreep(2) dCreep(1) dCreep(2)    eCreep')")
           DO i = 0, LRn
               IF (LR_is_defined(i).AND.LR_is_used(i)) THEN
                   IF(Verbose) WRITE(iUnitVerb, "(' ', I8, F6.3, F6.3, F6.3, F7.3, ES10.2, ES10.2, F10.0, F10.0, F10.4, F10.4, ES10.2, ES10.2, F10.5)") &
                         & i, LR_set_fFric(i), LR_set_cFric(i), LR_set_Biot(i), LR_set_Byerly(i), &
                         & LR_set_aCreep(1:2, i), LR_set_bCreep(1:2, i), LR_set_cCreep(1:2, i), LR_set_dCreep(1:2, i), &
                         & LR_set_eCreep(i)
               END IF
           END DO
           IF(Verbose) WRITE(iUnitVerb, "('===========================================================================================================================')")
           IF(Verbose) WRITE(iUnitVerb, *)
       END IF ! LRn > 0
! ---------------------------------------------------------------------

!   Check topology, and compute geometric properties:

       log_strike_adjustments = .FALSE.
       skipBC = .TRUE.
       CALL Square (brief, fDip, iUnitVerb, &
     &              log_strike_adjustments, &
     &              mxBn, mxEl, mxFEl, mxNode, &
     &              mxStar, nFl, nodeF, nodes, &
     &              numEl, numNod, skipBC, radius, wedge, &     ! INTENT(IN)
     &              xNode, yNode, &                             ! INTENT(INOUT)
     &              area, detJ, &
     &              dxs, dys, dxsp, dysp, edgeFS, &
     &              edgeTs, fLen, fPFlt, fpsfer, &
     &              fArg, nCond, nodCon, sita, &        ! INTENT(OUT)
     &              checkN, list)                       ! WORKSPACE
       sphere = (nCond == 0)

!   Input velocities at node points:

      IF(Verbose) WRITE(iUnitVerb, 10) iUnitV
   10 FORMAT(/' ATTEMPTING TO READ VELOCITIES OF NODES ON UNIT', I3/)
      CALL OldVel (iUnitVerb, iUnitV, mxNode, numNod, &  ! INTENT(IN)
     &             haveNV, title1, title2, title3, v) ! INTENT(OUT)
      IF (.NOT.haveNV) THEN
	     write(ErrorMsg,'(A,I2)') "OrbScore ERROR: No nodal velocities found on unit ",iUnitV
		 call FatalError(ErrorMsg,ThID)
      END IF
!----------------------------------------------------------------

!   COMPUTE RMS VALUE OF VELOCITIES OF NODES
!   (Perhaps RMS value of surface integral would be more
!     elegant, but this is much easier and almost the same.)

!   First, it is necessary to remove any net rotation:

       numAdj = numNod
       IF (mxAdj < numAdj) THEN
	     write(ErrorMsg,'(A/,A/,A/,A/,A/,A/,A,I7,A)') "In order to remove net rotation from model velocities", &
     &              " (so as to accurately assess RMS velocity),",      &
     &              " it is necessary to send a dummy set of (0, 0)",   &
     &              " benchmark velocities to subprogram Adjust.",      &
     &              " This in turn requires that parameter maxAdj",     &
     &              " must be at least equal to the number of",         &
     &              " nodes, which currently is: ", numNod, ". Please fix!"
		 call FatalError(ErrorMsg,ThID)
       END IF
       DO 250 i = 1, numNod
            nLat(i) = 90.0D0 - xNode(i) * oezOPi
            eLon(i) = yNode(i) * oezOPi
            predic(1, i) = v(1, i)
            predic(2, i) = v(2, i)
            data(1, i) = 0.0D0
            data(2, i) = 0.0D0
  250  CONTINUE
       IF (floats) THEN
            ALLOCATE ( weights(numNod) )
            weights = 1.0D0 ! whole array
            IF(Verbose) WRITE(iUnitVerb, 251)
  251       FORMAT (/' Calling Adjust with equal weight on each node:')
            CALL Adjust (data, eLon, iUnitVerb, mxAdj, numNod, nLat, & ! INTENT(IN)
     &                   radius, weights, &
     &                   predic, &                                  ! INTENT(INOUT)
     &                   eLonP, nLatP, rate)                        ! INTENT(OUT)
            DEALLOCATE (weights)
       END IF
       sum = 0.0D0
       DO 290 i = 1, numNod
            sum = sum + predic(1, i)**2 + predic(2, i)**2
  290  CONTINUE
       rmsv = DSQRT(sum / numNod)
       t2 = 1000.0D0 * secPYr * rmsv
       IF (floats) THEN
            rate = -rate
            t = ABS(rate) * 1000.0D0 * secPYr * radius
            IF(Verbose) WRITE(iUnitVerb, 295) rate, t, eLonP, nLatP, rmsv, t2
  295       FORMAT (/' After removal of net rotation of ', ES10.2, &
     &               ' radians/second (', F5.1, ' mm/yr at equator)' &
     &              /' about a pole at ', F7.2, ' degrees East, ', F6.2, &
     &               ' degrees North,' &
     &              /' the RMS nodal velocity for this model is ', &
     &                 ES10.2, &
     &               ' (', F7.2, ' mm/yr).')
       ELSE
            IF(Verbose) WRITE(iUnitVerb, 296) rmsv, t2
  296       FORMAT (/' The RMS nodal velocity for this model is ', &
     &                 ES10.2, &
     &               ' (', F7.2, ' mm/yr).')
       END IF

!-------------------------------------------------------------------

!           ********* GEODESY SCORING SECTION **********

!   Read geodetic data for scoring:

!  (Dislocation-in-elastic-halfspace corrections for
!   temporary locking of brittle parts of faults within
!   the model were added in Summer 2000 by Zhen Liu and Peter Bird.)
!   Note that this does NOT include corrections for
!   temporary locking of any subduction zones which might lie
!   along the boundary of a model, UNLESS fault elements
!   are used to represent those subduction shear zones.)

       CALL GetGEO (iUnitB, iUnitVerb, iUnitY, iUnitC, mxGeo, &
     &              pltGeo, pltCos, pltInt, iUnitI, &                 ! INTENT(IN)
     &              geoTag, geoPhi, geoThe, &                ! INTENT(OUT)
     &              geoVel, geoSig, geoAzi, &
     &              gpsFMT, numGeo)

!   Flag "explod" indicates that number of benchmarks
!   is equal to number of nodes.  This happens only when
!   the "geodesy" (.gps) file was created by Explode3, and
!   has "benchmark" locations for every node in the original
!   .feg file, but slightly displaced away from nodes.

!   In this case, processing is identical, but an extra
!   output file is produced to facilitate plotting
!   the smoothed (interseismic) velocities.
       explod = (numGeo == numNod)

       IF (numGeo <= 0) THEN
            geodes = 0.0D0
            CLOSE (UNIT = iUnitY, DISP = "DELETE")
            GO TO 400
       END IF
       IF(Verbose) WRITE(iUnitVerb, 301) numGeo, iUnitB
  301  FORMAT (/' ', I6, ' Geodetic velocity data were read from unit ', I2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Find weights for geodetic benchmarks, based on associated nearest-neighbor surface areas:

      IF (DIs_Initialization_Needed()) THEN
          subdivision = 4
          IF(Verbose) WRITE(iUnitVerb, "(/' Initializing uniform global grid at subdivision ', &
     &           I1,' for area-weighting of data...')") subdivision
          CALL DInitialize_Weighting (subdivision)
      END IF

      IF(Verbose) WRITE(iUnitVerb, *)
      IF(Verbose) WRITE(iUnitVerb, "(' Computing area-weights for geodetic benchmarks...')")
      ALLOCATE ( weights(numGeo) )
      CALL DPerform_Weighting (number_of_data = numGeo, &
                             & theta_radians  = geoThe, &
                             & phi_radians    = geoPhi, & ! INTENT(IN)
                             & weights        = weights)  ! INTENT(OUT)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Find predicted velocities at benchmarks
!    (A) Steady-state components:

      IF (explod) THEN
           IF(Verbose) WRITE(iUnitVerb, 302)
  302      FORMAT(/' Calling FindIt to determine element containing' &
     &            /' each displaced node; be patient!....')
      END IF

       i = 0
  303  i = i + 1
            theta = geoThe(i)
            phi = geoPhi(i)
            CALL FindIt (mxEl, mxNode, nodes, numEl, xNode, yNode, & ! INTENT(IN)
     &                   theta, phi, &
     &                   iEle, s1, s2, s3)                           ! INTENT(OUT)
            IF (iEle == 0) THEN
!                Bemchmark is outside of grid; delete it.
                 IF(Verbose) WRITE(iUnitVerb, 305) geoTag(i)
  305            FORMAT (' Benchmark ', A20, ' was outside the grid: deleted.')
                 DO 310 j = i, numGeo - 1
                      geoTag(j) = geoTag(j + 1)
                      geoPhi(j) = geoPhi(j + 1)
                      geoThe(j) = geoThe(j + 1)
                      geoVel(j) = geoVel(j + 1)
                      geoSig(j) = geoSig(j + 1)
                      geoAzi(j) = geoAzi(j + 1)
  310            CONTINUE
                 numGeo = numGeo - 1
                 i = i - 1
            ELSE
                 theLat = 90.0D0 - theta * oezOPi
                 theLon = phi * oezOPi
                 IF (theLon > +180.0D0) theLon = theLon - 360.0D0
                 IF (theLon > +180.0D0) theLon = theLon - 360.0D0
                 IF (theLon < -180.0D0) theLon = theLon + 360.0D0
                 IF (theLon < -180.0D0) theLon = theLon + 360.0D0
                 IF (.NOT.explod) THEN
                     showis = .FALSE.
                     IF (showis) THEN
                         IF(Verbose) WRITE(iUnitVerb, 315)geoTag(i), theLon, theLat, iEle, s1, s2, s3
  315                    FORMAT(' ', A20,' (', F7.2, 'E, ', F6.2, 'N) = #', I5, ':(', F5.3, ',', F5.3, ',', F5.3,')')
                     END IF
                 END IF
!                Interpolate v at location (s1, s2, s3) in element #iEle:
                 CALL Projec (iEle, &
     &                        mxEl, mxNode, nodes, &
     &                        s1, s2, s3, &
     &                        xNode, yNode, v, &     ! INTENT(IN)
     &                        vTheta, vPhi)          ! INTENT(OUT)
                 predic(1, i) = vTheta
                 predic(2, i) = vPhi
                 longterm_at_GPS(1, i) = vTheta ! This duplicate copy of the long-term velocities
                 longterm_at_GPS(2, i) = vPhi   ! at GPS benchmarks (predicted by Shells) will NOT
                                                ! be modified by CALL Adjust (or in any other way),
                                                ! so that it will remain available IF (pltInt).
                 data(1, i) = geoVel(i) * COS(geoAzi(i))
                 data(2, i) = geoVel(i) * SIN(geoAzi(i))
                 nLat(i) = 90.0D0 - geoThe(i) * oezOPi
                 eLon(i) = geoPhi(i) * oezOPi
            END IF
!      End of variable-length loop:
       IF (i < numGeo) GO TO 303

       IF (explod.AND.(numGeo < numNod)) THEN
!        Condition (numGeo == numNod) that made "explod" .TRUE. is no longer valid:
            explod = .FALSE.
			write(ErrorMsg,'(A/,A/,A/,A/,A/,A)') "This should not happen OFTEN when benchmark",                         &
				&             " locations are really node locations, slightly displaced by Explode3.",                        &
				&             " Probably the bad benchmark(s) is/are only slightly outside the domain of the F-E model.",     &
				&             " Use OrbWin (or OrbWeaver) to view the Explode3d grid and see why some nodes/benchmarks moved",&
				&             " out of bounds.  Then adjust the corresponding benchmark coordinates in the false-benchmark",  &
				&             " data file, and re-run OrbScore2."
			call FatalError(ErrorMsg,ThID)
       END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Find predicted velocities at benchmarks;
!    (B) Coseismic displacement rate due to brittle
!        parts of faults within the model domain.

!   Find thicknesses of crust and mantle-lithosphere:

       CALL Interp (zMNode, mxEl, mxNode, nodes, numEl, & ! INTENT(IN)
     &              zMoho)                                ! INTENT(OUT)

       CALL Interp (tLNode, mxEl, mxNode, nodes, numEl, & ! INTENT(IN)
     &              tLInt)                                ! INTENT(OUT)

!   Compute tactical values of limits on viscosity, and weights for
!   imposition of constraints in linear systems:

       CALL Limits (area, detJ, iUnitVerb, mxEl, numEl, &
     &              OKDelV, radius, refStr, sphere, tLInt, &
     &              trHMax, zMoho, &                       ! INTENT(IN)
     &              constr, etaMax, fMuMax, visMax)        ! INTENT(OUT)

!   Locate the brittle/ductile transition in each fault:

!      Begin with very crude initialization:
       DO i = 1, nFl
            zTranF(1, i) = MAX(zMNode(nodeF(1, i)) / 2.0D0, oneKm)
            zTranF(2, i) = MAX(tLNode(nodeF(1, i)) / 2.0D0, oneKm)
       END DO
!      Then search for actual brittle/ductile transition, iteratively:
       DO iMohr = 1, 3
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
       END DO
!    NEW FEATURE (2018.11.26):
!    In all subduction zones (oceanic or continental),
!    override the brittle/ductile transition(s)
!    determined (above) by Mohr, and replace them with a single 26-km-thick
!    frictional layer, following Oleskevich et al. [1999, J. Geophys. Res.]
!    and Bird & Kagan [2004, Bull. Seismol. Soc. Am.]:
     DO i = 1, nFl
         dip = 0.5D0 * (fDip(1, i) + fDip(2, i)) ! averaging the dips at the 2 ends of fault trace.
         IF ((dip <= slide).OR.((Pi - dip) <= slide)) THEN ! N.B. slide is limiting (maximum) value for SUBs, from PARAMETER subDip, after conversion to radians.
             !This fault is a SUB:
             zTranF(1, i) = 26.0D3 ! 26 km = average depth range of seismogenic patches in SUBs, from Oleskevich et al. [1999, JGR].
             zTranF(2, i) =  0.0D3 ! (Not adding any other frictional layer at the top of the mantle, which is probably serpentinized.)
         END IF
     END DO

!    Compute benchmark velocities due to elastic dislocations
!    from earthquakes in the brittle part of each fault, at the slip-rate
!    predicted by the -Shells- finite-element model:

      IF (explod) THEN
           IF(Verbose) WRITE(iUnitVerb, 340)
  340      FORMAT(/' Calling Coseis to evaluate long-term-average coseismic' &
     &            /' velocity of each node; be patient!....')
      ELSE
           IF(Verbose) WRITE(iUnitVerb, 341)
  341      FORMAT(/' Calling Coseis to evaluate long-term-average coseismic' &
     &            /' velocity of each benchmark....')
      END IF

      CALL DCoseis (fArg, fDip, geoThe, geoPhi, mxNode, &
     &              mxFEl, numGeo, nFl, nodeF, radius, &
     &              v, wedge, xNode, yNode, zMNode, zTranF, &     ! INTENT(IN)
     &              geoUTh, geoUPh)                               ! INTENT(OUT)

!    For a typical year with no earthquakes,
!    combine the predicted benchmark velocities
!    by subtracting the coseismic part:

      DO 350 i = 1, numGeo
           predic(1, i) = predic(1, i) - geoUTh(i)
           predic(2, i) = predic(2, i) - geoUPh(i)
  350 CONTINUE

       IF (explod) THEN
            IF(Verbose) WRITE(iUnitVerb, 351)
  351       FORMAT(//' Creating file with smoothed (interseismic)' &
     &             /' velocities in v_____.out format, to be used' &
     &             /' as input to FiniteMap for visualizing them.')
            IF(Verbose) WRITE(iUnitVerb, 352) iUnitX
  352       FORMAT(/' Attempting to create output v____.out file for' &
     &             /' smoothed (interseismic) velocities at benchmarks' &
     &             /' (which were slightly displaced from nodes by' &
     &             /' Explode3) using unit ',I3/)
            OPEN (UNIT = iUnitX, FILE = '')
            line = "Explode3 of : "//TRIM(title1(1:66))
            WRITE(iUnitX, "(A)") TRIM(line)
            line = "Interseismic: "//TRIM(title2(1:66))
            WRITE(iUnitX, "(A)") TRIM(line)
            line = "velocities  : "//TRIM(title3(1:66))
            WRITE(iUnitX, "(A)") TRIM(line)
            WRITE(iUnitX, "(4ES18.9)")((predic(i, j), i = 1, 2), j = 1, numGeo)
!           Note: numGeo == numNod, or we wouldn't be here.
            CLOSE(iUnitX)

			write(ErrorMsg,'(A/,A/,A/,A/,A/,A)') "Smoothed velocities have been written.", &
			&             "This run of OrbScore2 will now end, as",                           &
			&             "it was not intended for actual scoring,",                          &
			&             "but just to evaluate interseismic",                                &
			&             "velocities at ""benchmarks"" which are",                           &
			&             "really slightly-displaced nodes."
			call FatalError(ErrorMsg,ThID)
       END IF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Adjust predictions by rigid rotation to best-fit geodesy:

       IF(Verbose) WRITE(iUnitVerb, 362)
  362  FORMAT (/' Calling Adjust with area-weight on each benchmark' &
     &         /' to minimize global area-integral of velocity-error-squared:')
       IF (floats) THEN
            CALL Adjust (data, eLon, iUnitVerb, mxAdj, numGeo, nLat, &
     &                   radius, weights, &     ! INTENT(IN)
     &                   predic, &              ! INTENT(INOUT)
     &                   eLonP, nLatP, rate)    ! INTENT(OUT)
            t = ABS(rate) * 1000.0D0 * secPYr * radius
            IF(Verbose) WRITE(iUnitVerb, 370) rate, t, eLonP, nLatP
  370       FORMAT (/' Model predictions adjusted by rotation of ', ES10.2, &
     &               ' radians/second (', F7.2,' mm/yr at equator)' &
     &              /' about a pole at ', F7.2, ' degrees East, ', F6.2, &
     &               ' degrees North,' &
     &              /' before comparison to geodetic data.')
       END IF

       IF(Verbose) WRITE(iUnitVerb, 380)
  380  FORMAT (// &
     &' Geodetic Benchmark Velocities Versus Model Predictions:'// &
     &' Benchmark            Velocity (mm/yr) Azimuth   Model.v (mm/yr)' &
     &,' Model.Az     Error (mm/yr) (Sigmas) area-weight'/ &
     &' ---------            --------  -----  -------   -------  -----', &
     &'  --------     -----  -----   ------  ----------')
  381  FORMAT (' ', A20, ES9.2, ' (', F5.1, ') ', F7.1, ES10.2, &
              &' (', F5.1, ') ', F8.1, ES10.2, ' (', F5.1, ')  (', F5.1, ') ', F10.5)

       sum0s = 0.0D0
       sum1 = 0.0D0
       sum1s = 0.0D0
       sum2 = 0.0D0
       sum2s = 0.0D0
       sumN = 0.0D0
       sumNs = 0.0D0
       DO 390 i = 1, numGeo
            gVMma = geoVel(i) * secPYr * 1000.0D0
            azim = 180.0D0 - geoAzi(i) * oezOPi
            IF (azim < 0.0D0) azim = azim + 360.0D0
            IF (azim > 360.0D0) azim = azim - 360.0D0
            vM = DSQRT(predic(1, i)**2 + predic(2, i)**2)
            vMMma = vM * secPYr * 1000.0D0
            pAzim = 180.0D0 - ATan2F(predic(2, i), predic(1, i)) * oezOPi
            IF (pAzim < 0.0D0) pAzim = pAzim + 360.0D0
            IF (pAzim > 360.0D0) pAzim = pAzim - 360.0D0
            bad = DSQRT((predic(1, i) - data(1, i))**2 + &
     &                  (predic(2, i) - data(2, i))**2)
            badMma = bad * secPYr * 1000.0D0
            badSig = bad / geoSig(i)
            IF(Verbose) WRITE(iUnitVerb, 381) geoTag(i), geoVel(i), gVMma, azim, &
     &                          vM, vMMma, pAzim, bad, badMma, badSig, &
     &                          weights(i)
            IF (pltGeo) THEN
                 eLonDe = geoPhi(i) / piO180
                 nLatDe = 90.0D0 - (geoThe(i) / piO180)
                 nLatDe = MIN(90.0D0, MAX(-90.0D0, nLatDe))
                 vEMmpa = +(predic(2, i) - data(2, i)) * 1000.0D0 * secPYr ! comparing (reference-frame-adjusted Shells predictions
                 vNMmpa = -(predic(1, i) - data(1, i)) * 1000.0D0 * secPYr ! of interseismic velocities) at benchmarks to (data)
                 vESigm = 0.1D0
                 vNSigm = 0.1D0
                 correl = 0.0D0
                 frame = "[same]"
                 tag = geoTag(i)
                 WRITE (iUnitY, gpsFMT) eLonDe, nLatDe, vEMmpa, vNMmpa, &
     &                                  vESigm, vNSigm, correl, &
     &                                  TRIM(frame), TRIM(tag)
            END IF ! pltGeo
            IF (pltCos) THEN ! WRITE out the Coseismic part of the Shells-model predicted geodetic velocities:
                 eLonDe = geoPhi(i) / piO180
                 nLatDe = 90.0D0 - (geoThe(i) / piO180)
                 nLatDe = MIN(90.0D0, MAX(-90.0D0, nLatDe))
                 vEMmpa = +geoUPh(i) * 1000.0D0 * secPYr
                 vNMmpa = -geoUTh(i) * 1000.0D0 * secPYr
                 vESigm = 0.1D0
                 vNSigm = 0.1D0
                 correl = 0.0D0
                 frame = "[same]"
                 tag = geoTag(i)
                 WRITE (iUnitC, gpsFMT) eLonDe, nLatDe, vEMmpa, vNMmpa, &
     &                                  vESigm, vNSigm, correl, &
     &                                  TRIM(frame), TRIM(tag)
            END IF ! pltCos
            IF (pltInt) THEN ! WRITE out the (NON-frame-adjusted) Interseismic part of the Shells-model predicted geodetic velocities:
                 eLonDe = geoPhi(i) / piO180
                 nLatDe = 90.0D0 - (geoThe(i) / piO180)
                 nLatDe = MIN(90.0D0, MAX(-90.0D0, nLatDe))
                 vEMmpa = +(longterm_at_GPS(2, i) - geoUPh(i)) * 1000.0D0 * secPYr ! long-term Shells velocity minus mean co-seismic velocity,
                 vNMmpa = -(longterm_at_GPS(1, i) - geoUTh(i)) * 1000.0D0 * secPYr ! both evaluated at geodetic benchmarks
                 vESigm = 0.1D0
                 vNSigm = 0.1D0
                 correl = 0.0D0
                 frame = "[same]"
                 tag = geoTag(i)
                 WRITE (iUnitI, gpsFMT) eLonDe, nLatDe, vEMmpa, vNMmpa, &
     &                                  vESigm, vNSigm, correl, &
     &                                  TRIM(frame), TRIM(tag)
            END IF ! pltInt
            IF (badSig > 2.0D0) sum0s = sum0s + weights(i)
            sum1 = sum1 + weights(i) * bad
            sum1s = sum1s + weights(i) * badSig
            sum2 = sum2 + weights(i) * bad**2
            sum2s = sum2s + weights(i) * badSig**2
            sumN = MAX(sumN, bad)
            sumNs = MAX(sumNs, badSig)
  390  CONTINUE
       sum0s = sum0s / numGeo
       sum1 = sum1 / numGeo
       sum1mm = sum1 * 1000.0D0 * secPYr
       sum1s = sum1s / numGeo
       sum2 = SQRT(sum2 / numGeo)
       sum2mm = sum2 * secPYr * 1000.0D0
       sum2s = SQRT(sum2s / numGeo)
       sumNmm = sumN * 1000.0D0 * secPYr
       IF(Verbose) WRITE(iUnitVerb, 397) sum0s, sum1, sum1mm, sum1s, &
     &                    sum2, sum2mm
  397  FORMAT(/ &
     &        /' Summary of geodetic prediction errors:' &
     &        /' Fraction of predictions wrong by over 2-sigma: ', F5.3 &
     &        /' Mean velocity error:  ', ES10.2, ' (', F7.2, ' mm/yr)' &
     &        /' Mean number of sigmas in error: ', F5.2, &
     &        /' RMS velocity error:   ', ES10.2, ' (', F7.2, ' mm/yr)' &
     &)
       IF (floats) THEN
            IF(Verbose) WRITE(iUnitVerb, 398)
  398       FORMAT(' *** Note: Line above was optimized by Adjust. ***')
       END IF
       IF(Verbose) WRITE(iUnitVerb, 399) sum2s, sumN, sumNmm, sumNs
  399  FORMAT(' RMS number of sigmas in error: ', F5.2 &
     &        /' Worst velocity error :', ES10.2, ' (', F7.2, ' mm/yr)' &
     &        /' Largest error, in sigmas: ', F6.2)

!    Select one indicator of error (RMS error, after frame-change):
       geodes = sum2mm

       IF (pltGeo) THEN
            CLOSE(iUnitY)
            IF(Verbose) WRITE(iUnitVerb, "(/' Wrote ERRORS.gps, for plotting geodetic errors with -FiniteMap-.')")
       END IF
       IF (pltCos) THEN
            CLOSE(iUnitC)
            IF(Verbose) WRITE(iUnitVerb, "(/' Wrote COSEISMIC.gps, for plotting model coseismic velocities with -FiniteMap-.')")
       END IF
       IF (pltInt) THEN
            CLOSE(iUnitI)
            IF(Verbose) WRITE(iUnitVerb, "(/' Wrote INTERSEISMIC.gps, for plotting model interseismic velocities with -FiniteMap-.')")
       END IF
       DEALLOCATE (weights)

!           ********* END SPACE-GEODESY SCORING SECTION **********

!-------------------------------------------------------------------

!             ****** STRESS-DIRECTION SCORING SECTION ********

!   Read most-compressive horizontal principal stress azimuths for scoring.

  400 CALL GetSTR (iUnitS, iUnitVerb, mxStr, &        ! INTENT(IN)
     &             strTag, strThe, strPhi, &
     &             strArg, strQua, strReg, numStr) ! INTENT(OUT)
      IF (numStr <= 0) THEN
           stress = 0.0D0
           regime = 0.0D0
           GO TO 1301
      END IF
      IF(Verbose) WRITE(iUnitVerb, 401) numStr, iUnitS
  401 FORMAT (/' ',I6,' Stress direction data were read from unit ', I2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (DIs_Initialization_Needed()) THEN
          subdivision = 4
          IF(Verbose) WRITE(iUnitVerb, "(/' Initializing uniform global grid at subdivision ', &
     &              I1, ' for area-weighting of data...')") subdivision
          CALL DInitialize_Weighting (subdivision)
      END IF

      IF(Verbose) WRITE(iUnitVerb, "(' Computing area-weights for stress direction/regime data...')")
      ALLOCATE ( weights(numStr) )
      weights = 1.0D0 ! whole array
!      CALL DPerform_Weighting (number_of_data = numStr, &
!                             & theta_radians  = strThe, &
!                             & phi_radians    = strPhi, & ! INTENT(IN)
!                             & weights = weights)         ! INTENT(OUT)
!  N. B. Area-weighting turned off here when we switched to use of
!        robust_interpolated_stress_for_OrbScore2, 2007.03.12 in the Earth5 project.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Compute the (theta, phi) and then Cartesian (x, y, z) of integration points:

       DO 450 i = 1, numEl
            DO 420 k = 1, 3
                 CartVs(1, k) = DCOS(yNode(nodes(k, i))) * &
     &                          DSIN(xNode(nodes(k, i)))
                 CartVs(2, k) = DSIN(yNode(nodes(k, i))) * &
     &                          DSIN(xNode(nodes(k, i)))
                 CartVs(3, k) = DCOS(xNode(nodes(k, i)))
  420       CONTINUE
            DO 440 m = 1, 7
                 DO 430 j = 1, 3
                      tempV(j) = 0.0D0
                      DO 425 k = 1, 3
                           tempV(j) = tempV(j) + CartVs(j, k) * points(k, m)
  425                 CONTINUE
  430            CONTINUE
                 CALL Unit(tempV) ! INTENT(INOUT)
                 IF (ABS(tempV(3)) <= 0.5D0) THEN
                      xIP(m, i) = DACOS(tempV(3))
                 ELSE
                      equPar = DSQRT(tempV(1)**2 + tempV(2)**2)
                      xIP(m, i) = ATan2F(equPar, tempV(3))
                 END IF
                 yIP(m, i) = ATan2F(tempV(2), tempV(1))
                 CartR(1, m, i) = DSIN(xIP(m, i)) * DCOS(yIP(m, i))
                 CartR(2, m, i) = DSIN(xIP(m, i)) * DSIN(yIP(m, i))
                 CartR(3, m, i) = DCOS(xIP(m, i))
  440       CONTINUE
  450  CONTINUE

!   Strain-rates at integration points:

       CALL EDot (dxs, dys, &
     &            fpsfer, mxEl, &
     &            mxNode, nodes, numEl, radius, sita, v, & ! INTENT(IN)
     &            eRate)                                   ! INTENT(OUT)

       IF(Verbose) WRITE(iUnitVerb, 460)
  460  FORMAT (// &
     &' Horizontal most-compressive stress directions versus ', &
     &   'model predictions:'/ &
     &/' Datum E.Long. N.Latt. Element Point', &
     & ' Quality Azimuth Model.Az Error', &
     & ' Regime       E1       E2       E3      Err Match? area-weight' &
     &/' ----- ------- ------- ------- ----- ------- ', &
     &  '------- -------- -----', &
     & ' ------ -------- -------- -------- -------- ------ -----------')
  461  FORMAT (' ', A5, 1X, F7.2, 1X, F7.2, 1X, I7, 1X, I5, 1X, F7.0, 1X, &
     &         F7.0, 1X, F8.0, 1X, F5.0, &
     &         5X, A2, 4ES9.1, 1X, A4, 2X, F10.5)
!: Accumulator of bad stress regimes:
       nBadRe = 0     ! simple integer count of bad regimes
       aBadRe = 0.0D0 ! area-weighted error count
       aRegime_denominator = 0.0D0
!: Accumulators for general data (all qualities):
       sum1 = 0.0D0
       sum1d = 0.0D0
       sum2 = 0.0D0
       sum2d = 0.0D0
!: Accumulators for A,B-quality data only:
       n_AB_data = 0
       aSum1 = 0.0D0
       aSum1D = 0.0D0
       aSum2 = 0.0D0
       aSum2D = 0.0D0

!   NOTE: Integration point must be within one-half element-
!         width of datum, or datum is not used.
!         Here, the typical element side is deterMined by
!         averaging together all element sides,
!         and expressed in radians:

       sumSid = 0.0D0
       DO 471 i = 1, numEl
            DO 470 k = 1, 3
                 kp1 = 1 + MOD(k, 3)
                 n1 = nodes(k, i)
                 n2 = nodes(kp1, i)
                 x1 = xNode(n1)
                 x2 = xNode(n2)
                 y1 = yNode(n1)
                 y2 = yNode(n2)
                 rA(1) = DSIN(x1) * DCOS(y1)
                 rA(2) = DSIN(x1) * DSIN(y1)
                 rA(3) = DCOS(x1)
                 rB(1) = DSIN(x2) * DCOS(y2)
                 rB(2) = DSIN(x2) * DSIN(y2)
                 rB(3) = DCOS(x2)
                 dot = rA(1) * rB(1) + rA(2) * rB(2) + rA(3) * rB(3)
                 CALL DCross (rA, rB, & ! INTENT(IN)
     &                        cUvec)    ! INTENT(OUT)
                 croSiz = DSQRT(cUvec(1)**2 + cUvec(2)**2 + cUvec(3)**2)
                 side = ATan2F(croSiz, dot)
                 sumSid = sumSid + side
  470       CONTINUE
  471  CONTINUE
       toler = sumSid / (6.0D0 * numEl)

       numDen = numStr
       DO 490 i = 1, numStr
!           Next stress direction datum:
            theLon = oezOPi * strPhi(i)
            IF (theLon > 180.0D0) theLon = theLon - 360.0D0
            IF (theLon < -180.0D0) theLon = theLon + 360.0D0
            theLat = 90.0D0 - oezOPi * strThe(i)
            rS(1) = DSIN(strThe(i)) * DCOS(strPhi(i))
            rS(2) = DSIN(strThe(i)) * DSIN(strPhi(i))
            rS(3) = DCOS(strThe(i))
!           Find closest integration point in the grid:
            r2min = 9.99D29
            DO 478 m = 1, 7
                 DO 476 j = 1, numEl
                      r2 = (rS(1) - CartR(1, m, j))**2 + &
     &                     (rS(2) - CartR(2, m, j))**2 + &
     &                     (rS(3) - CartR(3, m, j))**2
                      IF (r2 < r2min) THEN
                           iEle = j
                           mIP = m
                           r2min = r2
                      END IF
  476            CONTINUE
  478       CONTINUE
            rMin = DSQRT(r2min)
            IF (rMin > toler) THEN
                 IF(Verbose) WRITE(iUnitVerb, 480) strtag(i), theLon, theLat, toler
  480            FORMAT (' ',A5,1X,F7.2,1X,F7.2, &
     &                   ' Stress datum is more than ', &
     &                     F6.4, ' radians from nearest' &
     &                  ,' integration point: Ignored.')
                 numDen = numDen - 1
            ELSE
!                Get sigma_1H direction from E_1 direction:
                 e11h = eRate(1, mIP, iEle)
                 e22h = eRate(2, mIP, iEle)
                 e12h = eRate(3, mIP, iEle)
                 CALL Prince (e11h, e22h, e12h, &           ! INTENT(IN)
     &                        e1h, e2h, u1x, u1y, u2x, u2y) ! INTENT(OUT)
                 err = -(e11h + e22h)
!               (or, = -(e1h + e2h) )
                 e1 = MIN(e1h, e2h, err)
                 e3 = MAX(e1h, e2h, err)
                 IF ((err > e1).AND.(err < e3)) THEN
                      e2 = err
                 ELSE IF ((e1h > e1).AND.(e1h < e3)) THEN
                      e2 = e1h
                 ELSE
                      e2 = e2h
                 END IF
                 pAzim = 180.0D0 - oezOPi * ATan2F(u1y, u1x)
                 IF (pAzim < 0.0D0) pAzim = pAzim + 180.0D0
                 IF (pAzim < 0.0D0) pAzim = pAzim + 180.0D0
                 IF (pAzim >= 180.0D0) pAzim = pAzim - 180.0D0
                 IF (pAzim >= 180.0D0) pAzim = pAzim - 180.0D0
                 azim = 180.0D0 - oezOPi * strArg(i)
                 IF (azim < 0.0D0) azim = azim + 180.0D0
                 IF (azim < 0.0D0) azim = azim + 180.0D0
                 IF (azim >= 180.0D0) azim = azim - 180.0D0
                 IF (azim >= 180.0D0) azim = azim - 180.0D0
                 bad = ABS(pAzim - azim)
                 IF (bad > 90.0D0) bad = ABS(bad - 180.0D0)
                 IF (bad > 90.0D0) bad = ABS(bad - 180.0D0)
                 IF (bad > 90.0D0) bad = ABS(bad - 180.0D0)
                 IF (bad > 90.0D0) bad = ABS(bad - 180.0D0)
                 sum1 = sum1 + bad * strQua(i) *weights(i)
                 sum1d = sum1d + strQua(i) * weights(i)
                 sum2 = sum2 + strQua(i) * bad**2 * weights(i)
                 sum2d = sum2d + strQua(i) * weights(i)
!                Maintain separate totals for A,B-quality only:
                 IF (strQua(i) > 2.9D0) THEN
                      n_AB_data = n_AB_data + 1
                      aSum1 = aSum1 + bad * weights(i)
                      aSum1D = aSum1D + weights(i)
                      aSum2 = aSum2 + bad**2 * weights(i)
                      aSum2D = aSum2D + weights(i)
                 END IF

!                Check whether stress regime is correct:
                 reg = strreg(i)
                 outcom = 'OK? '
                 IF (reg == 'NF') THEN
!                     Normal faulting;
!                     e1 should be vertical
                      IF (e1 /= err) THEN
                           nBadRe = nBadRe + 1
                           aBadRe = aBadRe + weights(i)
                           outcom = 'BAD '
                      ELSE
                           outcom = 'GOOD'
                      END IF
                 ELSE IF (reg == 'SS') THEN
!                     Strike-slip faulting;
!                     e2 should be vertical
                      IF (e2 /= err) THEN
                           nBadRe = nBadRe + 1
                           aBadRe = aBadRe + weights(i)
                           outcom = 'BAD '
                      ELSE
                           outcom = 'GOOD'
                      END IF
                 ELSE IF (reg == 'TF') THEN
!                     Thrust-faulting;
!                     e3 should be vertical
                      IF (e3 /= err) THEN
                           nBadRe = nBadRe + 1
                           aBadRe = aBadRe + weights(i)
                           outcom = 'BAD '
                      ELSE
                           outcom = 'GOOD'
                      END IF
                 ELSE IF ((reg == 'ST').OR.(reg == 'TS')) THEN
!                     Transpression;
!                     err should be positive, and
!                     the two horizontal principal stresses
!                     should have opposite signs.
                      IF ((err <= 0.0D0).OR.((e1h * e2h) >= 0.0D0)) THEN
                           nBadRe = nBadRe + 1
                           aBadRe = aBadRe + weights(i)
                           outcom = 'BAD '
                      ELSE
                           outcom = 'GOOD'
                      END IF
                 ELSE IF ((reg == 'NS').OR.(reg == 'SN')) THEN
!                     Transtension;
!                     err should be negative, and
!                     the two horizontal principal stresses
!                     should have opposite signs.
                      IF ((err >= 0.0D0).OR.((e1h * e2h) >= 0.0D0)) THEN
                           nBadRe = nBadRe + 1
                           aBadRe = aBadRe + weights(i)
                           outcom = 'BAD '
                      ELSE
                           outcom = 'GOOD'
                      END IF
                 END IF
                 aRegime_denominator = aRegime_denominator + weights(i)
                 IF(Verbose) WRITE(iUnitVerb, 461) strtag(i), theLon, theLat, iEle, mIP, &
     &                               strQua(i), azim, pAzim, bad, &
     &                               strreg(i), e1, e2, e3, ERR, outcom, weights(i)
            END IF
  490  CONTINUE

       IF (sum1d > 0.0D0) THEN
            sum1 = sum1 / sum1d
       ELSE
            sum1 = 0.0D0
       END IF
       IF (sum2d > 0.0D0) THEN
            sum2 = DSQRT(sum2 / sum2d)
       ELSE
            sum2 = 0.0D0
       END IF
       IF(Verbose) WRITE(iUnitVerb, 496) numDen, sum1, sum2
  496  FORMAT (/' Summary of Stress-Direction Errors:' &
     &         /' Number of stress directions used for scoring: ', I6, '.' &
     &         /' Area- & quality-weighted MEAN error: ', F5.2, ' degrees.' &
     &         /' Area- & quality-weighted RMS  error: ', F5.2, ' degrees.')
       IF ((aSum1D > 0.0D0).AND.(aSum2D > 0.0D0)) THEN
            aSum1 = aSum1 / aSum1D
            aSum2 = DSQRT(aSum2 / aSum2D)
            IF(Verbose) WRITE(iUnitVerb, 497) n_AB_data, aSum1, aSum2
  497  FORMAT (/' OR, using only the ',I6,' A,B-quality data:' &
     &         /' Area-weighted MEAN error: ', F5.2, ' degrees.' &
     &         /' Area-weighted RMS  error: ', F5.2, ' degrees.')
       END IF

       IF (numDen > 0.0D0) THEN
            perBad = (100.0D0 * nBadRe) / (1.0D0 * numDen)
       ELSE
            perBad = 0.0D0
       END IF
       IF(Verbose) WRITE(iUnitVerb, 498) nBadRe, numDen, perBad
  498  FORMAT (' Regime:', I6, ' bad stress regimes out of ', I6, &
     &         ' data = ', F6.2, ' % bad.')

       IF (aRegime_denominator > 0.0D0) THEN
            perBad = (100.0D0 * aBadRe) / aRegime_denominator
       ELSE
            perBad = 0.0D0
       END IF
       IF(Verbose) WRITE(iUnitVerb, 499) aBadRe, aRegime_denominator, perBad
  499  FORMAT (' Regime:', F10.3, ' area-weighted bad stress regimes out of ', F10.3, &
     &         ' area-weighted tests = ', F6.2, ' % bad.')

!   Arbitrarily choose measure of error (mean error, in degrees):
       stress = sum1

       DEALLOCATE ( weights )

!           ********* END STRESS SCORING SECTION **********

!-------------------------------------------------------------------

!           ********* SLIP-RATE SCORING SECTION ***********

!  Read geological slip rates on faults which are represented by
!  fault elements in the .FEG file.
!  Score both components: strike-slip, and relative-vertical (throw) data
!  where available.
!  This part of the code was modified from the -PLATES- scoring program,
!  -Score_AK.for.
!  The data format is given in: Slip_rate_format.txt.

!   Read test dataset of geologic slip-rate limits for faults:

 1301  CALL GetSLF(iUnitD, iUnitVerb, mxGSR, &                    ! INTENT(IN)
     &             fName, rLat, delVZ, fltLon, fltLat, nFData) ! INTENT(OUT)

       IF (nFData <= 0) THEN
           slpErr = 0.0D0
           GO TO 501
       END IF

!   Compute slip rates (in mm/year) at specified points
!   and assign scores and probabilities:

       IF(Verbose) WRITE(iUnitVerb, 1310)
 1310  FORMAT (/' Geologic sliprate predictions and their errors (in mm/year):')
       IF(Verbose) WRITE(iUnitVerb, 1320)
 1320  FORMAT ( &
     & '                                Minimum Maximum   Model' &
     &,' Dextral        Minimum Maximum   Model   Throw       ' &
     &/' Fault name (first 30 bytes)    Dextral Dextral Dextral' &
     &,'   Error Chance   Throw   Throw   Throw   Error Chance' &
     &/' ------------------------------ ------- ------- -------' &
     &,' ------- ------ ------- ------- ------- ------- ------')

       shrink = 1.0D0
       nGSDat = 0
       nB = 0
       sumB = 0.0D0
       sumB2 = 0.0D0
       bigB = 0.0D0

       DO 1390 j = 1, nFData

!        Find fault element closest to (fltLon, fltLat):

!           xDatum is Theta, S from North pole, in radians:
            xDatum = (90.0D0 - fltLat(j)) * piO180
!           yDatum is Phi, E from Greenwich meridan, in radians:
            yDatum = fltLon(j) * piO180
            IF (nFl == 0) THEN
			  write(ErrorMsg,'(A,I8)') "ERR0R: NO FAULT ELEMENTS IN THIS .FEG FILE TO MATCH FAULT SLIP-RATE DATUM #",j
			  call FatalError(ErrorMsg,ThID)
            END IF
            i = 0
            s = 0.0D0
            r2min = 9.99D+29
            DO 1340 k = 1, nFl
                 n1 = nodeF(1, k)
                 n2 = nodeF(2, k)
                 n3 = nodeF(3, k)
                 n4 = nodeF(4, k)
                 x1 = xNode(n1)
                 x2 = xNode(n2)
                 x3 = xNode(n3)
                 x4 = xNode(n4)
                 y1 = yNode(n1)
                 y2 = yNode(n2)
                 y3 = yNode(n3)
                 y4 = yNode(n4)
                 factor = SIN(0.5D0 * (x1 + x2))
                 r2flt = (x2 - x1)**2 + (factor * (y1 - y2))**2
                 r2n1 = (xDatum - x1)**2 + (factor * (yDatum - y1))**2
                 r2n2 = (xDatum - x2)**2 + (factor * (yDatum - y2))**2
                 IF ((r2n1 <= r2flt).OR.(r2n2 <= r2flt)) THEN
!                     Determine internal coordinates for this fault:
                      IF ((r2n1 <= r2flt).AND.(r2n2 <= r2flt)) THEN
                           st = r2n1 / (r2n1 + r2n2)
                           smallc = SQRT(r2flt)
                           smallb = SQRT(r2n1)
                           smalla = SQRT(r2n2)
!                          a = angle at the node-1 end:
                           cosa = (r2n1 + r2flt - r2n2) / (2.0 * smallb * smallc)
                           a = ACOS(cosa)
!                          b = angle at the node-2 end:
                           cosb = (r2n2 + r2flt - r2n1) / (2.0 * smalla * smallc)
                           b = ACOS(cosb)
!                          r = offset from the fault line:
                           r = smallc * (SIN(a) * SIN(b) / SIN(a + b))
                           r2 = r * r
                      ELSE IF (r2n1 <= r2flt) THEN
                           st = 0.0
                           r2 = r2n1
                      ELSE IF (r2n2 <= r2flt) THEN
                           st = 1.0
                           r2 = r2n2
                      END IF
!                     Replace current best solution?
                      IF (r2 <= r2min) THEN
                           i = k
                           s = st
                           r2min = r2
                      END IF
                 END IF
 1340       CONTINUE
            IF (i == 0) THEN
    	      write(ErrorMsg,'(A,I4,A,A,A,F10.4,A,F10.4,A)') "ERR0R: NO FAULT ELEMENT NEAR TO DATUM ",j,' ',TRIM(fName(j))," AT (",fltLon(j),",",fltLat(j),")"
	    	  call FatalError(ErrorMsg,ThID)
            END IF

!         Assuming that i = fault-element was found, continue:

            n1 = nodeF(1, i)
            n2 = nodeF(2, i)
            n3 = nodeF(3, i)
            n4 = nodeF(4, i)
            x1 = xNode(n1)
            x2 = xNode(n2)
            x3 = xNode(n3)
            x4 = xNode(n4)
            y1 = yNode(n1)
            y2 = yNode(n2)
            y3 = yNode(n3)
            y4 = yNode(n4)
            vx1 = v(1, n1)
            vx2 = v(1, n2)
            vx3 = v(1, n3)
            vx4 = v(1, n4)
            vy1 = v(2, n1)
            vy2 = v(2, n2)
            vy3 = v(2, n3)
            vy4 = v(2, n4)
            delVX = fhival(s, vx1, vx2) - fhival(s, vx4, vx3)
            delVY = fhival(s, vy1, vy2) - fhival(s, vy4, vy3)
            s_dp = s
            angle = chord(fArg(1, i), s_dp, fArg(2, i))
            unitX = COS(angle)
            unitY = SIN(angle)
            unitBX = SIN(angle)
            unitBY = -COS(angle)
            spread = delVX * unitBX + delVY * unitBY
            sinist = delVX * unitX + delVY * unitY
            dExtra = -sinist

!   Note conversion of units to mm/year in next line:
            rLt = dExtra * unitL / unitT

            opens = spread
            closes = -spread

            dip = fhival(s, fDip(1, i), fDip(2, i))
            vertic = (ABS(dip - 1.57079632679490D0) < wedge)

            IF (vertic) THEN
                 relVZ = 0.0D0
                 probVZ = 1.0D0
            ELSE

!                Note conversion of units to mm/year in next line:
                 relVZ = closes * ABS(TAN(dip)) * unitL / unitT

                 CALL Likely (delVZ(1, j), delVZ(2, j), relVZ, & ! INTENT(IN)
     &                        probVZ)                            ! INTENT(OUT)
            END IF

            CALL Likely (rLat(1, j), rLat(2, j), rLt, & ! INTENT(IN)
     &                   probRL)                        ! INTENT(OUT)

            shrink = shrink * probRL * probVZ

            rErr = 0.0
            IF ((rLat(1, j) /= 0.0).OR.(rLat(2, j) /= 0.0)) THEN
                 nGSDat = nGSDat + 1
                 IF ((rLt < rLat(1, j)).AND.(rLat(1, j) /= 0.)) THEN
                      nB = nB + 1
                      rErr = rLat(1, j) - rLt
                      sumB = sumB + rErr
                      sumB2 = sumB2 + rErr**2
                      bigB = MAX(bigB, rErr)
                 END IF
                 IF ((rLt > rLat(2, j)).AND.(rLat(2, j) /= 0.)) THEN
                      nB = nB + 1
                      rErr = rLt - rLat(2, j)
                      sumB = sumB + rErr
                      sumB2 = sumB2 + rErr**2
                      bigB = MAX(bigB, rErr)
                 END IF
            END IF

            zErr = 0.0D0
            IF ((delVZ(1, j) /= 0.0D0).OR.(delVZ(2, j) /= 0.0D0)) THEN
                  nGSDat = nGSDat + 1
                  IF ((relVZ < delVZ(1, j)).AND.(delVZ(1, j) /= 0.0D0)) THEN
                       nB = nB + 1
                       zErr = delVZ(1, j) - relVZ
                       sumB = sumB + zErr
                       sumB2 = sumB2 + zErr**2
                       bigB = MAX(bigB, zErr)
                  END IF
                  IF ((relVZ > delVZ(2, j)).AND.(delVZ(2, j) /= 0.0D0)) THEN
                       nB = nB + 1
                       zErr = relVZ - delVZ(2, j)
                       sumB = sumB + zErr
                       sumB2 = sumB2 + zErr**2
                       bigB = MAX(bigB, zErr)
                  END IF
            END IF

!      Convert some numbers to '?' or "---" or "agrees" :

            IF (rLat(1, j) == 0.0D0) THEN
                 c8r1 = '       ?'
            ELSE
                 WRITE(c8r1, "(F8.3)") rLat(1, j)
            END IF

            IF (rLat(2, j) == 0.0D0) THEN
                 c8r2 = '       ?'
            ELSE
                 WRITE(c8r2, "(F8.3)") rLat(2, j)
            END IF

            IF (delVZ(1, j) == 0.0D0) THEN
                 c8v1 = '       ?'
            ELSE
                 WRITE(c8v1, "(F8.3)") delVZ(1, j)
            END IF

            IF (delVZ(2, j) == 0.0D0) THEN
                 c8v2 = '       ?'
            ELSE
                 WRITE(c8v2, "(F8.3)") delVZ(2, j)
            END IF

            IF ((rLat(1, j) == 0.0D0).AND.(rLat(2, j) == 0.0D0)) THEN
                 c8rerr = '       ?'
            ELSE IF (rErr == 0.0D0) THEN
                 c8rerr = '  agrees'
            ELSE
                 WRITE (c8rerr, "(F8.3)") rErr
            END IF

            IF ((delVZ(1, j) == 0.0D0).AND.(delVZ(2, j) == 0.0D0)) THEN
                 c8verr = '       ?'
            ELSE IF (zErr == 0.0D0) THEN
                 c8verr = '  agrees'
            ELSE
                 WRITE (c8verr, "(F8.3)") zErr
            END IF

            IF ((rLat(1, j) == 0.0D0).AND.(rLat(2, j) == 0.0D0)) THEN
                 c7rcha = ' ------'
            ELSE
                 WRITE (c7rcha, "(F7.4)") probRL
            END IF

            IF ((delVZ(1, j) == 0.0D0).AND.(delVZ(2, j) == 0.0D0)) THEN
                 c7vcha = ' ------'
            ELSE
                 WRITE (c7vcha, "(F7.4)") probVZ
            END IF


            IF (vertic) THEN
                 IF(Verbose) WRITE(iUnitVerb, 1370) fName(j)(1:30), &
     &                               c8r1, c8r2, &
     &                               rLt, c8rerr, c7rcha
 1370            FORMAT(' ', A, 2A8, F8.3, A8, A7, &
     &                  ' Throw not used with strike-slip fault.')
            ELSE
                 IF(Verbose) WRITE(iUnitVerb, 1380) fName(j)(1:30), &
     &                               c8r1, c8r2, &
     &                               rLt, c8rerr, c7rcha, &
     &                               c8v1, c8v2, &
     &                               relVZ, c8verr, c7vcha
 1380            FORMAT(' ', A, 2A8, F8.3, A8, A7, 2A8, F8.3, A8, A7)
            END IF
 1390  CONTINUE

!   Summary geologic slip-rate statistics for this model:

       IF(Verbose) WRITE(iUnitVerb, 1392) TRIM(title1), TRIM(title2), TRIM(title3)
 1392  FORMAT (/' FOR MODEL WITH TITLES:'/ &
     &         ' ', A / ' ', A / ' ', A /)
       enB = FLOAT(nB) / FLOAT(nGSDat)
       sumB = sumB / FLOAT(nGSDat)
       sumB2 = SQRT(sumB2 / FLOAT(nGSDat))
       IF(Verbose) WRITE(iUnitVerb, 1394) enB, sumB, sumB2, bigB, shrink
 1394  FORMAT ( &
     & /' FRACTION of predictions in error =', F6.3,'.' &
     & /' AVERAGE error of all predictions =         ', F10.3, ' mm/year.' &
     & /' Root-Mean-Square error of all predictions =', F10.3, ' mm/year.' &
     & /' LARGEST error of any prediction =          ', F10.3, ' mm/year.' &
     & /' Probability that model input parameters and physical equations are consistent' &
     & /' with the true geology represented by this dataset is', D10.3, '.')

       chance = shrink**(1.0D0 / nGSDat)
       IF(Verbose) WRITE(iUnitVerb, 1398) chance
 1398  FORMAT( &
     &        ' ',11('====')/ &
     &  ' #  Estimate, based on this dataset, of     #'/ &
     &  ' #  probability that, in a given case,      #'/ &
     &  ' #  these physical equations and parameter  #'/ &
     &  ' #  values lead to a predicted slip rate    #'/ &
     &  ' #  consistent with geologic data (within   #'/ &
     &  ' #  its uncertainties) is ',F8.5,'.         #'/ &
     &  ' ', 11('====')/)

!      Use root-mean-square error (in mm/year) as summary statistic:
       slpErr = sumB2

!           ****** END GEOLOGIC SLIP-RATE SCORING ********

!-------------------------------------------------------------------

!           ****** SPREADING-RATE SCORING SECTION ********

!   Read seafloor-spreading rates on mid-ocean ridges for scoring:
  501 CALL GetMOR (iUnitM, iUnitVerb, mxMOR, &   ! INTENT(IN)
     &             MORTag, MORPhi, MORThe, &  ! INTENT(OUT)
     &             MORVel, MORSig, numMOR)
       IF (numMOR <= 0) THEN
            spread = 0.0D0
            GO TO 600
       END IF


       toler = 200.0D0 / 6371.0D0
       IF(Verbose) WRITE(iUnitVerb, 530)
  530  FORMAT (// &
     &' Seafloor-spreading velocities versus model predictions:'// &
     &' Ridge E.Lon. N.Lat. Velocity (mm/yr)   Model.V (mm/yr)', &
     &'    Error (mm/yr) (Sigmas)'/ &
     &' ----- ------ ------ --------  -----    -------  -----', &
     &'    -----  -----   ------')
  531  FORMAT (' ', A5, F7.2, F7.2, ES9.2, ' (', F5.1, ')', ES10.2, &
              &' (', F5.1, ')', ES9.2, ' (', F5.1,')  (', F4.1, ')')

       sum0s = 0.0D0
       sum1 = 0.0D0
       sum1s = 0.0D0
       sum2 = 0.0D0
       sum2s = 0.0D0
       sumN = 0.0D0
       sumNs = 0.0D0
       numDen = numMOR
       DO 590 i = 1, numMOR
!          For each mid-ocean ridge spreading-rate,
            theLon = oezOPi * MORPhi(i)
            theLat = 90.0 - oezOPi * MORThe(i)
!          Find closest normal fault:
            rS(1) = SIN(MORThe(i)) * COS(MORPhi(i))
            rS(2) = SIN(MORThe(i)) * SIN(MORPhi(i))
            rS(3) = COS(MORThe(i))
            rMin = 9.99E29
            DO 550 j = 1, nFl
                 dip = 0.5D0 * (fDip(1, j) + fDip(2, j))
!                Dip must be greater than 37.5 degrees:
                 IF (ABS(dip - 1.57079632679490D0) < 0.916294D0) THEN
!                    Dip must be less than 64.0 degrees:
                     IF (ABS(dip - 1.57079632679490D0) > 0.448585D0) THEN
!                         Positions of end-nodes of fault:
                          rA(1) = SIN(xNode(nodeF(1, j))) * &
         &                        COS(yNode(nodeF(1, j)))
                          rA(2) = SIN(xNode(nodeF(1, j))) * &
         &                        SIN(yNode(nodeF(1, j)))
                          rA(3) = COS(xNode(nodeF(1, j)))
                          rB(1) = SIN(xNode(nodeF(2, j))) * &
         &                        COS(yNode(nodeF(2, j)))
                          rB(2) = SIN(xNode(nodeF(2, j))) * &
         &                        SIN(yNode(nodeF(2, j)))
                          rB(3) = COS(xNode(nodeF(2, j)))
!                         Vectpr along fault from node 1 to 2:
                          vf(1) = rB(1) - rA(1)
                          vf(2) = rB(2) - rA(2)
                          vf(3) = rB(3) - rA(3)
                          fLong = SQRT(vf(1)**2 + vf(2)**2 + vf(3)**2)
!                         Vector from node 1 to datum:
                          voff(1) = rS(1) - rA(1)
                          voff(2) = rS(2) - rA(2)
                          voff(3) = rS(3) - rA(3)
!                         Determine internal coordinate (0 ... 1) of datum:
                          frac = (voff(1) * vf(1) + voff(2) * vf(2) + voff(3) * vf(3)) / (fLong**2)
                          IF (frac < 0.0D0) THEN
                               frac = 0.0D0
                               r = SQRT(voff(1)**2 + voff(2)**2 + voff(3)**2)
                          ELSE IF (frac > 1.0D0) THEN
                               frac = 1.0D0
                               r = SQRT((rS(1) - rB(1))**2 + &
         &                              (rS(2) - rB(2))**2 + &
         &                              (rS(3) - rB(3))**2)
                          ELSE
                               CALL DCross (vf, voff, vt)
                               vtl = SQRT(vt(1)**2 + vt(2)**2 + vt(3)**2)
                               r = vtl / fLong
                          END IF
                          IF (r < rMin) THEN
                               rMin = r
                               iEle = j
                               fSave = frac
                          END IF
                     END IF
                 END IF
  550       CONTINUE
            IF (rMin > toler) THEN
                 IF(Verbose) WRITE(iUnitVerb, 551) i, MORTag(i), theLon, theLat, toler
  551            FORMAT (' #', I5, ' ON RISE ', A5,' AT ',F7.2, &
     &                   'E, ',F6.2,'N IS', &
     &                   ' MORE THAN ',F5.3,' RADIANS', &
     &                   ' FROM NORMAL FAULT: IGNORED.')
                 numDen = numDen - 1
            ELSE
!               Compute spreading-rate (iEle, fSave)
                dvat = v(1, nodeF(4, iEle)) - v(1, nodeF(1, iEle))
                dvap = v(2, nodeF(4, iEle)) - v(2, nodeF(1, iEle))
                dvbt = v(1, nodeF(3, iEle)) - v(1, nodeF(2, iEle))
                dvbp = v(2, nodeF(3, iEle)) - v(2, nodeF(2, iEle))
                dvt = dvat + fSave * (dvbt - dvat)
                dvp = dvap + fSave * (dvbp - dvap)
!               Find theta and phi unit-vectors at this point:
                uP(1) = -rS(2)
                uP(2) = rS(1)
                uP(3) = 0.0D0
                CALL Unit (uP) ! INTENT(INOUT)
                CALL DCross (uP, rS, uT)
                voff(1) = dvt * uT(1) + dvp * uP(1)
                voff(2) = dvt * uT(2) + dvp * uP(2)
                voff(3) = dvt * uT(3) + dvp * uP(3)
                rA(1) = SIN(xNode(nodeF(1, iEle))) * &
     &                  COS(yNode(nodeF(1, iEle)))
                rA(2) = SIN(xNode(nodeF(1, iEle))) * &
     &                  SIN(yNode(nodeF(1, iEle)))
                rA(3) = COS(xNode(nodeF(1, iEle)))
                rB(1) = SIN(xNode(nodeF(2, iEle))) * &
     &                  COS(yNode(nodeF(2, iEle)))
                rB(2) = SIN(xNode(nodeF(2, iEle))) * &
     &                  SIN(yNode(nodeF(2, iEle)))
                rB(3) = COS(xNode(nodeF(2, iEle)))
!               Vector along fault from node 1 to 2:
                vf(1) = rB(1) - rA(1)
                vf(2) = rB(2) - rA(2)
                vf(3) = rB(3) - rA(3)
                fLong = SQRT(vf(1)**2 + vf(2)**2 + vf(3)**2)
                CALL DCross (vf, voff, vt)
                dot = vt(1) * rS(1) + vt(2) * rS(2) + vt(3) * rS(3)
                spread = dot / fLong
                bad = ABS(spread - MORVel(i))
                badSig = bad / MORSig(i)
                mvmma = MORVel(i) * secPYr * 1000.0D0
                sprmma = spread * secPYr * 1000.0D0
                badMma = bad * secPYr * 1000.0D0
                IF(Verbose) WRITE(iUnitVerb, 531) MORTag(i), theLon, theLat, &
     &                              MORVel(i), mvmma, &
     &                              spread, sprmma, bad, badMma, badSig
                IF (badSig > 2.0D0) sum0s = sum0s + 1
                sum1 = sum1 + bad
                sum1s = sum1s + badSig
                sum2 = sum2 + bad**2
                sum2s = sum2s + badSig**2
                sumN = MAX(sumN, bad)
                sumNs = MAX(sumNs, badSig)
            END IF
  590  CONTINUE
       IF (numDen > 0) THEN
            sum0s = sum0s / numDen
            sum1 = sum1 / numDen
       ELSE
            sum0s = 0.0D0
            sum1 = 0.0D0
       END IF
       sum1mm = sum1 * 1000.0D0 * secPYr
       IF (numDen > 0) THEN
            sum1s = sum1s / numDen
            sum2 = SQRT(sum2 / numDen)
       ELSE
            sum1s = 0.0D0
            sum2 = 0.0D0
       END IF
       sum2mm = sum2 * secPYr * 1000.0D0
       IF (numDen > 0) THEN
            sum2s = SQRT(sum2s / numDen)
       ELSE
            sum2s = 0.0D0
       END IF
       sumNmm = 1000. * secPYr * sumN
       IF(Verbose) WRITE(iUnitVerb, 599) numDen, sum0s, sum1, sum1mm, sum1s, &
     &                     sum2, sum2mm, sum2s, sumN, sumNmm, sumNs
  599  FORMAT(/' Summary of spreading prediction errors:' &
     &        /' Number of spreading-rates used for scoring: ', I6 &
     &        /' Fraction of predictions wrong by over 2 sigma: ', F5.3 &
     &        /' Mean velocity error:  ', ES10.2, ' (', F7.2, ' mm/yr)' &
     &        /' Mean number of sigmas in error: ', F7.2, &
     &        /' RMS velocity error :  ', ES10.2, ' (', F7.2, ' mm/yr)' &
     &        /' RMS number of sigmas in error: ', F5.2, &
     &        /' Worst velocity error :', ES10.2, ' (', F7.2, ' mm/yr)' &
     &        /' Worst error, in sigmas: ', F7.2)

!   Arbitrarily select an error measure (mean error, in mm/year):
       spread = sum1mm
!          *********** END SPREADING-RATE SCORING ************

!----------------------------------------------------------------

!          *********** BEGIN SEISMICITY SCORING ************

  600  numEqs = 0
       time1 = 9999.0D0 * secPYr
       time2 = 0.0D0

!      Attempt to read catalog file (which may not be present)
!      in format produced by my Seismicity.f90 program:

       IF(Verbose) WRITE(iUnitVerb, 601) iUnitE
  601  FORMAT(/' ATTEMPTING TO READ SEISMIC CATALOG FILE ON UNIT', I3/)

  602  READ(iUnitE, 603, IOSTAT = ios) flag10, &
     &      jYear, month, kDay, &
     &      iHour, minute, iSecon, iTenth, &
     &      eqELon(numEqs + 1), eqNLat(numEqs + 1), &
     &      iDepth, eqMag(numEqs + 1)

  603  FORMAT(A10, &
     &        I4, 1X, I2, 1X, I2, 1X, &
     &        I2, 1X, I2, 1X, I2, 1X, I1, 1X, &
     &        F8.3, 1X, F7.3, 1X, &
     &        I3, F6.2)
       IF (ios /= 0) THEN
!           Problem with last READ (e.g., no file, or end-of-file)
            GO TO 610
       ELSE
            numEqs = numEqs + 1
            IF (numEqs > maxEqs) THEN
	          write(ErrorMsg,'(A)') "INCREASE PARAMETER maxEqs TO EQUAL THE NUMBER OF EARTHQUAKES IN INPUT FILE."
		      call FatalError(ErrorMsg,ThID)
            END IF
            tSec = jYear * secPYr + (month - 1) * 2.6298D6 + (kDay - 1) * 86400.0D0 + &
     &            (iHour - 1) * 3600.0D0 + (minute - 1) * 60.0D0 + iSecon + iTenth / 10.0D0
!           Time1 and time2 are earliest and last earthquake times,
!           in seconds.
            time1 = MIN(time1, tSec)
            time2 = MAX(time2, tSec)
            GO TO 602
       END IF
  610  seismi = 0.0D0
!     (Initializing error measure in case of no catalog)
       IF (numEqs <= 0) GO TO 700

       IF(Verbose) WRITE(iUnitVerb, 611) numEqs
  611  FORMAT(//' -----------------------------------------------' &
     &         /' Beginning seismicity scoring section;' &
     &        //' Read catalog file, found ', I8, ' earthquakes.')

!   Smooth delta-function seismicity with a Gaussian spatial filter,
!   and scale to convert to a scalar strain-rate;
!   evaluate this strain-rate at each node (for plotting):
!   =======================
       sigma = 60.0D0 * 1000.0D0
!   =======================
!      Sigma is the smoothing distance, in meters.
!      Note that it should allow for epicenter mislocation and
!      also for finite source dimensions.
!      However, it is not used to achieve a regionally-smooth
!      scalar strain-rate field; this will be done below.

       sigma2 = sigma**2

       vP = 8100.0D0
       rho = 3300.0D0
!     (SI values appropriate for upper mantle/oceanic lithosphere)
       vS = vP / SQRT(3.)
       SIMu = rho * vS**2
!      SIMu is the elastic shear modulus, in Pa:
       thick = 20000.0D0
!      Thick is the thickness of the seismogenic layer, in m.
       deltaT = time2 - time1
!      deltaT is the length of the catalog, in seconds.
       factor = 1.0D0 / (3.14159265358979D0 * SIMu * sigma2 * thick * deltaT)

       DO 615 i = 1, numNod
            eDotNC(i) = 0.0D0
            uvecN(1, i) = SIN(xNode(i)) * COS(yNode(i))
            uvecN(2, i) = SIN(xNode(i)) * SIN(yNode(i))
            uvecN(3, i) = COS(xNode(i))
  615  CONTINUE

       DO 630 j = 1, numEqs
            IF ((flag10(1:3) == 'ISC').OR.(flag10(1:3) == 'IGN')) THEN
!                ISC or IGN catalog, so magnitudes are body-wave
!                   magnitudes or m_b:
                 eqMB = eqMag(j)

                 eqMS = (eqMB - 2.5D0) / 0.63D0
!                Conversion from m_b to surface-wave magnitude (eqMS)
!                from Richter's Elementary Seismology textbook.

                 eqM0 = 10.0D0**(1.5D0 * eqMS + 9.14D0)
!                Conversion from M_s to moment (eqM0) from
!                Ekstrom & Dziewonski (1998, Nature, vol. 332, pp. 319-323),
!                which updates formula of
!                Kanamori and Anderson (BSSA, vol. 65, p. 1073),
!                in which the constant was +8.5.
!                The formala gives eqM0 in Newton-meters (SI).

            ELSE IF ((flag10(1:10) == 'HarvardCMT').OR.(flag10(1:9) == 'GlobalCMT')) THEN
!                Magnitudes are moment-magnitudes already.
                 eqMW = eqMag(j)
                 eqM0 = 10.0D0**(1.5D0 * eqMW + 9.1D0)
!                per definition of moment-magnitude in
!                Hanks and Kanamori (1979, J. Geophys. Res., pp. 2348-2350).

            ELSE IF (flag10(1:2) == 'NZ') THEN

                ! CHANGED BY Z.L.: NEW ZEALAND LOCAL MAGNITUDE ML IS CALCULATED BASED ON
                ! FOMULA GIVEN BY HAINES IN BSSA,71,275-294,1981. SYSTEMATIC DIFFERENCE
                ! BETWEEN ML AND BODY MAGNITUDE mb EXISTS AND IS DEPTH-DEPENDENT. SEE AL
                ! DAVID HARTE AND DAVID VERE-JONES, N.Z.J.GEOL.GEOPHYS.,42,237-253,1999
                ! VERY ROUGHLY, mb-ML ~0.2 FOR EVENTS SHALLOWER THAN 40KM BASED ON HISTO
                ! GIVEN IN HARTE AND VERE-JONES' PAPER. SO WE SIMPLY ADD 0.2 FROM ML TO
                ! IT INTO mb. THEN CONVERT mb INTO MS AND MW USING ABOVE FORMULA

                 eqMB = eqMag(j) + 0.2D0
                 eqMS = (eqMB - 2.5D0) / 0.63D0
                 eqM0 = 10.0D0**(1.5D0 * eqMS + 9.14D0)

            ELSE
			  write(ErrorMsg,'(A)') "ILLEGAL CATALOG FLAG: ",TRIM(flag10)
			  call FatalError(ErrorMsg,ThID)
            END IF

            equVec(1) = COS(eqNLat(j) * piO180) * &
     &                  COS(eqELon(j) * piO180)
            equVec(2) = COS(eqNLat(j) * piO180) * &
     &                  SIN(eqELon(j) * piO180)
            equVec(3) = SIN(eqNLat(j) * piO180)

            DO 620 i = 1, numNod
                 dot = equVec(1) * uvecN(1, i) + &
     &                 equVec(2) * uvecN(2, i) + &
     &                 equVec(3) * uvecN(3, i)
                 CALL DCross (equVec, uvecN(1, i), cUvec)
                 croSiz = SQRT(cUvec(1)**2 + cUvec(2)**2 + cUvec(3)**2)
                 arc = ATAN2(croSiz, dot)
                 r2 = (arc * radius)**2
                 arg = -r2 / sigma2
                 IF (arg > -80.0D0) THEN
!                     (Avoid underflows for distant nodes)
                      eDotNC(i) = eDotNC(i) + factor * eqM0 * EXP(arg)
                 END IF
  620       CONTINUE
  630  CONTINUE

!      Average catalog strain-rate to element centers
!      and impose reasonable minimum of 0.1% / 5 Ga = 6.0D-21 /s
!     (to prevent correlation test of logarithms from
!      being dominated by artificial, uncertain low-end values!)

       DO 640 i = 1, numEl
            eDotEC(i) = (eDotNC(nodes(1, i)) + &
     &                   eDotNC(nodes(2, i)) + &
     &                   eDotNC(nodes(3, i))) / 3.0D0
            eDotEC(i) = MAX(eDotEC(i), 6.0D-21)
  640  CONTINUE

!   END SEISMIC-CATALOG SCALAR STRAIN-RATE FIELD (ROUGH!)
!   ------------------------------------------------------------
!   BEGIN FINITE-ELEMENT-MODEL SCALAR STRAIN-RATE FIELD

!   Modify velocity solution by averaging together the velocities
!   of all nodes that share a common location, thus forcing the
!   velocity field to be continuous (and eliminating infinite
!   strain-rates in fault elements):

       DO 642 i = 1, numNod
            checkN(i) = .FALSE.
!          (which means "not yet involved in averaging')
  642  CONTINUE
       DO 670 i = 1, nFl
            DO 660 j1 = 1, 2
                 nj1 = nodeF(j1, i)
!               (Fault ends are the only places that can have delta.V's)
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
                      IF (dy >  3.14159265358979D0) dy = dy - 6.28318530717959D0
                      IF (dy < -3.14159265358979D0) dy = dy + 6.28318530717959D0
                      dy = dy * SIN(xNode(nj1))
                      short = SQRT(dx**2 + dy**2)
                      DO 646 k = 1, nFl
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
                                IF (dy >  3.14159265358979D0) dy = dy - 6.28318530717959D0
                                IF (dy < -3.14159265358979D0) dy = dy + 6.28318530717959D0
                                dy = dy * SIN(xNode(nl1))
                                test = SQRT(dx**2 + dy**2)
                                short = MIN(short, test)
                           END IF
  646                 CONTINUE
!                    Collect all corner nodes within 10% of this distance:
                      toler = short / 10.0D0
                      t2 = toler**2
                      DO 650 k = 1, numNod
                           IF (.NOT.checkN(k)) THEN
                                dx = xNode(nj1) - xNode(k)
                                dy = yNode(nj1) - yNode(k)
                                IF (dy >  3.14159265358979D0) dy = dy - 6.28318530717959D0
                                IF (dy < -3.14159265358979D0) dy = dy + 6.28318530717959D0
                                dy = dy * SIN(xNode(nj1))
                                r2 = dx**2 + dy**2
                                IF (r2 < t2) THEN
                                     nInSum = nInSum + 1
                                          IF (nInSum > mxStar) THEN
											write(ErrorMsg,'(A)') "INCREASE VALUE OF PARAMETER maxAtP."
											call FatalError(ErrorMsg,ThID)
                                          END IF
                                     list(nInSum) = k
                                     checkN(k) = .TRUE.
                                END IF
                           END IF
  650                 CONTINUE
!                    Averaging actually begins here, based on list:
                      xSum = 0.0D0
                      ySum = 0.0D0
                      DO 652 k = 1, nInSum
                           xSum = xSum + v(1, list(k))
                           ySum = ySum + v(2, list(k))
  652                 CONTINUE
                      xMean = xSum / nInSum
                      yMean = ySum / nInSum
                      DO 656 k = 1, nInSum
                           v(1, list(k)) = xMean
                           v(2, list(k)) = yMean
  656                 CONTINUE
                 END IF
  660       CONTINUE
  670  CONTINUE

!   Compute model strain-rates after removing delta.V's:

       CALL EDot (dxs, dys, &
     &            fpsfer, mxEl, &
     &            mxNode, nodes, numEl, radius, sita, v, & ! INTENT(IN)
     &            eRate)                                  ! INTENT(OUT)

!   Convert element-center values to scalars:
!   ((e3*-e1*)/2) at integration point 1
!  (Note: * means strain-rate partitioned(?) if E2 /= 0)

       DO 672 i = 1, numEl
            exx = eRate(1, 1, i)
            eyy = eRate(2, 1, i)
            exy = eRate(3, 1, i)
            CALL Prince (exx, eyy, exy, &            ! INTENT(IN)
     &                   e1, e2, u1x, u1y, u2x, u2y) ! INTENT(OUT)
            ez = -(exx + eyy)
            IF ((e2 == 0.0D0).AND.(e1 == 0.0D0)) THEN
                 eDotEM(i) = 1.0D-20
            ELSE
                 IF ((e2 * ez) > 0.0D0) THEN
!                e1 has the unique sign and is partitioned:
                      e1Part = .TRUE.
                      e2Part = .FALSE.
                      eZPart = .FALSE.
                 ELSE IF ((e1 * ez) > 0.0D0) THEN
!                e2 has the unique sign and is partitioned:
                      e1Part = .FALSE.
                      e2Part = .TRUE.
                      eZPart = .FALSE.
                 ELSE
!                ezz has the unique sign and is partitioned:
                      e1Part = .FALSE.
                      e2Part = .FALSE.
                      eZPart = .TRUE.
                 END IF
!                strike-slip rate (e2me1)
                 IF (e1Part) THEN
                      e2me1 = 2.0D0 * ABS(e2)
                 ELSE IF (e2Part) THEN
                      e2me1 = 2.0D0 * ABS(e1)
                 ELSE
                      e2me1 = ABS(e2 - e1)
                 END IF
!                thrust-faulting rate (ezme1)
                 IF (e1Part) THEN
                      ezme1 = 2.0D0 * ABS(ez)
                 ELSE IF (eZPart) THEN
                      ezme1 = 2.0D0 * ABS(e1)
                 ELSE
                      ezme1 = ABS(ez - e1)
                 END IF
!                normal-faulting rate (e2mez)
                 IF (e2Part) THEN
                      e2mez = 2.0D0 * ABS(ez)
                 ELSE IF (eZPart) THEN
                      e2mez = 2.0D0 * ABS(e2)
                 ELSE
                      e2mez = ABS(e2 - ez)
                 END IF
                 eLarge = MAX(e2me1, ezme1, e2mez)
                 eDotEM(i) = MAX(eLarge / 2.0D0, 1.0D-20)
            END IF
  672  CONTINUE

!       Repeatedly smooth scalar strain-rate fields
!       by averaging element values onto nodes
!      (treating faults as healed),
!       and then interpolating from nodes to element centers.

!       The same smoothing operator is applied to both
!       catalog and finite-element-model scalar strain-rate fields.

!   ====================================
       nBland = 16
!   ====================================

       nBland = MAX(nBland, 1)
       IF(Verbose) WRITE(iUnitVerb, 673) nBland
  673  FORMAT(/' Using ', I2, ' smoothing cycles on both strain-rates.')
       DO 684 nB = 1, nBland

!        (A) Extrapolate from elements to nodes (across faults):

            DO 674 i = 1, numNod
                 eDotNC(i) = 0.0D0
                 eDotNM(i) = 0.0D0
                 atnode(i) = 0.0D0
  674       CONTINUE
            DO 676 i = 1, numEl
                 DO 675 j = 1, 3
                      node = nodes(j, i)
                      eDotNC(node) = eDotNC(node) + eDotEC(i)
                      eDotNM(node) = eDotNM(node) + eDotEM(i)
                      atnode(node) = atnode(node) + 1.0D0
  675            CONTINUE
  676       CONTINUE
            DO 678 i = 1, nFl
                 DO 677 j = 1, 2
                      node1 = nodeF(j, i)
                      node2 = nodeF(5 - j, i)
                      t1 = eDotNC(node1)
                      t2 = eDotNC(node2)
                      eDotNC(node1) = eDotNC(node1) + t2
                      eDotNC(node2) = eDotNC(node2) + t1
                      t1 = eDotNM(node1)
                      t2 = eDotNM(node2)
                      eDotNM(node1) = eDotNM(node1) + t2
                      eDotNM(node2) = eDotNM(node2) + t1
                      t1 = atnode(node1)
                      t2 = atnode(node2)
                      atnode(node1) = atnode(node1) + t2
                      atnode(node2) = atnode(node2) + t1
  677            CONTINUE
  678       CONTINUE
            DO 679 i = 1, numNod
                 eDotNC(i) = eDotNC(i) / MAX(atnode(i), 1.0D0)
                 eDotNM(i) = eDotNM(i) / MAX(atnode(i), 1.0D0)
  679       CONTINUE

!        (B) Interpolate from nodes to element centers:

            DO 680 i = 1, numEl
                 eDotEC(i) = (eDotNC(nodes(1, i)) + &
     &                        eDotNC(nodes(2, i)) + &
     &                        eDotNC(nodes(3, i))) / 3.0D0
                 eDotEM(i) = (eDotNM(nodes(1, i)) + &
     &                        eDotNM(nodes(2, i)) + &
     &                        eDotNM(nodes(3, i))) / 3.0D0
  680       CONTINUE
  684  CONTINUE

!   ----- Output smoothed scalar strain-rate fields
!                     to two .FEG files: ----------------------

       WRITE(c3, "(I3)") nBland
       IF(Verbose) WRITE(iUnitVerb, 686) sigma, iUnitO
  686  FORMAT(/' EQ magnitudes (m_b) have been converted to moments,' &
     &        /' the moments distributed with a Gaussian filter of' &
     &        /' characteristic distance ', ES10.3,' m;' &
     &        /' then the scalar strain-rate was diffusively smoothed;' &
     &        /' Use the .feg file written to unit ', I3 &
     &        ,' to plot the results.'/)
       enctit = 'STRAIN-RATES FROM SEISMIC CATALOG, ' // 'diffusively smoothed (n =' // c3 // ')'
       CALL PutNet (iUnitO, &
     &              .TRUE., eDotNC, eDotNC, fDip, &
     &              mxEl, mxFEl, mxNode, n1000, &
     &              nFakeN, nFl, nodeF, nodes, &
     &              nRealN, numEl, numNod, offset, &
     &              enctit, eDotNC, xNode, yNode, eDotNC) ! INTENT(IN)

       IF(Verbose) WRITE(iUnitVerb, 688) iUnitZ
  688  FORMAT(/' Finite-element model velocities were made continuous' &
     &        /' by healing all active faults;' &
     &        /' then the scalar strain-rate was diffusively smoothed;' &
     &        /' Use the .feg file written to unit ', I3 &
     &        ,' to plot the results.'/)
       enctit = 'STRAIN-RATES FROM HEALED-FAULT FE MODEL, ' // 'diffusively smoothed (n =' // c3 // ')'
       CALL PutNet (iUnitZ, &
     &              .TRUE., eDotNM, eDotNM, fDip, &
     &              mxEl, mxFEl, mxNode, n1000, &
     &              nFakeN, nFl, nodeF, nodes, &
     &              nRealN, numEl, numNod, offset, &
     &              enctit, eDotNM, xNode, yNode, eDotNM) ! INTENT(IN)

!   -----------------------------------------------

!   Compute correlation coefficient between:
!    *log10(smoothed scalar strain-rate from seismic catalog)
!    *log10(smoothed scalar strain-rate from fault-healed model)

!   Convert to common logarithms:

       DO 690 i = 1, numEl
            eDotEC(i) = LOG10(MAX(eDotEC(i), 1.0D-30))
            eDotEM(i) = LOG10(MAX(eDotEM(i), 1.0D-30))
  690  CONTINUE

!   Find averages of log10's of scalar strain-rates; catalog & model:

       aveC = 0.0D0
       aveM = 0.0D0
       sum = 0.0D0
       DO 693 i = 1, numEl
            aveC = aveC + area(i) * eDotEC(i)
            aveM = aveM + area(i) * eDotEM(i)
            sum = sum + area(i)
  693  CONTINUE
       aveC = aveC / sum
       aveM = aveM / sum

       varC = 0.0D0
       varM = 0.0D0
       crossP = 0.0D0
       DO 694 i = 1, numEl
            varC = varC + area(i) * (eDotEC(i) - aveC)**2
            varM = varM + area(i) * (eDotEM(i) - aveM)**2
            crossP = crossP + area(i) * (eDotEC(i) - aveC) * (eDotEM(i) - aveM)
  694  CONTINUE
       correl = crossP / SQRT(varC * varM)
       IF(Verbose) WRITE(iUnitVerb, 695) correl
  695  FORMAT(/' Correlation coefficient between:' &
     &        /'   * log10(smoothed scalar strain-rate from seismic catalog)' &
     &        /'   * log10(smoothed scalar strain-rate from fault-headed model)' &
     &        /' is: ', F7.4)

!   Arbitrarily select an error measure (correlation coefficient)
       seismi = correl

!          *********** END SEISMICITY SCORING ************

!----------------------------------------------------------------

!  ****** UPPER-MANTLE SEISMIC ANISOTROPY SCORING SECTION ********

! N.B. The basic assumption is that upper-mantle seismic anisotropy
!     (expressed as SKS fast-polarization {"phi"} azimuths and
!      splitting times in s} is due to simple shear in the
!      asthenosphere, and thus is comparable to directions of
!      basal shear tractions inferred from a torque-report file
!     (but only for slab-less plates, with slab_q(iPlate) = .FALSE.).

!      Bear in mind that anisotropy could also be due to several other causes:
!  (1) Most-extensional direction of Neogene finite strain in lithosphere.
!  (2) Most-extensional direction of ancient finite strain in lithosphere.
!  (3) Azimuth of aligned magma-filled dikes (or water-filled cracks).
!      and in these cases the data would NOT be comparable as assumed!

!   READ fast-polarization azimuths ("phi") and splitting times
!      of SKS waves (e.g., compilation of Matt Fouch).

  700  CALL GetSKS (iUnitK, iUnitVerb, mxSKS,  &
     &              numSKS, SKS_tag, SKS_theta, SKS_phi, & ! INTENT(IN)
     &              SKS_argument, SKS_delay)               ! INTENT(OUT)
       IF (numSKS <= 0) THEN
            anisotropy = 0.0D0 ! summary error measure
            GO TO 900
       END IF
       IF(Verbose) WRITE(iUnitVerb, 701) numSKS, iUnitK
  701  FORMAT (/' ',I6,' Fast SKS azimuths and delay times were read from unit ', I2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (DIs_Initialization_Needed()) THEN
          subdivision = 4
          IF(Verbose) WRITE(iUnitVerb, "(/' Initializing uniform global grid at subdivision ', &
     &                I1, ' for area-weighting of data...')") subdivision
          CALL DInitialize_Weighting (subdivision)
       END IF

       IF(Verbose) WRITE(iUnitVerb, "(' Computing area-weights for SKS azimuth/delay data...')")
       IF (ALLOCATED(weights)) DEALLOCATE ( weights )
       ALLOCATE ( weights(numSKS) )
       CALL DPerform_Weighting (number_of_data = numSKS, &
                              & theta_radians  = SKS_theta, &
                              & phi_radians    = SKS_phi, & ! INTENT(IN)
                              & weights = weights)          ! INTENT(OUT)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!   Read plate outlines, for assigning sites to the correct plate:
       IF(Verbose) WRITE(iUnitVerb, 2) iUnitF
    2  FORMAT (/' ATTEMPTING TO READ OUTLINES OF PLATES FROM UNIT', I3/)
       CALL GetPBx (iUnitF, iUnitVerb, names, nPBnd, nPlate, & ! INTENT(IN)
     &              nBoundaryPoints, pLat, pLon)            ! INTENT(OUT)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       IF(Verbose) WRITE(iUnitVerb, "(' Assigning SKS azimuth/delay data to plates...')")
       DO 720 i = 1, numSKS
            CALL ThetaPhi_2_pLate (iUnitVerb, &
     &                             nPBnd, nBoundaryPoints, &
     &                             nPlate, pLat, pLon, &
     &                             SKS_theta(i), SKS_phi(i),  & ! INTENT(IN)
     &                             iPlate)                      ! INTENT(OUT)
            SKS_iPlate(i) = iPlate
  720  CONTINUE

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       IF(Verbose) WRITE(iUnitVerb, "(' Finding basal shear azimuths to compare with SKS azimuths...')")
       CALL Tractor(iunitQ, iUnitVerb, nPlate, numSKS, &
     &              slab_q, SKS_phi, SKS_theta, SKS_iPlate, &  ! INTENT(IN)
     &              basal_shear_tractions)                     ! INTENT(OUT)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    Prepare to compute azimuths of fastest-horizontal-extension principal strain-rate
!    axes, as alternate predictor of fast-polarization azimuths (phi) of SKS data:

!   Compute the (theta, phi) and then (x, y, z) of integration points:

       DO 750 i = 1, numEl
            DO 725 k = 1, 3
                 CartVs(1, k) = COS(yNode(nodes(k, i))) * &
     &                          SIN(xNode(nodes(k, i)))
                 CartVs(2, k) = SIN(yNode(nodes(k, i))) * &
     &                          SIN(xNode(nodes(k, i)))
                 CartVs(3, k) = COS(xNode(nodes(k, i)))
  725       CONTINUE
            DO 740 m = 1, 7
                 DO 730 j = 1, 3
                      tempV(j) = 0.0
                      DO 728 k = 1, 3
                           tempV(j) = tempV(j) + CartVs(j, k) * points(k, m)
  728                 CONTINUE
  730            CONTINUE
                 CALL Unit (tempV) ! INTENT(INOUT)
                 IF (ABS(tempV(3)) <= 0.5D0) THEN
                      xIP(m, i) = ACOS(tempV(3))
                 ELSE
                      equPar = SQRT(tempV(1)**2 + tempV(2)**2)
                      xIP(m, i) = ATan2F(equPar, tempV(3))
                 END IF
                 yIP(m, i) = ATan2F(tempV(2), tempV(1))
                 CartR(1, m, i) = SIN(xIP(m, i)) * COS(yIP(m, i))
                 CartR(2, m, i) = SIN(xIP(m, i)) * SIN(yIP(m, i))
                 CartR(3, m, i) = COS(xIP(m, i))
  740       CONTINUE
  750  CONTINUE

!   Strain rates at integration points:

       CALL EDot (dxs, dys, &
     &            fpsfer, mxEl, &
     &            mxNode, nodes, numEl, radius, sita, v, & ! INTENT(IN)
     &            eRate)                                   ! INTENT(OUT)

!   Note: Integration point must be within one-half-element-
!         width of datum, or datum is not used.
!         Here, the typical element side is determined by
!         averaging together all element sides,
!         and expressed in radians:

       sumSid = 0.0D0
       DO 756 i = 1, numEl
            DO 753 k = 1, 3
                 kp1 = 1 + MOD(k, 3)
                 n1 = nodes(k, i)
                 n2 = nodes(kp1, i)
                 x1 = xNode(n1)
                 x2 = xNode(n2)
                 y1 = yNode(n1)
                 y2 = yNode(n2)
                 rA(1) = SIN(x1) * COS(y1)
                 rA(2) = SIN(x1) * SIN(y1)
                 rA(3) = COS(x1)
                 rB(1) = SIN(x2) * COS(y2)
                 rB(2) = SIN(x2) * SIN(y2)
                 rB(3) = COS(x2)
                 dot = rA(1) * rB(1) + rA(2) * rB(2) + rA(3) * rB(3)
                 CALL DCross (rA, rB, cUvec)
                 croSiz = SQRT(cUvec(1)**2 + cUvec(2)**2 + cUvec(3)**2)
                 side = ATan2F(croSiz, dot)
                 sumSid = sumSid + side
  753       CONTINUE
  756  CONTINUE
       toler = sumSid / (6.0D0 * numEl)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       IF(Verbose) WRITE(iUnitVerb, 760)
  760  FORMAT (//' Fast SKS-polarization azimuths versus model predictions (in slab-less plates):'/  &
               &/' Datum E.Long. N.Latt. Plate Delay,s Azimuth Model.Az1,2 Error Area-weight Slabless plate?' &
               &/' ----- ------- ------- ----- ------- ------- ----------- ----- ----------- ---------------')
  761  FORMAT (  ' ', A5, 1X, F7.2, 1X, F7.2, 4X, A2, 1X, F7.2, 1X, F7.0, 1X, F5.0, ',', F5.0, 1X, F5.0, 1X, F11.4, 1X, '    .TRUE.     ')
  762  FORMAT (  ' ', A5, 1X, F7.2, 1X, F7.2, 4X, A2, 1X, F7.2, 1X, F7.0, 1X, '   ??,   ??', 1X, '   --', 9X, '---', 1X, '    .FALSE.    ')

!: Initialize accumulators:
       nSKS_bad = 0  ! SKS data misfit by more than 20 degrees (L0 numerator)
       nSKS_used = 0 ! Count of data in slabless plates (L0 denominator)
       sum1 = 0.0D0  ! Numerator for L1 (mean size of error)
       sum1d = 0.0D0 ! Denominator for L1
       sum2 = 0.0D0  ! Numerator for L2 (RMS size of error)
       sum2d = 0.0D0 ! Denominator for L2

       DO 790 i = 1, numSKS

!           Next SKS datum:
            theLon = oezOPi * SKS_phi(i)
            IF (theLon >  180.0D0) theLon = theLon - 360.0D0
            IF (theLon < -180.0D0) theLon = theLon + 360.0D0
            theLat = 90.0D0 - oezOPi * SKS_theta(i)
           !Cartesian unit vector:
            rS(1) = SIN(SKS_theta(i)) * COS(SKS_phi(i))
            rS(2) = SIN(SKS_theta(i)) * SIN(SKS_phi(i))
            rS(3) = COS(SKS_theta(i))

           !Datum: Azimuth of fast polarization of SKS, in degrees:
            azimuth = 180.0D0 - oezOPi * SKS_argument(i)
            IF (azimuth <    0.0D0) azimuth = azimuth + 180.0D0
            IF (azimuth <    0.0D0) azimuth = azimuth + 180.0D0
            IF (azimuth >= 180.0D0) azimuth = azimuth - 180.0D0
            IF (azimuth >= 180.0D0) azimuth = azimuth - 180.0D0

            IF (slab_q(SKS_iPlate(i))) THEN

                !No scoring; just echo the datum, without providing any model predictions...
                 IF(Verbose) WRITE(iUnitVerb, 762) SKS_tag(i), theLon, theLat, names(SKS_iPlate(i)), SKS_delay(i), azimuth

            ELSE ! slabless-plate; proceed with scoring!

                !Predictor1: Azimuth of simple shear in asthenosphere, in degrees:
                 bAzim = 180.0D0 - oezOPi * ATan2F(basal_shear_tractions(2, i), basal_shear_tractions(1, i)) ! (y = phi = E, x = theta = S)
                 IF (bAzim <    0.0D0) bAzim = bAzim + 180.0D0
                 IF (bAzim <    0.0D0) bAzim = bAzim + 180.0D0
                 IF (bAzim >= 180.0D0) bAzim = bAzim - 180.0D0
                 IF (bAzim >= 180.0D0) bAzim = bAzim - 180.0D0

                !Predictor2: Azimuth of e3 (fastest-extension) horizontal principal strain-rate in lithosphere, in degrees:
                !(a) Find closest integration point in the grid...
                 r2min = 9.99D29
                 DO 770 m = 1, 7
                      DO 765 j = 1, numEl
                           r2 = (rS(1) - CartR(1, m, j))**2 + &
     &                          (rS(2) - CartR(2, m, j))**2 + &
     &                          (rS(3) - CartR(3, m, j))**2
                           IF (r2 < r2min) THEN
                                iEle = j
                                mIP = m
                                r2min = r2
                           END IF
  765                 CONTINUE
  770            CONTINUE
                 rMin = SQRT(r2min)
                 IF (rMin > toler) THEN
                      IF(Verbose) WRITE(iUnitVerb, 780) strtag(i), theLon, theLat, toler
  780                 FORMAT (' ', A5, 1X, F7.2, 1X, F7.2, &
     &                        ' SKS DATUM IS MORE THAN ', &
     &                          F6.4, ' RADIANS FROM NEAREST' &
     &                       ,' INTEGRATION POINT: IGNORED.')
                      e3Azim = bAzim ! so prediction errors will be the same, and bAzim error won't be underbid.
                 ELSE
                    !(b) Get e_1H (fastest-horizontal-shortening) direction:
                      e11h = eRate(1, mIP, iEle)
                      e22h = eRate(2, mIP, iEle)
                      e12h = eRate(3, mIP, iEle)
                      CALL Prince (e11h, e22h, e12h, & ! INTENT(IN)
     &                             e1h, e2h, u1x, u1y, u2x, u2y) ! INTENT(OUT)
                      err = -(e11h + e22h)
                     !(or, = -(e1h + e2h) )
                      e1 = MIN(e1h, e2h, err)
                      e3 = MAX(e1h, e2h, err)
                      IF ((err > e1).AND.(err < e3)) THEN
                           e2 = err
                      ELSE IF ((e1h > e1).AND.(e1h < e3)) THEN
                           e2 = e1h
                      ELSE
                           e2 = e2h
                      END IF
                      pAzim = 180.0D0 - oezOPi * ATan2F(u1y, u1x)
                      IF (pAzim <    0.0D0) pAzim = pAzim + 180.0D0
                      IF (pAzim <    0.0D0) pAzim = pAzim + 180.0D0
                      IF (pAzim >= 180.0D0) pAzim = pAzim - 180.0D0
                      IF (pAzim >= 180.0D0) pAzim = pAzim - 180.0D0
                      e3Azim = pAzim + 90.0D0
                 END IF
                 IF (e3Azim <    0.0D0) e3Azim = e3Azim + 180.0D0
                 IF (e3Azim <    0.0D0) e3Azim = e3Azim + 180.0D0
                 IF (e3Azim >= 180.0D0) e3Azim = e3Azim - 180.0D0
                 IF (e3Azim >= 180.0D0) e3Azim = e3Azim - 180.0D0

                !Prediction error of Predictor1 (simple shear in asthenosphere):
                 simple_bad = ABS(bAzim - azimuth)
                 IF (simple_bad > 90.0D0) simple_bad = ABS(simple_bad - 180.0D0)
                 IF (simple_bad > 90.0D0) simple_bad = ABS(simple_bad - 180.0D0)
                 IF (simple_bad > 90.0D0) simple_bad = ABS(simple_bad - 180.0D0)
                 IF (simple_bad > 90.0D0) simple_bad = ABS(simple_bad - 180.0D0)

                !Prediction error of Predictor2 (pure-shear in lithosphere):
                 pure_bad = ABS(e3Azim - azimuth)
                 IF (pure_bad > 90.0D0) pure_bad = ABS(pure_bad - 180.0D0)
                 IF (pure_bad > 90.0D0) pure_bad = ABS(pure_bad - 180.0D0)
                 IF (pure_bad > 90.0D0) pure_bad = ABS(pure_bad - 180.0D0)
                 IF (pure_bad > 90.0D0) pure_bad = ABS(pure_bad - 180.0D0)

                !Select smaller of two prediction errors:
                 bad = MIN(simple_bad, pure_bad)

                 IF (bad > 20.0D0) nSKS_bad = nSKS_bad + 1        ! SKS data misfit by more than 20 degrees (L0 numerator)
                 nSKS_used = nSKS_used + 1                        ! Count of data in slabless plates (L0 denominator)
                 sum1 = sum1 + bad * SKS_delay(i) * weights(i)    ! Numerator for L1 (mean size of error)
                 sum1d = sum1d + SKS_delay(i) * weights(i)        ! Denominator for L1
                 sum2 = sum2 + bad**2 * SKS_delay(i) * weights(i) ! Numerator for L2 (RMS size of error)
                 sum2d = sum2d + SKS_delay(i) * weights(i)        ! Denominator for L2

                 IF(Verbose) WRITE(iUnitVerb, 761) SKS_tag(i), theLon, theLat, names(SKS_iPlate(i)), SKS_delay(i), azimuth, bAzim, e3Azim, bad, weights(i)

            END IF ! slab_q(SKS_iPlate(i)), or NOT <== proceed with scoring
  790  CONTINUE

       IF (nSKS_used > 0) THEN
            perBad = (100.0D0 * nSKS_bad) / (1.0D0 * nSKS_used)
       ELSE
            perBad = 0.0D0
       END IF
       IF (sum1d > 0.0D0) THEN
            sum1 = sum1 / sum1d
       ELSE
            sum1 = 0.0D0
       END IF
       IF (sum2d > 0.0D0) THEN
            sum2 = SQRT(sum2 / sum2d)
       ELSE
            sum2 = 0.0D0
       END IF
       IF(Verbose) WRITE(iUnitVerb, 796) nSKS_used, perBad, sum1, sum2
  796  FORMAT (/' SUMMARY OF SEISMIC ANISOTROPY ERR0RS:' &
     &         /' Number of SKS data in slabless plates, used for scoring: ', I6, '.' &
     &         /' Percentage of directions mis-predicted by more than 20 degrees: ', F6.2 &
     &         /' Area- & delay-weighted MEAN error: ', F5.2, ' degrees.' &
     &         /' Area- & delay-weighted RMS  error: ', F5.2, ' degrees.')

!   Arbitrarily choose measure of error (mean error, in degrees):
       anisotropy = sum1

       DEALLOCATE (weights)

!           ********* END ANISOTROPY SCORING SECTION **********

!----------------------------------------------------------------

!             ******   SUMMARY OF SCORES ********

  900 IF(Verbose) then
        WRITE(iUnitVerb, 901)
        IF (geodes     /= 0.0D0) WRITE(iUnitVerb, 910) geodes
        IF (spread     /= 0.0D0) WRITE(iUnitVerb, 920) spread
        IF (stress     /= 0.0D0) WRITE(iUnitVerb, 930) stress
	    IF (slpErr     /= 0.0D0) WRITE(iUnitVerb, 940) slpErr
	    IF (seismi     /= 0.0D0) WRITE(iUnitVerb, 950) seismi
	    IF (anisotropy /= 0.0D0) WRITE(iUnitVerb, 970) anisotropy
	  END IF
    
	inquire(unit=iUnitOrbDat, opened=WScores)
	if(WScores) then
      write(iUnitOrbDat, 901)
      if (geodes     /= 0.0D0) write(iUnitOrbDat, 910) geodes
      if (spread     /= 0.0D0) write(iUnitOrbDat, 920) spread
      if (stress     /= 0.0D0) write(iUnitOrbDat, 930) stress
	  if (slpErr     /= 0.0D0) write(iUnitOrbDat, 940) slpErr
	  if (seismi     /= 0.0D0) write(iUnitOrbDat, 950) seismi
	  if (anisotropy /= 0.0D0) write(iUnitOrbDat, 970) anisotropy
	  write(iUnitOrbDat, *) ''
	end if


  901 FORMAT (' ------------------------------------------------' &
     &      //' REVIEW OF SUMMARY MEASURES OF ERROR OR QUALITY:')
  910 FORMAT ('   GEODETIC VELOCITY       (ERROR): ', F7.2)
  920 FORMAT ('   SEAFLOOR SPREADING RATE (ERROR): ', F7.2)
  930 FORMAT ('   STRESS DIRECTION        (ERROR): ', F7.2)
  940 FORMAT ('   FAULT SLIP RATE         (ERROR): ', F7.2)
  950 FORMAT ('   (SMOOTHED) SEISMICITY CORRELATION (QUALITY): ', F7.4)
  970 FORMAT ('   SEISMIC ANISOTROPY      (ERROR): ', F7.2)

      IF(Verbose) WRITE(iUnitVerb, 999) TRIM(title1), TRIM(title2), TRIM(title3)
  999 FORMAT (/' FOR MODEL:' &
     &        /' ', A &
     &        /' ', A &
     &        /' ', A/)

    misfits    = 0.0D0
	GeoMisfit  = geodes
	SprdMisfit = spread
    StrsMisfit = stress
    SlpMisfit  = slpErr
    SeisMisfit = seismi
    AniMisfit  = anisotropy
	misfits    = (/GeoMisfit,SprdMisfit,StrsMisfit,SlpMisfit,SeisMisfit,AniMisfit,0.0D0/)

END SUBROUTINE

end module
