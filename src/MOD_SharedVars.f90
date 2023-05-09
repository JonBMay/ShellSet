!*******************************************************************************
! Contains shared global variables and block data
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

module SharedVars

implicit none


logical,dimension(3),parameter :: plt = (/.True.,.True.,.True./)
character(len=200) :: ErrorMsg
real,dimension(:),allocatable :: ErrorArray
character(len=2),dimension(:),allocatable :: ErrorArrayChar
integer,dimension(:),allocatable :: ErrorArrayInt
integer :: ThID
character(len=100) :: title1, title2, title3
logical :: Verbose, everyP
integer,parameter :: iUnitVerb = 6 ! Unit number for optional verbose.txt file
integer,parameter :: iUnitOrbDat = 66 ! Unit number for OrbScore misfits scores during misfit convergence

integer,parameter :: nPlate = 52 ! Number of pLates in PB2002 model of Bird [2003]
integer :: iConve,iPVRef,maxItr
character(len=2) :: pltRef

real*8 :: alphaT, conduc, constr, &
        & fFric, cFric, Biot, Byerly, aCreep, bCreep, cCreep, dCreep, eCreep, & ! default d_XXXX = LR_set_XXXX(0)
        & dipMax, etaMax, fMuMax, gMean, gradie, &
        & offMax, okDelV, okToQt, omega, oneKm, &
        & radio, radius, refStr, rhoAst, rhoBar, rhoH2O, &
        & slide, subDip, tAdiab, tauMax, temLim, trHMax, tSurf, &
        & vTimes, visMax, wedge, zBAsth

dimension alphaT(2), conduc(2), &
        & aCreep(2), bCreep(2), cCreep(2), dCreep(2), & ! default d_XXXX(1:2) = LR_set_XXXX(1:2, 0)
        & radio(2),  rhoBar(2), tauMax(2), temLim(2), &
		& omega(3, nPlate)

integer :: i,j

character*2  :: names
dimension names(nPlate)


!  Following rotation vectors in Cartesian (x,y,z) components,
!  with units of radians per million years:
!  [Bird, 2003, G**3, Table 1]
       data ((omega(i, j), i = 1, 3), j = 1, nPlate) / &
     &    0.002401, -0.007939,  0.013892, &
     &    0.000949, -0.008643,  0.013725, &
     &    0.000689, -0.006541,  0.013676, &
     &    0.002042, -0.013153,  0.008856, &
     &    0.008570, -0.005607,  0.017497, &
     &    0.000148, -0.003070,  0.010915, &
     &    0.015696,  0.002467,  0.023809, &
     &    0.009349,  0.000284,  0.016252, &
     &    0.000184,  0.005157,  0.001150, &
     &   -0.000871, -0.002268,  0.002507, &
     &   -0.019124,  0.030087,  0.010227, &
     &    0.011506, -0.044526,  0.007197, &
     &    0.001688, -0.009048,  0.012815, &
     &    0.003716, -0.003791,  0.000949, &
     &   -0.008915, -0.026445,  0.020895, &
     &   -0.061175,  0.005216, -0.013755, &
     &    0.070136,  0.160534,  0.094328, &
     &    0.000529, -0.007235,  0.013123, &
     &   -0.083251, -0.002464, -0.014923, &
     &    0.016256,  0.089364,  0.015035, &
     &    0.008180, -0.004800,  0.016760, &
     &    0.006512,  0.003176,  0.005073, &
     &    0.108013,  0.299461,  0.230528, &
     &    0.033318, -0.001813,  0.036441, &
     &   -0.013835,  0.008245,  0.015432, &
     &   -0.777844,  0.440872, -0.047437, &
     &    0.001521,  0.007739,  0.013437, &
     &    0.038223, -0.058291,  0.013679, &
     &    0.001768, -0.008439,  0.009817, &
     &   -0.004336,  0.003769, -0.000402, &
     &    0.000111, -0.006361,  0.010449, &
     &    0.044913, -0.009546,  0.010601, &
     &   -0.055342, -0.010890,  0.006794, &
     &   -0.000022, -0.013417,  0.019579, &
     &    0.001041, -0.008305,  0.012143, &
     &   -0.026223,  0.020184,  0.037208, &
     &    0.000000,  0.000000,  0.000000, &
     &   -0.000040, -0.009291,  0.012815, &
     &    0.012165, -0.012510, -0.000366, &
     &   -0.019183, -0.070604,  0.036798, &
     &    0.000472, -0.006355,  0.009100, &
     &    0.121443, -0.078836,  0.027122, &
     &    0.001117, -0.007434,  0.008534, &
     &   -0.000833, -0.006701,  0.013323, &
     &    0.001287, -0.008754,  0.014603, &
     &   -0.017196,  0.017186,  0.008623, &
     &    0.003201, -0.010440,  0.015854, &
     &    0.023380, -0.019369, -0.010465, &
     &   -0.009400,  0.023063,  0.008831, &
     &    0.142118,  0.005616,  0.078214, &
     &   -0.016831,  0.018478,  0.010166, &
     &   -0.000836, -0.006169,  0.016274 /

!--------------------------------------------------------------------
!                           DATA statements

!  "subDip" is the maximum dip (from horizontal, in degrees) for a
!     fault in a whole-Earth model (sphere = .TRUE.) to be treated as
!     a subduction zone (in which case, the footwall nodes require
!     boundary conditions).
!   In all models, faults with less than this dip have the down-dip
!     integral of traction limited to tauMax.  "tauMax" is an array
!     of two values, for oceanic and continental subduction zones,
!     respectively.  If such limits are not wanted, then
!     the tauMax vaules can be set to very large numbers (e.g.,9.99E29).
       data subDip / 19.D0 /

!  "dipMax" is the maximum dip (from horizontal, in degrees) for a
!   fault element to be treated as a dip-slip fault, with two degrees
!   of freedom per node-pair.  At steeper dips, the degree of freedom
!   corresponding to opening or convergence of the opposite sides is
!   eliminated by a constraint equation, and the fault is treated as
!   a vertical strike-slip fault.  This arbitrary limit is necessary
!   because the equations for dip-slip faults become singular as the
!   dip approaches 90 degrees.  In practice, it is best to specify dips
!   as either (1) vertical, or (2) clearly less than dipMax, within
!   each fault element.  If the dip varies within an element in such a
!   way that it passes through this limit within one element, then
!   the representation of that fault element in the equations may
!   be inaccurate.
       data dipMax / 75.0D0 /

! PB2002 plate names [Bird, 2003, G**3, Table 1]:
       data names / 'AF', 'AM', 'AN', &
     &              'AP', 'AR', 'AS', &
     &              'AT', 'AU', 'BH', &
     &              'BR', 'BS', 'BU', &
     &              'CA', 'CL', 'CO', &
     &              'CR', 'EA', 'EU', &
     &              'FT', 'GP', 'IN', &
     &              'JF', 'JZ', 'KE', &
     &              'MA', 'MN', 'MO', &
     &              'MS', 'NA', 'NB', &
     &              'ND', 'NH', 'NI', &
     &              'NZ', 'OK', 'ON', &
     &              'PA', 'PM', 'PS', &
     &              'RI', 'SA', 'SB', &
     &              'SC', 'SL', 'SO', &
     &              'SS', 'SU', 'SW', &
     &              'TI', 'TO', 'WL', &
     &              'YA' /
! Index number of Africa plate in this model.
integer, parameter :: iPAfri = 1

end module

BLOCK data BD1

!   Define "weight" (Gaussian integration weights) of the
!   seven integration points in each element, defined by internal
!   coordinates points(3, 7), where points(1-3, m) holds the s1-s3 of
!   integration point number m.
!   Because all of these arrays are functions of internal
!   coordinates, they are not affected by scaling or shape of
!   particular elements.

implicit none
double precision points, weight
common / S1S2S3 / points
common / WgtVec / weight
dimension points(3, 7), weight(7)

!  "points" contains the internal coordinates (s1, s2, s3) of the 7
!   Gaussian integration points (for area integrals) of the
!   triangular elements.  "points" is also the set of nodal functions
!   for unprojected scalar quantities within an element:
data points / &
& 0.3333333333333333D0, 0.3333333333333333D0, 0.3333333333333333D0, &
& 0.0597158733333333D0, 0.4701420633333333D0, 0.4701420633333333D0, &
& 0.4701420633333333D0, 0.0597158733333333D0, 0.4701420633333333D0, &
& 0.4701420633333333D0, 0.4701420633333333D0, 0.0597158733333333D0, &
& 0.7974269866666667D0, 0.1012865066666667D0, 0.1012865066666667D0, &
& 0.1012865066666667D0, 0.7974269866666667D0, 0.1012865066666667D0, &
& 0.1012865066666667D0, 0.1012865066666667D0, 0.7974269866666667D0 /

!  "weight" is the Gaussian weight (for area integrals) of the 7
!   integration points in each triangular element:
data weight / 0.2250000000000000D0, &
& 0.1323941500000000D0, 0.1323941500000000D0, 0.1323941500000000D0, &
& 0.1259391833333333D0, 0.1259391833333333D0, 0.1259391833333333D0 /

end BLOCK data BD1

BLOCK data BD2

!   Define fPhi (nodal functions) and fGauss (Gaussian integration
!   weights) at the 7 integration points in each fault element,
!   defined by internal coordinate fPoint(m = 1:7),
!   which contains the relative position
!   (fractional length) of the integration points.
!   Because all of these arrays are functions of internal
!   coordinates, they are not affected by length or orientation of
!   particular fault elements.

implicit none
double precision fPhi, fPoint, fGauss
common / SFault / fPoint
common / FPhis /  fPhi
common / FGList / fGauss
dimension fPhi(4, 7), fPoint(7), fGauss(7)

!   fPoint contains the seven integration point locations for the fault
!   elements.  Each value gives a position as a fraction of total length
!   measured from node1 to node2 (of array nodeF):
data fPoint / &
&       0.0254461D0, &
&       0.1292344D0, &
&       0.2970774D0, &
&       0.5000000D0, &
&       0.7029226D0, &
&       0.8707656D0, &
&       0.9745539D0 /

!   fGauss contains the seven corresponding weight factors for use in
!   line integrals:
data fGauss / &
&        0.0647425D0, &
&        0.1398527D0, &
&        0.1909150D0, &
&        0.2089796D0, &
&        0.1909150D0, &
&        0.1398527D0, &
&        0.0647425D0 /

!   fPhi contains the values of the 4 nodal functions (one per node)
!   at each of these 7 integration points in the fault element.
!   A special convention is that the nodal function of node 3
!   is the negative of that for node 2, while the nodal function
!   for node 4 is the negative of that for node 1.  This simplifies
!   many expressions in which we would otherwise have to have
!   a separate factor of +1 or -1 for the two sides of the fault.
data fPhi / &
&      0.9745539D0,  0.0254461D0, -0.0254461D0, -0.9745539D0, &
&      0.8707656D0,  0.1292344D0, -0.1292344D0, -0.8707656D0, &
&      0.7029226D0,  0.2970774D0, -0.2970774D0, -0.7029226D0, &
&      0.5000000D0,  0.5000000D0, -0.5000000D0, -0.5000000D0, &
&      0.2970774D0,  0.7029226D0, -0.7029226D0, -0.2970774D0, &
&      0.1292344D0,  0.8707656D0, -0.8707656D0, -0.1292344D0, &
&      0.0254461D0,  0.9745539D0, -0.9745539D0, -0.0254461D0 /

end BLOCK data BD2