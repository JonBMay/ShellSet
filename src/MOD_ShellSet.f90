!*******************************************************************************
! Module containing ShellSet subroutines
!
! Copyright (C) 2002 Jon May, Peter Bird, Michele Carafa
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

module ShellSetSubs

use SharedVars

implicit none

contains


!------------------------------------------------------------------------------
! Read Input
!------------------------------------------------------------------------------

subroutine InputFileCheck(OData,Abort) ! Check required OrbData, Shells, OrbScore files exist

logical,intent(in) :: OData
logical,intent(out) :: Abort

character(len=100),dimension(6) :: DataFiles
character(len=100),dimension(9) :: ScoreFiles
logical,dimension(7) :: ScoreFilesL
character(len=100),dimension(9) :: ShellsFiles
logical,dimension(4) :: ShellsFilesL
character(len=6) :: dir
integer :: i


Abort = .False.
write(dir,"(A)") "INPUT/"

if(OData) then
  call readDATA(DataFiles)
  do i = 1,size(DataFiles)
    if(.not. FileExist(trim(dir)//DataFiles(i))) then
      Abort = .True.
      print*,"OrbData required file not found: ",trim(dir)//trim(DataFiles(i))
    end if
  end do
end if

call readSHELLS(ShellsFiles,ShellsFilesL)
do i = 1,5
  if(.not. FileExist(dir//ShellsFiles(i))) then
    Abort = .True.
    print*,"Shells required file not found: ",dir//trim(ShellsFiles(i))
  end if
end do

call readSHELLSFINAL(ShellsFiles,ShellsFilesL)
do i = 1,5
  if(.not. FileExist(dir//ShellsFiles(i))) then
    Abort = .True.
    print*,"Shells Final required file not found: ",dir//trim(ShellsFiles(i))
  end if
end do

call readSCORE(ScoreFiles,ScoreFilesL)
do i = 1,2
  if(.not. FileExist(trim(dir)//ScoreFiles(i))) then
    Abort = .True.
    print*,"OrbScore required file not found: ",trim(dir)//trim(ScoreFiles(i))
  end if
end do

end subroutine


subroutine readDATA(DataFiles) ! Read OrbData part of InputFiles.in

character(len=100),dimension(6),intent(out) :: DataFiles

character(len=100) :: Grid,Para,ETOPO,Age,CrustThick,DeltaTS


open(unit=111,file='INPUT/InputFiles.in',status='old',action='read')

  read(111,*) ! Skip first line
  read(111,*) Para
  read(111,*) Grid
  read(111,*) ETOPO
  read(111,*) Age
  read(111,*) CrustThick
  read(111,*) DeltaTS

close(111)

DataFiles = (/Grid,Para,ETOPO,Age,CrustThick,DeltaTS/)

end subroutine


subroutine ReadGridInput(VarNames,VarValsMin,VarValsMax,NumMods,Cells,Levels,TotNumMods,ModRpt,RptMod,Abort) ! Read GridInput.in

integer :: NumVars
character(len=10),dimension(:),allocatable,intent(out) :: VarNames
real*8, dimension(:),allocatable,intent(out) :: VarValsMin
real*8, dimension(:),allocatable,intent(out) :: VarValsMax
integer,dimension(:),allocatable,intent(out) :: NumMods
integer,intent(out) :: Cells, Levels
integer,intent(out) :: TotNumMods
logical,intent(out) :: ModRpt
integer,intent(out) :: RptMod
logical,intent(out) :: Abort

integer :: i,ierr


TotNumMods = 1
ModRpt = .True.
Abort = .False.

open(unit=1,file='INPUT/GridInput.in')
  read(1,*,iostat=ierr) NumVars,Cells,Levels
  if(ierr /= 0) Abort = .True.
  allocate(VarNames(NumVars),VarValsMin(NumVars),VarValsMax(NumVars),NumMods(NumVars))
  read(1,*)
  if(.not. Abort) then
    do i = 1,NumVars
      read(1,*,iostat=ierr) VarNames(i)
      if(ierr /= 0) then
        Abort = .True.
        exit
      end if
      read(1,*,iostat=ierr) VarValsMin(i)
      if(ierr /= 0) then
        Abort = .True.
        exit
      end if
      read(1,*,iostat=ierr) VarValsMax(i)
      if(ierr /= 0) then
        Abort = .True.
        exit
      end if
      read(1,*,iostat=ierr) NumMods(i)
      if(ierr /= 0) then
        Abort = .True.
        exit
      end if
      TotNumMods = TotNumMods*NumMods(i)
      if(mod(NumMods(i),2)==0) ModRpt = .False.
      read(1,*) ! skip -----
    end do
  end if
close(1)

if(.not. Abort) then
! Find repeated model - will always be central model
  if(ModRpt) then
    RptMod = (TotNumMods+1)/2
  else
    RptMod = 0
  end if
end if

do i = 1,size(VarNames)
  if(trim(VarNames(i)) /= 'fFric' .AND. trim(VarNames(i)) /= 'cFric' .AND. trim(VarNames(i)) /= 'Biot' .AND. &
  &  trim(VarNames(i)) /= 'Byerly' .AND. trim(VarNames(i)) /= 'aCreep_C' .AND. trim(VarNames(i)) /= 'aCreep_M' .AND. &
  & trim(VarNames(i)) /= 'bCreep_C' .AND. trim(VarNames(i)) /= 'bCreep_M' .AND. trim(VarNames(i)) /= 'cCreep_C' .AND. &
  & trim(VarNames(i)) /= 'cCreep_M' .AND. trim(VarNames(i)) /= 'dCreep_C' .AND. trim(VarNames(i)) /= 'dCreep_M' .AND. &
  & trim(VarNames(i)) /= 'eCreep' .AND. trim(VarNames(i)) /= 'tAdiab' .AND. trim(VarNames(i)) /= 'gradie' .AND. &
  & trim(VarNames(i)) /= 'zBAsth' .AND. trim(VarNames(i)) /= 'pltRef' .AND. trim(VarNames(i)) /= 'iConve' .AND. &
  & trim(VarNames(i)) /= 'trHMax' .AND. trim(VarNames(i)) /= 'tauMax' .AND. trim(VarNames(i)) /= 'tauMax_S' .AND. &
  & trim(VarNames(i)) /= 'tauMax_L' .AND. trim(VarNames(i)) /= 'rhoH2O' .AND. trim(VarNames(i)) /= 'rhoBar_C' .AND. &
  & trim(VarNames(i)) /= 'rhoBar_M' .AND. trim(VarNames(i)) /= 'rhoAst' .AND. trim(VarNames(i)) /= 'gMean' .AND. &
  & trim(VarNames(i)) /= 'oneKm' .AND. trim(VarNames(i)) /= 'radius' .AND. trim(VarNames(i)) /= 'alphaT_C' .AND. &
  & trim(VarNames(i)) /= 'alphaT_M' .AND. trim(VarNames(i)) /= 'conduc_C' .AND. trim(VarNames(i)) /= 'conduc_M' .AND. &
  & trim(VarNames(i)) /= 'radio_C' .AND. trim(VarNames(i)) /= 'radio_M' .AND. trim(VarNames(i)) /= 'tSurf' .AND. &
  & trim(VarNames(i)) /= 'temLim_C' .AND. trim(VarNames(i)) /= 'temLim_M') then
    print*,'One or more Variable names not recognised'
    Abort = .True.
  end if
end do

if(Abort) print"(A)","Error reading GridInput.in - program will terminate"

end subroutine


subroutine ReadListInput(VarNames,VarValues,Abort) ! Read ListInput.in

character(len=10),dimension(:),allocatable,intent(out) :: VarNames
real*8,dimension(:,:),allocatable,intent(out) :: VarValues
logical,intent(out) :: Abort

integer :: NumVar,NumModels
integer :: i,j,ierr


Abort = .False.
open(unit=1,file="INPUT/ListInput.in")

  read(1,*,iostat=ierr) NumVar,NumModels
  if(ierr /= 0) then
    Abort = .True.
  end if
  allocate(VarNames(NumVar),VarValues(NumModels,NumVar))

! Read Variable names from ListInput.in
  if(.not. Abort) then
    do i = 1,NumVar
      read(1,*,iostat=ierr) VarNames(i)
      if(ierr /= 0) then
        Abort = .True.
        exit
      end if
    end do
    read(1,*) ! to skip ---- in input file
  end if

! Read Variable values from ListInput.in in order
  if(.not. Abort) then
    do i = 1,NumVar
    do j = 1,NumModels
      read(1,*,iostat=ierr) VarValues(j,i)
      if(ierr /= 0) then
        Abort = .True.
        exit
      end if
      if(j==NumModels) read(1,*) ! to skip ---- in input file
    end do
    if(Abort) exit
    end do
  end if
close(1)


do i = 1,size(VarNames)
  if(trim(VarNames(i)) /= 'fFric'    .AND. trim(VarNames(i)) /= 'cFric'    .AND. trim(VarNames(i)) /= 'Biot'     .AND. &
  &  trim(VarNames(i)) /= 'Byerly'   .AND. trim(VarNames(i)) /= 'aCreep_C' .AND. trim(VarNames(i)) /= 'aCreep_M' .AND. &
  &  trim(VarNames(i)) /= 'bCreep_C' .AND. trim(VarNames(i)) /= 'bCreep_M' .AND. trim(VarNames(i)) /= 'cCreep_C' .AND. &
  &  trim(VarNames(i)) /= 'cCreep_M' .AND. trim(VarNames(i)) /= 'dCreep_C' .AND. trim(VarNames(i)) /= 'dCreep_M' .AND. &
  &  trim(VarNames(i)) /= 'eCreep'   .AND. trim(VarNames(i)) /= 'tAdiab'   .AND. trim(VarNames(i)) /= 'gradie'   .AND. &
  &  trim(VarNames(i)) /= 'zBAsth'   .AND. trim(VarNames(i)) /= 'pltRef'   .AND. trim(VarNames(i)) /= 'iConve'   .AND. &
  &  trim(VarNames(i)) /= 'trHMax'   .AND. trim(VarNames(i)) /= 'tauMax'   .AND. trim(VarNames(i)) /= 'tauMax_S' .AND. &
  &  trim(VarNames(i)) /= 'tauMax_L' .AND. trim(VarNames(i)) /= 'rhoH2O'   .AND. trim(VarNames(i)) /= 'rhoBar_C' .AND. &
  &  trim(VarNames(i)) /= 'rhoBar_M' .AND. trim(VarNames(i)) /= 'rhoAst'   .AND. trim(VarNames(i)) /= 'gMean'    .AND. &
  &  trim(VarNames(i)) /= 'oneKm'    .AND. trim(VarNames(i)) /= 'radius'   .AND. trim(VarNames(i)) /= 'alphaT_C' .AND. &
  &  trim(VarNames(i)) /= 'alphaT_M' .AND. trim(VarNames(i)) /= 'conduc_C' .AND. trim(VarNames(i)) /= 'conduc_M' .AND. &
  &  trim(VarNames(i)) /= 'radio_C'  .AND. trim(VarNames(i)) /= 'radio_M'  .AND. trim(VarNames(i)) /= 'tSurf'    .AND. &
  &  trim(VarNames(i)) /= 'temLim_C' .AND. trim(VarNames(i)) /= 'temLim_M') then
    print*,'One or more Variable names not recognised'
    Abort = .True.
  end if
end do

if(Abort) print"(A)","Error reading ListInput.in - program will terminate"

end subroutine


subroutine readSCORE(ScoreFiles,ScoreFilesL,ModNum,MisRun) ! Read OrbScore part of InputFiles.in

character(len=100),dimension(9),intent(out) :: ScoreFiles
logical,dimension(7),intent(out) :: ScoreFilesL
integer,intent(in),optional :: ModNum
logical,dimension(6),intent(out),optional :: MisRun

character(len=100) :: Grid,Para
character(len=100) :: Geo,SeaSprd,StresDir,Fslip,SeisCorre,SeisAniso,Plates
logical :: exists,GeoL,SeaSprdL,StresDirL,FslipL,SeisCorreL,SeisAnisoL,PlatesL
character(len=6) :: dir
logical,save :: called = .False.


GeoL = .False.
SeaSprdL = .False.
StresDirL = .False.
FslipL = .False.
SeisCorreL = .False.
SeisAnisoL = .False.
PlatesL = .False.

dir = 'INPUT/'

open(unit=111,file='INPUT/InputFiles.in',status='old',action='read')

  do i = 1,36
    read(111,*) ! Skip first 36 lines
  end do
  read(111,*) Grid
  read(111,*) Para
  read(111,*) ! ----------
  read(111,*) Geo
  read(111,*) StresDir
  read(111,*) Fslip
  read(111,*) SeaSprd
  read(111,*) SeisCorre
  read(111,*) SeisAniso
  read(111,*) Plates

close(111)


if((Geo .NE. 'X')       .and. (FileExist(trim(dir)//trim(Geo))))       GeoL = .TRUE.
if((StresDir .NE. 'X')  .and. (FileExist(trim(dir)//trim(StresDir))))  StresDirL = .TRUE.
if((Fslip .NE. 'X')     .and. (FileExist(trim(dir)//trim(Fslip))) )    FslipL = .TRUE.
if((SeaSprd .NE. 'X')   .and. (FileExist(trim(dir)//trim(SeaSprd))))   SeaSprdL = .TRUE.
if((SeisCorre .NE. 'X') .and. (FileExist(trim(dir)//trim(SeisCorre)))) SeisCorreL = .TRUE.
if((SeisAniso .NE. 'X') .and. (FileExist(trim(dir)//trim(SeisAniso)))) SeisAnisoL = .TRUE.
if((Plates .NE. 'X')    .and. (FileExist(trim(dir)//trim(Plates))))    PlatesL = .TRUE.


ScoreFiles = (/Grid,Para,Geo,StresDir,Fslip,SeaSprd,SeisCorre,SeisAniso,Plates/)
ScoreFilesL = (/GeoL,StresDirL,FslipL,SeaSprdL,SeisCorreL,SeisAnisoL,PlatesL/)
if(present(MisRun)) MisRun = (/GeoL,SeaSprdL .and. StresDirL,StresDirL,FslipL,SeisCorreL,SeisCorreL .and. SeisAnisoL .and. PlatesL/) ! Must match order in misfit vector output of OrbScore

if((.not. called) .and. present(ModNum)) then
  called = .True.
  do i =1,7
    if((.not. ScoreFilesL(i)) .and. ScoreFiles(i+2) .ne. 'X') then
      write(ErrorMsg,"(3A)") "file not found: ",trim(dir)//trim(ScoreFiles(i+2))
      call NonFatalError(ErrorMsg,ModNum)
    end if
  end do

end if


end subroutine


subroutine readSHELLS(ShellsFiles,ShellsFilesL,ModNum) ! Read Shells part of InputFiles.in

character(len=100),dimension(9),intent(out) :: ShellsFiles
logical,dimension(4),intent(out) :: ShellsFilesL
integer,intent(in),optional :: ModNum

character(len=100) :: Grid,Bound,Para,Plate,PlateB
character(len=100) :: Vel,Mantle,Torque,Litho
logical :: exists,VelL,MantleL,TorqueL,LithoL
character(len=6) :: dir
logical,save :: called = .False.


VelL = .False.
MantleL = .False.
TorqueL = .False.
LithoL = .False.

dir = 'INPUT/'

open(unit=111,file='INPUT/InputFiles.in',status='old',action='read')

  do i = 1,10
    read(111,*) ! Skip first 10 lines
  end do
  read(111,*) Grid
  read(111,*) Bound
  read(111,*) Para
  read(111,*) Plate
  read(111,*) PlateB
  read(111,*) ! ----------
  read(111,*) Vel
  read(111,*) Mantle
  read(111,*) Torque
  read(111,*) Litho

close(111)


if((Vel .NE. 'X')    .and. (FileExist(trim(dir)//trim(Vel))))    VelL = .TRUE.
if((Mantle .NE. 'X') .and. (FileExist(trim(dir)//trim(Mantle)))) MantleL = .TRUE.
if((Torque .NE. 'X') .and. (FileExist(trim(dir)//trim(Torque)))) TorqueL = .TRUE.
if((Litho .NE. 'X')  .and. (FileExist(trim(dir)//trim(Litho))))  LithoL = .TRUE.


ShellsFiles = (/Grid,Bound,Para,Plate,PlateB,Vel,Mantle,Torque,Litho/)
ShellsFilesL = (/VelL,MantleL,TorqueL,LithoL/)

if((.not. called) .and. present(ModNum)) then
  called = .True.
  do i =1,4
    if((.not. ShellsFilesL(i)) .and. ShellsFiles(i+5) .ne. 'X') then
      write(ErrorMsg,"(3A)") "file not found: ",trim(dir)//trim(ShellsFiles(i+5))
      call NonFatalError(ErrorMsg,ModNum)
    end if
  end do

end if

end subroutine


subroutine readSHELLSFINAL(ShellsFiles,ShellsFilesL,ModNum) ! Read Shells final run part of InputFiles.in

character(len=100),dimension(9),intent(out) :: ShellsFiles
logical,dimension(4),intent(out) :: ShellsFilesL
integer,intent(in),optional :: ModNum

character(len=100) :: Grid,Bound,Para,Plate,PlateB
character(len=100) :: Vel,Mantle,Torque,Litho
logical :: exists,VelL,MantleL,TorqueL,LithoL
character(len=6) :: dir
logical,save :: called = .False.


VelL = .False.
MantleL = .False.
TorqueL = .False.
LithoL = .False.

dir = 'INPUT/'

open(unit=111,file='INPUT/InputFiles.in',status='old',action='read')

  do i = 1,23
    read(111,*) ! Skip first 23 lines
  end do
  read(111,*) Grid
  read(111,*) Bound
  read(111,*) Para
  read(111,*) Plate
  read(111,*) PlateB
  read(111,*) ! ----------
  read(111,*) Vel
  read(111,*) Mantle
  read(111,*) Torque
  read(111,*) Litho

close(111)


if((Vel .NE. 'X')    .and. (FileExist(trim(dir)//trim(Vel))))    VelL = .TRUE.
if((Mantle .NE. 'X') .and. (FileExist(trim(dir)//trim(Mantle)))) MantleL = .TRUE.
if((Torque .NE. 'X') .and. (FileExist(trim(dir)//trim(Torque)))) TorqueL = .TRUE.
if((Litho .NE. 'X')  .and. (FileExist(trim(dir)//trim(Litho))))  LithoL = .TRUE.


ShellsFiles = (/Grid,Bound,Para,Plate,PlateB,Vel,Mantle,Torque,Litho/)
ShellsFilesL = (/VelL,MantleL,TorqueL,LithoL/)

if((.not. called) .and. present(ModNum)) then
  called = .True.
  do i =1,4
    if((.not. ShellsFilesL(i)) .and. ShellsFiles(i+5) .ne. 'X') then
      write(ErrorMsg,"(3A)") "file not found: ",trim(dir)//trim(ShellsFiles(i+5))
      call NonFatalError(ErrorMsg,ModNum)
    end if
  end do

end if

end subroutine



!------------------------------------------------------------------------------
! Output
!------------------------------------------------------------------------------

subroutine FormatStrings(Opt,VarNames,NumMis,frmt) ! Create output format

character(len=4),intent(in) :: Opt
character(len=10),dimension(:),intent(in) :: VarNames
integer,intent(in) :: NumMis

character(len=500),intent(out) :: frmt

integer :: i


if(Opt == 'List') frmt = "(I0,X,I0," ! ThID, Model Numbers
if(Opt == 'Grid') frmt = "(I0,X,I0,X,I0,X,I0," ! ModNum, ThID, Level, Cell Numbers
if(Opt == 'MisL') frmt = "(" ! ModNum, ThID, Level, Cell Numbers

do i = 1,size(VarNames)

  select case(VarNames(i))
    case('fFric')
      frmt = trim(frmt)//"X,F7.5,"
    case('cFric')
      frmt = trim(frmt)//"X,F7.5,"
    case('Biot')
      frmt = trim(frmt)//"X,F7.5,"
    case('Byerly')
      frmt = trim(frmt)//"X,F7.5,"
    case('aCreep_C')
      frmt = trim(frmt)//"X,ES12.5,"
    case('aCreep_M')
      frmt = trim(frmt)//"X,ES12.5,"
    case('bCreep_C')
      frmt = trim(frmt)//"X,F11.5,"
    case('bCreep_M')
      frmt = trim(frmt)//"X,F11.5,"
    case('cCreep_C')
      frmt = trim(frmt)//"X,F11.5,"
    case('cCreep_M')
      frmt = trim(frmt)//"X,F11.5,"
    case('dCreep_C')
      frmt = trim(frmt)//"X,ES12.5,"
    case('dCreep_M')
      frmt = trim(frmt)//"X,ES12.5,"
    case('eCreep')
      frmt = trim(frmt)//"X,F12.6,"
    case('tAdiab')
      frmt = trim(frmt)//"X,F12.6,"
    case('gradie')
      frmt = trim(frmt)//"X,ES12.5,"
    case('zBAsth')
      frmt = trim(frmt)//"X,ES12.5,"
    case('trHMax')
      frmt = trim(frmt)//"X,ES12.5,"
    case('tauMax')
      frmt = trim(frmt)//"X,ES12.5,"
    case('tauMax_S')
      frmt = trim(frmt)//"X,ES12.5,"
    case('tauMax_L')
      frmt = trim(frmt)//"X,ES12.5,"
    case('rhoH2O')
      frmt = trim(frmt)//"X,F11.5,"
    case('rhoBar_C')
      frmt = trim(frmt)//"X,F11.5,"
    case('rhoBar_M')
      frmt = trim(frmt)//"X,F11.5,"
    case('rhoAst')
      frmt = trim(frmt)//"X,F11.5,"
    case('gMean')
      frmt = trim(frmt)//"X,F10.5,"
    case('oneKm')
      frmt = trim(frmt)//"X,F10.5,"
    case('radius')
      frmt = trim(frmt)//"X,F15.5,"
    case('alphaT_C')
      frmt = trim(frmt)//"X,ES12.5,"
    case('alphaT_M')
      frmt = trim(frmt)//"X,ES12.5,"
    case('conduc_C')
      frmt = trim(frmt)//"X,F11.5,"
    case('conduc_M')
      frmt = trim(frmt)//"X,F11.5,"
    case('radio_C')
      frmt = trim(frmt)//"X,ES12.5,"
    case('radio_M')
      frmt = trim(frmt)//"X,ES12.5,"
    case('tSurf')
      frmt = trim(frmt)//"X,F11.5,"
    case('temLim_C')
      frmt = trim(frmt)//"X,F11.5,"
    case('temLim_M')
      frmt = trim(frmt)//"X,F11.5,"
    case default
      print*,'One or more Variable Names not recognised inside FormatStrings'
  end select

end do


select case(NumMis)
  case(0)
    frmt = trim(frmt)//"X,A"
  case(1)
    frmt = trim(frmt)//"X,F12.5"
  case(2)
    frmt = trim(frmt)//"X,F12.5,X,F12.5"
  case(3)
    frmt = trim(frmt)//"X,F12.5,X,F12.5,X,F12.5"
  case(4)
    frmt = trim(frmt)//"X,F12.5,X,F12.5,X,F12.5,X,F12.5"
  case(5)
    frmt = trim(frmt)//"X,F12.5,X,F12.5,X,F12.5,X,F12.5,X,F12.5"
  case(6)
    frmt = trim(frmt)//"X,F12.5,X,F12.5,X,F12.5,X,F12.5,X,F12.5,X,F12.5"
  case default
    print*,'ERROR in number of misfits'
end select


frmt = trim(frmt)//")"

end subroutine


subroutine WriteMisfitTarget(AllModels,VarNames,TargMisM,LvlNumMods,Levels,Step,DirName,MisType,iLevel) ! Write to Models_Lim.txt

real*8,dimension(:,:),intent(in) :: AllModels
character(len=10),dimension(:),intent(in) :: VarNames
real,intent(in) :: TargMisM
integer,intent(in):: LvlNumMods,Levels
real*8,dimension(:),intent(in) :: Step
character(len=100),intent(in) :: DirName
character(len=3),intent(in) :: MisType
integer,intent(in),optional :: iLevel

character(len=500) :: frmt
character(len=2) :: SizeStep
integer :: i,j
logical,save :: Called = .False.

if(.not. Called) then
  open(1,file=trim(DirName)//'/Models_Lim.txt',status='new')
  Called = .True.
else
  open(1,file=trim(DirName)//'/Models_Lim.txt',access='append',status='old')
end if

frmt = "(A,I0,A,I0,A,I0,A,"
write(SizeStep,'(I0)'), size(Step)
frmt = trim(frmt)//trim(SizeStep)//"ES12.5"
frmt = trim(frmt)//")"
write(1,trim(frmt)) 'Level ', iLevel ,' of ',Levels,' with ',LvlNumMods,' models and Cell sizes +-:',Step/2.0

call FormatStrings('MisL',VarNames,1,frmt)

if(MisType == 'SC') then
  do i = 1,size(AllModels,1)
    if(AllModels(i,size(AllModels,2)) >= TargMisM) write(1,frmt) (AllModels(i,j), j=1,size(AllModels,2))
  end do
else
  do i = 1,size(AllModels,1)
    if(AllModels(i,size(AllModels,2)) <= TargMisM) write(1,frmt) (AllModels(i,j), j=1,size(AllModels,2))
  end do
end if

write(1,*)

close(1)

end subroutine


subroutine WriteOut(Opt,VarNames,msg_tag,msg_src,iLevel,iCell,VarValues,misfits,MisRun,MisType,Failed) ! Write to Models.txt

character(len=4),intent(in) :: Opt
character(len=10),dimension(:),intent(in) :: VarNames
integer,intent(in) :: msg_tag,msg_src,iLevel,iCell
real*8,dimension(:),intent(in) :: VarValues,misfits
logical,dimension(:),intent(in) :: MisRun
character(len=3),intent(in),optional :: MisType
integer,intent(in),optional :: Failed

real*8,dimension(:),allocatable :: Non0Mis
character(len=500) :: frmt
integer :: i,j
logical,dimension(6) :: mask
character(len=100) :: outhdr
integer :: NumMis
logical,save :: Called = .False.


if(.not. Called) then
  if(Opt == 'Grid') then
    write(outhdr,'(A,I0,A)') "global model, ThID, Level, Cell, VarValues(",size(VarNames),")"
  elseif(Opt == 'List') then
    write(outhdr,'(A,I0,A)') "global model, ThID, VarValues(",size(VarNames),")"
  end if
end if

NumMis = count(misfits/=0)
call FormatStrings(Opt,VarNames,NumMis,frmt)
allocate(Non0Mis(NumMis))

if(opt == 'Grid') then

  mask = .True.

  if(trim(MisType) == 'GV') then
    Non0Mis(1) = misfits(1)
    mask(1) = .False.
    outhdr = trim(outhdr)//", GV"
  elseif(trim(MisType) == 'SSR') then
    Non0Mis(1) = misfits(2)
    mask(2) = .False.
    outhdr = trim(outhdr)//", SSR"
  elseif(trim(MisType) == 'SD') then
    Non0Mis(1) = misfits(3)
    mask(3) = .False.
    outhdr = trim(outhdr)//", SD"
  elseif(trim(MisType) == 'FSR') then
    Non0Mis(1) = misfits(4)
    mask(4) = .False.
    outhdr = trim(outhdr)//", FSR"
  elseif(trim(MisType) == 'SC') then
    Non0Mis(1) = misfits(5)
    mask(5) = .False.
    outhdr = trim(outhdr)//", SC"
  elseif(trim(MisType) == 'SA') then
    Non0Mis(1) = misfits(6)
    mask(6) = .False.
    outhdr = trim(outhdr)//", SA"
  end if


  j = 2
  do i = 1,size(misfits)
    if(misfits(i) > 0.0 .and. MisRun(i) .and. mask(i)) then
      Non0Mis(j) = misfits(i)
      mask(i) = .False.
      j=j+1
      if(.not. Called) then
        if(i == 1) outhdr = trim(outhdr)//", GV"
        if(i == 2) outhdr = trim(outhdr)//", SSR"
        if(i == 3) outhdr = trim(outhdr)//", SD"
        if(i == 4) outhdr = trim(outhdr)//", FSR"
        if(i == 5) outhdr = trim(outhdr)//", SC"
        if(i == 6) outhdr = trim(outhdr)//", SA"
      end if
    end if
  end do

elseif(opt == 'List') then

  j = 1
  do i = 1,size(misfits)
    if(misfits(i) > 0.0 .and. MisRun(i)) then
      Non0Mis(j) = misfits(i)
      j=j+1
      if(.not. Called) then
        if(i == 1) outhdr = trim(outhdr)//", GV"
        if(i == 2) outhdr = trim(outhdr)//", SSR"
        if(i == 3) outhdr = trim(outhdr)//", SD"
        if(i == 4) outhdr = trim(outhdr)//", FSR"
        if(i == 5) outhdr = trim(outhdr)//", SC"
        if(i == 6) outhdr = trim(outhdr)//", SA"
      end if
    end if
  end do
end if

if(.not. Called) then
  write(101,'(A)') outhdr
end if

if(.not. present(Failed)) then
  if(Opt == 'Grid') then
    write(101,trim(frmt)) msg_tag,msg_src,iLevel,iCell,VarValues,Non0Mis
  elseif(Opt == 'List') then
    write(101,trim(frmt)) msg_tag,msg_src,VarValues,Non0Mis
  end if

elseif(present(Failed)) then

  call FormatStrings(Opt,VarNames,0,frmt)

  if(Failed == 100000) then! Shells Conv
    if(Opt == 'Grid') then
      write(*,trim(frmt)) msg_tag,msg_src,iLevel,iCell,VarValues,'Failed to converge in required MKL iterations'
      write(101,trim(frmt)) msg_tag,msg_src,iLevel,iCell,VarValues,'Failed to converge in required MKL iterations'
    elseif(Opt == 'List') then
      write(*,trim(frmt)) msg_tag,msg_src,VarValues,'Failed to converge in required MKL iterations'
      write(101,trim(frmt)) msg_tag,msg_src,VarValues,'Failed to converge in required MKL iterations'
    end if

  elseif(Failed == 100001) then! Var check
    if(Opt == 'Grid') then
      write(101,trim(frmt)) msg_tag,msg_src,iLevel,iCell,VarValues,'Failed input variable check'
    elseif(Opt == 'List') then
      write(101,trim(frmt)) msg_tag,msg_src,VarValues,'Failed input variable check'
    end if
  end if

end if

Called = .True.

end subroutine



!------------------------------------------------------------------------------
! Error Handling
!------------------------------------------------------------------------------

subroutine abort(errnum) ! Error messages

integer,intent(in) :: errnum


select case(errnum)

  case(0)
    print*,"Called with -np 1, program must be invoked with >1 threads"

  case(1)
    print'(A)',"Problem with 1 or more CLA inputs - these are case sensitive"
    print"(A/,A/,A/,A/,A/,A/)","The following CLAs print information to the terminal without running the program:", &
    & "-help  -> print information currently shown", &
    & "-info  -> print general information about the program", &
    & "-cite  -> print information about citing this program", &
    & "-copy  -> print information about program copyright", &
    & "-abort -> print information about MPI_Abort error numbers"
    print"(A/,A/,A/)","The following CLAs are required for the program to function:", &
    & "-Iter  -> number of iterations of 'main wirk loop'", &
    & "-InOpt -> either 'List' or 'Grid' to determine model input generation"
    print"(A/,A/,A/,A/,A/,A/,A/,A/)","The following CLAs are optional and add functionality:", &
    & "-MC  -> information for misfit convergence", &
    & "-MT  -> misfit type used to drive grid search", &
    & "-Dir -> new working directory", &
    & "-AEF -> make all errors fatal", &
    & "-V   -> create Verbose file containing extra informational output", &
    & "-ML  -> misfit limit beneath which all models are stored in an extra output file", &
    & "-KL  -> misfit limit at which the program is ended"

  case(2)
    print*,"Problem reading ListInput.in or GridInput.in"

  case(3)
    print*,"FatalError.txt file found"

  case(4)
    print*,"All errors are fatal and either ModelError.txt or FatalError.txt file found"

  case(5)
    print*,"Program initiated with more threads than ever needed"

  case(6)
    print*,"Grid search driver or misfit convergence type not run,"
	print*,"check flag -MT or -MC and ensure file exists with correct name"

  case(7)
    print*,"One or more required input files were not found"

  case(8)
    print*,"Non-running CLA"

  case(9)
    print*,"Error when reading CLAs, see previous error message"

  case default
    print"(A/,A/,A/,A/,A/,A/,A/,A/,A/,A/,A/)","The following list details the error numbers given inside an MPI_Abort program failure:", &
    &   "0 -> Called with -np 1, program must be involed with >1 threads", &
    &   "1  -> CLA problem", &
    &   "2  -> Problem reading ListInput.in or GridInput.in", &
    &   "3  -> FatalError.txt file found", &
    &   "4  -> All errors are fatal and either ModelError.txt or FatalError.txt file found", &
    &   "5  -> Program initiated with more threads than ever needed", &
    &   "6  -> Grid search program driving misfit is either not recognised or will not be run by OrbScore", &
    &   "7  -> One or more required input files were not found", &
    &   "8  -> Non-running CLA performed", &
    &   "9  -> Error when reading CLAs"

end select

end subroutine


subroutine FatalError(msg,ThID,ErrArr,ErrArrCh,ErrArrInt) ! Report fatal errors

character(len=*),intent(in) :: msg
integer,intent(in) :: ThID
real,dimension(:),intent(in),optional :: ErrArr
character(len=2),dimension(:),intent(in),optional :: ErrArrCh
integer,dimension(:),intent(in),optional :: ErrArrInt

character(len=100) :: Filename
integer :: i


call execute_command_line('mkdir -p Error')

! Personal Fatal Error file
write(Filename,'(A,I0,A)') 'Error/FatalError_',ThID,'.txt'
if(FileExist(Filename)) then
  open(unit=1,file=trim(Filename),status='old',access='append')
else
  open(unit=1,file=trim(Filename),status='new',action='write')
end if

    write(1,'(A,A)') ' Model failed with the following message: ',trim(msg)
    if(present(ErrArr)) then
      write(1,*) (ErrArr(i), i=1,size(ErrArr))
    elseif(present(ErrArrCh)) then
      write(1,*) (ErrArrCh(i), i=1,size(ErrArrCh))
    elseif(present(ErrArrInt)) then
      write(1,*) (ErrArrInt(i), i=1,size(ErrArrInt))
    end if
close(1)

! Global Fatal Error file
Filename = 'FatalError.txt'
if(.not. FileExist(Filename)) then
  open(unit=1,file=trim(Filename))
  close(1)
end if

end subroutine


subroutine ModelError(msg,ModNum) ! Report model errors

character(len=*),intent(in) :: msg
integer,intent(in) :: ModNum

character(len=100) :: Filename


call execute_command_line('mkdir -p Error')

! Personal Model Error file
write(Filename,'(A,I0,A)') 'Error/ModelError_',ModNum,'.txt'
if(FileExist(Filename)) then
  open(unit=1,file=trim(Filename),status='old',access='append')
else
  open(unit=1,file=trim(Filename),status='new',action='write')
end if

write(1,'(A)') trim(msg)
close(1)

! Global Model Error file
Filename = 'ModelError.txt'
if(.not. FileExist(Filename)) then
  open(unit=1,file=trim(Filename))
  close(1)
end if

end subroutine


subroutine NonFatalError(msg,ModNum) ! Report non-fatal errors

character(len=*),intent(in) :: msg
integer,intent(in) :: ModNum
character(len=100) :: Filename


call execute_command_line('mkdir -p Error')

! Personal Non-Fatal Error file
write(Filename,'(A,I0,A)') 'Error/NonFatalError_',ThID,'.txt'
if(FileExist(Filename)) then
  open(unit=1,file=trim(Filename),status='old',access='append')
else
  open(unit=1,file=trim(Filename),status='new',action='write')
end if

write(1,'(A,I0,A/,A)') 'Model ',ModNum,' gave the following message: ',trim(msg)
close(1)


! Global NonFatalError Error file
Filename = 'NonFatalError.txt'
if(.not. FileExist(Filename)) then
  open(unit=1,file=trim(Filename))
  close(1)
end if


end subroutine



!------------------------------------------------------------------------------
! Other
!------------------------------------------------------------------------------

subroutine cpuInfo(cpus) ! Find number of cpus

integer,intent(out) :: cpus
integer :: iostat
character(len=100) :: line
integer :: num


open(unit=1,file='/proc/cpuinfo',iostat=iostat)

  do while(iostat==0)
    read(1,"(A)",iostat=iostat) line
    if(line(1:8) == 'siblings') then
      read(line(index(line, ':')+1 : index(line, ':')+4),*) num
      cpus = num
      exit
    elseif(line(1:9) == 'cpu cores') then
      read(line(index(line, ':')+1 : index(line, ':')+4),*) num
      cpus = num*2
      exit
    end if
  end do
close(1)

end subroutine


logical function FileExist(Filename) ! Check if a file exists

character(len=*),intent(in) :: Filename

inquire(file=trim(Filename),exist=FileExist)

end function


subroutine OrbDataCheck(VarNames,OData) ! Check if OrbData should run

character(len=10),dimension(:),intent(in) :: VarNames
logical,intent(out) :: OData
integer :: i


i=1
OData = .False.

do while(.not. OData .and. i <= size(VarNames))
  if(trim(VarNames(i)) == 'tAdiab')   OData = .True.
  if(trim(VarNames(i)) == 'gradie')   OData = .True.
  if(trim(VarNames(i)) == 'zBAsth')   OData = .True.
  if(trim(VarNames(i)) == 'trHMax')   OData = .True.
  if(trim(VarNames(i)) == 'rhoH2O')   OData = .True.
  if(trim(VarNames(i)) == 'rhoBar_C') OData = .True.
  if(trim(VarNames(i)) == 'rhoBar_M') OData = .True.
  if(trim(VarNames(i)) == 'rhoAst')   OData = .True.
  if(trim(VarNames(i)) == 'alphaT_C') OData = .True.
  if(trim(VarNames(i)) == 'alphaT_M') OData = .True.
  if(trim(VarNames(i)) == 'conduc_C') OData = .True.
  if(trim(VarNames(i)) == 'conduc_M') OData = .True.
  if(trim(VarNames(i)) == 'radio_C')  OData = .True.
  if(trim(VarNames(i)) == 'radio_M')  OData = .True.
  if(trim(VarNames(i)) == 'tSurf')    OData = .True.
  if(trim(VarNames(i)) == 'temLim_C') OData = .True.
  if(trim(VarNames(i)) == 'temLim_M') OData = .True.
  i = i+1
end do

end subroutine


subroutine Variable_Update(fFric   , cFric , Biot  , Byerly, aCreep, & ! Update variable values with user input values
                        &  bCreep  , cCreep, dCreep, eCreep, tAdiab, &
                        &  gradie  , zBAsth, trHMax, tauMax, rhoH2O, &
                        &  rhoBar  , rhoAst, gMean , oneKm , radius, &
                        &  alphaT  , conduc, radio , tSurf , temLim, &
                        &  VarNames,VarValues)

real*8,intent(inout) :: alphaT(2) , conduc(2) , fFric  , cFric  , Biot   , Byerly , &
                     &  aCreep(2) , bCreep(2) , eCreep , gMean  , gradie , oneKm   , &
                     &  cCreep(2) , dCreep(2) , radius , rhoAst , rhoH2O , tAdiab  , &
                     &  rhoBar(2) , temLim(2) , trHMax , tSurf  , zBAsth  , &
                     &  tauMax(2) , radio(2)

character(len=10),dimension(:),intent(in) :: VarNames
real*8,dimension(:),intent(in) :: VarValues

character(len=500) :: frmt


do i = 1,size(VarNames)

  select case(VarNames(i))

    case('fFric')
      fFric = VarValues(i)
    case('cFric')
      cFric = VarValues(i)
    case('Biot')
      Biot = VarValues(i)
    case('Byerly')
      Byerly = VarValues(i)
    case('aCreep_C')
      aCreep(1) = VarValues(i)
    case('aCreep_M')
      aCreep(2) = VarValues(i)
    case('bCreep_C')
      bCreep(1) = VarValues(i)
    case('bCreep_M')
      bCreep(2) = VarValues(i)
    case('cCreep_C')
      cCreep(1) = VarValues(i)
    case('cCreep_M')
      cCreep(2) = VarValues(i)
    case('dCreep_C')
      dCreep(1) = VarValues(i)
    case('dCreep_M')
      dCreep(2) = VarValues(i)
    case('eCreep')
      eCreep = VarValues(i)
    case('tAdiab')
      tAdiab = VarValues(i)
    case('gradie')
      gradie = VarValues(i)
    case('zBAsth')
      zBAsth = VarValues(i)
    case('trHMax')
      trHMax = VarValues(i)
    case('tauMax')
      tauMax(1) = VarValues(i)
      tauMax(2) = VarValues(i)
    case('tauMax_S')
      tauMax(1) = VarValues(i)
    case('tauMax_L')
      tauMax(2) = VarValues(i)
    case('rhoH2O')
      rhoH2O = VarValues(i)
    case('rhoBar_C')
      rhoBar(1) = VarValues(i)
    case('rhoBar_M')
      rhoBar(2) = VarValues(i)
    case('rhoAst')
      rhoAst = VarValues(i)
    case('gMean')
      gMean  = VarValues(i)
    case('oneKm')
      oneKm  = VarValues(i)
    case('radius')
      radius = VarValues(i)
    case('alphaT_C')
      alphaT(1) = VarValues(i)
    case('alphaT_M')
      alphaT(2) = VarValues(i)
    case('conduc_C')
      conduc(1) = VarValues(i)
    case('conduc_M')
      conduc(2) = VarValues(i)
    case('radio_C')
      radio(1) = VarValues(i)
    case('radio_M')
      radio(2) = VarValues(i)
    case('tSurf')
      tSurf = VarValues(i)
    case('temLim_C')
      temLim(1) = VarValues(i)
    case('temLim_M')
      temLim(2) = VarValues(i)
    case default
      print*,'One or more Variable Names not recognised inside Variable_Update'
  end select
end do

if(Verbose) then  
  write(iUnitVerb,"(/A)")"The following variables were updated before being used by any of OrbData, Shells or OrbScore"
  write(frmt,"(A,I0,A)") "(",size(VarNames),"(A,X),/)"
  write(iUnitVerb,frmt), (trim(VarNames(i)), i =1,size(VarNames))

  write(iUnitVerb,"(A)") "With the following values:"
  call FormatStrings('MisL',VarNames,0,frmt)
  write(iUnitVerb,frmt), (VarValues(i), i =1,size(VarValues)),' '
  write(iUnitVerb,"(A/)"), "-------------------------------------------------------------------------------"
end if

end subroutine



!------------------------------------------------------------------------------
! File Handling
!------------------------------------------------------------------------------

subroutine CloseInput(prog) ! Close OrbData, Shells, OrbScore input files

character(len=2),intent(in) :: prog ! SH=Shells, OD=Data, OS=Score


select case(prog)

  case("OD")

    close(1)
    close(2)
    close(3)
    close(7)
    close(11)
    close(12)


  case("SH")

    close(1)
    close(2)
    close(3)
    close(8)
    close(9)
    close(13)


  case("OS")

    close(1)
    close(2)
    close(3)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(33)
    close(34)
    close(35)

end select

end subroutine


subroutine CloseOutput(ThID,ModNum,prog,ListDIR,rpeat) ! Close OrbData, Shells, OrbScore output files

character(len=2),intent(in) :: prog ! SH=Shells, OD=Data, OS=Score
character(len=100),intent(in) :: ListDIR
integer,intent(in),optional :: ThID,ModNum
character(len=200) :: outfile,infile
character(len=100) :: outdir,indir
integer,intent(in),optional :: rpeat


select case(prog)

  case("OD")

    close(13)
    close(14)
! OrbData needs to copy the new FEG file for SHELLS & OrbScore to use as input
    write(outdir,"(A,'/','ThID_',I0,'_Data_output/')") trim(ListDIR),ThID
    write(outfile,"(A,'FEG_',i0,'.14')") trim(outdir),ModNum
! FEG for SHELLS
    write(indir,"(A,'/','ThID_',I0,'_Shells_input/')") trim(ListDIR),ThID
    write(infile,"(A,'fort_',i0,'.1')") trim(indir),ModNum
    call execute_command_line( "cp "//trim(outfile)//" "//trim(infile) )
! FEG for OrbScore
    write(indir,"(A,'/','ThID_',I0,'_Score_input/')") trim(ListDIR),ThID
    write(infile,"(A,'fort_',i0,'.1')") trim(indir),ModNum
    call execute_command_line( "cp "//trim(outfile)//" "//trim(infile) )


  case("SH")

    close(21)
    close(22)
    close(23)
    close(24)
! SHELLS needs to copy the generated vEarth file for OrbScore to use as input
    write(outdir,"(A,'/','ThID_',I0,'_Shells_output/')") trim(ListDIR),ThID
    write(indir,"(A,'/','ThID_',I0,'_Score_input/')") trim(ListDIR),ThID
    write(outfile,"(A,'vEarth_',i0,'_',i0,'.22')") trim(outdir),ModNum,rpeat
    write(infile,"(A,'fort_',i0,'.3')") trim(indir),ModNum
    call execute_command_line("cp "//trim(outfile)//" "//trim(infile))


  case("OS")

    close(31)
    close(32)
    close(37)
    close(21)
    close(22)

end select

end subroutine


subroutine InputSetup(ThID,ModNum,prog,ListDIR,MisRun,Verbose) ! Copy OrbData, Shells, OrbScore input files

integer,intent(in) :: ThID,ModNum
character(len=2),intent(in) :: prog ! SH=Shells, OD=Data, OS=Score
character(len=100),intent(in) :: ListDIR
logical,dimension(6),intent(inout),optional :: MisRun
logical,intent(in),optional :: Verbose

integer,parameter :: iUnitVerb = 6 ! Unit number for optional verbose.txt file
character(len=100):: dir,filename
logical :: exists
character(len=100) :: ShellsFiles(9),ShellsFilesFinal(9),ScoreFiles(9),DataFiles(6)
logical :: ShellsFilesL(4),ShellsFilesFinalL(4),ScoreFilesL(7)
logical,save :: callOD=.False.,callSH=.False.,callOS=.False.


select case(prog)

  case("OD")

    call readDATA(DataFiles) ! (/Grid,Para,ETOPO,Age,CrustThick,DeltaTS/)

    write(dir,"(A,'/','ThID_',I0,'_Data_input')") trim(ListDIR),ThID
    if(.not. callOD) then
      call execute_command_line('mkdir '//dir)
      callOD = .True.
    end if
    write(filename,"('fort_',I0,'.')") ModNum

    call execute_command_line('cp INPUT/'//trim(DataFiles(2)) //' '//trim(dir)//'/'//trim(filename)//'1')
    call execute_command_line('cp INPUT/'//trim(DataFiles(1)) //' '//trim(dir)//'/'//trim(filename)//'2')
    call execute_command_line('cp INPUT/'//trim(DataFiles(4)) //' '//trim(dir)//'/'//trim(filename)//'7')
    call execute_command_line('cp INPUT/'//trim(DataFiles(5)) //' '//trim(dir)//'/'//trim(filename)//'11')
    call execute_command_line('cp INPUT/'//trim(DataFiles(6)) //' '//trim(dir)//'/'//trim(filename)//'12')


  case("SH")

    call readSHELLS(ShellsFiles,ShellsFilesL,ModNum=ModNum) ! (/Grid,Bound,Para,Plate,PlateB,Vel,Mantle,Torque,Litho/)

    write(dir,"(A,'/','ThID_',I0,'_Shells_input')") trim(ListDIR),ThID
    if(.not. callSH) then
      call execute_command_line('mkdir '//dir)
      callSH = .True.
    end if
    write(filename,"('fort_',I0,'.')") ModNum

    call execute_command_line('cp INPUT/'//trim(ShellsFiles(1)) //' '//trim(dir)//'/'//trim(filename)//'1')
    call execute_command_line('cp INPUT/'//trim(ShellsFiles(2)) //' '//trim(dir)//'/'//trim(filename)//'2')
    call execute_command_line('cp INPUT/'//trim(ShellsFiles(3)) //' '//trim(dir)//'/'//trim(filename)//'3')
    call execute_command_line('cp INPUT/'//trim(ShellsFiles(4)) //' '//trim(dir)//'/'//trim(filename)//'8')
    call execute_command_line('cp INPUT/'//trim(ShellsFiles(5)) //' '//trim(dir)//'/'//trim(filename)//'9')

    if(ShellsFilesL(1)) call execute_command_line('cp INPUT/'//trim(ShellsFiles(6)) //' '//trim(dir)//'/'//trim(filename)//'11')
    if(ShellsFilesL(2)) call execute_command_line('cp INPUT/'//trim(ShellsFiles(7)) //' '//trim(dir)//'/'//trim(filename)//'12')
    if(ShellsFilesL(3)) call execute_command_line('cp INPUT/'//trim(ShellsFiles(8)) //' '//trim(dir)//'/'//trim(filename)//'13')
    if(ShellsFilesL(4)) call execute_command_line('cp INPUT/'//trim(ShellsFiles(9)) //' '//trim(dir)//'/'//trim(filename)//'14')


  case("OS")

    call readSCORE(ScoreFiles,ScoreFilesL,ModNum=ModNum,MisRun=MisRun) ! (/Grid,Para,Geo,StresDir,Fslip,SeaSprd,SeisCorre,SeisAniso/)

    write(dir,"(A,'/','ThID_',I0,'_Score_input')") trim(ListDIR),ThID
    if(.not. callOS) then
      call execute_command_line('mkdir '//dir)
      callOS = .True.
    end if
    write(filename,"('fort_',I0,'.')") ModNum

    call execute_command_line('cp INPUT/'//trim(ScoreFiles(1))  //' '//trim(dir)//'/'//trim(filename)//'1')
    call execute_command_line('cp INPUT/'//trim(ScoreFiles(2))  //' '//trim(dir)//'/'//trim(filename)//'2')

    if(ScoreFilesL(1)) call execute_command_line('cp INPUT/'//trim(ScoreFiles(3))  //' '//trim(dir)//'/'//trim(filename)//'11')
    if(ScoreFilesL(2)) call execute_command_line('cp INPUT/'//trim(ScoreFiles(4))  //' '//trim(dir)//'/'//trim(filename)//'12')
    if(ScoreFilesL(3)) call execute_command_line('cp INPUT/'//trim(ScoreFiles(5))  //' '//trim(dir)//'/'//trim(filename)//'13')
    if(ScoreFilesL(4)) call execute_command_line('cp INPUT/'//trim(ScoreFiles(6))  //' '//trim(dir)//'/'//trim(filename)//'14')
    if(ScoreFilesL(5)) call execute_command_line('cp INPUT/'//trim(ScoreFiles(7))  //' '//trim(dir)//'/'//trim(filename)//'15')
    if(ScoreFilesL(5) .AND. ScoreFilesL(6) .AND. ScoreFilesL(7)) then
      call execute_command_line('cp INPUT/'//trim(ScoreFiles(8))  //' '//trim(dir)//'/'//trim(filename)//'33')
      call execute_command_line('cp INPUT/'//trim(ScoreFiles(9)) //' '//trim(dir)//'/'//trim(filename)//'35')
    end if


  case("BF")

    call readSHELLSFINAL(ShellsFilesFinal,ShellsFilesFinalL,ModNum=ModNum) ! (/Grid,Bound,Para,Plate,PlateB,Vel,Mantle,Torque,Litho/)
    call readSHELLS(ShellsFiles,ShellsFilesL)

    write(dir,"(A,'/','ThID_',I0,'_Shells_input')") trim(ListDIR),ThID
    write(filename,"('fort_',I0,'.')") ModNum

    if(ShellsFiles(1) /= ShellsFilesFinal(1)) then
      if(Verbose) write(iUnitVerb,*) 'Unit 1 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(1)) //' '//trim(dir)//'/'//trim(filename)//'1')
    end if
    if(ShellsFiles(2) /= ShellsFilesFinal(2)) then
      if(Verbose) write(iUnitVerb,*) 'Unit 2 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(2)) //' '//trim(dir)//'/'//trim(filename)//'2')
    end if
    if(ShellsFiles(3) /= ShellsFilesFinal(3)) then
      if(Verbose) write(iUnitVerb,*) 'Unit 3 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(3)) //' '//trim(dir)//'/'//trim(filename)//'3')
    end if
    if(ShellsFiles(4) /= ShellsFilesFinal(4)) then
      if(Verbose) write(iUnitVerb,*) 'Unit 8 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(4)) //' '//trim(dir)//'/'//trim(filename)//'8')
    end if
    if(ShellsFiles(5) /= ShellsFilesFinal(5)) then
      if(Verbose) write(iUnitVerb,*) 'Unit 9 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(5)) //' '//trim(dir)//'/'//trim(filename)//'9')
    end if

    if(ShellsFilesFinalL(1) .AND. (ShellsFiles(6) /= ShellsFilesFinal(6))) then
      if(Verbose) write(iUnitVerb,*) 'Unit 11 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(6)) //' '//trim(dir)//'/'//trim(filename)//'11')
    end if
    if(ShellsFilesFinalL(2) .AND. (ShellsFiles(7) /= ShellsFilesFinal(7))) then
      if(Verbose) write(iUnitVerb,*) 'Unit 12 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(7)) //' '//trim(dir)//'/'//trim(filename)//'12')
    end if
    if(ShellsFilesFinalL(3) .AND. (ShellsFiles(8) /= ShellsFilesFinal(8))) then
      if(Verbose) write(iUnitVerb,*) 'Unit 13 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(8)) //' '//trim(dir)//'/'//trim(filename)//'13')
    end if
    if(ShellsFilesFinalL(4) .AND. (ShellsFiles(9) /= ShellsFilesFinal(9))) then
      if(Verbose) write(iUnitVerb,*) 'Unit 14 for SHELLS and SHELLS final are not the same'
      call execute_command_line('cp INPUT/'//trim(ShellsFilesFinal(9)) //' '//trim(dir)//'/'//trim(filename)//'14')
    end if

end select

end subroutine


subroutine OpenInput(ThID,ModNum,prog,ListDIR,rpeat) ! Open OrbData, Shells, OrbScore input files

integer,intent(in) :: ThID,ModNum
character(len=2),intent(in) :: prog ! SH=Shells, OD=Data, OS=Score
character(len=100),intent(in) :: ListDIR

character(len=200) :: filename
character(len=100) :: dir
integer,intent(in),optional :: rpeat
logical :: exists


select case(prog)

  case("OD")

    write(dir,"(A,'/','ThID_',I0,'_Data_input/')") trim(ListDIR),ThID
    write(filename,"('fort_',I0,'.')") ModNum
    open(unit=1,file=trim(dir)//trim(filename)//'1')
    open(unit=2,file=trim(dir)//trim(filename)//'2')
    open(unit=3,file=trim(dir)//trim(filename)//'3')
    open(unit=7,file=trim(dir)//trim(filename)//'7')
    open(unit=11,file=trim(dir)//trim(filename)//'11')
    open(unit=12,file=trim(dir)//trim(filename)//'12')


  case("SH")
    write(dir,"(A,'/','ThID_',I0,'_Shells_input/')") trim(ListDIR),ThID
    write(filename,"('fort_',I0,'.')") ModNum
    open(unit=1,file=trim(dir)//trim(filename)//'1')
    open(unit=2,file=trim(dir)//trim(filename)//'2')
    open(unit=3,file=trim(dir)//trim(filename)//'3')
    open(unit=8,file=trim(dir)//trim(filename)//'8')
    open(unit=9,file=trim(dir)//trim(filename)//'9')

    if(FileExist(trim(dir)//trim(filename)//'11')) open(unit=11,file=trim(dir)//trim(filename)//'11')
    if(FileExist(trim(dir)//trim(filename)//'12')) open(unit=12,file=trim(dir)//trim(filename)//'12')
    if(FileExist(trim(dir)//trim(filename)//'13')) open(unit=13,file=trim(dir)//trim(filename)//'13')
    if(FileExist(trim(dir)//trim(filename)//'14')) open(unit=14,file=trim(dir)//trim(filename)//'14')

    if(rpeat>1) then ! use old qEarth file as input
      write(dir,"(A,'/','ThID_',I0,'_Shells_output/')") trim(ListDIR),ThID
      write(filename,"('qEarth_',I0,'_',I0,'.')") ModNum,rpeat-1
      open(unit=13,file=trim(dir)//trim(filename)//'24')
  end if


  case("OS")

    write(dir,"(A,'/','ThID_',I0,'_Score_input/')") trim(ListDIR),ThID
    write(filename,"('fort_',I0,'.')") ModNum
    open(unit=1,file=trim(dir)//trim(filename)//'1')
    open(unit=2,file=trim(dir)//trim(filename)//'2')
    open(unit=3,file=trim(dir)//trim(filename)//'3')

    if(FileExist(trim(dir)//trim(filename)//'11')) open(unit=11,file=trim(dir)//trim(filename)//'11')
    if(FileExist(trim(dir)//trim(filename)//'12')) open(unit=12,file=trim(dir)//trim(filename)//'12')
    if(FileExist(trim(dir)//trim(filename)//'13')) open(unit=13,file=trim(dir)//trim(filename)//'13')
    if(FileExist(trim(dir)//trim(filename)//'14')) open(unit=14,file=trim(dir)//trim(filename)//'14')
    if(FileExist(trim(dir)//trim(filename)//'15')) open(unit=15,file=trim(dir)//trim(filename)//'15')

    if(FileExist(trim(dir)//trim(filename)//'33')) then
      open(unit=33,file=trim(dir)//trim(filename)//'33')
      open(unit=35,file=trim(dir)//trim(filename)//'35')

      write(dir,"(A,'/','ThID_',I0,'_Shells_output/')") trim(ListDIR),ThID
      write(filename,"('qEarth_',I0,'_',I0,'.24')") ModNum,rpeat
      open(unit=34,file=trim(dir)//trim(filename))
    end if

end select

end subroutine


subroutine OpenOutput(ThID,ModNum,prog,ListDIR,plt,rpeat) ! open OrbData, Shells, OrbScore output files

integer,intent(in) :: ThID,ModNum
character(len=2),intent(in) :: prog ! SH=Shells, OD=Data, OS=Score
character(len=100),intent(in) :: ListDIR
logical,dimension(3),intent(in),optional :: plt
integer,intent(in),optional :: rpeat

character(len=200) :: filename
character(len=100) :: dir
logical,save :: callOD=.False.,callSH=.False.,callOS=.False.
character(len=100) :: ScoreFiles(9)
logical :: ScoreFilesL(7)


select case(prog)

  case("OD")

    write(dir,"(A,'/','ThID_',I0,'_Data_output/')") trim(ListDIR),ThID
    if(.not. callOD) then
      call execute_command_line('mkdir '//dir)
      callOD = .True.
    end if

    write(filename,"(A,'FEG_',i0,'.14')") trim(dir),ModNum
    open(unit=14,file=trim(filename),status='replace',action='write')

    write(filename,"(A,'Assign_Log_',i0,'.13')") trim(dir),ModNum
    open(unit=13,file=trim(filename),status='replace',action='write')


  case("SH")

    write(dir,"(A,'/','ThID_',I0,'_Shells_output/')") trim(ListDIR),ThID
    if(.not. callSH) then
      call execute_command_line('mkdir '//dir)
      callSH = .True.
    end if

    write(filename,"(A,'vEarth_',i0,'_',i0,'.22')") trim(dir),ModNum,rpeat
    open(unit=22,file=trim(filename),action='write',status='replace')

    write(filename,"(A,'fEarth_',i0,'_',i0,'.23')") trim(dir),ModNum,rpeat
    open(unit=23,file=trim(filename),action='write',status='replace')

    write(filename,"(A,'qEarth_',i0,'_',i0,'.24')") trim(dir),ModNum,rpeat
    open(unit=24,file=trim(filename),action='write',status='replace')

    write(filename,"(A,'Log_',i0'.21')") trim(dir),ModNum
    if(FileExist(trim(filename))) then
      open(unit=21,file=trim(filename),action='write',status='old',access='append')
    else
      open(unit=21,file=trim(filename),action='write',status='replace')
    end if


  case("OS")

    write(dir,"(A,'/','ThID_',I0,'_Score_output/')") trim(ListDIR),ThID
    if(.not. callOS) then
      call execute_command_line('mkdir '//dir)
      callOS = .True.
    end if

    if(plt(1)) then ! pltGeo
      write(filename,"('ERRORS_',I0,'_',I0,'.gps')") ModNum,rpeat
      open(unit=31,file=trim(dir)//trim(filename))
    end if

    if(plt(2)) then ! pltCos
      write(filename,"('COSEISMIC_',I0,'_',I0,'.gps')") ModNum,rpeat
      open(unit = 32,file=trim(dir)//trim(filename))
    end if

    if(plt(3)) then ! pltInt
      write(filename,"('INTERSEISMIC_',I0,'_',I0,'.gps')") ModNum,rpeat
      open(unit=37,file=trim(dir)//trim(filename))
    end if

    call readSCORE(ScoreFiles,ScoreFilesL)
    if(ScoreFilesL(5)) then ! Seismicity error(s)
      write(filename,"('fort_',I0,'_',I0,'.21')") ModNum,rpeat
      open(unit=21,file=trim(dir)//trim(filename))
      write(filename,"('fort_',I0,'_',I0,'.22')") ModNum,rpeat
      open(unit=22,file=trim(dir)//trim(filename))
    end if

end select

end subroutine



!------------------------------------------------------------------------------
! Grid
!------------------------------------------------------------------------------

subroutine CalcGrid(VarValsMin,VarValsMax,NumMods,Grid,Step) ! Calculate variable values in grid search level

real*8, dimension(:),intent(in) :: VarValsMin
real*8, dimension(:),intent(in) :: VarValsMax
integer,dimension(:),intent(in) :: NumMods
real*8,dimension(:,:),allocatable,intent(out) :: Grid ! 2D grid of vectors, max(NumMods)*NumVariables
real*8,dimension(:),intent(out) :: Step

integer :: i,j


allocate(Grid(maxval(NumMods),size(NumMods)))

do i = 1,size(NumMods)
  Step(i) = (VarValsMax(i) - VarValsMin(i) )/(NumMods(i))
  Grid(1,i) = VarValsMin(i) + (Step(i)/2.0)
  Grid(NumMods(i),i) = VarValsMax(i) - (Step(i)/2.0)
    do j = 2,NumMods(i)-1
      Grid(j,i) = (VarValsMin(i)+(Step(i)/2.0)) + (Step(i)*(j-1))
    end do
end do

end subroutine


subroutine DriverChk(MisType,Abort) ! Check that the driving misfit is actually run

character(len=3),intent(in) :: MisType
logical,intent(out) :: Abort

logical,dimension(6) :: MisRun


Abort = .True.
call MisfitRun(MisRun)

select case(MisType)

  case('GV')
    if(MisRun(1)) Abort = .False.

  case('SSR')
    if(MisRun(2)) Abort = .False.

  case('SD')
    if(MisRun(3)) Abort = .False.

  case('FSR')
    if(MisRun(4)) Abort = .False.

  case('SC')
    if(MisRun(5)) Abort = .False.

  case('SA')
    if(MisRun(6)) Abort = .False.

  case default
    print '(2A, /)', 'Unrecognised Misfit Type: ', MisType
    Abort = .True.

end select

end subroutine


subroutine GenerateModels(NumMods,TotNumMods,Grid,Models) ! Calculate models in grid search level

integer,dimension(:),intent(in)   :: NumMods
integer,intent(in) :: TotNumMods
real*8, dimension(:,:),intent(in) :: Grid
real*8,dimension(:,:),intent(out) :: Models

integer :: i,j,k,ModsBefore,ModsAfter,jump


do k = size(Grid,2),1,-1 ! k controls accessed variable (starting from final)
  ModsBefore = 1

  do i = 1,size(NumMods) - (1 + size(Grid,2) - k)
    ModsBefore = ModsBefore*NumMods(i)
  end do

  ModsAfter = (TotNumMods/ModsBefore)/NumMods(i) ! ModsAfter is repeat number

  do i = 1,ModsBefore
    jump = (i-1)*NumMods(k)*ModsAfter
    do j = 1,NumMods(k)
      Models(((j-1)*ModsAfter)+1+jump:(j*ModsAfter)+jump,k) = Grid(j,k) ! create all models in Models
    end do
  end do
end do

end subroutine


subroutine MisfitAnalysis(misfits,Cells,MisType,Keep) ! Find best models

real*8,dimension(:),intent(in) :: misfits
integer,intent(in) :: Cells
character(len=3),intent(in) :: MisType
integer,dimension(:),allocatable,intent(out) :: Keep

logical,dimension(:),allocatable :: mask
integer :: i


allocate(Keep(Cells),mask(size(misfits,1)))

mask = .True.

if(MisType == 'SC') then
  do i = 1,Cells
    Keep(i) = maxloc(misfits,dim=1,mask=mask)
    mask(maxloc(misfits,mask=mask)) = .False.
  end do
else
  do i = 1,Cells
    Keep(i) = minloc(misfits,dim=1,mask=mask)
    mask(minloc(misfits,mask)) = .False.
  end do
end if

end subroutine


subroutine MisfitRun(MisRun) ! Find which misfits are run (used in ShellSetMain & DriverChk)

logical,dimension(6),intent(out) :: MisRun

character(len=100),dimension(9):: ScoreFiles
logical,dimension(7) :: ScoreFilesL


call readSCORE(ScoreFiles,ScoreFilesL,MisRun=MisRun)

end subroutine


subroutine UpdateSetup(Models,Step,VarValsMin,VarValsMax) ! Update setup variables to new grid region

real*8,dimension(:),intent(in) :: Models
real*8,dimension(:),intent(in) :: Step
real*8,dimension(:),intent(out) :: VarValsMin
real*8,dimension(:),intent(out) :: VarValsMax

integer :: i


do i = 1,size(Models,1)
  VarValsMin(i) = Models(i) - Step(i)/2.0
  VarValsMax(i) = Models(i) + Step(i)/2.0
end do

end subroutine



!------------------------------------------------------------------------------
! CLA
!------------------------------------------------------------------------------

subroutine CLA_check(CLA_STOP) ! Check number of and for non-running CLAs

logical,intent(out) :: CLA_STOP

integer :: CLA_Count,i
character(len=10) :: arg


CLA_STOP = .False.
CLA_Count = command_argument_count()

if(CLA_Count == 0) then
  CLA_STOP = .True.
else
  do i = 1,CLA_Count
    call get_command_argument(i,arg)
    if(trim(arg) == '-help') then
      call help()
      CLA_STOP = .True.
    elseif(trim(arg) == '-info') then
      call info() ;
      CLA_STOP = .True.
    elseif(trim(arg) == '-cite') then
      call cite() ;
      CLA_STOP = .True.
    end if
  end do
end if

end subroutine


subroutine ReadCLA(MisConv,MisType,MisRange,MinIter,MaxIter,DirName,Verbose,AllErrorFatal,InOpt,TargMisM,TargMisK,Abort) ! Read CLAs

integer,intent(out) :: MinIter,MaxIter
logical,intent(out) :: MisConv,Verbose,AllErrorFatal,Abort
real,intent(out) :: MisRange
character(len=3),intent(out) :: MisType
character(len=100),intent(out) :: DirName
character(len=4),intent(out) :: InOpt
real,intent(out) :: TargMisM,TargMisK

integer :: i,CLA_Count,ierr,c1,c2
logical :: exists
character(len=100) :: arg


TargMisM = 0.0
TargMisK = 0.0
MinIter = 0
MaxIter = 0
MisConv = .false.
Verbose = .false.
AllErrorFatal = .false.
InOpt = '    '
MisRange = 0.0
MisType = 'null'
DirName = 'test'
Abort = .False.
arg = ' '
i=1

do while(i<=command_argument_count())
  call get_command_argument(i,arg)

  if(trim(arg) == '-ML') then
    call get_command_argument(i+1,arg)
    read(arg,'(F8.4)',iostat=ierr) TargMisM
    if(ierr /= 0) then
      print"(A)","Error with TargMisM - program will terminate"
      Abort = .True.
    end if
    i = i+2

  elseif(trim(arg) == '-MC') then
    MisConv = .True.
    call get_command_argument(i+1,arg)
    c1 = index(arg, ',')
    c2 = index(arg(c1+1:), ',')
    if(trim(arg(:c1-1)) == 'GV') then
      MisType = 'GV'
    elseif(trim(arg(:c1-1)) == 'SSR') then
      MisType = 'SSR'
    elseif(trim(arg(:c1-1)) == 'SD') then
      MisType = 'SD'
    elseif(trim(arg(:c1-1)) == 'FSR') then
      MisType = 'FSR'
    elseif(trim(arg(:c1-1)) == 'SC') then
      MisType = 'SC'
    elseif(trim(arg(:c1-1)) == 'SA') then
      MisType = 'SA'
    else
      print"(A)","No valid option for Misfit type was entered - program will terminate"
      Abort = .True.
    end if
    read(arg(c1+1:c1+c2-1),'(I4)',iostat=ierr) MinIter
    read(arg(c1+c2+1:),'(F8.4)',iostat=ierr) MisRange
    i = i+2

  elseif(trim(arg) == '-V') then
    Verbose = .True.
    i = i+1

  elseif(trim(arg) == '-AEF') then
    AllErrorFatal = .True.
    i = i+1

  elseif(trim(arg) == '-MT') then
    call get_command_argument(i+1,arg)
    if(trim(arg) == 'GV') then
      MisType = 'GV'
    elseif(trim(arg) == 'SSR') then
      MisType = 'SSR'
    elseif(trim(arg) == 'SD') then
      MisType = 'SD'
    elseif(trim(arg) == 'FSR') then
      MisType = 'FSR'
    elseif(trim(arg) == 'SC') then
      MisType = 'SC'
    elseif(trim(arg) == 'SA') then
      MisType = 'SA'
    else
      print"(A)","No valid option for Misfit type was entered - program will terminate"
      Abort = .True.
    end if
    i = i+2

  elseif(trim(arg) == '-KL') then
    call get_command_argument(i+1,arg)
    read(arg,'(F8.4)',iostat=ierr) TargMisK
    if(ierr /= 0) then
      print"(A)","Error with TargMisK - program will terminate"
      Abort = .True.
    end if
    i = i+2

  elseif(trim(arg) == '-Iter') then
    call get_command_argument(i+1,arg)
    read(arg,'(I4)',iostat=ierr) MaxIter
    if(ierr /= 0) then
      print"(A)","Error with Iter - program will terminate"
      Abort = .True.
    end if
    i = i+2

  elseif(trim(arg) == '-Dir') then
    call get_command_argument(i+1,arg)
    read(arg,'(A)',iostat=ierr) DirName
    if(ierr /= 0) then
      print"(A)","Error with DirName - program will terminate"
      Abort = .True.
    else
      call execute_command_line('rm -rf '//trim(DirName))
    end if
    i = i+2

  elseif(trim(arg) == '-InOpt') then
    call get_command_argument(i+1,arg)
    read(arg,'(A4)',iostat=ierr) InOpt
    if(ierr /= 0) then
      print"(A)","Error with InOpt - program will terminate"
      Abort = .True.
    end if
    i = i+2

  else
    print '(2A, /)', "Unrecognised command-line option: ", arg
    Abort = .True.

  end if

end do

! Check everything required is known
if(MisConv) then ! MisConv requires: MisType, MisRange, MinIter
  if((MisType /= 'GV'  .and. MisType /= 'SSR' .and. MisType /= 'SD' .and. MisType /= 'FSR' .and. &
  &  MisType /= 'SC' .and. MisType /= 'SA')  .or. MisRange == 0.0 .or. MinIter == 0.0) then
    print*,"Not everything required for Misfit convergence was entered - program will terminate"
    print*,"MisType=  ",trim(MisType)
    print*,"MisRange= ",MisRange
    print*,"MinIter=  ",MinIter
    Abort = .True.
  end if

end if

! Remove existing Models_Lim.txt file
if(TargMisM /= 0.0 .and. FileExist(trim(DirName)//"Models_Lim.txt")) call execute_command_line('rm '//trim(DirName)//'Models_Lim.txt')

end subroutine


! Non-running CLA subroutines
! need to check and finish
subroutine help() ! Print help information

print"(/A)","The program requires certain input in the form of command line arguments (CLAs)"
print"(A)","Here are two example invokation lines:"
print"(A/,A/,A)","mpirun -np 2 ./ShellSet -Iter 5 -Dir grid -InOpt Grid -MT SD",&
   & "    This runs a Grid search with 5 iterations of the main work loop (see figure .... of guide)",&
   & "    using the Stress Direction misfit to select best models, the working directory is 'grid'"
print"(A/,A/,A/,A/)","mpirun -np 2 ./ShellSet -Iter 7 -Dir list -InOpt List -MC GV,2,0.2",&
   & "    This runs List input models with 7 iterations of the main work loop (see figure .... of guide)",&
   & "    using the Geodetic Velocity misfit to perform misfit convergance tests after the 2nd work loop",&
   & "    iteration, the working directory is 'list'"

print"(A)","An explanaition of the CLAs follows:"

print"(A/,A/,A/,A/)","The following CLAs print information to the terminal without running the program:", &
   & "-help  -> print information currently shown", &
   & "-info  -> print general information about the program", &
   & "-cite  -> print information about citing this program"

print"(A/,A/,A/)","The following CLAs are required for the program to function:", &
   & "-Iter  -> number of iterations of 'main loop'", &
   & "-InOpt -> either 'List' or 'Grid' to determine model input generation"

print"(A/,A/,A/,A/,A/,A/,A/,A/)","The following CLAs are optional and add functionality:", &
   & "-MC  -> information for misfit convergence", &
   & "-MT  -> misfit type used to drive grid search", &
   & "-Dir -> new working directory", &
   & "-AEF -> make all errors fatal", &
   & "-V   -> create Verbose file containing extra informational output", &
   & "-ML  -> misfit limit beneath which all models are stored in an extra output file", &
   & "-KL  -> misfit limit at which the program is ended"

call abort(100)

call moreinfo()

end subroutine


subroutine cite() ! Print citation information

print"(/A/,A/,A/)","The most up-to date version of ShellSet can be found on GitHub, with a complete", &
   &   "program version history available on Zenodo. If you use any version of ShellSet in your work", &
   &   "please cite it using the program DOI"

print"(A/,A/)","A publication related to ShellSet is under review, upon acceptance and publication", &
   &   "information will be added allowing the citation of the paper."

call moreinfo()

end subroutine


subroutine info() ! Print program information

print"(/A/,A/,A/)","This program is a parallel combination of OrbData5, SHELLSV5.0 and OrbScore2 which were", &
   &    "developed by Prof. Peter Bird UCLA. The original programs, input files, detailed explanaitions and", &
   &    "examples are available from Prof. Birds website: http://peterbird.name"

print"(A/,A/)","The program allows for the testing of multiple input variables inside the SHELLS simulator", &
   &    "in parallel using either a list of variable name & values or a grid search option."

print"(A/,A/)","The program was developed inside a WSL2 environment running Ubuntu, using Intel OneAPI", &
   &    "Basekit & HPCkit. It has been tested both in the development environment and Linux OS machines"

call moreinfo()

end subroutine


subroutine moreinfo() ! Print extra information

print"(A/,A/)","For more information please check the README file and guide that should have been", &
   &    "included within the download package."

print "(A/,A/,A/,A/,A/,A/,A/)", "**********************************************************************************", &
   &    "ShellSet Copyright (C) 2002 Jon May, Peter Bird, Michele Carafa", &
   &    "This program comes with ABSOLUTELY NO WARRANTY.", &
   &    "This is free software, you are welcome to redistribute it under conditions", &
   &    "set out in GNU General Public License version 3, or any later version.", &
   &    "See the included license file, or https://www.gnu.org/licenses/, for details.", &
   &    "**********************************************************************************"

print"(A/)","If you use any version of ShellSet in your work please cite it using the DOI: "

end subroutine



!------------------------------------------------------------------------------
! From Original Program Units
!------------------------------------------------------------------------------

SUBROUTINE ReadPm (iUnit7, iUnitT, names , numPlt, offMax, & ! Read parameter input file
&                    aCreep, alphaT, bCreep, Biot  , &         ! output
&                    Byerly, cCreep, cFric , conduc, &
&                    dCreep, eCreep, everyP, fFric , gMean , &
&                    gradie, iConve, iPVRef, &
&                    maxItr, okDelV, okToQt, oneKm,  radio,  &
&                    radius, refStr, rhoAst, rhoBar, rhoH2O, &
&                    tAdiab, tauMax, temLim, title3, &
&                    trHMax, tSurf,  vTimes, zBAsth, pltRef)



!   Reads strategic and tactical input parameters from device iUnit7,
!   and echoes them on device iUnitT with annotations.

!      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       INTEGER, INTENT(IN) :: iUnit7, iUnitT                                                  ! input
       CHARACTER*2, INTENT(IN) :: names                                                       ! input
       INTEGER, INTENT(IN) :: numPlt                                                          ! input
       REAL*8, INTENT(IN) :: offMax                                                           ! input
       REAL*8, INTENT(OUT) :: aCreep, alphaT, bCreep, Biot, Byerly, cCreep, cFric , conduc, & ! output
          & dCreep, eCreep                                                                    ! output
       LOGICAL, INTENT(OUT) :: everyP                                                         ! output
       REAL*8, INTENT(OUT) :: fFric , gMean , gradie                                          ! output
       INTEGER, INTENT(OUT) :: iConve, iPVRef, maxItr                                         ! output
       REAL*8, INTENT(OUT) :: okDelV, okToQt, oneKm,  radio, radius, refStr, &                ! output
          & rhoAst, rhoBar, rhoH2O, tAdiab, tauMax, temLim                                    ! output
       CHARACTER*100, INTENT(OUT) :: title3                                                    ! output
       REAL*8, INTENT(OUT) :: trHMax, tSurf,  vTimes, zBAsth                                  ! output
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       CHARACTER*2,intent(out) :: pltRef
       INTEGER i, ios
       REAL*8 tempV, vector
       DIMENSION aCreep(2), alphaT(2), bCreep(2), cCreep(2), conduc(2), &
     &           dCreep(2), names(numplt), radio(2), &
     &           rhoBar(2), tauMax(2), temLim(2), tempv(2), vector(2)

       IF(Verbose) WRITE(iUnitT,1) iunit7
    1  FORMAT(//' Attempting to read input parameter file from unit ', I3/)
       title3 = ' '
       READ (iunit7, 2, IOSTAT = ios) title3
       IF (ios /= 0) THEN
          write(ErrorMsg,'(A)') "File not found, or file is empty, or file is too short."
          call FatalError(ErrorMsg,ThID)
       END IF
    2  FORMAT (A80)
       IF(Verbose) WRITE (iUnitT, 3) title3
    3  FORMAT (/' [OMIT from iXXX.in] Title of parameter set ='/' ',A80)

       IF(Verbose) WRITE (iUnitT, 4)
    4  FORMAT (' [OMIT from iXXX.in]' &
     &        /' [OMIT from iXXX.in]', &
     &         ' **************************************************' &
     &        /' [OMIT from iXXX.in]', &
     &         ' It is the user''s responsibility to input all of the' &
     &        /' [OMIT from iXXX.in]', &
     &         ' following numerical quantities in consistent units,' &
     &        /' [OMIT from iXXX.in]', &
     &         ' such as Systeme International (SI) or cm-g-s (cgs).' &
     &        /' [OMIT from iXXX.in]', &
     &         ' Note that time unit must be the second (hard-coded).' &
     &        /' [OMIT from iXXX.in]', &
     &         ' **************************************************' &
     &        /' [OMIT from iXXX.in]' &
     &        /' [OMIT from iXXX.in]', &
     &         ' ========== Strategic parameters (define the real', &
     &                                    '-Earth problem) ======' &
     &        /' [OMIT from iXXX.in]')

       READ (iunit7, * , IOSTAT = ios) fFric
       IF (ios /= 0) THEN
          write(ErrorMsg,'(A)') "File not found, or file is empty, or file is too short."
          call FatalError(ErrorMsg,ThID)
       END IF
       IF(Verbose) WRITE (iUnitT, 20) fFric
   20  FORMAT (' ',11X,F10.3,' fFric  = coefficient of friction on faults')
       IF (fFric < 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: negative fault friction fFric is not physical."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) cFric
       IF(Verbose) WRITE (iUnitT, 30) cFric
   30  FORMAT (' ',11X,F10.3,' cFric  = coefficient of friction within blocks')
       IF (cFric <= 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: continuum friction cFric must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) Biot
       IF(Verbose) WRITE (iUnitT, 40) Biot
   40  FORMAT (' ',11X,F10.4,' Biot   = effective-pressure coefficient', &
     &          ' of Biot: 0. (dry) to 1. (wet)')
       IF ((Biot < 0.0D0).OR.(Biot > 1.0D0)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: Biot coefficient must be in range 0.0 to 1.0."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) Byerly
       IF (offMax > 0.0D0) THEN
            IF(Verbose) WRITE (iUnitT, 43) Byerly
   43       FORMAT (' ',11X,F10.4,' Byerlee coefficient (0. to 0.999) ='/ &
     &              11X,' fractional friction reduction on master fault'/ &
     &              11X,' (Other faults have less reduction, in proportion to'/ &
     &              11X,' their total past offsets)')
            IF ((Byerly < 0.0D0).OR.(Byerly > 1.0D0)) THEN
              write(ErrorMsg,'(A)') "ERROR in parameter input file: Byerlee coefficient must be in range 0.0 to 1.0."
              call FatalError(ErrorMsg,ThID)
            END IF
       ELSE
            IF(Verbose) WRITE (iUnitT, 46) Byerly
   46       FORMAT (' ',11X,F10.4,' Byerlee coefficient (not used in', &
     &                  ' this run, as all fault offsets are zero).')
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             aCreep)              ! output
       IF(Verbose) WRITE (iUnitT, 50) aCreep(1), aCreep(2)
   50  FORMAT (' ',1P, E10.2,' ',E10.2,' aCreep = A for creep = ', &
     &            'pre-exponential shear', &
     &        ' stress constant for creep. (crust/mantle)')
       IF ((aCreep(1) <= 0.).OR.(aCreep(2) <= 0.0D0)) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: aCreep must be positive in each layer."
         call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             bCreep)              ! output
       IF(Verbose) WRITE (iUnitT, 60) bCreep(1), bCreep(2)
   60  FORMAT (' ',F10.0,' ',F10.0,' bCreep = B for creep =', &
     &             ' activation_energy/R/n', &
     &           ' (in K). (crust/mantle)')
       IF ((bCreep(1) < 0.0D0).OR.(bCreep(2) < 0.0D0)) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: Negative bCreep in either layer is unphysical."
         call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             cCreep)              ! output
       IF(Verbose) WRITE (iUnitT, 70) cCreep(1), cCreep(2)
   70  FORMAT (' ',1P, E10.2,' ',E10.2,' cCreep = C for creep = rho*', &
     &            'g*V_star*eCreep/R (in K/m). (crust/mantle)')
       IF ((cCreep(1) < 0.0D0).OR.(cCreep(2) < 0.0D0)) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: Negative cCreep in either layer is unphysical."
         call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             dCreep)              ! output
       IF(Verbose) WRITE (iUnitT, 80) dCreep(1), dCreep(2)
   80  FORMAT (' ',1P,E10.2,' ',E10.2,' dCreep = D for creep = max', &
     &         'imum shear stress under any conditions. (crust/mantle)')
       IF ((dCreep(1) <= 0.0D0).OR.(dCreep(2) <= 0.0D0)) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: dCreep must be positive in each layer."
         call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) eCreep
       IF(Verbose) WRITE (iUnitT, 90) eCreep
   90  FORMAT (' ',11X,F10.6,' eCreep = E for creep = strain-rate expo', &
     &           'nent for creep (1/n).  (Same for crust and mantle!)')
       IF (eCreep <= 0.0D0) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: eCreep must be positive."
         call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) tAdiab, gradie
       IF(Verbose) WRITE (iUnitT, 92) tAdiab, gradie
   92  FORMAT (' ',F10.0,' ',1P,E10.2,' tAdiab, GRADIE = intercept and ' &
     &        ,'slope of upper mantle adiabat below plate (K, K/m)')
       IF ((tAdiab < 0.0D0).OR.(gradie < 0.0D0)) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: Negative Kelvin temperature and/or negative adiabatic gradient is/are unphysical."
         call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) zBAsth
       IF(Verbose) WRITE (iUnitT, 94) zBAsth
   94  FORMAT (' ',11X,1P,E10.2,' zBAsth = depth of base of', &
     &                          ' asthenosphere')
       IF (zBAsth <= 0.0D0) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: zBAsth must be positive."
         call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, 952) pltRef
  952  FORMAT(A2)
       IF(Verbose) WRITE (iUnitT, 954) pltRef
  954  FORMAT(' ',A2,'<==================', &
     &        ' pltRef = plate defining velocity ', &
     &        'reference frame (AF, NA, EU, ...)')
       iPVRef = 0
       DO 956 i = 1, numPlt
            IF (names(i) == pltRef) iPVRef = i
  956  CONTINUE
       IF (iPVRef == 0) THEN
          write(ErrorMsg,'(A/,A/,A)') "ERROR in parameter input file: ",&
          &  "In line 13 (after zBAsth, before iConve), in the first two columns of the line,",&
          &  "define the velocity reference frame by entering one of the following plate names:"
          allocate(ErrorArrayChar(numPlt))
          ErrorArrayChar = names
          call FatalError(ErrorMsg,ThID,ErrArrCh=ErrorArrayChar)
       END IF

       READ (iunit7, * ) iConve
       IF(Verbose) WRITE (iUnitT, 96) iConve
   96  FORMAT (' ',11X,I10,' iConve = code for mantle flow below the ', &
     &                           'asthenosphere:' &
     &/' ','[OMIT from iXXX.in]  ',' 0 = static (with respect to AF)' &
     &/' ','[OMIT from iXXX.in]  ',' 1 = Hager and O''Connell (1979)', &
     &                                 ' Model II' &
     &/' ','[OMIT from iXXX.in]  ',' 2 = Baumgardner (1988) Figure', &
     &                                                    ' 7A-F, *10.' &
     &/' ','[OMIT from iXXX.in]  ',' 3 = PB2002 (Bird, 2003)' &
     &/' ','[OMIT from iXXX.in]  ',' 4 = PB2002 drags continents;', &
     &                                 ' no ocean drag' &
     &/' ','[OMIT from iXXX.in]  ',' 5 = drag on base of subduction', &
     &                                 ' forearc only' &
     &/' ','[OMIT from iXXX.in]  ',' 6 = sense & traction from trac', &
     &                                  'tion pole vectors' &
     & )
       IF ((iConve < 0).OR.(iConve > 6)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: iConve is not one of the values listed."
          call FatalError(ErrorMsg,ThID)
       END IF
       IF (iConve > 0) THEN
            BACKSPACE iunit7
            CALL ReadN (iunit7, iUnitT, 2, & ! input
     &                  tempv)               ! output
            IF (tempv(2) >= 0) THEN
                 vTimes = tempv(2)
                 IF(Verbose) WRITE (iUnitT, 98) vTimes
   98            FORMAT (' ',11X,F10.4,' vTimes = speed factor for con', &
     &                   'vection model identified above')
            ELSE
              write(ErrorMsg,'(A)') "ERROR in parameter input file: Uninterpretable value for vTimes."
              call FatalError(ErrorMsg,ThID)
            END IF
       ELSE
            vTimes = 1.0D0
       END IF

       READ (iunit7, * ) trHMax
       IF(Verbose) WRITE (iUnitT, 101) trHMax
  101  FORMAT (' ',11X,1P,E10.2,' trHMax = limit on horizontal', &
     &                ' tractions applied to base of plate')
       IF (trHMax < 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: trHMax may not be negative."
          call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             vector)              ! output
       tauMax(1) = vector(1)
       tauMax(2) = vector(2)
!      Provide for old parameter files with only one tauMax:
       IF (tauMax(2) <= 0.0D0) tauMax(2) = tauMax(1)
       IF(Verbose) WRITE (iUnitT, 106) tauMax(1), tauMax(2)
  106  FORMAT (' ',1P, E10.2,' ',E10.2, &
     &         ' tauMax = upper limit(s) on subduction zone shear', &
     &         ' coupling, integrated down-dip (N/m).  One value:', &
     &         ' universal. Two values: sea, land.')
       IF ((tauMax(1) < 0.0D0).OR.(tauMax(2) < 0.0D0)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: tauMax may not be negative."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) rhoH2O
       IF(Verbose) WRITE (iUnitT, 110) rhoH2O
  110  FORMAT (' ',11X,1P,E10.3,' rhoH2O = density of groundwater,', &
     &                          ' lakes, & oceans')
       IF (rhoH2O < 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: rhoH2O may not be negative."
          call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             rhoBar)              ! output
       IF(Verbose) WRITE (iUnitT, 120) rhoBar(1), rhoBar(2)
  120  FORMAT (' ',1P,E10.3,' ',E10.3,' rhoBar = mean density,', &
     &        ' corrected to 0 degrees Kelvin. (crust/mantle)')
       IF ((rhoBar(1) <= 0.0D0).OR.(rhoBar(2) <= 0.0D0)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: rhoBar must be positive in each layer."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) rhoAst
       IF(Verbose) WRITE (iUnitT, 130) rhoAst
  130  FORMAT (' ',11X,1P,E10.3,' rhoAst = density of asthenosphere')
       IF (rhoAst <= 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: rhoAst must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) gMean
       IF(Verbose) WRITE (iUnitT, 140) gMean
  140  FORMAT (' ',11X,1P,E10.3,' gMean  = mean gravitational', &
     &                          ' acceleration', &
     &       ' (length/s**2)')
       IF (gMean <= 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: gMean must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) oneKm
       IF(Verbose) WRITE (iUnitT, 150) oneKm
  150  FORMAT (' ',11X,1P,E10.3,' oneKm  = number of length units', &
     &  ' needed to make 1 kilometer (e.g., 1000. in SI, 1.D5 in cgs)')
       IF (oneKm <= 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: oneKm must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) radius
       IF(Verbose) WRITE (iUnitT, 155) radius
  155  FORMAT (' ',11X,1P,E10.3,' radius = radius of the planet')
       IF (radius <= 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: radius must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             alphaT)              ! output
       IF(Verbose) WRITE (iUnitT, 160) alphaT(1), alphaT(2)
  160  FORMAT (' ',1P,E10.2,' ',E10.2,' alphaT = volumetric thermal', &
     &                                ' expansion', &
     &            ' (1/V)*(dV/dT). (crust/mantle)')
       IF ((alphaT(1) < 0.0D0).OR.(alphaT(2) < 0.0D0)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: Negative alphaT in either layer is unphysical."
          call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             conduc)              ! output
       IF(Verbose) WRITE (iUnitT, 170) conduc(1), conduc(2)
  170  FORMAT (' ',1P,E10.2,' ',E10.2,' conduc = thermal conductivity,', &
     &        ' energy/length/s/deg. (crust/mantle)')
       IF ((conduc(1) <= 0.0D0).OR.(conduc(2) <= 0.0D0)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: conduc must be positive in each layer."
          call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             radio)               ! output
       IF(Verbose) WRITE (iUnitT, 180) radio(1), radio(2)
  180  FORMAT (' ',1P,E10.2,' ',E10.2,' radio  = radioactive heat ', &
     &        'production, energy/volume/s. (crust/mantle)')
       IF ((radio(1) < 0.0D0).OR.(radio(2) < 0.0D0)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: Negative radio in either layer is unphysical."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) tSurf
       IF(Verbose) WRITE (iUnitT, 185) tSurf
  185  FORMAT (' ',11X,F10.0,' tSurf  = surface temperature, on', &
     &        ' absolute scale (deg. K)')
       IF (tSurf <= 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: tSurf must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       CALL ReadN (iunit7, iUnitT, 2, & ! input
     &             temLim)              ! output
       IF(Verbose) WRITE (iUnitT, 190) temLim(1), temLim(2)
  190  FORMAT (' ',F10.0,' ',F10.0,' temLim = convecting', &
     &       ' temperature (Tmax), on absolute scale. (crust/mantle)')
       IF ((temLim(1) <= 0.0D0).OR.(temLim(2) <= 0.0D0)) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: temLim must be positive in each layer."
          call FatalError(ErrorMsg,ThID)
       END IF

       IF(Verbose) WRITE (iUnitT, 199)
  199  FORMAT (' [OMIT from iXXX.in]' &
     &        /' [OMIT from iXXX.in]', &
     &         ' ========== Tactical parameters (How to reach ', &
     &                             'the solution?) ==========' &
     &        /' [OMIT from iXXX.in]')

       READ (iunit7, * ) maxItr
       IF(Verbose) WRITE (iUnitT, 200) maxItr
  200  FORMAT (' ',11X,I10,' maxItr = maximum iterations within', &
     &           ' velocity solution')
       IF (maxItr < 1) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: maxItr must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) okToQt
       IF(Verbose) WRITE (iUnitT, 210) okToQt
  210  FORMAT (' ',11X,F10.6,' okToQt = acceptable fractional change', &
     &                  ' in velocity (stops iteration early)')
       IF (okToQt < 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: okToQt may not be negative."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) refStr
       IF(Verbose) WRITE (iUnitT, 220) refStr
  220  FORMAT (' ',11X,1P,E10.2,' refStr = expected mean value of', &
     &        ' shear stress in plate (used to set stiffness limits)')
       IF (refStr <= 0.0D0) THEN
          write(ErrorMsg,'(A)') "ERROR in parameter input file: refStr must be positive."
          call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) okDelV
       IF(Verbose) WRITE (iUnitT, 230) okDelV
  230  FORMAT (' ',11X,1P,E10.2,' okDelV = magnitude of velocity', &
     &                  ' errors allowed due to finite stiffness' &
     &        /' [OMIT from iXXX.in]  ', &
     & ' (Such errors may appear in such forms as:' &
     &        /' [OMIT from iXXX.in]  ', &
     & '  1. fictitious basal slip of plate over asthenosphere' &
     &        /' [OMIT from iXXX.in]  ', &
     & '  2. erroneous convergence/divergence at vertical faults' &
     &        /' [OMIT from iXXX.in]  ', &
     & '  3. velocity effect of fictitious viscous compliances' &
     &        /' [OMIT from iXXX.in]  ', &
     & ' HOWEVER, VALUES WHICH ARE TOO SMALL WILL CAUSE ILL-CONDITIONED' &
     &        /' [OMIT from iXXX.in]  ', &
     & ' LINEAR SYSTEMS AND STRESS ERR0RS, ', &
     &                                'AND MAY PREVENT CONVERGENCE!)' &
     &)
       IF (okDelV <= 0.0D0) THEN
         write(ErrorMsg,'(A)') "ERROR in parameter input file: okDelV must be positive."
         call FatalError(ErrorMsg,ThID)
       END IF

       READ (iunit7, * ) everyP
       IF(Verbose) WRITE (iUnitT, 240) everyP
  240  FORMAT (' ',11X,L10,' everyP = should nodal velocities be', &
     &        ' output in every iteration? (for convergence studies)')

       IF(Verbose) WRITE (iUnitT, 999)
  999  FORMAT (' --------------------------------------------------', &
     &          '-----------------------------')
END SUBROUTINE ReadPm


SUBROUTINE ReadN (iUnitP, iUnitT, n, & ! Split input lines (used by ReadPm)
&                   vector)              ! INTENT(OUT)

!   A utility routine designed to permit -Faults- input files
!      to also be used by -Plates-, which expects more numbers
!      in some records.
!   This routine attempts to READ 'n' floating-point values
!      (using * FORMAT) from the next record on device 'iUnitP'.
!   If the input line is too short, the missing values are set to zero.
!   Note that this routine can also be used for reading INTEGER input,
!      as long as the calling program takes responsibility for
!      converting REAL*8 --> INTEGER after each CALL.

       INTEGER, INTENT(IN) :: iUnitP, iUnitT, n
       REAL*8, DIMENSION(:), INTENT(OUT) :: vector

       CHARACTER*1   :: c
       CHARACTER*132 :: line
       INTEGER :: i, number
       LOGICAL       :: anyIn, dotted, expon, signed

       READ (iUnitP, "(A)") line

       number = 0
       anyIn = .FALSE.
       expon = .FALSE.
       signed = .FALSE.
       dotted = .FALSE.
       DO 10 i = 1, 132
            c = line(i:i)
            IF ((c == ' ').OR.(c == ',').OR.(c == '/')) THEN
                 signed = .FALSE.
                 expon = .FALSE.
                 dotted = .FALSE.
                 IF (anyIn) THEN
                      number = number + 1
                      anyIn = .FALSE.
                 END IF
            ELSE IF ((c == '+').OR.(c == '-')) THEN
                 IF (signed) THEN
                      GO TO 50
                 ELSE
                      signed = .TRUE.
                 END IF
            ELSE IF ((c == 'E').OR.(c == 'D').OR. &
     &               (c == 'e').OR.(c == 'd')) THEN
                 IF (expon) THEN
                      GO TO 50
                 ELSE
                      expon = .TRUE.
                      signed = .FALSE.
                      dotted = .TRUE.
                 END IF
            ELSE IF (c == '.') THEN
                 IF (dotted) THEN
                      GO TO 50
                 ELSE
                      dotted = .TRUE.
                 END IF
            ELSE IF ((c == '0').OR.(c == '1').OR.(c == '2').OR. &
     &               (c == '3').OR.(c == '4').OR.(c == '5').OR. &
     &               (c == '6').OR.(c == '7').OR.(c == '8').OR. &
     &               (c == '9')) THEN
                 signed = .TRUE.
                 anyIn = .TRUE.
            ELSE
                 GO TO 50
            END IF
   10  CONTINUE
       IF (anyIn) number = number + 1

   50  IF (number == 0) THEN
         write(ErrorMsg,'(A,I0,A/,A)') "A LINE OF PARAMETER INPUT WHICH WAS SUPPOSED TO CONTAIN 1~",n," NUMBERS",&
            & "COULD NOT BE INTERPRETED. LINE FOLLOWS: ",line
         call FatalError(ErrorMsg,ThID)
       ELSE
            number = MIN(number, n)
            BACKSPACE iUnitP
            READ (iUnitP, * ) (vector(i), i = 1, number)
            IF (number < n) THEN
                 DO 99 i = number + 1, n
                      vector(i) = 0.0D0
   99            CONTINUE
            END IF
       END IF
END SUBROUTINE ReadN


end module
