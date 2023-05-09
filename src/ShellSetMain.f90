!*******************************************************************************
! ShellSet program
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


program ShellSet

use mkl_service
use ShellSetSubs
use SHELLS
use OrbScore
use OrbData
use SharedVars
use VariableCheck


implicit none

include "mpif.h"


! List input variables
character(len=10),dimension(:),allocatable :: ListVarNames
real*8,dimension(:,:),allocatable :: ListVarValues
integer,dimension(:),allocatable :: lengths
integer :: NumVar

! Grid input variables
integer,dimension(:),  allocatable :: NumMods,Keep
integer,dimension(:,:),allocatable :: WorkList
real*8, dimension(:),  allocatable :: VarValsMin,VarValsMax,Step,StepOld
real*8, dimension(:,:),allocatable :: Grid,AllModels,KeepCells
integer :: CellNumMods,Cells,Levels,RptMod,iLevel,iCell,N,Tag,NumMis,ModWrt=0,TotNumMods
logical :: ModRpt,Kill
real :: TargMisK,TargMisM
logical,dimension(:),allocatable :: mask

! MKL setup variables
integer :: cpus,MKLthreads

! MPI specific variables
integer :: np,ierr,msg_stat(MPI_STATUS_size)

! Command line argument Variables
integer :: MinIter,MaxIter
logical :: MisConv,AEF,MPIAbort,CLA_STOP=.False.,ListIn ! True if list input, False if grid search
real :: MisRange
character(len=3) :: MisRank
character(len=100) :: VerbFile,DirName
character(len=4) :: InOpt
character(len=3),dimension(:),allocatable :: GM_Types
integer :: NumGM

! Work iteration & control
integer :: ModNum,rpeat,recvd,repeated
real*8 :: misfitMin,misfitMax
real*8,dimension(:),allocatable :: Misfits
real*8,dimension(:),allocatable :: allmisfits
logical :: OData,Run,Work,SHELLSconv,FErrChk=.False.

! Updating variables inside iteration
character(len=100) :: Iterfilename
character(len=100),dimension(9) :: ShellsFiles
logical,dimension(4) :: ShellsFilesL

character(len=500) :: frmt,filename
character(len=200) :: arg,outmsg,cmd
logical,dimension(7) :: MisRun
character(len=100),dimension(6) :: DataFiles
real*8 :: start,finish,ProgStart,ProgFinish

!------------------------------------------ Executing code ----------------------------------------------
! ----- Start MPI -----
call MPI_INIT (ierr)
call MPI_COMM_size (MPI_COMM_WORLD, np, ierr)
call MPI_COMM_RANK (MPI_COMM_WORLD, ThID, ierr)


if(ThID == 0) then
  call cpu_time(ProgStart)
  if(np == 1) then
    call abort(0)
    call MPI_Abort(MPI_COMM_WORLD,0)
  end if

! ----- CLA -----
  call CLA_check(CLA_STOP)
  if(CLA_STOP) then
    if(command_argument_count() == 0) then
      call abort(1)
      call MPI_Abort(MPI_COMM_WORLD,1)
    else
      call abort(8)
      call MPI_Abort(MPI_COMM_WORLD,8)
    end if
  end if

  call ReadCLA(MisConv,MisRank,MisRange,MinIter,MaxIter,DirName,Verbose,AEF,InOpt,TargMisM,TargMisK,GM_Types,MPIAbort)
  NumGM = size(GM_Types)

  if(MPIAbort) then
    call abort(9)
    call MPI_Abort(MPI_COMM_WORLD,9)
  end if
  if(InOpt=='Grid') then
    ListIN = .False.
    call DriverChk(MisRank,MPIAbort,GM_Types=GM_Types)
    if(MPIAbort) then
      call abort(6)
      call MPI_Abort(MPI_COMM_WORLD,6)
    end if
  elseif(InOpt=='List') then
    ListIN = .True.
  else
    MPIAbort = .True.
  end if
  if(MPIAbort) then
    call abort(1)
    call MPI_Abort(MPI_COMM_WORLD,1)
  end if

! ----- Input files -----
  if(ListIn) then
    call ReadListInput(ListVarNames,ListVarValues,MPIAbort)
    if(MPIAbort) then
      call abort(2)
      call MPI_Abort(MPI_COMM_WORLD,2)
    end if
    NumVar = size(ListVarValues,2)
    CellNumMods = size(ListVarValues,1)

    if(np-1 > CellNumMods) then ! Kill program if too many threads
      call abort(5)
      call MPI_Abort(MPI_COMM_WORLD,5)
    end if

  else
    call ReadGridInput(ListVarNames,VarValsMin,VarValsMax,NumMods,Cells,Levels,CellNumMods,ModRpt,RptMod,MPIAbort)
	TotNumMods = CellNumMods+(CellNumMods*Cells*(Levels-1))
    if(MPIAbort) then
      call abort(2)
      call MPI_Abort(MPI_COMM_WORLD,2)
    end if

    if(np-1 > CellNumMods) then ! Kill program if too many threads
      call abort(5)
      call MPI_Abort(MPI_COMM_WORLD,5)
    end if

    NumVar = size(ListVarNames)
    allocate(ListVarValues(CellNumMods,size(NumMods)))
    allocate(Misfits(CellNumMods))
    allocate(AllModels(CellNumMods,1+size(NumMods)))
    allocate(mask(1+size(NumMods)))
    allocate(KeepCells(Cells,1+size(NumMods)))
    allocate(StepOld(size(NumMods)))
    allocate(Step(size(NumMods)))
  end if

  call OrbDataCheck(ListVarNames,OData)
  call InputFileCheck(OData,MPIAbort)
  if(MPIAbort) then
    call abort(7)
    call MPI_Abort(MPI_COMM_WORLD,7)
  end if

end if ! (ThID == 0)

! ----- CLA Broadcasting -----
call MPI_Bcast(MisConv,     1, MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(MisRank,   3, MPI_CHARACTER,0, MPI_COMM_WORLD, ierr)
if(MisConv) call MPI_Bcast(MisRange,  1, MPI_REAL,     0, MPI_COMM_WORLD, ierr)

call MPI_Bcast(NumGM, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if(ThID /=0) allocate(GM_Types(NumGM))
do i = 1,NumGM
  call MPI_Bcast(GM_Types(i), 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
end do

call MPI_Bcast(MinIter,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(MaxIter,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(DirName,   100, MPI_CHARACTER,0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(Verbose,     1, MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(AEF,         1, MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(ListIn,      1, MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(CellNumMods,  1, MPI_integer,  0, MPI_COMM_WORLD, ierr)

! ----- Input Broadcasting -----
call MPI_Bcast(NumVar, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

if(ThID /=0) allocate(ListVarNames(NumVar),ListVarValues(1,NumVar))

do i = 1,NumVar
  call MPI_Bcast(ListVarNames(i), 10, MPI_character, 0, MPI_COMM_WORLD, ierr)
end do
allocate(allmisfits(7))
! -------------------- Boss --------------------
if(ThID == 0) then

! ---------- Grid model selector ----------
  if(.not. ListIN) then

! Make model output directory
    call execute_command_line('mkdir '//trim(DirName))
    call execute_command_line('cp INPUT/InputFiles.in '//trim(DirName))
    call execute_command_line('cp INPUT/GridInput.in '//trim(DirName))

    allocate(WorkList(np-1,2))

    Tag = 0
    ModWrt = 0
	
    do iLevel = 1,Levels
      print'(2(A,I0))','Starting level ',iLevel,' of ',Levels	
	  call cpu_time(start)
      
	  if(iLevel == 1) then
        N = 1
      else
        N = Cells
      end if

      if(iLevel == 2) deallocate(AllModels,mask)
      if(iLevel > 1) StepOld = Step
      do iCell = 1,N
        if(iLevel > 1) call UpdateSetup(KeepCells(iCell,1:size(KeepCells,2)-1),StepOld,VarValsMin,VarValsMax)
        call CalcGrid(VarValsMin,VarValsMax,NumMods,Grid,Step) ! Calculate the Grid of values
        call GenerateModels(NumMods,CellNumMods,Grid,ListVarValues) ! Generate Models from Grid
        recvd = 0
        do ModNum = 1,np-1
		  if(iLevel > 1 .and. ModRpt .and. ModNum == RptMod) then
            Tag = Tag + 1
			repeated = Tag	
			recvd = recvd + 1
            if(iLevel == 1) WorkList(ModNum,1) = ModNum
            WorkList(ModNum,2) = Tag-ModWrt	
			print'(A,I0,A)','Model ',repeated,' is being skipped as a repeated model'
            cycle
		  end if		
          call MPI_Send(.True., 1, MPI_LOGICAL, ModNum, ModNum, MPI_COMM_WORLD, ierr)
          Tag = Tag + 1
          if(iLevel == 1) WorkList(ModNum,1) = ModNum
          WorkList(ModNum,2) = Tag-ModWrt
          call MPI_Send(ListVarValues(ModNum,:), NumVar, MPI_double_precision, ModNum, Tag, MPI_COMM_WORLD, ierr)
          print'(2(A,I0))','Starting model: ',Tag,' of ',TotNumMods
		end do

! Write to output file
        if(iLevel==1) then
          open(unit=101,file=trim(DirName)//'/Models.txt',action='write')
          call get_command(arg)
          write(outmsg,'(A,A)') 'Program invoked with: ',trim(arg)
          write(101,'(A)') trim(outmsg)
        end if

! Find which Misfits are run
		call MisfitRun(MisRun)

! Send remaining models to free MPI threads as they become free
        if(np-1 <= size(ListVarValues,1)) then
          do ModNum = np,size(ListVarValues,1)
		    if(iLevel > 1 .and. ModRpt .and. ModNum == RptMod) then
              Tag = Tag + 1
			  repeated = Tag		  
			  recvd = recvd + 1
			  print'(A,I0,A)','Model ',repeated,' is being skipped as a repeated model'
              cycle
		    end if	  
            call MPI_Recv(allmisfits, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
            print'(A,I0)','Finished model: ',msg_stat(MPI_TAG)
			recvd = recvd + 1
            if(any(allmisfits == 100000)) then
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100000)
            elseif(any(allmisfits == 100001)) then
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100001)
            else
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank)
            end if
            if(MisRank == 'GV')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(1) ! Geodetic Velocity
            if(MisRank == 'SSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(2) ! Seafloor Spreading Rates
            if(MisRank == 'SD')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(3) ! Stress Direction
            if(MisRank == 'FSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(4) ! Fault Slip Rate
            if(MisRank == 'SC')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(5) ! Smoothed Seismicity Correlation
            if(MisRank == 'SA')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(6) ! Seismic Anisotropy
            if(MisRank == 'GM')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(7) ! Geometric Mean
            allmisfits = 0.0
            if(AEF) then
              FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
              if(FErrChk) then
                call abort(4)
                if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
                if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
                call MPI_Abort(MPI_COMM_WORLD,4)
              end if
            else
              FErrChk = FileExist('FatalError.txt')
              if(FErrChk) then
                call abort(3)
                if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
                call MPI_Abort(MPI_COMM_WORLD,3)
              end if
            end if

            Tag = Tag + 1			
            call MPI_Send(.True., 1, MPI_LOGICAL, msg_stat(MPI_SOURCE), Tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(ListVarValues(ModNum,:), NumVar, MPI_double_precision, msg_stat(MPI_SOURCE), Tag, MPI_COMM_WORLD, ierr) ! updated Tag
			print'(2(A,I0))','Starting model: ',Tag,' of ',TotNumMods
			WorkList(msg_stat(MPI_SOURCE),2) = Tag-ModWrt
          end do
! Receive last models
          do while(recvd/=size(ListVarValues,1))! recvd counts messages back
            call MPI_Recv(allmisfits, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
            print'(A,I0)','Finished model: ',msg_stat(MPI_TAG)
			recvd = recvd + 1
            if(any(allmisfits == 100000)) then
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100000)
            elseif(any(allmisfits == 100001)) then
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100001)
            else
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank)
            end if
            if(MisRank == 'GV')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(1) ! Geodetic Velocity
            if(MisRank == 'SSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(2) ! Seafloor Spreading Rates
            if(MisRank == 'SD')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(3) ! Stress Direction
            if(MisRank == 'FSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(4) ! Fault Slip Rate
            if(MisRank == 'SC')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(5) ! Smoothed Seismicity Correlation
            if(MisRank == 'SA')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(6) ! Seismic Anisotropy
            if(MisRank == 'GM')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(7) ! Geometric Mean
            allmisfits = 0.0
            if(AEF) then
              FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
              if(FErrChk) then
                call abort(4)
                if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
                if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
                call MPI_Abort(MPI_COMM_WORLD,4)
              end if
            else
              FErrChk = FileExist('FatalError.txt')
              if(FErrChk) then
                call abort(3)
                if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
                call MPI_Abort(MPI_COMM_WORLD,3)
              end if
            end if
          end do
        end if
        if(iLevel==2 .and. iCell==1) allocate(AllModels(Cells*CellNumMods,1+size(ListVarValues,2)),mask(1+size(ListVarValues,2)))
        AllModels(1+(iCell-1)*CellNumMods : iCell*CellNumMods, 1:size(AllModels,2)-1) = ListVarValues
        AllModels(1+(iCell-1)*CellNumMods : iCell*CellNumMods,   size(AllModels,2))   = Misfits
! Copy repeated model data
        if(iLevel > 1 .and. ModRpt) then
          AllModels(RptMod+(iCell-1)*CellNumMods,1:size(AllModels,2)-1) = KeepCells(iCell,1:size(KeepCells,2)-1)
          AllModels(RptMod+(iCell-1)*CellNumMods,  size(AllModels,2))   = KeepCells(iCell,  size(KeepCells,2))
          call writeRepeated(repeated,ThID,iLevel,iCell,ListVarNames,KeepCells(iCell,:))        
        end if
        ModWrt = ModWrt + size(ListVarValues,1)
      end do ! cells
! Decide which cells to keep
      call MisfitAnalysis(AllModels(:,size(AllModels,2)),Cells,MisRank,Keep)
      do i = 1,size(Keep)
        KeepCells(i,1:size(KeepCells,2)-1) = AllModels(Keep(i),1:size(AllModels,2)-1)
        KeepCells(i,  size(KeepCells,2))   = AllModels(Keep(i),  size(AllModels,2))
      end do
! Write target misfits to file
      if(TargMisM /= 0.0) then
        if(MisRank == 'SC') then
          if(iLevel == 1 .and. maxval(KeepCells(:,size(KeepCells,2)))>TargMisM) then
            call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,CellNumMods,Levels,Step,DirName,MisRank,iLevel)
          elseif(iLevel > 1 .and. maxval(KeepCells(:,size(KeepCells,2)))>TargMisM) then
            call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,Cells*CellNumMods,Levels,Step,DirName,MisRank,iLevel)
          end if
        else
          if(iLevel == 1 .and. minval(KeepCells(:,size(KeepCells,2)))<TargMisM) then
            call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,CellNumMods,Levels,Step,DirName,MisRank,iLevel)
          elseif(iLevel > 1 .and. minval(KeepCells(:,size(KeepCells,2)))<TargMisM) then
            call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,Cells*CellNumMods,Levels,Step,DirName,MisRank,iLevel)
          end if
        end if
      end if

      if(TargMisK /= 0.0) then
        if(MisRank == 'SC') then
          if(maxval(AllModels(:,size(AllModels,2)),mask=AllModels(:,size(AllModels,2))>0)>TargMisK) then
            write(101,'(A,F11.5,A)') "Misfit kill value reached (",TargMisK,") - program exiting early"
            exit
          end if
        else
          if(minval(AllModels(:,size(AllModels,2)),mask=AllModels(:,size(AllModels,2))>0)<TargMisK) then
            write(101,'(A,F11.5,A)') "Misfit kill value reached (",TargMisK,") - program exiting early"
            exit
          end if
        end if
      end if

      call cpu_time(finish)
	  print'(A,I0,A,I0,A,I2.2,F0.5/)','Level ',iLevel,' finished in (mm:ss): ',int((finish-start)/60.0),':',int(mod((finish-start),60.0)),mod((finish-start),60.0)-int(mod((finish-start),60.0))  
    end do ! Levels

    do ModNum = 1,np-1
      call MPI_Send(.False., 1, MPI_LOGICAL, ModNum, ModNum, MPI_COMM_WORLD, ierr)
    end do



  elseif(ListIN) then ! ---------- List model selector ----------

! Make model output directory
    call execute_command_line('mkdir '//trim(DirName))
    call execute_command_line('cp INPUT/InputFiles.in '//trim(DirName))
    call execute_command_line('cp INPUT/ListInput.in '//trim(DirName))

    if(np-1 <= size(ListVarValues,1)) then
! Send first np models to np MPI threads
      do ModNum = 1,np-1
        call MPI_Send(.True., 1, MPI_LOGICAL, ModNum, ModNum, MPI_COMM_WORLD, ierr)
        call MPI_Send(ListVarValues(ModNum,:), NumVar, MPI_double_precision,  ModNum, ModNum, MPI_COMM_WORLD, ierr)
        print'(2(A,I0))','Starting model: ',ModNum,' of ',size(ListVarValues,1)
	  end do

! Open output file
      open(unit=101,file=trim(DirName)//'/Models.txt',action='write')
      call get_command(arg)
      write(outmsg,'(A,A)') 'Program invoked with: ',trim(arg)
      write(101,'(A)') trim(outmsg)

! Find which Misfits are run
		call MisfitRun(MisRun)

      recvd = 0
! Send remaining models to MPI threads as they become free
      if(np-1 < size(ListVarValues,1)) then
        do ModNum = np,size(ListVarValues,1)
          call MPI_Recv(allmisfits, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
          print'(A,I0)','Finished model: ',msg_stat(MPI_TAG)
		  recvd = recvd+1
          if(any(allmisfits == 100000)) then
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100000)
          elseif(any(allmisfits == 100001)) then
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100001)
          else
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank)
          end if
          allmisfits = 0.0
          if(AEF) then
            FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
            if(FErrChk) then
              call abort(4)
              if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,4)
            end if
          else
            FErrChk = FileExist('FatalError.txt')
            if(FErrChk) then
              call abort(3)
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,3)
            end if
          end if
          call MPI_Send(.True., 1, MPI_LOGICAL, msg_stat(MPI_SOURCE), ModNum, MPI_COMM_WORLD, ierr)
          call MPI_Send(ListVarValues(ModNum,:), NumVar, MPI_double_precision,  msg_stat(MPI_SOURCE), ModNum, MPI_COMM_WORLD, ierr)
		  print'(2(A,I0))','Starting model: ',ModNum,' of ',size(ListVarValues,1)
		end do

        do while(recvd/=size(ListVarValues,1))! recvd counts messages back
          call MPI_Recv(allmisfits, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
          print'(A,I0)','Finished model: ',msg_stat(MPI_TAG)
		  recvd = recvd+1
          if(any(allmisfits == 100000)) then
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100000)
          elseif(any(allmisfits == 100001)) then
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100001)
          else
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank)
          end if
          allmisfits = 0.0
          if(AEF) then
            FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
            if(FErrChk) then
              call abort(4)
              if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,4)
            end if
          else
            FErrChk = FileExist('FatalError.txt')
            if(FErrChk) then
              call abort(3)
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,3)
            end if
          end if
          call MPI_Send(.False., 1, MPI_LOGICAL, msg_stat(MPI_SOURCE), 0, MPI_COMM_WORLD,ierr) ! Tell workers to stop
        end do

      elseif(np-1 == size(ListVarValues,1)) then
        do ModNum = 1,np-1
          call MPI_Recv(allmisfits, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
          print'(A,I0)','Finished model: ',msg_stat(MPI_TAG)
		  if(any(allmisfits == 100000)) then
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100000)
          elseif(any(allmisfits == 100001)) then
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank,Failed=100001)
          else
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisRank=MisRank)
          end if
          allmisfits = 0.0
          if(AEF) then
            FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
            if(FErrChk) then
              call abort(4)
              if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,4)
            end if
          else
            FErrChk = FileExist('FatalError.txt')
            if(FErrChk) then
              call abort(3)
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,3)
            end if
          end if
          call MPI_Send(.False., 1, MPI_LOGICAL, msg_stat(MPI_SOURCE), 0, MPI_COMM_WORLD,ierr) ! Tell workers to stop
        end do
      end if

    end if

! Close output file
    close(101)
  end if
end if




! -------------------- Worker --------------------
if(ThID/=0)then

! Worker array allocation
  if(MisConv) allocate(Misfits(MaxIter))

! MKL setup
  if(np-1 <= CellNumMods) then
    call cpuInfo(cpus)
    MKLthreads = nint( real(cpus)/real(np) ) 
  end if

! https://www.intel.com/content/www/us/en/develop/documentation/onemkl-windows-developer-guide/top/managing-performance-and-memory/improving-performance-with-threading.html
  call mkl_set_num_threads(MKLthreads) ! This is a suggestion, MKL may use fewer  (https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/support-functions/threading-control/mkl-set-num-threads.html)
  call mkl_set_dynamic(1) ! Dynamic true to allow fewer threads to be called      (https://www.intel.com/content/www/us/en/develop/documentation/onemkl-windows-developer-guide/top/managing-performance-and-memory/improving-performance-with-threading/using-additional-threading-control/mkl-dynamic.html)

  call MPI_Recv(Work, 1, MPI_LOGICAL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)! Wait for work instruction
  do while(Work)
    call MPI_Recv(ListVarValues, NumVar, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
    if(Verbose) then
      write(VerbFile,"(A,'/Verbose',I0,'.txt')") trim(DirName),ThID
      open(unit=iUnitVerb,file=trim(VerbFile))
    end if

    ModNum = msg_stat(MPI_TAG)

    call readDATA(DataFiles)
    write(filename,'(2A)') "INPUT/"//trim(DataFiles(2))

    open(unit=1,file=trim(filename))
    call ReadPm (     1, iUnitVerb, names,  nPlate, offMax, & ! INTENT(IN)
      &          aCreep, alphaT,    bCreep, Biot  , Byerly, & ! INTENT(OUT)
      &          cCreep, cFric,     conduc, dCreep, eCreep, &
      &          everyP, fFric,     gMean , gradie, iConve, &
      &          iPVRef, maxItr,    OKDelV, OKToQt, oneKm,  &
      &          radio , radius,    refStr, rhoAst, rhoBar, &
      &          rhoH2O, TAdiab,    tauMax, temLim, title3, &
      &          trHMax, TSurf,     vTimes, zBAsth, pltRef)
    close(1)

    Run = .True.
! Check iteration number & iConve value are valid mix
    if(MaxIter /= 0 .and. (iConve == 1 .or. iConve == 2 .or. iConve == 3 .or. iConve == 4 .or. iConve == 5)) then
      write(ErrorMsg,'(A,I0,A)') "The iConve value ",iConve," expects -Iter to be 0 or not entered"
      call FatalError(ErrorMsg,ModNum)
      call abort(10)
      if(FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
      call MPI_Abort(MPI_COMM_WORLD,10)
    elseif(MaxIter == 0 .and. (iConve == 0 .or. iConve == 6)) then
      write(ErrorMsg,'(A,I0,A)') "The iConve value ",iConve," expects -Iter to be >0"
      call FatalError(ErrorMsg,ModNum)
      call abort(10)
      if(FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
      call MPI_Abort(MPI_COMM_WORLD,10)
    end if

    call Variable_Update( fFric, cFric,  Biot,   Byerly, aCreep, &
      &                  bCreep, cCreep, dCreep, eCreep, tAdiab, &
      &                  gradie, zBAsth, trHMax, tauMax, rhoH2O, &
      &                  rhoBar, rhoAst, gMean,  oneKm,  radius, &
      &                  alphaT, conduc, radio,  tSurf,  temLim, &
      &                  ListVarNames, ListVarValues(1,:)) ! ListVarValues(1,:) because only relevant values sent to worker
    if(Verbose) then
      write(iUnitVerb,"(/A)")"The following variables were updated before being used by any of OrbData, Shells or OrbScore"
      write(frmt,"(A,I0,A)") "(",size(ListVarNames),"(A,X),/)"
      write(iUnitVerb,frmt), (trim(ListVarNames(i)), i =1,size(ListVarNames))

      write(iUnitVerb,"(A)") "With the following values:"
      call FormatStrings('MisL',ListVarNames,0,frmt)
      write(iUnitVerb,frmt), (ListVarValues(1,i), i =1,size(ListVarValues,2)),' '
      write(iUnitVerb,"(A/)"), "-------------------------------------------------------------------------------"
    end if

    call EarthChk(ModNum,ThID,ListVarNames,ListVarValues(1,:),Run)

    if(Run) then
      SHELLSconv = .True.
      if(MaxIter>0) call InputSetup(ThID,ModNum,"SH",DirName,iConve=iConve) ! Call before Data as OUT-IN handled by Data
      call InputSetup(ThID,ModNum,"SF",DirName,iConve=iConve) ! Call before Data as OUT-IN handled by Data
      call InputSetup(ThID,ModNum,"OS",DirName) ! Call before Data as OUT-IN handled by Data
      if(AEF) then
        FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
        if(FErrChk) then
          call abort(4)
          if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
          if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
          call MPI_Abort(MPI_COMM_WORLD,4)
        end if
      else
        FErrChk = FileExist('FatalError.txt')
        if(FErrChk) then
          call abort(3)
          if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
          call MPI_Abort(MPI_COMM_WORLD,3)
        end if
      end if

	  call OrbDataCheck(ListVarNames,OData)
      if(OData) then
        call InputSetup(ThID,ModNum,"OD",DirName)
        call OpenInput(ThID,ModNum,"OD",DirName)
        call OpenOutput(ThID,ModNum,"OD",DirName)
        call OrbData5(ModNum,ListVarNames,ListVarValues(1,:),DirName)
        if(AEF) then
          FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
          if(FErrChk) then
            call abort(4)
            if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
            if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
            call MPI_Abort(MPI_COMM_WORLD,4)
          end if
        else
          FErrChk = FileExist('FatalError.txt')
          if(FErrChk) then
            call abort(3)
            if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
            call MPI_Abort(MPI_COMM_WORLD,3)
          end if
        end if
	    call CloseInput("OD")
        call CloseOutput(ThID,ModNum,"OD",DirName,Iter=MaxIter) ! Copies FEG file to SHELLS & OrbScore input
      end if

      rpeat=0
      do while(rpeat<MaxIter .and. SHELLSconv)
        rpeat = rpeat+1
! Shells
        call OpenInput(ThID,ModNum,"SH",DirName,rpeat=rpeat)
        call OpenOutput(ThID,ModNum,"SH",DirName,rpeat=rpeat)
        call Shells_v5p0(ListVarNames,ListVarValues(1,:),rpeat,SHELLSconv)
        if(AEF) then
          FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
          if(FErrChk) then
            call abort(4)
            if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
            if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
            call MPI_Abort(MPI_COMM_WORLD,4)
          end if
        else
          FErrChk = FileExist('FatalError.txt')
          if(FErrChk) then
            call abort(3)
            if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
            call MPI_Abort(MPI_COMM_WORLD,3)
          end if
        end if
		call CloseInput("SH")
        call CloseOutput(ThID,ModNum,"SH",DirName,rpeat=rpeat)
! Misfit Convergence
        if(rpeat>=MinIter .and. SHELLSconv .and. MisConv) then
          call OpenInput(ThID,ModNum,"OS",DirName,rpeat=rpeat,Final=.False.)
          call OpenOutput(ThID,ModNum,"OS",DirName,plt=plt,rpeat=rpeat,MC=MisConv)
          call OrbScore2(ListVarNames,ListVarValues(1,:),allmisfits)
          if(AEF) then
            FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
            if(FErrChk) then
              call abort(4)
              if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,4)
            end if
          else
            FErrChk = FileExist('FatalError.txt')
            if(FErrChk) then
              call abort(3)
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,3)
            end if
          end if

		  call Closeinput("OS")
          call CloseOutput(ThID,ModNum,"OS",DirName)

! Calculate Geometric mean
          if(trim(MisRank) == 'GM') then
            allmisfits(7) = GeoMean(allmisfits,GM_Types=GM_Types)
          else	  
            allmisfits(7) = GeoMean(allmisfits)
          end if

		  if(MisRank == 'GV')  Misfits(rpeat) = allmisfits(1) ! Geodetic Velocity
          if(MisRank == 'SSR') Misfits(rpeat) = allmisfits(2) ! Seafloor Spreading Rates
          if(MisRank == 'SD')  Misfits(rpeat) = allmisfits(3) ! Stress Direction
          if(MisRank == 'FSR') Misfits(rpeat) = allmisfits(4) ! Fault Slip Rate
          if(MisRank == 'SC')  Misfits(rpeat) = allmisfits(5) ! Smoothed Seismicity Correlation
          if(MisRank == 'SA')  Misfits(rpeat) = allmisfits(6) ! Seismic Anisotropy
          if(MisRank == 'GM')  Misfits(rpeat) = allmisfits(7) ! Geometric Mean

          if(rpeat>MinIter) then
            misfitMin = (1.0 - MisRange) * Misfits(rpeat-1)
            misfitMax = (1.0 + MisRange) * Misfits(rpeat-1)
            if(Misfits(rpeat)>misfitMin .and. Misfits(rpeat)<misfitMax) then
			  if(Verbose) then
			    write(iUnitVerb,'(A,I0,3(A,f10.4))') 'Model ',ModNum,' Misfit convergence satisfied: ',Misfits(rpeat),' is within', 100*MisRange,'% of ',Misfits(rpeat-1)
			    print'(2(A,I0))', 'Model ',ModNum,' Misfit convergence satisfied at iteration ',rpeat
			  end if
			  exit ! Misfit convergence satisfied - exit main work loop
			end if
          end if
        end if

! Update variables between 1st & 2nd iteration?
        if(rpeat==1 .and. FileExist('INPUT/UpVar.in')) then
          if(ThID==1 .and. Verbose) write(iUnitVerb,'(A)') 'UpVar.in file detected, updating listed variables'
          call IterVar(fFric, cFric,  Biot,   Byerly, aCreep, &
          &           bCreep, cCreep, dCreep, eCreep, tAdiab, &
          &           gradie, zBAsth, trHMax, tauMax, rhoH2O, &
          &           rhoBar, rhoAst, gMean,  oneKm,  radius, &
          &           alphaT, conduc, radio,  tSurf,  temLim)
        end if


      end do ! do while(rpeat<MaxIter .and. SHELLSconv)

! Final or non-iterating Shells call
      if(SHELLSconv) then ! Final call with (optional) different files
        if(MaxIter /= 0) rpeat = rpeat+1 ! necessary since rpeat updated at beginning of main loop
        call OpenInput(ThID,ModNum,"SF",DirName,rpeat=rpeat)
        call OpenOutput(ThID,ModNum,"SF",DirName,rpeat=rpeat)
        call Shells_v5p0(ListVarNames,ListVarValues(1,:),rpeat,SHELLSconv)
        if(AEF) then
          FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
          if(FErrChk) then
            call abort(4)
            if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
            if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
            call MPI_Abort(MPI_COMM_WORLD,4)
          end if
        else
          FErrChk = FileExist('FatalError.txt')
          if(FErrChk) then
            call abort(3)
            if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
            call MPI_Abort(MPI_COMM_WORLD,3)
          end if
        end if
        call CloseInput("SF")
        call CloseOutput(ThID,ModNum,"SF",DirName,rpeat=rpeat)

        if(SHELLSconv) then
          call OpenInput(ThID,ModNum,"OS",DirName,rpeat=rpeat,Final=.True.)
          call OpenOutput(ThID,ModNum,"OS",DirName,plt=plt,rpeat=rpeat,MC=MisConv)
          call OrbScore2(ListVarNames,ListVarValues(1,:),allmisfits)
          if(AEF) then
            FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
            if(FErrChk) then
              call abort(4)
              if (FileExist('ModelError.txt')) call execute_command_line('rm ModelError.txt')
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,4)
            end if
          else
            FErrChk = FileExist('FatalError.txt')
            if(FErrChk) then
              call abort(3)
              if (FileExist('FatalError.txt')) call execute_command_line('rm FatalError.txt')
              call MPI_Abort(MPI_COMM_WORLD,3)
            end if
          end if
          call Closeinput("OS")
          call CloseOutput(ThID,ModNum,"OS",DirName)
        end if
      end if

      if(.not. SHELLSconv) then
        allmisfits = 100000 ! SHELLSconv False
      end if
    else ! Run False
      allmisfits = 100001 ! Differentiate between failures
    end if

! Calculate Geometric mean
      if(trim(MisRank) == 'GM') then
        allmisfits(7) =  GeoMean(allmisfits,GM_Types=GM_Types)
      else
        allmisfits(7) = GeoMean(allmisfits)
      end if

    call MPI_Send(allmisfits, 7, MPI_DOUBLE_PRECISION, 0, ModNum, MPI_COMM_WORLD,ierr)
    call MPI_Recv(Work, 1, MPI_LOGICAL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)

  end do ! do while(Work) loop

end if ! Boss/Worker selection

call MPI_Barrier(MPI_COMM_WORLD, ierr)

! Tidy up
if(ThID /=0) then 
  if(Verbose) close(iUnitVerb)

else
  write(cmd,'(3A,I0,A)') "rm -r ",trim(DirName),"/ThID_*_input"
  call execute_command_line(trim(cmd))

  call cpu_time(ProgFinish)
  print "(A,A,A,I0,A,I2.2,F0.5/)", "ShellSet ",trim(InOpt)," test finished in (mm:ss): ",int((ProgFinish-ProgStart)/60.0),':',int(mod((ProgFinish-ProgStart),60.0)),mod((ProgFinish-ProgStart),60.0)-int(mod((ProgFinish-ProgStart),60.0))
  print "(A/,A/,A/,A/,A/,A/,A/)", "****************************************************************************", &
    &    "ShellSet Copyright (C) 2023 Jon May, Peter Bird, Michele Carafa", &
    &    "This program comes with ABSOLUTELY NO WARRANTY.", &
    &    "This is free software, you are welcome to redistribute it under conditions", &
    &    "set out in GNU General Public License version 3, or any later version.", &
    &    "See the included license file, or https://www.gnu.org/licenses/, for details", &
    &    "****************************************************************************"

end if

call MPI_FINALIZE (ierr) ! Finalize communication

end program