!*******************************************************************************
! ShellSet program
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

program ShellSet

use mkl_service
use ShellSetSubs
use SHELLS
use OrbScore
use OrbData
use SharedVars
use VariableCheck


implicit none

include "mpif.h"  ! MPI library


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
integer :: TotNumMods,Cells,Levels,RptMod,iLevel,iCell,N,Tag,NumMis,ModWrt = 0
logical :: ModRpt,Kill
real :: TargMisK,TargMisM
logical,dimension(:),allocatable :: mask

!MKL setup variables
integer :: cpus,MKLthreads

!MPI specific variables
integer :: np,ierr,msg_stat(MPI_STATUS_size)

! Command line argument Variables
integer :: MinIter,MaxIter
logical :: MisConv,AEF,MPIAbort,CLA_STOP=.False.,ListIn ! True if list input, False if grid search
real :: MisRange
character(len=3) :: MisType
character(len=100) :: VerbFile,DirName
character(len=4) :: InOpt

! Work iteration & control
integer :: ModNum,rpeat,recvd
real*8 :: misfitMin,misfitMax
real*8,dimension(:),allocatable :: Misfits
real*8,dimension(6) :: allmisfits
logical :: OData,Run,Work,SHELLSconv,FErrChk=.False.

character(len=100) :: frmt,filename
character(len=200) :: arg,outmsg,cmd
logical,dimension(6) :: MisRun
character(len=100),dimension(6) :: DataFiles
real*8 :: start,finish

!------------------------------------------ Executing code ----------------------------------------------
! ----- Start MPI -----
call MPI_INIT (ierr)
call MPI_COMM_size (MPI_COMM_WORLD, np, ierr)
call MPI_COMM_RANK (MPI_COMM_WORLD, ThID, ierr)


if(ThID == 0) then
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

  call ReadCLA(MisConv,MisType,MisRange,MinIter,MaxIter,DirName,Verbose,AEF,InOpt,TargMisM,TargMisK,MPIAbort)
  if(MPIAbort) then
    call abort(9)
    call MPI_Abort(MPI_COMM_WORLD,9)
  end if  
  if(InOpt=='Grid') then
    ListIN = .False.
    call DriverChk(MisType,MPIAbort)
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
	TotNumMods = size(ListVarValues,1)

    if(np-1 > TotNumMods) then ! Kill program if too many threads
      call abort(5)
      call MPI_Abort(MPI_COMM_WORLD,5)
    end if

  else
    call ReadGridInput(ListVarNames,VarValsMin,VarValsMax,NumMods,Cells,Levels,TotNumMods,ModRpt,RptMod,MPIAbort)
	if(MPIAbort) then
	  call abort(2)
	  call MPI_Abort(MPI_COMM_WORLD,2)
	end if

    if(np-1 > TotNumMods) then ! Kill program if too many threads
      call abort(5)
      call MPI_Abort(MPI_COMM_WORLD,5)
    end if

	NumVar = size(ListVarNames)
    allocate(ListVarValues(TotNumMods,size(NumMods)))
    allocate(Misfits(TotNumMods))
    allocate(AllModels(TotNumMods,1+size(NumMods)))
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
if(MisConv) then
  call MPI_Bcast(MisType,   3, MPI_CHARACTER,0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(MisRange,  1, MPI_REAL,     0, MPI_COMM_WORLD, ierr)
end if
call MPI_Bcast(MinIter,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(MaxIter,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(DirName,   100, MPI_CHARACTER,0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(Verbose,     1, MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(AEF,         1, MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(ListIn,      1, MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(TotNumMods,  1, MPI_integer,  0, MPI_COMM_WORLD, ierr)

! ----- Input Broadcasting -----
call MPI_Bcast(NumVar, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
if(ThID /=0) allocate(ListVarNames(NumVar),ListVarValues(1,NumVar))

do i = 1,NumVar
  call MPI_Bcast(ListVarNames(i), 10, MPI_character, 0, MPI_COMM_WORLD, ierr)
end do

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
        call GenerateModels(NumMods,TotNumMods,Grid,ListVarValues) ! Generate Models from Grid

        do ModNum = 1,np-1
          call MPI_Send(.True., 1, MPI_LOGICAL, ModNum, ModNum, MPI_COMM_WORLD, ierr)
	      Tag = Tag + 1
	      if(iLevel == 1) WorkList(ModNum,1) = ModNum
	      WorkList(ModNum,2) = Tag-ModWrt
	      call MPI_Send(ListVarValues(ModNum,:), NumVar, MPI_double_precision, ModNum, Tag, MPI_COMM_WORLD, ierr)
        end do

! Write to output file
        if(iLevel==1) then
          open(unit=101,file=trim(DirName)//'/Models.txt',action='write')
          call get_command(arg)
          write(outmsg,'(A,A)') 'Program invoked with: ',trim(arg)
          write(101,'(A)') trim(outmsg)
        end if

        recvd = 0
! Send remaining models to free MPI threads as they become free
        if(np-1 <= size(ListVarValues,1)) then
          do ModNum = np,size(ListVarValues,1)
            call MPI_Recv(allmisfits, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
			call MPI_Recv(MisRun, 6, MPI_LOGICAL, msg_stat(MPI_SOURCE), MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
			recvd = recvd+1
			if(MisType == 'GV')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(1) ! Geodetic Velocity
	        if(MisType == 'SSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(2) ! Seafloor Spreading Rates
	        if(MisType == 'SD')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(3) ! Stress Direction
	        if(MisType == 'FSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(4) ! Fault Slip Rate
	        if(MisType == 'SC')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(5) ! Smoothed Seismicity Correlation
	        if(MisType == 'SA')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(6) ! Seismic Anisotropy
		    if(any(allmisfits == 100000)) then
			  call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100000)
		    elseif(any(allmisfits == 100001)) then
			  call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100001)
		    else
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType)
	        end if
			allmisfits = 0.0
			if(AEF) then
	          FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		      if(FErrChk) then
		        call abort(4)
		        call MPI_Abort(MPI_COMM_WORLD,4)
	          end if
	        else
	          FErrChk = FileExist('FatalError.txt')
		      if(FErrChk) then
			    call abort(3)
			    call MPI_Abort(MPI_COMM_WORLD,3)
	          end if
	        end if

		    Tag = Tag + 1
	        call MPI_Send(.True., 1, MPI_LOGICAL, msg_stat(MPI_SOURCE), Tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(ListVarValues(ModNum,:), NumVar, MPI_double_precision, msg_stat(MPI_SOURCE), Tag, MPI_COMM_WORLD, ierr) ! updated tag
		    WorkList(msg_stat(MPI_SOURCE),2) = Tag-ModWrt
	      end do

          do while(recvd/=size(ListVarValues,1))! recvd counts messages back
            call MPI_Recv(allmisfits, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
			call MPI_Recv(MisRun, 6, MPI_LOGICAL, msg_stat(MPI_SOURCE), MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
	        recvd = recvd+1
            if(MisType == 'GV')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(1) ! Geodetic Velocity
	        if(MisType == 'SSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(2) ! Seafloor Spreading Rates
	        if(MisType == 'SD')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(3) ! Stress Direction
	        if(MisType == 'FSR') Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(4) ! Fault Slip Rate
	        if(MisType == 'SC')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(5) ! Smoothed Seismicity Correlation
	        if(MisType == 'SA')  Misfits(WorkList(msg_stat(MPI_SOURCE),2)) = allmisfits(6) ! Seismic Anisotropy

		    if(any(allmisfits == 100000)) then
			  call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100000)
		    elseif(any(allmisfits == 100001)) then
  			  call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100001)
		    else
              call WriteOut('Grid',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType)
	        end if
			allmisfits = 0.0
			  if(AEF) then
	            FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		        if(FErrChk) then
			      call abort(4)
			      call MPI_Abort(MPI_COMM_WORLD,4)
	            end if
	          else
  	            FErrChk = FileExist('FatalError.txt')
    		    if(FErrChk) then
	              call abort(3)
			      call MPI_Abort(MPI_COMM_WORLD,3)
	            end if
	          end if
	      end do
        end if

        if(iLevel==2 .and. iCell==1) allocate(AllModels(Cells*TotNumMods,1+size(ListVarValues,2)),mask(1+size(ListVarValues,2)))
        AllModels(1+(iCell-1)*TotNumMods : iCell*TotNumMods, 1:size(AllModels,2)-1) = ListVarValues
  	    AllModels(1+(iCell-1)*TotNumMods : iCell*TotNumMods,   size(AllModels,2))   = Misfits
! Copy repeated model data
	    if(iLevel > 1 .and. ModRpt) then
          AllModels(RptMod+(iCell-1)*TotNumMods,1:size(AllModels,2)-1) = KeepCells(iCell,1:size(KeepCells,2)-1)
   	      AllModels(RptMod+(iCell-1)*TotNumMods,  size(AllModels,2))   = KeepCells(iCell,  size(KeepCells,2))
		end if
	    ModWrt = ModWrt + size(ListVarValues,1)
      end do ! cells

! Decide which cells to keep
      call MisfitAnalysis(AllModels(:,size(AllModels,2)),Cells,MisType,Keep)
      do i = 1,size(Keep)
        KeepCells(i,1:size(KeepCells,2)-1) = AllModels(Keep(i),1:size(AllModels,2)-1)
        KeepCells(i,  size(KeepCells,2))   = AllModels(Keep(i),  size(AllModels,2))
      end do
! Write target misfits to file
      if(TargMisM /= 0.0) then
        if(MisType == 'SC') then
          if(iLevel == 1 .and. maxval(KeepCells(:,size(KeepCells,2)))>TargMisM) then
	        call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,TotNumMods,Levels,Step,DirName,MisType,iLevel)
          elseif(iLevel > 1 .and. maxval(KeepCells(:,size(KeepCells,2)))>TargMisM) then
	        call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,Cells*TotNumMods,Levels,Step,DirName,MisType,iLevel)
          end if
        else
          if(iLevel == 1 .and. minval(KeepCells(:,size(KeepCells,2)))<TargMisM) then
	        call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,TotNumMods,Levels,Step,DirName,MisType,iLevel)
          elseif(iLevel > 1 .and. minval(KeepCells(:,size(KeepCells,2)))<TargMisM) then
	        call WriteMisfitTarget(AllModels,ListVarNames,TargMisM,Cells*TotNumMods,Levels,Step,DirName,MisType,iLevel)
          end if
        end if
      end if

! Check Kill misfit value at end of each level
      if(TargMisK /= 0) then
        if(minval(AllModels(:,size(AllModels,2)),mask=AllModels(:,size(AllModels,2))>0)<TargMisK) then
	      write(101,'(A,F11.5,A)') "Misfit kill value reached (",TargMisK,") - program exiting early"
	      exit
	    end if
      end if

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
      end do

! Open output file
      open(unit=101,file=trim(DirName)//'/Models.txt',action='write')
      call get_command(arg)
      write(outmsg,'(A,A)') 'Program invoked with: ',trim(arg)
      write(101,'(A)') trim(outmsg)

      recvd = 0
! Send remaining models to free MPI threads as they become free
      if(np-1 < size(ListVarValues,1)) then
        do ModNum = np,size(ListVarValues,1)
          call MPI_Recv(allmisfits, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
		  call MPI_Recv(MisRun, 6, MPI_LOGICAL, msg_stat(MPI_SOURCE), MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
          recvd = recvd+1
		  if(any(allmisfits == 100000)) then
			call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100000)
		  elseif(any(allmisfits == 100001)) then
			call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100001)
		  else
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType)
	      end if
	      allmisfits = 0.0
		  if(AEF) then
	        FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		    if(FErrChk) then
			  call abort(4)
			  call MPI_Abort(MPI_COMM_WORLD,4)
	        end if
	      else
	        FErrChk = FileExist('FatalError.txt')
		    if(FErrChk) then
			  call abort(3)
			  call MPI_Abort(MPI_COMM_WORLD,3)
	        end if
	      end if
	      call MPI_Send(.True., 1, MPI_LOGICAL, msg_stat(MPI_SOURCE), ModNum, MPI_COMM_WORLD, ierr)
          call MPI_Send(ListVarValues(ModNum,:), NumVar, MPI_double_precision,  msg_stat(MPI_SOURCE), ModNum, MPI_COMM_WORLD, ierr)
	    end do

        do while(recvd/=size(ListVarValues,1))! recvd counts messages back
          call MPI_Recv(allmisfits, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
          call MPI_Recv(MisRun, 6, MPI_LOGICAL, msg_stat(MPI_SOURCE), MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
          recvd = recvd+1
		  if(any(allmisfits == 100000)) then
			call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100000)
		  elseif(any(allmisfits == 100001)) then
			call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100001)
		  else
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType)
	      end if
		  allmisfits = 0.0
		  if(AEF) then
	        FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		    if(FErrChk) then
			  call abort(4)
			  call MPI_Abort(MPI_COMM_WORLD,4)
	        end if
	      else
	        FErrChk = FileExist('FatalError.txt')
		    if(FErrChk) then
			  call abort(3)
			  call MPI_Abort(MPI_COMM_WORLD,3)
	        end if
	      end if
	  call MPI_Send(.False., 1, MPI_LOGICAL, msg_stat(MPI_SOURCE), 0, MPI_COMM_WORLD,ierr) ! Tell workers to stop
    end do

      elseif(np-1 == size(ListVarValues,1)) then
        do ModNum = 1,np-1
          call MPI_Recv(allmisfits, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
          call MPI_Recv(MisRun, 6, MPI_LOGICAL, msg_stat(MPI_SOURCE), MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
		  if(any(allmisfits == 100000)) then
			call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100000)
		  elseif(any(allmisfits == 100001)) then
			call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType,Failed=100001)
		  else
            call WriteOut('List',ListVarNames,msg_stat(MPI_TAG),msg_stat(MPI_SOURCE),iLevel,iCell,ListVarValues(msg_stat(MPI_TAG)-ModWrt,:),allmisfits,MisRun,MisType=MisType)
	      end if
	      allmisfits = 0.0
		  if(AEF) then
	        FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		    if(FErrChk) then
		      call abort(4)
			  call MPI_Abort(MPI_COMM_WORLD,4)
	        end if
	      else
	        FErrChk = FileExist('FatalError.txt')
		    if(FErrChk) then
			  call abort(3)
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
  if(np-1 <= TotNumMods) then ! use np-1 to calculate MKL thread numbers
    call cpuInfo(cpus)
    MKLthreads = nint(real( (real(cpus)/max( real((np-1)),real(size(ListVarValues,1)) )) ))
  end if

  call mkl_set_num_threads(MKLthreads) ! This is a suggestion, MKL may use less  (https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/support-functions/threading-control/mkl-set-num-threads.html)
  call mkl_set_dynamic(1) ! Dynamic true to allow less threads to be called      (https://www.intel.com/content/www/us/en/develop/documentation/onemkl-windows-developer-guide/top/managing-performance-and-memory/improving-performance-with-threading/using-additional-threading-control/mkl-dynamic.html)

  call MPI_Recv(Work, 1, MPI_LOGICAL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)! Wait for work instruction
  do while(Work)
    call MPI_Recv(ListVarValues, NumVar, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)
	if(Verbose) then
      write(VerbFile,"(A,'/Verbose',I0,'.txt')") trim(DirName),ThID
      open(unit=iUnitVerb,file=trim(VerbFile))
    end if
    rpeat=0
	ModNum = msg_stat(MPI_TAG)
    Run = .False.

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

	call Variable_Update( fFric, cFric,  Biot,   Byerly, aCreep, &
      &                  bCreep, cCreep, dCreep, eCreep, tAdiab, &
      &                  gradie, zBAsth, trHMax, tauMax, rhoH2O, &
      &                  rhoBar, rhoAst, gMean,  oneKm,  radius, &
      &                  alphaT, conduc, radio,  tSurf,  temLim, &
      &                  ListVarNames, ListVarValues(1,:))

    call EarthChk(ModNum,ThID,ListVarNames,ListVarValues(1,:),Run)

    if(Run) then
	  SHELLSconv = .True.
	  do while(rpeat<MaxIter .and. SHELLSconv)
	    rpeat = rpeat+1
        if(rpeat==1) then
		  call InputSetup(ThID,ModNum,"SH",DirName) ! Call before Data as OUT-IN handled by Data
		  call InputSetup(ThID,ModNum,"OS",DirName,MisRun=MisRun) ! Call before Data as OUT-IN handled by Data
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
			    call MPI_Abort(MPI_COMM_WORLD,4)
	          end if
	        else
	          FErrChk = FileExist('FatalError.txt')
		      if(FErrChk) then
			    call abort(3)
			    call MPI_Abort(MPI_COMM_WORLD,3)
	          end if
	        end if
            call CloseInput("OD")
            call CloseOutput(ThID,ModNum,"OD",DirName) ! Copies FEG file to SHELLS & OrbScore input
          end if
        end if

    	call OpenInput(ThID,ModNum,"SH",DirName,rpeat=rpeat)
		call OpenOutput(ThID,ModNum,"SH",DirName,rpeat=rpeat)
		call Shells_v5p0(ListVarNames,ListVarValues(1,:),rpeat,SHELLSconv,OData)
	    if(AEF) then
	      FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		  if(FErrChk) then
		    call abort(4)
		    call MPI_Abort(MPI_COMM_WORLD,4)
	      end if
	    else
	      FErrChk = FileExist('FatalError.txt')
		  if(FErrChk) then
		    call abort(3)
			call MPI_Abort(MPI_COMM_WORLD,3)
	      end if
	    end if
		call CloseInput("SH")
		call CloseOutput(ThID,ModNum,"SH",DirName,rpeat=rpeat) ! Copy vEarth to OrbScore input

    	if(rpeat>=MinIter .and. SHELLSconv .and. MisConv) then
          call OpenInput(ThID,ModNum,"OS",DirName,rpeat=rpeat)
  		  call OpenOutput(ThID,ModNum,"OS",DirName,plt=(/.True.,.True.,.True./),rpeat=rpeat)
		  call OrbScore2(ListVarNames,ListVarValues(1,:),allmisfits)
	      if(AEF) then
	        FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		    if(FErrChk) then
			  call abort(4)
			  call MPI_Abort(MPI_COMM_WORLD,4)
	        end if
	      else
	        FErrChk = FileExist('FatalError.txt')
		    if(FErrChk) then
			  call abort(3)
			  call MPI_Abort(MPI_COMM_WORLD,3)
	        end if
	      end if
		  call Closeinput("OS")
          call CloseOutput(ThID,ModNum,"OS",DirName)

	      if(MisType == 'GV')  Misfits(rpeat) = allmisfits(1) ! Geodetic Velocity
	      if(MisType == 'SSR') Misfits(rpeat) = allmisfits(2) ! Seafloor Spreading Rates
	      if(MisType == 'SD')  Misfits(rpeat) = allmisfits(3) ! Stress Direction
	      if(MisType == 'FSR') Misfits(rpeat) = allmisfits(4) ! Fault Slip Rate
	      if(MisType == 'SC')  Misfits(rpeat) = allmisfits(5) ! Smoothed Seismicity Correlation
	      if(MisType == 'SA')  Misfits(rpeat) = allmisfits(6) ! Seismic Anisotropy

          if(rpeat>MinIter) then
	        misfitMin = (1.0 - MisRange) * Misfits(rpeat-1)
		    misfitMax = (1.0 + MisRange) * Misfits(rpeat-1)
		    if(Misfits(rpeat)>misfitMin .and. Misfits(rpeat)<misfitMax) exit ! Misfit convergence satisfied - exiting main loop
          end if
	    end if
	  end do ! do while(rpeat<MaxIter .and. SHELLSconv)

	  if(SHELLSconv) then ! Final call with (optional) different files
	    rpeat = rpeat+1 ! necessary since rpeat updated at beginning of main loop
	    call InputSetup(ThID,ModNum,"BF",DirName,Verbose=Verbose) ! Updates files in SHELLS input dir
        call OpenInput(ThID,ModNum,"SH",DirName,rpeat=rpeat)
        call OpenOutput(ThID,ModNum,"SH",DirName,rpeat=rpeat)
	    call Shells_v5p0(ListVarNames,ListVarValues(1,:),rpeat,SHELLSconv,OData)
	    if(AEF) then
	      FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		  if(FErrChk) then
		    call abort(4)
		    call MPI_Abort(MPI_COMM_WORLD,4)
		  end if
	    else
    	  FErrChk = FileExist('FatalError.txt')
		  if(FErrChk) then
		    call abort(3)
		    call MPI_Abort(MPI_COMM_WORLD,3)
	      end if
	    end if
        call CloseInput("SH")
        call CloseOutput(ThID,ModNum,"SH",DirName,rpeat=rpeat)

        if(SHELLSconv) then
          call OpenInput(ThID,ModNum,"OS",DirName,rpeat=rpeat)
          call OpenOutput(ThID,ModNum,"OS",DirName,plt=(/.True.,.True.,.True./),rpeat=rpeat)
		  call OrbScore2(ListVarNames,ListVarValues(1,:),allmisfits)
	      if(AEF) then
	        FErrChk = (FileExist('FatalError.txt') .or. FileExist('ModelError.txt'))
		    if(FErrChk) then
			  call abort(4)
			  call MPI_Abort(MPI_COMM_WORLD,4)
			end if
	      else
	        FErrChk = FileExist('FatalError.txt')
		    if(FErrChk) then
			  call abort(3)
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
	  call MisfitRun(MisRun)
	  allmisfits = 100001 ! SHELLS fails use 100000 to differentiate between failures
	end if

	call MPI_Send(allmisfits, 6, MPI_DOUBLE_PRECISION, 0, ModNum, MPI_COMM_WORLD,ierr)
	call MPI_Send(MisRun, 6, MPI_LOGICAL, 0, ModNum, MPI_COMM_WORLD,ierr)
	call MPI_Recv(Work, 1, MPI_LOGICAL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, msg_stat, ierr)

  end do ! do while(Work) loop

end if ! Boss/Worker selection

call MPI_Barrier(MPI_COMM_WORLD, ierr)
if(ThID /=0) then ! Tidy up
  if(Verbose) close(6)
  if(.not. Work) then
	if(OData) then
	  write(cmd,'(3A,I0,3A,I0,3A,I0,A)') "rm -r ",trim(DirName),"/ThID_",ThID,"_Data_input"
	  call execute_command_line(trim(cmd))
	end if
	write(cmd,'(3A,I0,3A,I0,3A,I0,A)') "rm -r ",trim(DirName),"/ThID_",ThID,"_Shells_input ",trim(DirName),"/ThID_",ThID,"_Score_input"
	call execute_command_line(trim(cmd))
  end if

else
  print "(A,A,A/)", "ShellSet ",trim(InOpt)," test finished."
  write(frmt,"(A,I0,A)") "(A,",size(ListVarNames),"(A,X),/)"
  print(frmt), "Variables tested: ", (trim(ListVarNames(i)), i =1,size(ListVarNames))
  print "(A/,A/,A/,A/,A/,A/,A/)", "**********************************************************************************", &
	&    "ShellSet Copyright (C) 2022 Jon May, Peter Bird, Michele Carafa", &
	&    "This program comes with ABSOLUTELY NO WARRANTY.", &
	&    "This is free software, you are welcome to redistribute it under conditions", &
	&    "set out in GNU General Public License version 3, or any later version.", &
	&    "See the included license file, or https://www.gnu.org/licenses/, for details.", &
	&    "**********************************************************************************"

end if

call MPI_FINALIZE (ierr) ! Finalize communication


end program
