!*******************************************************************************
! Module containing Variable check subroutines
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

module VariableCheck

use ShellSetSubs
use SharedVars

contains


subroutine EarthChk(ModNum,ThID,VarNames,VarValues,Run)

implicit none
integer,intent(in) ::  ModNum,ThID
character(len=10),dimension(:),intent(in) :: VarNames
real*8,dimension(:),intent(in) :: VarValues
logical,intent(out) :: Run

integer :: i,j
logical :: tested
logical,dimension(:),allocatable :: mask


allocate(mask(size(VarValues)))

mask = .False.
Run = .True.

do i = 1,size(VarNames)
! rhoBar_M > rhoBar_C?
  if(trim(VarNames(i)) == 'rhoBar_C' .and. (.not. mask(i))) then
    tested = .False.
    mask(i) = .True.
	do j = i+1,size(VarNames)
	  if(trim(VarNames(j)) == 'rhoBar_M') then
	    tested = .True.
	    mask(j) = .True.
	    if(VarValues(j) <= VarValues(i)) then ! if(rhoBar_M <= rhoBar_C) then
		  Run = .False.
		  write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with rhoBar values for model: ",ModNum, VarValues(i), VarValues(j)
		  call NonFatalError(ErrorMsg,ModNum)
		end if
	  end if
	end do

	if(.not. tested) then
	  if(rhoBar(2) <= VarValues(i)) then ! if(rhoBar_M <= rhoBar_C) then
		Run = .False.
		write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with rhoBar values for model: ",ModNum, VarValues(i), rhoBar(2)
		call NonFatalError(ErrorMsg,ModNum)
	  end if
	end if


  elseif(trim(VarNames(i)) == 'rhoBar_M' .and. (.not. mask(i))) then
    tested = .False.
    mask(i) = .True.
	do j = i+1,size(VarNames)
	  if(trim(VarNames(j))=='rhoBar_C') then
	    tested = .True.
	    mask(j) = .True.
	    if(VarValues(i) <= VarValues(j)) then ! if(rhoBar_M <= rhoBar_C) then
		  Run = .False.
		  write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with rhoBar values for model: ",ModNum, VarValues(j), VarValues(i)
		  call NonFatalError(ErrorMsg,ModNum)
		end if
	  end if
	end do

	if(.not. tested) then
	  if(VarValues(i) <= rhoBar(1)) then ! if(rhoBar_M <= rhoBar_C) then
		Run = .False.
		write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with rhoBar values for model: ",ModNum, rhoBar(1), VarValues(i)
		call NonFatalError(ErrorMsg,ModNum)
	  end if
	end if
  end if



! fFric < cFric?
  if(trim(VarNames(i)) == 'fFric' .and. (.not. mask(i))) then
    tested = .False.
    mask(i) = .True.
	do j = i+1,size(VarNames)
	  if(trim(VarNames(j)) == 'cFric') then
	    tested = .True.
	    mask(j) = .True.
	    if(VarValues(j) <= VarValues(i)) then ! if(cFric <= fFric) then
		  Run = .False.
		  write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with Fric values for model: ",ModNum, VarValues(i), VarValues(j)
		  call NonFatalError(ErrorMsg,ModNum)
		end if
	  end if
	end do

	if(.not. tested) then
	  if(cFric <= VarValues(i)) then ! if(cFric <= fFric) then
		Run = .False.
		write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with Fric values for model: ",ModNum, VarValues(i), cFric
		call NonFatalError(ErrorMsg,ModNum)
	  end if
	end if


  elseif(trim(VarNames(i)) == 'cFric' .and. (.not. mask(i))) then
    tested = .False.
    mask(i) = .True.
	do j = i+1,size(VarNames)
	  if(trim(VarNames(j)) == 'fFric') then
	    tested = .True.
	    mask(j) = .True.
	    if(VarValues(i) <= VarValues(j)) then ! if(cFric <= fFric) then
		  Run = .False.
		  write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with Fric values for model: ",ModNum, VarValues(j), VarValues(i)
		  call NonFatalError(ErrorMsg,ModNum)
		end if
	  end if
	end do

	if(.not. tested) then
	  if(VarValues(i) <= fFric) then ! if(cFric <= fFric) then
		Run = .False.
		write(ErrorMsg,"(A,I0,F10.3,F10.3)") "Model skipped due to problem with Fric values for model: ",ModNum, fFric, VarValues(i)
		call NonFatalError(ErrorMsg,ModNum)
	  end if
	end if
  end if

end do

end subroutine


end module