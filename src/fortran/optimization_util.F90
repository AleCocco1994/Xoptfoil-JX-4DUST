!  This file is part of XOPTFOIL.

!  XOPTFOIL is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  XOPTFOIL is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with XOPTFOIL.  If not, see <http://www.gnu.org/licenses/>.

!  Copyright (C) 2017-2019 Daniel Prosser

module optimization_util

  use os_util

! Module containing optimization routines

  implicit none

  contains


  subroutine init_random_seed()

  !! Initializes a random seed (subroutine from gcc.gnu.org)

  ! For ifort compatibility
#ifdef intel_compilers
  use ifport, only : getpid  
#endif

  integer, dimension(:), allocatable :: myseed
  integer :: i, n, un, istat, dt(8), pid, t(2), s
  integer(8) :: count, tms
  
  call random_seed(size = n)
  allocate(myseed(n))

  ! First try if the OS provides a random number generator

  un = 18
  open(newunit=un, file="/dev/urandom", access="stream",                       &
       form="unformatted", action="read", status="old", iostat=istat)

  if (istat == 0) then

     read(un) myseed
     close(un)

  else

     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.

     call system_clock(count)
     if (count /= 0) then
        t = transfer(count, t)
     else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = transfer(tms, t)
     end if

     s = ieor(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = ieor(s, pid)
     if (n >= 3) then
        myseed(1) = t(1) + 36269
        myseed(2) = t(2) + 72551
        myseed(3) = pid
        if (n > 3) then
           myseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
     else
        myseed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     end if

  end if

  call random_seed(put=myseed)
  deallocate(myseed)

end subroutine init_random_seed



subroutine initial_designs(dv, objval, fevals, objfunc, xmin, xmax, &
                           initial_x0_based, x0, &
                           feasible_init, feasible_limit, attempts)

 !! Creates initial designs and tries to make them feasible, if desired
 !    With 'initial_x0_based' the designs are close to x0, which is good for Bezier
 !    whereas  Hicks-Henne and Thickness-Camber get designs within xmin and xmax                        
 

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval
  integer, intent(out) :: fevals
  double precision, dimension(:), intent(in) :: xmin, xmax, x0
  logical, intent(in) :: initial_x0_based, feasible_init
  double precision, intent(in) :: feasible_limit
  integer, intent(in) :: attempts

  interface
    double precision function objfunc(x, evaluate_only_geometry)
      double precision, dimension(:), intent(in) :: x
      logical, intent(in), optional :: evaluate_only_geometry
    end function
  end interface

  integer :: i, j, pop, nvars, initcount
  double precision, dimension(:), allocatable :: randvec1, designstore
  double precision, dimension(:,:), allocatable :: dv_delta
  double precision :: minstore

  ! Initial settings and memory allocation

  nvars = size(dv,1)
  pop = size(dv,2)
  allocate(randvec1(nvars))
  allocate(designstore(nvars))

  fevals = pop

! Set up matrix of random numbers for initial designs
  
  call random_number(dv)

  dv_delta = (dv - 0.5d0) * 0.05d0          ! deviation factor from x0  

! Initial population of designs set between xmin and xmax

  write(*,'(/," - ",A)', advance ='no') 'Generating and evaluating '//stri(pop)//' initial designs: '

  ! take x0 as initial for the first particle
  dv(:,1) = x0
  objval(1) = objfunc(x0, .true.)         ! should always be 1.0 ...


!$OMP parallel do private(j)
  do i = 2, pop

    if (initial_x0_based) then
      ! init values will be random delta to x0        
      do j = 1, nvars
        dv(j,i) = x0(j) + dv_delta(j,i) * (xmax(j)-xmin(j))
        dv(j,i) = max (dv(j,i), xmin(j))
        dv(j,i) = min (dv(j,i), xmax(j))
      end do 
    else 
      ! freestyle init values between min and max       
      dv(:,i) = maxval(xmax - xmin)*dv(:,i) + xmin
    end if 

    objval(i) = objfunc(dv(:,i), .true.)  ! evaluate only geometry
    
  end do
!$omp end parallel do

  call show_design_info (feasible_limit, .true., objval)
  write(*,*) 


  ! Enforce initially feasible designs

  if (feasible_init) then

    write(*,'(" - ",A)', advance ='no') 'Trying to improve with max '//stri(attempts)//' init attempts: '

!$OMP parallel do private(j, initcount, minstore, designstore, randvec1)
    do i = 1, pop

      initcount = 0
      minstore = objval(i)
      designstore = dv(:,i)

      ! Take a number of tries to fix infeasible designs

      do while ((initcount <= attempts) .and.                                  &
               (objval(i) >= feasible_limit))

        call random_number(randvec1)

        if (initial_x0_based) then
          ! init values will be random delta to x0 
          randvec1 = (randvec1 - 0.5d0) * 0.2d0        
          do j = 1, nvars
            dv(j,i) = x0(j) + randvec1 (j) * (xmax(j)-xmin(j))
            dv(j,i) = max (dv(j,i), xmin(j))
            dv(j,i) = min (dv(j,i), xmax(j))
          end do 
        else 
          ! freestyle init values between min and max       
          dv(:,i) = maxval(xmax - xmin)*randvec1 + xmin
        end if 
    
        objval(i) = objfunc(dv(:,i), .true.)  ! evaluate only geometry
    
        if (objval(i) < minstore) then
          minstore = objval(i)
          designstore = dv(:,i)
        end if
        initcount = initcount + 1
!$omp critical
        fevals = fevals + 1
!$omp end critical
      end do

      ! Pick the best design tested if a feasible design was not found

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        dv(:,i) = designstore
        objval(i) = minstore
      end if

    end do

!$omp end parallel do

    call show_design_info (feasible_limit, .false., objval)
    write(*,*) 
  
  end if

end subroutine initial_designs




subroutine  show_design_info (feasible_limit, first_run, objval)
  !! Shows user info about result of inital design 

  double precision, intent(in)  :: feasible_limit
  logical                       :: first_run
  double precision, dimension (:), intent(in)  :: objval
  integer       :: color, i
  Character (1) :: sign 

  call print_colored (COLOR_NOTE, ' ')

  do i = 1, size(objval)
    if (objval(i) == 1d0) then 
      color = COLOR_NOTE                          
      ! color = COLOR_GOOD                            
      sign  = '+'
    elseif (objval(i) <= feasible_limit) then 
      color = COLOR_NOTE                           
      sign  = 'o'
    else if (objval(i) > feasible_limit) then   
      if (first_run) then 
        color = COLOR_NOTE                         
      else
        color = COLOR_ERROR                           
      end if  
      sign  = 'x'
    else  
      color = COLOR_NOTE                         
      sign  = '-'
    end if 

    call print_colored (color, sign)            
  end do 
  
end subroutine show_design_info






!=============================================================================80
!
! Computes max radius of designs (used for evaluating convergence
!
!=============================================================================80
function design_radius(dv)

  use math_deps, only : norm_2

  double precision, dimension(:,:), intent(in) :: dv
  double precision design_radius

  integer :: i, ndesigns
  double precision, dimension(size(dv,1)) :: design_centroid
  double precision :: radius

  ! Compute centroid of designs

  ndesigns = size(dv,2)
  design_centroid(:) = 0.d0
  do i = 1, ndesigns
    design_centroid = design_centroid + dv(:,i)
  end do
  design_centroid = design_centroid / dble(ndesigns)

  ! Compute max design radius

  design_radius = 0.d0
  do i = 1, ndesigns
    radius = norm_2(dv(:,i) - design_centroid)
    if (radius > design_radius) design_radius = radius
  end do

end function

!=============================================================================80
!
! Sorts a set of designs according to their objective function value
!
!=============================================================================80
subroutine bubble_sort(dv, objvals)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals

  double precision, dimension(size(dv,1),size(dv,2)) :: tempdv
  double precision, dimension(size(dv,2)) :: tempvals
  integer, dimension(size(dv,2)) :: finalorder, temporder
  integer :: nvars, ndesigns, i, sortcounter
  logical :: sorted

  nvars = size(dv,1)
  ndesigns = size(dv,2)

! Set up indexing array

  do i = 1, ndesigns
    finalorder(i) = i
  end do
  temporder = finalorder

! Bubble sorting algorithm

  sorted = .false.
  tempvals = objvals
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, ndesigns - 1
      if (objvals(i+1) < objvals(i)) then

!       Flip the order of these elements. temp arrays are to preserve values.

        tempvals(i) = objvals(i+1)
        tempvals(i+1) = objvals(i)
        temporder(i) = finalorder(i+1)
        temporder(i+1) = finalorder(i)
        finalorder(i) = temporder(i)
        finalorder(i+1) = temporder(i+1)
        objvals(i) = tempvals(i)
        objvals(i+1) = tempvals(i+1)
        sortcounter = sortcounter + 1

      end if
    end do
    if (sortcounter == 0) sorted = .true.
    
  end do

! Use indexing array to rearrange order of designs

  do i = 1, ndesigns
    tempdv(:,i) = dv(:,finalorder(i))
  end do
  dv = tempdv

end subroutine bubble_sort

!=============================================================================80
!
! Pops item out of a vetor.  Note: doesn't actually change size of vector, just
! shuffles data so that the first nitems-1 entries represent the new vector.
!
!=============================================================================80
subroutine pop_double_vector(vector, nitems, popidx)

  double precision, dimension(:), intent(inout) :: vector
  integer, intent(in) :: nitems, popidx

  integer :: i
  double precision, dimension(size(vector,1)) :: tempvector

! Copy input vector

  tempvector = vector

! Populate output vector

  do i = 1, nitems-1
    if (i < popidx) then
      vector(i) = tempvector(i)
    else
      vector(i) = tempvector(i+1)
    end if
  end do

end subroutine pop_double_vector

!=============================================================================80
!
! Pops item out of a vetor.  Note: doesn't actually change size of vector, just
! shuffles data so that the first nitems-1 entries represent the new vector.
!
!=============================================================================80
subroutine pop_integer_vector(vector, nitems, popidx)

  integer, dimension(:), intent(inout) :: vector
  integer, intent(in) :: nitems, popidx

  integer :: i
  integer, dimension(size(vector,1)) :: tempvector

! Copy input vector

  tempvector = vector

! Populate output vector

  do i = 1, nitems-1
    if (i < popidx) then
      vector(i) = tempvector(i)
    else
      vector(i) = tempvector(i+1)
    end if
  end do

end subroutine pop_integer_vector

!=============================================================================80
!
! Writes design variables to file
!
!=============================================================================80
subroutine write_design(filename, filestat, variables, counter)

  character(*), intent(in) :: filename, filestat
  double precision, dimension(:), intent(in) :: variables
  integer, intent(in) :: counter

  integer, save :: iunit
  integer :: nvars, i
  character(30) :: text

  nvars = size(variables,1)
  iunit = 17

! Open the file and write to it if requested

  if (trim(filestat) == 'new') then
    open(unit=iunit, file=filename, status='replace')
    write(text,*) nvars
    text = adjustl(text)
    write(iunit,'(A)') 'Number of variables: '//trim(text)
  else
    open(unit=iunit, file=filename, status='old', position='append')
  end if

! Write iteration number and the design variables to file

  write(text,*) counter
  text = adjustl(text)
  write(iunit,'(A)') 'Design number '//trim(text)
  do i = 1, nvars
    write(iunit,'(es25.16)') variables(i)
  end do

! Close the file 

  close(iunit)

end subroutine write_design

!=============================================================================80
!
! Reads commands from run_control file and clears any unrecognized commands
!
!=============================================================================80
subroutine read_run_control(commands, ncommands)

  character(80), dimension(:), intent(inout) :: commands
  integer, intent(out) :: ncommands

  character(80) :: buffer
  integer :: rcunit, ioerr
  
  commands(:) = ""
  ncommands = 0

  rcunit = 18
  open(unit=rcunit, file='run_control', status='old', iostat=ioerr, err=501)
  if (ioerr /= 0) then
    return
  else
    do while (1 .eq. 1)
      read(rcunit,'(A)',end=501) buffer
      if ( (trim(buffer) == "stop") .or.                                       &
           (trim(buffer) == "stop_monitoring") ) then
        commands(ncommands+1) = buffer
        ncommands = ncommands + 1
      else
        write(*,*) 'Warning: unrecognized command "'//trim(buffer)//           &
                   '" in run_control.' 
      end if
    end do
  end if
   
  write(*,*) "Warning: error encountered while reading run_control. Skipping."
  return
501 close(rcunit)

!   open(unit=rcunit, file='run_control', status='replace', err=502)
!   do i = 1, ncommands
!     write(rcunit,'(A)') commands(ncommands)
!   end do
!   close(rcunit) 

!   return

! 502 write(*,*) "Warning: error encountered while reading run_control. Skipping."

end subroutine read_run_control



subroutine write_history_header (iunit)
  !! write csv header of op points data 

  integer, intent(in) :: iunit
  write(iunit,'(A)') 'Iteration; Design; Objective function; % Improvement; Design radius'
  flush(iunit)

end subroutine


subroutine  write_history (iunit, step, new_design, designcounter, radius, fmin, f0)
  !! write iteration result to history file 'iunit' during optimization

  integer, intent(in)           :: iunit, step, designcounter 
  logical, intent(in)           :: new_design
  double precision, intent(in)  :: radius ,fmin, f0 
  double precision :: relfmin

  relfmin = (f0 - fmin)/f0 * 100.d0
  if (new_design) then 
    write(iunit,'(I6,";",I6, 3(";", F11.7))') step, designcounter, fmin, relfmin, radius
  else 
    write(iunit,'(I6,";",A6, 3(";", F11.7))') step, ''           , fmin, relfmin, radius
  end if 
  flush(iunit)

end subroutine  write_history



end module optimization_util
