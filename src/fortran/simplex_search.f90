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

module simplex_search

! Module containing simplex search optimization routine

  implicit none

! Options type for direct searches

  type ds_options_type
    double precision :: tol       ! tolerance in simplex radius before
                                  !   triggering a stop condition
    integer :: maxit              ! Max steps allowed before stopping
  end type ds_options_type

  contains

!=============================================================================80
!
! Nelder-Mead simplex search algorithm
!
!=============================================================================80
subroutine simplexsearch(xopt, fmin, step, fevals, objfunc, x0, given_f0_ref,  &
                         f0_ref, ds_options)

  !! xopt          out: designvars result 
  !! fmin          out: smallest value of objective function     
  !! step          out: iteration steps needed
  !! fevals        out: number of evaluation of objective function 
  !! objfunc       interface objective function 
  !! x0            start values of designvars 
  !! given_f0_ref  is there a reference reference start value of objective function 
  !! f0_ref        inout: reference start value of objective function 

  use optimization_util, only : bubble_sort, design_radius, write_design

  double precision, dimension(:), intent(inout) :: xopt
  double precision, intent(out) :: fmin
  integer, intent(out) :: step, fevals

  interface
    double precision function objfunc(x, evaluate_only_geometry)
      double precision, dimension(:), intent(in) :: x
      logical, intent(in), optional :: evaluate_only_geometry
    end function
  end interface

  double precision, dimension(:), intent(in) :: x0
  double precision, intent(inout) :: f0_ref
  logical, intent(in) :: given_f0_ref
  type (ds_options_type), intent(in) :: ds_options

  double precision, dimension(size(x0,1),size(x0,1)+1) :: dv
  double precision, dimension(size(x0,1)+1) :: objvals
  double precision, dimension(size(x0,1)) :: xcen, xr, xe, xc

  double precision :: rho, xi, gam, sigma, fr, fe, fc, f0, mincurr, radius
  integer :: i, j, nvars, designcounter
  logical :: converged, needshrink, signal_progress

! Standard Nelder-Mead constants

  rho = 1.d0
  xi = 2.d0
  gam = 0.5d0
  sigma = 0.5d0

! Set up or read initialzation data

  nvars = size(x0,1)


!   Get f0 (reference seed design objective function)

  if (given_f0_ref) then
    f0 = f0_ref
  else 
    f0 = objfunc(x0)
    f0_ref = f0
  end if

!   Set up initial simplex

  fevals = 0
  do j = 1, nvars
    do i = 1, nvars
      if (i == j) then
        if (x0(i) == 0.d0) then
          dv(i,j) = 0.00025d0
        else
          if (x0(i) > 0) then 
            dv(i,j) = x0(i) + 0.05d0 !0.1d0 ! 1.05d0*x0(i)
          else
            dv(i,j) = x0(i) - 0.05d0 ! 0.1d0 ! 1.05d0*x0(i)
          end if 
        end if
      else
        dv(i,j) = x0(i)
      end if
    end do
    objvals(j) = objfunc(dv(:,j))
    fevals = fevals + 1
  end do

  dv(:,nvars+1) = x0
  objvals(nvars+1) = objfunc(x0)
  fevals = fevals + 1

! Counters

  step = 0
  designcounter = 0

! Initial minimum value

  fmin = minval(objvals)
  mincurr = fmin

! Iterative procedure for optimization
 
  needshrink = .false.
  converged = .false.

  main_loop: do while (.not. converged)

    step = step + 1
    if (step == ds_options%maxit) converged = .true.
    
!   Sort according to ascending objective function value

    call bubble_sort(dv, objvals)
    mincurr = objvals(1)

!   Update fmin if appropriate

    if (mincurr < fmin) then
      fmin = mincurr
      signal_progress = .true.
    else
      signal_progress = .false.
    end if

!   Check for convergence

    radius = design_radius(dv)
    if (radius < ds_options%tol) converged = .true.

!   Compute the centroid of the best nvals designs

    xcen(:) = 0.d0
    do i = 1, nvars
      xcen = xcen + dv(:,i)
    end do
    xcen = xcen/dble(nvars)

!   Compute the reflection point and evaluate its objective function value

    xr = (1.d0 + rho)*xcen - rho*dv(:,nvars+1)
    fr = objfunc(xr)
    fevals = fevals + 1

    expand_or_contract: if (objvals(1) <= fr .and. fr < objvals(nvars)) then

!      Accept reflection point

       dv(:,nvars+1) = xr
       objvals(nvars+1) = fr
       cycle

    elseif (fr < objvals(1)) then

!     Expand

      xe = (1.d0 + rho*xi)*xcen - rho*xi*dv(:,nvars+1)
      fe = objfunc(xe)
      fevals = fevals + 1
      if (fe < fr) then
        dv(:,nvars+1) = xe
        objvals(nvars+1) = fe
      else
        dv(:,nvars+1) = xr
        objvals(nvars+1) = fr
      end if
      cycle

    elseif (fr >= objvals(nvars)) then

!     Outside contraction

        contraction: if (fr < objvals(nvars+1)) then

          xc = (1.d0 + rho*gam)*xcen - rho*gam*dv(:,nvars+1)
          fc = objfunc(xc)
          fevals = fevals + 1

          if (fc < fr) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            needshrink = .false.
          else
            needshrink = .true.
          end if

!       Inside contraction

        else 

          xc = (1.d0 - gam)*xcen + gam*dv(:,nvars+1)
          fc = objfunc(xc)
          fevals = fevals + 1
          
          if (fc < objvals(nvars+1) ) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            needshrink = .false.
          else
            needshrink = .true.
          end if

        end if contraction

!       Shrink

        shrink: if (needshrink) then

          do i = 2, nvars + 1
            dv(:,i) = dv(:,1) + sigma*(dv(:,i) - dv(:,1))
            objvals(i) = objfunc(dv(:,i))
            fevals = fevals + 1
          end do
          cycle

        else

          cycle

        end if shrink

    end if expand_or_contract

  end do main_loop

! Sort one more time according to ascending objective function value

  call bubble_sort(dv, objvals)
  xopt = dv(:,1)
  fmin = objvals(1)

! Check for convergence one more time

  radius = design_radius(dv)

! Display warning if max iterations are reached
  
  if (step == ds_options%maxit .and. (radius >= ds_options%tol)) then
    write(*,*) 'Warning: Simplex optimizer forced to exit due to the max number'
    write(*,*) '         of iterations being reached.'
  end if

end subroutine simplexsearch


end module simplex_search
