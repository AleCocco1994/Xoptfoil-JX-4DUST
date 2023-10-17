!===============================================================================
!
! Runs Xfoil at a specified op_point which is either
!  - at an angle of attack
!  - at an lift coefficient
!
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
!
!===============================================================================

subroutine run_op_point (op_point_spec,        &
                          viscous_mode, maxit, show_details,  &
                          op_point_result)

  use xfoil_inc

  type(op_point_specification_type), intent(in)  :: op_point_spec
  logical,                           intent(in)  :: viscous_mode, show_details
  integer,                           intent(in)  :: maxit
  type(op_point_result_type),        intent(out) :: op_point_result

  integer         :: niter_needed
  doubleprecision :: save_ACRIT


  op_point_result%cl    = 0.d0
  op_point_result%cd    = 0.d0
  op_point_result%alpha = 0.d0
  op_point_result%cm    = 0.d0
  op_point_result%xtrt  = 0.d0
  op_point_result%xtrb  = 0.d0
  op_point_result%converged = .true.

! Support Type 1 and 2 re numbers  
  REINF1 = op_point_spec%re%number
  RETYP  = op_point_spec%re%type 
  MATYP  = op_point_spec%ma%type 
  call MINFSET(op_point_spec%ma%number)

! Set compressibility parameters from MINF
  CALL COMSET

! Set ncrit per point
  save_ACRIT = ACRIT
  if (op_point_spec%ncrit /= -1d0) ACRIT = op_point_spec%ncrit

! Inviscid calculations for specified cl or alpha
  if (op_point_spec%spec_cl) then
    LALFA = .FALSE.
    ALFA = 0.d0
    CLSPEC = op_point_spec%value
    call SPECCL
  else 
    LALFA = .TRUE.
    ALFA = op_point_spec%value * DTOR
    call SPECAL
  end if

  if (abs(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0E-5) LVCONV = .false.

  ! Viscous calculations (if requested)

  op_point_result%converged = .true. 

  if (viscous_mode) then 
    
    call VISCAL(maxit, niter_needed)

    ! coverged? 

    if (niter_needed > maxit) then 
      op_point_result%converged = .false.
    ! RMSBL equals to viscrms() formerly used...
    else if (.not. LVCONV .or. (RMSBL > 1.D-4)) then 
      op_point_result%converged = .false.
    else 
      op_point_result%converged = .true.
    end if

  end if

  ! Restore default ncrit 

  ACRIT = save_ACRIT

  ! Outputs

  op_point_result%cl    = CL
  op_point_result%cm    = CM
  op_point_result%alpha = ALFA/DTOR

  if (viscous_mode) then 
    op_point_result%cd   = CD 
    op_point_result %xtrt = XOCTR(1)
    op_point_result%xtrb = XOCTR(2)
    if (op_point_result%converged) then 
      ! call detect_bubble (op_point_result%bubblet, op_point_result%bubbleb)
      op_point_result%bubblet%found = .false.
      op_point_result%bubbleb%found = .false.
    end if
  else
    op_point_result%cd   = CDP
    op_point_result%xtrt = 0.d0
    op_point_result%xtrb = 0.d0
    op_point_result%bubblet%found = .false.
    op_point_result%bubbleb%found = .false.
  end if

! Final check for NaNs

  if (isnan(op_point_result%cl)) then
    op_point_result%cl = -1.D+08
    op_point_result%converged = .false.
  end if
  if (isnan(op_point_result%cd)) then
    op_point_result%cd = 1.D+08
    op_point_result%converged = .false.
  end if
  if (isnan(op_point_result%cm)) then
    op_point_result%cm = -1.D+08
    op_point_result%converged = .false.
  end if

  if(show_details) then 
    if (op_point_result%converged) then
      !call print_colored (COLOR_NORMAL,  '.' // stri(niter_needed))
      call print_colored (COLOR_NORMAL,  '.')
    else
      !call print_colored (COLOR_WARNING, 'x' // stri(niter_needed))
      call print_colored (COLOR_WARNING, 'x')
    end if
  end if

! when xfoil doesn't converge, it tends to write invalid values in the original field ...
  if (.not. op_point_result%converged) then
    if (op_point_spec%spec_cl) then
      op_point_result%cl    = op_point_spec%value
    else 
      op_point_result%alpha = op_point_spec%value 
    end if
  end if 


end subroutine run_op_point