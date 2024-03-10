module mini_ch_ce_interp
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  logical :: first_call = .True.

  ! T, P and vmr data arrays from table
  integer :: nT, nP, npoint, nmol
  character(len=10), allocatable, dimension(:) :: m_name
  real(dp), allocatable, dimension(:) :: T_t, P_t, lP_t, lT_t
  real(dp), allocatable, dimension(:,:) :: mu_t
  real(dp), allocatable, dimension(:,:,:) :: vmr_t

  public :: interp_ce_table
  private :: first_call, read_ce_table, Bezier_interp, locate

  contains

  subroutine interp_ce_table(n_sp, T_in, P_in, VMR, mu, table)
    implicit none

    integer, intent(in) :: n_sp
    real(dp), intent(in) :: T_in, P_in
    character(len=200), intent(in) :: table

    real(dp), intent(out) :: mu
    real(dp), dimension(n_sp), intent(out) :: VMR

    integer :: m, i_p1, i_p2, i_p3, i_t1, i_t2, i_t3
    real(dp) :: lP_in, lT_in, y_out
    real(dp), dimension(3) :: lTa, lPa, mua, xa, xa_out, mua_out

    if (first_call .eqv. .True.) then
      call read_ce_table(table)
      first_call = .False.
    end if

    ! log input T and P
    lP_in = log10(P_in/1e5_dp)
    lT_in = log10(T_in)

    ! Find upper and lower T and P triplet indexes
    call locate(P_t, nP, P_in, i_p2)
    i_p1 = i_p2 - 1
    i_p3 = i_p2 + 1

    if (i_p1 <= 0) then
      i_p1 = 1
      i_p2 = 2
      i_p3 = 3
    else if (i_p3 > nP) then
      i_p1 = nP - 2
      i_p2 = nP - 1
      i_p3 = nP
    end if

    ! Check if input temperature is within table range
    if (T_in <= T_t(1)) then

      ! Perform Bezier interpolation at minimum table temperature
      lPa(1) = lP_t(i_p1)
      lPa(2) = lP_t(i_p2)
      lPa(3) = lP_t(i_p3)
      mua(1) = mu_t(i_p1,1)
      mua(2) = mu_t(i_p2,1)
      mua(3) = mu_t(i_p3,1)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mu)

      do m = 1, n_sp
        xa(1) = vmr_t(m,i_p1,1)
        xa(2) = vmr_t(m,i_p2,1)
        xa(3) = vmr_t(m,i_p3,1)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, y_out)
        VMR(m) = 10.0_dp**y_out
      end do

    else if (T_in >= T_t(nT)) then

      ! Perform Bezier interpolation at maximum table temperature
      lPa(1) = lP_t(i_p1)
      lPa(2) = lP_t(i_p2)
      lPa(3) = lP_t(i_p3)
      mua(1) = mu_t(i_p1,nT)
      mua(2) = mu_t(i_p2,nT)
      mua(3) = mu_t(i_p3,nT)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mu)

      do m = 1, n_sp
        xa(1) = vmr_t(m,i_p1,nT)
        xa(2) = vmr_t(m,i_p2,nT)
        xa(3) = vmr_t(m,i_p3,nT)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, y_out)
        VMR(m) = 10.0_dp**y_out
      end do

    else

      ! Perform 2D Bezier interpolation by performing interpolation 4 times

      ! Find temperature index triplet
      call locate(T_t, nT, T_in, i_t2)
      i_t1 = i_t2 - 1
      i_t3 = i_t2 + 1

      lPa(1) = lP_t(i_p1)
      lPa(2) = lP_t(i_p2)
      lPa(3) = lP_t(i_p3)

      mua(1) = mu_t(i_p1,i_t1)
      mua(2) = mu_t(i_p2,i_t1)
      mua(3) = mu_t(i_p3,i_t1)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mua_out(1)) ! Result at T1, P_in
      mua(1) = mu_t(i_p1,i_t2)
      mua(2) = mu_t(i_p2,i_t2)
      mua(3) = mu_t(i_p3,i_t2)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mua_out(2)) ! Result at T2, P_in
      mua(1) = mu_t(i_p1,i_t3)
      mua(2) = mu_t(i_p2,i_t3)
      mua(3) = mu_t(i_p3,i_t3)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mua_out(3)) ! Result at T3, P_in
      lTa(1) = lT_t(i_t1)
      lTa(2) = lT_t(i_t2)
      lTa(3) = lT_t(i_t3)
      call Bezier_interp(lTa(:), mua_out(:), 3, lT_in, mu) ! Result at T_in, P_in

      do m = 1, n_sp
        xa(1) = vmr_t(m,i_p1,i_t1)
        xa(2) = vmr_t(m,i_p2,i_t1)
        xa(3) = vmr_t(m,i_p3,i_t1)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, xa_out(1)) ! Result at T1, P_in
        xa(1) = vmr_t(m,i_p1,i_t2)
        xa(2) = vmr_t(m,i_p2,i_t2)
        xa(3) = vmr_t(m,i_p3,i_t2)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, xa_out(2)) ! Result at T2, P_in
        xa(1) = vmr_t(m,i_p1,i_t3)
        xa(2) = vmr_t(m,i_p2,i_t3)
        xa(3) = vmr_t(m,i_p3,i_t3)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, xa_out(3)) ! Result at T3, P_in
        lTa(1) = lT_t(i_t1)
        lTa(2) = lT_t(i_t2)
        lTa(3) = lT_t(i_t3)
        call Bezier_interp(lTa(:), xa_out(:), 3, lT_in, y_out) ! Result at T_in, P_in
        VMR(m) = 10.0_dp**y_out
      end do

    end if

  end subroutine interp_ce_table

  subroutine read_ce_table(table)
    implicit none

    character(len=200), intent(in) :: table

    integer :: i, j, u

    !! Read T and P grid from file + mu and VMR data

    open(newunit=u, file=trim(table),action='read',form='formatted',status='old')

    read(u,*) nT, nP, npoint, nmol

    !print*, np, nmol

    allocate(m_name(nmol))
    read(u,*) (m_name(i),i=1,nmol)

    allocate(T_t(nT))
    read(u,*) (T_t(i),i=1,nT)

    allocate(P_t(nP))
    read(u,*) (P_t(i),i=1,nP)

    allocate(mu_t(nP,nT), vmr_t(nmol,nP,nT))
    do i = 1, nT
      do j = 1, nP
        read(u,*) mu_t(j,i), vmr_t(1:nmol,j,i)
        vmr_t(:,j,i) = log10(max(vmr_t(:,j,i),1e-99_dp))
      end do
    end do

    ! Log10 arrays of T-p grid
    allocate(lT_t(nT),lP_t(nP))
    lT_t = log10(T_t)
    lP_t = log10(P_t)

    close(u)

  end subroutine read_ce_table

  
  ! Perform Bezier interpolation
  subroutine Bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, wh, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if ((x > xi(1)) .and. (x < xi(2))) then
      ! left hand side interpolation
      !print*,'left'
      wh = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if ((wh <= min(wlim,wlim1)) .or. (wh >= max(wlim,wlim1))) then
        wh = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (wh*dy/dx + (1.0_dp - wh)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      wh = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if ((wh <= min(wlim,wlim1)) .or. (wh >= max(wlim,wlim1))) then
        wh = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (wh*dy1/dx1 + (1.0_dp - wh)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine Bezier_interp


  subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real(dp), dimension(n), intent(in) :: arr
    real(dp), intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = n+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

end module mini_ch_ce_interp