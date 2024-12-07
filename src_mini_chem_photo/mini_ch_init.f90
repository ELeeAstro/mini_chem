module mini_ch_init
  use mini_ch_precision
  use mini_ch_class
  implicit none

  public :: read_react_list, init_photochem
  !private ::

contains

  subroutine read_react_list(data_file, sp_file, net_dir, met)
    implicit none

    character(len=200), intent(in) :: data_file, sp_file, net_dir, met
    integer :: i, u, j, k, u2

    !! Read in the reaction list
    open(newunit=u,file=trim(data_file),status='old',action='read',form='formatted')

    do i = 1, 7
      read(u,*)
    end do

    read(u,*) n_reac

    allocate(re(n_reac))

    do i = 1, n_reac
      re(i)%id = i
      read(u,*)
      read(u,*) re(i)%re_t, re(i)%Tmin, re(i)%Tmax
      if (re(i)%re_t == 2) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%A, re(i)%B, re(i)%C
      else if (re(i)%re_t == 3) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%A0, re(i)%B0, re(i)%C0, re(i)%Ainf, re(i)%Binf, re(i)%Cinf
      else if (re(i)%re_t == 4) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%fname
      else if (re(i)%re_t == 5) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%br_idx
      end if

      allocate(re(i)%stoi_re(re(i)%n_re), re(i)%stoi_pr(re(i)%n_pr))
      read(u,*) (re(i)%stoi_re(j),j=1,re(i)%n_re), (re(i)%stoi_pr(k),k=1,re(i)%n_pr)

      allocate(re(i)%c_re(re(i)%n_re),re(i)%c_pr(re(i)%n_pr))
      read(u,*) (re(i)%c_re(j),j=1,re(i)%n_re), (re(i)%c_pr(k),k=1,re(i)%n_pr)

      ! Read reaction rate table if from table
      if (re(i)%re_t == 4) then
        print*, 'Reading: ', trim(met)//'_'//trim(re(i)%fname)
        open(newunit=u2,file=trim(net_dir)//trim(met)//'_'//trim(re(i)%fname),status='old',action='read',form='formatted')
        read(u2,*)
        read(u2,*) re(i)%nT, re(i)%nP, re(i)%nkf
        allocate(re(i)%T(re(i)%nT), re(i)%P(re(i)%nP), re(i)%kf(re(i)%nT,re(i)%nP))
        allocate(re(i)%lT(re(i)%nT), re(i)%lP(re(i)%nP), re(i)%lkf(re(i)%nT,re(i)%nP))
        read(u2,*)
        read(u2,*) (re(i)%T(j), j = 1, re(i)%nT)
        re(i)%lT(:) = log10(re(i)%T(:))
        read(u2,*)
        read(u2,*) (re(i)%P(k), k = 1, re(i)%nP)
        ! Convert P from bar to dyne
        re(i)%P(:) = re(i)%P(:) * 1.0e6_dp
        re(i)%lP(:) = log10(re(i)%P(:))
        read(u2,*)
        do j = 1,  re(i)%nT
          do k = 1, re(i)%nP
            read(u2,*) re(i)%kf(j,k)
            re(i)%kf(j,k) = max(re(i)%kf(j,k),1.0e-199_dp)
          end do
        end do
        re(i)%lkf(:,:) = log10(re(i)%kf(:,:))
        close(u2)
      end if

      ! print*, (re(i)%T(j), j = 1, re(i)%nT)
      ! print*, (re(i)%P(k), k = 1, re(i)%nP)

      ! print*, re(i)%id, re(i)%re_t, re(i)%Tmin, re(i)%Tmax
      ! if (re(i)%re_t == 2) then
      !   print*, re(i)%n_re, re(i)%n_pr, re(i)%A, re(i)%B, re(i)%C
      ! else if (re(i)%re_t == 3) then
      !   print*, re(i)%n_re, re(i)%n_pr, re(i)%A0, re(i)%B0, re(i)%C0, re(i)%Ainf, re(i)%Binf, re(i)%Cinf
      ! else if (re(i)%re_t == 4) then
      !   print*, re(i)%n_re, re(i)%n_pr, re(i)%fname
      ! end if
      ! print*, (re(i)%stoi_re(j),j=1,re(i)%n_re), (re(i)%stoi_pr(k),k=1,re(i)%n_pr)
      ! print*, (re(i)%c_re(j),j=1,re(i)%n_re), ' -> ', (re(i)%c_pr(k),k=1,re(i)%n_pr)

      ! We now find the delta mu.the difference in the number of products - reactants
      re(i)%dmu = real(sum(re(i)%stoi_pr(:)) - sum(re(i)%stoi_re(:)),dp)

      !print*, i, re(i)%dmu

    end do

    close(u)

    !! Get the species list from the species input file
    open(newunit=u,file=trim(sp_file),status='old',action='read',form='formatted')
    do i = 1, 6
      read(u,*)
    end do
    read(u,*) n_sp
    read(u,*)

    allocate(g_sp(n_sp))

    do i = 1, n_sp
      g_sp(i)%id = i
      read(u,*) g_sp(i)%c, g_sp(i)%mw, g_sp(i)%n_a
      allocate(g_sp(i)%a_l(g_sp(i)%n_a))
      read(u,*) g_sp(i)%a_l(:)
      allocate(g_sp(i)%a_h(g_sp(i)%n_a))
      read(u,*) g_sp(i)%a_h(:)
      ! print*, g_sp(i)%id, g_sp(i)%c, g_sp(i)%mw,  g_sp(i)%n_a
      ! print*, g_sp(i)%a_l(:)
      ! print*, g_sp(i)%a_h(:)
    end do

    close(u)

    !! For each reaction, we find the id number for each of the reacting species
    do i = 1, n_reac
      ! Allocate the gas index arrays
      allocate(re(i)%gi_re(re(i)%n_re), re(i)%gi_pr(re(i)%n_pr))
      ! Find the reactant index
      do k = 1, re(i)%n_re
        do j = 1, n_sp
          if (g_sp(j)%c == re(i)%c_re(k)) then
            re(i)%gi_re(k) = j
            exit
          end if
        end do
      end do

      ! Find the product index
      do k = 1, re(i)%n_pr
        do j = 1, n_sp
          if (g_sp(j)%c == re(i)%c_pr(k)) then
            re(i)%gi_pr(k) = j
            exit
          end if
        end do
      end do

      !print*, i, re(i)%c_re(:), re(i)%gi_re(:)
      !print*, i, re(i)%c_pr(:), re(i)%gi_pr(:)
    end do

  end subroutine read_react_list

  subroutine init_photochem(nlay, dbin1, dbin2, dbin_12trans, wl_s, wl_e, stellar_file)
    implicit none

    integer, intent(in) :: nlay
    real(dp), intent(in) :: dbin1, dbin2, dbin_12trans, wl_s, wl_e
    character(len=200), intent(in) :: stellar_file

    integer :: i, n1, n2, u, nlines, io
    integer :: idx, idx1
    real(dp) :: wl_now
    real(dp), allocatable, dimension(:) :: wl_f, flx_f

    !! First calculate wavelength grid for photochemistry calculations
    !! We follow the vulcan method, splitting the range into 2 parts, one at width dbin1 and dbin2
    !! that transitions at dbin_12trans - then output number of wavelengths and generate wavelength grid

    !! Count number of intervals in first length
    n1 = 1
    wl_now = wl_s
    do while (wl_now < dbin_12trans)
      wl_now = wl_now + dbin1
      n1 = n1 + 1
      !print*, n1, wl_now,  dbin_12trans
    end do

    !! Remove one in case of overshoot
    n1 = n1 - 1

    !! Now count number of intervals in second length
    n2 = 1
    wl_now = dbin_12trans
    do while (wl_now < wl_e)
      wl_now = wl_now + dbin2
      n2 = n2 + 1
      !print*, n2, wl_now, wl_e
    end do

    !! Remove one in case of overshoot
    n2 = n2 - 1

    !! Generate grid of wavelength according to total number of intervals
    nwl = n1 + n2
    allocate(wl_grid(nwl)) ! Global variable

    wl_grid(1) = wl_s
    do i = 2, nwl-1
      if (i <= n1) then 
        wl_grid(i) = wl_grid(i-1) + dbin1
      else 
        wl_grid(i) = wl_grid(i-1) + dbin2
      end if
    end do
    wl_grid(nwl) = wl_e

    !! Convert nm to cm for integration purposes
    wl_grid(:) = wl_grid(:)*1e-7_dp

    !! Allocate actinic flux array and stellar flux
    allocate(a_flux(nlay,nwl), s_flux(nwl))

    !! Read in stellar flux file and interpolate to wavelength grid
    !! Scale flux to distance and radius of planet

    open(newunit=u,file=trim(stellar_file),status='old',action='read',form='formatted')

    ! Read header
    read(u,*)

    ! Find number of lines in file
    nlines = 0
    do
      read(u,*,iostat=io)
      if (io /= 0) then 
        exit
      end if
      nlines = nlines + 1
    end do

    ! Allocate values for stellar flux file
    allocate(wl_f(nlines),flx_f(nlines)) 
    
    ! Rewind file
    rewind(u)
    ! Read header again
    read(u,*)

    ! Read file data
    do i = 1, nlines
      read(u,*) wl_f(i), flx_f(i)
      !print*, i, wl_f(i), flx_f(i)
    end do
    wl_f(:) = wl_f(:)*1e-7_dp

    close(u)


    ! Now we must interpolate the file to the decided wavelength grid
    do i = 1, nwl

      call locate(wl_f(:), nlines, wl_grid(i), idx)

      if (idx < 1) then
        s_flux(i) = 0.0_dp
      else if (idx >= nwl) then
        s_flux(i) = 0.0_dp
      else
        idx1 = idx + 1
        call linear_log_interp(wl_grid(i), wl_f(idx), wl_f(idx1), flx_f(idx), flx_f(idx1), s_flux(i))
      end if

      !print*, i, wl_grid(i)*1e7_dp, s_flux(i)
    end do


    !! Now we must read in photochemical cross sections for each species and interpolate as
    !! well as calculate Rayleigh scattering for each species

  end subroutine init_photochem

end module mini_ch_init
