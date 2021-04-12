module mini_ch_read_reac_list
  use mini_ch_precision
  use mini_ch_class
  implicit none

contains

  subroutine read_react_list()
    implicit none

    integer :: i, u, j, k, u2

    !! Read in the reaction list
    open(newunit=u,file='mini_chem_data_HO.txt',status='old',action='read',form='formatted')

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
      end if

      allocate(re(i)%stoi_re(re(i)%n_re), re(i)%stoi_pr(re(i)%n_pr))
      read(u,*) (re(i)%stoi_re(j),j=1,re(i)%n_re), (re(i)%stoi_pr(k),k=1,re(i)%n_pr)

      allocate(re(i)%c_re(re(i)%n_re),re(i)%c_pr(re(i)%n_pr))
      read(u,*) (re(i)%c_re(j),j=1,re(i)%n_re), (re(i)%c_pr(k),k=1,re(i)%n_pr)

      ! Read reaction rate table if from table
      if (re(i)%re_t == 4) then
        print*, 'Reading: ', trim(re(i)%fname)
        open(newunit=u2,file=trim(re(i)%fname),status='old',action='read',form='formatted')
        read(u2,*)
        read(u2,*) re(i)%nT, re(i)%nP, re(i)%nkf
        allocate(re(i)%T(re(i)%nT), re(i)%P(re(i)%nP), re(i)%kf(re(i)%nT,re(i)%nP))
        read(u2,*)
        read(u2,*) (re(i)%T(j), j = 1, re(i)%nT)
        read(u2,*)
        read(u2,*) (re(i)%P(k), k = 1, re(i)%nP)
        ! Convert P from bar to dyne
        re(i)%P(:) = re(i)%P(:) * 1.0e6_dp
        read(u2,*)
        do j = 1,  re(i)%nT
          do k = 1, re(i)%nP
            read(u2,*) re(i)%kf(j,k)
          end do
        end do
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
    open(newunit=u,file='mini_chem_sp_HO.txt',status='old',action='read',form='formatted')
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

      ! print*, i, re(i)%c_re(:), re(i)%gi_re(:)
      ! print*, i, re(i)%c_pr(:), re(i)%gi_pr(:)
    end do

  end subroutine read_react_list


end module mini_ch_read_reac_list
