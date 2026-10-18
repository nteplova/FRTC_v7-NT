!=======================================================================
! Модуль FP2D_bridge_mod
! Связка между 1D кодом FRTC и 2D решателем FP2D_solver_mod.
!
! Режимы (fp_solver):
!   2 — 2D нерелятивистский
!   3 — 2D релятивистский
!
! Отображение k ↔ знак p_par:
!   В 1D-коде (см. 3_FokkerPlanck.f90, fokkerplanck_compute):
!     k=1 → dij(:,:,1) (положительный спектр, ispectr=+1)
!           znak = -1 (обратное электрическое поле в 1D-уравнении)
!     k=2 → dij(:,:,2) (отрицательный спектр, ispectr=-1)
!           znak = +1
!   В 2D-коде мы работаем с полным распределением на p_par ∈ [-p_max,+p_max],
!   поэтому:
!     k=1 → D_QL прикладывается к p_par < 0
!     k=2 → D_QL прикладывается к p_par > 0
!
! ВАЖНО: fij0 в этом модуле — это маргинал НАЧАЛЬНОГО максвеллиана
! (без E-поля и без волн). Это отличается от 1D-кода, где fij0 — решение
! FP-уравнения только с E-полем. См. подпрограмму fp2d_extract_1d.
!=======================================================================
module FP2D_bridge_mod
    use kind_module,     only: wp
    use FP2D_solver_mod
    use nr_grid,         only: MAX_NR, vij, fij, fij0, dfij, dij
    use maxwell,         only: i0
    use plasma,          only: fvt, ft, zefff, enorm, fst
    use rt_parameters,   only: nr
    implicit none

    real(wp), allocatable, save :: f2D(:,:,:,:)      ! (nx, ny, nr, 2)
    real(wp), allocatable, save :: D_QL2D(:,:,:,:)   ! (nx, ny, nr, 2)
    real(wp), save :: p_par_grid(nx), p_perp_grid(ny)
    logical,  save :: fp2d_initialized = .FALSE.
    logical,  save :: fp2d_rel_saved   = .FALSE.

    real(wp), parameter :: c_light = 3.0e10_wp
    real(wp), parameter :: m_e_c2  = 8.187e-7_wp

    ! ---- Множитель нормировки маргинала ----
    ! Если f2D содержит (2π)^{-3/2} exp(-p²/2), то marginal = ∫ f 2π p_⊥ dp_⊥
    ! Если f2D содержит exp(-p²/2) без префактора, то нормировку нужно
    ! задавать так, чтобы результат совпадал с 1D-кодом:
    !   f_1D(v) = exp(-v²/2) / sqrt(2π)
    real(wp), parameter :: norm_2pi32 = 15.7496099457224_wp   ! (2π)^{3/2}

    ! ---- Масштаб для D_QL ----
    ! Исторически в bridge стоял множитель 0.1. Оставлен как именованная
    ! константа, чтобы его можно было обнулить или откалибровать.
    ! Должен быть 1.0 для прямого соответствия 1D-коду!
    real(wp), parameter :: D_QL_scale = 1.0_wp

    ! ---- Диагностика ----
    logical, save :: fp2d_verbose = .TRUE.

   ! print *, '=== FP2D 15:09 ==='

contains

!=======================================================================
subroutine fp2d_reset()
    if (allocated(f2D))    deallocate(f2D)
    if (allocated(D_QL2D)) deallocate(D_QL2D)
    fp2d_initialized = .FALSE.
    fp2d_rel_saved   = .FALSE.
end subroutine fp2d_reset

!=======================================================================
subroutine fp2d_init(use_rel)
    logical, intent(in) :: use_rel
    integer :: jr

    call fp2d_reset()

    allocate(f2D(nx, ny, nr, 2))
    allocate(D_QL2D(nx, ny, nr, 2))
    f2D    = 0.0_wp
    D_QL2D = 0.0_wp

    call get_2d_grid(p_par_grid, p_perp_grid)

    do jr = 1, nr
        call set_physical_params(jr, use_rel)
        call init_maxwell_juttner(f2D(:,:,jr,1), theta_e, vte_over_c)
        call init_maxwell_juttner(f2D(:,:,jr,2), theta_e, vte_over_c)
    end do

    fp2d_initialized = .TRUE.
    fp2d_rel_saved   = use_rel

    if (fp2d_verbose) then
        print *, '=== FP2D init ==='
        print *, '  use_relativistic = ', use_rel
        print *, '  nx, ny, nr       = ', nx, ny, nr
        print *, '  p_max            = ', p_max
        print *, '  dp_par, dp_perp  = ', dp_par, dp_perp
        print *, '  maxval(f2D)      = ', maxval(f2D)
        print *, '  f2D(33,1,1,1)    = ', f2D(33,1,1,1)
        print *, '  f2D(33,3,1,1)    = ', f2D(33,3,1,1)
        print *, '================='
    end if
end subroutine fp2d_init

!=======================================================================
subroutine set_physical_params(jr, use_rel)
    integer,  intent(in) :: jr
    logical,  intent(in) :: use_rel
    real(wp) :: r, vte

    r   = dble(jr)/dble(nr+1)
    vte = fvt(r)
    theta_e    = ft(r)/m_e_c2
    vte_over_c = vte/c_light
    z_eff      = zefff(r)
    e_field    = enorm(jr)
    gamma_coll = 1.0_wp

    call set_solver_params(theta_e, z_eff, e_field, vte_over_c, &
                            gamma_coll, use_rel)
end subroutine set_physical_params

!=======================================================================
subroutine fp2d_driver(dt, nt, use_rel)
    real(wp), intent(in) :: dt
    integer,  intent(in) :: nt
    logical,  intent(in) :: use_rel
    integer  :: jr, k, nt_local, jr_mid
    real(wp) :: f2d_max_before, f2d_max_after, DQL_max
    real(wp) :: dtau_total, dt_local
    real(wp), parameter :: dt0 = 0.01_wp

    if (.not. fp2d_initialized) then
        call fp2d_init(use_rel)
    else if (use_rel .neqv. fp2d_rel_saved) then
        call fp2d_init(use_rel)
    end if

    f2d_max_before = maxval(f2D)

    do jr = 1, nr
        call set_physical_params(jr, use_rel)

        dtau_total = dt * dble(nt) * fst(jr)

        nt_local = 1
        if (dtau_total > dt0) nt_local = 1 + int(dtau_total / dt0)
        dt_local = dtau_total / dble(nt_local)

        do k = 1, 2
        call interp_dql_2d(jr, k)

        if (maxval(D_QL2D(:,:,jr,k)) < 1.0e-20_wp) then
            ! волновой диффузии нет — принудительно возвращаем максвеллиан,
            ! чтобы кинетика по p_par не дрейфовала от начального состояния
            call set_physical_params(jr, use_rel)
            call init_maxwell_juttner(f2D(:,:,jr,k), theta_e, vte_over_c)
            cycle
        end if

        call solve_fokker_planck_2d(f2D(:,:,jr,k), D_QL2D(:,:,jr,k), &
                                     dt_local, nt_local)
        end do
    end do

    f2d_max_after = maxval(f2D)
    DQL_max       = maxval(D_QL2D)

    if (fp2d_verbose) then
        jr_mid = nr/2
        if (jr_mid < 1) jr_mid = 1
        print *, '=== FP2D diagnostics ==='
        print *, '  maxval(D_QL2D)         = ', DQL_max
        print *, '  maxval(f2D) до solve   = ', f2d_max_before
        print *, '  maxval(f2D) после      = ', f2d_max_after
        print *, '  f2D(33,1,jr_mid,1)     = ', f2D(33,1,jr_mid,1)
        print *, '  f2D(33,1,jr_mid,2)     = ', f2D(33,1,jr_mid,2)
        print *, '  asym(33,1) = f2 - f1   = ', &
                    f2D(33,1,jr_mid,2) - f2D(33,1,jr_mid,1)
        print *, '  dtau, nt_local, dt     = ', &
                    dtau_total, nt_local, dt_local
        print *, '==================================='
    end if

    call fp2d_extract_1d(use_rel)
end subroutine fp2d_driver

!=======================================================================
subroutine interp_dql_2d(jr, k)
    integer, intent(in) :: jr, k
    real(wp) :: p_par, D_val
    integer  :: i, j

    D_QL2D(:,:,jr,k) = 0.0_wp

    do j = 1, ny
        do i = 1, nx
            p_par = -p_max + dble(i-1)*dp_par

            if (k == 1 .and. p_par < 0.0_wp) cycle
            if (k == 2 .and. p_par > 0.0_wp) cycle

            D_val = interp_1d(vij(:,jr), dij(:,jr,k), i0, abs(p_par))
            if (D_val > 0.0_wp) then
                D_QL2D(i,j,jr,k) = D_QL_scale * D_val
            end if
        end do
    end do
end subroutine interp_dql_2d

!=======================================================================
function interp_1d(xs, ys, n, x) result(y)
    real(wp), intent(in) :: xs(:), ys(:), x
    integer,  intent(in) :: n
    real(wp) :: y
    integer  :: k

    y = 0.0_wp
    if (x <= xs(1)) then
        y = ys(1); return
    end if
    if (x >= xs(n)) then
        y = ys(n); return
    end if

    do k = 1, n-1
        if (x >= xs(k) .and. x <= xs(k+1)) then
            if (abs(xs(k+1)-xs(k)) > 1.0e-30_wp) then
                y = ys(k) + (ys(k+1)-ys(k))*(x-xs(k))/(xs(k+1)-xs(k))
            else
                y = 0.5_wp*(ys(k)+ys(k+1))
            end if
            return
        end if
    end do
end function interp_1d

!=======================================================================
subroutine fp2d_extract_1d(use_rel)
    logical, intent(in) :: use_rel
    integer  :: jr, k, i
    real(wp) :: p, p_par, f_marg, f_marg_init

    do jr = 1, nr
        do i = 1, i0
            p = vij(i, jr)
            do k = 1, 2
                if (p < 1.0e-30_wp) then
                    p_par = 0.0_wp
                else if (k == 1) then
                    p_par = p
                else
                    p_par = - p 
                end if

                call marginal_p_perp(f2D(:,:,jr,k), p_par, f_marg)
               ! call marginal_maxwell(p_par, jr, use_rel, f_marg_init)
                
                f_marg_init = exp(-0.5_wp * p**2) / sqrt(2.0_wp * pi)

                fij (i, jr, k) = f_marg
                fij0(i, jr, k) = f_marg_init
            end do
        end do
        call compute_dfij(jr)
    end do

    if (fp2d_verbose) then
        print *, '=== FP2D extract (jr_mid) ==='
        print *, '  fij (33,k=1), fij (33,k=2) = ', &
                    fij (33,nr/2,1), fij (33,nr/2,2)
        print *, '  fij0(33,k=1), fij0(33,k=2) = ', &
                    fij0(33,nr/2,1), fij0(33,nr/2,2)
        print *, '  diff k=1, diff k=2         = ', &
                    fij(33,nr/2,1)-fij0(33,nr/2,1), &
                    fij(33,nr/2,2)-fij0(33,nr/2,2)
        print *, '  asym = fij(k=2) - fij(k=1) = ', &
                    fij(33,nr/2,2) - fij(33,nr/2,1)
        print *, '============================='
    end if
end subroutine fp2d_extract_1d

!=======================================================================
subroutine marginal_p_perp(f2, p_par_target, f_marg)
    real(wp), intent(in)  :: f2(nx, ny)
    real(wp), intent(in)  :: p_par_target
    real(wp), intent(out) :: f_marg
    real(wp) :: p_par_i, p_perp, f_interp, integral
    integer  :: i, j
    logical  :: found

    integral = 0.0_wp
    found = .FALSE.

    do i = 1, nx-1
        p_par_i = p_par_grid(i)
        if (p_par_target >= p_par_i .and. &
            p_par_target <= p_par_grid(i+1)) then

            do j = 1, ny
                p_perp = p_perp_grid(j)
                f_interp = f2(i,j) + (f2(i+1,j)-f2(i,j)) * &
                            (p_par_target - p_par_i)/dp_par
                if (j == 1 .or. j == ny) then
                    integral = integral + 0.5_wp * f_interp * 2.0_wp*pi * &
                               p_perp * dp_perp
                else
                    integral = integral + f_interp * 2.0_wp*pi * &
                               p_perp * dp_perp
                end if
            end do
            found = .TRUE.
            exit
        end if
    end do

    if (found) then
        f_marg = integral / norm_2pi32
    else
        f_marg = 0.0_wp
    end if
end subroutine marginal_p_perp

!=======================================================================
subroutine marginal_maxwell(p_par_target, jr, use_rel, f_marg)
    real(wp), intent(in)  :: p_par_target
    integer,  intent(in)  :: jr
    logical,  intent(in)  :: use_rel
    real(wp), intent(out) :: f_marg
    real(wp) :: p_perp, p, gamma, f_val, w, integral
    real(wp), parameter :: dpp = 0.1_wp
    integer  :: j

    call set_physical_params(jr, use_rel)

    integral = 0.0_wp
    do j = 1, 200
        p_perp = dble(j-1)*dpp
        p = sqrt(p_par_target**2 + p_perp**2)

        if (use_rel) then
            gamma = sqrt(1.0_wp + (p*vte_over_c)**2)
            f_val = exp(-(gamma - 1.0_wp)/theta_e)
        else
            f_val = exp(-0.5_wp * p**2)
        end if

        if (j == 1 .or. j == 200) then
            w = 0.5_wp
        else
            w = 1.0_wp
        end if

        integral = integral + w * f_val * 2.0_wp*pi * p_perp * dpp
    end do

    f_marg = integral / norm_2pi32
end subroutine marginal_maxwell

!=======================================================================
subroutine compute_dfij(jr)
    integer, intent(in) :: jr
    real(wp) :: dv
    integer  :: i, k

    do k = 1, 2
        do i = 2, i0-1
            dv = vij(i+1,jr) - vij(i-1,jr)
            if (abs(dv) > 1.0e-30_wp) then
                dfij(i,jr,k) = (fij(i+1,jr,k) - fij(i-1,jr,k))/dv
            end if
        end do
        dfij(1, jr, k) = dfij(2, jr, k)
        dfij(i0,jr, k) = dfij(i0-1, jr, k)
    end do
end subroutine compute_dfij

!=======================================================================
subroutine fp2d_save(time)
    real(wp), intent(in) :: time
    character(len=120)   :: fname
    character(len=8)     :: tag
    integer, parameter   :: iu = 22
    integer              :: jr, k
    logical              :: exist

    if (.not. fp2d_initialized) return

    if (fp2d_rel_saved) then
        tag = 'rel     '
    else
        tag = 'nonrel  '
    end if

    write(fname,'("lhcd/f2d/f2d_", A, "_", f9.7, ".bin")') &
          trim(tag), time

    inquire(file=fname, exist=exist)
    if (exist) then
        !print *, 'fp2d_save: file exists, skip ', trim(fname)
        return
    end if

    open(iu, file=fname, status='new', action='write', &
         access='stream', form='unformatted')
    write(iu) nx, ny, nr
    write(iu) p_par_grid
    write(iu) p_perp_grid
    do jr = 1, nr
        do k = 1, 2
            write(iu) f2D(:,:,jr,k)
        end do
    end do
    close(iu)

    print *, 'fp2d_save: written ', trim(fname)
end subroutine fp2d_save

end module FP2D_bridge_mod