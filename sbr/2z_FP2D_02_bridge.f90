!=======================================================================
! Модуль FP2D_bridge_mod
! Связка между 1D кодом FRTC и 2D решателем FP2D_solver_mod.
!
! Поддерживает два режима (fp_solver):
!   2 — 2D нерелятивистский (классический Ландау)
!   3 — 2D релятивистский (Беляев-Будкер)
!
! Учитывает раздельно положительный и отрицательный спектры:
!   k = 1 → p_par < 0  (отрицательный N_par)
!   k = 2 → p_par > 0  (положительный N_par)
!
! ВАЖНО: для интерполяции D_QL используется dij(:,jr,k) из nr_grid —
! это уже пересчитанный в физические единицы Фоккера-Планка коэффициент.
! Сетка vij(:,jr) и 2D-сетка (p_par, p_perp) — в одинаковых единицах v_te(r).
!=======================================================================
module FP2D_bridge_mod
    use kind_module,     only: wp
    use FP2D_solver_mod
    use nr_grid,         only: MAX_NR, vij, fij, fij0, dfij, dij
    use maxwell,         only: i0
    use plasma,          only: fvt, ft, zefff, enorm
    use rt_parameters,   only: nr
    implicit none

    ! --- Состояние модуля ---
    real(wp), allocatable, save :: f2D(:,:,:,:)      ! (nx, ny, nr, 2)
    real(wp), allocatable, save :: D_QL2D(:,:,:,:)   ! (nx, ny, nr, 2)
    real(wp), save :: p_par_grid(nx), p_perp_grid(ny)
    logical,  save :: fp2d_initialized = .FALSE.
    logical,  save :: fp2d_rel_saved   = .FALSE.

    ! --- Константы ---
    real(wp), parameter :: c_light = 3.0e10_wp
    real(wp), parameter :: m_e_c2  = 8.187e-7_wp
    real(wp), parameter :: norm_2pi32 = 15.7496099457224_wp  ! (2π)^{3/2}

    ! --- Диагностика ---
    logical,  save :: fp2d_verbose = .TRUE.

contains

!=======================================================================
! Сброс состояния модуля
!=======================================================================
subroutine fp2d_reset()
    if (allocated(f2D))    deallocate(f2D)
    if (allocated(D_QL2D)) deallocate(D_QL2D)
    fp2d_initialized = .FALSE.
    fp2d_rel_saved   = .FALSE.
end subroutine fp2d_reset

!=======================================================================
! Инициализация 2D массивов
!=======================================================================
subroutine fp2d_init(use_rel)
    logical, intent(in) :: use_rel
    integer :: jr

    call fp2d_reset()

    allocate(f2D(nx, ny, nr, 2))
    allocate(D_QL2D(nx, ny, nr, 2))
    f2D    = 0.0_wp
    D_QL2D = 0.0_wp

    ! Сетка импульсов (один раз)
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
        print *, 'use_relativistic = ', use_rel
        print *, 'nx, ny, nr = ', nx, ny, nr
        print *, 'p_max = ', p_max, '  dp_par = ', dp_par, '  dp_perp = ', dp_perp
        print *, 'maxval(f2D) после init = ', maxval(f2D)
        print *, '================='
    end if
end subroutine fp2d_init

!=======================================================================
! Установка физических параметров в точке jr
!=======================================================================
subroutine set_physical_params(jr, use_rel)
    integer,  intent(in) :: jr
    logical,  intent(in) :: use_rel
    real(wp) :: r, vte

    r   = dble(jr)/dble(nr+1)
    vte = fvt(r)                    ! см/с
    theta_e    = ft(r)/m_e_c2       ! T_e/(m_e c^2)
    vte_over_c = vte/c_light
    z_eff      = zefff(r)
    e_field    = enorm(jr)
    gamma_coll = 1.0_wp

    call set_solver_params(theta_e, z_eff, e_field, vte_over_c, &
                            gamma_coll, use_rel)
end subroutine set_physical_params

!=======================================================================
! Основной драйвер: прогон 2D решателя для всех (jr, k)
!=======================================================================
subroutine fp2d_driver(dt, nt, use_rel)
    real(wp), intent(in) :: dt
    integer,  intent(in) :: nt
    logical,  intent(in) :: use_rel
    integer  :: jr, k
    real(wp) :: f2d_max_before, f2d_max_after
    real(wp) :: DQL_max
    integer  :: jr_mid

    ! Переинициализация при смене режима
    if (.not. fp2d_initialized) then
        call fp2d_init(use_rel)
    else if (use_rel .neqv. fp2d_rel_saved) then
        call fp2d_init(use_rel)
    end if

    ! --- Снимок состояния ДО solve ---
    f2d_max_before = maxval(f2D)

    do jr = 1, nr
        call set_physical_params(jr, use_rel)
        do k = 1, 2
            call interp_dql_2d(jr, k, use_rel)
            call solve_fokker_planck_2d(f2D(:,:,jr,k), D_QL2D(:,:,jr,k), dt, nt)
        end do
    end do

    ! --- Снимок состояния ПОСЛЕ solve ---
    f2d_max_after = maxval(f2D)
    DQL_max       = maxval(D_QL2D)

    if (fp2d_verbose) then
        jr_mid = nr/2
        if (jr_mid < 1) jr_mid = 1

        print *, '=== FP2D diagnostics ==='
        print *, 'use_relativistic     = ', use_rel
        print *, 'maxval(D_QL2D)       = ', DQL_max
        print *, 'maxval(f2D) до solve = ', f2d_max_before
        print *, 'maxval(f2D) после    = ', f2d_max_after
        print *, 'f2D(32, 1, jr_mid, 1) = ', f2D(32, 1, jr_mid, 1)
        print *, 'f2D(32, 1, jr_mid, 2) = ', f2D(32, 1, jr_mid, 2)
        print *, 'f2D(48, 1, jr_mid, 1) = ', f2D(48, 1, jr_mid, 1)
        print *, 'f2D(48, 1, jr_mid, 2) = ', f2D(48, 1, jr_mid, 2)
        print *, 'f2D(32,10, jr_mid, 1) = ', f2D(32,10, jr_mid, 1)
        print *, 'f2D(32,10, jr_mid, 2) = ', f2D(32,10, jr_mid, 2)

        ! Симметрия между k=1 и k=2
        print *, 'asym (32,:) = ', f2D(32,1,jr_mid,2) - f2D(32,1,jr_mid,1)
        print *, 'asym (48,:) = ', f2D(48,1,jr_mid,2) - f2D(48,1,jr_mid,1)

        print *, '==================================='
    end if

    call fp2d_extract_1d(use_rel)
end subroutine fp2d_driver

!=======================================================================
! Интерполяция dij (1D, пересчитанный в физические единицы ФП)
! на 2D сетку (p_par, p_perp).
!
! Единицы vij(:,jr) и 2D-сетки совпадают — обе в единицах v_te(r).
! Интерполяция прямая, без пересчёта масштаба.
!=======================================================================
subroutine interp_dql_2d(jr, k, use_rel)
    integer, intent(in) :: jr, k
    logical, intent(in) :: use_rel
    real(wp) :: p_par, D_val
    integer  :: i, j

    D_QL2D(:,:,jr,k) = 0.0_wp

    do j = 1, ny
        do i = 1, nx
            p_par = -p_max + dble(i-1)*dp_par

            ! --- Разделение по знаку p_par ---
            ! k=1 (отрицательный спектр) → взаимодействие при p_par < 0
            ! k=2 (положительный спектр) → взаимодействие при p_par > 0
            if (k == 1 .and. p_par > 0.0_wp) cycle
            if (k == 2 .and. p_par < 0.0_wp) cycle

            D_val = interp_1d(vij(:,jr), dij(:,jr,k), i0, abs(p_par))
            if (D_val > 0.0_wp) then
                D_QL2D(i,j,jr,k) = D_val
            end if
        end do
    end do
end subroutine interp_dql_2d

!=======================================================================
! Линейная интерполяция y(xs) в точку x
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
! Извлечение 1D маргиналов и запись в fij/fij0/dfij
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
                    p_par = -p
                else
                    p_par = p
                end if

                call marginal_p_perp(f2D(:,:,jr,k), p_par, f_marg)
                call marginal_maxwell(p_par, jr, use_rel, f_marg_init)

                fij (i, jr, k) = f_marg
                fij0(i, jr, k) = f_marg_init
            end do
        end do
        call compute_dfij(jr)
    end do

    if (fp2d_verbose) then
        print *, '=== FP2D extract ==='
        print *, 'fij(32, jr_mid, 1)  = ', fij(32, nr/2, 1)
        print *, 'fij(32, jr_mid, 2)  = ', fij(32, nr/2, 2)
        print *, 'fij0(32, jr_mid, 1) = ', fij0(32, nr/2, 1)
        print *, 'fij0(32, jr_mid, 2) = ', fij0(32, nr/2, 2)
        print *, 'fij(48, jr_mid, 1)  = ', fij(48, nr/2, 1)
        print *, 'fij(48, jr_mid, 2)  = ', fij(48, nr/2, 2)
        print *, 'fij0(48, jr_mid, 1) = ', fij0(48, nr/2, 1)
        print *, 'fij0(48, jr_mid, 2) = ', fij0(48, nr/2, 2)
        print *, 'diff fij - fij0 (32) = ', &
                    fij(32,nr/2,2) - fij0(32,nr/2,2)
        print *, 'diff fij - fij0 (48) = ', &
                    fij(48,nr/2,2) - fij0(48,nr/2,2)
        print *, '================='
    end if
end subroutine fp2d_extract_1d

!=======================================================================
! Маргинал 2D распределения по p_perp:
! f_1D(p_par) = 2π ∫ f_2D(p_par, p_perp) p_perp dp_perp / (2π)^{3/2}
!=======================================================================
subroutine marginal_p_perp(f2, p_par_target, f_marg)
    real(wp), intent(in)  :: f2(nx, ny)
    real(wp), intent(in)  :: p_par_target
    real(wp), intent(out) :: f_marg
    real(wp) :: p_par_i, p_perp, f_interp, integral
    integer  :: i, j

    integral = 0.0_wp

    ! Линейная интерполяция по p_par
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

            f_marg = integral / norm_2pi32
            return
        end if
    end do

    f_marg = 0.0_wp
end subroutine marginal_p_perp

!=======================================================================
! Маргинал начального распределения
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
! Производная ∂f/∂v по 1D массивам
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
! Сохранение 2D распределения
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
        print *, 'fp2d_save: file exists, skip ', trim(fname)
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

!=======================================================================
! Диагностика
!=======================================================================
subroutine fp2d_diagnostics(jr, k, J_par, E_kin)
    integer,  intent(in)  :: jr, k
    real(wp), intent(out) :: J_par, E_kin
    real(wp) :: p_par, p_perp, p, gamma, v_par
    integer  :: i, j

    J_par = 0.0_wp
    E_kin = 0.0_wp

    do i = 1, nx
        p_par = p_par_grid(i)
        do j = 1, ny
            p_perp = p_perp_grid(j)
            p = sqrt(p_par**2 + p_perp**2)
            if (p < eps_reg) cycle

            if (fp2d_rel_saved) then
                gamma = sqrt(1.0_wp + (p*vte_over_c)**2)
                v_par = p_par/gamma
            else
                gamma = 1.0_wp
                v_par = p_par
            end if

            J_par = J_par + f2D(i,j,jr,k) * v_par * p_perp * dp_par * dp_perp
            E_kin = E_kin + f2D(i,j,jr,k) * (gamma - 1.0_wp) * &
                    p_perp * dp_par * dp_perp
        end do
    end do

    J_par = J_par * 2.0_wp*pi
    E_kin = E_kin * 2.0_wp*pi
end subroutine fp2d_diagnostics

end module FP2D_bridge_mod