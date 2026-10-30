!=======================================================================
! Модуль FP2D_bridge_mod
! Связка между 1D кодом FRTC и 2D решателем FP2D_solver_mod.
!
! Режимы (fp_solver):
!   2 — 2D нерелятивистский
!   3 — 2D релятивистский
!
! Отображение k ↔ знак p_par:
!   k=1 (ispectr=+1, положительный спектр) → D_QL прикладывается к p_par > 0
!   k=2 (ispectr=-1, отрицательный спектр) → D_QL прикладывается к p_par < 0
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
    real(wp), allocatable, save :: f2D0(:,:,:,:)     ! (nx, ny, nr, 2) — с D_QL=0
    real(wp), allocatable, save :: D_QL2D(:,:,:,:)   ! (nx, ny, nr, 2)
    real(wp), save :: p_par_grid(nx), p_perp_grid(ny)
    logical,  save :: fp2d_initialized = .FALSE.
    logical,  save :: fp2d_rel_saved   = .FALSE.

    real(wp), parameter :: c_light = 3.0e10_wp
    real(wp), parameter :: m_e_c2  = 8.187e-7_wp

    ! ---- Множитель нормировки маргинала ----
    ! f_1D(v) = exp(-v²/2) / sqrt(2π)
    real(wp), parameter :: norm_2pi32 = 15.7496099457224_wp   ! (2π)^{3/2}

    ! ---- Масштаб для D_QL ----
    ! Должен быть 1.0 для прямого соответствия 1D-коду.
    real(wp), parameter :: D_QL_scale = 1.0_wp

    ! ---- Диагностика: печать ОДИН раз за прогон ----
    logical, save :: fp2d_verbose     = .FALSE.
    integer, save :: fp2d_diag_level  = 0
    logical, save :: fp2d_diag_done   = .FALSE.

contains

!=======================================================================
subroutine fp2d_reset()
    if (allocated(f2D))    deallocate(f2D)
    if (allocated(f2D0))   deallocate(f2D0)
    if (allocated(D_QL2D)) deallocate(D_QL2D)
    fp2d_initialized = .FALSE.
    fp2d_rel_saved   = .FALSE.
    fp2d_verbose     = .FALSE.
    fp2d_diag_level  = 0
    fp2d_diag_done   = .FALSE.
end subroutine fp2d_reset

!=======================================================================
subroutine fp2d_init(use_rel)
    logical, intent(in) :: use_rel
    integer :: jr

    call fp2d_reset()

    allocate(f2D(nx, ny, nr, 2))
    allocate(f2D0(nx, ny, nr, 2))
    allocate(D_QL2D(nx, ny, nr, 2))
    f2D    = 0.0_wp
    f2D0   = 0.0_wp
    D_QL2D = 0.0_wp

    call get_2d_grid(p_par_grid, p_perp_grid)

    do jr = 1, nr
        call set_physical_params(jr, use_rel)
        call init_maxwell_juttner(f2D(:,:,jr,1), theta_e, vte_over_c)
        call init_maxwell_juttner(f2D(:,:,jr,2), theta_e, vte_over_c)
        call init_maxwell_juttner(f2D0(:,:,jr,1), theta_e, vte_over_c)
        call init_maxwell_juttner(f2D0(:,:,jr,2), theta_e, vte_over_c)
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
        print *, '  maxval(f2D0)     = ', maxval(f2D0)
        print *, '  diag_level       = ', fp2d_diag_level
        call fp2d_solver_stats(f2D(:,:,1,1),  'f2D  jr=1 k=1')
        call fp2d_solver_stats(f2D0(:,:,1,1), 'f2D0 jr=1 k=1')
        print *, '================='
    end if

    if (fp2d_verbose .and. fp2d_diag_level >= 4) then
        print *, '=== running clean self-test (D_QL=0, E=0) ==='
        call fp2d_conservation_test(nsteps=20, dt_step=0.01_wp, &
                                    use_rel=use_rel, tol_rel=1.0e-3_wp)
        print *, '=== end of self-test ==='
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
    integer  :: jr, k, nt_local, jr_res, i_res_l
    real(wp) :: f2d_max_before, f2d_max_after, DQL_max
    real(wp) :: dtau_total, dt_local
    real(wp), parameter :: dt0     = 1.0_wp
    real(wp), parameter :: rho_res = 0.05_wp   ! ρ≈0.05 — резонанс в 1D

    ! --- Инициализация; печатаем init-диагностику один раз ---
    if (.not. fp2d_initialized) then
        fp2d_verbose    = .TRUE.
        fp2d_diag_level = 2
        call fp2d_init(use_rel)
        fp2d_diag_done  = .FALSE.   ! разрешить один блок print_all
    else if (use_rel .neqv. fp2d_rel_saved) then
        call fp2d_init(use_rel)
    end if

    f2d_max_before = maxval(f2D)

    block
    real(wp) :: D_zero(nx, ny)
    D_zero = 0.0_wp

    do jr = 1, nr
        call set_physical_params(jr, use_rel)

        dtau_total = dt * dble(nt) * fst(jr)

        nt_local = 1
        if (dtau_total > dt0) nt_local = 1 + int(dtau_total / dt0)
        dt_local = dtau_total / dble(nt_local)

        do k = 1, 2
            ! --- 1. f2D0: только столкновения ---
            call solve_fokker_planck_2d(f2D0(:,:,jr,k), D_zero, &
                                         dt_local, nt_local)

            ! --- 2. Основное распределение с волнами ---
            call interp_dql_2d(jr, k)

            if (maxval(D_QL2D(:,:,jr,k)) < 1.0e-20_wp) then
                f2D(:,:,jr,k) = f2D0(:,:,jr,k)
                cycle
            end if

            call solve_fokker_planck_2d(f2D(:,:,jr,k), D_QL2D(:,:,jr,k), &
                                         dt_local, nt_local)
        end do
    end do
    end block

    f2d_max_after = maxval(f2D)
    DQL_max       = maxval(D_QL2D)

    jr_res = max(1, min(nr, nint(rho_res * dble(nr+1))))

    ! --- Сводка по слою (один раз за прогон) ---
    if (fp2d_verbose) then
        i_res_l = nint((5.0_wp + p_max) / dp_par) + 1
        if (i_res_l < 1)  i_res_l = 1
        if (i_res_l > nx) i_res_l = nx
        print *, '=== FP2D diagnostics ==='
        print *, '  jr_res, rho_res        = ', jr_res, &
                 dble(jr_res)/dble(nr+1)
        print *, '  p_par at i_res         = ', p_par_grid(i_res_l)
        print *, '  maxval(D_QL2D)         = ', DQL_max
        print *, '  maxval(f2D) до solve   = ', f2d_max_before
        print *, '  maxval(f2D) после      = ', f2d_max_after
        print *, '  f2D (i_res,1,jr_res,1) = ', f2D (i_res_l,1,jr_res,1)
        print *, '  f2D0(i_res,1,jr_res,1) = ', f2D0(i_res_l,1,jr_res,1)
        print *, '  f2D (i_res,1,jr_res,2) = ', f2D (i_res_l,1,jr_res,2)
        print *, '  f2D0(i_res,1,jr_res,2) = ', f2D0(i_res_l,1,jr_res,2)
        print *, '  asym k=1 (+p): f2D - f2D0 = ', &
                    f2D(i_res_l,1,jr_res,1) - f2D0(i_res_l,1,jr_res,1)
        print *, '  asym k=2 (-p): f2D - f2D0 = ', &
                    f2D(nx+1-i_res_l,1,jr_res,2) - &
                    f2D0(nx+1-i_res_l,1,jr_res,2)
        print *, '  dtau, nt_local, dt     = ', &
                    dtau_total, nt_local, dt_local
        print *, '==================================='
    end if

    call fp2d_extract_1d(use_rel)

    ! --- Развёрнутая диагностика резонансного слоя — один раз ---
    if (fp2d_diag_level >= 2 .and. .not. fp2d_diag_done) then
        call fp2d_print_all(jr_res, 1, 'jr_res k=1 (pos)')
        call fp2d_print_all(jr_res, 2, 'jr_res k=2 (neg)')
        if (fp2d_diag_level >= 3) then
            call fp2d_print_2d_slice(jr_res, 1, 'jr_res k=1')
            call fp2d_print_2d_slice(jr_res, 2, 'jr_res k=2')
        end if
        fp2d_diag_done = .TRUE.
    end if

    ! --- Выключаем диагностику после первого прохода ---
    fp2d_verbose    = .FALSE.
    fp2d_diag_level = 0

end subroutine fp2d_driver

!=======================================================================
subroutine interp_dql_2d(jr, k)
    integer, intent(in) :: jr, k
    real(wp) :: p_par, D_val, dij_max
    integer  :: i, j, n_nonzero

    D_QL2D(:,:,jr,k) = 0.0_wp
    dij_max = maxval(dij(:,jr,k))
    if (dij_max <= 0.0_wp) return

    n_nonzero = 0
    do j = 1, ny
        do i = 1, nx
            p_par = -p_max + dble(i-1)*dp_par

            ! k=1 → p_par > 0 ;  k=2 → p_par < 0
            if (k == 1 .and. p_par < 0.0_wp) cycle
            if (k == 2 .and. p_par > 0.0_wp) cycle

            ! Усреднение dij по 2D-ячейке — интеграл сохраняется точно
            D_val = integrate_dij_cell(vij(:,jr), dij(:,jr,k), i0, &
                                        max(abs(p_par)-0.5_wp*dp_par, 0.0_wp), &
                                        abs(p_par)+0.5_wp*dp_par)

            if (.not. (D_val >= 0.0_wp)) D_val = 0.0_wp

            if (D_val > 0.0_wp) then
                D_QL2D(i,j,jr,k) = D_QL_scale * D_val
                n_nonzero = n_nonzero + 1
            end if
        end do
    end do

    if (fp2d_verbose) then
        print '(A,I0,A,I0,A,I0,A,ES12.4)', &
              '  interp_dql_2d: jr=', jr, '  k=', k, &
              '  nnz=', n_nonzero, '  max=', maxval(D_QL2D(:,:,jr,k))
    end if
end subroutine

!=======================================================================
function integrate_dij_cell(xs, ys, n, x_lo, x_hi) result(y)
    !! Среднее dij по ячейке [x_lo, x_hi]. Интеграл сохраняется точно.
    real(wp), intent(in) :: xs(:), ys(:), x_lo, x_hi
    integer,  intent(in) :: n
    real(wp) :: y
    integer  :: k
    real(wp) :: x1, x2, y1, y2, dx, xl, xr

    y  = 0.0_wp
    xl = max(x_lo, xs(1))
    xr = min(x_hi, xs(n))
    if (xr <= xl) return

    do k = 1, n-1
        x1 = max(xl, xs(k))
        x2 = min(xr, xs(k+1))
        if (x2 <= x1) cycle

        dx = xs(k+1) - xs(k)
        if (dx > 1.0e-30_wp) then
            y1 = ys(k) + (ys(k+1)-ys(k))*(x1-xs(k))/dx
            y2 = ys(k) + (ys(k+1)-ys(k))*(x2-xs(k))/dx
        else
            y1 = 0.5_wp*(ys(k)+ys(k+1))
            y2 = y1
        end if
        y = y + 0.5_wp*(y1+y2)*(x2-x1)
    end do

    y = y / (xr - xl)
end function integrate_dij_cell

!=======================================================================
subroutine fp2d_extract_1d(use_rel)
    logical, intent(in) :: use_rel
    integer  :: jr, k, i
    real(wp) :: p, p_par, f_marg, f_marg_init

    do jr = 1, nr
        if (fp2d_verbose) then
            do k = 1, 2
                print *, 'jr=', jr, ' k=', k, &
                    ' min dfij =', minval(dfij(:,jr,k)), &
                    ' max|dfij| =', maxval(abs(dfij(:,jr,k))), &
                    ' max dij  =', maxval(dij(:,jr,k)), &
                    ' max D_QL2D =', maxval(D_QL2D(:,:,jr,k))
            end do
        end if
        do i = 1, i0
            p = vij(i, jr)
            do k = 1, 2
                if (p < 1.0e-30_wp) then
                    p_par = 0.0_wp
                else if (k == 1) then
                    p_par = p
                else
                    p_par = -p
                end if

                call marginal_p_perp(f2D (:,:,jr,k), p_par, f_marg)
                call marginal_p_perp(f2D0(:,:,jr,k), p_par, f_marg_init)

                fij (i, jr, k) = f_marg
                fij0(i, jr, k) = f_marg_init
            end do
        end do
        ! call compute_dfij(jr)
        ! ВАЖНО: compute_dfij НЕ вызывается намеренно.
        ! В 1D-коде dfij тоже не обновляется — остаётся начальным максвелловским.
        ! Чтобы 1D и 2D были в одинаковых условиях при вычислении dfundv
        ! (plasma::calculate_dfundv), оставляем dfij фиксированным.
    end do
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
                if (j == 1) then
                    integral = integral + f_interp * 0.25_wp * pi * dp_perp**2
                else if (j == ny) then
                    integral = integral + 0.5_wp * f_interp * 2.0_wp*pi * p_perp * dp_perp
                else
                    integral = integral + f_interp * 2.0_wp*pi * p_perp * dp_perp
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

    write(fname,'("lhcd/fp_2d/fp2d_", A, "_", f9.7, ".bin")') &
          trim(tag), time

    inquire(file=fname, exist=exist)
    if (exist) return

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
subroutine fp2d_find_i(p_target, i_out)
    real(wp), intent(in)  :: p_target
    integer,  intent(out) :: i_out
    i_out = nint((p_target + p_max) / dp_par) + 1
    if (i_out < 1)  i_out = 1
    if (i_out > nx) i_out = nx
end subroutine fp2d_find_i

!=======================================================================
subroutine fp2d_print_all(jr, k, label)
    integer,          intent(in) :: jr, k
    character(len=*), intent(in) :: label
    call fp2d_print_marginal_1d(jr, k, label)
    call fp2d_print_profiles (jr, k, label)
    call fp2d_print_moments  (jr, k, label)
    call fp2d_print_dql      (jr, k, label)
end subroutine fp2d_print_all

!=======================================================================
subroutine fp2d_print_marginal_1d(jr, k, label)
    integer,          intent(in) :: jr, k
    character(len=*), intent(in) :: label
    real(wp) :: v, f_marg, f_marg0, f_1d, f0_1d, diff
    integer  :: idx, i
    real(wp), parameter :: vs(8) = (/ 0.5_wp, 1.0_wp, 2.0_wp, 3.0_wp, &
                                       5.0_wp, 7.0_wp, 10.0_wp, 15.0_wp /)
    real(wp) :: vmax_v, dv_v

    vmax_v = vij(i0, jr)
    dv_v   = vmax_v / dble(i0-1)

    print *, '--- marginal vs 1D (jr=', jr, ', k=', k, ', ', trim(label), ') ---'
    print '(A12,5A15)', '    v', 'f_2D_marg', 'f_2D0_marg', &
                        'f_1D', 'f0_1D', 'd(2D-1D)'
    do idx = 1, 8
        v = vs(idx)
        if (k == 1) then
            call marginal_p_perp(f2D (:,:,jr,k),  v, f_marg)
            call marginal_p_perp(f2D0(:,:,jr,k),  v, f_marg0)
        else
            call marginal_p_perp(f2D (:,:,jr,k), -v, f_marg)
            call marginal_p_perp(f2D0(:,:,jr,k), -v, f_marg0)
        end if
        i = nint(v / dv_v) + 1
        if (i < 1)  i = 1
        if (i > i0) i = i0
        f_1d  = fij (i, jr, k)
        f0_1d = fij0(i, jr, k)
        diff  = f_marg - f_1d
        print '(F12.4,5ES15.5E3)', v, f_marg, f_marg0, f_1d, f0_1d, diff
    end do
end subroutine fp2d_print_marginal_1d

!=======================================================================
subroutine fp2d_print_profiles(jr, k, label)
    integer,          intent(in) :: jr, k
    character(len=*), intent(in) :: label
    real(wp) :: p_samp(10), p, f_val, f0_val
    integer  :: idx, i
    real(wp), parameter :: ps(10) = (/ 0.0_wp, 0.5_wp, 1.0_wp, 2.0_wp, 3.0_wp, &
                                       5.0_wp, 7.0_wp, 10.0_wp, 20.0_wp, 40.0_wp /)
    if (k == 1) then
        p_samp = ps
    else
        p_samp = -ps
    end if

    print *, '--- f2D profile (j_p⊥=1, jr=', jr, ', k=', k, ', ', &
             trim(label), ') ---'
    print '(A12,3A15)', '    p_par', 'f2D', 'f2D0', 'f2D-f2D0'
    do idx = 1, 10
        p = p_samp(idx)
        call fp2d_find_i(p, i)
        f_val  = f2D (i, 1, jr, k)
        f0_val = f2D0(i, 1, jr, k)
        print '(F12.4,3ES15.5E3)', p, f_val, f0_val, f_val - f0_val
    end do
end subroutine fp2d_print_profiles

!=======================================================================
subroutine fp2d_compute_moments(f2, n_val, mom_par, mom_sq)
    real(wp), intent(in)  :: f2(nx, ny)
    real(wp), intent(out) :: n_val, mom_par, mom_sq
    real(wp) :: p_par, p_perp, w
    integer  :: i, j

    n_val = 0.0_wp; mom_par = 0.0_wp; mom_sq = 0.0_wp

    do j = 1, ny
        p_perp = dble(j-1)*dp_perp
        if (j == 1) then
            w = 0.25_wp * pi * dp_perp**2
        else
            w = 2.0_wp*pi * p_perp * dp_perp
            if (j == ny) w = 0.5_wp * w
        end if
        w = w * dp_par

        do i = 1, nx
            p_par = -p_max + dble(i-1)*dp_par
            if (i == 1 .or. i == nx) then
                n_val   = n_val   + 0.5_wp * w * f2(i,j)
                mom_par = mom_par + 0.5_wp * w * p_par * f2(i,j)
                mom_sq  = mom_sq  + 0.5_wp * w * (p_par**2 + p_perp**2) * f2(i,j)
            else
                n_val   = n_val   + w * f2(i,j)
                mom_par = mom_par + w * p_par * f2(i,j)
                mom_sq  = mom_sq  + w * (p_par**2 + p_perp**2) * f2(i,j)
            end if
        end do
    end do

    n_val   = n_val   / norm_2pi32
    mom_par = mom_par / norm_2pi32
    mom_sq  = mom_sq  / norm_2pi32
end subroutine fp2d_compute_moments

subroutine fp2d_print_moments(jr, k, label)
    integer,          intent(in) :: jr, k
    character(len=*), intent(in) :: label
    real(wp) :: n_val, n0_val, mp, mp0, mq, mq0
    real(wp) :: n_1d, mp_1d, mq_1d, n0_1d, mp0_1d, mq0_1d
    real(wp) :: v, dv
    integer  :: i

    call fp2d_compute_moments(f2D (:,:,jr,k), n_val, mp, mq)
    call fp2d_compute_moments(f2D0(:,:,jr,k), n0_val, mp0, mq0)

    n_1d = 0.0_wp; mp_1d = 0.0_wp; mq_1d = 0.0_wp
    n0_1d = 0.0_wp; mp0_1d = 0.0_wp; mq0_1d = 0.0_wp
    do i = 1, i0-1
        dv = vij(i+1,jr) - vij(i,jr)
        v  = 0.5_wp*(vij(i+1,jr) + vij(i,jr))
        n_1d   = n_1d   + fij (i,jr,k) * dv
        mp_1d  = mp_1d  + fij (i,jr,k) * v * dv
        mq_1d  = mq_1d  + fij (i,jr,k) * v**2 * dv
        n0_1d  = n0_1d  + fij0(i,jr,k) * dv
        mp0_1d = mp0_1d + fij0(i,jr,k) * v * dv
        mq0_1d = mq0_1d + fij0(i,jr,k) * v**2 * dv
    end do

    print *, '--- moments (jr=', jr, ', k=', k, ', ', trim(label), ') ---'
    print '(A12,4A15)', '    ', '2D', '2D0', '1D', '1D0'
    print '(A,4ES15.5E3)', '  ∫f    : ', n_val, n0_val, n_1d, n0_1d
    print '(A,4ES15.5E3)', '  ∫vf   : ', mp,    mp0,    mp_1d, mp0_1d
    print '(A,4ES15.5E3)', '  ∫v²f  : ', mq,    mq0,    mq_1d, mq0_1d
    print '(A,ES15.5E3)', '  Δmom(2D) = ', mp - mp0
    print '(A,ES15.5E3)', '  Δmom(1D) = ', mp_1d - mp0_1d
end subroutine fp2d_print_moments

!=======================================================================
subroutine fp2d_print_dql(jr, k, label)
    integer,          intent(in) :: jr, k
    character(len=*), intent(in) :: label
    real(wp) :: dql_max, dql_int, p_perp, w
    real(wp) :: dij_max, dij_int, dv
    integer  :: i, j, i_max, j_max, n_above

    dql_max = 0.0_wp
    dql_int = 0.0_wp
    i_max = 1; j_max = 1

    ! Проход 1: max, интеграл, координаты пика
    do j = 1, ny
        do i = 1, nx
            p_perp = dble(j-1)*dp_perp
            w = 2.0_wp*pi*p_perp*dp_par*dp_perp
            dql_int = dql_int + w * D_QL2D(i,j,jr,k)
            if (D_QL2D(i,j,jr,k) > dql_max) then
                dql_max = D_QL2D(i,j,jr,k)
                i_max = i; j_max = j
            end if
        end do
    end do

    ! Проход 2: количество точек > 0.1·max (порог теперь известен)
    n_above = 0
    do j = 1, ny
        do i = 1, nx
            if (D_QL2D(i,j,jr,k) > 0.1_wp * dql_max) &
                n_above = n_above + 1
        end do
    end do

    dij_max = maxval(dij(:,jr,k))
    dij_int = 0.0_wp
    do i = 1, i0-1
        dv = vij(i+1,jr) - vij(i,jr)
        dij_int = dij_int + dij(i,jr,k) * dv
    end do

    print *, '--- D_QL (jr=', jr, ', k=', k, ', ', trim(label), ') ---'
    print '(A,ES15.5E3,A,2F10.3)', &
          '  max D_QL2D        : ', dql_max, &
          '   at (p_par,p_perp) = ', &
          -p_max + dble(i_max-1)*dp_par, dble(j_max-1)*dp_perp
    print '(A,ES15.5E3)', '  ∫ D_QL2D 2πp⊥dp  : ', dql_int
    print '(A,I6)',       '  точек >0.1·max    : ', n_above
    print '(A,ES15.5E3)', '  max dij (1D)      : ', dij_max
    print '(A,ES15.5E3)', '  ∫ dij dv (1D)     : ', dij_int
    print '(A,ES15.5E3)', '  dql_int/dij_int   : ', &
                          dql_int / max(dij_int, 1.0e-30_wp)
end subroutine fp2d_print_dql

!=======================================================================
subroutine fp2d_print_2d_slice(jr, k, label)
    integer,          intent(in) :: jr, k
    character(len=*), intent(in) :: label
    real(wp) :: p_par_s
    real(wp) :: p_perp_v(5)
    integer  :: ii, jj, i
    real(wp), parameter :: pp = 5.0_wp, pl = 2.0_wp

    do jj = 1, 5
        p_perp_v(jj) = pl * dble(jj-1)
    end do

    print *, '--- f2D slice (jr=', jr, ', k=', k, ', ', &
             trim(label), ') ---'
    print '(A12,5A14)', 'p_par\p_perp', '0', '2', '4', '6', '8'
    do ii = 1, 5
        if (k == 1) then
            p_par_s = pp * ii
        else
            p_par_s = -pp * ii
        end if
        call fp2d_find_i(p_par_s, i)
        write(*, '(F12.2,5ES14.4)') p_par_s, &
            (f2D(i, fp2d_find_j_val(p_perp_v(jj)), jr, k), jj=1, 5)
    end do
end subroutine fp2d_print_2d_slice

!=======================================================================
integer function fp2d_find_j_val(p_target) result(j_out)
    real(wp), intent(in) :: p_target
    j_out = nint(p_target / dp_perp) + 1
    if (j_out < 1)  j_out = 1
    if (j_out > ny) j_out = ny
end function fp2d_find_j_val

!=======================================================================
subroutine fp2d_conservation_test(nsteps, dt_step, use_rel, tol_rel)
    integer,  intent(in) :: nsteps
    real(wp), intent(in) :: dt_step
    logical,  intent(in) :: use_rel
    real(wp), intent(in), optional :: tol_rel

    real(wp) :: f_test(nx,ny)
    real(wp) :: D_zero(nx,ny)
    real(wp) :: fmax0, fmax1, I0, I1, mp, mq
    real(wp) :: rel_fmax, rel_I, tol
    integer  :: istep

    tol = 1.0e-3_wp
    if (present(tol_rel)) tol = tol_rel

    call set_physical_params(1, use_rel)
    call set_solver_params(theta_e, z_eff, 0.0_wp, vte_over_c, &
                            gamma_coll, use_rel)
    call init_maxwell_juttner(f_test, theta_e, vte_over_c)
    D_zero = 0.0_wp

    call compute_2d_integrals(f_test, I0, mp, mq)
    fmax0 = maxval(f_test)

    print *, '=== FP2D conservation test (E=0, D_QL=0) ==='
    print '(A,I4,A,ES15.5E3,A,ES15.5E3)', &
          '  step=', 0, '   fmax=', fmax0, '   I=', I0
    call fp2d_boundary_flux_diagnostic(f_test, 't=0')

    do istep = 1, nsteps
        call solve_fokker_planck_2d(f_test, D_zero, dt_step, 1)
        call compute_2d_integrals(f_test, I1, mp, mq)
        fmax1 = maxval(f_test)
        rel_fmax = abs(fmax1 - fmax0)/max(abs(fmax0), 1.0e-30_wp)
        rel_I    = abs(I1    - I0   )/max(abs(I0   ), 1.0e-30_wp)
        print '(A,I4,2ES15.5E3,A,ES10.3,A,ES10.3)', &
              '  step=', istep, fmax1, I1, &
              '   dmax=', rel_fmax, '   dI=', rel_I
        fmax0 = fmax1; I0 = I1
    end do

    call fp2d_boundary_flux_diagnostic(f_test, 't=end')
    print *, '=== end of conservation test ==='
    print '(A,ES15.5E3)', '  total rel. drift in maxval = ', rel_fmax
    print '(A,ES15.5E3)', '  total rel. drift in I      = ', rel_I
    if (rel_fmax > tol .or. rel_I > tol) then
        print *, '*** WARNING: drift exceeds tolerance ', tol
    end if
end subroutine fp2d_conservation_test

!=======================================================================
subroutine compute_2d_integrals(f2, N, Mp, Mq)
    real(wp), intent(in)  :: f2(nx, ny)
    real(wp), intent(out) :: N, Mp, Mq
    real(wp) :: p_par, p_perp, w
    integer  :: i, j

    N = 0.0_wp; Mp = 0.0_wp; Mq = 0.0_wp
    do j = 1, ny
        p_perp = dble(j-1)*dp_perp
        if (j == 1) then
            w = 0.25_wp * pi * dp_perp**2
        else
            w = 2.0_wp*pi * p_perp * dp_perp
            if (j == ny) w = 0.5_wp * w
        end if
        w = w * dp_par

        do i = 1, nx
            p_par = -p_max + dble(i-1)*dp_par
            if (i == 1 .or. i == nx) then
                N  = N  + 0.5_wp * w * f2(i,j)
                Mp = Mp + 0.5_wp * w * p_par      * f2(i,j)
                Mq = Mq + 0.5_wp * w * (p_par**2 + p_perp**2) * f2(i,j)
            else
                N  = N  + w * f2(i,j)
                Mp = Mp + w * p_par      * f2(i,j)
                Mq = Mq + w * (p_par**2 + p_perp**2) * f2(i,j)
            end if
        end do
    end do
end subroutine compute_2d_integrals

end module FP2D_bridge_mod