!=======================================================================
! Модуль FP2D_bridge_mod
! Связка между 1D кодом FRTC и 2D решателем FP2D_solver_mod.
!
! Поддерживает два режима:
!   fp_solver = 2 — 2D нерелятивистский (классический Ландау)
!   fp_solver = 3 — 2D релятивистский (Беляев-Будкер)
!
! Процедуры:
!   fp2d_driver      — основной прогон 2D решателя
!   fp2d_extract_1d  — извлечение 1D маргиналов в fij/fij0/dfij
!   fp2d_save        — сохранение 2D распределения
!   fp2d_reset       — сброс состояния (при смене режима)
!=======================================================================
module FP2D_bridge_mod
    use kind_module,     only: wp
    use FP2D_solver_mod
    use nr_grid,         only: MAX_NR, vij, fij, fij0, dfij, dij
    use maxwell,         only: i0
    use small_vgrid,     only: vgrid, dql, nvpt, MAX_PT
    use plasma,          only: fvt, ft, zefff, enorm
    use rt_parameters,   only: nr
    implicit none

    ! --- Состояние модуля ---
    real(wp), allocatable, save :: f2D(:,:,:)      ! (nx, ny, nr)
    real(wp), allocatable, save :: D_QL2D(:,:,:)   ! (nx, ny, nr)
    real(wp), save :: p_par_grid(nx), p_perp_grid(ny)
    logical,  save :: fp2d_initialized = .FALSE.
    logical,  save :: fp2d_rel_saved   = .FALSE.

    ! --- Константы ---
    real(wp), parameter :: c_light = 3.0e10_wp
    real(wp), parameter :: m_e_c2  = 8.187e-7_wp

contains

!=======================================================================
! Сброс состояния
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

    allocate(f2D(nx, ny, nr))
    allocate(D_QL2D(nx, ny, nr))
    f2D    = 0.0_wp
    D_QL2D = 0.0_wp

    ! Сетка (один раз)
    call get_2d_grid(p_par_grid, p_perp_grid)

    do jr = 1, nr
        call set_physical_params(jr, use_rel)
        call init_maxwell_juttner(f2D(:,:,jr), theta_e, vte_over_c)
    end do

    fp2d_initialized = .TRUE.
    fp2d_rel_saved   = use_rel
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
! Основной драйвер
!=======================================================================
subroutine fp2d_driver(dt, nt, use_rel)
    real(wp), intent(in) :: dt
    integer,  intent(in) :: nt
    logical,  intent(in) :: use_rel
    integer  :: jr

    ! Переинициализация при смене режима
    if (.not. fp2d_initialized) then
        call fp2d_init(use_rel)
    else if (use_rel .neqv. fp2d_rel_saved) then
        call fp2d_init(use_rel)
    end if

    do jr = 1, nr
        call set_physical_params(jr, use_rel)
        call interp_dql_2d(jr, use_rel)
        call solve_fokker_planck_2d(f2D(:,:,jr), D_QL2D(:,:,jr), dt, nt)
    end do

    call fp2d_extract_1d(use_rel)
end subroutine fp2d_driver

!=======================================================================
! Интерполяция dql из 1D сетки small_vgrid в 2D сетку
!=======================================================================
subroutine interp_dql_2d(jr, use_rel)
    integer, intent(in) :: jr
    logical, intent(in) :: use_rel
    real(wp) :: vto, p_par, p_vte0, D_val
    integer  :: i, j

    ! Пересчёт масштаба: small_vgrid в единицах v_te(0),
    ! 2D сетка — в единицах v_te(r)
    if (use_rel) then
        ! в релятивистском режиме шкала та же
        vto = fvt(dble(jr)/dble(nr+1)) / fvt(0.0_wp)
    else
        vto = fvt(dble(jr)/dble(nr+1)) / fvt(0.0_wp)
    end if

    D_QL2D(:,:,jr) = 0.0_wp

    do j = 1, ny
        do i = 1, nx
            p_par = -p_max + dble(i-1)*dp_par
            ! D_QL локализована при фиксированном |p_par|
            p_vte0 = abs(p_par) * vto
            D_val  = interp_1d(vgrid(:,jr), dql(:,jr), nvpt, p_vte0)
            if (D_val > 0.0_wp) then
                D_QL2D(i,j,jr) = D_val
            end if
        end do
    end do
end subroutine interp_dql_2d

!=======================================================================
! Линейная интерполяция на 1D сетке
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
            if (p < 1.0e-30_wp) then
                ! граничная точка: f(0)
                do k = 1, 2
                    call marginal_p_perp(f2D(:,:,jr), 0.0_wp, f_marg)
                    call marginal_maxwell(0.0_wp, jr, use_rel, f_marg_init)
                    fij (i, jr, k) = f_marg
                    fij0(i, jr, k) = f_marg_init
                end do
            else
                do k = 1, 2
                    if (k == 1) then
                        p_par = -p
                    else
                        p_par = p
                    end if
                    call marginal_p_perp(f2D(:,:,jr), p_par, f_marg)
                    call marginal_maxwell(p_par, jr, use_rel, f_marg_init)
                    fij (i, jr, k) = f_marg
                    fij0(i, jr, k) = f_marg_init
                end do
            end if
        end do
        call compute_dfij(jr)
    end do
end subroutine fp2d_extract_1d

!=======================================================================
! Маргинал 2D распределения по p_perp
! f_1D(p_par) = 2π ∫ f_2D(p_par, p_perp) p_perp dp_perp
!=======================================================================
subroutine marginal_p_perp(f2, p_par_target, f_marg)
    real(wp), intent(in)  :: f2(nx, ny)
    real(wp), intent(in)  :: p_par_target
    real(wp), intent(out) :: f_marg
    real(wp) :: p_par_i, p_perp, f_interp
    integer  :: i, j

    f_marg = 0.0_wp

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
                    f_marg = f_marg + 0.5_wp * f_interp * 2.0_wp*pi * &
                             p_perp * dp_perp
                else
                    f_marg = f_marg + f_interp * 2.0_wp*pi * &
                             p_perp * dp_perp
                end if
            end do
            return
        end if
    end do
end subroutine marginal_p_perp

!=======================================================================
! Маргинал начального распределения (Максвелл или Максвелл-Юттнер)
!=======================================================================
subroutine marginal_maxwell(p_par_target, jr, use_rel, f_marg)
    real(wp), intent(in)  :: p_par_target
    integer,  intent(in)  :: jr
    logical,  intent(in)  :: use_rel
    real(wp), intent(out) :: f_marg
    real(wp) :: p_perp, p, gamma, f_val, w
    real(wp), parameter :: dpp = 0.1_wp
    integer  :: j

    ! Обновить параметры для точки jr (нужно для theta_e, vte_over_c)
    call set_physical_params(jr, use_rel)

    f_marg = 0.0_wp
    do j = 1, 200
        p_perp = dble(j-1)*dpp
        p = sqrt(p_par_target**2 + p_perp**2)

        if (use_rel) then
            gamma = sqrt(1.0_wp + (p*vte_over_c)**2)
            f_val = exp(-(gamma - 1.0_wp)/theta_e)
        else
            f_val = exp(-0.5_wp * p**2 / theta_e)
        end if

        if (j == 1 .or. j == 200) then
            w = 0.5_wp
        else
            w = 1.0_wp
        end if

        f_marg = f_marg + w * f_val * 2.0_wp*pi * p_perp * dpp
    end do
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
! Формат: lhcd/f2d/f2d_<tag>_<time>.bin
!   integer  nx, ny, nr
!   real(wp) f2D(:,:,jr)  для jr = 1..nr
!=======================================================================
subroutine fp2d_save(time)
    real(wp), intent(in) :: time
    character(len=120)   :: fname
    character(len=8)     :: tag
    integer, parameter   :: iu = 22
    integer              :: jr
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
        write(iu) f2D(:,:,jr)
    end do
    close(iu)

    print *, 'fp2d_save: written ', trim(fname)
end subroutine fp2d_save

!=======================================================================
! Диагностика: полный ток и энергия в 2D
!=======================================================================
subroutine fp2d_diagnostics(jr, J_par, E_kin)
    integer,  intent(in)  :: jr
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

            J_par = J_par + f2D(i,j,jr) * v_par * p_perp * dp_par * dp_perp
            E_kin = E_kin + f2D(i,j,jr) * (gamma - 1.0_wp) * &
                    p_perp * dp_par * dp_perp
        end do
    end do

    J_par = J_par * 2.0_wp*pi
    E_kin = E_kin * 2.0_wp*pi
end subroutine fp2d_diagnostics

end module FP2D_bridge_mod