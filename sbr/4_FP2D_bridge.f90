!=======================================================================
! Модуль FP2D_bridge_mod: связь между 1D кодом FRTC и 2D решателем
!=======================================================================
module FP2D_bridge_mod
    use kind_module
    use FP2D_Solver
    use maxwell,      only: i0, vij, fij, fij0, dfij, dij
    use current,      only: dql
    use iterator_mod, only: vgrid, nvpt
    use plasma,       only: fvt, ft, zefff, enorm
    use rt_parameters, only: nr
    implicit none

    real(wp), parameter :: c_light = 3.0e10_wp
    real(wp), parameter :: m_e_c2  = 8.187e-7_wp   ! m_e c^2 в эрг

    ! --- 2D состояния ---
    real(wp) :: f2D(nx, ny, 100, 2)
    real(wp) :: f2D_init(nx, ny, 100, 2)
    real(wp) :: D_QL_2D(nx, ny, 100, 2)
    logical, save :: fp2d_initialized = .FALSE.

contains

!=======================================================================
! Инициализация 2D массивов
!=======================================================================
subroutine fp2d_init_all()
    integer :: jr, k
    real(wp) :: theta

    do jr = 1, nr
        theta = ft(dble(jr)/dble(nr+1)) / m_e_c2
        do k = 1, 2
            call init_maxwell_juttner_2D(f2D_init(:,:,jr,k), theta)
            f2D(:,:,jr,k) = f2D_init(:,:,jr,k)
        end do
    end do
    fp2d_initialized = .TRUE.
end subroutine fp2d_init_all

!=======================================================================
! Драйвер: решает 2D уравнение для всех радиальных точек
!=======================================================================
subroutine fokkerplanck_2D_driver(dt, n_steps)
    real(wp), intent(in) :: dt
    integer,  intent(in) :: n_steps
    integer  :: jr, k
    real(wp) :: theta, vte, r_loc

    if (.not. fp2d_initialized) call fp2d_init_all()

    do jr = 1, nr
        r_loc = dble(jr) / dble(nr+1)
        theta = ft(r_loc) / m_e_c2
        vte   = fvt(r_loc)

        T_e        = theta
        Z_eff      = zefff(r_loc)
        e_field    = enorm(jr)
        Gamma_coll = 1.0_wp

        do k = 1, 2
            call interpolate_dql_2D(jr, vte, D_QL_2D(:,:,jr,k))
            call solve_fokker_planck_2D(f2D(:,:,jr,k), D_QL_2D(:,:,jr,k), &
                                        dt, n_steps)
            call extract_1D_from_2D(jr, k, vte)
        end do
    end do
end subroutine fokkerplanck_2D_driver

!=======================================================================
! Интерполяция dql (в v_te единицах) на 2D сетку (в m_e c)
!=======================================================================
subroutine interpolate_dql_2D(jr, vte, D_2D)
    integer,  intent(in)  :: jr
    real(wp), intent(in)  :: vte
    real(wp), intent(out) :: D_2D(nx, ny)
    real(wp) :: p_par, v_par, D_val
    integer  :: i, j, klo, khi

    D_2D = 0.0_wp
    do j = 1, ny
        do i = 1, nx
            p_par = -p_par_max + dble(i-1) * dp_par
            v_par = abs(p_par) * c_light / vte
            D_val = interp_1D_grid(vgrid(:,jr), nvpt, dql(:,jr), v_par, klo, khi)
            D_2D(i,j) = D_val
        end do
    end do
end subroutine interpolate_dql_2D

!=======================================================================
! Извлечение 1D маргинала из 2D
!=======================================================================
subroutine extract_1D_from_2D(jr, k, vte)
    integer,  intent(in)  :: jr, k
    real(wp), intent(in)  :: vte
    real(wp) :: p_target
    integer  :: i

    do i = 1, i0
        p_target = vij(i,jr) * vte / c_light
        if (k == 1) p_target = -p_target
        fij(i,jr,k) = interp_f2D_axis(f2D(:,:,jr,k), p_target)
    end do
end subroutine extract_1D_from_2D

!=======================================================================
! Вспомогательные функции интерполяции
!=======================================================================
real(wp) function interp_1D_grid(xs, n, ys, x, klo, khi) result(y)
    real(wp), intent(in)  :: xs(:), ys(:), x
    integer,  intent(in)  :: n
    integer,  intent(out) :: klo, khi