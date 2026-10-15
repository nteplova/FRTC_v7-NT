!=======================================================================
! Модуль FP2D_solver_mod
! 2D решатель уравнения Фоккера-Планка для кода FRTC.
!
! Схема: Ченг-Купер + 9-точечная неявная схема с AF-факторизацией.
! Координаты: (p_par, p_perp) в единицах v_te(r).
!
! Граничное условие при p_perp = 0:
!   (1/p) ∂/∂p (p D ∂f/∂p) → 4 D ∂²f/∂p²
!=======================================================================
module FP2D_solver_mod
    use kind_module, only: wp
    implicit none

    integer, parameter :: nx = 64
    integer, parameter :: ny = 64

    real(wp), parameter :: p_max   = 100.0_wp
    real(wp), parameter :: dp_par  = 2.0_wp*p_max/dble(nx-1)
    real(wp), parameter :: dp_perp = p_max/dble(ny-1)

    integer,  parameter :: max_iter = 100
    real(wp), parameter :: tol_conv = 1.0e-6_wp
    real(wp), parameter :: eps_reg  = 1.0e-12_wp
    real(wp), parameter :: pi       = 3.14159265358979323846_wp

    real(wp), save :: theta_e    = 1.0e-3_wp
    real(wp), save :: z_eff      = 1.0_wp
    real(wp), save :: e_field    = 0.0_wp
    real(wp), save :: vte_over_c = 0.05_wp
    real(wp), save :: gamma_coll = 1.0_wp

    logical, save :: use_relativistic = .TRUE.

    ! --- Флаг диагностики ---
    logical, save :: debug_init = .TRUE.

contains

!=======================================================================
! Установка физических параметров
!=======================================================================
subroutine set_solver_params(theta, z, ef, vc, gc, rel)
    real(wp), intent(in) :: theta, z, ef, vc, gc
    logical,  intent(in) :: rel
    theta_e          = theta
    z_eff            = z
    e_field          = ef
    vte_over_c       = vc
    gamma_coll       = gc
    use_relativistic = rel
end subroutine set_solver_params

!=======================================================================
! Начальное распределение с диагностикой
!=======================================================================
subroutine init_maxwell_juttner(f, theta, vc)
    real(wp), intent(out) :: f(nx, ny)
    real(wp), intent(in)  :: theta, vc
    real(wp) :: p_par, p_perp, p, gamma
    integer  :: i, j

    ! --- Жёсткая диагностика ---
    print *, '>>> init_maxwell_juttner CALLED <<<'
    print *, '    use_relativistic = ', use_relativistic
    print *, '    theta            = ', theta
    print *, '    vc               = ', vc
    print *, '    p_max            = ', p_max
    print *, '    dp_par           = ', dp_par
    ! ------------------------------------

    do j = 1, ny
        p_perp = dble(j-1) * dp_perp
        do i = 1, nx
            p_par = -p_max + dble(i-1) * dp_par
            p = sqrt(p_par**2 + p_perp**2)

            if (use_relativistic) then
                gamma = sqrt(1.0_wp + (p*vc)**2)
                f(i,j) = exp(-(gamma - 1.0_wp)/theta)
            else
                f(i,j) = exp(-0.5_wp * p**2)
            end if

            ! --- Печать в точке (i=33, j=1) ---
            if (i == 33 .and. j == 1) then
                print *, '>>> (i=33,j=1) p = ', p
                print *, '    use_rel = ', use_relativistic
                print *, '    exp(-0.5 p^2)         = ', exp(-0.5_wp * p**2)
                print *, '    exp(-0.5 p^2/theta)   = ', exp(-0.5_wp * p**2 / theta)
                print *, '    f(i,j)                = ', f(i,j)
            end if
        end do
    end do
end subroutine init_maxwell_juttner

!=======================================================================
! Решатель
!=======================================================================
subroutine solve_fokker_planck_2d(f, D_QL, dt, nt)
    real(wp), intent(inout) :: f(nx, ny)
    real(wp), intent(in)    :: D_QL(nx, ny)
    real(wp), intent(in)    :: dt
    integer,  intent(in)    :: nt

    real(wp) :: f_new(nx, ny), A(nx, ny, 9), RHS(nx, ny)
    real(wp) :: residual, max_val
    integer  :: step, iter

    do step = 1, nt
        do iter = 1, max_iter
            call build_operator(f, D_QL, dt, A, RHS)
            call approximate_factorization(A, RHS, f_new)

            residual = maxval(abs(f_new - f))
            max_val  = maxval(abs(f_new)) + 1.0e-30_wp
            residual = residual / max_val
            f = f_new
            if (residual < tol_conv) exit
        end do
    end do
end subroutine solve_fokker_planck_2d

!=======================================================================
! Построение 9-точечного оператора
!=======================================================================
subroutine build_operator(f, D_QL, dt, A, RHS)
    real(wp), intent(in)  :: f(nx, ny), D_QL(nx, ny), dt
    real(wp), intent(out) :: A(nx, ny, 9), RHS(nx, ny)
    real(wp) :: p_par, p_perp, p, gamma, v
    real(wp) :: D_pp, D_tt, F_v, D_par, D_perp, F_par, F_perp
    real(wp) :: delta_p, delta_t, w_p, w_t
    real(wp) :: ap, bp, am, bm, pp, pm
    real(wp) :: D_perp_eff
    integer  :: i, j

    A   = 0.0_wp
    RHS = 0.0_wp

    ! --- Внутренние точки: j = 2..ny-1, i = 2..nx-1 ---
    do j = 2, ny-1
        p_perp = dble(j-1)*dp_perp
        do i = 2, nx-1
            p_par = -p_max + dble(i-1)*dp_par
            p = sqrt(p_par**2 + p_perp**2)
            if (p < eps_reg) p = eps_reg

            if (use_relativistic) then
                gamma = sqrt(1.0_wp + (p*vte_over_c)**2)
                v     = p / gamma
            else
                gamma = 1.0_wp
                v     = p
            end if

            call collision_coeffs(p, v, gamma, D_pp, D_tt, F_v)

            D_perp = (p_perp/p)**2*D_pp + (p_par/p)**2*D_tt
            D_par  = (p_par/p)**2 *D_pp + (p_perp/p)**2*D_tt
            F_perp = (p_perp/p)*F_v
            F_par  = (p_par/p) *F_v

            call cc_weight(F_perp*dp_perp, D_perp, delta_p, w_p)
            call cc_weight(F_par *dp_par,  D_par,  delta_t, w_t)

            ! --- Поток по p_perp ---
            pp = p_perp + 0.5_wp*dp_perp
            pm = p_perp - 0.5_wp*dp_perp
            if (pm < eps_reg) pm = eps_reg

            ap = -D_perp/dp_perp + F_perp*delta_p
            bp =  D_perp/dp_perp + F_perp*(1.0_wp-delta_p)
            am = -D_perp/dp_perp + F_perp*delta_p
            bm =  D_perp/dp_perp + F_perp*(1.0_wp-delta_p)

            A(i,j,2) = A(i,j,2) + dt/(p_perp*dp_perp) * pm * bm
            A(i,j,8) = A(i,j,8) - dt/(p_perp*dp_perp) * pp * ap
            A(i,j,5) = A(i,j,5) - dt/(p_perp*dp_perp)*pp*bp &
                                 + dt/(p_perp*dp_perp)*pm*am

            ! --- Поток по p_par ---
            ap = -D_par/dp_par + F_par*delta_t
            bp =  D_par/dp_par + F_par*(1.0_wp-delta_t)
            am = -D_par/dp_par + F_par*delta_t
            bm =  D_par/dp_par + F_par*(1.0_wp-delta_t)

            A(i,j,4) = A(i,j,4) + dt/dp_par * bm
            A(i,j,6) = A(i,j,6) - dt/dp_par * ap
            A(i,j,5) = A(i,j,5) - dt/dp_par*bp + dt/dp_par*am

            ! --- Квазилинейная диффузия ---
            if (D_QL(i,j) > 0.0_wp) then
                A(i,j,4) = A(i,j,4) - dt*D_QL(i,j)/dp_par**2
                A(i,j,6) = A(i,j,6) - dt*D_QL(i,j)/dp_par**2
                A(i,j,5) = A(i,j,5) + 2.0_wp*dt*D_QL(i,j)/dp_par**2
            end if

            ! --- Электрическое поле ---
            if (e_field > 0.0_wp) then
                A(i,j,4) = A(i,j,4) - dt*e_field/dp_par
                A(i,j,5) = A(i,j,5) + dt*e_field/dp_par
            else if (e_field < 0.0_wp) then
                A(i,j,6) = A(i,j,6) + dt*e_field/dp_par
                A(i,j,5) = A(i,j,5) - dt*e_field/dp_par
            end if

            A(i,j,5) = A(i,j,5) + 1.0_wp/dt
            RHS(i,j) = f(i,j)/dt
        end do
    end do

    ! --- Особая точка j = 1: p_perp = 0 ---
    do i = 2, nx-1
        p_par  = -p_max + dble(i-1)*dp_par
        p = abs(p_par)
        if (p < eps_reg) p = eps_reg

        if (use_relativistic) then
            gamma = sqrt(1.0_wp + (p*vte_over_c)**2)
            v     = p / gamma
        else
            gamma = 1.0_wp
            v     = p
        end if

        call collision_coeffs(p, v, gamma, D_pp, D_tt, F_v)

        D_par = D_pp
        F_par = 0.0_wp
        if (p_par >  eps_reg) F_par =  F_v
        if (p_par < -eps_reg) F_par = -F_v

        call cc_weight(F_par*dp_par, D_par, delta_t, w_t)

        ap = -D_par/dp_par + F_par*delta_t
        bp =  D_par/dp_par + F_par*(1.0_wp-delta_t)
        am = -D_par/dp_par + F_par*delta_t
        bm =  D_par/dp_par + F_par*(1.0_wp-delta_t)

        A(i,1,4) = A(i,1,4) + dt/dp_par * bm
        A(i,1,6) = A(i,1,6) - dt/dp_par * ap
        A(i,1,5) = A(i,1,5) - dt/dp_par*bp + dt/dp_par*am

        D_perp_eff = 4.0_wp * D_tt / dp_perp**2
        A(i,1,5) = A(i,1,5) + dt * D_perp_eff
        A(i,1,8) = A(i,1,8) - dt * D_perp_eff

        if (D_QL(i,1) > 0.0_wp) then
            A(i,1,4) = A(i,1,4) - dt*D_QL(i,1)/dp_par**2
            A(i,1,6) = A(i,1,6) - dt*D_QL(i,1)/dp_par**2
            A(i,1,5) = A(i,1,5) + 2.0_wp*dt*D_QL(i,1)/dp_par**2
        end if

        if (e_field > 0.0_wp) then
            A(i,1,4) = A(i,1,4) - dt*e_field/dp_par
            A(i,1,5) = A(i,1,5) + dt*e_field/dp_par
        else if (e_field < 0.0_wp) then
            A(i,1,6) = A(i,1,6) + dt*e_field/dp_par
            A(i,1,5) = A(i,1,5) - dt*e_field/dp_par
        end if

        A(i,1,5) = A(i,1,5) + 1.0_wp/dt
        RHS(i,1) = f(i,1)/dt
    end do

    ! --- Граничные условия (Dirichlet) ---
    do j = 1, ny
        A(1,j,:)  = 0.0_wp; A(1,j,5)  = 1.0_wp; RHS(1,j)  = 0.0_wp
        A(nx,j,:) = 0.0_wp; A(nx,j,5) = 1.0_wp; RHS(nx,j) = 0.0_wp
    end do
    do i = 1, nx
        A(i,ny,:) = 0.0_wp; A(i,ny,5) = 1.0_wp; RHS(i,ny) = 0.0_wp
    end do

end subroutine build_operator

!=======================================================================
! Коэффициенты столкновений
!=======================================================================
subroutine collision_coeffs(p, v, gamma, D_pp, D_tt, F_v)
    real(wp), intent(in)  :: p, v, gamma
    real(wp), intent(out) :: D_pp, D_tt, F_v
    real(wp) :: vth2, G, vs, gg

    vth2 = theta_e
    G    = gamma_coll
    vs   = max(v, eps_reg)

    if (use_relativistic) then
        gg = max(gamma, 1.0_wp)
        D_pp = G * vth2 / (vs**3 * gg)
        D_tt = G * (z_eff + 1.0_wp) / (2.0_wp * vs * gg)
        F_v  = -G * vth2 / (theta_e * vs**2 * gg)
    else
        D_pp = G * vth2 / vs**3
        D_tt = G * (z_eff + 1.0_wp) / (2.0_wp * vs)
        F_v  = -G * vth2 / (theta_e * vs**2)
    end if
end subroutine collision_coeffs

!=======================================================================
! Вес Ченга-Купера
!=======================================================================
subroutine cc_weight(C_dp, D, delta, w)
    real(wp), intent(in)  :: C_dp, D
    real(wp), intent(out) :: delta, w

    if (D < eps_reg) then
        delta = 0.5_wp; w = 0.0_wp; return
    end if
    w = C_dp / D

    if (abs(w) < 1.0e-6_wp) then
        delta = 0.5_wp - w/12.0_wp + w**3/720.0_wp
    else if (w > 50.0_wp) then
        delta = 1.0_wp - 1.0_wp/w
    else if (w < -50.0_wp) then
        delta = -1.0_wp/w
    else
        delta = 1.0_wp/w - 1.0_wp/(exp(w) - 1.0_wp)
    end if
end subroutine cc_weight

!=======================================================================
! Приближённая факторизация
!=======================================================================
subroutine approximate_factorization(A, RHS, f_new)
    real(wp), intent(in)  :: A(nx, ny, 9), RHS(nx, ny)
    real(wp), intent(out) :: f_new(nx, ny)
    real(wp) :: g(nx, ny)
    real(wp) :: ar(nx), br(nx), cr(nx), dr(nx), xr(nx)
    integer  :: i, j

    do j = 1, ny
        do i = 1, nx
            ar(i) = A(i,j,4); br(i) = A(i,j,5)
            cr(i) = A(i,j,6); dr(i) = RHS(i,j)
        end do
        call tridiag(nx, ar, br, cr, dr, xr)
        do i = 1, nx; g(i,j) = xr(i); end do
    end do

    do i = 1, nx
        do j = 1, ny
            ar(j) = A(i,j,2); br(j) = A(i,j,5)
            cr(j) = A(i,j,8); dr(j) = g(i,j)
        end do
        call tridiag(ny, ar, br, cr, dr, xr)
        do j = 1, ny; f_new(i,j) = xr(j); end do
    end do
end subroutine approximate_factorization

!=======================================================================
! Метод прогонки Томаса
!=======================================================================
subroutine tridiag(n, a, b, c, d, x)
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: a(n), b(n), c(n), d(n)
    real(wp), intent(out) :: x(n)
    real(wp) :: cp(n), dp(n), m
    integer  :: i

    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
    do i = 2, n
        m = b(i) - a(i)*cp(i-1)
        if (abs(m) < 1.0e-30_wp) m = 1.0e-30_wp
        cp(i) = c(i)/m
        dp(i) = (d(i) - a(i)*dp(i-1))/m
    end do
    x(n) = dp(n)
    do i = n-1, 1, -1
        x(i) = dp(i) - cp(i)*x(i+1)
    end do
end subroutine tridiag

!=======================================================================
! Доступ к сетке
!=======================================================================
subroutine get_2d_grid(p_par_out, p_perp_out)
    real(wp), intent(out) :: p_par_out(nx), p_perp_out(ny)
    integer :: i, j
    do i = 1, nx
        p_par_out(i) = -p_max + dble(i-1)*dp_par
    end do
    do j = 1, ny
        p_perp_out(j) = dble(j-1)*dp_perp
    end do
end subroutine get_2d_grid

end module FP2D_solver_mod