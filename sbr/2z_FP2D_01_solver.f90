!=======================================================================
! Модуль FP2D_solver_mod
! 2D решатель уравнения Фоккера-Планка для кода FRTC.
!
! Численная схема: Douglas ADI (2 шага на dt).
!   Шаг 1: (I - dt*Lx) f*      = f^n + dt*Ly f^n
!   Шаг 2: (I - dt*Ly) f^{n+1} = f*  - dt*Ly f^n
!
! Дискретизация потоков - Chang-Cooper.
! Сила трения F = -D*p (well-balanced форма).
!
! ВАЖНО: D и F вычисляются В ЦЕНТРЕ ячейки (не на грани!) --
! это стабильная форма, проверенная на тестах. Попытка вычислить
! D и F на грани (face-based) приводит к неустойчивости.
!=======================================================================
module FP2D_solver_mod
    use kind_module, only: wp
    implicit none

    integer, parameter :: nx = 129
    integer, parameter :: ny = 129

    real(wp), parameter :: p_max   = 120.0_wp
    real(wp), parameter :: dp_par  = 2.0_wp*p_max/dble(nx-1)
    real(wp), parameter :: dp_perp = p_max/dble(ny-1)

    real(wp), parameter :: eps_reg  = 1.0e-4_wp
    real(wp), parameter :: pi       = 3.14159265358979323846_wp

    real(wp), save :: theta_e    = 1.0e-3_wp
    real(wp), save :: z_eff      = 1.0_wp
    real(wp), save :: e_field    = 0.0_wp
    real(wp), save :: vte_over_c = 0.05_wp
    real(wp), save :: gamma_coll = 1.0_wp

    logical, save :: use_relativistic = .TRUE.

contains

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
subroutine init_maxwell_juttner(f, theta, vc)
    real(wp), intent(out) :: f(nx, ny)
    real(wp), intent(in)  :: theta, vc
    real(wp) :: p_par, p_perp, p, gamma
    integer  :: i, j

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
        end do
    end do
end subroutine init_maxwell_juttner

!=======================================================================
subroutine solve_fokker_planck_2d(f, D_QL, dt_in, nt_in)
    real(wp), intent(inout) :: f(nx, ny)
    real(wp), intent(in)    :: D_QL(nx, ny)
    real(wp), intent(in)    :: dt_in
    integer,  intent(in)    :: nt_in

    real(wp) :: ax(nx), bx(nx), cx(nx)
    real(wp) :: ay(ny), by(ny), cy(ny)
    real(wp) :: rhs(nx, ny), fstar(nx, ny), Lyf(nx, ny)
    real(wp) :: dt
    integer  :: step, i, j, nt

    call compute_stable_dt(D_QL, dt_in, dt, nt)

    do step = 1, nt
        call apply_Ly(f, Lyf)

        ! --- Шаг 1: (I - dt*Lx) f* = f + dt*Ly f ---
        rhs = f + dt * Lyf
        do j = 1, ny
            rhs(1,  j) = 0.0_wp
            rhs(nx, j) = 0.0_wp
        end do
        do j = 1, ny
            call build_x_operator(j, D_QL, dt, ax, bx, cx)
            call tridiag(nx, ax, bx, cx, rhs(:,j), fstar(:,j))
        end do

        ! --- Шаг 2: (I - dt*Ly) f^{n+1} = f* - dt*Ly f ---
        rhs = fstar - dt * Lyf
        do i = 1, nx
            rhs(i, ny) = 0.0_wp
        end do
        do i = 1, nx
            call build_y_operator(i, dt, ay, by, cy)
            call tridiag(ny, ay, by, cy, rhs(i,:), f(i,:))
        end do

        ! Защита от отрицательных значений и выбросов
        do j = 1, ny
            do i = 1, nx
                if (f(i,j) < 0.0_wp) f(i,j) = 0.0_wp
                if (f(i,j) > 1.0e10_wp) f(i,j) = 1.0e10_wp
            end do
        end do
    end do

    do j = 1, ny
        f(1,  j) = 0.0_wp
        f(nx, j) = 0.0_wp
    end do
    do i = 1, nx
        f(i, ny) = 0.0_wp
    end do
end subroutine solve_fokker_planck_2d

!=======================================================================
subroutine compute_stable_dt(D_QL, dt_in, dt, nt)
    real(wp), intent(in)  :: D_QL(nx, ny), dt_in
    real(wp), intent(out) :: dt
    integer,  intent(out) :: nt

    real(wp) :: D_max, dt_max, p, v, gamma, D_tt
    integer  :: i, j

    D_max = 0.0_wp
    do j = 1, ny
        do i = 1, nx
            p = sqrt((-p_max + dble(i-1)*dp_par)**2 &
                   + (dble(j-1)*dp_perp)**2)
            if (p < 0.1_wp) p = 0.1_wp
            if (use_relativistic) then
                gamma = sqrt(1.0_wp + (p*vte_over_c)**2)
                v = p/gamma
            else
                gamma = 1.0_wp; v = p
            end if
            D_tt = gamma_coll*(z_eff+1.0_wp) &
                 / (2.0_wp * max(v, eps_reg) * max(gamma, 1.0_wp))
            if (D_tt > D_max) D_max = D_tt
        end do
    end do

    dt_max = 5.0_wp * dp_perp**2 / max(4.0_wp*D_max, eps_reg)
    nt     = max(1, int(dt_in / dt_max) + 1)
    dt     = dt_in / dble(nt)
end subroutine compute_stable_dt

!=======================================================================
! Столкновительные коэффициенты в точке
!=======================================================================
subroutine collision_at_point(p_par, p_perp, D_par, F_par, D_perp, F_perp)
    real(wp), intent(in)  :: p_par, p_perp
    real(wp), intent(out) :: D_par, F_par, D_perp, F_perp
    real(wp) :: p, gamma, v, D_pp, D_tt, F_v

    p = sqrt(p_par**2 + p_perp**2)
    if (p < eps_reg) p = eps_reg
    if (use_relativistic) then
        gamma = sqrt(1.0_wp + (p*vte_over_c)**2)
        v = p/gamma
    else
        gamma = 1.0_wp
        v = p
    end if

    call collision_coeffs(p, v, gamma, D_pp, D_tt, F_v)

    D_par  = (p_par /p)**2 * D_pp + (p_perp/p)**2 * D_tt
    D_perp = (p_perp/p)**2 * D_pp + (p_par /p)**2 * D_tt

    ! Well-balanced формы
    F_par  = -D_pp  * p_par    ! как в формуле 1Д
    F_perp = -D_pp * p_perp
end subroutine collision_at_point

!=======================================================================
! Стандартный вес Chang-Cooper
!=======================================================================
subroutine cc_weight(F_dp, D, delta, w)
    real(wp), intent(in)  :: F_dp, D
    real(wp), intent(out) :: delta, w

    w = F_dp / max(D, eps_reg)

    if (abs(w) < 1.0e-6_wp) then
        delta = 0.5_wp - w/12.0_wp
    else if (w > 50.0_wp) then
        delta = 1.0_wp - 1.0_wp/w
    else if (w < -50.0_wp) then
        delta = -1.0_wp/w
    else
        delta = 1.0_wp/w - 1.0_wp/(exp(w) - 1.0_wp)
    end if
end subroutine cc_weight

!=======================================================================
subroutine collision_coeffs(p, v, gamma, D_pp, D_tt, F_v)
    real(wp), intent(in)  :: p, v, gamma
    real(wp), intent(out) :: D_pp, D_tt, F_v
    real(wp) :: G, vs, gg
    real(wp), parameter :: v_min_reg = 0.05_wp   ! <-- новое

    G  = gamma_coll
    ! было:  vs = max(v, eps_reg)
    vs = max(v, v_min_reg)

    if (use_relativistic) then
        gg = max(gamma, 1.0_wp)
        D_pp = G / (vs**3 * gg)
        D_tt = G * (z_eff + 1.0_wp) / (2.0_wp * vs * gg)
        F_v  = -G / (vs**2 * gg)
    else
        D_pp = G / vs**3
        D_tt = G * (z_eff + 1.0_wp) / (2.0_wp * vs)
        F_v  = -G / vs**2
    end if
end subroutine collision_coeffs

!=======================================================================
! Явное применение Ly (перпендикулярное направление)
!=======================================================================
subroutine apply_Ly(f, Lyf)
    real(wp), intent(in)  :: f(nx, ny)
    real(wp), intent(out) :: Lyf(nx, ny)
    real(wp) :: p_par, p_perp
    real(wp) :: D_dum, F_dum, D_perp_p, F_perp_p
    real(wp) :: delta, w
    real(wp) :: ap, bp, am, bm, pp, pm
    real(wp) :: S_p, S_m
    real(wp) :: D_pp_tmp, D_tt_tmp, F_v_tmp, p_, v_, g_
    integer  :: i, j

    Lyf = 0.0_wp

    do i = 1, nx
        p_par = -p_max + dble(i-1)*dp_par

        do j = 2, ny-1
            p_perp = dble(j-1)*dp_perp

            call collision_at_point(p_par, p_perp, &
                                    D_dum, F_dum, D_perp_p, F_perp_p)
            call cc_weight(F_perp_p*dp_perp, D_perp_p, delta, w)

            ap =  D_perp_p/dp_perp + F_perp_p*delta
            bp = -D_perp_p/dp_perp + F_perp_p*(1.0_wp - delta)
            am =  D_perp_p/dp_perp + F_perp_p*delta
            bm = -D_perp_p/dp_perp + F_perp_p*(1.0_wp - delta)

            pp = p_perp + 0.5_wp*dp_perp
            pm = max(p_perp - 0.5_wp*dp_perp, 0.0_wp)

            S_p = ap*f(i,j) + bp*f(i,j+1)
            S_m = am*f(i,j-1) + bm*f(i,j)

            Lyf(i,j) = -(pp*S_p - pm*S_m)/(p_perp*dp_perp)
        end do

        ! Ось j=1: простая формула 4·D·(f_2-f_1)/dp²
        p_ = abs(p_par)
        if (p_ < 0.05_wp) p_ = 0.05_wp
        if (use_relativistic) then
            g_ = sqrt(1.0_wp + (p_*vte_over_c)**2)
            v_ = p_/g_
        else
            g_ = 1.0_wp; v_ = p_
        end if
        call collision_coeffs(p_, v_, g_, D_pp_tmp, D_tt_tmp, F_v_tmp)
        Lyf(i,1) = 4.0_wp*D_tt_tmp*(f(i,2) - f(i,1))/dp_perp**2

        Lyf(i,ny) = 0.0_wp
    end do
end subroutine apply_Ly

!=======================================================================
! x-оператор: (I - dt*Lx)
!=======================================================================
subroutine build_x_operator(j, D_QL, dt, a, b, c)
    integer,  intent(in)  :: j
    real(wp), intent(in)  :: D_QL(nx,ny), dt
    real(wp), intent(out) :: a(nx), b(nx), c(nx)

    real(wp) :: p_par, p_perp, dp
    real(wp) :: D_par, F_par, D_dum, F_dum
    real(wp) :: delta, w
    real(wp) :: ap, bp, am, bm
    real(wp) :: alpha, beta, gamma_, ql_coef
    integer  :: i

    a = 0.0_wp; b = 0.0_wp; c = 0.0_wp
    p_perp = dble(j-1)*dp_perp
    dp = dp_par

    do i = 2, nx-1
        p_par = -p_max + dble(i-1)*dp

        call collision_at_point(p_par, p_perp, &
                                D_par, F_par, D_dum, F_dum)
        call cc_weight(F_par*dp, D_par, delta, w)

        ap =  D_par/dp + F_par*delta
        bp = -D_par/dp + F_par*(1.0_wp - delta)
        am =  D_par/dp + F_par*delta
        bm = -D_par/dp + F_par*(1.0_wp - delta)

        alpha  =  am/dp
        beta   = -(ap - bm)/dp
        gamma_ = -bp/dp

        ql_coef = D_QL(i,j)
        alpha  = alpha  + ql_coef/dp**2
        beta   = beta   - 2.0_wp*ql_coef/dp**2
        gamma_ = gamma_ + ql_coef/dp**2

        if (e_field > 0.0_wp) then
            alpha = alpha - e_field/dp
            beta  = beta  + e_field/dp
        else if (e_field < 0.0_wp) then
            beta  = beta  - e_field/dp
            gamma_= gamma_+ e_field/dp
        end if

        a(i) = -dt*alpha
        b(i) =  1.0_wp - dt*beta
        c(i) = -dt*gamma_
    end do

    a(1)  = 0.0_wp;  b(1)  = 1.0_wp;  c(1)  = 0.0_wp
    a(nx) = 0.0_wp;  b(nx) = 1.0_wp;  c(nx) = 0.0_wp
end subroutine build_x_operator

!=======================================================================
! y-оператор: (I - dt*Ly)
!=======================================================================
subroutine build_y_operator(i, dt, a, b, c)
    integer,  intent(in)  :: i
    real(wp), intent(in)  :: dt
    real(wp), intent(out) :: a(ny), b(ny), c(ny)

    real(wp) :: p_par, p_perp
    real(wp) :: D_perp, F_perp, D_dum, F_dum
    real(wp) :: delta, w
    real(wp) :: ap, bp, am, bm, pp, pm
    real(wp) :: alpha, beta, gamma_
    real(wp) :: D_pp_tmp, D_tt_tmp, F_v_tmp, p_, v_, g_
    integer  :: j

    a = 0.0_wp; b = 0.0_wp; c = 0.0_wp
    p_par = -p_max + dble(i-1)*dp_par

    do j = 2, ny-1
        p_perp = dble(j-1)*dp_perp

        call collision_at_point(p_par, p_perp, &
                                D_dum, F_dum, D_perp, F_perp)
        call cc_weight(F_perp*dp_perp, D_perp, delta, w)

        ap =  D_perp/dp_perp + F_perp*delta
        bp = -D_perp/dp_perp + F_perp*(1.0_wp - delta)
        am =  D_perp/dp_perp + F_perp*delta
        bm = -D_perp/dp_perp + F_perp*(1.0_wp - delta)

        pp = p_perp + 0.5_wp*dp_perp
        pm = max(p_perp - 0.5_wp*dp_perp, 0.0_wp)

        alpha  =  pm*am/(p_perp*dp_perp)
        beta   = -(pp*ap - pm*bm)/(p_perp*dp_perp)
        gamma_ = -pp*bp/(p_perp*dp_perp)

        a(j) = -dt*alpha
        b(j) =  1.0_wp - dt*beta
        c(j) = -dt*gamma_
    end do

    ! Ось p_perp=0: простая формула
    block
        p_ = abs(p_par)
        if (p_ < 0.05_wp) p_ = 0.05_wp
        if (use_relativistic) then
            g_ = sqrt(1.0_wp + (p_*vte_over_c)**2)
            v_ = p_/g_
        else
            g_ = 1.0_wp; v_ = p_
        end if
        call collision_coeffs(p_, v_, g_, D_pp_tmp, D_tt_tmp, F_v_tmp)
        a(1) = 0.0_wp
        b(1) = 1.0_wp + 4.0_wp*dt*D_tt_tmp/dp_perp**2
        c(1) = -4.0_wp*dt*D_tt_tmp/dp_perp**2
    end block

    a(ny) = 0.0_wp; b(ny) = 1.0_wp; c(ny) = 0.0_wp
end subroutine build_y_operator

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

!=======================================================================
! Диагностика состояния 2D распределения: интеграл, max, min, позиция пика
!=======================================================================
subroutine fp2d_solver_stats(f, label)
    real(wp), intent(in) :: f(nx, ny)
    character(len=*), intent(in) :: label
    real(wp) :: total, fmax, fmin
    integer  :: i, j, i_max, j_max
    logical  :: has_nan

    total = 0.0_wp
    fmax  = 0.0_wp
    fmin  = 1.0e30_wp
    i_max = 1; j_max = 1

    do j = 1, ny
        do i = 1, nx
            if (f(i,j) > fmax) then
                fmax = f(i,j); i_max = i; j_max = j
            end if
            if (f(i,j) < fmin) fmin = f(i,j)
            total = total + f(i,j) * 2.0_wp*pi*dble(j-1)*dp_perp &
                                * dp_par * dp_perp
        end do
    end do

    has_nan = .FALSE.
    if (fmax /= fmax) has_nan = .TRUE.
    if (fmin /= fmin) has_nan = .TRUE.
    if (total /= total) has_nan = .TRUE.

    print '(A,A,A)',      '  [', trim(label), ']'
    print '(A,ES15.5E3,A,2F10.3)', '    ∫f·2πp⊥dp = ', total, &
          '   peak at (p_par,p_perp) = ', &
          -p_max + dble(i_max-1)*dp_par, dble(j_max-1)*dp_perp
    print '(A,ES15.5E3,A,ES15.5E3)', '    fmax = ', fmax, &
          '   fmin = ', fmin
    if (has_nan) print *, '    *** NaN/Inf detected! ***'
end subroutine fp2d_solver_stats

end module FP2D_solver_mod