!=======================================================================
! Модуль FokkerPlanck_module
! Драйвер решения уравнения Фоккера-Планка.
!
! Поддерживает четыре режима через fp_solver:
!   0 — 1D Савельев    (нерелятивистский)
!   1 — 1D Ченг-Купер  (нерелятивистский)
!   2 — 2D нерелятивистский  (классический Ландау + Максвелл)
!   3 — 2D релятивистский    (Беляев-Будкер + Максвелл-Юттнер)
!
! Вызов: fokkerplanck_compute(time, TAU)
!   time — текущее время (для имён файлов)
!   TAU  — шаг по времени
!=======================================================================
module FokkerPlanck_module
    use kind_module
    implicit none

contains

!=======================================================================
! Главная процедура: расчёт функции распределения на шаге TAU
!=======================================================================
subroutine fokkerplanck_compute(time, TAU)
    use FokkerPlanck1D_mod, only: FokkerPlanck1D
    use FP2D_bridge_mod,    only: fp2d_driver, fp2d_save
    use rt_parameters,      only: nr, save_interval, fp_solver
    use writer_module,      only: write_v_array, binary_write_array
    use maxwell,            only: jindex, kindex, flag_d0
    use nr_grid,            only: vij, fij, fij0, dij, dfij
    use plasma,             only: fvt, enorm, fst
    implicit none

    type(FokkerPlanck1D) :: fokker_planck
    real(wp), intent(in) :: time, TAU

    real(wp) :: dtstep, dtau, dt
    real(wp) :: xend, h, r
    real(wp) :: znak, dt0
    integer, parameter :: ntau = 10
    integer :: i, j, it, nt, k

    real(wp), save :: last_time_saving = 0.0_wp

    dt0 = 0.1_wp
    dtstep = TAU / dble(ntau)

    print *, 'fokkerplanck_compute: time =', time, ' fp_solver =', fp_solver

    ! =================================================================
    ! ВЫБОР РЕШАТЕЛЯ
    ! =================================================================
    select case (fp_solver)

    ! -----------------------------------------------------------------
    case (0, 1)
    ! ---- 1D решатель (Савельев = 0, Ченг-Купер = 1) ----
    ! -----------------------------------------------------------------
        do j = 1, nr
            jindex = j
            dtau = dtstep * fst(j)
            nt = 1
            if (dtau > dt0) nt = 1 + int(dtau / dt0)
            dt = dtau / nt
            r = dble(j) / dble(nr+1)
            xend = 3.0e10_wp / fvt(r)

            do k = 1, 2
                kindex = k
                znak = 2.0_wp * dble(k) - 3.0_wp  ! k=1 → -1, k=2 → +1

                ! --- 1. Начальное (максвелловское) распределение ---
                flag_d0 = .TRUE.
                fokker_planck = FokkerPlanck1D(znak*enorm(j), xend, &
                                                vij(:,j), fij0(:,j,k))
                call fokker_planck%init_zero_diffusion
                do i = 1, ntau
                    call fokker_planck%solve_time_step(dt, nt)
                end do
                fij0(:,j,k) = fokker_planck%f

                ! --- 2. Основное распределение с диффузией ---
                flag_d0 = .FALSE.
                fokker_planck = FokkerPlanck1D(znak*enorm(j), xend, &
                                                vij(:,j), fij(:,j,k))
                call fokker_planck%init_diffusion(dij(:,j,k))
                do i = 1, ntau
                    call fokker_planck%solve_time_step(dt, nt)
                end do
                fij(:,j,k) = fokker_planck%f
            end do
        end do

    ! -----------------------------------------------------------------
    case (2)
    ! ---- 2D нерелятивистский ----
    ! -----------------------------------------------------------------
        call fp2d_driver(dtstep, ntau, use_rel=.FALSE.)

    ! -----------------------------------------------------------------
    case (3)
    ! ---- 2D релятивистский ----
    ! -----------------------------------------------------------------
        call fp2d_driver(dtstep, ntau, use_rel=.TRUE.)

    ! -----------------------------------------------------------------
    case DEFAULT
    ! -----------------------------------------------------------------
        write(*,'(A,I0)') 'ERROR: unknown fp_solver = ', fp_solver
        error stop 1

    end select

    ! =================================================================
    ! СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
    ! =================================================================
    if ((time - last_time_saving) >= save_interval) then

        ! --- 1D массивы (для всех режимов) ---
        call binary_write_array(vij, fij0(:,1:nr,:), time, 'maxwell_fij0')
        call write_v_array(vij, fij (:,1:nr,:), time, 'maxwell')
        call write_v_array(vij, dfij(:,1:nr,:), time, 'f_derivative')
        call write_v_array(vij,  dij(:,1:nr,:), time, 'diffusion')

        ! --- 2D массивы (только для 2D режимов) ---
        if (fp_solver == 2 .or. fp_solver == 3) then
            call fp2d_save(time)
        end if

        last_time_saving = time
    end if

end subroutine fokkerplanck_compute

!=======================================================================
! Инициализация диффузии для схемы Савельева
! (вспомогательная процедура, вызывается из FokkerPlanck1D)
!=======================================================================
subroutine init_diffusion(h, n, vj, dj, d1, d2, d3)
    use lock_module
    implicit none
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: h
    real(wp), dimension(:), intent(in)  :: vj, dj
    real(wp), dimension(:), intent(out) :: d1, d2, d3
    real(wp), dimension(:), allocatable :: xx
    integer :: i0
    integer :: i, klo, khi, ierr, klo1, khi1
    integer :: klo2, klo3, khi2, khi3, ierr1, ierr2, ierr3

    i0 = size(vj)

    allocate(xx(n+1))
    do i = 1, n+1
        xx(i) = h/2.0_wp + h*dble(i-1)
    end do

    do i = 1, n+1
        call lock(vj, i0, xx(i), klo1, khi1, ierr1)
        call lock(vj, i0, xx(i)-h/2.0_wp, klo2, khi2, ierr2)
        call lock(vj, i0, xx(i)+h/2.0_wp, klo3, khi3, ierr3)

        if (ierr1 == 1 .or. ierr2 == 1 .or. ierr3 == 1) then
            write(*,*) 'lock error in init_diffusion'
            write(*,*) 'i=', i, ' xx(i)=', xx(i)
            write(*,*) 'vj(1)=', vj(1), ' vj(i0)=', vj(i0)
            error stop 2
        end if

        d1(i) = dj(klo1)
        d2(i) = dj(klo2)
        d3(i) = dj(klo3)
    end do

    deallocate(xx)

end subroutine init_diffusion

end module FokkerPlanck_module