module small_vgrid
    use kind_module   
    use nr_grid, only: MAX_NR
    implicit none
    integer, parameter :: MAX_PT= 1000 ! размер температурной сетки
    real(wp) :: vmid(MAX_PT),vz1(MAX_PT),vz2(MAX_PT)
    !integer  :: ibeg(100),iend(100)

    real(wp) :: vrj(MAX_PT), dj(MAX_PT), djnew(MAX_PT)
    real(wp) :: dj2(MAX_PT),d2j(MAX_PT)

    real(wp) :: vgrid(MAX_PT,MAX_NR), dfundv(MAX_PT,MAX_NR)
    !!common/gridv/vgrid(101,100),dfundv(101,100)

    real(wp) :: dql(MAX_PT, MAX_NR)
    !!
    !real(wp) :: dq1(101,MAX_NR)
    !real(wp) :: dq2(101,MAX_NR)
    real(wp) :: dqi0(50,MAX_NR) 
    !common /alph/ dqi0(50,100)    

    real(wp) :: dncount(MAX_PT, MAX_NR)
    !common/findsigma/dncount(101,100)    

    integer  :: nvpt
    !!common/gridv/nvpt
    integer :: ipt1, ipt2, ipt

    real(wp) :: vlf,vrt,dflf,dfrt
    !common /a0ghp/ vlf,vrt,dflf,dfrt
        
    integer, parameter :: kpt1=20, kpt3=20

contains

    subroutine init_small_vgrid_arrays
        use constants, only: zero
        dql=zero
        !dq1=zero
        !dq2=zero
        dqi0=zero
        dncount=zero  
    end subroutine init_small_vgrid_arrays

    subroutine find_achieved_radial_points(nvpt)
        !!  find achieved radial points jbeg-jend
        use rt_parameters, only : nr
        implicit none
        integer, intent(in) :: nvpt
        integer i, j, jbeg, jend, nvmin, nvach

        nvmin=1 !minimum counted events at a given radius rho
        jbeg=1
        jend=0
        do j=1,nr
            nvach=0
            do i=1,nvpt
                nvach = nvach + dncount(i,j)
            end do
            if (nvach.lt.nvmin) then
                if (jend.eq.0) jbeg = jbeg + 1
            else
                jend=j
            end if
        end do
        if (jend.eq.0.or.jbeg.ge.jend) then
            write(*,*)'failure: jbeg=',jbeg,' jend=',jend 
            pause
            stop
        end if
    end subroutine  

    subroutine distr(vz,j,ifound,fder)
        use lock_module      
        use rt_parameters, only: nr
        implicit none
        integer, intent(in) :: j
        integer, intent(inout) :: ifound
        real*8 vz,fder
        integer i,klo,khi,ierr,nvp
        real*8,dimension(:),allocatable:: vzj,dfdvj
        real(wp) :: dfout
        !real*8 vlf,vrt,dflf,dfrt
        !common /a0ghp/ vlf,vrt,dflf,dfrt
        !common/gridv/vgrid(101,100),dfundv(101,100),nvpt
        nvp=nvpt
        allocate(vzj(nvp),dfdvj(nvp))
        do i=1, nvp
            vzj(i)=vgrid(i,j)
            dfdvj(i)=dfundv(i,j)
        end do
        call lock2(vzj,nvp,vz,klo,khi,ierr)
        if(ierr.eq.0) then !vgrid(1,j) <= vz <= vgrid(nvpt,j)
            call linf(vzj,dfdvj,vz,dfout,klo,khi)
            ifound=klo
            vlf=vzj(klo)
            vrt=vzj(khi)
            fder=dfout
            dflf=dfdvj(klo)
            dfrt=dfdvj(khi)
        else if(ierr.eq.1) then !vz < vgrid(1,j)
            write(*,*)'exception: ierr=1 in distr()'
            pause'next key = stop'
            stop
        else if(ierr.eq.2) then !vz > vgrid(nvpt,j)
            !write(*,*)'exception: ierr=2 in distr()'
            print *, 'vz > vgrid(nvpt,j)'
            ifound=nvp
            vlf=vzj(nvp)
            vrt=vz
            fder=0
            dflf=0
            dfrt=0
            !pause'next key = stop'
            !stop
        else if(ierr.eq.3) then
            write(*,*)'exception in distr, klo=khi=',klo,' j=',j,' nvp=',nvp
            write(*,*)'vz=',vz,' v1=',vzj(1),' v2=',vzj(nvp)
            pause'next key = stop'
            stop
        end if
        deallocate(vzj,dfdvj)
    end   

    subroutine gridvel(v1,v2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
        implicit none
        real(wp), intent(in) :: v1, v2, vmax, cdel
        integer,  intent(in) :: ni1, ni2, ipt1, kpt3
        real(wp), intent(inout) :: vrj(:)
        integer kpt1, kpt2, k
        real(wp) :: v12
        kpt1=ipt1-1
        kpt2=ni1+ni2+1
        do k=1,kpt1  !0<=v<v1
            vrj(k)=dble(k-1)*v1/dble(kpt1)
        end do
        v12=v1+(v2-v1)*cdel
        do k=1,ni1+1 !v1<=v<=v12
            vrj(k+kpt1)=v1+dble(k-1)*(v12-v1)/dble(ni1)
        end do
        do k=2,ni2+1 !!v12<v<=v2
            vrj(k+kpt1+ni1)=v12+dble(k-1)*(v2-v12)/dble(ni2)
        end do     
        do k=1,kpt3  !v2<v<=vmax
            vrj(k+kpt1+kpt2)=v2+dble(k)*(vmax-v2)/dble(kpt3)
        end do
    end    

end module small_vgrid