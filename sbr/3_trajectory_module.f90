module trajectory_module
    use kind_module
    use trajectory_data
    implicit none

    !integer, parameter :: max_num_trajectories = 30000
    !type(Trajectory), target ::  trajectories(max_num_trajectories)
    integer :: number_of_trajectories
    type(Trajectory), allocatable, target ::  trajectories(:)
contains



    !real(wp) function rini(xm, tet, point, ifail) !sav2009
    subroutine rini(traj, point, iw0)
        !! вычисление начальной точки входа луча
        use constants, only : zero
        use rt_parameters, only : inew, nr, iw, spectrum_coordinate_system
        use spectrum_mod, only : SpectrumPoint
        use decrements, only : dhdnr 
        use dispersion_module, only: ivar, yn3, izn, znakstart
        use dispersion_module, only: find_all_roots, dhdomega
        use metrics, only: g22, g33, co, si
        use metrics, only: calculate_metrics
        use driver_module, only: irs
        use trajectory_data
        implicit none
        
        class(Trajectory), intent(inout) :: traj
        type(SpectrumPoint), intent(in)  :: point
        integer, intent(in)              :: iw0
        real(wp)  :: xm
        real(wp)  :: tet 

        integer :: ntry
        real(wp) :: pa, prt, prm, hr 

        real(wp),  parameter :: rhostart=1.d0
        integer,   parameter :: ntry_max=5
        integer :: num_roots
        real(wp) :: xnr_root(4)
        real(wp) :: znak
        irs = 1
        iw = iw0
        izn = 1
        hr = 1.d0/dble(nr+1)
        tet = traj%tetin

        ntry = 0
        pa = rhostart
        do while (ntry.lt.ntry_max.and.pa.ge.2d0*hr)
            pa = rhostart-hr*dble(ntry)-1.d-4
            ntry = ntry+1

            ! вычисление g22 и g33
            call calculate_metrics(pa, tet)

            select case(spectrum_coordinate_system)
            case(0) ! toroidal coordinates
                yn3 = point%Ntor*dsqrt(g33)
                xm  = point%Npol*dsqrt(g22)
            case(1) ! magnetic coordinates
                yn3 = point%Ntor*dsqrt(g33)/co 
                xm  = point%Npol*dsqrt(g22)/si
            case DEFAULT
                print *, 'bad spectrum coordinate system'
                stop
            end select  

            num_roots = find_all_roots(pa,xm,tet,xnr_root)
             
            if (num_roots>0) then
                ! определение znakstart
                znak = dhdomega(pa,tet,xnr_root(1),xm)
                izn = 1
                if(-znak*dhdnr.gt.zero) then
                    izn=-1
                    znak = dhdomega(pa,tet,xnr_root(2),xm)
                    if(-znak*dhdnr.gt.zero) then
                        write(*,*)'Exception: both modes go outward !!'
                        stop
                    end if
                endif
                traj%rin = pa
                traj%tetzap = tet
                traj%xmzap = xm
                traj%rzap = pa
                traj%yn3zap = yn3 
                traj%irszap = irs 
                traj%iwzap = iw
                traj%iznzap = izn
                traj%znakstart = znak ! оказалось,что нужно запоминать. znakstart используется в ext4
                return
            end if
        end do
        !print *, 'error: no roots'
        traj%mbad = 1 ! плохоая траектория
    end      


subroutine init_trajectories(iw0, spectr)
    use constants
    use rt_parameters, only: max_number_of_traj
    use rt_parameters, only: nr, kv, ntet
    use driver_module
    use plasma, only: tet1, tet2, gap_tet_minus, gap_tet_plus
    use spectrum_mod
    implicit none
    integer,intent(in)         :: iw0
    type (Spectrum),intent(in) :: spectr
    type (SpectrumPoint) point
    integer itet, inz, tet_count
    real(wp) htet
    real(wp) tetin
    real(wp) power_coeff
    if (allocated(trajectories)) deallocate(trajectories)
    allocate(trajectories(max_number_of_traj))

    if (ntet.ne.1) htet = (tet2-tet1)/(ntet-1)
    !-----------------------------------------
    !    find initial radius for a trajectory
    !    on the 1th iteration
    !-----------------------------------------
    number_of_trajectories = 0 
    tet_count = 0 
    do itet = 1,ntet
        tetin = tet1+htet*(itet-1)
        !if (abs(tetin)<abs(tet1/4)) cycle
        if ((gap_tet_minus<tetin).and.(tetin<gap_tet_plus)) cycle
        tet_count = tet_count + 1
        do inz = 1, spectr%size
            number_of_trajectories = number_of_trajectories+1
            if (number_of_trajectories>max_number_of_traj) then
                print *, "number_of_trajectories = ", number_of_trajectories
                print *, "The limit on the number of trajectories has been exceeded"
                stop
            endif     
            current_trajectory => trajectories(number_of_trajectories)
            point = spectr%data(inz)
            call current_trajectory%init(tetin, inz)
            call rini(current_trajectory, point, iw0)
        enddo
    enddo
    power_coeff = float(ntet)/float(tet_count)
    do inz = 1, spectr%size
        point = spectr%data(inz)
        point%power = point%power * power_coeff
    enddo

end subroutine 

subroutine write_trajectories(tview, ispectr,nnz,ntet) !sav2008
!!!writing trajectories into a file
    use constants
    use approximation
    use plasma
    use decrements, only: pdec1,pdec2,pdec3,pdecv,pdecal,dfdv
    use decrements, only: zatukh
    use rt_parameters, only :  nr, itend0, kv, nmaxm, traj_len_seved, save_interval
    use small_vgrid, only : dflf, dfrt, distr
    use driver_module !, only: jrad, iww, izz, length
    use trajectory_data
    implicit none
    
    real(wp), intent(in) :: tview

    integer, intent(in) :: ispectr, nnz, ntet  !sav#

    type(Trajectory) :: traj
    type(TrajectoryPoint) ::tp
    !common /bcef/ ynz,ynpopq
    !common /vth/ vthc(length),poloidn(length)
    real(wp) vthcg,npoli
    !common /a0ghp/ vlf,vrt,dflf,dfrt
    
    integer i, n, itr, ntraj
    integer jrc,nturn,ib,ie,jr,ifast,idir,iv
    integer jznak,jdlt,mn,mm,jchek,itet,inz
    integer traj_size
    integer, parameter :: unit_bias = 10
    integer, parameter :: m=7
    real(wp), parameter :: pleft=1.d-10 !m may be chaged together with name(m)

    real(wp) :: h, xr, xdl, xdlp, xly, xlyp, xgm, xgmp, th
    real(wp) :: x, xx, z, zz, pl, pc, pa
    real(wp) :: pdec1z, pdec3z, pintld, pintal
    real(wp) :: cotet, sitet
    real(wp) :: v, refr, dek3, parn, argum
    real(wp) :: df, powpr, powd, powal, pil, pic
    real(wp) :: powcol, pia, pt, denom, powdamped, domin, fff
    real(wp) :: pow 
    real(wp), save  :: last_time_saving_pos = 0
    real(wp), save  :: last_time_saving_neg = 0
    real(wp) :: delta_time
    character(14) folder
    character(40) ver_fn
    character(40) fname

    print *, 'view_time=',tview
    !print *, name(m)
    if (ispectr>0) then
        folder = "lhcd/traj/pos/"
        delta_time = tview - last_time_saving_pos
    else
        folder = "lhcd/traj/neg/"
        delta_time = tview - last_time_saving_neg
    endif

    write(ver_fn,'(A, "v2")') folder
    print *, ver_fn
    open(1,file=ver_fn)
    write (1,*) 'version 2'
    close(1)

    write(fname,'(A, f9.7,".dat")') folder, tview
    
    print *, 'traj file:', fname, delta_time
    if (delta_time<save_interval) then 
        print *,  "skip saving"
        return
    endif
    if (ispectr>0) then
        last_time_saving_pos = tview
    else
        last_time_saving_neg = tview 
    endif
    if (traj_len_seved == 0) then
        return
    endif

    h=1d0/dble(nr+1)

    open(1,file=fname)

    ntraj=0 !sav2008
    !do itr=1,nnz*ntet !sav2008
    do itr=1, number_of_trajectories
        pow=1.d0
        pl=zero
        pc=zero
        pa=zero
        pdec1=zero
        pdec1z=zero
        pdec3=zero
        pdec3z=zero
        pdecv=zero
        pintld=zero
        pintal=zero
        jrc=nr+1
        jznak=-1
        nturn=1

        traj = trajectories(itr)
        call traj%write_info(1)

        if(traj%mbad.eq.0) then 
            ntraj=ntraj+1
            write(1,3) !write header 
            if (traj_len_seved< 0) then
                traj_size= traj%size
            else
                traj_size = min(traj%size, traj_len_seved)
            endif
            do i=1, traj_size
                tp = traj%points(i)
                v  = tp%vel
                jr = tp%jrad
                refr  = tp%perpn
                npoli = tp%poloidn
                ifast = tp%iww
                vthcg = tp%vthc
                idir  = tp%izz
                dek3  = zero
                th = tp%tetai
                parn = tp%xnpar
                if(itend0.gt.0) then
                    argum=clt/(refr*valfa)
                    dek3=zatukh(argum,abs(jr),vperp,kv)
                end if

                call distr(v,abs(jr),iv,df)
                if(jr.lt.0) then    !case of turn
                    jr=-jr
                    !variant          pintld=-dland(i)*df
                    !!          pintld=-dland(i)*(dflf+dfrt)/2d0
                    pintld = dabs(tp%dland*(dflf+dfrt)/2d0)
                    pdec2  = dexp(-2d0*tp%dcoll)
                    pintal = dabs(tp%dalf*dek3)
                else
                    pdec2 = tp%dcoll
                    pdecv = tp%dland
                    !!          pdec1=-pdecv*df
                    pdec1 = dabs(pdecv*df)
                    pdec3 = dabs(tp%dalf*dek3)
                    pintld = (pdec1+pdec1z)/2d0*h
                    pintal = (pdec3+pdec3z)/2d0*h
                    pdec1z = pdec1
                    pdec3z = pdec3
                end if
                powpr=pow
                powd=pow*dexp(-2d0*pintld)
                powcol=powd*pdec2
                powal=powcol*dexp(-2d0*pintal)
                pow=powal
                pil=pintld
                pic=.5d0*dabs(dlog(pdec2))
                pia=pintal
                pt=1.d0-pow  !total absorbed power
                denom=pil+pic+pia
                powdamped=1.d0-dexp(-2.d0*denom)
                domin=powpr*powdamped
                if(denom.ne.zero) then
                    fff=domin/denom
                    pl=pl+dabs(pil*fff)  !el. Landau absorbed power
                    pc=pc+dabs(pic*fff)  !el. collisions absorbed power
                    pa=pa+dabs(pia*fff)  !alpha Landau absorbed power
                end if
                xr= tp%rho !h*dble(jr)
                cotet=dcos(th)
                sitet=dsin(th)
                xdl=fdf(xr,cdl,ncoef,xdlp)
                xly=fdf(xr,cly,ncoef,xlyp)
                xgm=fdf(xr,cgm,ncoef,xgmp)
                xx=-xdl+xr*cotet-xgm*sitet**2
                zz=xr*xly*sitet
                x=(r0+rm*xx)/1d2
                z=(z0+rm*zz)/1d2 
                jdlt=jr-jrc
                jrc=jr
                if(jdlt*jznak.lt.0.and.nturn.lt.m-1) then
                    nturn=nturn+1
                    jznak=-jznak
                end if

                        !   R, Z, rho, theta, N_par, N_pol, P_tot,P_land, P_coll, vth,  'slow=1','idir=1','ipow', driver, N_traj 
                write(1, 7) x, z, xr,  th,    parn,  npoli, pt,   pl,     pc,     vthcg, ifast,   idir, tp%ipow,  tp%driver,  itr
                
                if(pt.ge.1d0-pleft) go to 11 !maximal absorbed power along a ray
            end do
11          continue
            write (1,*)
        end if
    end do


    close(1)


1     format(2x,'N_traj',3x,'mbad',6x,'theta',9x,'Npar',9x,'rho_start')
2     format('R_pass',4x,'Ptot',6x,'Pland',6x,'Pcoll',8x,'Pa',7x,'dPtot',6x,'dPland',5x,'dPcoll',6x,'dPa')
3     format(5x,'R',10x,'Z',11x,'rho',8x,'theta',7x,'N_par',7x,'N_pol',6x,'P_tot',7x,'P_land',6x,'P_coll',6x,'vth',4x,'slow=1',4x,'out=1',4x,'ipow',4x,'driver',4x'N_traj',6x)
4     format(i3,5x,8(f6.3,5x))
5     format(6(e13.6,3x))
6     format(2(i6,2x),4(e13.6,1x))
7     format(10(e11.4,1x),i5,2x,i5,2x,i3,2x,i5,2x,i5)
8     format('after radial pass=',i3,2x,' P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
9     format('Total passes:           P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
20    format('written time slice (seconds) =',f9.3)
    end

    !subroutine traj(xm0, tet0, xbeg, nmax, nb1, nb2, pabs) 
    subroutine tracing(traj, nmax, nb1, nb2, pow, pabs) 
        use constants, only : tiny1
        use rt_parameters, only: eps, rrange, hdrob, nr, ipri, iw
        use dispersion_module, only: izn, yn3
        use dispersion_module, only: extd4, disp2
        use dispersion_module, only: find_all_roots, find_all_roots_simple
        use driver_module, only: im4, hrad, irs, iabsorp, iznzz, iwzz, irszz, rzz
        use driver_module, only: tetzz, xmzz
        use driver_module, only: driver2, driver4
        use dispersion_equation, only: ynz
        use trajectory_data

        implicit none
        class(Trajectory), intent(inout) :: traj
        real(wp), intent(inout)    :: pow
        real(wp), intent(in)    :: pabs
        integer,  intent(inout) :: nmax        
        integer,  intent(inout) :: nb1, nb2        
        integer :: num_roots
        !integer,  intent(in)    :: nomth, nomnz


        real(wp) :: xm0
        real(wp) :: tet0
        real(wp) :: xbeg
        integer :: nrefl
        integer :: irep
        integer :: irf, irf1
        integer :: ib2
        integer :: irs0        
        integer :: inak_saved 
        real(wp), parameter :: pgdop=0.02d0
        real(wp), parameter :: hmin=0.d-7 !sav2008, old hmin=1.d-7
        real(wp) :: eps0
        real(wp) :: rrange0, hdrob0, tet, xm, hr
        real(wp) :: xsav, xend,hsav, h1
        real(wp) :: ystart(2),yy(4)

        real(wp) :: xnr
        real(wp) :: ynz0, x1, x2, rexi, tetnew
        real(wp) :: xmnew, rnew, xnrnew
        real(wp) :: xnr_root(4)
        real(wp) :: pg1, pg2, pg3, pg4, pg

        ! copy initial parameters for a trajectory
        xm0  = traj%xmzap
        tet0 = traj%tetzap
        xbeg = traj%rzap
        yn3  = traj%yn3zap
        irs  = traj%irszap
        iw   = traj%iwzap
        izn  = traj%iznzap

        eps0=eps
        rrange0=rrange
        hdrob0=hdrob
        nrefl=0
        im4=0
        nb1=0
        nb2=0
        irep=0
        tet=tet0
        xm=xm0
        hr=1.d0/dble(nr+1) !sav2008
        hrad=hr
        !---------------------------------------
        ! find saving point and define
        ! parameters, depending on direction
        !---------------------------------------
  
  10    irf1=idnint(xbeg/hr)
        if (dabs(irf1*hr-xbeg).lt.tiny1)  then
            xsav=hr*irf1
        else
            irf=int(xbeg/hr)
            if (irs.eq.1)  xsav=hr*irf
            if (irs.eq.-1) xsav=hr*(irf+1)
        end if
        xend= 0.5d0 - 0.5d0*irs + tiny1*irs
        if (ipri.gt.2) write (*,*) 'xbeg-xend',xbeg,xend
        hsav = -hr*irs
        h1 = hsav
        !---------------------------------------
        ! solve eqs. starting from xbeg
        !---------------------------------------
        ystart(1) = tet
        ystart(2) = xm
        call driver2(ystart,xbeg,xend,xsav,hmin,h1, pow, pabs)
        tet = ystart(1)
        xm = ystart(2)
        ib2 = 0

        !---------------------------------------
        ! absorption
        !---------------------------------------
        if(iabsorp.ne.0) then
            if(ipri.gt.2) write (*,*)'in traj() iabsorp=',iabsorp
            nmax=nrefl
            !return
            goto 90
        end if
        if (xend.eq.xbeg) nb1=nb1+1
        !sav2008 20    continue
  
        !--------------------------------------------------------
        !  pass turning point
        !--------------------------------------------------------
        irs0=irs
        num_roots= find_all_roots(xend, xm, tet, xnr_root)
        !num_roots= find_all_roots_simple(xend, xm, tet, xnr_root)

        xnr = xnr_root(1)
        if (num_roots == 0) then
            print *,'no roots'
            !pause
        endif

        ynz0 = ynz
  40    yy(1)=tet 
        yy(2)=xm
        yy(3)=xend
        yy(4)=xnr
        x1=0d0
        x2=1d+10
        rexi=xend
        inak_saved = current_trajectory%size !inak
        call driver4(yy,x1,x2,rexi,hmin, extd4)
        if(iabsorp.eq.-1) goto 90 !return !failed to turn

        tetnew = yy(1)
        xmnew  = yy(2)
        rnew   = yy(3)
        xnrnew = yy(4)

        if(ipri.gt.2) write (*,*) 'from r=',rexi,'to r=',rnew
 
        !---------------------------------------
        ! find mode
        !---------------------------------------

        !call disp2_iroot3(rnew, xmnew, tetnew, xnr_root)
        !print *, xnrnew
        num_roots= find_all_roots(rnew, xmnew, tetnew, xnr_root)
        !num_roots= find_all_roots_simple(rnew, xmnew, tetnew, xnr_root)

        pg1 = abs(xnrnew-xnr_root(1))
        pg2 = abs(xnrnew-xnr_root(2))
        pg3 = abs(xnrnew-xnr_root(3))
        pg4 = abs(xnrnew-xnr_root(4))

        pg = dmin1(pg1,pg2,pg3,pg4)
        if (dabs(pg/xnrnew).gt.pgdop) then
            !---------------------------------------------
            ! bad accuracy, continue with 4 equations
            !--------------------------------------------
            ib2=ib2+1
            nb2=nb2+1
            if (ib2.gt.4) then
                if (ipri.gt.1) write (*,*) 'error: cant leave 4 eqs'
                iabsorp=-1
                print *,'exit ib2.gt.4'
                !return
                goto 90
            end if
            eps=eps/5d0
            rrange=rrange*2d0
            hdrob=hdrob*2d0
            call current_trajectory%reset(inak_saved) ! костыль - восстанавливаю занчение счетчика точек траектории
            goto 40
        end if
        !-------------------------------------
        !          change wave type
        !-------------------------------------
        if (pg.ne.pg1) then
            if (pg.eq.pg2) izn=-izn
            if (pg.eq.pg3) iw=-iw
            if (pg.eq.pg4) iw=-iw
            if (pg.eq.pg4) izn=-izn
        end if
        if (irs0.ne.irs) nrefl=nrefl+1
        xbeg=rnew
        tet=tetnew
        xm=xmnew
        im4=1
        eps=eps0
        rrange=rrange0
        hdrob=hdrob0
        if(nrefl.lt.nmax) goto 10
        rzz=xbeg
        tetzz=tet
        xmzz=xm
        iznzz=izn
        iwzz=iw
        irszz=irs

        !---------------------------------------
        ! remember end point of trajectory
        !---------------------------------------
90      traj%rzap   = rzz
        traj%tetzap = tetzz
        traj%xmzap  = xmzz
        traj%yn3zap = yn3
        traj%iznzap = iznzz
        traj%iwzap  = iwzz
        traj%irszap = irszz
    end  


    

end module trajectory_module