module manager_mod
    !! модуль определяет начальные значения лучей и запускает трассировку
    use kind_module
    implicit none

contains

    subroutine manager(iterat,iw0, ntet, spectr)
        use constants            
        use plasma
        use rt_parameters, only : nr, ipri, iw, nmaxm, pabs0, eps, eps_const            
        use trajectory_module
        use spectrum_mod
        use power, only: plost, pnab
        use dispersion_module, only: yn3, ivar, izn, znakstart
        use driver_module, only: irs, iabsorp
        use trajectory_data
        implicit none
        type (Spectrum) spectr
        type (SpectrumPoint) point
        real(wp) pabs
        real(wp) pow
        integer ntet, iout, itr,  nnj,  n_it
        integer maxref, iterat, nmax0, ibad
        !integer inz, itet
        integer iw0, ifail, iabsirp, inak0,ib,ie
        integer nmax, i, nb1,nb2
        real(wp) htet, hr, rin, xmin!, rstart
        real(wp) powexit, dltpow,  pow1, pgamma !, xm
        real(wp) tetin0

        pabs = spectr%max_power*pabs0/1.d2
        !print *, 'pabs =', pabs, spectr%max_power, pabs0
        htet = zero
        hr = 1.d0/dble(nr+1) !sav2008
        if (ntet.ne.1) htet = (tet2-tet1)/(ntet-1)
        irs = 1
        iout = 0
        itr = 0
        nnj = 0
        do n_it = 0,3
            nnj = nnj+nmaxm(n_it+1)
        end do
        maxref = nnj
        if (iterat.lt.3) nmax0=nmaxm(iterat+1)
        if (iterat.ge.3) nmax0=nmaxm(4)
        if (ipri.gt.1) then
            write(*,1001) iterat+1
            write(*,1002)
        end if

        ibad = 0
        itr = 0 
        !--------------------------------------
        ! begin outer loop on teta
        !--------------------------------------
        !do itet = 1,ntet
             !--------------------------------------
            ! begin inner loop on nz
            !--------------------------------------
        !    do inz = 1, spectr%size
        do itr=1, number_of_trajectories
        !        itr = itr+1
                current_trajectory => trajectories(itr)
                znakstart = current_trajectory%znakstart !  znakstart используется в ext4 
                point = spectr%data(current_trajectory%spectrum_point_index)

                if (current_trajectory%mbad.ne.0) then
                    plost = plost+point%power
                    go to 31
                end if

                powexit = point%power
                dltpow = pabs
                call dqliter(dltpow,current_trajectory,hr,powexit,iout)
                if (iout.eq.0) then
                    go to 30
                end if

                pow = powexit 
                ! продолжение траектории 
                ! initial parameters for a trajectory

                nmax = nmax0
                iabsorp = 0
                !-------------------------------------
                ! call ray tracing
                !-------------------------------------
                call tracing(current_trajectory, nmax, nb1, nb2, pow, pabs)
                eps = eps_const 
                !nbad1 = nbad1+nb1
                !nbad2 = nbad2+nb2
                current_trajectory%nrefj = current_trajectory%nrefj + nmax
                powexit = pow
                if (iabsorp.lt.0) then
                    !-------------------------------------
                    !    encounted problems
                    !-------------------------------------

                    if (ipri.gt.1) then
                        tetin0=current_trajectory%tetin   !tet1+htet*(itet-1)
                        write (*,111) tetin0, point%Ntor
111                     format(1x,'traj. with tet0=',f10.5,1x,', Ninput=',f10.5,1x,'failed')
                    end if
                    current_trajectory%mbad = 1 ! плохоая траектория
                    plost= plost+pow
                    goto 30
                end if

20              continue

                if(current_trajectory%nrefj.gt.maxref.and.pow.gt.pabs) then !forced absorp
                    if (current_trajectory%nrefj.gt.maxref+1000) then 
                        print *, 'exceeded limit of reflection number'
                        print *, 'max_refl= ', maxref
                        print *, 'nrefj=', current_trajectory%nrefj
                        print *, 'size=',current_trajectory%size
                        print *, 'pow=', pow, pabs
                        print *, 'point%power=', point%power
                        print *,'------------'
                    endif
                    if(pow.ge.point%power) go to 30 !sav2008
                    pow1 = pow
                    pgamma = 1.d0-pow1/point%power
                    powexit = pow1/pgamma
                    dltpow = powexit-pow1+pabs
                    call dqliter(dltpow, current_trajectory, hr,powexit,iout)
                    powexit = powexit-dltpow+pabs
                    if(powexit.lt.zero) powexit=zero
                end if
30              continue
                pnab = pnab+powexit
31              continue
        !    end do
            !if(ipri.gt.1) write(*,1003)itet,icall1,icall2,current_trajectory%nrefj,nbad1,nbad2
        end do
1001    format (30x,i4,' iteration')
1002    format (6x,'n',5x,'call2',6x,'call4',6x,'nrefl',4x,'last',5x,'bad2',5x,'bad4')
1003    format (3x,i4,2(1x,i10),2x,i7,2x,i8,2(1x,i7),2(2x,i7))
1004    format(1x,i8)
1005    format(1x,i5)
1006    format (e14.7)
    end    

    subroutine dqliter(dltpow, traj, h, powexit, iout) !sav2008
        !! вычисление поглощенной мощности вдоль траектории
        use constants, only: clt, zero
        use rt_parameters, only: itend0, kv
        use small_vgrid, only: vlf, vrt, dflf, dfrt
        use small_vgrid, only: distr
        use decrements, only: pdec1, pdec2, pdec3, pdecv
        use decrements, only: zatukh
        use power,  only: psum4, dfind
        use plasma, only: vperp
        use trajectory_data
        implicit none

        type(Trajectory), pointer, intent(in) :: traj
        real(wp), intent(in)   :: dltpow
        real(wp), intent(in)   :: h
        real(wp), intent(inout)  :: powexit
        integer, intent(inout) :: iout

        type(TrajectoryPoint) :: tp
        integer  :: i, iv,  jr, ifast, jchek
        real(wp) :: pdec1z, pdec3z, pintld, pintal
        real(wp) :: v, refr, dek3, argum, valfa
        real(wp) :: df, dfsr, vsr, pcurr, dcv
        real(wp) :: powpr, powd, powcol, powal
        real(wp) :: pil, pic, pia

        real(wp)    :: pow

        pow=powexit
        pdec1=zero
        pdec1z=zero
        pdec3=zero
        pdec3z=zero
        pdecv=zero
        pintld=zero
        pintal=zero
        iout=0
        do i = 1, traj%size
            !-----------------------------------
            ! restore memorized decrements and
            ! integrate power equation
            !------------------------------------
            tp = traj%points(i)
            v=tp%vel
            jr=tp%jrad
            refr=tp%perpn
            ifast=tp%iww
            dek3=zero
            if(itend0.gt.0) then
                argum=clt/(refr*valfa)
                dek3=zatukh(argum,abs(jr),vperp,kv)
            end if
            !!!old variant
            !!!       call raspr(v,abs(jr),iv,df)
            !!!       if(iv.eq.0) iv=1
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            call distr(v,abs(jr),iv,df)
            !!       dfsr=v*df*(vrt-vlf)
            !!       vsr=v*(vrt-vlf)
            dfsr=(vlf*dflf+vrt*dfrt)/2d0*(vrt-vlf) !sav2008
            vsr=(vrt+vlf)*(vrt-vlf)/2d0 !sav2008
            if(jr.lt.0) then !case of turn
                jr=-jr
                !variant        pintld=-dland(i)*df
                !!        pintld=-dland(i)*(dflf+dfrt)/2d0
                pintld=dabs(tp%dland*(dflf+dfrt)/2d0)
                pdec2=dexp(-2d0*tp%dcoll)
                pintal=dabs(tp%dalf*dek3)
                pcurr=pdec2*dexp(-2d0*pintld-2d0*pintal)
                psum4=psum4+pow*(1d0-pcurr)
                dcv=tp%dland/vsr
            else
                pdec2=tp%dcoll
                pdecv=tp%dland
                !!        pdec1=-pdecv*df
                pdec1=dabs(pdecv*df)
                pdec3=dabs(tp%dalf*dek3)
                pintld=(pdec1+pdec1z)/2d0*h
                pintal=(pdec3+pdec3z)/2d0*h
                pdec1z=pdec1
                pdec3z=pdec3
                dcv=pdecv*h/vsr
            end if
            powpr=pow
            if(dltpow.ne.zero) then
                powd=pow*dexp(-2d0*pintld)
                powcol=powd*pdec2
                powal=powcol*dexp(-2d0*pintal)
                pow=powal
            end if
            pil=pintld
            pic=.5d0*dabs(dlog(pdec2))
            pia=pintal
            call dfind(jr,iv,v,powpr,pil,pic,pia,dfsr,dcv, &
                    refr,vlf,vrt,ifast)
            if(pow.lt.dltpow) then
                powexit=pow
                return
            end if
        end do

    iout=1
    powexit=pow

    end    

end module manager_mod
