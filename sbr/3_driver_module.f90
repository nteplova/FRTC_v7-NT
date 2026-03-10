module driver_module
    use kind_module
    implicit none
    integer, parameter :: length = 5000000

    !real(wp) dland(length),dcoll(length),perpn(length),dalf(length)
    !real(wp) vel(length),tetai(length)
    !real(wp) xnpar(length)
    !real(wp) :: rho(length)
    !integer izz(length),iww(length),jrad(length)
    !!common/agh/xnpar,vel,dland,dcoll,dalf,perpn,tetai,jrad,iww,izz

    integer     :: irs
    !common /abcd/ irs
    integer     :: iabsorp
    !common /abcdg/ iabsorp

    real(wp)    :: rzz, tetzz, xmzz
    !!common /abc/ rzz,tetzz,xmzz    
    integer     :: iznzz, iwzz, irszz
    !!common /abc/ iznzz,iwzz,irszz

    real(wp)    :: hrad
    !common /bcg/ hrad    

    integer     :: ind
    !common /cmn/ ind    

    integer     :: im4
    !common /bg/ im4    

    !real(wp)    :: pow
    !!common /acg/ pow

    !integer     :: inak, lenstor, lfree
    !common /ag/ inak,lenstor,lfree

    real(wp)    :: pintld4,pintcl4,pintal4
    !common /dg/ pintld4,pintcl4,pintal4

    !real(wp)    ::  vthc(length),poloidn(length)
    !common /vth/ vthc(length),poloidn(length)
contains

    subroutine memorize_trajectory_point4(ro, theta)
        use decrements, only: vfound
        use decrements, only: cf2
        use rt_parameters, only : nr
        implicit none
        real(wp), intent(in) :: ro, theta
        integer jf
        jf=idnint(ro/hrad)
        cf2 = theta
        if(jf.le.0) jf=1
        if(jf.gt.nr) jf=nr
        call memorize_trajectory_point(vfound, jf, ro, 1d0, 4)
    end subroutine

    subroutine memorize_trajectory_point(vz, j, ro, powccc, driver)
        !!  memorize trajectory point
        use rt_parameters, only : nr
        use plasma, only: fvt
        use decrements, only: cf1, cf2, cf3, cf6
        use decrements, only: icf1, icf2
        use dispersion_module, only: ipow
        use decrements, only: pdecv, pdecal
        use trajectory_data

        implicit none
        type(TrajectoryPoint) :: tp
        real(wp), intent(in)     :: vz, ro, powccc
        integer, intent(in)      :: j, driver
        real(wp)    :: radth
        !inak=inak+1
        !if(inak.eq.lenstor) then
            !write(*,*)'storage capacity exceeded !'
            !iabsorp=-1
            !inak=lenstor-1
            !return
        !end if
        tp%vel = vz
        tp%rho =ro
        tp%perpn = cf1 !refr
        tp%poloidn = cf6 !npoloid
        tp%tetai = cf2 ! tet_i
        radth=dble(j)/dble(nr) !!! 
        tp%vthc = 3.d10/fvt(radth)
        tp%iww = icf1 ! было ifast 
        tp%izz = icf2 ! было idir
        tp%xnpar = cf3 !было xparn
        tp%driver = driver
        tp%ipow = ipow
        !print *,'driver=', driver, current_trajectory%size
        if(im4.eq.1) then
            tp%jrad = -j
            tp%dland = pintld4
            tp%dcoll = pintcl4
            tp%dalf  = pintal4
            im4=0
            call current_trajectory%add_point(tp)
            return
        end if
        tp%jrad  = j
        tp%dland = pdecv
        tp%dalf  = pdecal
        if(ipow.ne.1) tp%dcoll = powccc
        if(ipow.eq.1) tp%dcoll = 1d0
        call current_trajectory%add_point(tp)
    end subroutine

    subroutine driver2(ystart,x1,x2,xsav,hmin,h1, pow, pabs) !sav2008
        !! solve eqs. starting from xbeg
        !! ystart(1) = tet
        !! ystart(2) = xm
        !! x1 = xbeg rini
        !! x2 = xend 
        use constants, only: zero, tiny1
        use rt_parameters, only : nr, ipri, rbord, maxstep2, hmin1, iw, eps
        use decrements, only: ifound, vfound
        use dispersion_module, only: ipow,  jfoundr, iconv, irefl, izn, ivar
        use dispersion_equation, only: ynz
        implicit none
        real(wp), intent(inout) :: ystart(2)
        real(wp), intent(in)    :: x1
        real(wp), intent(inout) :: x2, xsav
        real(wp), intent(in)    :: hmin,h1
        real(wp), intent(inout) :: pow
        real(wp), intent(in)    :: pabs
        !common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
        !common /abcd/ irs
        !common /abcdg/ iabsorp
        !common /bcef/ ynz,ynpopq
        !common /bcg/ hrad
        !integer :: ind
        !common /cmn/ ind
        integer, parameter :: nvar=2
        real(wp) :: yscal(nvar),y(nvar),dydx(nvar),yold(nvar),dyold(nvar)
        real(wp) :: x, xold
        real(wp) :: h, hsav, hdid, hnext
        real(wp) :: dstsav, dyd, dst3, dst2, dst1
        real(wp) :: ynz0 !!!!!! проверить значение
        real(wp) :: powccc
        integer  :: i, ii, irep, nstp
        x=x1
        h=dsign(h1,x2-x1)
        ind=0
        ipow=-1
        xold=x
        hsav=hrad*irs
        hdid=zero
        do i=1,nvar
            y(i)=ystart(i)
            yold(i)=y(i)
        end do
        !-------------- start moving -------------
        do nstp=1, maxstep2
            !--------- netpoint control -----------
            dstsav= dabs(x-xsav)
            if(dstsav.lt.tiny1) then
                ipow=ipow+2
                jfoundr=idnint(x/hrad)
                if(jfoundr.le.0) jfoundr=1
                if(jfoundr.gt.nr) jfoundr=nr
            end if
            !print *, 'call extd2 ivar=', ivar
            call extd2(x,y,dydx)
            irep=0
            if(iconv+irefl.ne.0) then
                ! if(iconv+irefl.ne.0.or.ynz.lt.0.d0) then
                !---------------------------------------------
                ! made step to nontransparent zone-return back
                !----------------------------------------------
                x=xold
                do ii=1,nvar
                    y(ii)=yold(ii)
                    dydx(ii)=dyold(ii)
                end do
                irep=1
                h=hdid/2
                hdid=h
                ipow=0
                if(dabs(h).lt.hmin1) then
                    ind=3
                    go to 20
                end if
                ynz=ynz0
                go to 10
            end if            
            !c--------------------------------------
            !c memorize step data
            !c--------------------------------------
            xold=x
            do i=1,nvar
                dyd=dabs(dydx(i))
                yscal(i)=dabs(y(i))+dabs(h*dyd)+1.d-30/(1.d0+dyd)
                yold(i)=y(i)
                dyold(i)=dydx(i)
            end do
            ynz0=ynz
            if(ipow.gt.0) then !integrate power equation
                powccc = dql1(ifound, jfoundr, pow, pabs)
                ! -----------------------------------
                !      memorize trajectory
                ! ----------------------------------
                call memorize_trajectory_point(vfound, jfoundr, x, powccc, 2)
                if(iabsorp.eq.1) then !absorption
                    rzz=x
                    tetzz=y(1)
                    xmzz=y(2)
                    iznzz=izn
                    iwzz=iw
                    irszz=irs
                    return
                end if
                if(iabsorp.eq.-1) return !problem
                ipow=0
                xsav=xsav-hsav
            end if
            !c--------------------------------------
            !c choose step size
            !c--------------------------------------
            dst3=(x-xsav)*(x+h-xsav)
            if(dst3.lt.zero.and.irep.eq.0) h=xsav-x
            if(x.gt.rbord.and.h.gt.zero) then
                ind=2
                go to 20
            end if
10          dst1=(x-rbord)*(x+h-rbord)
            dst2=x*(x+h)
            if((dst1.lt.zero.and.irs.eq.-1).or.dst2.lt.zero) then
                h=h/2.d0
                if(dabs(h).lt.hmin1) then
                    ind=4
                    go to 20
                end if
                go to 10
            end if
            !c--------------------------------------
            !c find solution at x=x+hdid
            !c---------------------------------------
            ynz0=ynz
            call difeq(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,extd2)
20          continue
            if(ind.ne.0) then !exit
                xsav=xsav+hsav
                x2=x
                do i=1,nvar
                    ystart(i)=y(i)
                end do
                ynz=ynz0
                return
            end if
            !c---------------------------------------
            if(dabs(hnext).lt.hmin) then
                if(ipri.gt.1) write(*,*) 'exit driver2: step is too small'
                go to 40
            end if
            h=hnext
        end do
            !c---------------------------------------
        if (ipri.gt.1) write (*,*) 'error in driver2: too many steps'
40      iabsorp=-1
        return
1001    format (10(e14.7,1x))

    end

    subroutine extd2(x, y, dydx)
        use dispersion_module, only: disp2
        implicit none
        real(wp), intent(in)    :: x      ! ro
        real(wp), intent(in)    :: y(:)   ! theta, yn2
        real(wp), intent(inout) :: dydx(:)
        real(wp) :: tt, xm
        real(wp) :: xnr, prt, prm
        tt=y(1)
        xm=y(2)
        !print *,'extd2'
        call disp2(x, xm, tt, xnr, prt, prm)
        !print *,'extd2 after'
        !pause
        dydx(1)=-prm
        dydx(2)=prt
    end

    function dql1(ifound, jfoundr, pow, pabs) result(powccc)
        use constants, only: clt, zero
        use rt_parameters, only: itend0, kv
        use plasma, only: fvt, vperp
        use power,  only: dfind
        use dispersion_module, only: ipow
        use decrements, only: vfound
        use decrements, only: cf1, cf2, cf3, cf4, cf5, cf6
        use decrements, only: icf1,icf2
        use small_vgrid, only: vlf,vrt,dflf,dfrt
        use decrements, only: zatukh ! function zatukh(psy,j,u,n)
        use decrements, only: pdec1,pdec2,pdec3,pdecv,pdecal,dfdv
        implicit none
        integer, intent(in)  :: ifound, jfoundr
        real(wp), intent(in) :: pabs
        real(wp), intent(inout) :: pow
        real(wp)    :: radth
        integer     :: i, j, ifast, idir
        real(wp)    :: powpr,  hdis, vz, refr
        real(wp)    :: dek3, dfsr
        real(wp)    :: vsr, pintld, pintcl, argum, valfa
        real(wp)    :: pintal, dcv, powd, powccc, powcol, powal
        real(wp)    :: pil, pic, pia
        
        powpr=pow
        iabsorp=0
        hdis=hrad
        vz=vfound
        i=ifound
        if(i.eq.0) i=1
        j=jfoundr
        refr=cf1
        !tet_i=cf2
        !npoloid=cf6
        !xparn=cf3

        ifast=icf1
        idir=icf2
        dek3=zero
        dfsr=(vlf*dflf+vrt*dfrt)/2d0*(vrt-vlf)
        vsr=(vrt+vlf)*(vrt-vlf)/2d0
        !c--------------------------------------
        !c   find power
        !c--------------------------------------
        if(im4.eq.1) then
            !!       pintld=-pintld4*dfdv
            pintld=dabs(pintld4*dfdv)
            pintcl=dabs(pintcl4)
            if(itend0.gt.0) then
                argum=clt/(refr*valfa)
                dek3=zatukh(argum,j,vperp,kv)
            end if
            pintal=dabs(pintal4*dek3)
            dcv=pintld4/vsr
        else
            pintld=dabs(pdec1*hdis)
            pintcl=dabs(pdec2*hdis)
            pintal=dabs(pdec3*hdis)
            dcv=pdecv*hdis/vsr
        end if
        if(pabs.ne.zero) then
            powd=pow*dexp(-2d0*pintld)
            powccc=dexp(-2d0*pintcl)
            powcol=powd*powccc
            powal=powcol*dexp(-2d0*pintal)
            pow=powal
        end if
        if(pow.le.pabs) iabsorp=1
        pil=pintld
        pic=pintcl
        pia=pintal
        call dfind(j,i,vz,powpr,pil,pic,pia,dfsr,dcv &
                                ,refr,vlf,vrt,ifast)

    end function

    subroutine difeq(y, dydx, nv,x, htry, eps, yscal, hdid, hnext, derivs)
        use rt_parameters, only : hmin1
        use runge_kutta_module, only: Iderivs_func
        implicit none
        real(wp), intent(inout) :: y(nv)
        real(wp), intent(in)    :: dydx(nv)
        integer,  intent(in)    :: nv
        real(wp), intent(inout) :: x
        real(wp), intent(in)    :: htry
        real(wp), intent(in)    :: eps
        real(wp), intent(inout) :: yscal(nv)
        real(wp), intent(inout) :: hdid
        real(wp), intent(inout) :: hnext
        !external derivs
        procedure(Iderivs_func) :: derivs

        real(wp)            :: dysav(nv) !sav#

        integer, parameter  :: nmax=50,kmaxx=8,imax=kmaxx+1
        real(wp),parameter  :: safe1=.25d0, safe2=.7d0 
        real(wp),parameter  :: redmax=1.d-5, redmin=.7d0
        real(wp),parameter  :: tiny=1.d-30, scalmx=.1d0

        !cu    uses derivs,mmid,pzextr
        integer  :: i,iq,k,kk,km,kmax,kopt,nseq(imax)
        real(wp) :: eps1,epsold,errmax,fact,h,red,scale,work,wrkmin
        real(wp) :: xest, xnew
        real(wp) :: a(imax),alf(kmaxx,kmaxx),err(kmaxx)
        real(wp) :: yerr(nmax),ysav(nmax),yseq(nmax)
        logical  :: first, reduct

        save a,alf,epsold,first,kmax,kopt,nseq,xnew !!! зачем save ?????
        
        real(wp) :: dyd
        integer  :: ii

        !common /cmn/ ind
        data first/.true./,epsold/-1.d0/
        data nseq /2,4,6,8,10,12,14,16,18/

        if(eps.ne.epsold)then
            hnext=-1.d29
            xnew=-1.d29
            eps1=safe1*eps
            a(1)=nseq(1)+1
            do k=1,kmaxx
                a(k+1)=a(k)+nseq(k+1)
            enddo
            do iq=2,kmaxx
                do k=1,iq-1
                    alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+1)))
                enddo
            enddo
            epsold=eps
            do kopt=2,kmaxx-1
                if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt)) goto 1
            enddo
  1         continue
            kmax=kopt
        endif
        h=htry
        do i=1,nv
            ysav(i)=y(i)
            dysav(i)=dydx(i)
        enddo
        if(h.ne.hnext.or.x.ne.xnew)then
            first=.true.
            kopt=kmax
        endif
        reduct=.false.
  2     do k=1,kmax
            xnew=x+h
            if(xnew.eq.x) then
                write(*,*) 'step size underflow in difeq'
               ! pause
            end if
            call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
            !sav#
            if(ind.eq.1) then
                h=h/2d0
                if (dabs(h).lt.hmin1) then
                    do ii=1,nv
                        y(ii)=ysav(ii)
                    end do
                    hnext=h
                    return
                end if
                do ii=1,nv
                        dyd=dabs(dysav(ii))
                        yscal(ii)=dabs(ysav(ii))+dabs(h*dyd)+1.d-30/(1d0+dyd)
                        y(ii)=ysav(ii)
                end do
                goto 2
            end if
            !sav#
            xest=(h/nseq(k))**2
            !var        call pzextr(k,xest,yseq,y,yerr,nv)  !polynomial extrapolation
            call rzextr(k,xest,yseq,y,yerr,nv) !rational extrapolation
            if (k.ne.1) then
                errmax=tiny
                do i=1,nv
                    errmax=max(errmax,abs(yerr(i)/yscal(i)))
                enddo
                errmax=errmax/eps
                km=k-1
                err(km)=(errmax/safe1)**(1.d0/(2*km+1))
            endif
            if (k.ne.1.and.(k.ge.kopt-1.or.first))then
                if (errmax.lt.1.d0) goto 4
                if (k.eq.kmax.or.k.eq.kopt+1) then
                    red=safe2/err(km)
                    goto 3
                else if(k.eq.kopt) then
                if (alf(kopt-1,kopt).lt.err(km)) then
                    red=1.d0/err(km)
                    goto 3
                endif
                else if (kopt.eq.kmax) then
                if (alf(km,kmax-1).lt.err(km)) then
                    red=alf(km,kmax-1)*safe2/err(km)
                    goto 3
                endif
                else if (alf(km,kopt).lt.err(km)) then
                    red=alf(km,kopt-1)/err(km)
                    goto 3
                endif
            endif
        enddo
  3     red=min(red,redmin)
        red=max(red,redmax)
        h=h*red
        reduct=.true.
        goto 2
  4     x=xnew
        hdid=h
        first=.false.
        wrkmin=1.d35
        do kk=1,km
            fact=max(err(kk),scalmx)
            work=fact*a(kk+1)
            if(work.lt.wrkmin)then
                scale=fact
                wrkmin=work
                kopt=kk+1
            endif
        enddo
        hnext=h/scale
        if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
            fact=max(scale/alf(kopt-1,kopt),scalmx)
            if(a(kopt+1)*fact.le.wrkmin)then
                hnext=h/fact
                kopt=kopt+1
            endif
        endif
        return
    end

    subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
        use dispersion_module, only: iconv,irefl
        use runge_kutta_module, only: Iderivs_func
        implicit none
        real(wp), intent(in)    :: y(nvar)
        real(wp), intent(in)    :: dydx(nvar)
        integer,  intent(in)    :: nvar
        real(wp), intent(in)    :: xs
        real(wp), intent(in)    :: htot
        integer,  intent(in)    :: nstep
        real(wp), intent(inout) :: yout(nvar)
        !external derivs
        procedure(Iderivs_func) :: derivs

        integer, parameter :: nmax=50
        integer  :: i,n
        real(wp) :: h,h2,swap,x,ym(nmax),yn(nmax)
        real(wp) :: yz1,yz2
        !integer iconv,irefl
        !common /cefn/ iconv,irefl
        !integer ind
        !common /cmn/ ind
        h=htot/nstep
        yz1=y(1) !sav#
        yz2=y(2) !sav#
        do i=1,nvar
            ym(i)=y(i)
            yn(i)=y(i)+h*dydx(i)
        enddo
        x=xs+h
        call derivs(x,yn,yout)
        if (iconv+irefl.ne.0) goto 10 !sav#
        h2=2.d0*h
        do n=2,nstep
            do i=1,nvar
                swap=ym(i)+h2*yout(i)
                ym(i)=yn(i)
                yn(i)=swap
            enddo
            x=x+h
            call derivs(x,yn,yout)
            if (iconv+irefl.ne.0) goto 10 !sav#
        enddo
        do i=1,nvar
            yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
        enddo
        ind=0 !sav#
        return
  10    ind=1 !sav#
        yout(1)=yz1 !sav#
        yout(2)=yz2 !sav#
        return !sav#
    end    

    subroutine rzextr(iest,xest,yest,yz,dy,nv)
        !! rational extrapolation
        integer,  intent(in)    :: iest, nv
        real(wp), intent(in)    :: xest
        real(wp), intent(in)    :: yest(nv)
        real(wp), intent(inout) :: yz(nv)
        real(wp), intent(inout) :: dy(nv)
        
     
        integer, parameter :: imax=13, nmax=50
        integer  :: j,k
        real(wp) :: b,b1,c,ddy,v,yy
        real(wp) :: d(nmax,imax),fx(imax),x(imax)
        save d,x !! зачем save ????
        x(iest)=xest
        if(iest.eq.1) then
            do j=1,nv
                yz(j)=yest(j)
                d(j,1)=yest(j)
                dy(j)=yest(j)
            enddo
        else
            do k=1,iest-1
                fx(k+1)=x(iest-k)/xest
            enddo
            do j=1,nv
                yy=yest(j)
                v=d(j,1)
                c=yy
                d(j,1)=yy
                do k=2,iest
                    b1=fx(k)*v
                    b=b1-c
                    if(b.ne.0.d0) then
                        b=(c-v)/b
                        ddy=c*b
                        c=b1*b
                    else
                        ddy=v
                    endif
                    if (k.ne.iest) v=d(j,k)
                    d(j,k)=ddy
                    yy=yy+ddy
                enddo
                dy(j)=ddy
                yz(j)=yy
            enddo
        endif
        return
    end

    
    subroutine pzextr(iest,xest,yest,yz,dy,nv)
        !! polynomial extrapolation
        integer iest,nv,imax,nmax
        double precision xest,dy(nv),yest(nv),yz(nv)
        parameter (imax=13,nmax=50)
        integer j,k1
        double precision delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax)
        save qcol,x
        x(iest)=xest
        do j=1,nv
            dy(j)=yest(j)
            yz(j)=yest(j)
        enddo
        if(iest.eq.1) then
            do j=1,nv
                qcol(j,1)=yest(j)
            enddo
        else
            do j=1,nv
                d(j)=yest(j)
            enddo
            do k1=1,iest-1
                delta=1.d0/(x(iest-k1)-xest)
                f1=xest*delta
                f2=x(iest-k1)*delta
                do j=1,nv
                    q=qcol(j,k1)
                    qcol(j,k1)=dy(j)
                    delta=d(j)-q
                    dy(j)=f1*delta
                    d(j)=f2*delta
                    yz(j)=yz(j)+dy(j)
                enddo
            enddo
            do j=1,nv
                qcol(j,iest)=dy(j)
            enddo
        endif
        return
    end

!------------------------------------------------------------------------------------------------
    subroutine driver4(ystart,x1,x2,rexi,hmin, derivs)
        use constants, only : zero
        use runge_kutta_module
        use rt_parameters, only: ipri, eps, hdrob, rbord, maxstep4, rrange
        use dispersion_module, only: idec, ivar, izn
        use dispersion_module, only: pdec14, pdec24, pdec34
        use dispersion_module, only: disp2, disp2_ivar3
        implicit none
        real(wp), intent(inout)  :: ystart(:)
        real(wp), intent(inout)  :: x1,x2
        real(wp), intent(in)     :: rexi, hmin
        procedure (Iderivs_func) :: derivs 
        !external derivs
        !common /abcd/ irs
        !common /abcde/ izn!,iw
        !common /abcdg/ iabsorp
        !common /bdeo/ ivar
        !common /bcef/ ynz,ynpopq
        !common /df/ pdec14,pdec24,pdec34,idec
        !real(wp) pintld4,pintcl4,pintal4
        !common /dg/ pintld4,pintcl4,pintal4
        integer,  parameter :: iturns=1, maxat=3, nvar=4
        real(wp), parameter :: hbeg=1.d-4 !sav2008
        real(wp)  :: x, xnr, dyd, hnext
        !real(wp)  :: prt, prm
        real(wp)  :: yscal(nvar),y(nvar),dydx(nvar),yold(nvar)
        real(wp)  :: eps1, rbord1, hdid, xold, h ,rmm
        real(wp)  :: hdrob1, pdec14zz, pdec24zz, pdec34zz
        integer   :: ipr1, iat, i, ii, nstp
        ipr1=0
        iat=0
        x=zero
        eps1=eps
        hdrob1=hdrob
        rbord1=rbord
        hdid=zero
        pintld4=zero
        pintcl4=zero
        pintal4=zero
        pdec14zz=zero
        pdec24zz=zero
        pdec34zz=zero
        xold=x
        do i=1,nvar
            y(i)=ystart(i)
            yold(i)=y(i)
        end do
        rmm=1d+10*irs
        !sav2008
        !old      rexi1=rexi+rrange
        !old      rexi2=rexi-rrange
        !old      if(rexi1.gt.0.95d0) rexi1=1.d10
        !old      if(rexi2.lt.0.05d0) rexi2=-1.d10
        !est      if(rexi1.gt.0.9d0) rexi1=1.1d0
    10  continue
        !c--------------------------------------
        !c start integration
        !c--------------------------------------
        do nstp=1,maxstep4
            idec=iturns
            call derivs(x,y,dydx)
            idec=0
            pintld4 = pintld4 + abs((pdec14+pdec14zz)/2d0*hdid)
            pintcl4 = pintcl4 + abs((pdec24+pdec24zz)/2d0*hdid)
            pintal4 = pintal4 + abs((pdec34+pdec34zz)/2d0*hdid)
            pdec14zz=pdec14
            pdec24zz=pdec24
            pdec34zz=pdec34
            if (nstp.eq.1) then
                h=hbeg
                !!var        if(dabs(dydx(3)).ne.zero) h=dabs(hmin1/dydx(3))/hdrob1
                if(dabs(dydx(3)).ne.zero) h=0.5d0*dabs(rrange/dydx(3))/hdrob1
            end if
    20      continue
            if(y(3).ge.rbord1.and.dydx(3).gt.zero) then
                !c--------------------------------------
                !c forced reflection from periphery
                !c--------------------------------------
                ivar=3
                izn=-izn
                call disp2_ivar3(y(3),y(2),y(1),xnr)
                if(ivar.eq.-1) then !out of dispersion curve - restart
                    print *, 'out of dispersion curve - restart'
                    do i=1,nvar
                        y(i)=ystart(i)
                    end do
                    x=zero
                    iat=iat+1
                    if(iat.gt.maxat) then
                        if(ipri.gt.1) write (*,*)'turn in driver4 failed'
                        ivar=0
                        goto 40
                    end if
                    eps1=eps1/2.d0
                    hdrob1=hdrob1*2.d0
                    ivar=0
                    goto 10
                end if
                irs=-irs
                y(4)=xnr
                call derivs(x,y,dydx)
                if(dydx(3).gt.zero.and.ipri.gt.1) then
                    write(*,*)'Unsuccesful turn: r, drds=',y(3),dydx(3)
                end if
                ivar=0
                iat=0
            end if
            !sav2008       if((y(3).gt.rexi1.or.y(3).lt.rexi2)) then  ! exit
            !!    if(dabs(y(3)-rexi).gt.rrange.or.nstp.eq.maxstep4) then  ! exit !sav2008
            !call memorize_trajectory_point4(y(3), y(1))
            if(dabs(y(3)-rexi).gt.rrange) then  ! exit !sav2008
                if(dydx(3).gt.zero) irs=-1
                if(dydx(3).lt.zero) irs=1
                if(dydx(3).eq.zero) then !sav2008
                    write(*,*)'exception dr/ds=0 in driver4'
                    pause 'zmi na pedal'
                    go to 1
                end if
                x2=x
                x1=rmm
                do i=1,nvar
                    ystart(i)=y(i)
                end do
                return
            end if
    1       continue
            !c---------------------------------------
            !c remember old values
            !c---------------------------------------
            xold=x
            do i=1,nvar
                dyd=dabs(dydx(i))
                yscal(i)=dabs(y(i))+dabs(h*dyd)+1.d-30/(1.d0+dyd)+1.d-30
                yold(i)=y(i)
            end do
            if (y(3)*irs.lt.rmm*irs) rmm=y(3)
    30      continue
            
            call runge_kutta_qs(y, dydx, nvar, x, h, eps1, yscal, hdid, hnext, derivs)
            if(y(3).ge.1.d0) then  ! crossed plasma boundary
                do ii=1,nvar
                    y(ii)=yold(ii)
                end do
                x=xold
                ipr1=ipr1+1
                if (ipr1.lt.maxat) then
                    h=h/3.d0
                    goto 30
                end if
                rbord1=y(3)-1.d-4
                goto 20
            end if
            ipr1=0
            if(dabs(hnext).lt.hmin) then
                if(ipri.gt.1) write(*,*)'error in dr4: step is too small'
                goto 40
            end if
            h=hnext
        end do
        if(ipri.gt.1) write(*,*)'error in dr4: too many steps.'
        if(ipri.gt.1) write(*,*)'tet=',y(1),'xm=',y(2),'xend=',y(3)
    40  iabsorp=-1
    end    
end module driver_module