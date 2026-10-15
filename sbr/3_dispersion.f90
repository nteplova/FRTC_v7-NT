module metrics
    use kind_module    
    implicit none

    real(wp) :: g11,g12,g22,g33,gg,g,si,co
    !!common/metrika/g11,g12,g22,g33,gg,g,si,co    
    real(wp) :: g2v1, g2jq, g3v
    real(wp) :: b, bp, bt

    real(wp) :: x0, x0t, xj, xjt
    real(wp) :: dxdr, dxdt, dzdr, dzdt
    real(wp) :: dxdrdt, dxdtdt, dzdrdt, dzdtdt
    real(wp) :: bat, bpt,  btt
    real(wp) :: g2jqt
    real(wp) :: g11t, g22t, g33t, g12t, gprt
contains        
    subroutine calculate_metrics(pa,ptet)
        use constants
        use approximation
        use plasma
        use rt_parameters
        use xpoint
        implicit none
        real(wp), intent(in) :: pa
        real(wp), intent(in) :: ptet
        real(wp) :: xdl, xdlp, xly, xlyp
        real(wp) :: xgm, xgmp, xmy, xmyp, xlyv, cotet, sitet


        xdl=fdf(pa,cdl,ncoef,xdlp)
        xly=fdf(pa,cly,ncoef,xlyp)
        xgm=fdf(pa,cgm,ncoef,xgmp)
        xmy=fdf(pa,cmy,ncoef,xmyp)
        xlyv=xlyp*pa+xly
        cotet=dcos(ptet)
        sitet=dsin(ptet)
        !----------------
        dxdr=-xdlp+cotet-xgmp*sitet**2
        dxdt=-(pa+two*xgm*cotet)*sitet
        dzdr=xlyv*sitet
        dzdt=xly*pa*cotet
        x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
        dxdrdt=-sitet-two*xgmp*sitet*cotet
        dzdrdt=xlyv*cotet
        dxdtdt=-pa*cotet-two*xgm*(cotet**2-sitet**2)
        dzdtdt=-xly*pa*sitet
        x0t=dxdt
        if(-5.gt.0) then !!!!!!!!!!!xtochk
            call xpnt(pa,ptet-1.5d0*pi,dxdr,dxdt,dzdr,dzdt,x0,dxdrdt,dxdtdt,dzdrdt,dzdtdt,x0t)
        end if
        !--------------------------------------
        ! components of metric tensor
        !--------------------------------------
        g11=dxdr**2+dzdr**2
        g22=dxdt**2+dzdt**2
        g12=dxdr*dxdt+dzdr*dzdt
        g33=x0**2
        xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
        gg=xj
        g=xj*g33
        g2v1=one/dsqrt(g22)
        g2jq=dsqrt(g22/xj)
        g3v=one/dsqrt(g33)
        !--------------------------------------
        !  magnetic field
        !--------------------------------------
        bt=b_tor*(r0/rm)/x0
        bp=g2jq*g3v*xmy ! dsqrt(g22/xj)/dsqrt(g33)
        b=dsqrt(bp*bp+bt*bt)
        si=bp/b
        co=bt/b

        !------ calculation of derivatives ------

        g11t=two*(dxdr*dxdrdt+dzdr*dzdrdt)
        g22t=two*(dxdt*dxdtdt+dzdt*dzdtdt)
        g33t=two*x0*(-pa*sitet-two*xgm*sitet*cotet)
        g12t=dxdrdt*dxdt+dxdr*dxdtdt+dzdrdt*dzdt+dzdr*dzdtdt
        xjt=g11t*g22+g22t*g11-two*g12*g12t
        btt=-b_tor*(r0/rm)/x0**2*x0t
        g2jqt=(g22t/xj-g22/xj**2*xjt)/(g2jq*two)
        bpt=xmy*(g2jqt*g3v-.5d0*g2jq*g3v/g33*g33t)
        bat=one/b*(bp*bpt+bt*btt)

    end subroutine calculate_metrics
end module metrics

module dielectric_tensor
    !! components of dielectric tensor
    use kind_module    
    use metrics
    implicit none
    real(wp) :: e1, e2, e3
    real(wp) :: pn, fnr, wpq, whe, v, u1, u
contains
    subroutine calculate_dielectric_tensor(pa)
        !! calculate components of dielectric tensor
        
        use constants, only: zero, one, two
        use constants, only: c0, c1, pi        
        use rt_parameters, only: inew
        use plasma, only: fn1, fn2
        use plasma, only: ww, xmi
        implicit none
        real(wp), intent(in) :: pa

        real(wp) :: fnrr
        pn = 0
        fnr = 0
        fnrr = 0
        if(inew.eq.0) then !vardens
            pn=fn1(pa,fnr)
        else
            pn=fn2(pa,fnr,fnrr)
        end if
        wpq=c0**2*pn ! 
        whe=b*c1     !
        v=wpq/ww**2  ! квадрат плазменной частоты на частоту волны
        u1=whe/ww    ! 
        u=u1**2
        
        e1=one-v*(one/xmi-one/u)
        e2=v/u1
        e3=one-v

    end subroutine
end module dielectric_tensor


module dispersion_equation
    !! dispersion equation
    use kind_module    
    use metrics
    implicit none
    real(wp) :: ynz, ynzq
    real(wp) :: as, bs, cs
    real(wp) :: pnew, yny, gpr, dls
contains
    subroutine calculate_dispersion_equation(yn2 , yn3)
        use constants, only: zero, one, two
        use constants, only: c0, c1
        use dielectric_tensor
        use rt_parameters, only: inew
        use plasma, only: ww, xsz
        implicit none
        real(wp), intent(in) :: yn2 , yn3
        !real(wp) :: 
       
        ! dispersion equation
        !--------------------------------------
        !sav2008      if(ivar.eq.2) yn2=(ynz-yn3*co*g3v)/(si*g2v1)
        ynz=yn2*si*g2v1+yn3*co*g3v  
        ynzq=ynz**2
        as=e1
        bs=-(e1**2-e2**2+e1*e3-(e1+e3)*ynzq)
        cs=e3*(e1**2-e2**2-two*e1*ynzq+ynzq**2)
        !-----------------------------
        !est !sav2009
        pnew=zero
        yny= - (yn2*g2v1*co-yn3*g3v*si)
        if(inew.gt.0) then
            if(inew.eq.1) then
                    yny= - (yn2*g2v1*co-yn3*g3v*si)
                elseif(inew.eq.2) then
                    yny= - g2jq*(yn2*g2v1*co-yn3*g3v*si)
            end if
            gpr=c0**2/ww**2/u1*fnr*xsz
            pnew=yny*gpr
            bs=bs+pnew
            cs=cs+pnew*(ynzq-e3)
        end if
        !------------------------------------
        dls=bs*bs-4d0*as*cs


        !c      write(*,*)'rho=',pa,' teta=',ptet
        !c      write(*,*)'N2=',yn2,' N3=',yn3
        !c      write(*,*)'Npar=',ynz,' e1=',e1
        !c      write(*,*)'v=',v,' u=',u
        !c      write(*,*)'whe=',whe,' ww=',ww
        !c      write(*,*)'e2=',e2,' e3=',e3
        !c      write(*,*)'bs=',bs,' as=',as
        !c      write(*,*)'cs=',cs,' dls=',dls
        !c      pause
    end subroutine
end module dispersion_equation


module partial_derivatives
    !! calculation of partial derivatives
    use kind_module    
    
    implicit none
contains
    subroutine calculate_partial_derivatives(yn2 , yn3, dl1, ynpopq, al, bl, dl2, prt, prm)
        use constants, only: zero, one, two
        use constants, only: c0, c1
        use metrics
        use dielectric_tensor
        use dispersion_equation
        use rt_parameters, only: inew
        use plasma, only: ww, xsz
        implicit none
        real(wp), intent(in) :: yn2 , yn3, dl1, ynpopq, al, bl, dl2
        real(wp), intent(out) :: prt, prm

        real(wp) :: s1, p1, p2, p3, ynzt, e1t, e2t, u1t, cot, sit
        real(wp) :: s2, dnm, v1, v2, vvt, vvm
        real(wp) :: ynyt, dnym
        real(wp) :: pnewt
        real(wp) :: s21, sjg, s23, s24, s22
        !--------------------------------------
        !   calculation of derivatives
        !--------------------------------------

        !g11t=two*(dxdr*dxdrdt+dzdr*dzdrdt)
        !g22t=two*(dxdt*dxdtdt+dzdt*dzdtdt)
        !g33t=two*x0*(-pa*sitet-two*xgm*sitet*cotet)
        !g12t=dxdrdt*dxdt+dxdr*dxdtdt+dzdrdt*dzdt+dzdr*dzdtdt
        !xjt=g11t*g22+g22t*g11-two*g12*g12t
        !btt=-b_tor*(r0/rm)/x0**2*x0t
        !g2jqt=(g22t/xj-g22/xj**2*xjt)/(g2jq*two)
        !bpt=xmy*(g2jqt*g3v-.5d0*g2jq*g3v/g33*g33t)
        !bat=one/b*(bp*bpt+bt*btt)
        !-----
        sit=bpt/b-bp/b**2*bat
        cot=btt/b-bt/b**2*bat
        u1t=c1/ww*bat
        e1t=-v/u**2*two*u1*u1t
        e2t=-v/u*u1t
        ynzt=yn2*sit*g2v1-yn2*si*g2v1**3/two*g22t + yn3*cot*g3v-yn3*co*g3v**3/two*g33t
        p1=two*ynz*ynzt
        p2=e1t
        p3=(e2*e2t)/e1-e2**2/(two*e1**2)*e1t

        s1=-p2/(two*e1**2)*e3*(ynzq-e1) + (e3+e1)/(two*e1)*(p1-p2)+p3
        s2=two*e3/e1*(ynzq-e1)*(p1-p2) - p2*e3/(e1**2)*(ynzq-e1)**2-two*e3*p3
        dnm=two*ynz*si*g2v1
        v1=(e3+e1)/(two*e1)*dnm
        v2=two*e3/e1*(ynzq-e1)*dnm
        !-----------------------------------
        !est !sav2009
        if(inew.gt.0) then
            gprt=-c0**2/ww**2*fnr/u*u1t*xsz
            if(inew.eq.1) then
                ynyt= - (yn2*cot*g2v1-yn2*co*g2v1**3/two*g22t - (yn3*sit*g3v-yn3*si*g3v**3/two*g33t))
                dnym= - co*g2v1
            else if(inew.eq.2) then
                ynyt= - g2jq*(yn2*cot*g2v1-yn2*co*g2v1**3/two*g22t - (yn3*sit*g3v-yn3*si*g3v**3/two*g33t))
                ynyt=ynyt - g2jqt*(yn2*g2v1*co-yn3*g3v*si)
                dnym= - g2jq*co*g2v1
            end if
            pnewt=(ynyt*gpr+yny*gprt)
            s1=s1+pnewt/(two*e1)-pnew/(two*e1**2)*e1t
            s2=s2+(pnewt*(ynzq-e3)+pnew*p1)/e1-pnew*(ynzq-e3)/e1**2*e1t
            v1=v1+dnym*gpr/(two*e1)
            v2=v2+gpr*(dnym*(ynzq-e3)+yny*dnm)/e1
        end if
        !---------------------------------------------
        vvt=-s1+(bs/as*s1-s2)/(two*dl1)
        vvm=-v1+(bs/as*v1-v2)/(two*dl1)
        s1=-yn2*(g12t/g22-g12/g22**2*g22t)
        s21=yn2**2*(g11t/g22-g11/g22**2*g22t)
        s22=yn3**2*( xjt/(g33*g22)-xj/(g33*g22)**2*(g33t*g22+g22t*g33) )
        sjg=(xjt*g22-xj*g22t)/g22**2
        s23=two*ynz*ynzt*xj/g22+sjg*ynzq
        s24=vvt*xj/g22+ynpopq*sjg
        s2=s21+s22-s23-s24
        prt=-s1+(two*(bl/al)*s1-s2)/(two*dl2)
        s1=-g12/g22
        s21=two*yn2*g11/g22
        s22=dnm*xj/g22
        s23=vvm*xj/g22
        s2=s21-s22-s23
        prm=-s1+(two*(bl/al)*s1-s2)/(two*dl2)

    end subroutine
end module partial_derivatives

module decrements
    use kind_module
    use nr_grid, only: MAX_NR
    implicit none
    real(wp) :: dnx
    real(wp) :: dhdnr

    integer  :: ifound
    real(wp) :: vfound
    !!common /eg1/ vfound,ifound
    integer  :: icf1,icf2
    real(wp) :: pdec1,pdec2,pdec3,pdecv,pdecal,dfdv
    !!common /eg2/ pdec1,pdec2,pdec3,pdecv,pdecal,dfdv,icf1,icf2
    real(wp) :: cf1,cf2,cf3,cf4,cf5,cf6
    !!common /eg3/ cf1,cf2,cf3,cf4,cf5,cf6

    real(wp) :: dgdu(50,MAX_NR)
    integer  :: kzero(MAX_NR)
    !!common /arr/ dgdu(50,100),kzero(100)
    !! используется в zatukh, ourlhcd2017 и alphas

contains
    subroutine calculate_decrements(pa, yn2, ptet, yn3, ynpopq,  xnr, jr, izn)
        use constants, only: zero, one, two
        use constants, only: c0, c1, pi
        use constants, only: zalfa, xmalfa, xlog, clt
        use plasma, only: fn1, fn2, fvt, ft, zefff
        use plasma, only: ww, xmi,xsz, cltn, cnye, cnyi
        use plasma, only: cnstal, valfa, vperp
        use rt_parameters, only: inew, iw, itend0, kv        
        use metrics
        use dielectric_tensor
        use dispersion_equation, only : as, bs, yny, ynz, ynzq
        use small_vgrid, only: distr
        use source_new_mod, only: source_new
        implicit none
        real(wp), intent(in) :: pa, yn2, ptet, yn3, ynpopq, xnr
        integer, intent(in)  :: jr, izn


        real(wp) :: sl1, vz, vt, fder, aimh, pnye, pnyi
        real(wp) :: tmp, fcoll, source, argum
        real(wp) :: dek1, dek2, dek3
        
        !--------------------------------------
        !  calculation of decrements
        !--------------------------------------
        dnx=two*as*ynpopq+bs
        dhdnr=dnx*(two*g22*xnr-two*g12*yn2)/xj
        sl1=(ynzq-e1)*(ynzq+ynpopq-e1)-e2**2
        cf3=ynz
        cf4=xnr
        cf5=yn2
        vz=cltn/dabs(ynz)
        !if(vz.gt.cltn) vz=cltn !sav2010
        vt=fvt(pa)
        !jr=jfoundr
        icf1=iw
        icf2=izn
        call distr(vz,jr,ifound,fder)
        dfdv=fder
        vfound=vz
        cf2=ptet
        cf6=yny
        aimh=wpq/ww**2*pi*sl1*cltn**2/ynzq
        pdecv=dabs(aimh/dhdnr/xsz)
        !!        pdec1=-pdecv*dfdv
        pdec1=dabs(pdecv*dfdv)
        pnye=cnye*wpq**2/(pn*vt**3)
        pnyi=cnyi*pnye*zefff(pa)
        pdec2=dabs(pnyi/ww*(wpq/whe**2*ynpopq+wpq/ww**2*ynzq)*ynpopq/dhdnr/xsz)
        cf1=dsqrt(ynpopq)
        if(itend0.gt.0) then
            tmp=ft(pa)/0.16d-8
            fcoll=.5d-13*pn*zalfa**2*xlog/xmalfa/tmp**1.5d0
            !cc          ddens=dn1*pn
            !cc          tdens=dn2*pn
            !cc          tt=fti(pa)**0.33333d0    ! (ti, kev)^1/3
            !cc          source=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
            call source_new(pa,source)
            dek1=cnstal*pdecv*(1.d0-e3/ynpopq)**2/cf1
            dek2=source/(fcoll*pn)
            pdecal=dek1*dek2
            pdec3=zero
            if(itend0.gt.0) then
                argum=clt/(cf1*valfa)
                dek3=zatukh(argum,jr,vperp,kv)
                pdec3=pdecal*dek3
            end if
        end if
    end subroutine

    real(wp) function zatukh(psy,j,u,n)
        use constants
        !implicit real*8 (a-h,o-z)
        implicit none
        real(wp), intent(in) :: psy
        real(wp), intent(in) :: u(:,:)
        integer,  intent(in) :: j, n
        real(wp) :: x(50),y(50),a(50),b(50)
        real(wp) :: um, s1, s2, ss1, ss2, sum
        !common /a0befr/ pi,pi2
        !common /arr/ dgdu(50,100),kzero(100)
        integer :: km, k, i, l
        km = kzero(j)
        um = u(km,j)
        if(um.ge.one) then
            zatukh = zero
            if(psy.lt.one) zatukh = .5d0*pi/psy**3
            return
        end if
        if(psy-um.le.zero.or.u(n,j)-psy.le.zero) then
            zatukh = zero
            return
        end if
        do k=1,n
            x(k) = u(k,j)
            y(k) = dgdu(k,j)
        end do
        i=n-1
        do l=1,n-1
            if (x(l+1)-psy.gt.zero.and.psy-x(l).ge.zero) i=l
        end do
        do k=i,n-1
            b(k) = (y(k+1)-y(k))/(x(k+1)-x(k))
            a(k) = y(k) - b(k)*x(k)
        end do
        s2  = sqrt((x(i+1)-psy)*(x(i+1)+psy))
        ss2 = x(i+1) + s2
        sum = a(i)*log(psy/ss2) - b(i)*s2
        do k=2,n-i
            s1  = sqrt((x(i+k-1)-psy)*(x(i+k-1)+psy))
            ss1 = x(i+k-1) + s1
            s2  = sqrt((x(i+k)-psy)*(x(i+k)+psy))
            ss2 = x(i+k) + s2
            sum = sum + a(i+k-1)*log(ss1/ss2) + b(i+k-1)*(s1-s2)
        end do
        zatukh=sum
        return
    end   

end module decrements

module dispersion_module
    use kind_module    
    implicit none

    real(wp) :: yn3
    !!common /abefo/ yn3

    integer :: ivar
    !!common /bdeo/ ivar   

    integer  :: izn
    !!common /abcde/ izn

    real(wp) :: ynpopq
    !real(wp) ::ynz, ynpopq
    !!common /bcef/ ynz,ynpopq

    integer iconv, irefl
    !!common /cefn/ iconv,irefl

    integer ipow, jfoundr
    !!common /ceg/ ipow,jfoundr

    real(wp) :: dhdm,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
    !!common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u

    real(wp) :: znakstart
    !!common/direct/znakstart

    real(wp) :: ham
    !!common/fjham/ham

    integer  :: idec
    real(wp) :: pdec14,pdec24,pdec34
    !!common /df/ pdec14,pdec24,pdec34,idec

contains


    subroutine disp2_ivar3(pa,yn2,ptet,xnro)
        ! case iroot == 1 ivar= 3
        use constants, only: zero, one, two
        use rt_parameters, only: iw
        use metrics
        use dielectric_tensor
        use dispersion_equation
        use partial_derivatives
        use decrements
        implicit none
        real(wp), intent(in) :: pa      ! ro
        real(wp), intent(in) :: yn2     ! ???
        real(wp), intent(in) :: ptet    ! theta
        real(wp), intent(out) :: xnro ! ???

        integer  :: jr

        real(wp) :: dl1, ynpopq1, al, bl, cl, cl1, dll
        real(wp) :: dl2, xnr 

        iconv=0
        irefl=0
        if(pa.ge.one.or.pa.le.zero) then
            print *, 'pa=',pa
            pause
            !goto 70
        endif
          
        call calculate_metrics(pa, ptet)

        call calculate_dielectric_tensor(pa)

        call calculate_dispersion_equation(yn2 , yn3)
        
        if(dls.lt.zero) then
            ! conversion
            print *, 'conversion'
            !pause
            iconv=1
            ivar=-1
            return
        end if
300     continue
        dl1=dfloat(iw)*dsqrt(dls)/two/as
        if(iw.eq.-1) ynpopq=-bs/(two*as)+dl1        ! = (-bs + sqrt(dls)) / (2*as)
        if(iw.eq.1)  ynpopq=two*cs/(-bs-two*as*dl1) ! = (-bs - sqrt(dls)) * (2*cs)??
        

        !cc      write(*,*)'iw=',iw,' izn=',izn,' Nperp=',dsqrt(ynpopq)
        !cc      write(*,*)'Nperp2=',ynpopq,' ynpopq1=',-bs/(two*as)-dl1
        !cc      pause

        !if (ynpopq.lt.zero) goto 70
        al=g22/xj
        bl=-yn2*g12/xj
        cl=g11*yn2**2/xj+yn3**2/g33-ynzq-ynpopq

        dll=bl*bl-al*cl

        !if ((dll.lt.zero).or.(ynpopq.lt.zero)) then
        if (dll.lt.zero) then
        !  70 reflection 
            irefl=1
            if (ivar.eq.10) then 
                ivar=-1
                return
            end if
            !print *, 'reflection'
            !pause
            iw=-iw
            ivar=10
            goto 300
            

        endif

400     dl2=-dfloat(izn)*dsqrt(dll)/al
        if(izn.eq.1) xnr=-bl/al+dl2         ! =  (-bl - sqrt(dll)) / al
        if(izn.eq.-1) xnr=cl/(-bl-al*dl2)   ! =  (-bl + sqrt(dll)) / al
        xnro=xnr
 
        !  find Nr of reflected wave
        dnx=two*as*ynpopq+bs
        dhdnr=dnx*(two*g22*xnr-two*g12*yn2)/xj
        if(-znakstart*dhdnr.gt.zero) then
            izn=-izn
            goto 400
        end if

        return


    end

    subroutine disp2(pa,yn2,ptet,xnro,prt,prm)
        ! case iroot == 1 ivar= 0 
        use constants, only: zero, one, two
        use rt_parameters, only: iw
        use metrics
        use dielectric_tensor
        use dispersion_equation
        use partial_derivatives
        use decrements
        implicit none
        real(wp), intent(in) :: pa      ! ro
        real(wp), intent(in) :: yn2     ! ???
        real(wp), intent(in) :: ptet    ! theta
        real(wp), intent(out) :: xnro ! ???
        real(wp), intent(out) :: prt  ! ???
        real(wp), intent(out) :: prm  ! ???       

        integer  :: jr

        real(wp) :: dl1, ynpopq1, al, bl, cl, cl1, dll
        real(wp) :: dl2, xnr 

       
        if (ivar/=0) then
            print *, 'disp2 ivar=', ivar
            print *, "disp2"
            pause
        endif
        iconv=0
        irefl=0
        if(pa.ge.one.or.pa.le.zero) goto 70

        call calculate_metrics(pa, ptet)

        call calculate_dielectric_tensor(pa)

        call calculate_dispersion_equation(yn2 , yn3)
        
        if(dls.lt.zero) then
            ! conversion
            iconv=1
            return
        end if
30      continue
        dl1=dfloat(iw)*dsqrt(dls)/two/as
        if(iw.eq.-1) ynpopq=-bs/(two*as)+dl1        ! = (-bs + sqrt(dls)) / (2*as)
        if(iw.eq.1)  ynpopq=two*cs/(-bs-two*as*dl1) ! = (-bs - sqrt(dls)) * (2*cs)??
        

        !cc      write(*,*)'iw=',iw,' izn=',izn,' Nperp=',dsqrt(ynpopq)
        !cc      write(*,*)'Nperp2=',ynpopq,' ynpopq1=',-bs/(two*as)-dl1
        !cc      pause

        !if (ynpopq.lt.zero) goto 70
        al=g22/xj
        bl=-yn2*g12/xj
        cl=g11*yn2**2/xj+yn3**2/g33-ynzq-ynpopq

        dll=bl*bl-al*cl

        if(dll.lt.zero) goto 70

40      dl2=-dfloat(izn)*dsqrt(dll)/al
        if(izn.eq.1) xnr=-bl/al+dl2         ! =  (-bl - sqrt(dll)) / al
        if(izn.eq.-1) xnr=cl/(-bl-al*dl2)   ! =  (-bl + sqrt(dll)) / al
        xnro=xnr

        !--------------------------------------
        !   calculation of derivatives
        !--------------------------------------
        call calculate_partial_derivatives(yn2, yn3, dl1, ynpopq, al, bl, dl2, prt, prm)

        if(ipow.gt.0) then
            !--------------------------------------
            !  calculation of decrements
            !--------------------------------------
            call calculate_decrements(pa, yn2, ptet, yn3, ynpopq, xnr, jfoundr, izn)

        end if
        return

        !  reflection
70      irefl=1
        return
    end

    function find_all_roots_simple(pa, yn2, ptet, root) result(num_roots)
        !! find all roots of dispersion equation
        use constants, only: zero, one, two
        use rt_parameters, only: iw
        use metrics
        use dielectric_tensor
        use dispersion_equation
        implicit none
        real(wp), intent(in) :: pa      ! ro
        real(wp), intent(in) :: yn2     ! ???
        real(wp), intent(in) :: ptet    ! theta
        real(wp), intent(inout) :: root(:)
    
        integer  :: num_roots
        real(wp) :: ynpopq1, al, bl, cl, cl1
        real(wp) :: dll1, dll
        real(wp) :: dl1,dl2
        real(wp) :: xnr

        iconv=0
        irefl=0
        if(pa.ge.one.or.pa.le.zero) then
            print *, ' find_all_roots_simple'
            pause
            return
        endif

        
        call calculate_metrics(pa, ptet)

        call calculate_dielectric_tensor(pa)

        call calculate_dispersion_equation(yn2 , yn3)

        num_roots = 0
        root(:)=1d+10

        if(dls.lt.zero) return

        ynpopq = (-bs + iw*sqrt(dls))/(2*as)  ! = - two*cs/ (bs + iw*sqrt(dls))
        ynpopq1= (-bs - iw*sqrt(dls))/(2*as)

        al=g22/xj
        bl=-yn2*g12/xj
        cl= g11*yn2**2/xj + yn3**2/g33 - ynzq - ynpopq
        dll=bl*bl-al*cl

        if (dll.ge.zero) then
            root(1)= (-bl - izn*sqrt(dll))/al
            root(2)= (-bl + izn*sqrt(dll))/al ! = cl/(-bl - sqrt(dll))
            num_roots = num_roots + 2
        end if
        
        cl1  = g11*yn2**2/xj + yn3**2/g33 - ynzq - ynpopq1
        dll1 = bl**2 - al*cl1

        if (dll1.ge.zero) then
            root(3)= (-bl - izn*sqrt(dll1))/al
            root(4)= (-bl + izn*sqrt(dll1))/al
            num_roots = num_roots + 2
        end if
    end

    function square_solver(a,b,c, root) result(num_roots)
        implicit none
        real(wp), intent(in) :: a
        real(wp), intent(in) :: b
        real(wp), intent(in) :: c
        real(wp), intent(inout) :: root(:)      
        integer  :: num_roots
        real(wp) :: epsilon, dec

        num_roots = 0

        epsilon = 1
        if (b<0) epsilon = -1

        dec = b**2 - 4*a*c
        if (dec>=0) then
            root(1) =  -(b + epsilon*sqrt(dec))/(2.0*a)        
            root(2) = -2.0*c/ (b + epsilon*sqrt(dec)) ! = (-bs + iw*sqrt(dls))/(2*as)  
            num_roots = num_roots + 2
        endif    
    end function

    function square_solver2(a,b,c, root) result(num_roots)
        implicit none
        real(wp), intent(in) :: a
        real(wp), intent(in) :: b
        real(wp), intent(in) :: c
        real(wp), intent(inout) :: root(:)      
        integer  :: num_roots, r1, r2
        real(wp) :: epsilon, dec

        num_roots = 0

        r1 = 1
        r2 = 2
        epsilon = 1
        if (b<0) then 
            epsilon = -1
            r1 = 2
            r2 = 1
        endif

        dec = b**2 - 4*a*c
        if (dec>=0) then
            root(r1) =  -(b + epsilon*sqrt(dec))/(2.0*a)        
            root(r2) = -2.0*c/ (b + epsilon*sqrt(dec)) ! = (-bs + iw*sqrt(dls))/(2*as)  
            num_roots = num_roots + 2
        endif    
    end function



    function find_all_roots(pa, yn2, ptet, root) result(num_roots)
        !! find all roots of dispersion equation
        !! уточненный вариант вычисления корней
        use constants, only: zero, one, two
        use rt_parameters, only: iw
        use metrics
        use dielectric_tensor
        use dispersion_equation
        implicit none
        real(wp), intent(in) :: pa      ! ro
        real(wp), intent(in) :: yn2     ! ???
        real(wp), intent(in) :: ptet    ! theta
        real(wp), intent(inout) :: root(:)
    
        integer  :: num_roots

        real(wp) :: al, bl, cl, cl1, ynpopq1
        real(wp) :: ynpopq_(2), rt(2)
        real(wp) :: tmp(4)
        real(wp) :: dll1, dll
        real(wp) :: dl1,dl2
        real(wp) :: xnr
        real(wp) :: epsilon
        real(wp) :: err_1, err_2
        iconv=0
        irefl=0
        if(pa.ge.one.or.pa.le.zero) then
            print *, 'find_all_roots'
            pause
            return
        endif

        
        call calculate_metrics(pa, ptet)

        call calculate_dielectric_tensor(pa)

        call calculate_dispersion_equation(yn2 , yn3)

        num_roots = 0
        root(:)=1d+10

        if (square_solver(as, bs, cs, ynpopq_)==0) return
        
        if (iw<0) then
            ynpopq  = ynpopq_(2) 
            ynpopq1 = ynpopq_(1) 
        else
            ynpopq1 = ynpopq_(2) 
            ynpopq =  ynpopq_(1) 
        endif

        al=g22/xj
        bl=-yn2*g12/xj
        cl= g11*yn2**2/xj + yn3**2/g33 - ynzq 

        if (square_solver2(al, 2*bl, cl - ynpopq, rt)>0) then
            if (izn == 1) then
                root(1) = rt(1)
                root(2) = rt(2)
            else 
                root(1) = rt(2)
                root(2) = rt(1)
            endif
            num_roots = num_roots + 2
        end if

        if (square_solver2(al, 2*bl, cl- ynpopq1, rt)>0) then
            if (izn == 1) then
                root(3)= rt(1) !(-bl - sqrt(dll))/al
                root(4)= rt(2) !cl/(-bl - sqrt(dll))
            else 
                root(3)= rt(2) !cl/(-bl - sqrt(dll))
                root(4)= rt(1) !(-bl + izn*sqrt(dll))/al ! = cl/(-bl - sqrt(dll))
            endif
        end if
    end




    function find_all_roots_test(pa, yn2, ptet, root) result(num_roots)
        !! find all roots of dispersion equation
        !! уточненный вариант вычисления корней
        use constants, only: zero, one, two
        use rt_parameters, only: iw
        use metrics
        use dielectric_tensor
        use dispersion_equation
        implicit none
        real(wp), intent(in) :: pa      ! ro
        real(wp), intent(in) :: yn2     ! ???
        real(wp), intent(in) :: ptet    ! theta
        real(wp), intent(inout) :: root(:)
    
        integer  :: num_roots

        real(wp) :: al, bl, cl, cl1, ynpopq1
        real(wp) :: ynpopq_(2), rt(2)
        real(wp) :: tmp(4)
        real(wp) :: dll1, dll
        real(wp) :: dl1,dl2
        real(wp) :: xnr
        real(wp) :: epsilon
        real(wp) :: err_1, err_2
        iconv=0
        irefl=0
        if(pa.ge.one.or.pa.le.zero) then
            print *, ' ivar=', ivar
            pause
            return
        endif

        call calculate_metrics(pa, ptet)

        call calculate_dielectric_tensor(pa)

        call calculate_dispersion_equation(yn2 , yn3)

        num_roots = 0
        root(:)=1d+10

        !if(dls.lt.zero) return
        if (square_solver(as, bs, cs, ynpopq_)==0) return

        epsilon = 1
        if (bs<0) epsilon = -1
        
        if (iw<0) then
            ynpopq  = ynpopq_(2) !-two*cs/ (bs + epsilon*sqrt(dls)) ! = (-bs + iw*sqrt(dls))/(2*as)  
            !ynpopq = (-bs + sqrt(dls))/(2.0*as)  ! = - two*cs/ (bs + iw*sqrt(dls))
            ynpopq1 = ynpopq_(1) !-(bs + epsilon*sqrt(dls))/(2.0*as)
        else
            ynpopq1 = ynpopq_(2) !-two*cs/ (bs + epsilon*sqrt(dls)) ! = (-bs + iw*sqrt(dls))/(2*as)  
            ynpopq =  ynpopq_(1) !-(bs + epsilon*sqrt(dls))/(2.0*as)
            !ynpopq = two*cs/ (-bs - iw*sqrt(dls)) ! = (-bs + iw*sqrt(dls))/(2*as)  
            !ynpopq1= two*cs/ (-bs + iw*sqrt(dls)) ! = (-bs - iw*sqrt(dls))/(2*as)
        endif

        print *, 'ynpopq', '   ynpopq1'
        !print *, ynpopq_
        print *, ynpopq, ynpopq1
        print *, two*cs/ (-bs - iw*sqrt(dls)) ! = (-bs + iw*sqrt(dls))/(2*as)  
        print *, two*cs/ (-bs + iw*sqrt(dls)) ! = (-bs - iw*sqrt(dls))/(2*as)
        !print *, ynpopq, ynpopq1
        !print *, ' --------------'
        !print *, as, bs + sqrt(dls)
        !print *, bs, cs
        

        al=g22/xj
        bl=-yn2*g12/xj
        cl= g11*yn2**2/xj + yn3**2/g33 - ynzq - ynpopq
        dll=bl*bl-al*cl
        print *,'dll=', dll
        if (dll.ge.zero) then
            if (izn == 1) then
                tmp(1)= (-bl - sqrt(dll))/al
                tmp(2)= cl/(-bl - sqrt(dll))
            else 
                tmp(1)= cl/(-bl - sqrt(dll))
                tmp(2)= (-bl + izn*sqrt(dll))/al ! = cl/(-bl - sqrt(dll))
            endif
            !num_roots = num_roots + 2
        end if


        if (square_solver2(al, 2*bl, cl, rt)>0) then
            if (izn == 1) then
                root(1) = rt(1)
                root(2) = rt(2)
            else 
                root(1) = rt(2)
                root(2) = rt(1)
            endif

            err_1 = abs(root(1) - tmp(1)) !(-bl - sqrt(dll))/al
            err_2 = abs(root(2) - tmp(2)) !cl/(-bl - sqrt(dll))

            if (err_1+err_2> 1d-9) then 
                print *, 'errror - 1'
                print *, err_1, err_2
                print *, tmp(1:2)                
                print *, root(1:2)
                print *, izn
                pause
            else
                !print *, 'errror good'
                !print *, err_1, err_2
            endif

            num_roots = num_roots + 2
        end if
        
        cl1  = g11*yn2**2/xj + yn3**2/g33 - ynzq - ynpopq1
        dll1 = bl**2 - al*cl1
        print *, '--- izn=', izn
        if (dll1.ge.zero) then
            tmp(3)= (-bl - izn*sqrt(dll1))/al
            tmp(4)= (-bl + izn*sqrt(dll1))/al
            !num_roots = num_roots + 2
            !print *, root(3), root(4)
        end if
        !return

        if (square_solver2(al, 2*bl, cl1, rt)>0) then
            if (izn == 1) then
                root(3)= rt(1) !(-bl - sqrt(dll))/al
                root(4)= rt(2) !cl/(-bl - sqrt(dll))
            else 
                root(3)= rt(2) !cl/(-bl - sqrt(dll))
                root(4)= rt(1) !(-bl + izn*sqrt(dll))/al ! = cl/(-bl - sqrt(dll))
            endif
            err_1 = abs(tmp(3) - root(3))
            err_2 = abs(tmp(4) - root(4))               
            num_roots = num_roots + 2
            if (err_1+err_2 > 1d-9) then 
                print *, 'errror - 2'
                print *, err_1, err_2
                print *, tmp(3:4)
                print *, root(3:4)
                print *, izn
                pause
            else
                !print *, 'errror good'
                !print *, err_1, err_2
            endif
        end if

    end

    subroutine disp4(in_pa,ptet,xnr,yn2)
        use constants, only: zero, one, two, c0, c1,pi, zalfa, xlog, xmalfa
        use approximation
        use plasma, only: fn1, fn2, fvt, zefff, ft
        use plasma, only: cdl, cly, cgm, cmy, ncoef
        use plasma, only: r0, rm, b_tor, ww, xmi, xsz
        use plasma, only: cltn, cnye, cnyi, cnstal
        use rt_parameters, only: inew, itend0
        use metrics
        use dispersion_equation, only: ynz
        use decrements, only : dhdnr !!!!!
        use source_new_mod, only: source_new
        use xpoint
        implicit none
        real(wp), intent(in)    :: in_pa   ! ro
        real(wp), intent(in)    :: yn2     ! ???
        real(wp), intent(in)    :: ptet    ! theta
        real(wp), intent(inout) :: xnr     ! ???
            
        real(wp) :: pa, pn
        real(wp) :: fnr, fnrr
        real(wp) :: wpq
        real(wp) :: xdl, xdlp, xdlpp
        real(wp) :: xly, xlyp, xlypp, xlyv
        real(wp) :: xgm, xgmp, xgmpp
        real(wp) :: xmy, xmyp
        real(wp) :: dxrt
        real(wp) :: dxdrdr, dxdtdr
        real(wp) :: dzdrdr, dzdtdr
        real(wp) :: dvdr, du1dr, dudr, du1dt, dudt
        real(wp) :: cotet, sitet
        real(wp) :: x0r, xjr
        !real(wp) :: dxdt, dzdrdt
        real(wp) :: g2gq, g22q, g33q
        real(wp) :: bpr, btr
        real(wp) :: bar
        real(wp) :: sit, cot, sir, cor
        real(wp) :: v, vt, vpop, vpopr, vpopt, u, u1
        real(wp) :: sl1, pnye, pnyi, tmp
        real(wp) :: fcoll, dek1, dek2
        real(wp) :: e1, e2, e3
        real(wp) :: e1t, e1r, e2t, e2r, e3r
        real(wp) :: ynzr, ynzt, ynzqr, ynzqt, ynpopqr, ynpopqt
        real(wp) :: whe, wde3dw, wdbsdw, wdcsdw, wdhdw
        real(wp) :: yny, ynzq, ynyr, ynyt
        real(wp) :: gpr, dgpr, gdop, gprr, gdopr, gdopt
        real(wp) :: as, bs, cs 
        real(wp) :: asr, bsr, csr, ast, bst, cst
        real(wp) :: dhdv, dhdu
        real(wp) :: dnx, dny, dnz, ddn2
        real(wp) :: gr, gt, g11r, g22r, g12r
        real(wp) :: g22qr, g22qt, g33qr, g33qt
        real(wp) :: g2jqr
        real(wp) :: g33r
        real(wp) :: g2gqt, g2gqr
        real(wp) :: source
        real(wp) :: aimh
        !common /bcef/ ynz,ynpopq

        !integer :: irefl, iconv
        !common /cefn/ iconv,irefl

        !common /df/ pdec14,pdec24,pdec34,idec
        !common/metrika/g11,g12,g22,g33,gg,g,si,co
        !common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
        !common/fjham/ham
        pa = in_pa
        irefl=0
        iconv=0
        if(pa.eq.zero) pa=1.d-7 ! TODO не хорошо, что изменяется pa
        if(pa.lt.zero) pa=dabs(pa)
        !sav2008      if(pa.gt.one) then
        !sav2008       dhdm=666d0
        !sav2008       dhdtet=-666d0
        !sav2008       dhdnr=666d0
        !sav2008       dhdr=-666d0
        !sav2008       irefl=1
        !sav2008       return
        !sav2008      end if
  
          !!      pn=fn1(pa,fnr)
        !!      pn=fn2(pa,fnr,fnrr)
        if(inew.eq.0) then !vardens
         pn=fn1(pa,fnr)
        else
         pn=fn2(pa,fnr,fnrr)
        end if
  
        !cc        hstp=1.d-7
        !cc        pplus=fn2(pa+hstp,fnr2,fnrr2)
        !cc        pminus=fn2(pa-hstp,fnr1,fnrr1)
        !cc        fnr=0.5d0*(pplus-pminus)/hstp
        !cc        fnrr=0.5d0*(fnr2-fnr1)/hstp
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        wpq=c0**2*pn
        xdl=fdfddf(pa,cdl,ncoef,xdlp,xdlpp)
        xly=fdfddf(pa,cly,ncoef,xlyp,xlypp)
        xgm=fdfddf(pa,cgm,ncoef,xgmp,xgmpp)
        xmy=fdf(pa,cmy,ncoef,xmyp)
        cotet=dcos(ptet)
        sitet=dsin(ptet)
        xlyv=xly+xlyp*pa
        !--------------------------------------
        ! components of metric tensor
        !--------------------------------------
        dxdr=-xdlp+cotet-xgmp*sitet**2
        dxdt=-(pa+two*xgm*cotet)*sitet
        dzdr=xlyv*sitet
        dzdt=xly*pa*cotet
        x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
        dxdrdr=-xdlpp-xgmpp*sitet**2
        dxdtdt=-cotet*(pa+two*xgm*cotet)+sitet**2*two*xgm
        dxdtdr=-sitet*(one+two*xgmp*cotet)
        dxdrdt=dxdtdr
        dzdrdr=(two*xlyp+pa*xlypp)*sitet
        dzdtdt=-xly*pa*sitet
        dzdtdr=xlyv*cotet
        dzdrdt=dzdtdr
        x0t=dxdt
        x0r=dxdr
        if(-5.gt.0) then !!!!!!!!!!!xtochk
            call xpnt2(pa,ptet-1.5d0*pi,dxdr,dxdt,dzdr,dzdt,x0,dxdrdt,dxdtdt,dzdrdt,dzdtdt,x0t,dxdrdr,dxdtdr,dzdrdr,dzdtdr,x0r)
        end if
        g11=dxdr**2+dzdr**2
        g22=dxdt**2+dzdt**2
        g12=dxdr*dxdt+dzdr*dzdt
        g33=x0**2
        xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
        gg=xj
        g=xj*g33
        g2jq=dsqrt(g22/xj)
        g2gq=dsqrt(g22/g)
        g22q=dsqrt(g22)
        g33q=dsqrt(g33)
        !c--------------------------------------
        !c  magnetic field
        !c--------------------------------------
        bt=b_tor*(r0/rm)/x0
        bp=g2gq*xmy ! dsqrt(g22/(xj*g33))
        b=dsqrt(bp*bp+bt*bt)
        whe=b*c1
        si=bp/b
        co=bt/b
        !c---------------------------------------
        !c components of dielectric tensor
        !c---------------------------------------
        v=wpq/ww**2
        u1=whe/ww
        u=u1**2
        e1=one-v*(one/xmi-one/u)
        e2=v/u1
        e3=one-v
        !c-------------------------------------
        !c dispersion equation
        !c--------------------------------------
        ynz=yn2*si/g22q+yn3*co/g33q
        ynzq=ynz**2
        vpop=xnr**2*g22-two*xnr*yn2*g12+g11*yn2**2
        ynpopq=vpop/xj+yn3**2/g33-ynzq
        as=e1
        bs=-(e1**2-e2**2+e1*e3-(e1+e3)*ynzq)
        cs=e3*(e1**2-e2**2-two*e1*ynzq+ynzq**2)
        !c----------------------------------------------------
        !sav2009
        dhdv=(1.d0/(u**2*xmi**2))*((2.d0-2.d0*ynpopq-3.d0*v)*v*xmi**2-&
            u**2*(3.d0*v**2+2.d0*v*(-1.d0+ynpopq*(1.d0+xmi)&
            +2.d0*xmi*(-1.d0+ynzq))+&
            xmi*(-2.d0+ynpopq+xmi*(-1.d0+ynzq))*(-1.d0+ynpopq+ynzq))&
            +u*xmi*(3.d0*v**2*(2.d0+xmi)&
            +(-2.d0+ynpopq)*xmi*(-1.d0+ynpopq+ynzq)&
            +v*(-4.d0+4.d0*ynpopq*(1.d0+xmi)+xmi*(-6.d0+4.d0*ynzq))))

        dhdu=-(1.d0/(u**3*xmi))*(v*(-1.d0+ynpopq+v)*(2.d0*u*v-2.d0*v*xmi &
            +u*(-2.d0+ynpopq+v)*xmi)+u*v*(-2.d0+ynpopq+2.d0*v)*xmi*ynzq)
        dhdv2v=2.d0*v*dhdv !w*d(-H)/dv
        dhdu2u=2.d0*u*dhdu !w*d(-H)/du
        !c----------------------------------------------------
        !est !sav2009
        if(inew.gt.0) then
            if(inew.eq.1) then
                yny= - (yn2*co/g22q-yn3*si/g33q)
            else if(inew.eq.2) then
                yny= - g2jq*(yn2*co/g22q-yn3*si/g33q)
            end if
            gpr=c0**2/ww**2/u1*fnr*xsz
            gdop=yny*gpr
            bs=bs+gdop
            cs=cs+gdop*(ynzq-e3)
            !sav2009:
            dgpr=-gpr        !w*d(gpr)/dw
            wde3dw=2.d0*v    !w*d(e3)/dw
            wdbsdw=yny*dgpr  !w*d(bs)/dw
            wdcsdw=(ynzq-e3)*wdbsdw-gdop*wde3dw   !w*d(cs)/dw
            wdhdw=wdbsdw*ynpopq+wdcsdw   !w*d(H1)/dw
            dhdv2v=dhdv2v-wdhdw !correction to dhdv2v: w*d(-H)/dv+w*d(-H1)/dw
        end if
        ham=as*ynpopq**2+bs*ynpopq+cs !sav2009
        !--------------------------------------------------------
        !!      dl=bs**2-4d0*as*bs
        !c--------------------------------------
        !c   calculation of derivatives
        !c--------------------------------------
        g11r=two*dxdr*dxdrdr+two*dzdr*dzdrdr
        g22r=two*dxdt*dxdtdr+two*dzdt*dzdtdr
        g11t=two*dxdr*dxdrdt+two*dzdr*dzdrdt
        g22t=two*dxdt*dxdtdt+two*dzdt*dzdtdt
        g12r=dxdrdr*dxdt+dxdr*dxdtdr+dzdrdr*dzdt+dzdr*dzdtdr
        g12t=dxdrdt*dxdt+dxdr*dxdtdt+dzdrdt*dzdt+dzdr*dzdtdt
        g33r=two*x0*x0r
        g33t=two*x0*x0t
        g22qr=g22r/(g22q*two)
        g22qt=g22t/(g22q*two)
        g33qr=g33r/(g33q*two)
        g33qt=g33t/(g33q*two)
        xjr=g11r*g22+g22r*g11-two*g12*g12r
        xjt=g11t*g22+g22t*g11-two*g12*g12t
        g2jqr=(g22r/xj-g22/xj**2*xjr)/(g2jq*two) !sav2009
        g2jqt=(g22t/xj-g22/xj**2*xjt)/(g2jq*two) !sav2009
        gr=xjr*g33+g33r*xj
        gt=xjt*g33+g33t*xj
        g2gqt=(g22t/g-g22/g**2*gt)/(g2gq*two)
        g2gqr=(g22r/g-g22/g**2*gr)/(g2gq*two)
        bpt=xmy*g2gqt
        bpr=g2gqr*xmy+g2gq*xmyp
        btr=-b_tor*(r0/rm)/x0**2*x0r
        btt=-b_tor*(r0/rm)/x0**2*x0t
        bat=one/b*(bp*bpt+bt*btt)
        bar=one/b*(bp*bpr+bt*btr)
        sit=bpt/b-bp/b**2*bat
        cot=btt/b-bt/b**2*bat
        sir=bpr/b-bp/b**2*bar
        cor=btr/b-bt/b**2*bar
        dvdr=fnr*c0**2/ww**2
        du1dr=c1*bar/ww
        dudr=two*u1*du1dr
        du1dt=c1*bat/ww
        dudt=two*u1*du1dt
        e1r=-dvdr*(one/xmi-one/u)-v*dudr/u**2
        e1t=-v*dudt/u**2
        e2r=dvdr/u1-v/u1**2*du1dr
        e2t=-v/u1**2*du1dt
        e3r=-dvdr
        ynzr=yn2*(sir/g22q-si/g22q**2*g22qr) + yn3*(cor/g33q-co/g33q**2*g33qr)
        ynzt=yn2*(sit/g22q-si/g22q**2*g22qt) + yn3*(cot/g33q-co/g33q**2*g33qt)
        ynzqr=two*ynz*ynzr
        ynzqt=two*ynz*ynzt
        vpopr=(xnr**2*g22r-two*xnr*yn2*g12r+yn2**2*g11r)
        vpopt=(xnr**2*g22t-two*xnr*yn2*g12t+yn2**2*g11t)
        ynpopqr=vpopr/xj-vpop/xj**2*xjr-yn3**2/g33**2*g33r-ynzqr
        ynpopqt=vpopt/xj-vpop/xj**2*xjt-yn3**2/g33**2*g33t-ynzqt
        asr=e1r
        bsr=(e3r+e1r)*(ynzq-e1)+(e3+e1)*(ynzqr-e1r)+two*e2*e2r
        csr=e3r*((ynzq-e1)**2-e2**2)+e3*(two*(ynzq-e1)*(ynzqr-e1r) - two*e2*e2r)
        ast=e1t
        bst=e1t*(ynzq-e1)+(e3+e1)*(ynzqt-e1t)+two*e2*e2t
        cst=e3*(two*(ynzq-e1)*(ynzqt-e1t)-two*e2*e2t)
        !---------------------------------------------------
        !est !sav2009
        if (inew.gt.0) then
            if(inew.eq.1) then
                ynyr= - (yn2*(cor/g22q-co/g22q**2*g22qr)&
                        -yn3*(sir/g33q-si/g33q**2*g33qr))
                ynyt= - (yn2*(cot/g22q-co/g22q**2*g22qt)-&
                        -yn3*(sit/g33q-si/g33q**2*g33qt))
            else if(inew.eq.2) then
                ynyr= - g2jq*(yn2*(cor/g22q-co/g22q**2*g22qr)&
                                -yn3*(sir/g33q-si/g33q**2*g33qr))
                ynyr=ynyr - g2jqr*(yn2*co/g22q-yn3*si/g33q)
                ynyt= - g2jq*(yn2*(cot/g22q-co/g22q**2*g22qt)-&
                                -yn3*(sit/g33q-si/g33q**2*g33qt))
                ynyt=ynyt - g2jqt*(yn2*co/g22q-yn3*si/g33q)
            end if
            gprr=c0**2/ww**2*(fnrr/u1-fnr/u1**2*du1dr)*xsz
            gprt=-c0**2/ww**2*fnr/u1**2*du1dt*xsz
            gdopr=ynyr*gpr+yny*gprr
            gdopt=ynyt*gpr+yny*gprt
            bsr=bsr+gdopr
            csr=csr+gdopr*(ynzq-e3)+gdop*(ynzqr-e3r)
            bst=bst+gdopt
            cst=cst+gdopt*(ynzq-e3)+gdop*ynzqt
        end if
        !c---------------------------------------------------
  
        dhdr=asr*ynpopq**2+bsr*ynpopq+as*two*ynpopq*ynpopqr+bs*ynpopqr+csr

        dhdtet=ast*ynpopq**2+bst*ynpopq+as*two*ynpopq*ynpopqt+bs*ynpopqt+cst

        dnx=two*as*ynpopq+bs
        dnz=ynpopq*(e1+e3)+two*(ynzq-e1)*e3
        dhdnr=dnx*two*(g22*xnr-g12*yn2)/xj
        dhdm=dnx*two*(yn2*g11-xnr*g12)/xj+(dnz-dnx)*two*ynz*si/g22q
        !sav2009
        dhdn3=two*((yn3-ynz*co*g33q)*dnx+ynz*co*g33q*dnz)/g33 !sav2009
        !c----------------------------------------------------------------
        !est !sav2009
        if (inew.gt.0) then
            if(inew.eq.1) then
                dny= - gpr*(ynpopq+ynzq-e3)
            else if(inew.eq.2) then
                dny= - g2jq*gpr*(ynpopq+ynzq-e3)
            end if
            dhdm=dhdm+dny*co/g22q+two*ynz*yny*gpr*si/g22q
            dhdn3=dhdn3-dny*si/g33q+two*ynz*yny*gpr*co/g33q !sav2009
        end if
        !c---------------------------------------------------
        !sav2009
        ddn2=g11*dhdnr**2+g22*dhdm**2+2.d0*g12*dhdnr*dhdm+g33*dhdn3**2
        ddn=dsqrt(ddn2)
        !
        if(idec.ne.0) then
            vt=fvt(pa)
            sl1=(ynzq-e1)*(ynzq+ynpopq-e1)-e2**2
            aimh=wpq/ww**2*pi*sl1*cltn**2/ynzq
            pdec14=dabs(aimh/xsz/ddn)
            pnye=cnye*wpq**2/(pn*vt**3)
            pnyi=cnyi*pnye*zefff(pa)
            pdec24=dabs(pnyi/ww*(wpq/whe**2*ynpopq+wpq/ww**2*ynzq)*ynpopq/&
            xsz/ddn)
            if(itend0.gt.0) then
                tmp=ft(pa)/0.16d-8
                fcoll=.5d-13*pn*zalfa**2*xlog/xmalfa/tmp**1.5d0
                !cc        ddens=dn1*pn
                !cc        tdens=dn2*pn
                !cc        tt=fti(pa)**0.33333d0    ! (ti, kev)^1/3
                !cc        source=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
                call source_new(pa, source)
                dek1=cnstal*pdec14*(1.d0-e3/ynpopq)**2/dsqrt(ynpopq)
                dek2=source/(fcoll*pn)
                pdec34=dek1*dek2
            end if
        end if
    end


    function dhdomega(rho,theta,yn1,yn2) result(znak) 
        !! вычисляет znakstart
        use decrements, only : dhdnr !!!!!
        implicit none
        real(wp), intent(in) :: rho
        real(wp), intent(in) :: theta
        real(wp), intent(inout) :: yn1
        real(wp), intent(in) :: yn2
        real(wp) :: wdhdw, znak
        !common /a0ef2/ ww
        !common /abefo/ yn3
        !common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
        !common/direct/znakstart
        !parameter(zero=0.d0,h=1.d-6)
        
        call disp4(rho, theta, yn1, yn2)
        
        !!w*dH/dw=wdhdw:
        wdhdw=-(yn1*dhdnr+yn2*dhdm+yn3*dhdn3+dhdv2v+dhdu2u)
        znak=dsign(1.d0,wdhdw)
        !znakstart=znak
        !c      write(*,*)'formula: znak=',znak
        !c      write(*,*)'wdhdw=',wdhdw,' H=',ham
        !c      write(*,*)'rho=',rho,' teta=',theta
        !c      write(*,*)'yn1=',yn1,' yn2=',yn2
        !c      write(*,*)'dhdnr=',dhdnr,' dhdm=',dhdm
        !c      write(*,*)'dhdr=',dhdr,' dhdtet=',dhdtet
        !c      write(*,*)'dhdn3=',dhdn3,' yn3=',yn3
        !c      write(*,*)'yn1*dhdnr=',yn1*dhdnr,' yn2*dhdm=',yn2*dhdm
        !c      write(*,*)'yn1*dhdnr+yn2*dhdm=',yn1*dhdnr+yn2*dhdm
        !cc      pause
  
        end    

    
    subroutine extd4(x,y,dydx)
        !use dispersion_module
        use decrements, only : dhdnr !!!!!
        implicit none
        !implicit real(wp) (a-h,o-z)
        real(wp), intent(in)    :: x
        real(wp), intent(in)    :: y(:)
        real(wp), intent(inout) :: dydx(:)
        !common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
        !common/direct/znakstart
        real(wp) :: znak, xxx, ptet, yn2, pa, yn1
        znak=znakstart
        xxx=x
        ptet=y(1)
        yn2=y(2)
        pa=y(3)
        yn1=y(4)
        call disp4(pa,ptet,yn1,yn2)
        !new variant
        dydx(1)=-znak*dhdm/ddn
        dydx(2)=znak*dhdtet/ddn
        dydx(3)=-znak*dhdnr/ddn
        dydx(4)=znak*dhdr/ddn
        !dydx(5)=-znak*dhdn3/ddn !! возможен выход за границы массива
  
        !c      dydx(1)=znak*dhdm/ddn
        !c      dydx(2)=-znak*dhdtet/ddn
        !c      dydx(3)=znak*dhdnr/ddn
        !c      dydx(4)=-znak*dhdr/ddn
        !c      dydx(5)=znak*dhdn3/ddn
        
        !old variant:
        !      dydx(1)=dhdm/ddn
        !      dydx(2)=-dhdtet/ddn
        !      dydx(3)=dhdnr/ddn
        !      dydx(4)=-dhdr/ddn
        !      dydx(5)=dhdn3/ddn
    end

 
end module dispersion_module