module maxwell
      !! Все что относится к распределению Максвелла
      use kind_module      
      use constants, only : zero, pisqrt, pi2sqrt, pqe
      implicit none
      integer, parameter :: i0 = 1002


      logical flag_d0
      !! бывший d0
      integer jindex, kindex
      !!common/dddql/ d0,jindex,kindex
      
contains

      function currlhcd(v,f) result(curs)
            implicit none
            real(wp), intent(in) :: v(:),f(:)
            real(wp) curs      
            integer i0,k

            real(wp) vl,vr,fl,fr
            curs=0.d0
            i0 = size(v)
            do k=1,i0-1
                  vl=v(k)
                  vr=v(k+1)
                  fl=f(k)
                  fr=f(k+1)
                  curs=curs+(fl*vl+fr*vr)/2d0*(vr-vl)
            end do
      end
      
      function create_vt_grid(vclt) result(vt_grid)
            !! создание сетки тепловых скоростей
            implicit none
            real(wp), intent(in) :: vclt
            real(wp) :: vt_grid(i0)            
            real(wp) vmax
            integer i
            !r=dble(j)/dble(nr+1)
            !vclt=3.d10/fvt(r)
            vmax = 2.d0*vclt
            do i=1,i0
                  vt_grid(i)=dble(i-1)*vmax/dble(i0-1)
            end do
      end function

      subroutine init_vi(vclt, vi)
            real(wp), intent(in) :: vclt
            real(wp), intent(out) :: vi(i0)            
            real(wp) vmax
            integer i
            vmax = 2.d0*vclt
            do i=1,i0
                  vi(i)=dble(i-1)*vmax/dble(i0-1)
            end do
      end subroutine      

      subroutine init_fmaxw_classic(vclt, enorm, fi, dfi)
            real(wp), intent(in) :: vclt, enorm
            real(wp), intent(out) :: fi(i0), dfi(i0)
            real(wp) vi, vmax
            integer i
            vmax = 2.d0*vclt
            do i=1,i0
                  vi = dble(i-1)*vmax/dble(i0-1)
                  if(vi < vclt) then
                        fi(i) = fmaxw_classic(vi, enorm, dfi(i))
                  else
                        fi(i) = zero
                        dfi(i) = zero
                  end if
            end do
      end subroutine      

      subroutine init_fmaxw_ext(vclt, enorm, fi, dfi)
            real(wp), intent(in) :: vclt, enorm
            real(wp), intent(out) :: fi(i0), dfi(i0)
            real(wp) vi, vmax
            integer i
            vmax=2.d0*vclt
            do i=1,i0
                  vi = dble(i-1)*vmax/dble(i0-1)
                  if(vi < vclt) then
                        fi(i) = fmaxw_ext(vi, enorm, dfi(i))
                  else
                        fi(i) = zero
                        dfi(i) = zero
                  end if
            end do
      end subroutine      

      real(wp) function funmaxwell(v,dfunmaxwell)
            !! распределение Максвелла 
            !! $$ f(v) = \frac{1}{\sqrt{2\pi}} \exp(-\frac{1}{2} v^2 ))$$
            !! и его производная
            !! $$ dfmaxw = - v \cdot f(v) $$
            implicit none
            real(wp) v, dfunmaxwell, arg

            arg = -0.5d0*v**2
            funmaxwell = exp(arg)/pi2sqrt
            dfunmaxwell = -v*funmaxwell
      end      

      real(wp) function fmaxw_classic(v,alfa2,dfmaxw)
            !! распределение Максвелла с альфа-частицами
            !! $$ f(v) = \frac{1}{\sqrt{2\pi}} \exp(-\frac{1}{2} v^2 (1.0 + \frac{1}{2} \cdot alfa_2 \cdot v^2))$$
            !! и его производная
            !! $$ dfmaxw = - v \cdot (1.0 + alfa_2 \cdot v^2) \cdot f(v) $$
            implicit none
            real(wp) v,alfa2,dfmaxw
            real(wp) arg,alfa,api,b,psiq,f,df
            
            arg=-0.5d0*v**2*(1.d0+0.5d0*alfa2*v**2)
            fmaxw_classic=dexp(arg)/pi2sqrt
            dfmaxw=-v*(1.d0+alfa2*v**2)*fmaxw_classic
      end

      real(wp) function fmaxw_ext(v,alfa2,dfmaxw)
            !! $$ alfa = \sqrt{alfa_2} $$
            !! $$ api = 2 \cdot alfa \cdot \exp({-\frac{1}{4 alfa_2}}) $$
            !! $$ b = 2 - erf(0.5/alfa) + api $$
            !! $$ f = psiq(v, alfa_2) $$
            !! $$ fmaxw = \frac{f+api}{b \sqrt{2\pi}} $$
            implicit none
            real(wp) v,alfa2,dfmaxw
            real(wp) arg,alfa,api,b,f,df

            alfa=dsqrt(alfa2)
            api=2.d0*alfa*dexp(-0.25d0/alfa2)/pisqrt
            b=2.d0-erfcc(0.5d0/alfa)+api
            f=psiq(v,alfa2)
            fmaxw_ext=(f+api)/b/pi2sqrt
            df=-v*((1.d0-alfa2*v**2)*f+api)
            dfmaxw=df/b/pi2sqrt
      end      

      real(wp) function fmaxw(v,alfa2,dfmaxw)
            implicit none
            real(wp) v,alfa2,dfmaxw
            real(wp) arg,alfa,api,b,f,df

            if(alfa2.le.zero) then
                  arg=-0.5d0*v**2*(1.d0-0.5d0*alfa2*v**2)
                  fmaxw=dexp(arg)/pi2sqrt
                  dfmaxw=-v*(1.d0-alfa2*v**2)*fmaxw
            else
                  alfa=dsqrt(alfa2)
                  api=2.d0*alfa*dexp(-0.25d0/alfa2)/pisqrt
                  b=2.d0-erfcc(0.5d0/alfa)+api
                  f=psiq(v,alfa2)
                  fmaxw=(f+api)/b/pi2sqrt
                  df=-v*((1.d0-alfa2*v**2)*f+api)
                  dfmaxw=df/b/pi2sqrt
            end if
      end      

      real(wp) function psiq(v,alfa2)
            !! $$ psiq=exp(ksiV**2)*erfcc(ksiV)*exp(-0.25/alfa2) $$
            implicit none
            real(wp) v,alfa2,df
            real(wp) x,t,z,f,asymp,alfa,q,u
            real(wp), parameter :: zmax = 10.d0
            
            alfa=dsqrt(alfa2)
            q=-0.25d0/alfa2
            x=0.5d0*(alfa*v**2-1.d0/alfa)
            z=abs(x)
            if(z.gt.zmax) then !asymptotics
                  f=dexp(q)*(1.d0-0.5d0/z**2+0.75d0/z**4-15.d0/8.d0/z**6)/z/pisqrt
            else
                  t=1.d0/(1.d0+0.5d0*z)
                  f=t*exp(q-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*&
                  &(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*&
                  &(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
            end if
            if(x.lt.zero) then
                  u=-0.5d0*v**2+0.25d0*alfa2*v**4 !u=x**2-0.25d0/alfa2
                  f=2.d0*dexp(u)-f
            end if
            psiq=f
            return
            end

      function erfcc(x)
            implicit none
            real(wp) erfcc,x
            real(wp) t,z
            real(wp), parameter :: zmax = 10.d0
            
            z=abs(x)
            if(z.gt.zmax) then !asymptotics
                  erfcc=(1.d0-0.5d0/z**2+0.75d0/z**4-15.d0/8.d0/z**6)/z/pisqrt
                  erfcc=exp(-z*z)*erfcc
            else
                  t=1.d0/(1.d0+0.5d0*z)
                  erfcc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*&
                  &(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*&
                  &(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
            end if
            if(x.lt.zero) erfcc=2.d0-erfcc
            return
            end      
end module maxwell
