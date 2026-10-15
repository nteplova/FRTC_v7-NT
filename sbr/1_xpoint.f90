module xpoint
      !! Все что относится к распределению Максвелла
      use kind_module      
      use constants, only : zero, pisqrt, pi2sqrt, pqe
      implicit none
      !integer, parameter :: i0 = 1002


      logical flag_d0
      !! бывший d0
      integer jindex, kindex
      !!common/dddql/ d0,jindex,kindex
      
contains

    subroutine readxpnt
      implicit real*8 (a-h,o-z)
      common/b/ apx,apz,bpx,bpz,cfx,cfz,alz
      common/arrays/ cx(10,100),cz(10,100),mp,mh
      integer  mh,mp,lunin,k,i
      


      lunin=41
      open(lunin,file='lhcd/input',status='old')
       read(lunin,*) apx
       read(lunin,*) apz
       read(lunin,*) bpx
       read(lunin,*) bpz
       read(lunin,*) cfx
       read(lunin,*) cfz
       read(lunin,*) alz
       read(lunin,*) mh
       read(lunin,*) mp
       do k=1,mh
        do i=1,mp
         read(lunin,*) cx(i,k),cz(i,k)
        end do
       end do

      close(lunin) 



    end

     


    subroutine xpnt(r,t,dxdr,dxdt,dzdr,dzdt,dxdf,dxdrdt,dxdtdt,dzdrdt,dzdtdt,dxdfdt)
      implicit real*8 (a-h,o-z)
      !external fbas,dfbas,ddfbas
      dimension wx(100),dwx(100),ddwx(100)
      dimension wz(100),dwz(100),ddwz(100)
      common/b/ apx,apz,bpx,bpz,cfx,cfz,alz
      common/arrays/ cx(10,100),cz(10,100),mp,mh
      parameter(alfa=5.d0,beta=3.d-5)
      integer  mh,mp
      !mp=11
      !mh=3
      call fcoeff(r,wx,dwx,ddwx,wz,dwz,ddwz)
      cotet=dcos(t)
      sitet=dsin(t)
      tt=1d0-cotet
      rr=1d0-r
      s=tt+rr
      sa=tt+alz*rr
      ss=dsqrt(sa)
      fx=dexp(-cfx*s)
      fz=dexp(-cfz*s)
      fxr=-fx
       fxrr=-fxr
        fxt=fx*sitet
         fxrt=fxr*sitet
          fxtt1=fx*cotet
           fxtt2=fx*sitet**2
      fzr=-fz
       fzrr=-fzr
        fzt=fz*sitet
         fzrt=fzr*sitet
          fztt1=fz*cotet
           fztt2=fz*sitet**2
      ux=cfx*(bpx*ss-apx)
      uz=cfz*(bpz*ss-apz)
      ox=bpx/(2d0*ss)
       oxx=ox/(2d0*sa)
        oz=bpz/(2d0*ss)
         ozz=oz/(2d0*sa)
      uxr=fxr*(ux-ox*alz)
      uzr=fzr*(uz-oz*alz)
      uxt=fxt*(ux-ox)
      uzt=fzt*(uz-oz)
      r2=r**2
      r4=r2**2
      fe=dexp(alfa*r2)
      f=r4*fe
      ff=(beta+f)
      f_new=f/ff
      f1=2.d0*r2*r*(2.d0+alfa*r2)*fe
      ff2=ff**2
      df_new=beta*f1/ff2
      addx=(apx-ss*bpx)*fx-summ(t,wx,mh,fbas)
      addz=(apz-ss*bpz)*fz-summ(t,wz,mh,fbas)
      adxdr=uxr-summ(t,dwx,mh,fbas)
      adzdr=uzr-summ(t,dwz,mh,fbas)
      adxdt=uxt-summ(t,wx,mh,dfbas)
      adzdt=uzt-summ(t,wz,mh,dfbas)
      dxdr=dxdr+f_new*adxdr+addx*df_new
      dzdr=dzdr+f_new*adzdr+addz*df_new
      dxdt=dxdt+f_new*adxdt
      dzdt=dzdt+f_new*adzdt
      dxdf=dxdf+f_new*addx
      uxrt=fxrt*(alz*oxx+cfx*(ox*(1d0+alz)-ux))
      uzrt=fzrt*(alz*ozz+cfz*(oz*(1d0+alz)-uz))
      uxtt=fxtt1*(ux-ox)+fxtt2*(oxx+cfx*(ox*2d0-ux))
      uztt=fztt1*(uz-oz)+fztt2*(ozz+cfz*(oz*2d0-uz))
      adxdrdt=uxrt-summ(t,dwx,mh,dfbas)
      adzdrdt=uzrt-summ(t,dwz,mh,dfbas)
      adxdtdt=uxtt-summ(t,wx,mh,ddfbas)
      adzdtdt=uztt-summ(t,wz,mh,ddfbas)
      dxdrdt=dxdrdt+f_new*adxdrdt+adxdt*df_new
      dzdrdt=dzdrdt+f_new*adzdrdt+adzdt*df_new
      dxdtdt=dxdtdt+f_new*adxdtdt
      dzdtdt=dzdtdt+f_new*adzdtdt
      dxdfdt=dxdfdt+f_new*adxdt
    end

    subroutine xpnt_coord(r,t,addx,addz)
      implicit real*8 (a-h,o-z)
      !external fbas,dfbas,ddfbas
      dimension wx(100),dwx(100),ddwx(100)
      dimension wz(100),dwz(100),ddwz(100)
      common/b/ apx,apz,bpx,bpz,cfx,cfz,alz
      common/arrays/ cx(10,100),cz(10,100),mp,mh
      parameter(alfa=5.d0,beta=3.d-5)
      integer  mh,mp
      !mp=11
      !mh=3
      call fcoeff(r,wx,dwx,ddwx,wz,dwz,ddwz)
      cotet=dcos(t)
      sitet=dsin(t)
      tt=1d0-cotet
      rr=1d0-r
      s=tt+rr
      sa=tt+alz*rr
      ss=dsqrt(sa)
      fx=dexp(-cfx*s)
      fz=dexp(-cfz*s)
      fxr=-fx
       fxrr=-fxr
        fxt=fx*sitet
         fxrt=fxr*sitet
          fxtt1=fx*cotet
           fxtt2=fx*sitet**2
      fzr=-fz
       fzrr=-fzr
        fzt=fz*sitet
         fzrt=fzr*sitet
          fztt1=fz*cotet
           fztt2=fz*sitet**2
      ux=cfx*(bpx*ss-apx)
      uz=cfz*(bpz*ss-apz)
      ox=bpx/(2d0*ss)
       oxx=ox/(2d0*sa)
        oz=bpz/(2d0*ss)
         ozz=oz/(2d0*sa)
      uxr=fxr*(ux-ox*alz)
      uzr=fzr*(uz-oz*alz)
      uxt=fxt*(ux-ox)
      uzt=fzt*(uz-oz)
      r2=r**2
      r4=r2**2
      fe=dexp(alfa*r2)
      f=r4*fe
      ff=(beta+f)
      f_new=f/ff
      f1=2.d0*r2*r*(2.d0+alfa*r2)*fe
      ff2=ff**2
      df_new=beta*f1/ff2
      addx=(apx-ss*bpx)*fx-summ(t,wx,mh,fbas)
      addz=(apz-ss*bpz)*fz-summ(t,wz,mh,fbas)
      adxdr=uxr-summ(t,dwx,mh,fbas)
      adzdr=uzr-summ(t,dwz,mh,fbas)
      adxdt=uxt-summ(t,wx,mh,dfbas)
      adzdt=uzt-summ(t,wz,mh,dfbas)
      dxdr=dxdr+f_new*adxdr+addx*df_new
      dzdr=dzdr+f_new*adzdr+addz*df_new
      dxdt=dxdt+f_new*adxdt
      dzdt=dzdt+f_new*adzdt
      dxdf=dxdf+f_new*addx
      uxrt=fxrt*(alz*oxx+cfx*(ox*(1d0+alz)-ux))
      uzrt=fzrt*(alz*ozz+cfz*(oz*(1d0+alz)-uz))
      uxtt=fxtt1*(ux-ox)+fxtt2*(oxx+cfx*(ox*2d0-ux))
      uztt=fztt1*(uz-oz)+fztt2*(ozz+cfz*(oz*2d0-uz))
      adxdrdt=uxrt-summ(t,dwx,mh,dfbas)
      adzdrdt=uzrt-summ(t,dwz,mh,dfbas)
      adxdtdt=uxtt-summ(t,wx,mh,ddfbas)
      adzdtdt=uztt-summ(t,wz,mh,ddfbas)
      dxdrdt=dxdrdt+f_new*adxdrdt+adxdt*df_new
      dzdrdt=dzdrdt+f_new*adzdrdt+adzdt*df_new
      dxdtdt=dxdtdt+f_new*adxdtdt
      dzdtdt=dzdtdt+f_new*adzdtdt
      dxdfdt=dxdfdt+f_new*adxdt
    end



    subroutine xpnt2(r,t,dxdr,dxdt,dzdr,dzdt,dxdf,dxdrdt,dxdtdt,dzdrdt,dzdtdt,dxdfdt,dxdrdr,dxdtdr,dzdrdr,dzdtdr,dxdfdr)
      implicit real*8 (a-h,o-z)
      !external fbas,dfbas,ddfbas
      dimension wx(100),dwx(100),ddwx(100)
      dimension wz(100),dwz(100),ddwz(100)
      common/b/ apx,apz,bpx,bpz,cfx,cfz,alz
      common/arrays/ cx(10,100),cz(10,100),mp,mh
      parameter(alfa=5.d0,beta=3.d-5)
      integer mh,mp
      call fcoeff(r,wx,dwx,ddwx,wz,dwz,ddwz)
      cotet=dcos(t)
      sitet=dsin(t)
      tt=1d0-cotet
      rr=1d0-r
      s=tt+rr
      sa=tt+alz*rr
      ss=dsqrt(sa)
      fx=dexp(-cfx*s)
      fz=dexp(-cfz*s)
      fxr=-fx
       fxrr=-fxr
        fxt=fx*sitet
         fxrt=fxr*sitet
          fxtt1=fx*cotet
           fxtt2=fx*sitet**2
      fzr=-fz
       fzrr=-fzr
        fzt=fz*sitet
         fzrt=fzr*sitet
          fztt1=fz*cotet
           fztt2=fz*sitet**2
      ux=cfx*(bpx*ss-apx)
      uz=cfz*(bpz*ss-apz)
      ox=bpx/(2d0*ss)
      oxx=ox/(2d0*sa)
      oz=bpz/(2d0*ss)
      ozz=oz/(2d0*sa)
      uxr=fxr*(ux-ox*alz)
      uzr=fzr*(uz-oz*alz)
      uxt=fxt*(ux-ox)
      uzt=fzt*(uz-oz)
      r2=r**2
      r4=r2**2
      fe=dexp(alfa*r2)
      f=r4*fe
      ff=(beta+f)
      f_new=f/ff
      f1=2.d0*r2*r*(2.d0+alfa*r2)*fe
      f2=4.d0*r2*(3.d0+4.5d0*alfa*r2+r4*alfa**2)*fe
      ff2=ff**2
      df_new=beta*f1/ff2
      ddf_new=beta*(f2-2.d0*f1**2/ff)/ff2
      addx=(apx-ss*bpx)*fx-summ(t,wx,mh,fbas)
      addz=(apz-ss*bpz)*fz-summ(t,wz,mh,fbas)
      adxdr=uxr-summ(t,dwx,mh,fbas)
      adzdr=uzr-summ(t,dwz,mh,fbas)
      adxdt=uxt-summ(t,wx,mh,dfbas)
      adzdt=uzt-summ(t,wz,mh,dfbas)
      dxdr=dxdr+f_new*adxdr+addx*df_new
      dzdr=dzdr+f_new*adzdr+addz*df_new
      dxdt=dxdt+f_new*adxdt
      dzdt=dzdt+f_new*adzdt
      dxdf=dxdf+f_new*addx
      uxrt=fxrt*(alz*oxx+cfx*(ox*(1d0+alz)-ux))
      uzrt=fzrt*(alz*ozz+cfz*(oz*(1d0+alz)-uz))
      uxtt=fxtt1*(ux-ox)+fxtt2*(oxx+cfx*(ox*2d0-ux))
      uztt=fztt1*(uz-oz)+fztt2*(ozz+cfz*(oz*2d0-uz))
      adxdrdt=uxrt-summ(t,dwx,mh,dfbas)
      adzdrdt=uzrt-summ(t,dwz,mh,dfbas)
      adxdtdt=uxtt-summ(t,wx,mh,ddfbas)
      adzdtdt=uztt-summ(t,wz,mh,ddfbas)
      dxdrdt=dxdrdt+f_new*adxdrdt+adxdt*df_new
      dzdrdt=dzdrdt+f_new*adzdrdt+adzdt*df_new
      dxdtdt=dxdtdt+f_new*adxdtdt
      dzdtdt=dzdtdt+f_new*adzdtdt
      dxdfdt=dxdfdt+f_new*adxdt
      dxdtdr=dxdrdt
      dzdtdr=dzdrdt
      uxrr=fxrr*(oxx*alz**2+cfx*(ox*2d0*alz-ux))
      uzrr=fzrr*(ozz*alz**2+cfz*(oz*2d0*alz-uz))
      adxdrdr=uxrr-summ(t,ddwx,mh,fbas)
      adzdrdr=uzrr-summ(t,ddwz,mh,fbas)
      dxdrdr=dxdrdr+f_new*adxdrdr+2.d0*df_new*adxdr+addx*ddf_new
      dzdrdr=dzdrdr+f_new*adzdrdr+2.d0*df_new*adzdr+addz*ddf_new
      dxdfdr=dxdfdr+f_new*adxdr+addx*df_new
    end




    double precision function fbas(k,x,m)
      implicit real*8 (a-h,o-z)
      integer k, m
      real*8 x
      real*8 mh, q1, q2
       mh=(m-1)/2
       q1=dble(k-1)
       q2=dble(k-mh-1)
         if(k.le.mh+1) then
           fbas=dcos(x*q1)
         else
           fbas=dsin(x*q2)
         end if
      return
    end


    double precision function dfbas(k,x,m)
      implicit real*8 (a-h,o-z)
      integer k, m
      real*8 x
      real*8 mh, q1, q2
       mh=(m-1)/2
       q1=dble(k-1)
       q2=dble(k-mh-1)
         if(k.le.mh+1) then
           dfbas=-q1*dsin(x*q1)
         else
           dfbas=q2*dcos(x*q2)
         end if
      return
    end


    double precision function ddfbas(k,x,m)
      implicit real*8 (a-h,o-z)
      integer k, m
      real*8 x
      real*8 mh, q1, q2
       mh=(m-1)/2
       q1=dble(k-1)
       q2=dble(k-mh-1)
         if(k.le.mh+1) then
           ddfbas=-q1**2*dcos(x*q1)
         else
           ddfbas=-q2**2*dsin(x*q2)
         end if
      return
    end


    subroutine fcoef(r,wx,wz)
      implicit real*8 (a-h,o-z)
      integer k,mh,mp
      integer i
      dimension cffx(10),cffz(10)
      dimension wx(100),wz(100)
      common/arrays/ cx(10,100),cz(10,100),mp,mh
          do k=1,mh
            do i=1,mp
             cffx(i)=cx(i,k)
             cffz(i)=cz(i,k)
            end do
           wx(k)=fdffdf(r,cffx,mp,df)
           wz(k)=fdffdf(r,cffz,mp,df)
          end do
      return
    end


    subroutine fcoeff(r,wx,dwx,ddwx,wz,dwz,ddwz)
      implicit real*8 (a-h,o-z)
      integer k,mh,mp
      integer i
      dimension cffx(10),cffz(10)
      dimension wx(100),dwx(100),ddwx(100)
      dimension wz(100),dwz(100),ddwz(100)
      common/arrays/ cx(10,100),cz(10,100),mp,mh
          do k=1,mh
            do i=1,mp
             cffx(i)=cx(i,k)
             cffz(i)=cz(i,k)
            end do
           wx(k)=fdfddff(r,cffx,mp,dwx(k),ddwx(k))
           wz(k)=fdfddff(r,cffz,mp,dwz(k),ddwz(k))
          end do
      return
    end


    double precision function summ(x,c,m,f)
      implicit real*8 (a-h,o-z)
      integer i,m
      dimension c(m)
       out=0d0
        do i=1,m
         out=out+c(i)*f(i,x,m)
        end do
       summ=out
    end

    double precision function fdffdf(x,c,n,df)
      implicit real*8 (a-h,o-z)
      dimension c(n)
      integer i,j,n
      p=c(n)
      dp=0.d0
      do j=n-1,1,-1
        dp=dp*x+p
        p=p*x+c(j)
      end do
      fdffdf=p
      df=dp
    end



    double precision function fdfddff(x,c,n,df,ddf)
      implicit real*8 (a-h,o-z)
      dimension c(n)
      integer i,j,n
      p=c(n)
      dp=0d0
      ddp=0d0
      do j=n-1,1,-1
        ddp=ddp*x+2d0*dp
        dp=dp*x+p
        p=p*x+c(j)
      end do
      fdfddff=p
      df=dp
      ddf=ddp
    end
       
end module xpoint
