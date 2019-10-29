Program stany
    implicit none
    real dr, dU,dr2,dx,dy,dxp,dyp,drp2
    real T,q,rad,drp
    integer d,i,j,k,o,m,licznik,z,licznik2
    integer N,L,Mcs,Rmax,ch
    integer check,time,ale
    parameter (T=0.7,L=35, Mcs=230000,dr=0.01,drp=0.01,q=0.03,rad=1.0,N=L*L*q,Rmax = int(L/2.0 - 1.0),check=30000,time=100,ale=(Rmax/drp))!dr-ruch atomu,drp-szerokosc pierscienia,q-gestosc,rad-radiusczast,Rmax-sprawdzamy dotad pierscien
    real x(0:N),y(0:N),xp(0:N),yp(0:N),ra!ra odleglosc do prawdo
    real*8 ran1,w,r
    real start, finish !time count
    real fun(0:ale)
    real a,dra
    !----------------------------------------------------
    d=-1
    licznik=0
    licznik2=0
    ch=0
    open(10,file="co.dat")
    open(12,file="prob.dat")
    !------------Polozenia-------------------------------
    do i=1,(Rmax/drp)
        fun(i)=0
    enddo    
    do i=0,N
100     continue   
        xp(i)=ran1(d)*L
200     continue   
        yp(i)=ran1(d)*L
        do j=0,N
            if(xp(i).eq.x(j))then
                goto 100
            endif
            if(yp(i).eq.y(j))then
                goto 200
            endif
            dy=yp(i)-y(j)
            dx=xp(i)-x(j)
            if(sqrt(dy**2+dx**2).le.rad)then
                goto 100
            endif
        enddo
    y(i)=yp(i)
    x(i)=xp(i)
    enddo        
    !------------algorytm MC----------------------------
        licznik=0
    do k=1, Mcs
        !call cpu_time(start)
        licznik2=licznik2+1
        write(*,*)licznik2
        do o=1,N
            write(10,*)x(o),y(o)
        enddo
        rewind(10)
        do i=0,N
            dU=0
            r=ran1(d)
300         continue
            xp(i)=x(i)+(r-0.5)*dr! nowa wartosc
400         continue
            r=ran1(d)
            yp(i)=y(i)+(r-0.5)*dr
            if(xp(i).gt.L)then!tak, aby nie wyszlo za przedzial L
                xp(i)=xp(i)-L
            endif
            if(yp(i).gt.L)then
                yp(i)=yp(i)-L
            endif
            if(xp(i).lt.0.)then
                xp(i)=xp(i)+L
            endif
            if(yp(i).lt.0.)then
                yp(i)=yp(i)+L
            endif
            !do m=1,N!sprawdzenie, czy wartosci nie powtarzaja sie
            !    if(x(m).eq.xp(i))then
            !        goto 300
            !    endif
            !    if(y(m).eq.yp(i))then
            !        goto 400
            !    endif
            !enddo
            do j=0,N
                if(i.ne.j) then
                    dxp=abs(x(j)-xp(i))                    
                    dyp=abs(y(j)-yp(i))
                    if(dxp.lt.2.5.and.dyp.lt.2.5)then
                        dy=abs(y(j)-y(i))
                        dx=abs(x(j)-x(i))
                        if(dx.gt.L/2) then
                            dx=L-dx
                        endif
                        if(dxp.gt.L/2) then
                            dxp=L-dxp
                        endif
                        if(dy.gt.L/2) then
                            dy=L-dy
                        endif
                        if(dyp.gt.L/2) then
                            dyp=L-dyp
                        endif  
                        dr2=dx**2+dy**2
                        drp2=dxp**2+dyp**2
                        if(sqrt(drp2).le.rad)then
                            goto 500
                        endif
                        !(drp2**(-6)-drp2**(-3))-(dr2**(-6)-dr2**(-3))
                        dU=dU+4*((drp2**(-6)-drp2**(-3))-(dr2**(-6)-dr2**(-3)))!(drp2**(-6)-dr2**(-6)))
                        endif
                endif    
            enddo     
            w=min(1.0,exp(-dU/T))
            if(r.le.w) then!akceptacja
                x(i)=xp(i)
                y(i)=yp(i)
            endif
500         continue 
            !call cpu_time(finish)
            !write(11,*)(finish-start)
        enddo!koniec algorytmu Metropolisa
        if(k.gt.check .and. mod(k,time)==0) then
            do z=0,ale
                ra=z*drp
                call Prob(ra,drp,N,Rmax,x(:),y(:),q,a,L,rad)
                fun(z)=a+fun(z)
            enddo  
            licznik=licznik+1
        endif
    enddo
    do z=0,ale
        fun(z)=fun(z)/licznik
        write(12,*)z*drp,fun(z)
    enddo
end
    
    
FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) then
          idum=idum+IM
      endif
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
    end
subroutine Prob(r,da,N,Rmax,x,y,q,probability,L,rad)
    implicit none
    real probability,P
    integer N,i,j,m,rmax,L
    real r,da,x(N),y(N),q,deltaR
    real pi,dr,dx,dy,rad
    parameter (pi=3.141592653)
    P=0
    dr=da*5
    probability=P
    if(r.le.dr.and.r.le.rad) return
    do j=0, N
        M = 0
        do i=0, N
            if(i.ne.j) then
                dx=abs(x(i)-x(j))
                dy=abs(y(i)-y(j))
                if(dx.gt.L/2) then
                    dx=L-dx
                endif
                if(dy.gt.L/2) then
                    dy=L-dy
                endif
                deltaR = sqrt( dx**2 + dy**2 )
                if(deltaR.gt.r-dr.and.deltaR.le.r.and.deltaR.le.rad) then
                    M = M + 1
                endif
            endif
        enddo
          P = P + M/(2.0*PI*r*dr*q)
          !ring width: r-dr to r
    enddo
    Probability = P/(N+0.0)
    return
end