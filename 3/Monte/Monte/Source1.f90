Program Monte
    implicit none
    real*8 ran1
    real T
    integer L,Mcs
    integer i,j,k,B,kroki
    real pocz,kon,krok
    parameter (T=3,L=150,Mcs=200000)
    integer S(L,L) 
    integer nx(L),ny(L),px(L),py(L)
    real e(5),m,mk
    real Binder,Tcrit
    real im
    !----------------------------------------------------
    m=0
    !--------pliki----------------------------------------
    open(10,file="wartT3.out")
    !open(101,file="Tcrit.out")
    !---------------S------------------------------------
    do i=1,L
        do j=1,L
            S(i,j)=1
        enddo
    enddo
    !---wypelnienie tab nastepnych polozen---------------
    do i=1,L
        nx(i)=i+1
        ny(i)=i+1
        px(i)=i-1
        py(i)=i-1
    enddo
    nx(L)=1
    ny(L)=1
    py(1)=L
    px(1)=L
    !------------tablica wart e--------------------------
    do i=1,5
        e(i)=exp(-(-12+4*i)/T)
    enddo
    pocz=1.8d0
    kon=3.0d0
    krok=0.01
    kroki=(kon-pocz)/krok
    do im=pocz,kon,krok!do srednich
        if(im.ge.(2.05)) then
            krok=0.002
        endif
        if(im.gt.2.4) then
            krok=0.01
        endif
        !call MonT(S(:,:),L,im,nx(:),ny(:),px(:),py(:),Mcs)
    enddo
    call MonM(S(:,:),L,e(:),nx(:),ny(:),px(:),py(:),Mcs)
    !Tcrit=Binder(kroki)
    !write(101,*)Tcrit
    !call Energy(kroki)
    !call logi(kroki)
end
!-------------------------------------------------------------------
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
!----------------------------------------------------------------------
Subroutine writeM(A,L)
    integer A(L,L)
    integer a1,a2
    !--------------------------------------------
    open(unit=2, file='graphT1.7.txt', ACTION="write")
    !--------------------------------------------    
    do a1=1,L
        write(2, '(*(I3))')((A(a1,a2)) ,a2=1,L)
    enddo
    write(2,*) "------------------------------------------------------------------------------"!wypisywanie macierzy
    return
end
!------------------------------------------
Subroutine MonM(S,L,e,nx,ny,px,py,Mcs)
    implicit none
    real*8 ran1,r
    real T,kb
    integer L,Mcs,d
    integer i,j,k,n,o,p
    integer S(L,L),dEn
    integer nx(L),ny(L),px(L),py(L)
    real e(5),w,m
    !----------------------------------------------------
    d=-1
    dEn=0
    m=0
    !------------------------------------------------
    !------------algorytm MC-----------------------------
    do i=1, Mcs
        !call writeM(S,L)
        do j=1,L
            do k=1,L
                S(j,k)=-S(j,k)
                dEn = (-2)*S(j,k)*(S(nx(j),k)+S(px(j),k)+S(j,ny(k))+S(j,py(k)))
                r=ran1(d)
                select case(dEn)
                    case (-8)
                        n=1
                    case (-4)
                        n=2
                    case (0)
                        n=3
                    case (4)
                        n=4
                    case (8)
                        n=5
                endselect
                w=min(1.0,e(n))
                if(r.le.w) then
                    continue
                else
                    S(j,k)=-S(j,k)
                endif
            enddo      
        enddo!koniec algorytmu Metropolisa
        if(i.gt.30000 .and. mod(i,500)==0) then
            do o=1,L
                do p=1,L
                    m=m+S(o,p)
                enddo
            enddo
            m=m/(L**2)
            write(10,'(F20.15,I12)') m,i
            m=0
        endif

        !------------------------------------------------------!
        !                                                      !
        !       miejsce na obliczenie potrzebnych wartosci     !
        !                                                      !        
        !------------------------------------------------------!        
    enddo
    call writeM(S,L)
    end
!----------------------------

Subroutine MonT(S,L,T,nx,ny,px,py,Mcs)
    implicit none
    real*8 ran1,r
    real T,kb
    integer L,Mcs,d
    integer i,j,k,n,o,p
    integer S(L,L),dEn
    integer nx(L),ny(L),px(L),py(L)
    real e(5),w,m,ms,M2,x,H,Hsr,H2,C !tablica exp, magnet,srednia mag,Namag^2,podatnosc,Hamilt,Sredni,Ham^2,pojemnosc cieplna
    integer licznik
    real Ul,M4,Md
    real*8 beta,ni
    !----------------------------------------------------
    d=-1
    beta=0.125
    ni=1.0d0
    dEn=0
    m=0
    ms=0
    H=0
    H2=0
    Hsr=0
    m2=0
    Ul=0
    m4=0
    md=0
    licznik=0
    !--------pliki----------------------------------------
    open(11,file="wartL40.out")
    open(12,file="UL40.out")
    !----------------------------------------------------
    do i=1,L    
        do j=1,L !poczatkowy stan spinów
            S(i,j)=1
        enddo
    enddo
    !------------tablica wart e--------------------------
    do i=1,5
        e(i)=exp(-(-12+4*i)/T)
    enddo
    !------------algorytm MC-----------------------------
    do i=1, Mcs
        !call writeM(S,L)
        do j=1,L
            do k=1,L
                S(j,k)=-S(j,k)
                dEn = (-2)*S(j,k)*(S(nx(j),k)+S(px(j),k)+S(j,ny(k))+S(j,py(k)))
                r=ran1(d)
                select case(dEn)
                    case (-8)
                        n=1
                    case (-4)
                        n=2
                    case (0)
                        n=3
                    case (4)
                        n=4
                    case (8)
                        n=5
                endselect
                w=min(1.0,e(n))
                if(r.le.w) then
                    continue
                else
                    S(j,k)=-S(j,k)
                endif
            enddo      
        enddo
        if(i.ge.30000 .and. mod(i,500)==0) then
            do o=1,L
                do p=1,L
                    m=m+S(o,p)!liczenie magnetyzacji
                    H=H+S(o,p)*(S(nx(o),p)+S(px(o),p)+S(o,nx(p))+S(o,px(p)))
                enddo
            enddo
            m2=m2+m**2!do sredniej z kwadr M
            m4=m4+m**4!do sredniej z 4 potegi M
            H2=H2+H**2!do sredniej z kwadr energii
            md=md+abs(m)!do sredniej z M
            m=m/(L**2)!do sredniej magnetyzacji
            Hsr=Hsr+H
            ms=ms+abs(m)
            m=0
            H=0
            licznik=licznik+1
        endif
        
        !------------------------------------------------------!
        !                                                      !
        !       miejsce na obliczenie potrzebnych wartosci     !
        !                                                      !        
        !------------------------------------------------------!        
    enddo
    ms=ms/licznik
    m2=m2/licznik
    m4=m4/licznik
    md=md/licznik
    x=(m2-md**2)/(L*L*T)!podatnosc 
    Hsr=Hsr/licznik
    H2=H2/licznik
    C=(H2-Hsr**2)/(L*L*T*T)!pojemnosc
    Ul=1-(m4/(3*m2*m2))
    write(11,'(5F20.4)')ms,x,C,Hsr,T
    write(12,'(2F25.8)')Ul,T
    end   
!----------------------------
    function Binder(mcs)
    implicit none
    integer mcs,i,c,j,k,a,b
    real U100(mcs),U40(mcs),U10(mcs)
    real error
    real T(mcs),Tk
    real Binder,xl,xr,xm,z
    !--------pliki----------------------------------------
    open(11,file="UL10.out")
    open(12,file="UL40.out")
    open(15,file="Temps.out")
    !open(13,file="UL100.out")
    do i=1,mcs
        read(11,*,iostat=c) U10(i),T(i)
        read(12,*,iostat=c) U40(i)
        if(c.ne.0) then
            close(11)
            close(12)
            EXIT
        endif
    enddo
    
    !-------------------------------------------
    error=0.02
    xl=T(1)
    a=1
    b=mcs
    
    k=1
    !-------------------------------------------
100 xr=T(b)
    xl=T(a)
    xm=abs((xr+xl)/2.d0)
    write(15,*)(U10(k)/U40(k)),k,(1-(U10(a)/U40(a)))*(1-(U10(k)/U40(k)))
    do i=1,mcs-1
        j=i+1
        if(xm.ge.T(i).and.xm.le.T(j))then
            if((T(i)-xm).lt.(T(j)-xm))then
                xm=T(i)
                k=i
                goto 300
            else
                xm=T(j)
                k=j
                goto 300
            endif
        endif
    enddo
300 continue
    if((xr-xl).le.error) then 
        goto 200
    endif
    if((1-(U10(a)/U40(a)))*(1-(U10(k)/U40(k))).le.(0))then
        xr=xm
        b=k        
        goto 100
    else
        xl=xm
        a=k
        goto 100
    endif
200 continue
    Binder=(T(k)+T(k+1))/2
    return
    end
!---------------------------------------------------------------------
    subroutine Energy(kroki)
    integer kroki
    real T10(kroki),T40(kroki),T100(kroki/2)
    real E10(kroki),E40(kroki),E100(kroki/2)
    
    open(20,file="EnergyT10.out")
    open(21,file="EnergyT40.out")
    open(22,file="EnergyT100.out")
    open(3,file="wartL100.out")
    open(4,file="wartL40.out")
    open(5,file="wartL10.out")
    do i=1,int(kroki/2)
        read(3,*)x,x,x,E100(i),T100(i)
        write(22,*)-E100(i)/10000,T100(i)
    enddo
    do i=1,kroki
        read(5,*)x,x,x,E10(i),T10(i)
        read(4,*)x,x,x,E40(i),T40(i)
        write(20,*)-E10(i)/100,T10(i)
        write(21,*)-E40(i)/1600,T40(i)
    enddo

    return
    end
!--------------------------------------------------------
    subroutine logi(kroki)
    integer kroki
    real T10(kroki),T40(kroki),T100(kroki/2)
    real m10(kroki),m40(kroki),m100(kroki/2)
    real v,b,Tcrit
    !--------------------------------------
    v=1
    b=0.125
    open(20,file="logT10.out")
    open(21,file="logT40.out")
    open(22,file="logT100.out")
    open(3,file="wartL100.out")
    open(4,file="wartL40.out")
    open(5,file="wartL10.out")
    open(1,file="Tcrit.out")
    read(1,*)Tcrit
    do i=1,int(kroki/2)
        read(3,*)m100(i),x,x,x,T100(i)
        write(22,*)log(m100(i)*100**(b/v)),log(abs(1-(T100(i)/Tcrit))*100**(1/v))
    enddo
    do i=1,kroki
        read(5,*)m10(i),x,x,x,T10(i)
        read(4,*)m40(i),x,x,x,T40(i)
        write(20,*)log(m10(i)*10**(b/v)),log(abs(1-(T10(i)/Tcrit))*10**(1/v))
        write(21,*)log(m40(i)*40**(b/v)),log(abs(1-(T40(i)/Tcrit))*40**(1/v))
    enddo

    return
    end
        
    
    
    
