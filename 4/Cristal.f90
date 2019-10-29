Program Monte
    implicit none
    real*8 ran1,r
    real pi,ne,no!wspolczynniki
    real P2
    integer r1,r2,r3,r4,r5,r6,r7,r8
    real n,w, E
    real dEn,Ptab(-180:180)
    integer Mcs,L,D
    integer i,j,da,k,licznik,o,p
    real ksi,dfi,psi,dU,zmiana!tablica katów, sta³a, sta³a, nowy k¹t,zmiana energii
    parameter (L=20,D=20,Mcs=230000,pi=3.141592653)
    parameter (ksi=20)
    parameter (no=1.5,ne=1.7,dfi=17)
    integer nx(L),px(L),ny(D),py(D),fi(L,D)
    !----------------------------------------------------
    da=-1
    E=0
    zmiana=0.0002
    n=0
    !--------pliki----------------------------------------
    open(10,file="n.out")
    !---------------S------------------------------------
    do i=1,L
        do j=1,D
            fi(i,j)=5
        enddo
    enddo
    !---wypelnienie tab nastepnych polozen---------------
    do i=1,L
        nx(i)=i+1
        px(i)=i-1
    enddo
    do i=1,D
        py(i)=i-1
        ny(i)=i+1
    enddo
    nx(L)=1
    ny(D)=D
    py(1)=1
    px(1)=L
    !------------tablica wart P2(x)--------------------------
    do i=-180,180
        Ptab(i)=P2(i*pi/180)
    enddo
    !----------------------------------------------------
    !------------algorytm MC-----------------------------
    do E=0.757,0.78,zmiana
        licznik=0
        do k=1, Mcs
            do i=1,L
                do j=2,D-1
                    r=ran1(da)
                    psi =int(fi(i,j) + (r-0.5)*Dfi)
                    if(psi.gt.90) then
                        psi = psi - 180
                    endif
                    if(psi.lt.-90) then
                        psi = psi + 180
                    endif
                    r1=int(fi(i,j)-fi(nx(i),j))
                    r2=int(fi(i,j)-fi(px(i),j))
                    r3=int(fi(i,j)-fi(i,py(j)))
                    r4=int(fi(i,j)-fi(i,ny(j)))
                    r5=-int(-psi+fi(i,py(j)))
                    r6=-int(-psi+fi(i,ny(j)))
                    r7=-int(-psi+fi(nx(i),j))
                    r8=-int(-psi+fi(px(i),j))
                    dEn = ksi*(Ptab(r1)+Ptab(r2)+Ptab(r3)+Ptab(r4)-Ptab(r5)-Ptab(r6)-Ptab(r7)-Ptab(r8))
                    dEn = dEn+(E**2)*((Ptab(int(90-fi(i,j)))-Ptab(90-int(psi))))
                    w=min(1.0,exp(-dEn))
                    !write(*,*) psi
                    if(ran1(da).le.w) then!akceptacja
                        fi(i,j)=psi
                        continue
                    else
                        continue
                    endif
                enddo      
            enddo!koniec algorytmu Metropolisa
            if(k.gt.30000 .and. mod(k,100)==0) then
                do o=1,L
                    do p=1,D
                        n=n+ne*no/sqrt(((sin(fi(o,p)*pi/180))**2)*(ne**2)+((cos(fi(o,p)*pi/180))**2)*(no**2))!efektywny n
                    enddo
                enddo
                licznik=licznik+1
            endif
        enddo
        n=n/(L*D)
        n=n/(licznik+0.0)
        write(10,*) n,E
        
        n=0
    enddo
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
Function P2(x)
real x,P2
P2=1.5*(cos(x))**2 - 0.5

return
end