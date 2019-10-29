Program Dyfusion
    implicit none
    integer L,N,K
    integer d,i,j,m!do petli itd.
    parameter (L=30,N=100,K=300) !K-liczba krokow, N-liczba czastek,L-dlugosc siatki
    integer nx(0:L),ny(0:L),px(0:L),py(0:L)!tablice sasiadow 
    integer x(0:N),y(0:N),A(L,L)!tablice polozen
    integer x0,y0!pomocnicze
    real przesun,difN !funkcje
    real*8 r,ran1,C
    external przesun
    !--------pliki----------------------------------------
    open(10,file="wart.out")
    open(11,file="wartL.out")
    d=-1
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
    !----------------zerowanie tablicy polozen-----------------------------
    do i=1,L
        do j=1,L
            A(i,j)=0
        enddo
    enddo
    !-------------wypelnianie poczatkowe tablicy polozen-------------------
    do i=1,N
18      continue
        r=ran1(d)
        x0=int(r*L)+1
        r=ran1(d)
        y0=int(r*L)+1
        if(A(x0,y0)==0) then
            x(i)=x0
            y(i)=y0
            A(x(i),y(i))=1
        else 
            goto 18
        endif
    enddo
    !----------dyfuzja dla roznych krokow---------
    do m=1,K,10
        write(10,'(F12.6,3I10)')  przesun(L,m,N,A(:,:),nx(:),ny(:),px(:),py(:),x(:),y(:)),m,L,N
        !call writeM(A(:,:),L)
    enddo
    !--------dyfuzja dla roznych macierzy----------
    do m=L**2-10,0,(-5)
        C=(m+0.0)/((L+0.0)**2)
        write(11,'(F12.6,F12.4,3I10)') difN(L,K,m,A(:,:),nx(:),ny(:),px(:),py(:)),C,m,L,N
    enddo      
    end
    !-||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||-!
    !----------------------------------------end--------------------------------------------------------
    
    function przesun(L,K,N,A,nx,ny,px,py,x,y)
    implicit none
    integer x0,y0
    integer L,K,N,d,B
    integer nx(0:L),ny(0:L),px(0:L),py(0:L) !x0,y0-pomocnicze
    integer dx(0:N),dy(0:N)    
    integer x(0:N),y(0:N),A(L,L)
    integer s,a1,a2
    real dR(0:N),przesun,sr,sr2!sr2 do sredniej
    integer i,j,z,o! do petli
    real*8 r,ran1
    !----------------------------------------------
    B=500
    d=-1
    sr2=0
    !---------------------------------------------
    do o=1,B
    !---------------------------------------------
        do i=1,N
            dR(i)=0
            dx(i)=0
            dy(i)=0      
        enddo
        !-----glowna petla-----------------------------------------
            do i=1,K
                do j=1,N        
                    r=ran1(d)
                    s=int(4*r+0.0)
                    x0=x(j)
                    y0=y(j)
                    if(s==0) then
                        y0=ny(y0)
                    endif
                    if(s==1) then
                        x0=nx(x0)
                    endif
                    if(s==2) then
                        y0=py(y0)
                    endif
                    if(s==3) then
                        x0=px(x0)
                    endif
                    if(A(x0,y0)==1) then
                        goto 19
                    else
                        A(x(j),y(j))=0
                        if(s==0) then 
                            dy(j)=dy(j)+1!(-y(j)+y0)
                        endif
                        if(s==1) then 
                            dx(j)=dx(j)+1!(-x(j)+x0)
                        endif
                        if(s==2) then 
                            dy(j)=dy(j)-1!(y(j)-y0)
                        endif
                        if(s==3) then 
                            dx(j)=dx(j)-1!(x(j)-x0)
                        endif
                        x(j)=x0
                        y(j)=y0
                        A(x(j),y(j))=1
                    endif
19                  continue
                enddo
            enddo
        sr=0
        do z=1,N
            dR(z)=(dx(z))**2+(dy(z))**2
            sr=sr+dR(z)        
        enddo
        sr=sr/N
        sr2=sr2+sr
    enddo
    sr2=sr2/B
    przesun=sr2/(4*K)
    return
    end
    !-------------------------------------------------
    Function difN(L,K,N,A,nx,ny,px,py)
    implicit none
    integer x0,y0
    integer L,K,N,d,B
    integer nx(0:L),ny(0:L),px(0:L),py(0:L) !x0,y0-pomocnicze
    integer dx(0:N),dy(0:N)    
    integer x(0:N),y(0:N),A(L,L)
    integer s,a1,a2
    real dR(0:N),difN,sr,sr2!sr2 do sredniej z calosci
    integer i,j,z,o! do petli
    real*8 r,ran1
    !--------------------------------------------    
    d=-1
    sr2=0
    B=500
    !----------------zerowanie tablicy polozen-----------------------------
    do i=1,L
        do j=1,L
            A(i,j)=0
        enddo
    enddo
    !-------------wypelnianie poczatkowe tablicy polozen od nowa, by zmieniac ny------------------------------
    do i=1,N
18      continue
        r=ran1(d)
        x0=int(r*L)+1
        r=ran1(d)
        y0=int(r*L)+1
        if(A(x0,y0)==0) then
            x(i)=x0
            y(i)=y0
            A(x(i),y(i))=1
        else 
            goto 18
        endif
    enddo
    !-----------------------petla-------------------------
    do o=1,B
    !---------------------------------------------
        do i=1,N
            dR(i)=0
            dx(i)=0
            dy(i)=0      
        enddo
    !-----glowna petla-----------------------------------------
            do i=1,K
                do j=1,N        
                    r=ran1(d)
                    s=int(4*r+0.0)
                    x0=x(j)
                    y0=y(j)
                    if(s==0) then
                        y0=ny(y0)
                    endif
                    if(s==1) then
                        x0=nx(x0)
                    endif
                    if(s==2) then
                        y0=py(y0)
                    endif
                    if(s==3) then
                        x0=px(x0)
                    endif
                    if(A(x0,y0)==1) then
                        goto 19
                    else
                        A(x(j),y(j))=0
                        if(s==0) then 
                            dy(j)=dy(j)+1!(-y(j)+y0)
                        endif
                        if(s==1) then 
                            dx(j)=dx(j)+1!(-x(j)+x0)
                        endif
                        if(s==2) then 
                            dy(j)=dy(j)-1!(y(j)-y0)
                        endif
                        if(s==3) then 
                            dx(j)=dx(j)-1!(x(j)-x0)
                        endif
                        x(j)=x0
                        y(j)=y0
                        A(x(j),y(j))=1
                    endif
19                  continue
                enddo
            enddo
        sr=0
        !liczenie sredniej====================================================
        do z=1,N
            dR(z)=(dx(z))**2+(dy(z))**2
            sr=sr+dR(z)        
        enddo
        sr=sr/N
        sr2=sr2+sr
    enddo
    !-----------------------koniec petli usredniaj¹cej--------------------------------------
    sr2=sr2/B
    difN=sr2/(4*K)
    return
    end
    !--------===================]]]]]]]]]]]]]]]]]]]]][[[[[[[[[[[[[[[[[[[[[[[[[[[==========================------------!
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
    !-----------------------------------------------------------------------------------------------------------------!
    Subroutine writeM(A,L)
    integer A(L,L)
    integer a1,a2
    !--------------------------------------------
    open(unit=2, file='graph1.txt', ACTION="write")
    !--------------------------------------------    
    do a1=1,L
        write(2, '(*(F3.0))')(real(A(a1,a2)) ,a2=1,L)
    enddo
    write(2,*) "------------------------------------------------------------------------------"!wypisywanie macierzy
    return
    end
    