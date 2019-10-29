Program Dyfusion
    integer x0,y0,uloz
    integer L,N,K,d,z!z do niesprawdzania ifów za kazdym razem, tylko 1
    parameter(L=30,N=3) !K-liczba krokow, N-liczba czastek,L-dlugosc siatki
    integer dx(0:N),dy(0:N),nx(0:L),ny(0:L),px(0:L),py(0:L) !x0,y0-pomocnicze
    integer x(0:N),y(0:N),A(L,L)
    integer s,a1,a2
    real dR(0:N)
    real*8 r,ran1,dyf,sr
    open(10,file="wart.out")
    open(unit=2, file='graph1.txt', ACTION="write")
    d=-1
    !---wypelnienie tab nastepnych---------------
    sr=0!srednia =0
    K=10000
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
    !---------------------------------------------
    !----------wypelnienie tablic polozen---------
    do i=1,L
        do j=1,L
            A(i,j)=0
        enddo
    enddo
    do i=1,N
        dx(i)=0
        dy(i)=0!by uzupelnic delty
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
    !----------------------------------------------
    do i=1,K
        do a1=1,L
             write(2, '(*(F3.0))')(real(A(a1,a2)) ,a2=1,L)
        enddo
        write(2,*) "------------------------------------------------------------------------------"!wypisywanie macierzy
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
                    dy(j)=dy(j)+(y(j)-y0)
                endif
                if(s==1) then 
                    dx(j)=dx(j)+(x(j)-x0)
                endif
                if(s==2) then 
                    dy(j)=dy(j)+(y(j)-y0)
                endif
                if(s==3) then 
                    dx(j)=dx(j)+(x(j)-x0)
                endif
                x(j)=x0
                y(j)=y0
                A(x(j),y(j))=1
            endif
19          continue
        enddo
    enddo
    !-------------------------------------------------
    sr=0
    do z=1,N
        dR(z)=dx(z)**2+dy(z)**2
        sr=sr+dR(z)        
    enddo
    sr=sr/N
    dyf=sr/4
    dyf=dyf/K
    !------------------------------------------------
    write(10,'(4F12.6)') dyf,K,N,L      
    
    
    
    end
    !--------------------------------------------------------------------------------------------------------------------------------- 
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
    !-------------------------------------------------------------------------------------------------------------------------------------
    