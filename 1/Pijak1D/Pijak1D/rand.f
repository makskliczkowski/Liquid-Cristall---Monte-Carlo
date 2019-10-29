      program RozkladN
      integer N,K,d,N2,N3,x
      parameter(N=10000,K=10000,N2=1000)
      real xn2,xn,r,var,ran1,logv(0:1000),logN(0:1000)
      open(10,file="wart.out")
      r=0
      d=-1
*     N3=N
      do a=0,N2
         xn=0
         xn2=0
         var=0
         N3=N*ran1(d)
*mnozenie razy losowa liczba tak o, mimo, ze wywolanie zabiera wi©cej t
*a)      N3=N3-1
             do i=1,K
                x=0
                do j=1,N3
                    r=ran1(d)
                    if(r.le.(0.5)) then
                        x=x+1
                    else
                        x=x-1
                    endif
                enddo
                xn2=xn2+x*x
                xn=xn+x
             enddo
          var=xn2/K - xn/K
          logv(a)=log(var)
          logn(a)=log(N3+0.0)
      enddo
      do i=0,N2
         write(10,*)logv(i),logn(i)
      enddo
      end
*-----------------------------------------------------------------------
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
      if (idum.lt.0) idum=idum+IM

      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      end

      
