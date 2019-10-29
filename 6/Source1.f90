program liczenie
    implicit none
    integer i
    real a,b
    open(11,file="088.dat")
    open(10,file="prob.dat")
    do i=1,1600
        read(10,*)a,b
        b=b/1188645
        write(11,*)a,b
    enddo
    end
    