        program Generate_xyzts

        implicit none

        character*80 ipath

        integer i,j,k,myrank,stat
        integer,parameter :: maxp=1000000

        character*50 :: filenum
        character*50 :: filenum1
        real*8 x,y,z,w,t,r,n
        integer m
        real*8 x1(1:maxp),y1(1:maxp),z1(1:maxp)

        real*8, allocatable :: xread(:), yread(:), zread(:)

        open(1989,file='myrank.txt')
        open(1990,file='xyzts_small.dat')
        open(1991,file='xyzts_big.dat')

        write(1990,*) 5827, 1, 1.0000e-6, 5, 14, 1
        write(1991,*) 5827, 1, 1.0000e-6, 5, 14, 1
        m=0
        do while(.true.)
           read(1989,*,iostat=stat) myrank
           if(stat==0) then
                 if( myrank.ge.0 .and. myrank.le.9)then
                    write(filenum1,'(i1.1)')myrank
                 else if(myrank.ge.10 .and. myrank.le.99)then
                    write(filenum1,'(i2.2)')myrank
                 else if(myrank.ge.100 .and. myrank.le.999)then
                    write(filenum1,'(i3.3)')myrank
                 else if(myrank.ge.1000 .and. myrank.le.9999)then
                    write(filenum1,'(i4.4)')myrank
                 else if(myrank.ge.10000 .and. myrank.le.99999)then
                    write(filenum1,'(i5.5)')myrank
                 else if(myrank.ge.100000 .and. myrank.le.999999)then
                    write(filenum1,'(i6.6)')myrank
                 else if(myrank.ge.1000000 .and. myrank.le.9999999)then
                    write(filenum1,'(i7.7)')myrank
                 else if(myrank.ge.10000000 .and. myrank.le.99999999)then
                    write(filenum1,'(i8.8)')myrank
                 else if(myrank.ge.100000000 .and. myrank.le.999999999)then
                    write(filenum1,'(i9.9)')myrank
                 else if(myrank.ge.1000000000 .and. myrank.le.9999999999)then
                    write(filenum1,'(i10.10)')myrank
                 end if

       
        open(1+myrank, file = 'bctin.dat.'//trim(filenum1)//'')

        i=0

        read(1+myrank,*) n
        m=m+n

        do i=1,n
        read(1+myrank,*) x,y,z
        if (x.gt.0.01) then
        write(1991,*) x-0.0001,y,z
        else
        write(1991,*) x,y,z
        endif
        write(1990,*) x,y,z
        end do

        close(1+myrank)

        else if (stat>0) then
        write(*,*) 'Something is wrong with reading operation'
        else 
        exit
        endif
!        write(*,*) 'hello1'
        end do

!        write(*,*) 'hello2'
        close(1990)
        close(1991)
        write(*,*) 'Total number of nodes:',m

        allocate(xread(m))
        allocate(yread(m))
        allocate(zread(m))

        open(1990,file = "xyzts_big.dat",status="old",action="read")
        read(1990,*)
        do i=1,m
           read(1990,*)xread(i),yread(i),zread(i)
        enddo
        close (1990)

        open (1991,file="xyzts.dat",status="old",action="write")
        write(1991,*)m,1, 1.0000e-6, 5, 14, 1
        do i=1,m
           write(1991,*)xread(i),yread(i),zread(i)
        enddo
        close(1991)

        deallocate(xread)
        deallocate(yread)
        deallocate(zread)
        end program
        
