c-------------------------------------------------------------

c     Author: nsaini
c     Created: 2018-07-17

c-------------------------------------------------------------

      program Generate_xyzts
      implicit none

c     --------------------------------------------------------
      real*8, allocatable :: x1(:), y1(:), z1(:)
      real*8, allocatable :: xread(:), yread(:), zread(:)
      real*8 epsd, dist
      real*8 x,y,z
      
      integer stat
      integer myrank
      integer filetotal
      integer i,m,j
      integer index
      integer flag

      character*80 infile
      character*80 tempchar
c     --------------------------------------------------------      

      allocate(x1(100000))
      allocate(y1(100000))
      allocate(z1(100000))
      open(1,file = "myrank.txt",status = "old",action="read")
      open(10,file = "xyzts_big.dat",status = "replace")
      write(10,*) 5827, 1, 1.0000e-6, 5, 14, 1

      epsd = 1.0E-9
      index = 0
      do while(.true.)
         read(1,*,iostat=stat)myrank
         write(*,*)"processing rank # ",myrank
         if(stat .eq. 0)then
            if( myrank.ge.0 .and. myrank.le.9)then
               write(tempchar,'(i1.1)')myrank
            else if(myrank.ge.10 .and. myrank.le.99)then
               write(tempchar,'(i2.2)')myrank
            else if(myrank.ge.100 .and. myrank.le.999)then
               write(tempchar,'(i3.3)')myrank
            else if(myrank.ge.1000 .and. myrank.le.9999)then
               write(tempchar,'(i4.4)')myrank
            else if(myrank.ge.10000 .and. myrank.le.99999)then
               write(tempchar,'(i5.5)')myrank
            else if(myrank.ge.100000 .and. myrank.le.999999)then
               write(tempchar,'(i6.6)')myrank
            else if(myrank.ge.1000000 .and. myrank.le.9999999)then
               write(tempchar,'(i7.7)')myrank
            else if(myrank.ge.10000000 .and. myrank.le.99999999)then
               write(tempchar,'(i8.8)')myrank
            else if(myrank.ge.100000000 .and. myrank.le.999999999)then
               write(tempchar,'(i9.9)')myrank
            else if(myrank.ge.1000000000 .and. myrank.le.9999999999)then
               write(tempchar,'(i10.10)')myrank
            end if

            open(2, file = "bctin.dat."//trim(tempchar)//"")
            read(2,*)filetotal

            do i=1,filetotal
               read(2,*)x,y,z
               if(index .eq. 0)then
                  index = index + 1
                  write(10,*)x,y,z
                  x1(index) = x
                  y1(index) = y
                  z1(index) = z
               else
                  flag = 0
                  do j=1,index
                     dist = sqrt((x-x1(j))**2.0 + (y-y1(j))**2.0
     $                    + (z-z1(j))**2.0)
                     if(dist .lt. epsd)then
                        flag = 1
                     endif
                  enddo
                  if(flag .eq. 0)then
                     index = index + 1
                     write(10,*)x,y,z
                     x1(index) = x
                     y1(index) = y
                     z1(index) = z
                  endif                  
               endif
            enddo
            close(2)
         else
            exit
         endif
      enddo
      close(1)
      close(10)

      m=index
      
      allocate(xread(m))
      allocate(yread(m))
      allocate(zread(m))
      
      open(1990,file = "xyzts_big.dat",status="old",action="read")
      read(1990,*)
      do i=1,m
         read(1990,*)xread(i),yread(i),zread(i)
      enddo
      close (1990)
      
      open (1991,file="xyzts.dat",status="replace",action="write")
      write(1991,*)m,1, 1.0000e-6, 5, 14, 1
      do i=1,m
         write(1991,*)xread(i),yread(i),zread(i)
      enddo
      close(1991)

      write(*,*) 'Total number of nodes:',m
      deallocate(xread)
      deallocate(yread)
      deallocate(zread)

      deallocate(x1)
      deallocate(y1)
      deallocate(z1)
      end program
