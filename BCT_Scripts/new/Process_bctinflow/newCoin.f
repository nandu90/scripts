c-------------------------------------------------------------

c     Author: nsaini
c     Created: 2018-07-17

c-------------------------------------------------------------

      program Convert_Inflow
      implicit none

c     --------------------------------------------------------
C     Variable declarations
      real*8, allocatable :: x1(:), y1(:), z1(:)
      real*8 x,y,z
      real*8 dist
      real*8 epsd
      
      integer myrankstat, bctstat
      integer myrank
      integer index
      integer i
      integer flag
      
      character*80 infile
      character*80 outfile
      character*80 tempchar
c     --------------------------------------------------------

      epsd = 1.0E-9

      allocate(x1(10000))
      allocate(y1(10000))
      allocate(z1(10000))
      
      open(1,file = "myrank.txt")
      do while(.true.)
         read(1,*,iostat=myrankstat)myrank
         if(myrankstat.eq.0)then
            infile = ""
            tempchar = ""
            outfile = ""
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
            endif

            infile = "bct/bct.dat."//trim(tempchar)//""

            open(2, file=infile,action="read",status = "old")
            index = 0
            
            do
               read(2,*,iostat=bctstat)myrank,x,y,z
               if(bctstat .eq. 0)then
                  if(index .eq. 0)then
                     index = 1
                     x1(index) = x
                     y1(index) = y
                     z1(index) = z
                  else
                     flag = 0
                     do i=1,index
                        dist = sqrt((x-x1(i))**2.0 + (y-y1(i))**2.0 +
     $                       (z-z1(i))**2.0)
                        if(dist .lt. epsd)then
                           flag = 1
                        endif
                     enddo
                     if(flag .eq. 0)then
                        index = index + 1
                        x1(index) = x
                        y1(index) = y
                        z1(index) = z
                     endif
                     
                  endif
               else
                  close(2)
                  write(*,*) index,
     $            'unique points loaded from the original PHASTA file'
                  exit
               endif               
            enddo

            outfile = "bctin.dat."//trim(tempchar)//""

            open(3,file=outfile,status="replace",action="write")
            write(3,*)index

            do i=1,index
               write(3,"(10E16.7)")x1(i),y1(i),z1(i)
            enddo
            close(3)
         else
            exit
         endif
      
      enddo
      close(1)
      deallocate(x1)
      deallocate(y1)
      deallocate(z1)
      end program
