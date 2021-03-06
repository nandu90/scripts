c-------------------------------------------------------------

c     Author: nsaini
c     Created: 2018-07-05

c-------------------------------------------------------------

      program varts_to_BCT

      implicit none
c     --------------------------------------------------------
C     Variable declarations
      integer np, nskip, nd1, nd2
      integer recordcode, iBCTgap, istart, Nfiles
      integer istep
      integer ts1, ts2
      integer i, j, k, m
      integer lstep
      integer inlength

      real*8 tol
      real*8 ctime
      real*8 varts(20)
      real*8 oldtime
      real*8, allocatable :: tstart(:), tstop(:), totalt(:),deltat(:)
      real*8, allocatable :: x(:), y(:), z(:)

      character*80 fileFMT
      character*80 tempchar
      character*80 infile
      character*80 outfile
c     --------------------------------------------------------

c     --------------------------------------------------------
C     Input: xyzts.dat
      open(2, file = '../Process_bctinflow/xyzts.dat')      
      read(2, *) np, nskip, tol, nd1, nd2
      allocate(x(np))
      allocate(y(np))
      allocate(z(np))
      do i=1,np
         read(2,*)x(i), y(i), z(i)
      enddo
      close(2)
      write(*,*)"Done reading xyzts"

C     Input: merge.inp
      open(1, file = 'merge.inp')
      do i = 1, 6		! skip the header
         read(1,*)
      end do      
      read(1,*) recordcode
      read(1,*)
      read(1,*) iBCTgap		
      read(1,*) 
      read(1,*) istart
      read(1,*) 
      read(1,*) Nfiles

      write(*,*)
      write(*,*)"Merge.inp parameters:"
      write(*,*)"Recordcode: ",recordcode
      write(*,*)"iBCTgap: ",iBCTgap
      write(*,*)"Starting time step: ",istart
      write(*,*)"Number of file to process: ",Nfiles
      write(*,*)
      close(1)

C     Input: vbct.inp
      allocate(tstart(Nfiles))
      allocate(tstop(Nfiles))
      allocate(totalt(Nfiles))
      allocate(deltat(Nfiles))
      open(3, file="vbct.inp")
      do i = 1, 6		! skip the header
         read(3,*)
      end do
      do i=1,Nfiles
         read(3,*) tstart(i)    ! Range of time in the result
         read(3,*) tstop(i)
      enddo
!      read(3,*)
!      do i=1,Nfiles
!         read(3,*) totalt(i)
!      enddo
!      read(3,*)
!      do i=1,Nfiles
!         read(3,*) deltat(i)
!      enddo

      write(*,*)
      write(*,*) "vbct input parameters:"
      do i=1,Nfiles
         write(*,*)"Starting time value: ",tstart(i)
         write(*,*)"End time value: ",tstop(i)
!         write(*,*)"total t value: ",totalt(i)
!         write(*,*)"deltat t value: ",deltat(i)
         write(*,*)
      enddo
      close(3)
c     --------------------------------------------------------

      
      call system ("mkdir bctFiles")
C     Loop over the number of files
      do i=1,Nfiles
         istep = istart + iBCTgap*(i-1)

         if(istep .lt. 10)then
            fileFMT = "(I1)"
         elseif(istep .ge. 10 .and. istep .lt. 100)then
            fileFMT = "(I2)"
         elseif(istep .ge. 100 .and. istep .lt. 1000)then
            fileFMT = "(I3)"
         elseif(istep .ge. 1000 .and. istep .lt. 10000)then
            fileFMT = "(I4)"
         else
            fileFMT = "(I5)"
         endif

         write(*,*) "Processing from time step # ",istep  
         write(tempchar,fileFMT)istep

         inlength = 1*8+4*15
         infile = "mergedFiles/varts_run_"//trim(tempchar)//".dat"
         open(10, file=infile, status="old", form ="formatted",
     &        action="read", access="direct", recl=inlength)
         
         outfile ="bctFiles/bct."//trim(tempchar)//".dat"
         open(20, file=outfile, status="new", form = "formatted",
     &        action="write")

         

         write(20,"(1x, 2I10)")np, iBCTgap

         
         do j=1,np
            ctime = tstart(i)
            write(20,"(1x, 3E15.7, I10)")x(j),y(j),z(j),iBCTgap
            do k=1,iBCTgap
               read(10, "(1x, I8, 4E15.7)", REC=np*(k-1)+j)lstep,
     &              (varts(m),m=2,5)
!               if(k .eq. 1)then
!                  oldtime = varts(5)
!               else
               if(k .gt. 1)then
                  ctime = ctime + varts(5)
!                  oldtime = varts(5)
               endif
               if(ctime .ge. tstart(i) .and. ctime .lt. tstop(i))then
                  write(20,"(1x, 4E15.7, 10I10)")varts(2:4),ctime
               endif                 
               
            enddo
         enddo
         
         close(10)
         close(20)
      enddo
      
      
      

      deallocate(x)
      deallocate(y)
      deallocate(z)
      end program
      
