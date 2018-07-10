c-------------------------------------------------------------

c     Author: nsaini
c     Created: 2018-07-09

c-------------------------------------------------------------

      program bct_partition

      implicit none
c     --------------------------------------------------------
C     Variable Declarations

      integer i, j, k, m
      integer recordcode
      integer iBCTgap
      integer istart
      integer Nfiles
      integer istep
      integer np, nts, isize
      integer myrank
      integer stat
      integer np_bct
      integer flag

      real*8, allocatable :: u(:,:), v(:,:), w(:,:), t(:,:)
      real*8, allocatable :: x(:), y(:), z(:)
      real*8, allocatable :: x1(:), y1(:), z1(:)
      real*8 x_p, y_p, z_p
      real*8 dist
      real*8 epsd
      
      character*80 infile
      character*80 fileFMT
      character*80 tempchar
      character*80 tempchar2
      character*80 outfile

      epsd = 1.0E-9
      
c     --------------------------------------------------------
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

c     --------------------------------------------------------

      call system("mkdir ../bctInput")
      
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

         write(tempchar,fileFMT)istep

         call system("mkdir ../bctInput/"//trim(tempchar)//"")
         
         infile = "bctFiles/bct."//trim(tempchar)//".dat"
         open(10, file=infile, status="old", form ="formatted",
     &        action="read")

         read(10,*)np,nts
         
         isize = np*nts

         write(*,*)"Processing time step number: ",istep
         write(*,*)"np, nts, isize: ",np, nts, isize

         allocate(u(np,nts))
         allocate(v(np,nts))
         allocate(w(np,nts))
         allocate(t(np,nts))

         allocate(x1(np))
         allocate(y1(np))
         allocate(z1(np))

         do j=1, np             !loop over the points
            read(10,*) x1(j), y1(j), z1(j)
            do k =1,nts
               read(10,*) u(j,k), v(j,k), w(j,k), t(j,k)
            enddo
         enddo
         close(10)

         open(20, file = "../Process_bctinflow/myrank.txt")

         do while(.true.)
            tempchar = ""
            infile = ""
            read(20,*,iostat=stat)myrank
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
               else if(myrank.ge.10000000 .and.myrank.le.99999999)then
                  write(tempchar,'(i8.8)')myrank
               else if(myrank.ge.100000000 .and.myrank.le.999999999)then
                  write(tempchar,'(i9.9)')myrank
               end if

               infile =
     &              "../Process_bctinflow/bctin.dat."//
     &              trim(tempchar)//''
               
               open(30,file = infile,action="read")

               read(30,*)np_bct

               allocate(x(np_bct))
               allocate(y(np_bct))
               allocate(z(np_bct))

               do j=1,np_bct
                  read(30,*)x(j), y(j), z(j)
               enddo              
               close(30)

               tempchar2 = ""
               write(tempchar2,fileFMT)istep
               outfile =
     &              "../bctInput/"//trim(tempchar2)//"/bct.dat."
     $              //trim(tempchar)//""

               open(40,file = outfile, status = "unknown")

               write(40,"(1x, 2I10)")np_bct, nts

               do j=1, np_bct
                  x_p = x(j)
                  y_p = y(j)
                  z_p = z(j)

                  write(40, "(1x, 3E15.7, I10)")x_p, y_p, z_p, nts
                  flag = 0
                  do k=1,np
                     dist = sqrt((x_p-x1(k))**2.0 + (y_p-y1(k))**2.0 +
     $                    (z_p-z1(k))**2.0)
                     
                     if(flag .eq. 0)then
                        if(dist .lt. epsd)then
                           do m=1,nts
                              write(40, "(1x, 4E15.7, 10I10)")
     $                             u(k,m),v(k,m),w(k,m),t(k,m)
                           enddo
                           flag = 1
                        else
                           continue
                        endif
                     endif
                  enddo
                  
               enddo

               close(40)
               
            elseif(stat .gt. 0)then
               write(*,*) 'something is wrong with reading operation'
            else
               exit
            endif

            deallocate(x)
            deallocate(y)
            deallocate(z)
            
         enddo
         
         close(20)

         deallocate(u)
         deallocate(v)
         deallocate(w)
         deallocate(t)
         
         deallocate(x1)
         deallocate(y1)
         deallocate(z1)
      enddo
      
      
      end program
