      program bct_partition

      implicit none
     
      character*80 filenum1
      character*80 filenum
      character*80 filenum2
      character*80 ipath

      integer np, nts, myrank, stat, np_bct, i, j, k, m, l, p, flag
      real*8, allocatable :: u(:,:), v(:,:), w(:,:), t(:,:)
      integer, parameter :: maxp=35159000 
      integer, parameter :: bctnum=1 !Nadish - This has to changed everytime

      real*8, allocatable :: x(:), y(:), z(:)
      real*8, allocatable :: x1(:), y1(:), z1(:)
      real*8 x_p, y_p, z_p
      real*8, parameter :: epsd=1.0e-9
      integer isize
      real*8 dist

      open(1,file='path.dat')
         read(1,50) ipath
50       format(A80)
      close(1)
!        write(*,*) bctnum
      do p=1, bctnum
         if( p.ge.0 .and. p.le.9)then
            write(filenum2,'(i1.1)')p
         else if(p.ge.10 .and.p .le.99)then
            write(filenum2,'(i2.2)')p
         else if(p.ge.100 .and.p.le.999)then
            write(filenum2,'(i3.3)')p
         else if(p.ge.1000 .and.p.le.9999)then
            write(filenum2,'(i4.4)')p
         end if
         open(p+1, file=trim(ipath)//'/bct.'//trim(filenum2)//'.dat')
         read(p+1,*) np, nts
         if(p .eq. bctnum)then  !Nadish
            isize = np*(nts)
         else
            isize = np*(nts+1)  !Nadish - This is the original value
         endif
         write(*,*) p,isize,bctnum !Nadish
         write(*,*) np, nts         !Numer of points and number of time steps
         
         allocate(u(np,nts))
         allocate(v(np,nts))
         allocate(w(np,nts))
         allocate(t(np,nts))

         allocate(x1(np))
         allocate(y1(np))
         allocate(z1(np))
    
         
         do m=1, np             !loop over the points
            read(p+1,*) x1(m), y1(m), z1(m)
            !write(*,*) x1(m), y1(m), z1(m)
            do k =1,nts
               read(p+1,*) u(m,k), v(m,k), w(m,k), t(m,k)
               !write(*,*) u(m,k), v(m,k), w(m,k), t(m,k)
            enddo
         enddo
         close(p+1)

        
         
         open(bctnum+1,file=trim(ipath)//'/myrank.txt')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do while(.true.)
            read(bctnum+1,*,iostat=stat) myrank
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
               else if(myrank.ge.10000000 .and.myrank.le.99999999)then
                  write(filenum1,'(i8.8)')myrank
               else if(myrank.ge.100000000 .and.myrank.le.999999999)then
                  write(filenum1,'(i9.9)')myrank
               else if(myrank.ge.1000000000 .and.myrank.le.9999999999)then
                  write(filenum1,'(i10.10)')myrank
               end if
               
               open(100000+myrank, file = trim(ipath)//'/bctin.dat.'//trim(filenum1)//'')
               read(100000+myrank, *) np_bct

               allocate(x(np_bct))
               allocate(y(np_bct))
               allocate(z(np_bct))
               
               do l=1, np_bct
                  read(100000+myrank, *) x(l), y(l), z(l)
               enddo
               close(100000+myrank)
               
               open(200000+myrank,file=trim(ipath)//'/bct_input/bct.dat.
     &'//trim(filenum2)//'.'//trim(filenum1)//'', status='unknown')
    
               write(*,*) 'Processing processor', myrank,'
     &           of time series', p
               write(200000+myrank,13) np_bct, nts
               
               do i=1, np_bct
                  x_p=x(i)
                  y_p=y(i)
                  z_p=z(i)
!                  if(myrank .eq. 5986)then
!                     write(*,*)x_p,y_p,z_p
!                  endif
                  write(200000+myrank,12) x_p, y_p, z_p, nts
                  flag=0
                  do j=1,np
                     dist = sqrt((x_p-x1(j))**2 +(y_p-y1(j))**2+
     &                    (z_p-z1(j))**2 )
!                     if(myrank .eq. 5986 .and.
!     &                     j.gt.23195 .and. j .lt.23197)then
!                        write(*,*)j,x1(j),y1(j),z1(j),dist
!                     endif
                     if (flag.eq.0) then
                        if (dist .lt. epsd) then
                           do k=1, nts
                              write(200000+myrank,11)u(j,k),v(j,k)
     &                             ,w(j,k),t(j,k)
                           enddo
                           flag=1
                        else
                           continue
                        endif
                     endif
                  enddo         !end of search identical points
                  !if(myrank .eq. 5986)call exit(1)
               enddo            !end of loop over points on processor
               
               close(200000+myrank)
            else if(stat>0) then
               write(*,*) 'something is wrong with reading operation'
            else
               exit
            endif
            deallocate(x)
            deallocate(y)
            deallocate(z)
            
         enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      close(bctnum+1)
      !   enddo !end of current processor
       !  close(bctnum+1)
! close(p) !end of current bct series
       deallocate(u)
       deallocate(v)
       deallocate(w)
       deallocate(t)

       deallocate(x1)
       deallocate(y1)
       deallocate(z1)
      enddo
13       format(1x, 2I10)
12       format(1x, 3E15.7, I10)   !Coordinate and time step
11       format(1x, 4E15.7, 10I10) !Velocity and time
         end program
