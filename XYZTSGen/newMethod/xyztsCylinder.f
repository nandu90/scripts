c-------------------------------------------------------------

c     Author: nsaini
c     Created: 2018-07-16

c-------------------------------------------------------------

c     --------------------------------------------------------
C     The philosophy i am following at the moment is to generate
C     mesh points such that there is a point in each cell. Thus,
C     mesh resolution is the determining parameter for the density
C     of xyzts points
c     --------------------------------------------------------      
      program xyztsgen
      implicit none

c     --------------------------------------------------------
C     Variable Declarations
      real*8 meshres
      real*8 outerR
      real*8 r
      real*8 arc
      real*8 random
      real*8 pi 
      real*8 deltathetha
      real*8 iangle, angle
      
      real*8,allocatable,dimension(:) :: x, y, z
      real*8,allocatable,dimension(:) :: specifiedx

      integer narc
      integer index
      integer i,j
      integer maxIndex
      integer nplanes
      integer totalpoints
      integer,allocatable,dimension(:) :: npoints
      
c     --------------------------------------------------------

c     --------------------------------------------------------      
C     Input Section
      pi = 3.14159265359
                  
      meshres = 2.667E-05*3.0

      outerR = 0.0025

      maxIndex = 100000

      nplanes = 3

      allocate(x(maxIndex))
      allocate(y(maxIndex))
      allocate(z(maxIndex))
      allocate(npoints(nplanes))
      allocate(specifiedx(nplanes))
      
      specifiedx(1) = 6.75E-3
      specifiedx(2) = 7.5E-3
      specifiedx(3) = 8.25e-3
c     --------------------------------------------------------

      call init_random_seed()
      
      call random_number(random)
      
      

      index = 1
      
      do j=1, nplanes
         r = 0.0d0
c     --------------------------------------------------------
C     Place the first point at the origin
         y(index) = 0.0d0
         z(index) = 0.0d0
c     --------------------------------------------------------

c     --------------------------------------------------------
C     Now move out until you hit outer R
         index = index + 1
         do while (r .le. outerR)
            r = r + meshres
            arc = 2.0*pi*r

C     Increment in angle
            deltathetha = meshres/r
         
C     Max Number of points that can fit on the circumference
            narc = floor(arc/meshres)
         
C     Choose an angle at random
            call random_number(random)
            iangle = random*2.0*pi
            angle = iangle
         
            do i = 1,narc
               y(index) = r*cos(angle)
               z(index) = r*sin(angle)
               angle = angle + deltathetha
               index = index + 1
               if(index .gt. maxIndex) then
                  exit
               endif
            enddo

            if(index .gt. maxIndex) then
               write(*,*) "Requires a larger array for xyzts"
               exit
            endif
         enddo

         
      
c     --------------------------------------------------------
C     Assign values to x coordinate
         if(j .eq. 1)then
            npoints(j) = index-1
            x(1:npoints(j)) = specifiedx(j)
         else
            npoints(j) = index-1-sum(npoints(1:j-1),1)
            x(sum(npoints(1:j-1),1)+1:index-1) = specifiedx(j)
         endif

         write(*,*)"Points on plane: ",j,npoints(j)
c     --------------------------------------------------------         
      enddo

      

      totalpoints = sum(npoints,1)
      write(*,*)"Generated xyzts points: ",totalpoints
c     --------------------------------------------------------
C     Write to file
      
      open(10, file="xyzts.dat",status="replace", action="write")

      write(10,*)totalpoints, 1, 1.0E-6, 5, 14, 1
      do i=1,totalpoints
         write(10,*)x(i),y(i),z(i)
      enddo
      
      close(10)
c     --------------------------------------------------------      
      end program


      
      subroutine init_random_seed()
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed
          
      call random_seed()
      allocate(seed(1))
      
      call system_clock(COUNT=clock)
      
      seed = clock + 37
      call random_seed(PUT = seed)
          
      deallocate(seed)
      end
