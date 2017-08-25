!The values that need to be changed are marked with nash
!Compile blmesher in this directory but run it from the directory outside this one
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
!     Last change:  AG   14 Jun 2004    1:35 am
        module coordinate_array
        implicit none
        INTEGER numnp,numbnp
        INTEGER n1m,n2m,n3m,nsd
        INTEGER numel,nenl,ien,ienb
!*	integer i,j,k
        REAL*8 x,z,delta_x, dytemp
        dimension x(:,:)
        DIMENSION ien(:,:),ienb(:,:)
        allocatable :: ien,ienb
        allocatable :: x
        end module

      program mesh_general
        USE coordinate_array

!       ***********Variables*****************************************
!* ***********Reading part(from the PLOT3D format)*******************
	INTEGER i,j,k
!********Output arrays(coordinates and connectivities)*********************

	integer itemp
	dimension itemp(:)
	allocatable :: itemp
        INTEGER l,n,ng,iel,incr,nn
	integer inc1,inc2,inc3,ince1,ince2,ince3

!***********switch for reading or generating mesh**************************
        integer mesh_import,extract_surf

!*this particular program reads nodal points for the SHIP simulations 
!* from the Finite Difference model mesh and create FE mesh

!***********Formats*********************************************
114   format(1(i9,2x),12(e22.15,2x))
115   format(9(i9,2x))

!***************************************************************
	call gener_coor!generates nodal coordinate's
!********************************************************************
!***********writing nodal coordinates********************************!********************************************************************
	open(unit=96,file='p.crd',status='replace')
        do i=1,numnp
            write(96,114)i,x(i,1),x(i,2),x(i,3)
	enddo


!********Connectivities part****************************
c        write(*,*) 'enter nenl'
c        read(*,*) nenl
        nenl=8  !nash - for hex
c       nenl=4  !nash - for tet
!	***********************************************
	 write(*,*)'nenl=',nenl
!*  create connectivity array over here (not in the genien)


!* global node numbers of the model corner's 
	
         if(nenl.eq.8) nelPbrick=1
         if(nenl.eq.4) nelPbrick=6
	ALLOCATE(itemp(8))
      itemp(1)=1
      itemp(2)=2
      itemp(3)=n3m+2
      itemp(4)=n3m+1
      itemp(5)=n3m*n2m+1
      itemp(6)=n3m*n2m+2
      itemp(7)=n3m*(n2m+1)+2
      itemp(8)=n3m*(n2m+1)+1
	n=1
	ng=1
        ince1=nelPbrick
	inc1=1
	ince2=(n3m-1)*nelPbrick
	inc2=n3m
	ince3=(n3m-1)*(n2m-1)*nelPbrick
	inc3=n3m*(n2m)
        if (ince2 .eq. 0) ince2 = n3m-1
        if (ince3 .eq. 0) ince3 = (n3m-1) * (n2m-1)

!* **********Calculating number of elements****************************************
        numel=(n1m-1)*(n2m-1)*(n3m-1)
        if(nenl.eq.4) numel=numel*6
!***********generate connectivities************************************************
	allocate(ien(nenl,numel))
        do 1800 i = 0, n1m-2  
        do 1800 j = 0, n2m-2  
        do 1800 k = 0, n3m-2
          iel  = n + k * ince1 + j * ince2 + i * ince3
          if (iel .gt. numel) stop
          incr = k * inc1 + j * inc2 + i * inc3
          if(nenl.eq.8) then
             do 1700 nn = 1, nenl
                ien(nn,iel) = itemp(nn) + incr
 1700        continue
          else
c tet 1 is  1 4 8 3 
             ien(1,iel)=itemp(1)+incr
             ien(2,iel)=itemp(4)+incr
             ien(3,iel)=itemp(8)+incr
             ien(4,iel)=itemp(3)+incr
             iel=iel+1
c tet 2 is  3 8 1 2 
             ien(1,iel)=itemp(3)+incr
             ien(2,iel)=itemp(8)+incr
             ien(3,iel)=itemp(1)+incr
             ien(4,iel)=itemp(2)+incr
             iel=iel+1
c tet 3 is  8 3 2 7 
             ien(1,iel)=itemp(8)+incr
             ien(2,iel)=itemp(3)+incr
             ien(3,iel)=itemp(2)+incr
             ien(4,iel)=itemp(7)+incr
             iel=iel+1
c tet 4 is  1 8 5 2 
             ien(1,iel)=itemp(1)+incr
             ien(2,iel)=itemp(8)+incr
             ien(3,iel)=itemp(5)+incr
             ien(4,iel)=itemp(2)+incr
             iel=iel+1
c tet 5 is  8 5 2 7 
             ien(1,iel)=itemp(8)+incr
             ien(2,iel)=itemp(5)+incr
             ien(3,iel)=itemp(2)+incr
             ien(4,iel)=itemp(7)+incr
             iel=iel+1
c tet 6 is  5 2 7 6 
             ien(1,iel)=itemp(5)+incr
             ien(2,iel)=itemp(2)+incr
             ien(3,iel)=itemp(7)+incr
             ien(4,iel)=itemp(6)+incr
             iel=iel+1
          endif

             
1800    continue

        open(unit=97,file='p.cnn',status='replace')
	do i=1,numel !* loop over elements
 	write(97,115) i,(ien(nn,i),nn=1,nenl) !* write element number and connectivities
	enddo !* end of the loop of the elements

!*******writing mesh in the binary format for use in E2sms******************************
	call write_binary

!**********closing file's**************************************
	close(96)
	close(97)
! ********deallocating ****************************************
	deallocate(x)
! *********writing***************************************************

	write(*,*)'Finished'
	stop
      end


!**********************************************************************

        subroutine gener_coor
        use coordinate_array
!*****This routine for generate coordinate's manually
        INTEGER counterI
	real*8 :: parameter, pi = 3.1415926535E0
        REAL*8 xmax,ymax,zmax,counter
        REAL*8 xtemp,ytemp,ztemp
        REAL*8 dx,dy,dz,nlx,nly,nlz
	integer Nble
        DIMENSION xtemp(:),ytemp(:,:),ztemp(:)
        Dimension dy1a(:),dy1c(:),alfaa(:),alfac(:)
        Dimension delta(:),cf(:),utau(:)
        ALLOCATABLE :: xtemp,ytemp,ztemp,dy1a,dy1c,alfaa,alfac
        ALLOCATABLE :: cf,utau,delta
!*******Declaring size's of arrays and allocate arrays
c        write(*,*) 'simple (1) or geometric growth (2)'
c        read(*,*) simple

c        write(*,*) 'nptsx,nptsy,nptsz'
c        read(*,*) n1m,n2m,n3m
c      simple=0
c        n1m=45 !*#of nodal points in the x direction
c        n2m=48 !*#of nodal points in the y direction
c        n3m=32 !*#of nodal points in the z direction
c        simple=1
c        n1m=2 !*#of nodal points in the x direction
c        n2m=3 !*#of nodal points in the y direction
c        n3m=3 !*#of nodal points in the z direction


c         simple=2
c         n1m=90
c         n2m=50 
c         n3m=40 
         simple=1
         n1m=127                !nash - Number of elments in x-direction
         n2m=200                 !nash - Number of elements in y-direction
         n3m=100                 !nash - Number of elements in z-direction
         nlx=0.02546D0            !nash - x length
         nly=0.04D0           !nash - y length
         nlz=0.02D0           !nash - z length     



        numel=(n1m-1)*(n2m-1)*(n3m-1)
        numnp=n1m*n2m*n3m
	numbnp=n2m*n3m
	nsd=3

!**********************************************************
        ALLOCATE(x(numnp,nsd))
        ALLOCATE(xtemp(n1m))
        ALLOCATE(ytemp(n1m,n2m))
        ALLOCATE(ztemp(n3m))
        ALLOCATE(dy1a(n1m))
        ALLOCATE(dy1c(n1m))
        ALLOCATE(alfaa(n1m))
        ALLOCATE(alfac(n1m))
        ALLOCATE(cf(n1m))
        ALLOCATE(utau(n1m))
        ALLOCATE(delta(n1m))

c     calculate the inflow BL thickness deltain
c
      Rezeta=300.0d0              !set inlet Re_zeta
      Rexin=(Rezeta/0.036)**1.25   !calculate inlet Rex
      df = 2.*0.0292*Rexin**(-0.2)    !calculate dimensionless shearing stress
      xnu = 0.00036574
      ucl = 1.0
!      xin=Rexin*xnu/ucl
      xin = 0.0E0
      deltain=1.0  !(df/2./0.0227)**(-4)*xnu/ucl
      
c      write(*,*) 'inlet info'
      write(*,*) 'deltain=', deltain
c      write(*,*) 'Rex=',Rexin
      write(*,*) 'xnu=',xnu
      write(*,*) 'xin=',xin
c      print*, 'df=',df      

c Length of domain (in bl heights)

      xl=nlx*deltain
      yl=nly*deltain
      zl=nlz*deltain

       print*, 'xl=',xl 
       print*, 'yl=',yl
       print*, 'zl=',zl    
 
c     x coodinates
c     
      do i=1,n1m
         xtemp(i)=xin+(i-1)*xl/(n1m-1)
!         cf(i)=2*0.0292*(ucl*xtemp(i)/xnu)**(-0.2)
         utau(i)= 6.671E-02   !ucl*(cf(i)/2.)**0.5
!         delta(i)=(cf(i)/2./0.0227)**(-4)*xnu/ucl
!         print*, 'cf', cf(i), 'xtemp(i)', xtemp(i)
      enddo
      write(*,*) 'delx_plus=',utau(1),xnu
      write(*,*) (xtemp(4)-xtemp(3))*utau(3)/xnu

c
c  inside BL y stretching
c
      if(simple.ne.1) then
      do i=1,n1m
         alfaa(i)=3.0
      enddo
      eka=n2m*(1.0/2.0) ! No of nodes innner layer
      do i=1,n1m
         dy1a(i)=1.2*xnu/utau(1)          
 599      f=dy1a(i)*(alfaa(i)**eka-1.0)/(alfaa(i)-1.0)-delta(i)
         df=dy1a(i)*(eka*alfaa(i)**(eka-1.0)/(alfaa(i)-1.0)
     &        +(1.0-alfaa(i)*eka)
     &        /(1.0-alfaa(i))**2)
         alfaa(i)=alfaa(i)-f/df
         error=f/df
         if(error.gt.0.0000001) then 
            go to 599 
         endif
         write(*,*)'alfaa',i,alfaa(i)
      enddo
c     
c   outside BL y stretching
c
      do i=1,n1m
         alfac(i)=2.0
      enddo

      ekc=n2m-1-eka
      do i=1,n1m
         dy1c(i)=dy1a(i)*alfaa(i)**(eka-1.0)
 100     f=dy1c(i)*(alfac(i)**ekc-1.0)
     &        /(alfac(i)-1.0)-(yl-delta(i))
         df=dy1c(i)*(ekc*alfac(i)**(ekc-1.0)/(alfac(i)-1.0)
     &      +(1.0-alfac(i)**ekc)/(1.0-alfac(i))**2)
         alfac(i)=alfac(i)-f/df
         error=abs(f/df)
c         write(*,*)'i',i,'f',f,'df',df,'alfac',alfac(i),'error',error
         if(error.gt.0.0000001) then
            go to 100
         else
            write(*,*)'i=',i,'alfac=',alfac(i)
         endif
      enddo
      endif
c     
c     y coodinates
c
	Nble = 1  ! 10
      if(simple.eq.1) then
	write(*,*) 'Start denerating the Y mesh'
         do i=1,n1m
! Update this one:
	    ytemp(i,1) = -0.5D0*nly  !-1.0
            do j=2,Nble
	     ytemp(i,j) = ytemp(i, j-1) + 0.5E0*1.2**real(j-2)*xnu/utau(1) 
	    end do
            do j = n2m-Nble, n2m
             ytemp(i,j) = -1.0*ytemp(i, n2m - j + 1)
!            if (i.eq.1) write(*,*) 'j, y = ', j, ytemp(i,j+1), ytemp(i,j)
            end do
         dytemp = 2.0*(ytemp(i, n2m-Nble) - ytemp(i, Nble))
     1  /(n2m-2*Nble+1)
          do j = Nble+1, n2m-Nble
            ytemp(i,j) = ytemp(i, j-1) + dytemp
          enddo
         enddo
      else
      do i=1,n1m
         ytemp(i,1)=0.0
         ytemp(i,2)=dy1a(i)     
         do j=3,eka+1
            ytemp(i,j)=alfaa(i)**(j-2)*dy1a(i)+ytemp(i,j-1)
         enddo
         do j=eka+2,n2m
            ytemp(i,j)=alfac(i)**(j-eka-2)*dy1c(i)+ytemp(i,j-1)
         enddo
      enddo
      endif
       
      do j=1,n2m-1
         write(*,*) j, ytemp(1,j), (ytemp(1,j+1)-ytemp(1,j))*utau(1)/xnu
         write(501,*) 0.0, ytemp(1,j), 0.0
      enddo   
	write(*,*) 'Look for y coordinates in fort.501'

c
c     z coordinates
c
      do 30 i=1,n3m
         ztemp(i)=(n3m-i)*zl/(n3m-1)
 30   continue
      write(*,*)(ztemp(1)-ztemp(2))*utau(1)/xnu


!*********************************************
        counterI=0
        do i=1,n1m !loop by the i-index
	 do j=1,n2m !loop by the j-index
           do k=1,n3m !loop by the k-index
            counterI=counterI+1
	    x(counterI,1)=xtemp(i)
            x(counterI,2)=ytemp(i,j)
            x(counterI,3)=ztemp(k)
           enddo
         enddo
        enddo
        DEALLOCATE(xtemp)
        DEALLOCATE(ytemp)
        DEALLOCATE(ztemp)



        end subroutine


        subroutine write_binary 
        use coordinate_array
	integer i,j,numelb,nfaces,numflx
                  

!************open the binary file***************************************
        open (unit=99, file='geom.dat.1', status='replace',
     &                    form='unformatted', err=999)


!*.... write the geometry file
!*
	nsd=3
	numelb=0
	if(nenl.eq.8) nfaces=6
        if(nenl.eq.4) nfaces=4
	numflx=0
        write (*,*) 'numnp, nsd, numel, numelb, nenl, nfaces, numflx'
        write (*,*) numnp, nsd, numel, numelb, nenl, nfaces, numflx
        write (99) numnp, nsd, numel, numelb, nenl, nfaces, numflx

	write (99) ((x(i,j), i=1,numnp),j=1,nsd)
        write (99) ((ien(i,j),j=1,numel),i=1,nenl)
	write(99)
	write(99)
!*
!*.... close geometry file
!*
        close (99)
!*
!*.... end of file error handling
!*
        return
999     write(*,*) 'wrgeom: error in opening' 
        end
!***************************************************************





































