c-------------------------------------------------------------

c     Author: nsaini
c     Created: 2018-07-05

c-------------------------------------------------------------
c     --------------------------------------------------------
C     Description of required input for merge.inp
C     PHASTA xyzts reocrd code: (same value as in solver.inp)
C     Number of time steps in a single file: (same value as in solver.inp)
C     Starting time step: (axiomatic)
C     Number of files to merge: (axiomatic)
c     --------------------------------------------------------      
      program merge
      
      implicit none
c     --------------------------------------------------------      
c     Variable declarations
      integer recordcode
      integer iBCTgap
      integer istart, Nfiles, istep
      integer np, nskip, nd1, nd2
C     integer inlength,
      integer outlength
      integer i, j, k, m
      integer lstep, jj, iphase
      integer icount
      
      real*8 tol
      real*8 varts(20)
      
      character*20 vartsFMT
      character*80 infile
      character*80 outfile
      character*80 tempchar
      character*20 fileFMT
c     --------------------------------------------------------

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

C     Input: xyzts.dat
      open(2, file = '../Process_bctinflow/xyzts.dat')      
      read(2, *) np, nskip, tol, nd1, nd2
      close(2)
c     --------------------------------------------------------      

c     The recordcode must match the value used in PHASTA, solver.inp
      
      if(recordcode .eq. 0)then !Default - Record everything
         vartsFMT='(2I8, I3, 19E15.7)'
         vartsFMT=trim(vartsFMT)
      else if(recordcode .eq. 1)then !Record only velocities
         vartsFMT='(2I8, I3, 4E15.7)'
         vartsFMT=trim(vartsFMT)
      else if(recordcode .eq. 2)then !Record only LS scalar
         vartsFMT='(2I8, I3, 1E15.7)'
         vartsFMT=trim(vartsFMT)
      else if(recordcode .eq. 3)then !Record velocities and LS scalar
         vartsFMT='(2I8, I3, 4E15.7)'
         vartsFMT=trim(vartsFMT)
      else
         write(*,*)"Specify record code between 0-3. Exiting"
         call exit(1)
      endif

      
c     Path to the directory  where all action takes place
           
      write(*,*) "Number of points per step = ",np

c$$$      if(recordcode .eq. 0)then
c$$$         inlength = 2*8+3+19*15
c$$$      elseif(recordcode .eq. 1)then
c$$$         inlength = 2*8+3+4*15
c$$$      elseif(recordcode .eq. 2)then
c$$$         inlength = 2*8+3+1*15
c$$$      elseif(recordcode .eq. 3)then
c$$$         inlength = 2*8+3+4*15
c$$$      endif
      outlength = 1*8 + 4*15
      
C     Open file for output
      call system ("mkdir mergedFiles")
!     Loop over the run files
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
         
         outfile = "mergedFiles/varts_run_"//trim(tempchar)//".dat"
         open(20, file=outfile, status="new", form = "formatted",
     &        action="write", recl=outlength, access="direct")

         infile = "varts/varts."//trim(tempchar)//".run.1.dat"
         open(10, file=infile, status="old", form ="formatted",
     &        action="read")

         icount = 0
         
         do j=1, iBCTgap
            write(*,*)"Reading step number: ",istep+j
            do k=1, np
               if(recordcode .eq. 0)then
                  read(10,vartsFMT,advance="no")lstep,jj,iphase
     $                 ,(varts(m), m=1,19)
               elseif(recordcode .eq. 1)then
                  read(10,vartsFMT,advance="no")lstep,jj,iphase
     $                 ,(varts(m), m=2,5)
               elseif(recordcode .eq. 2)then
                  read(10,vartsFMT,advance="no")lstep,jj,iphase
     $                 ,varts(15)
               elseif(recordcode .eq. 3)then
                  read(10,vartsFMT,advance="no")lstep,jj,iphase
     $                 ,(varts(m), m=2,5)
               endif
               icount = icount + 1
               write(20,"(1I8, 4E15.7)",REC=icount)
     &              lstep,(varts(m),m=2,5)
            enddo
         enddo
         close(10)
         close(20)
         
         
      enddo
        
      end program
    
