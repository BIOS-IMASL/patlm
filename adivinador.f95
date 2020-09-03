module guessing

  contains

    subroutine adivinador(guess,x,xg,constr,neq,pH,pHi)

      use parameters, only: iRpore,ihm2,dimr,dimz,xsolbulk,xb_I

      !#######################################################################
      !     Returns initial unknown guess (xg)
      !     guess=0 use bulk values (xsolbulk) as xg
      !     guess=1 use x as xg
      !     guess=2 read xg from file
      !#######################################################################
      implicit none
  
      integer,intent(in):: neq,guess
      real(8),intent(in):: x(neq),pH,pHi
      real(8),intent(out):: xg(neq),constr(neq)

      character(14):: auxcha,sim_id_in

      integer,parameter:: nguess=13

      integer:: i,nc_pore,nc_mem,nc_tot,nc_sol,nc_surf,neq_in

      real(8):: pH_in,csalt_in


      nc_pore=ihm2*iRpore 
      nc_mem=ihm2*dimr-nc_pore
      nc_tot=dimr*dimz
      nc_sol=nc_tot-nc_mem 
      nc_surf=(dimr-iRpore)+ihm2 ! # of surface cells

      ! remember that neq=nc_sol+nc_tot+2*nc_surf


      do i=1,neq
         xg(i)=x(i)

         constr(i)=0d0         ! constraint vector
      enddo



      do i=1,nc_sol
         !constr(i)=1d0
      enddo
      
      !do i=nc_sol+nc_tot+nc_surf+1,neq
      do i=nc_sol+nc_tot+1,nc_sol+nc_tot+nc_surf
         !constr(i)=1d0
      enddo




      if (guess==0) then

         write(*,'(/1x"using bulk guess...")')

         do i=1,nc_sol
            xg(i)=xsolbulk ! solvent volume fraction
         enddo
     
         do i=1,nc_tot
            xg(i+nc_sol)=0d0 ! electrostatic potential
         enddo
          
         do i=1,nc_surf
            xg(i+nc_tot+nc_sol)=1d0-xb_I ! lipid area fraction
            xg(i+nc_tot+nc_sol+nc_surf)=0d0 ! surface electrostatic potential
         enddo


      elseif (guess==2) then
     
         open(unit=nguess,file="guessx.in")


         read(nguess,*)auxcha,sim_id_in
         read(nguess,*)
         read(nguess,*)csalt_in,pH_in
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)
         read(nguess,*)neq_in
 
    
         if (neq_in /= neq) then
            write(*,*)"##### ERROR: wrong neq (L7) in file guessx.in!"
            stop
         endif

         
         do i=1,neq
            read(nguess,*)xg(i)
         enddo

         close(nguess)


         write(*,801)pH_in,csalt_in,sim_id_in
801      format(/1x"...Reading initial guess from file: guessx.in",1x/4x&
              &"with: pH:"1xf6.3,1x" -- csalt:"1xf6.4,1x" -- sim_id:"1xa14)

         
      else
         
         write(*,804)pHi
804      format(/1x"using previous solution pH:",1x,f6.3" as guess...")
     
      endif
  



  
    end subroutine adivinador




    !.......................................................................
    
    subroutine readsol(finum,x,neq,flname)

      use parameters      
      !#######################################################################
      !     Reads solution files
      !#######################################################################
      implicit none
  
      integer,intent(out):: neq
      real(8),intent(out),allocatable:: x(:)
      
      integer,intent(in):: finum

      character(50):: flname
      
      character(8):: date
      character(10):: time
      character(14):: auxcha,sim_id_in

      integer:: i

      integer,parameter:: solfile=21
      
      
      call date_and_time(DATE=date,TIME=time)
      
      sim_id=date(5:6)//"-"//date(7:8)//"_"//time(1:2)//"-"//time(3:4)//"-"//time(5:6)
      
      
      open(unit=solfile,file=flname)



      read(solfile,*)auxcha,sim_id_in
      read(solfile,*)
      read(solfile,*)csalt,pH
      read(solfile,*)dimr,deltar,dimz,deltaz
      read(solfile,*)iRpore,R_pore
      read(solfile,*)ihm2,h_mem
      read(solfile,*)xb_I,pK_I,z_I,a_I,a_N,a_l
      read(solfile,*)eps_m,lb_m,eps_w,lb_w
      read(solfile,*)temp,pKw
      read(solfile,*)vsol,vH,vOH,vpls,vmin
      read(solfile,*)zH,zOH,zpls,zmin
      read(solfile,*)rhopbulk
      read(solfile,*)
      
      rhopbulk=rhopbulk*Na
      
      vH=vH/vsol
      vOH=vOH/vsol
      vpls=vpls/vsol
      vmin=vmin/vsol

      a_I=a_I/a_l
      a_N=a_N/a_l


      read(solfile,*)neq

      allocate(x(neq))
      
      do i=1,neq
         read(solfile,*)x(i)
      enddo

      close(solfile)



      if (finum==1) then

         write(*,'(/"*****"43x"*****"/&
              &"***** Membrane with Pore - Peptide/Protein Sol. *****"&
              &   /"*****"43x"*****")')
         
      endif


      write(*,'(/1x,"Simulation ID: "a14)')sim_id

      write(*,'(/1x30("*")/)')
      write(*,'(1x"...Reading solutions from file: ",a50)')flname
      write(*,'(1x"...Input solution ID: ",a14,1x,"- please check!")')sim_id_in
      write(*,'(/1x"...post-processing,"/4x"pH:"1xf6.3/4x"csalt:"1xf6.4)')pH,csalt
  

    
    end subroutine readsol





  end module guessing
