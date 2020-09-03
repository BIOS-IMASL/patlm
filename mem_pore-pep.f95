
program mem_pore_pep
  
  use parameters
  use proteina
  use printout
  use guessing
  use call_solver

  implicit none

  integer:: guessx,ier,dimkrlv,mxiter,i,neq,nsol

  integer,parameter:: nsolfi=22      
    
  real(8):: pHmin,pHmax,dpHmin,dpHmax,dpHstp
  real(8):: dpH,pH0,pHvec(2),dpHvec(3)
  real(8):: errtol(2),fnorm
  real(8),allocatable:: x(:),xg(:),f(:),constr(:)
  
  real:: t_ini,t_fin,dt

  character(50):: solfiname

  logical:: loopH,goodx




  

  !..........read parameters

  call read_param(errtol,dimkrlv,mxiter,guessx,pHvec,dpHvec)


  

  !............................................................
  !..... solutions given in "sol-cases.in" file - only calculate free energy and stuff
  !............................................................
  ! Formato de sol-cases.in
  ! Entero = numero de soluciones a leer
  ! Nombres de los archivos de solucion uno debajo de otro
  ! Notar que el archivo de solucion contine los datos del archivo de entrada con que se corrio. 
  ! Ingonara lo que lee de param.in
  
  if (guessx==3) then
     
     open(nsolfi,file="sol-cases.in")
     
     read(nsolfi,*)nsol
              
     do i=1,nsol
                
        read(nsolfi,*)solfiname
           
        call readsol(i,x,neq,solfiname)

        if (i==1) call read_peptide(guessx,deltaz,deltar)
               
        ! Genera las cantidades de bulk. Ver notas, solo calcula parametros a partir de los inputs
        call bulkimeter(vsol,vH,zH,vOH,zOH,vpls,zpls,vmin,zmin,pKw,pH,csalt,&  
             &xsolbulk,xHbulk,xOHbulk,xplsbulk,xminbulk,rhopbulk,&
             &vaa,zaa,pKaa,naac,cnaa,ntyi)          
        
	! Procesa los datos.
        call freen(x,neq) ! Ahora esta subrutina no es funcional
        
     enddo

     close(nsolfi)

     write(*,'(/1x5("*")" Ciao!"/)')

     stop  
     
  endif

  call cpu_time(t_ini)

  pHmin=pHvec(1)
  pHmax=pHvec(2)
  
  dpHmax=dpHvec(1)
  dpHmin=dpHvec(2)
  dpHstp=dpHvec(3)
  
  if (pHmin>pHmax) dpHstp=-dpHstp


  !..........read/generate protein and rotate it
    

   ! Lee archivos de entrada del peptido y lo genera con RIS     
  call read_peptide(guessx,deltaz,deltar)
  
 
  call get_neq(neq)  ! calcula numero de ecuaciones a resolver
     
  allocate(f(neq))          ! f is dummy here
  allocate(x(neq))
  allocate(xg(neq))
  allocate(constr(neq))
  
  
  call cpu_time(t_fin)
  dt=t_fin-t_ini
  
          
  call print_screen(sim_id,neq,dt)
  


  !................
  !..... loop on pH
  !................
  

  ier=0

  loopH=.true.
  pH=pHmin
  pH0=0d0

  if (pHmin<=pHmax) then
     dpH=dpHmax
  else
     dpH=-dpHmax
  endif

  
  do while (loopH)        
	
     ! Genera las cantidades de bulk. Ver notas, solo calcula parametros a partir de los inputs	  
     call bulkimeter(vsol,vH,zH,vOH,zOH,vpls,zpls,vmin,zmin,pKw,pH,csalt,&  
          &xsolbulk,xHbulk,xOHbulk,xplsbulk,xminbulk,rhopbulk,&
          &vaa,zaa,pKaa,naac,cnaa,ntyi)         
     
     write(*,1401)pH,csalt
1401 format(/1x30("*")//1x"...starting calculation with:"/&
          &4x"pH:"1xf6.3/4x"csalt:"1xf6.4)
     
     ! Guess incial para el solver	          
     call adivinador(guessx,x,xg,constr,neq,pH,pH0) ! guess     
     guessx=1    ! A partir de aqui se toma la solucion al ph anterior para el nuevo calculo 
     
     call cpu_time(t_ini)
    
     ! subrutina que llama al solver, devuelve ier=0 si encontro solucion             
     ier=0

     !write(*,*)"STOP"
     !stop
     call llamador(x,xg,constr,neq,errtol,dimkrlv,mxiter,ier,fnorm) ! solve

     call cpu_time(t_fin)
     dt=t_fin-t_ini
     
     ! Imprime las soluciones (big files) que luego se pueden usar como inputs
     call print_xs(ier,fnorm,goodx,neq,x,constr,dt)
	          
     if (ier==0.and.goodx) then
        
        ! Se llama  a esta rutina si se quiere procesar la solucion 
        call freen(x,neq)

        pH0=pH
        pH=pH+dpH
        
     else
        
        write(*,'(/1x"...no solution found for pH:",1x,f6.3)')pH

        dpH=dpH-dpHstp
        pH=pH0+dpH

        x=xg

        if (pH0==0d0) exit ! firt pH calculation failed
        
     endif
     
     ! revisa que este todo bien con el loop en pH
     call StayOn_pHloop(pH,pH0,dpH,dpHmin,pHmin,pHmax,loopH)

     if (.not.loopH) exit
     
  enddo
  

  !..........
  

  write(*,*)
  
  
  
  
  
  
contains
  
  
  
    
  
  subroutine stayOn_pHloop(pH,pH0,dpH,dpHmin,pHmin,pHmax,loopH)
    implicit none
    
    real(8),intent(in):: pH0,dpH,dpHmin,pHmin,pHmax
    real(8),intent(inout):: pH
    
    logical,intent(out):: loopH
    
    
    
    loopH=.true.
    
    if (abs(dpH)<dpHmin) loopH=.false.
    
    if (dpH==0d0) then
       loopH=.false. ! do not do case again
       
       write(*,*)
       write(*,*)" ...pH_step=0 reached"

    elseif (dpH<0d0.and.pHmin<pHmax) then
       loopH=.false.       

       write(*,*)
       write(*,*)" ...minimum |pH_step| reached"

    elseif (dpH>=0.d0.and.pHmin>pHmax) then
       loopH=.false.

       write(*,*)
       write(*,*)" ...minimum |pH_step| reached"
    endif
    
    
    if (pH>pHmax.and.pHmin<=pHmax) then
       
       if (pH0==pHmax) then
          loopH=.false.
       else
          pH=pHmax
       endif
    endif
               
    if (pH<pHmax.and.pHmin>pHmax) then
                  
       if (pH0==pHmax) then
          loopH=.false.
       else
          pH=pHmax
       endif
    endif
    
    

  end subroutine stayOn_pHloop
  




end program mem_pore_pep
