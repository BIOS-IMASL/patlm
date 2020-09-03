module parameters

  real(8),parameter :: Na=0.602214,pi=acos(-1.0),cbjrrm=16710.075
  
  real(8):: zH,zOH,zpls,zmin
  real(8):: vsol,vH,vOH,vpls,vmin,rmin,rpls
  real(8):: temp,pKw,pOH,csalt,pH
  real(8):: xsolbulk,xHbulk,xOHbulk,xplsbulk,xminbulk,rhopbulk
  real(8):: eps_w,lb_w
  real(8):: xb_I,pK_I,z_I,a_I,a_N,a_l,eps_m,lb_m,epsunt
  real(8):: deltar,deltaz,R_pore,h_mem

  integer:: dimr,dimz,iRpore,ihm2,numeq

  character(14):: sim_id
  
contains

  subroutine read_param(errtol,dimkrlv,mxiter,guessx,pHvec,dpHvec) 
    
    !#######################################################################
    !     Reads (sets) all parameters of the calculation
    !#######################################################################

    implicit none
    integer:: i,guessx

    integer:: dimkrlv,mxiter
    real(8):: errtol(2)    
    real(8):: pHvec(2),dpHvec(3)
    
    character(10):: stread,title,time
    character(8):: date
    logical:: loop
    
    integer,parameter:: ninprm=1

    

    
    
    call date_and_time(DATE=date,TIME=time)
    
    sim_id=date(5:6)//"-"//date(7:8)//"_"//time(1:2)//"-"//time(3:4)//"-"//time(5:6)
    
    
    
    !...................
    !.....Solution......
    !...................
    
    open(unit=ninprm,file="param.in")
    
    title="Solution"

    loop=.true.
    do while(loop)
       read(ninprm,*,end=370)stread
       if (stread==title) loop=.false.
    enddo

370 if (loop) then
       write(*,'(/1x"Error in input file - No Solution found in param.in!"/)')
       stop
    endif
    
    !.....
    
    read(ninprm,*)
    read(ninprm,*)temp
    read(ninprm,*)
    read(ninprm,*)pKw
    
    !.....Volumes.....
    read(ninprm,*)
    read(ninprm,*)vsol
    read(ninprm,*)vH
    read(ninprm,*)vOH
    read(ninprm,*)rpls
    read(ninprm,*)rmin
    
    
    !.....Charges.....
    read(ninprm,*)
    read(ninprm,*)zH
    read(ninprm,*)zOH
    read(ninprm,*)zpls
    read(ninprm,*)zmin

    !.....Electrostatics.....
    read(ninprm,*)
    read(ninprm,*)eps_w



    !..................
    !.....Membrane.....
    !..................
    
    rewind(ninprm)

    title="Membrane"

    loop=.true.
    do while(loop)
       read(ninprm,*,end=371)stread
       if (stread==title) loop=.false.
    enddo
    
371 if (loop) then
       write(*,'(/1x"Error in input file - No Membrane found!"/)')
       stop
    endif
    
         
    read(ninprm,*)
    read(ninprm,*)xb_I
    read(ninprm,*)pK_I
    read(ninprm,*)z_I
    read(ninprm,*)
    read(ninprm,*)a_I
    read(ninprm,*)a_N
    read(ninprm,*)a_l
    read(ninprm,*)
    read(ninprm,*)eps_m
    read(ninprm,*)
    read(ninprm,*)h_mem


    !.................
    !.....Peptide.....
    !.................  
    
    rewind(ninprm)
    
    title="Peptide"
    
    loop=.true.
    do while(loop)
       read(ninprm,*,end=372)stread
       if (stread==title) loop=.false.
    enddo
    
372 if (loop) then
       write(*,'(/1x"Error in input file - No Peptide found!"/)')
       stop
    endif
         
         
    read(ninprm,*)
    read(ninprm,*)rhopbulk
    
    rhopbulk=rhopbulk*Na





    !...................
    !.....Resolution....
    !...................
    
    rewind(ninprm)

    title="Resolution"

    loop=.true.
    do while(loop)
       read(ninprm,*,end=373)stread
       if (stread==title) loop=.false.
    enddo

373 if (loop) then
       write(*,'(/1x"Error in input file - No Resolution found!"/)')
       stop
    endif

    read(ninprm,*)
    read(ninprm,*)guessx
    read(ninprm,*)
    read(ninprm,*)deltar,deltaz
    read(ninprm,*)dimr,dimz
    read(ninprm,*)
    read(ninprm,*)dimkrlv  ! maximum Krylov subspace dimesion, =0 means dimkrlv=neq
    read(ninprm,*)mxiter   ! max # of iterations, default=200
    read(ninprm,*)
    read(ninprm,*)errtol(1)
    read(ninprm,*)errtol(2)


    !...................
    !.....Environment....
    !...................
    
    rewind(ninprm)
    
    title="Environment"
    
    loop=.true.
    do while(loop)
       read(ninprm,*,end=374)stread
       if (stread==title) loop=.false.
    enddo
    
374 if (loop) then
       write(*,'(/1x"Error in input file - No Environment found!"/)')
       stop
    endif
    
    
    read(ninprm,*)         
    read(ninprm,*)pHvec(1),pHvec(2) ! pHs
    read(ninprm,*)dpHvec(1),dpHvec(2),dpHvec(3)
    
    read(ninprm,*)csalt







    !..............
    !.....Pore.....
    !..............
    
    rewind(ninprm)

    title="Pore"
    
    loop=.true.
    do while(loop)
       read(ninprm,*,end=375)stread
       if (stread==title) loop=.false.
    enddo
    
375 if (loop) then
       write(*,'(/1x"Error in input file - No Pore found!"/)')
       stop
    endif
    
         
    read(ninprm,*)
    read(ninprm,*)R_pore





    close(ninprm)
         





    !..........
    !..........
    !..........



    ! Pone cantidades en las unidades del programa: Volumen=volumen del solvente, excepto el del solvente que es nm^33, densidad= 1/nm^3, energias en kT, 
    iRpore=nint(R_pore/deltar)        ! dimensiones del poro en celdas
    R_pore=dble(iRpore)*deltar        ! corrige para que sea conmensurado con deltar


    ihm2=nint(h_mem/2d0/deltaz)     ! Lo mismo que antes para el espesor de la membrana
    h_mem=dble(ihm2)*deltaz*2d0


    a_I=a_I/a_l   ! divide por area de referencia
    a_N=a_N/a_l
      
    
    vpls=(4./3.)*pi*rpls**3/vsol    ! volumen sal
    vmin=(4./3.)*pi*rmin**3/vsol
    
    lb_w=cbjrrm/temp/eps_w   ! Bjerrum length in nm, necesaria para la ecuacion de Poisson 
    lb_m=cbjrrm/temp/eps_m   
    
    epsunt=4.d0*pi*cbjrrm/temp


   
    return
  end subroutine read_param
 




  !...................................................................
  !...................................................................
  !...................................................................
  
  ! Esta subrutina calcula el numero de ecuaciones.
  ! Recordar que para la solucion se resuelve Poisson y volumen, para la membrana solo Poisson
  ! y para la superficie area y condiciones de contorno. 
  
  subroutine get_neq(neq)
    implicit none
    
    integer,intent(out):: neq
    integer:: nc_pore,nc_mem,nc_tot,nc_sol,nc_surf
    

    nc_pore=ihm2*iRpore ! # volume cells inside pore, cantidad de celdas de volumen adentro del poro

    nc_mem=ihm2*dimr-nc_pore ! # volume cells inside membrane

    nc_tot=dimr*dimz ! # of cells including membrane+solution, numero total de celdas de volumen sin contar la superficie

    nc_sol=nc_tot-nc_mem ! # of solution cells

    nc_surf=(dimr-iRpore)+ihm2 ! # of surface cells   
    
    neq=nc_sol+nc_tot+2*nc_surf ! pack sol + poisson sol & mem + pack and electros surf   ! numero total de ecuaciones
    
    numeq=neq


  end subroutine get_neq






end module parameters
