module printout


contains



  
  subroutine print_screen(sim_id,neq,dt)
    implicit none
    !#######################################################################
    !     Prints inputs of the calculation
    !#######################################################################
    
    character(14),intent(in):: sim_id

    integer,intent(in):: neq
    real :: dt

    
    
    write(*,'(/"*****"43x"*****"/&
         &"***** Membrane with Pore - Peptide/Protein Sol. *****"&
         &   /"*****"43x"*****")')


    
    
    
    write(*,'(/1x,"You should place stuff you want to print to screen here!!!!")') 
    
    write(*,'(/1x,"# of equations:",i8)')neq
    
    write(*,*)"Setup took (seconds):",dt
    
    write(*,'(/1x,"Simulation ID: "a14)')sim_id



  
    
    
    
  end subroutine print_screen
  
  
  
  
  
  
  !.......................................................................
  
  subroutine print_xs(ier,fnorm,goodx,neq,x,constr,dt)
    
    !#######################################################################
    !     Prints solutions
    !#######################################################################
    use parameters
    implicit none
    
    integer,intent(in):: neq,ier
    
    real(8),intent(in):: fnorm,x(neq),constr(neq)
    real,intent(in):: dt
    
    logical,intent(out):: goodx
  
    character(6):: spH
    character(16):: outname

    integer:: i
    
    logical:: satcon,fnorm_is_finite,fnorm_is_good

    integer,parameter:: outfile=15



    ! first let's check this is a solution
  
    fnorm_is_finite=.true.
    if (fnorm/=fnorm) fnorm_is_finite=.false.
    fnorm_is_finite=(abs(fnorm)<huge(fnorm))
  
    fnorm_is_good=.true.
    if (fnorm_is_finite) then
       fnorm_is_good=(fnorm>0d0)
    else
       fnorm_is_good=.false.
    endif
    
    
    satcon=.true.    
    call satis_const(satcon,neq,x,constr)
    
    goodx=(fnorm_is_good.and.satcon.and.ier==0)
    
  
    !..........
  
  
    write(spH,"(f6.3)")pH
    if (pH<10d0.and.spH/="10.000") then
       write(spH,"(f5.3)")pH
       spH="0"//spH
    endif

  
    if (goodx) then
       outname="xs-"//spH//".dat"
    else
       outname="no_xs-"//spH//".bad"
    endif


    open(outfile,file=outname)


    !..........

  
    write(*,*)
    write(*,*)"FNORM=",fnorm
    write(*,*)
    write(*,*)"Time elapsed (seconds):",dt
  
  
    write(outfile,*)"Simulation_ID: ",sim_id
    write(outfile,*)ier,fnorm,dt
    write(outfile,*)csalt,pH
    write(outfile,*)dimr,deltar,dimz,deltaz
    write(outfile,*)iRpore,R_pore
    write(outfile,*)ihm2,h_mem
    write(outfile,*)xb_I,pK_I,z_I,a_I*a_l,a_N*a_l,a_l
    write(outfile,*)eps_m,lb_m,eps_w,lb_w
    write(outfile,*)temp,pKw
    write(outfile,*)vsol,vH*vsol,vOH*vsol,vpls*vsol,vmin*vsol
    write(outfile,*)zH,zOH,zpls,zmin
    write(outfile,*)rhopbulk/Na
    write(outfile,*)


    write(outfile,*)neq
    do i=1,neq
       write(outfile,*)x(i)
    enddo


    close(outfile)

  
    if (.not.goodx) write(*,'(/1x"##### ERROR: not a good solution!")')
  
  
    

  










  
  end subroutine print_xs






  !.............................................................................
  
  subroutine satis_const(satcon,neq,x,constr)
    
    !***************************************************************************
    !***************************************************************************
    
    implicit none
    logical,intent(out):: satcon
    
    integer,intent(in):: neq
    real(8),intent(in):: x(neq),constr(neq)
    
    integer:: i
    

    satcon=.true.

    do i=1,neq
       
       if (constr(i)==1d0) then
          satcon=(x(i)>=0d0)
          
       elseif(constr(i)==-1d0) then
          satcon=(x(i)<=0d0)
          
       elseif(constr(i)==2d0) then
          satcon=(x(i)>0d0)
          
       elseif (constr(i)==-2d0) then
          satcon=(x(i)<0d0)
          
       endif
       
       if (.not.satcon) return
       
    enddo
    
    
  end subroutine satis_const
  
  
  
  
  
  
  
end module printout
