module call_solver
  implicit none
  


contains
  
  
  
  subroutine llamador(x,xg,constr,neq_in,errtol,maxl_in,maxiter_in,ier,fnorm)
    !#######################################################################
    !     This subroutine calls kinsol
    !#######################################################################    
    implicit none
    
    integer,intent(in):: neq_in,maxl_in,maxiter_in

    real(8),intent(in):: errtol(2),xg(neq_in),constr(neq_in)
    
    integer,intent(out):: ier

    real(8),intent(out):: fnorm,x(neq_in)

    real(8):: fnormtol,scsteptol
    real(8):: scale(neq_in)
    real(8):: rout(2)        ! Kinsol additional output info
    
    integer :: i,globalstrat,maxl,maxlrst 
    
    integer*8:: maxiter,msbpre,iout(15),neq

    
    
    
    neq=neq_in
    
    

    
    !..........
    
    
    maxl=maxl_in              ! maximum Krylov subspace dimesion
    
    if (maxl_in==0) maxl=neq
    
    maxiter=maxiter_in            ! max # of iterations, default=200
    
    fnormtol=errtol(1)        ! function-norm stopping tolerance
    scsteptol=errtol(2)
    
    msbpre=10                 ! max # of iterations without prec setup    
    maxlrst=5                 ! maximum number of restarts
    globalstrat=0
    
    
    !..........


    do i=1,neq
       scale(i)=1d0          ! scaling vector 
       x(i)=xg(i)             ! initial guess
    enddo
      
      
    !..........


    call fnvinits(3,neq,ier)  ! inits NVECTOR module, 3 for kinsol
    if (ier/=0) then          ! ier error flag (0 is OK)
       write(*,'("SUNDIALS_ERROR: FNVINITS returned IER = ",i2)')ier
       return
    endif


    call fkinmalloc(iout,rout,ier) ! Allocates memory/output additional info
    if (ier/=0) then
       write(*,'("SUNDIALS_ERROR: FKINMALLOC returned IER = ",i2)')ier
       return
    endif
    
    
    call fkinsetiin('MAX_NITERS',maxiter,ier) ! Additional input info    
    call fkinsetiin('MAX_SETUPS',msbpre,ier) 
    call fkinsetrin('FNORM_TOL',fnormtol,ier)
    call fkinsetrin('SSTEP_TOL',scsteptol,ier)
    call fkinsetvin('CONSTR_VEC',constr,ier) ! constraint vector
    

    call fkinspgmr(maxl,maxlrst,ier)   
    if (ier/=0) then
       write(*,'("SUNDIALS_ERROR: FKINSPGMR returned IER = ",i2)')ier
       call fkinfree          ! free memory
       return
    endif
      

!    call fkinspilssetprec(1,ier) ! preconditions

    
    call fkinsol(x,globalstrat,scale,scale,ier) ! calls kinsol

    
    if (ier<0) then
       write(*,1240)ier,iout(9)
1240   format(/"SUNDIALS_ERROR: FKINSOL returned IER = ",i2,/,16x,&
            &"Linear Solver returned IER = ",i2)
       call fkinfree
       return
    endif


    write(*,'(/" FKINSOL return code is ",i3)')ier


    fnorm=rout(1)
    
    
    call fkinfree
    
    
    
  end subroutine llamador
  
  
 

 
end module call_solver


