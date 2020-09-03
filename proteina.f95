module proteina

  use parameters, only: vsol,sim_id,dimr,dimz,ihm2,iRpore
  implicit none
  
  real(8),parameter,private:: pi=acos(-1d0)
  
  integer,allocatable,private:: seed(:)

  integer:: ntyi,naac,ncfg

  integer,allocatable:: cnaa(:),typ(:),rhoq_i(:)

  
  real(8):: lsg_aa
  real(8),allocatable:: pKaa(:),zaa(:), vaa(:)
  real(8),allocatable,private:: rprot(:,:,:)
  integer,parameter:: ntht=36
  integer,parameter:: nr_int=60 ! must be even
  
  
  integer(kind=2), allocatable, dimension(:,:,:,:,:,:) :: naarz
  real(kind=4), allocatable, dimension(:,:,:,:,:,:) :: int_naarz
  

    
  real(8),allocatable,private:: tmp_rprot(:,:,:)


contains
  
  


  subroutine read_peptide(guess,deltaz,deltar)
    implicit none

    real(8),intent(in):: deltar,deltaz

    integer,intent(in):: guess

    integer:: i,icfg,iostat
    
    integer,parameter:: ncfgmx=105000 ! max # of (chain conformations)*rotations
    
    integer,parameter:: nrot=12   !numero de rotaciones al azar con angulos de Euler 

    integer,parameter:: nprotCG=83,nprot=63,nwrtpro=94

    character(5),allocatable:: aanm5(:)
    character(3),allocatable:: aanm(:)   

    logical:: writecfg

 




    if (.True.) then

       call with_solution(ncfgmx,deltar,deltaz)

       deallocate(rprot)

       return
    endif

    

    


    open(unit=nprotCG,file="protein_CG.mol")
    
    
    read(nprotCG,*)ntyi,lsg_aa

    
    allocate(pKaa(ntyi)) 
    allocate(zaa(ntyi))
    allocate(vaa(ntyi)) 
    allocate(aanm5(ntyi))
    allocate(aanm(ntyi))
    
    read(nprotCG,*)aanm5
    
    !aanm(0)="Neu"                ! AAC tipo cero es neutro !Usar 9 para neutro sin ads
    aanm(1:ntyi)=aanm5(:)(3:5)
    
    read(nprotCG,*)(pKaa(i),i=1,ntyi)
    read(nprotCG,*)(zaa(i),i=1,ntyi)
    read(nprotCG,*)(vaa(i),i=1,ntyi)

    vaa=vaa/vsol 
    
    close(nprotCG)
    
    
    allocate(cnaa(ntyi))
 
    
    open(unit=nprot,file="peptide.mol") ! this file contains peptide sequence
    
    read(nprot,*)naac
    
    allocate(typ(naac))
    
    read(nprot,*)(typ(i),i=1,naac)
    
    close(nprot)

    
    cnaa=0
    do i=1,naac
       if (typ(i)/=0) cnaa(typ(i))=cnaa(typ(i))+1
    enddo

    
    allocate(tmp_rprot(ncfgmx,naac,3)) 


    call creador_all(ncfg,nrot,naac,lsg_aa)



    allocate(rprot(ncfg,naac,3)) 

    rprot(:,:,:)=tmp_rprot(1:ncfg,:,:)

    deallocate(tmp_rprot)





    writecfg=.true.

    if (writecfg) then ! write xyz file with extra info

       open(unit=nwrtpro,file="peptide.xyz")

       write(nwrtpro,*)"Simulation ID:",sim_id

       write(nwrtpro,*)ntyi,lsg_aa
       write(nwrtpro,*)aanm5//"  "    
       write(nwrtpro,*)(pKaa(i),i=1,ntyi)
       write(nwrtpro,*)(zaa(i),i=1,ntyi)
       write(nwrtpro,*)(vaa(i)*vsol,i=1,ntyi)
       write(nwrtpro,*)

       write(nwrtpro,*)naac
       write(nwrtpro,*)(typ(i),i=1,naac)
       write(nwrtpro,*)

       do icfg=1,ncfg
          write(nwrtpro,*)naac ! xyz file
          write(nwrtpro,*)
          
          do i=1,naac            
             write(nwrtpro,*)aanm(typ(i)),rprot(icfg,i,1),rprot(icfg,i,2),rprot(icfg,i,3)
          enddo
       enddo
       
       close(nwrtpro)
       
    endif
    

    call discretize_rz(deltar,deltaz)


    deallocate(rprot)
    

    

    
  end subroutine read_peptide
  
  





  !####################################################################
  !....................................................................

  subroutine with_solution(ncfgmx,deltar,deltaz)
    implicit none

    integer,intent(in):: ncfgmx
    real(8),intent(in):: deltar,deltaz

    integer,parameter:: nrdpro=93

    integer:: i,iostat

    real(8),allocatable:: x(:),y(:),z(:)
    
    character(3):: auxcha3
    character(5),allocatable:: aanm5(:)
    character(3),allocatable:: aanm(:)   
    
    
    
    
    open(unit=nrdpro,file="protein_CG.mol")
    
    read(nrdpro,*)ntyi,lsg_aa
    
    allocate(aanm5(ntyi))
    allocate(aanm(ntyi))
    allocate(pKaa(ntyi)) 
    allocate(zaa(ntyi)) 
    allocate(vaa(ntyi))
    allocate(rhoq_i(ntyi))
    
    read(nrdpro,*)aanm5
    
    !aanm(0)="Neu"
    
    aanm(1:ntyi)=aanm5(:)(3:5)
    
    read(nrdpro,*)(pKaa(i),i=1,ntyi)
    read(nrdpro,*)(zaa(i),i=1,ntyi)
    read(nrdpro,*)(vaa(i),i=1,ntyi)
    read(nrdpro,*)(rhoq_i(i),i=1,ntyi)

    close(nrdpro)
       
    allocate(cnaa(ntyi)) 

    open(unit=99,file="input_1.xyz")

    read(99,*)naac
    
    allocate(typ(naac))
    allocate(x(naac),y(naac),z(naac))
    
    read(99,*)(typ(i),i=1,naac)
    read(99,*)
    
    vaa=vaa/vsol 
    
    cnaa=0
    do i=1,naac
       if (typ(i)/=0) cnaa(typ(i))=cnaa(typ(i))+1
    enddo
    
    
    ! now begin reading xyz file


    
    allocate(tmp_rprot(ncfgmx,naac,3)) 

    

    ncfg=0d0
    
    do 
       read(99,*,iostat=iostat)
       
       if (iostat<0) exit

       read(9,*,iostat=iostat)

       ncfg=ncfg+1
       
       do i=1,naac            
          read(99,*)auxcha3,x(i),y(i),z(i)

          tmp_rprot(ncfg,i,:)=(/x(i),y(i),z(i)/)
       enddo

    enddo

    close(99)
    


    
    allocate(rprot(ncfg,naac,3)) 
    
    rprot(:,:,:)=tmp_rprot(1:ncfg,:,:)
    
    deallocate(tmp_rprot)





    call discretize_rz(deltar,deltaz)

       
    
 




  end subroutine with_solution

     


  !####################################################################
  !....................................................................
  subroutine discretize_rz(deltar,deltaz)
    implicit none

    real(8),intent(in):: deltar,deltaz
    
    integer,parameter:: hmx=10 ! change if protein is larger

    
    
    integer:: icfg,ir,iz,hh,rr,r,itht,iaa,jr
    integer:: ir_max,r_min,r_max,h_min,h_max,iz_max,iz_min
    

    integer:: icfg_max, iaa_max


    

    real(8):: tr,tz,tx,ty,xxx,yyy,zzz,rrr,theta
    real(8):: vrr,vjr
    real(8):: zmax,zmin,rmax,rmin
    
    logical:: entra


    ir_max=dimr+hmx

    

    h_max=-1
    h_min=1
    
    r_max=-1
    r_min=1
    
    
    

    
    
    
    
    
    
    
    
    
    
!Calculate mins and maxs   

if (.false.) then
    ! Above Pore
    do iz=1,1

      do ir=1,ir_max
         
        do icfg=1,ncfg
            
          do itht=1,ntht 
            theta=2d0*pi*dble(itht-1)/dble(ntht)
            
              do jr=1,nr_int ! jr=nr_int/2 is redundant with previous loop...
                entra=.true.
                     
                tr=(ir-1+dble(2*jr-1)/dble(2*nr_int))*deltar ! posicion del centro de masa
                tz=(iz-0.5)*deltaz !Posicion CENTRO DE MASA en z

                tx=tr*cos(theta)
                ty=tr*sin(theta)


                ! Probar que la configuracion entre en el sistema
                do iaa=1,naac
                
                 xxx=rprot(icfg,iaa,1)+tx
                 yyy=rprot(icfg,iaa,2)+ty                 
                 zzz=rprot(icfg,iaa,3)+tz
                   
                 rrr=sqrt(xxx*xxx+yyy*yyy)


                 rr=int(rrr/deltar)+1-ir
                 hh=int(zzz/deltaz)+1-iz
                 if (zzz<0d0) hh=hh-1

    
                     
                enddo
              enddo
          enddo

        enddo

      enddo

    enddo
    
endif          

    
     
     
     
    
    zmax=maxval(abs(rprot(:,:,3)))
    rmax=maxval(sqrt(rprot(:,:,1)**2 + rprot(:,:,2)**2))
    
    
    
    !write(*,*) zmax, rmax
    
    r_max= int(rmax/deltar)+1
    r_min= -r_max
    h_max= nint(zmax/deltaz)
    h_min= -h_max
    
    
   

    
    ir_max=dimr-r_min ! asummes r_min<0

    iz_max=ihm2-h_min+1
    iz_min=1-h_max
    
    

     
    
    






    
    allocate(int_naarz(ncfg,1:ntyi,ir_max,r_min:r_max,iz_min:iz_max,h_min:h_max))
    int_naarz=0
    
    


    allocate(naarz(ncfg,1:ntyi,ir_max,r_min:r_max,iz_min:iz_max,h_min:h_max))
    naarz=0



!Fill naarz and int_naarz

   ! Inside Pore
    do iz=iz_min,ihm2

      do ir=1,iRpore
         
        do icfg=1,ncfg
            
          do itht=1,ntht
            entra=.true.
            theta=2d0*pi*dble(itht-1)/dble(ntht)

            tr=(ir-0.5)*deltar !Posicion CENTRO DE MASA en r
            tz=(iz-0.5)*deltaz !Posicion CENTRO DE MASA en z
                
            tx=tr*cos(theta)
            ty=tr*sin(theta)
            

            ! Probar que la configuracion entre en el sistema
            do iaa=1,naac
                
              xxx=rprot(icfg,iaa,1)+tx
              yyy=rprot(icfg,iaa,2)+ty                 
              zzz=rprot(icfg,iaa,3)+tz
                   
              rrr=sqrt(xxx*xxx+yyy*yyy)
              
              
              rr=int(rrr/deltar)+1
              hh=int(zzz/deltaz)+1
              if (zzz<0d0) hh=hh-1
                    
              if (hh<=ihm2 .and. rr>iRpore) then
                entra=.false.
              endif




            enddo

            

            if (entra) then
            
             do iaa=1,naac

               rr=int(rrr/deltar)+1-ir
               hh=int(zzz/deltaz)+1-iz
               if (zzz<0d0) hh=hh-1

                
               naarz(icfg,typ(iaa),ir,rr,iz,hh)=naarz(icfg,typ(iaa),ir,rr,iz,hh)+1
                  
            
            
              enddo
            endif
            
            
            
            
            

            do jr=1,nr_int ! jr=nr_int/2 is redundant with previous loop...
              entra=.true.
                     
              tr=(ir-1+dble(2*jr-1)/dble(2*nr_int))*deltar ! posicion del centro de masa
              tz=(iz-0.5)*deltaz !Posicion CENTRO DE MASA en z

              tx=tr*cos(theta)
              ty=tr*sin(theta)

              vjr=(2*(ir-1)*nr_int+2*jr-1)
              !vjr=2*jr-1


              ! Probar que la configuracion entre en el sistema
              do iaa=1,naac
                
                xxx=rprot(icfg,iaa,1)+tx
                yyy=rprot(icfg,iaa,2)+ty                 
                zzz=rprot(icfg,iaa,3)+tz

                   
                rrr=sqrt(xxx*xxx+yyy*yyy)
                
                
                
                rr=int(rrr/deltar)+1
                hh=int(zzz/deltaz)+1
                if (zzz<0d0) hh=hh-1
                    
                if (hh<=ihm2 .and. rr>iRpore) then
                  entra=.false.
                endif


              enddo
              
              if (entra) then
                do iaa=1,naac
                     
                  xxx=rprot(icfg,iaa,1)+tx
                  yyy=rprot(icfg,iaa,2)+ty                 
                  zzz=rprot(icfg,iaa,3)+tz
                     
                  rrr=sqrt(xxx*xxx+yyy*yyy)
                  
                  rr=int(rrr/deltar)+1-ir
                  hh=int(zzz/deltaz)+1-iz
                  if (zzz<0d0) hh=hh-1


                  vrr=2*(ir+rr)-1
                     
                  int_naarz(icfg,typ(iaa),ir,rr,iz,hh)=int_naarz(icfg,typ(iaa),ir,rr,iz,hh)+vjr/vrr

                     
                enddo
             
              endif
            
            enddo  
          enddo

        enddo

      enddo

    enddo


    ! Above Pore
    do iz=ihm2+1,iz_max

      do ir=1,ir_max
         
        do icfg=1,ncfg
            
          do itht=1,ntht 
            
            entra=.true.
            theta=2d0*pi*dble(itht-1)/dble(ntht)

            tr=(ir-0.5)*deltar !Posicion CENTRO DE MASA en r
            tz=(iz-0.5)*deltaz !Posicion CENTRO DE MASA en z
                
            tx=tr*cos(theta)
            ty=tr*sin(theta)
            
            ! Probar que la configuracion entre en el sistema
            do iaa=1,naac
                
              xxx=rprot(icfg,iaa,1)+tx
              yyy=rprot(icfg,iaa,2)+ty                 
              zzz=rprot(icfg,iaa,3)+tz
                   
              rrr=sqrt(xxx*xxx+yyy*yyy)
              
              rr=int(rrr/deltar)+1
              hh=int(zzz/deltaz)+1
              if (zzz<0d0) hh=hh-1
                    
              if (hh<=ihm2 .and. rr>iRpore) then
                entra=.false.
              endif




            enddo
            
            if (entra) then
                   
              do iaa=1,naac
                
                xxx=rprot(icfg,iaa,1)+tx
                yyy=rprot(icfg,iaa,2)+ty                 
                zzz=rprot(icfg,iaa,3)+tz
                   
                rrr=sqrt(xxx*xxx+yyy*yyy)

                rr=int(rrr/deltar)+1-ir
                hh=int(zzz/deltaz)+1-iz
  
                if (zzz<0d0) hh=hh-1

                
                naarz(icfg,typ(iaa),ir,rr,iz,hh)=naarz(icfg,typ(iaa),ir,rr,iz,hh)+1
                  
         
              enddo
            endif
            
            
            
            
            
           
            
            
                
            
            
            
            
            
              do jr=1,nr_int ! jr=nr_int/2 is redundant with previous loop...
                entra=.true.
                     
                tr=(ir-1+dble(2*jr-1)/dble(2*nr_int))*deltar ! posicion del centro de masa
                tz=(iz-0.5)*deltaz !Posicion CENTRO DE MASA en z

                tx=tr*cos(theta)
                ty=tr*sin(theta)

                vjr=(2*(ir-1)*nr_int+2*jr-1)
                !vjr=2*jr-1
              
                ! Probar que la configuracion entre en el sistema
                do iaa=1,naac
                
                 xxx=rprot(icfg,iaa,1)+tx
                 yyy=rprot(icfg,iaa,2)+ty                 
                 zzz=rprot(icfg,iaa,3)+tz
                   
                 rrr=sqrt(xxx*xxx+yyy*yyy)
                 
                 rr=int(rrr/deltar)+1
                 hh=int(zzz/deltaz)+1
                 if (zzz<0d0) hh=hh-1
                    
                 if (hh<=ihm2 .and. rr>iRpore) then
                   entra=.false.
                 endif




                enddo
              
            
                if (entra) then
                  
                  do iaa=1,naac
                     
                    xxx=rprot(icfg,iaa,1)+tx
                    yyy=rprot(icfg,iaa,2)+ty                 
                    zzz=rprot(icfg,iaa,3)+tz
                     
                    rrr=sqrt(xxx*xxx+yyy*yyy)
                  
                    rr=int(rrr/deltar)+1-ir
                    hh=int(zzz/deltaz)+1-iz
                    if (zzz<0d0) hh=hh-1

    
                    vrr=2*(ir+rr)-1
                       
                    int_naarz(icfg,typ(iaa),ir,rr,iz,hh)=int_naarz(icfg,typ(iaa),ir,rr,iz,hh)+vjr/vrr


                  enddo
                endif
              enddo
          enddo

        enddo

      enddo

    enddo
    
    


  
end subroutine discretize_rz


  






  
  !####################################################################
  !....................................................................

  subroutine creador_all(cuantas,nrot,nseg,lseg)
    implicit none
    
    integer,intent(in):: nseg,nrot
    integer,intent(inout):: cuantas

    real(8),intent(in):: lseg
  
    integer:: seed_size
    
    
    
    call random_seed() ! initialize with system generated seed
    call random_seed(size=seed_size) ! find out size of seed
    allocate(seed(seed_size))
    
    seed=1289
    
    call random_seed(put=seed) ! set current seed
    

    cuantas=0

    call sconf(3,nseg,lseg,nrot,cuantas)

    
  end subroutine creador_all
  
  



  !####################################################################
  !....................................................................
  
  subroutine cadena(nseg,xx,yy,zz,lseg,conf,selfav)
    implicit none
    
    integer,intent(in):: nseg
    real(8),intent(out):: xx(nseg),yy(nseg),zz(nseg)

    real(8),intent(in):: lseg

    character(1),intent(in):: conf(nseg)

    logical,intent(out):: selfav
    
    integer:: i,ii,j
    
    real(8):: sitheta,cotheta,siphip,cophip
    real(8):: m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3),I3(3,3)      
    real(8):: x,y,z,rn,dista,state
    
    
    
    sitheta=sin(68.0*pi/180.0)
    cotheta=cos(68.0*pi/180.0)
    siphip=sin(120.0*pi/180.0)
    cophip=cos(120.0*pi/180.0)
    
    
    tt(1,1)=cotheta
    tt(1,2)=sitheta
    tt(1,3)=0.0
    tt(2,1)=sitheta
    tt(2,2)=-cotheta
    tt(2,3)=0.0
    tt(3,1)=0.0
    tt(3,2)=0.0
    tt(3,3)=-1.0
    
    tp(1,1)=cotheta
    tp(1,2)=sitheta
    tp(1,3)=0.0
    tp(2,1)=sitheta*cophip
    tp(2,2)=-cotheta*cophip
    tp(2,3)=siphip
    tp(3,1)=sitheta*siphip
    tp(3,2)=-cotheta*siphip
    tp(3,3)=-cophip

    tm(1,1)=cotheta
    tm(1,2)=sitheta
    tm(1,3)=0.0
    tm(2,1)=sitheta*cophip
    tm(2,2)=-cotheta*cophip
    tm(2,3)=-siphip
    tm(3,1)=-sitheta*siphip
    tm(3,2)=cotheta*siphip
    tm(3,3)=-cophip
    
    I3=0d0
    forall(j=1:3) I3(j,j)=1d0

       
    !.....

    
    xx(1)=0d0
    yy(1)=0d0
    zz(1)=0d0

    m=I3
    
    do i=2,3 ! first monomers
       
       call mrrrr(m,tt,mm)
       m=mm
       
       z=m(1,1)*lseg
       x=m(2,1)*lseg
       y=m(3,1)*lseg
         
       zz(i)=zz(i-1)+z
       xx(i)=xx(i-1)+x
       yy(i)=yy(i-1)+y

    enddo

    do i=4,nseg
       
       select case(conf(i))
          
       case("t")          
          call mrrrr(m,tt,mm)
          
       case("+")
          call mrrrr(m,tp,mm)
          
       case("-")
          call mrrrr(m,tm,mm)
          
       end select

       m=mm
       
       z=m(1,1)*lseg
       x=m(2,1)*lseg
       y=m(3,1)*lseg
       
       zz(i)=zz(i-1)+z
       xx(i)=xx(i-1)+x
       yy(i)=yy(i-1)+y
       
    enddo
    
    
    selfav=.true.
    dista=0d0 ! check self-avoidance
     
    do i=4,nseg
       do j=1,i-3

          dista=sqrt((xx(j)-xx(i))**2+(yy(j)-yy(i))**2+(zz(j)-zz(i))**2)

          if (dista<lseg) selfav=.false. 

       enddo
    enddo
       

    
  end subroutine cadena





  !####################################################################
  !....................................................................
  
  subroutine mrrrr(a,b,c)
    implicit none
    
    integer:: i,j,k
    real(8),intent(in):: a(3,3),b(3,3)
    real(8),intent(out):: c(3,3)
    
    
    c=0d0
        
    do i=1,3
       do j=1,3
          do k=1,3
             c(i,j)=c(i,j)+a(i,k)*b(k,j)
          enddo
       enddo
    enddo
    
    
  end subroutine mrrrr
  
  
  
  
  
  !####################################################################
  !....................................................................
  
  subroutine rota(x,y,z,alpha,beta,gamma,xrot,yrot,zrot,nseg)
    implicit none
    
    real(8),intent(in):: x(nseg),y(nseg),z(nseg)
    real(8),intent(out):: xrot(nseg),yrot(nseg),zrot(nseg)
    
    real(8),intent(in):: alpha,beta,gamma
    
    integer,intent(in):: nseg
    integer:: i
    
    real(8):: a,b,c
    real(8):: sbe,cbe,sal,cal,sga,cga
    
    
    
    cbe=cos(beta)
    sbe=(1-cbe**2)**0.5
    
    cal=cos(alpha)
    sal=sin(alpha)
    
    cga=cos(gamma)
    sga=sin(gamma)
    
    do i=1,nseg
       
       a=z(i)
       b=x(i)
       c=y(i)
        
       zrot(i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
       xrot(i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
       yrot(i)=a*sbe*sga+b*sbe*cga+c*cbe
        
    enddo
     


  end subroutine rota
  
  
  
  
  
  
  !####################################################################
  !....................................................................
  
  subroutine rotaxy(x,y,theta,xrot,yrot,nseg)
    implicit none
    
    real(8),intent(in):: x(nseg),y(nseg)
    real(8),intent(out):: xrot(nseg),yrot(nseg)
    
    real(8),intent(in):: theta
    
    integer,intent(in):: nseg
    integer:: i
    
    real(8):: a,b
    real(8):: sithe,cothe
    
    
    sithe=sin(theta)
    cothe=cos(theta)
    
    
    do i=1,nseg
       
       a=x(i)
       b=y(i)

       xrot(i)=a*cothe-b*sithe
       yrot(i)=a*sithe+b*cothe
       
    enddo
    
    
    
  end subroutine rotaxy
  
  

  

  !####################################################################
  !....................................................................
  
  recursive subroutine sconf(iseg,nseg,lseg,nrot,cuantas)
    implicit none
     
    integer,intent(in):: iseg,nseg,nrot
    real(8),intent(in):: lseg

    integer,intent(inout):: cuantas
    
    integer:: i,j
    
    character(1),allocatable,save:: co(:)

    real(8):: alpha,beta,gamma,theta,rn    
    real(8):: xx(nseg),yy(nseg),zz(nseg)
    real(8):: xrot(nseg),yrot(nseg),zrot(nseg),xrtxy(nseg),yrtxy(nseg)
    
    logical:: selfav
    
    
    if (.not.allocated(co)) then
       allocate(co(nseg))
    endif
    
    
    if (iseg==nseg) then
       
       co(1:3)="0"
       call cadena(nseg,xx,yy,zz,lseg,co,selfav)

       if (.not.selfav) return
       
       cuantas=cuantas+1
       
       
       call align_z(xx,yy,zz,nseg)
          
       
       call save_conf(xx,yy,zz,nseg,cuantas)
       
       do i=1,nrot
          
          call random_number(rn)
          alpha=rn*2d0*pi
          
          call random_number(rn)
          beta=rn*pi
          
          call random_number(rn)
          gamma=rn*2d0*pi
          
          call rota(xx,yy,zz,alpha,beta,gamma,xrot,yrot,zrot,nseg)
          
          cuantas=cuantas+1
          
          call save_conf(xrot,yrot,zrot,nseg,cuantas)
          
       enddo
       
   
   
   
       
    else
       
       do i=1,3
          
          if (i==1) then
             co(iseg+1)="t"
          elseif (i==2) then
             co(iseg+1)="+"
          else 
             co(iseg+1)="-"
          endif
           

          call sconf(iseg+1,nseg,lseg,nrot,cuantas)
       enddo


    endif
 


    
  end subroutine sconf
  
  



  
  !####################################################################
  !....................................................................
  
  subroutine save_conf(x,y,z,nseg,icfg)
    implicit none
    
    real(8),intent(inout):: x(nseg),y(nseg),z(nseg)  !dimension numero de segmentos
    
    integer,intent(in):: nseg,icfg

    real(8):: xcm,ycm,zcm,rseg0(3)
    
    integer:: i,seg0,ncfgmx

    
    xcm=0.d0
    ycm=0.d0
    zcm=0.d0
    
    
    do i=1,nseg
       xcm=xcm+x(i)
       ycm=ycm+y(i)
       zcm=zcm+z(i)
    enddo
    
    xcm=xcm/dble(nseg)
    ycm=ycm/dble(nseg)
    zcm=zcm/dble(nseg)


    x=x-xcm
    y=y-ycm
    z=z-zcm


    do i=1,nseg
       tmp_rprot(icfg,i,1:3)=(/x(i),y(i),z(i)/)
    enddo
    
    
    ncfgmx=size(tmp_rprot,1)
    if (icfg>ncfgmx) then
       write(*,'(/1xa/)')"Error! Increase max # of chain conformations..."
       stop
    endif



  end subroutine save_conf
  
  
  !####################################################################
  !....................................................................
  
  subroutine align_z(x,y,z,nseg)
    implicit none
    integer:: i
    real(8),intent(inout):: x(nseg),y(nseg),z(nseg)
    integer,intent(in):: nseg
    real(8):: u,v,w, r31
    real(8):: tz(3,3)
    real(8):: xnew(nseg),ynew(nseg),znew(nseg)
    
    
    
    u=x(3)-x(1)
    v=y(3)-y(1)
    w=z(3)-z(1)
    
    r31=sqrt(u**2+v**2+w**2)
    
    tz(1,1)=w/r31
    tz(1,2)=0d0
    tz(1,3)=-sqrt(u**2+v**2)/r31
    tz(2,1)=0
    tz(2,2)=1
    tz(2,3)=0
    tz(3,1)=sqrt(u**2+v**2)/r31
    tz(3,2)=0
    tz(3,3)=w/r31
    
    do i=1,nseg
    
      xnew(i)=tz(1,1)*x(i)+tz(1,2)*y(i)+tz(1,3)*z(i)
      ynew(i)=tz(2,1)*x(i)+tz(2,2)*y(i)+tz(2,3)*z(i)
      znew(i)=tz(3,1)*x(i)+tz(3,2)*y(i)+tz(3,3)*z(i)
      
    enddo
    
    x=xnew
    y=ynew
    z=znew
    
    
    
  end subroutine align_z
    
  
  
  

  
end module proteina



























