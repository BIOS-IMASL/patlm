subroutine fkfun(x,f,ier)   !fkfunk nombre generico que usa el solver, x= solucion, f=0, ier =  mens. de salida del solver. ver manual
  
  !#######################################################################
  !     Calculates functions f to be solved (f(x)=0)
  !#######################################################################
  use parameters
  use proteina, only: zaa,pKaa,ntyi,naac,vaa,ncfg,cnaa,typ,lsg_aa,naarz,int_naarz, ntht, nr_int,rhoq_i

  implicit none

  integer:: ier,i,j,ii,iiaux
  integer:: iz,jz,kz,ir,jr,kr,hh,h,h_min,h_max,rr,r,r_min,r_max
  integer:: neq,nc_pore,nc_mem,nc_tot,nc_sol,nc_surf

  real(8):: f(numeq),x(numeq)
  real(8):: epsi,bpsi,xsol,ysol,xH,xOH,xmin,xpls
  real(8):: lnprz,prz,eDmup,fb_I,faux
  real(8):: dpsiz_s,dpsiz_m,dpsir_s,dpsir_m
  real(8):: d2psiz,d2psir,d2psi,constq
  real(8):: bchi,total_volume
  integer:: vrr

  integer:: ir_max,iz_min,iz_max, iz_mxx, ir_mxx

  real(8),allocatable,dimension(:,:),save:: psi,phiw,xtot,rhoq
  real(8),allocatable,dimension(:),save:: fbaa,faa_tmp
  real(8),allocatable,dimension(:),save:: psis1,xNs1,xIs1,f_Is1
  real(8),allocatable,dimension(:),save:: psis2,xNs2,xIs2,f_Is2
  real(8),allocatable,dimension(:,:,:),save:: avnaa,faa

  logical:: loopout

     
  !..........
  

   
  if (.not.allocated(psi)) allocate(psi(dimr,dimz))
  if (.not.allocated(phiw)) allocate(phiw(dimr,dimz))
  if (.not.allocated(xtot)) allocate(xtot(dimr,dimz))
  if (.not.allocated(rhoq)) allocate(rhoq(dimr,dimz))

  if (.not.allocated(fbaa)) allocate(fbaa(ntyi))
  if (.not.allocated(faa_tmp)) allocate(faa_tmp(ntyi))
  if (.not.allocated(faa)) allocate(faa(ntyi,dimr,dimz))
  if (.not.allocated(avnaa)) allocate(avnaa(1:ntyi,dimr,dimz))

  if (.not.allocated(psis1)) allocate(psis1(iRpore+1:dimr))
  if (.not.allocated(f_Is1)) allocate(f_Is1(iRpore+1:dimr)) !superficie de arriba
  if (.not.allocated(xIs1)) allocate(xIs1(iRpore+1:dimr))
  if (.not.allocated(xNs1)) allocate(xNs1(iRpore+1:dimr))

  if (.not.allocated(psis2)) allocate(psis2(ihm2))
  if (.not.allocated(f_Is2)) allocate(f_Is2(ihm2))
  if (.not.allocated(xIs2)) allocate(xIs2(ihm2))
  if (.not.allocated(xNs2)) allocate(xNs2(ihm2))


  psi=0d0
  phiw=0d0
  xtot=0d0
  rhoq=0d0

  fbaa=0d0
  faa_tmp=0d0

  faa=0d0
  avnaa=0d0

  psis1=0d0
  f_Is1=0d0
  xIs1=0d0
  xNs1=0d0

  psis2=0d0
  f_Is2=0d0
  xIs2=0d0
  xNs2=0d0

  bchi=3d0
  total_volume=0d0
      
  
  !.......... Ionizable lipid


  fb_I=1d0/(1d0+10d0**(pK_I-pH)) ! CAREFUL, THIS IS IN A DILUTE SOLUTION. Far from pore fbI depends on psis1(dimr)


  !.......... Protein
  
  do i=1,ntyi ! loop sobre tipos, numero de tipos de aac, excluyendo los descargados
     fbaa(i)=1d0/(1d0+10.d0**(sign(1d0,zaa(i))*(pH-pKaa(i))))
  enddo


  

  do j=1,ntyi
  	 total_volume=total_volume+dble(cnaa(j))*vaa(j)
  enddo

  eDmup=rhopbulk/xsolbulk**total_volume/dble(ncfg)

  do j=1,ntyi
     eDmup=eDmup*(1d0-fbaa(j))**cnaa(j)  
  enddo





!!!.............................

  
  nc_pore=ihm2*iRpore ! # volume cells inside pore
  
  nc_mem=ihm2*dimr-nc_pore ! # volume cells inside membrane
  
  nc_tot=dimr*dimz ! # of cells including membrane+solution
  
  nc_sol=nc_tot-nc_mem ! # of solution cells
  
  nc_surf=(dimr-iRpore)+ihm2 ! # of surface cells





  !.......... Assign xs ...........
    

  ! el vector x contiene las incognitas, al comienzo x contiene lo que viene del guess
  ! Inside Pore
  
  do iz=1,ihm2              ! Inside pore. Asigna soluciones previas a cantidades fisicas
     do ir=1,iRpore

        ii=ir+iRpore*(iz-1) ! busco indice de la solucion

        phiw(ir,iz)=x(ii)  ! volume fraction del H2O en ese punto.
         
        psi(ir,iz)=x(ii+nc_sol)    ! Potencial electrostatico en unidades de kT*electro-charge

        epsi=dexp(-psi(ir,iz))     ! exp de pot electrostatico en unidades de kT*electro-charge

        faa_tmp=fbaa/(1d0-fbaa)*epsi**zaa     ! auxiliar para grado de carga de aac en el punto f_tau(r,z) pag. 10 ! OJO operacion vetorial
        
        faa(:,ir,iz)=faa_tmp(:)/(1d0+faa_tmp(:))     ! Despeja f_tau    pag. 10
     enddo
  enddo
  
  
  iiaux=iRpore*ihm2
  
  do iz=ihm2+1,dimz ! Above pore en todas las direcciones = ariba de la membrana
     do ir=1,dimr
        
        jz=iz-ihm2
        ii=ir+dimr*(jz-1)+iiaux     ! busco indice de la solucion

        phiw(ir,iz)=x(ii)
        
        psi(ir,iz)=x(ii+nc_sol)
        
        epsi=dexp(-psi(ir,iz))
        
        faa_tmp=fbaa/(1d0-fbaa)*epsi**zaa
        
        faa(:,ir,iz)=faa_tmp(:)/(1d0+faa_tmp(:))

     enddo
  enddo
  


  
  do iz=1,ihm2              ! Inside membrane
     do ir=iRpore+1,dimr

        jr=ir-iRpore
        ii=jr+(dimr-iRpore)*(iz-1)+2*nc_sol    ! busco indice de la solucion

        psi(ir,iz)=x(ii)

     enddo
  enddo
  

  !celdas de superficie

  ii=nc_sol+nc_tot

  do ir=iRpore+1,dimr
     
     ii=ii+1

     xNs1(ir)=x(ii)            ! Fraccion de area que ocupa el lipido neutro  - equiva al volumen del agua en bulk 
     
     psis1(ir)=x(ii+nc_surf)   ! superficie de arriba. Potencial en esa superficie

     
     f_Is1(ir)=(fb_I*dexp(-psis1(ir)*z_I))/(1-fb_I+fb_I*dexp(-psis1(ir)*z_I))
     
   enddo
  

   do ir=iRpore+1,dimr

      
      faux=(1d0-fb_I+fb_I*dexp(-psis1(ir)*z_I))/(1d0-fb_I+fb_I*dexp(-psis1(dimr)*z_I))

      xIs1(ir)=xb_I*(xNs1(ir)/(1d0-xb_I))**(a_I/a_N)*faux

   enddo


  do iz=1,ihm2
     
     ii=ii+1
     
     xNs2(iz)=x(ii)    ! fraccion de area del lipido neutro en la superficie 2 , pared del cilindro
     
     psis2(iz)=x(ii+nc_surf)

     
     f_Is2(iz)=(fb_I*dexp(-psis2(iz)*z_I))/(1-fb_I+fb_I*dexp(-psis2(iz)*z_I))

     faux=(1d0-fb_I+fb_I*dexp(-psis2(iz)*z_I))/(1d0-fb_I+fb_I*dexp(-psis1(dimr)*z_I))

     xIs2(iz)=xb_I*(xNs2(iz)/(1d0-xb_I))**(a_I/a_N)*faux
     
  enddo







  r_min=lbound(naarz,4)
  r_max=ubound(naarz,4)

  h_min=lbound(naarz,6)
  h_max=ubound(naarz,6)

  iz_min=-h_max+1
  iz_max=dimz-h_min ! asummes h_min<0

  ir_max=dimr-r_min ! asummes r_min<0 
  
  
  iz_mxx=ubound(naarz,5)
  ir_mxx=ubound(naarz,3)




  avnaa=0d0

           


  !.......... calculate monomers vol fract ...........
  

  
  ! Inside Pore
  do iz=iz_min,ihm2 ! Inside pore
     do ir=1,iRpore ! should modify if pore is thin and peptide long
        
        if (iz>=1) then           
           epsi=dexp(-psi(ir,iz))   
           bpsi=psi(ir,iz)
           xsol=phiw(ir,iz)
        else
           epsi=dexp(-psi(ir,-iz+1))
           bpsi=psi(ir,-iz+1)
           xsol=phiw(ir,-iz+1)           
        endif
        
        ysol=xsol/xsolbulk
        
        xH=xHbulk*ysol**vH*epsi**zH          ! volume fraction H3O+ locales
        xOH=xOHbulk*ysol**vOH*epsi**zOH      ! volume fraction OH-   locales
        xmin=xminbulk*ysol**vmin*epsi**zmin  ! volume fraction anion  
        xpls=xplsbulk*ysol**vpls*epsi**zpls

  
        do i=1,ncfg ! protein 
           
           lnprz=0d0
           prz=0d0

           loopout=.false.

           do h=h_min,h_max
              do r=r_min,r_max ! Loop sobre x 
              
                 rr=ir+r 
                 if (rr<1) cycle
                 
                 hh=iz+h  
                 if (hh<1) then
                    hh=-hh+1  !uso de simetria en z
                 endif

                 if (hh<=ihm2.and.rr>iRpore) cycle              
                    


                 if (rr==iRpore .and. hh<ihm2+1) then
                 	lnprz=lnprz+naarz(i,8,ir,r,iz,h)*bchi/dble(ntht)
                 else if (rr>iRpore .and. hh==ihm2+1) then
                 	lnprz=lnprz+naarz(i,8,ir,r,iz,h)*bchi/dble(ntht)
                 endif
                 
                 ! sumo la contribucion de las cargas
                 ! ahora todos son j 
                 do j=1,ntyi
                 	lnprz=lnprz+naarz(i,j,ir,r,iz,h)*dlog(phiw(rr,hh))*vaa(j)/dble(ntht)
                    lnprz=lnprz-naarz(i,j,ir,r,iz,h)*dlog(1d0-faa(j,rr,hh))/dble(ntht)
                 enddo
                 
              enddo
              
           enddo

           

              prz=dexp(lnprz)

           
           do h=h_min,h_max
              do r=r_min,r_max ! Loop sobre x 
                 
                 hh=iz+h  
                 rr=ir+r
                 
                 if (rr<1) cycle
                 
                 
                 !vjr_aux=2*(ir-1)*nr_int
                 vrr=2*(ir+r)-1
                 if (hh>=1) then
                    avnaa(:,rr,hh)=avnaa(:,rr,hh)+prz*int_naarz(i,:,ir,r,iz,h)/dble(ntht*nr_int**2)
                 endif
                       
              enddo
           enddo
              
        enddo   ! Fin Loop sobre configuraciones
        
        if (iz<1) cycle
                   
        xtot(ir,iz)=xH+xOH+xpls+xmin+xsol      ! volume fraction de todo menos el peptido
        
        ! Densidad de carga de todo menos la contribucion del peptido
        rhoq(ir,iz)=xH*zH/vH+xOH*zOH/vOH+xpls*zpls/vpls+xmin*zmin/vmin 
        
     enddo ! Fin Loop sobre coordenada global (cm)
  enddo    ! Fin Loop sobre coordenada global (cm)
  
  


  do iz=ihm2+1,iz_max  ! Above pore
     do ir=1,ir_max

        if (iz>dimz) then           
           if (ir>dimr) then           
              epsi=dexp(-psi(dimr,dimz))
              bpsi=psi(dimr,dimz)
              xsol=phiw(dimr,dimz)            
           else    
              epsi=dexp(-psi(ir,dimz))
              bpsi=psi(ir,dimz)
              xsol=phiw(ir,dimz)     
           end if
        else
           if(ir>dimr)then
              epsi=dexp(-psi(dimr,iz))
              bpsi=psi(dimr,iz)
              xsol=phiw(dimr,iz)
           else
              epsi=dexp(-psi(ir,iz))
              bpsi=psi(ir,iz)
              xsol=phiw(ir,iz)    
           endif
        endif
           
        ysol=xsol/xsolbulk
        
        xH=xHbulk*ysol**vH*epsi**zH
        xOH=xOHbulk*ysol**vOH*epsi**zOH
        xmin=xminbulk*ysol**vmin*epsi**zmin
        xpls=xplsbulk*ysol**vpls*epsi**zpls
           
        do i=1,ncfg            ! protein 
              
           lnprz=0d0
           prz=0d0

           !loopout=.false.
           
           do h=h_min,h_max
              do r=r_min,r_max
              
                 rr=ir+r
                 if (rr<1) cycle
                 if (rr>dimr) rr=dimr

                 hh=iz+h
                 if (hh<1) hh=-hh+1
                 if (hh>dimz) hh=dimz
                 
                 if (hh<=ihm2.and.rr>iRpore) cycle
                    

                 jz=iz
                 jr=ir
                 
                 if (iz>iz_mxx) jz=iz_mxx
                 if (ir>ir_mxx) jr=ir_mxx
                    
                 

                 if (rr==iRpore .and. hh<ihm2+1) then
                 	lnprz=lnprz+naarz(i,8,jr,r,jz,h)*bchi/dble(ntht)
                 else if (rr>iRpore .and. hh==ihm2+1) then
                 	lnprz=lnprz+naarz(i,8,jr,r,jz,h)*bchi/dble(ntht)
                 endif

                 do j=1,ntyi
                 	lnprz=lnprz+naarz(i,j,jr,r,jz,h)*dlog(phiw(rr,hh))*vaa(j)/dble(ntht)
                    lnprz=lnprz-naarz(i,j,jr,r,jz,h)*dlog(1d0-faa(j,rr,hh))/dble(ntht)
                 enddo
                 
              enddo
              !if (loopout) exit
           enddo
           
           
           !if (loopout) then
            !  cycle
           !else
              prz=dexp(lnprz)
          ! endif
           

           do h=h_min,h_max
              do r=r_min,r_max ! Loop sobre x
                 
                 hh=iz+h  
                 rr=ir+r
                 if (rr<1) cycle
                 
                 
                 vrr=2*(ir+r)-1
                 if (hh>=1.and.hh<=dimz.and.rr<=dimr) then
                   avnaa(:,rr,hh)=avnaa(:,rr,hh)+prz*int_naarz(i,:,jr,r,jz,h)/dble(ntht*nr_int**2)
                 endif
                    
              enddo
           enddo

              
        enddo
        
        if (iz>dimz.or.ir>dimr) cycle
                   
        xtot(ir,iz)=xH+xOH+xpls+xmin+xsol
        
        rhoq(ir,iz)=xH*zH/vH+xOH*zOH/vOH+xpls*zpls/vpls+xmin*zmin/vmin 
        
     enddo
  enddo
  











 


  do j=1,ntyi
     avnaa(j,:,:)=avnaa(j,:,:)*vaa(j)*eDmup*vsol
  enddo
  
  



  
  
  if (.false.) then

     write(*,'(/1x,"ihm2=",i3,1x"-- iRpore=",i3)')ihm2,iRpore
     do ir=0,dimr-10,10
        write(*,'(1x"iz\ir",10(3x,i3,5x))')(i,i=ir+1,ir+10)
        do iz=1,20

        enddo
        write(*,*)
     enddo


     write(*,'(/1x"just printing initial volume fractions. Exiting."/)')
     stop

  endif
  


























  !..........................................
  !..... Now construct f(x)'s (eqns to solve)
  !..........................................
  

  constq=4.0*pi*lb_w/vsol
  
   
  f=0d0
  
  
  do iz=1,ihm2 ! Inside Pore
     do ir=1,iRpore

        kz=iz-1
        if (iz==1) kz=1
        
        ii=ir+iRpore*(iz-1)

        do j=1,ntyi
        	xtot(ir,iz)=xtot(ir,iz)+avnaa(j,ir,iz)
        end do
        
        f(ii)=xtot(ir,iz)-1d0


        d2psiz=psi(ir,iz+1)+psi(ir,kz)-2d0*psi(ir,iz)               
        d2psiz=d2psiz/deltaz**2
         

        if (ir==1) then

           d2psir=psi(2,iz)/(dble(ir)-0.5)/2d0        
           d2psir=d2psir+psi(2,iz)-2d0*psi(1,iz)


        elseif (ir==iRpore) then


           d2psir=(psi(ir,iz)-psi(ir-1,iz))/(dble(ir)-0.5)

           d2psir=d2psir+psis2(iz)+psi(ir-1,iz)-2d0*psi(ir,iz)


        else


           d2psir=(psi(ir,iz)-psi(ir-1,iz))/(dble(ir)-0.5)

           d2psir=d2psir+psi(ir+1,iz)+psi(ir-1,iz)-2d0*psi(ir,iz)
        
         
        endif


        d2psir=d2psir/deltar**2
        
        d2psi=d2psir+d2psiz
        
        do j=1,ntyi           
           rhoq(ir,iz)=rhoq(ir,iz)+faa(j,ir,iz)*zaa(j)*avnaa(rhoq_i(j),ir,iz)/vaa(rhoq_i(j)) ! vaa(:) DONE
        enddo
               
        f(ii+nc_sol)=d2psi+rhoq(ir,iz)*constq ! Poisson eqn
        
     enddo
  enddo
  
  
  


  iiaux=iRpore*ihm2
  
  do iz=ihm2+1,dimz ! Above pore
     do ir=1,dimr
        
        jz=iz-ihm2
        ii=ir+dimr*(jz-1)+iiaux
        
        do j=1,ntyi
        	xtot(ir,iz)=xtot(ir,iz)+avnaa(j,ir,iz)
        end do
        
        f(ii)=xtot(ir,iz)-1d0

        
        kz=iz+1
        if (iz==dimz) kz=dimz

        kr=ir+1
        if (ir==dimr) kr=dimr
        

        if (iz==ihm2+1.and.ir>iRpore) then

           d2psiz=psi(ir,iz+1)+psis1(ir)-2d0*psi(ir,iz)
           
        else

           d2psiz=psi(ir,kz)+psi(ir,iz-1)-2d0*psi(ir,iz)
               
        endif
 
        d2psiz=d2psiz/deltaz**2
        
        if (dimr>1) then

           if (ir==1) then
 
              d2psir=psi(2,iz)/(dble(ir)-0.5)/2d0        
              d2psir=d2psir+psi(2,iz)-2d0*psi(1,iz)

           
           else


              d2psir=(psi(ir,iz)-psi(ir-1,iz))/(dble(ir)-0.5)

              d2psir=d2psir+psi(kr,iz)+psi(ir-1,iz)-2d0*psi(ir,iz)


           endif

        else
           d2psir=0d0
        endif
        
        d2psir=d2psir/deltar**2
        
        d2psi=d2psir+d2psiz
        
        do j=1,ntyi
           rhoq(ir,iz)=rhoq(ir,iz)+faa(j,ir,iz)*zaa(j)*avnaa(rhoq_i(j),ir,iz)/vaa(rhoq_i(j))   !vaa(:) DONE
        enddo
        
        f(ii+nc_sol)=d2psi+rhoq(ir,iz)*constq ! Poisson eqn
        
     enddo
  enddo
  



  
  do iz=1,ihm2              ! Inside membrane
     do ir=iRpore+1,dimr
        
        jr=ir-iRpore
        ii=jr+(dimr-iRpore)*(iz-1)+2*nc_sol

        kr=ir+1
        if (ir==dimr) kr=dimr
        

        if (iz==1) then

           d2psiz=psi(ir,iz+1)-psi(ir,iz)


        elseif (iz==ihm2) then           
           
           d2psiz=psis1(ir)+psi(ir,iz-1)-2d0*psi(ir,iz)

                   
        else

           d2psiz=psi(ir,iz+1)+psi(ir,iz-1)-2d0*psi(ir,iz)        


        endif
        
        d2psiz=d2psiz/deltaz**2

                    


        if (ir>1) then

           d2psir=(psi(ir,iz)-psi(ir-1,iz))/(dble(ir)-0.5)
        
           d2psir=d2psir+psi(kr,iz)+psi(ir-1,iz)-2d0*psi(ir,iz)


        endif
        
        d2psir=d2psir/deltar**2
    


    
        d2psi=d2psir+d2psiz

        
        f(ii)=d2psi ! Laplace eqn







        
     enddo
  enddo
  
  

  


  
  ii=nc_sol+nc_tot  

  do ir=iRpore+1,dimr ! Planar surface
     
     ii=ii+1
     
     dpsiz_s=psi(ir,ihm2+1)-psis1(ir)


     dpsiz_m=psis1(ir)-psi(ir,ihm2)



     f(ii+nc_surf)=(-dpsiz_s/lb_w+dpsiz_m/lb_m)/deltaz&
          &-4d0*pi*xIs1(ir)*f_Is1(ir)*z_I/a_i/a_L


     f(ii)=xIs1(ir)+xNs1(ir)-1d0

  enddo
  

  





  ier=0
  
  
end subroutine fkfun
