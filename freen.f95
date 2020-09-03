 subroutine freen(x,neq)
  
  !#######################################################################
  !     Calculates stuff
  !#######################################################################
  use parameters
  use proteina

  implicit none

  integer,intent(in):: neq
  real(8),intent(in):: x(neq)

  integer:: ier

  integer:: i,j,ii,iiaux,iaa
  integer:: iz,jz,kz,ir,jr,hh,h,h_min,h_max,rr,r,r_min,r_max
  integer:: rxmax,rxmin,rymin,rymax,hmin,hmax,ry
  integer:: nc_pore,nc_mem,nc_tot,nc_sol,nc_surf

  real(8):: epsi,bpsi,xsol,ysol,xH,xOH,xmin,xpls
  real(8):: lnprz,prz,eDmup,fb_I,eDmui
  real(8):: dpsiz_s,dpsiz_m,dpsir_s,dpsir_m
  real(8):: constq
  real(8):: dpsir2, dpsiz2
  real(8):: xx,yy,rrr,faux,volir,voltot
  real(8):: ph0,q,qpore,qpore2,qlip,qlip2,sumavnaa,sumavnaabulk,gamma
  real(8):: Stf,Stp,bFcp,bFms,Ue,Ueb,FreeE,bFmup,bFmup2,bFmul2,bFchm,Ues
  real(8):: Stfb,Stpb,bFcpb,bFmsb,bFmupb,bFmup2b,bFmul2b,bFbchm,Uesb
  real(8):: bmuH,bmuOH,bmupls,bmumin
  real(8):: bchi, total_volume

  real(8):: q_sol,q_solRef,q_lip,q_lipRef

  

  integer:: vrr

  integer:: ir_max,iz_min,iz_max, iz_mxx, ir_mxx

  real(8),allocatable,dimension(:,:),save:: psi,phiw,xtot,rhoq
  real(8),allocatable,dimension(:),save:: fbaa,faa_tmp,confprz
  real(8),allocatable,dimension(:),save:: psis1,xNs1,xIs1,f_Is1
  real(8),allocatable,dimension(:),save:: psis2,xNs2,xIs2,f_Is2
  real(8),allocatable,dimension(:,:,:),save:: avnaa,faa

  real(8),allocatable,dimension(:,:),save:: lnrhop,rhop,umf,rhomin,rhopls,rhoH,rhoOH


  logical:: loopout

  character(6):: spH
  character(2):: sj
     
  !..........

  if (.not.allocated(psi)) allocate(psi(dimr,dimz))
  if (.not.allocated(phiw)) allocate(phiw(dimr,dimz))
  if (.not.allocated(xtot)) allocate(xtot(dimr,dimz))
  if (.not.allocated(rhoq)) allocate(rhoq(dimr,dimz))

  if (.not.allocated(fbaa)) allocate(fbaa(ntyi))
  if (.not.allocated(confprz)) allocate(confprz(ncfg))
  if (.not.allocated(faa_tmp)) allocate(faa_tmp(ntyi))
  if (.not.allocated(faa)) allocate(faa(ntyi,dimr,dimz))
  if (.not.allocated(avnaa)) allocate(avnaa(1:ntyi,dimr,dimz))

  if (.not.allocated(lnrhop)) allocate(lnrhop(dimr,dimz))
  if (.not.allocated(rhop)) allocate(rhop(dimr,dimz))
  if (.not.allocated(rhomin)) allocate(rhomin(dimr,dimz))
  if (.not.allocated(rhopls)) allocate(rhopls(dimr,dimz))
  if (.not.allocated(rhoH)) allocate(rhoH(dimr,dimz))
  if (.not.allocated(rhoOH)) allocate(rhoOH(dimr,dimz))



  if (.not.allocated(psis1)) allocate(psis1(iRpore+1:dimr))
  if (.not.allocated(f_Is1)) allocate(f_Is1(iRpore+1:dimr)) !superficie de arriba
  if (.not.allocated(xIs1)) allocate(xIs1(iRpore+1:dimr))
  if (.not.allocated(xNs1)) allocate(xNs1(iRpore+1:dimr))

  if (.not.allocated(psis2)) allocate(psis2(ihm2))
  if (.not.allocated(f_Is2)) allocate(f_Is2(ihm2))
  if (.not.allocated(xIs2)) allocate(xIs2(ihm2))
  if (.not.allocated(xNs2)) allocate(xNs2(ihm2))
  
  if (.not.allocated(umf)) allocate(umf(dimr,dimz))


  write(*,*)
  write(*,*)"...Process results here (replicate parts of fkfun subroutine)"

  !..........
   


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
  f_Is2=0
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
     eDmup=eDmup*(1d0-fbaa(j))**cnaa(j)  ! esto como estaba antes
  enddo










  bmuH=-dlog(xHbulk/vH/xsolbulk**vH) ! bmu^0_H
  bmuOH=-dlog(xOHbulk/vOH/xsolbulk**vOH) ! bmu^0_OH      
      
  bmupls=-dlog(xplsbulk/vpls/xsolbulk**vpls) ! bmu+
  bmumin=-dlog(xminbulk/vmin/xsolbulk**vmin) ! bmu-




  !.................

   

  nc_pore=ihm2*iRpore ! # volume cells inside pore
  
  nc_mem=ihm2*dimr-nc_pore ! # volume cells inside membrane
  
  nc_tot=dimr*dimz ! # of cells including membrane+solution
  
  nc_sol=nc_tot-nc_mem ! # of solution cells
  
  nc_surf=(dimr-iRpore)+ihm2 ! # of surface cells

  



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

     f_Is1(ir)=fb_I/(1d0-fb_I)*dexp(-psis1(ir)*z_I)      ! grado de disociacion del lipido ionizable
     f_Is1(ir)=f_Is1(ir)/(1d0+f_Is1(ir))             ! grado de disociacion del lipido ionizable
    
   enddo
  

   do ir=iRpore+1,dimr

      faux=fb_I/(1d0-fb_I)
      
      faux=(1d0+faux*dexp(-psis1(ir)*z_I))/(1d0+faux*dexp(-psis1(dimr)*z_I))

      xIs1(ir)=xb_I*(xNs1(ir)/(1d0-xb_I))**(a_I/a_N)*faux

   enddo


  do iz=1,ihm2
     
     ii=ii+1
     
     xNs2(iz)=x(ii)    ! fraccion de area del lipido neutro en la superficie 2 , pared del cilindro
     
     psis2(iz)=x(ii+nc_surf)

     f_Is2(iz)=fb_I/(1d0-fb_I)*dexp(-psis2(iz)*z_I)
     f_Is2(iz)=f_Is2(iz)/(1d0+f_Is2(iz))

     faux=fb_I/(1d0-fb_I)     
     faux=(1d0+faux*dexp(-psis2(iz)*z_I))/(1d0+faux*dexp(-psis1(dimr)*z_I))

     xIs2(iz)=xb_I*(xNs2(iz)/(1d0-xb_I))**(a_I/a_N)*faux
     
  enddo
  
  
  
  !!!!!!!!!!!!
  
  
  
  
  
  
  
  
  
  
  r_min=lbound(naarz,4)
  r_max=ubound(naarz,4)

  h_min=lbound(naarz,6)
  h_max=ubound(naarz,6)

  iz_min=-h_max+1
  iz_max=dimz-h_min ! asummes h_min<0

  ir_max=dimr-r_min ! asummes r_min<0 
  
  
  iz_mxx=ubound(naarz,5)
  ir_mxx=ubound(naarz,3)
  
  

           


  !.......... calculate monomers vol fract ...........
  !.......... version 1
  
  rhop=0d0
  umf=0d0
  rhomin=0d0
  rhopls=0d0
  rhoH=0d0
  rhoOH=0d0
  
  Stf=0d0
  Stfb=0d0
  Stp=0d0
  Stpb=0d0
  bFms=0d0
  bFmsb=0d0
  bFchm=0d0
  bFbchm=0d0
  sumavnaa=0d0
  sumavnaabulk=0d0
  bFmup=0d0
  bFmupb=0d0
  bFmup2=0d0
  bFmup2b=0d0
  bFmul2=0d0
  bFmul2b=0d0
  Ue=0d0
  Ueb=0d0
  Ues=0d0
  Uesb=0d0



  
  

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

           
           if (iz>0) then
            rhop(ir,iz)=rhop(ir,iz)+prz*eDmup
            if (eDmup==0d0) then
              Stp=Stp+0d0
            else
               Stp=Stp+prz*eDmup*(dlog(prz*eDmup*vsol)-1d0)*(2*ir-1)
            endif
           endif
           
           do h=h_min,h_max
              do r=r_min,r_max ! Loop sobre x 
                 
                 hh=iz+h  
                 rr=ir+r
                 
                 if (rr<1) cycle
                 
                 

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

        rhomin(ir,iz)=xmin
        rhopls(ir,iz)=xpls
        rhoH(ir,iz)=xH
        rhoOH(ir,iz)=xOH
        
        
          Stf=Stf+xH/vH*(dlog(xH/vH)-1d0)*(2*ir-1)
          Stf=Stf+xOH/vOH*(dlog(xOH/vOH)-1d0)*(2*ir-1)
          Stf=Stf+xmin/vmin*(dlog(xmin/vmin)-1d0)*(2*ir-1)
          Stf=Stf+xpls/vpls*(dlog(xpls/vpls)-1d0)*(2*ir-1)
          Stf=Stf+xsol*(dlog(xsol)-1d0)*(2*ir-1)
        
        
        
          bFchm=bFchm+xH/vH*bmuH*(2*ir-1)
          bFchm=bFchm+xOH/vOH*bmuOH*(2*ir-1)
          bFchm=bFchm+xpls/vpls*bmupls*(2*ir-1)
          bFchm=bFchm+xmin/vmin*bmumin*(2*ir-1)
        



     enddo ! Fin Loop sobre coordenada global (cm)
  enddo    ! Fin Loop sobre coordenada global (cm)
  
  

  do iz=ihm2+1,iz_max ! Above pore
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

           enddo
           
           

              prz=dexp(lnprz)

          
           if (iz<=dimz .and. ir<=dimr) then
            rhop(ir,iz)=rhop(ir,iz)+prz*eDmup
            
            if (eDmup==0d0) then
              Stp=Stp+0d0
            else
              Stp=Stp+prz*eDmup*(dlog(prz*eDmup*vsol)-1d0)*(2*ir-1)
            endif
           endif

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

        rhomin(ir,iz)=xmin
        rhopls(ir,iz)=xpls
        rhoH(ir,iz)=xH
        rhoOH(ir,iz)=xOH
        

        Stf=Stf+xH/vH*(dlog(xH/vH)-1d0)*(2*ir-1)
        Stf=Stf+xOH/vOH*(dlog(xOH/vOH)-1d0)*(2*ir-1)
        Stf=Stf+xmin/vmin*(dlog(xmin/vmin)-1d0)*(2*ir-1)
        Stf=Stf+xpls/vpls*(dlog(xpls/vpls)-1d0)*(2*ir-1)
        Stf=Stf+xsol*(dlog(xsol)-1d0)*(2*ir-1)
        
        
        
        bFchm=bFchm+xH/vH*bmuH*(2*ir-1)
        bFchm=bFchm+xOH/vOH*bmuOH*(2*ir-1)
        bFchm=bFchm+xpls/vpls*bmupls*(2*ir-1)
        bFchm=bFchm+xmin/vmin*bmumin*(2*ir-1)
        
        
     enddo
  enddo
  
  
  do j=1,ntyi
     avnaa(j,:,:)=avnaa(j,:,:)*vaa(j)*eDmup*vsol
  enddo
  
  
  


!!!--------------------------------------------------------------



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do iz=1,ihm2 ! Inside Pore
     do ir=1,iRpore
        do j=1,ntyi
          sumavnaa=sumavnaa+avnaa(rhoq_i(j),ir,iz)*faa(j,ir,iz)*dlog(faa(j,ir,iz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          sumavnaa=sumavnaa+avnaa(rhoq_i(j),ir,iz)*(1-faa(j,ir,iz))*dlog(1-faa(j,ir,iz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          
          
          if (sign(1d0,zaa(j))>0) then
          
            bFmup=bFmup+avnaa(rhoq_i(j),ir,iz)*faa(j,ir,iz)*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH)/vaa(rhoq_i(j))/vsol*(2*ir-1)
            
            if (eDmup==0d0) then

              bFmup2=bFmup2+0d0

            else

              bFmup2=bFmup2+rhop(ir,iz)*(-dlog(eDmup))*(2*ir-1)

            endif
            
          else
          
            bFmup=bFmup+avnaa(rhoq_i(j),ir,iz)*(1-faa(j,ir,iz))*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH)/vaa(rhoq_i(j))/vsol*(2*ir-1)
            
            if (eDmup==0d0) then

              bFmup2=bFmup2+0d0

            else

              bFmup2=bFmup2+rhop(ir,iz)*(-dlog(eDmup)+cnaa(j)*dlog((1-faa(j,ir,iz))/faa(j,ir,iz)))*(2*ir-1)

            endif
            
          endif
        enddo
     enddo
  enddo
        
        
  do iz=ihm2+1,dimz ! Above pore
     do ir=1,dimr
        do j=1,ntyi
          sumavnaa=sumavnaa+avnaa(rhoq_i(j),ir,iz)*faa(j,ir,iz)*dlog(faa(j,ir,iz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          sumavnaa=sumavnaa+avnaa(rhoq_i(j),ir,iz)*(1-faa(j,ir,iz))*dlog(1-faa(j,ir,iz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          
          
          if (sign(1d0,zaa(j))>0) then
          
            bFmup=bFmup+avnaa(rhoq_i(j),ir,iz)*faa(j,ir,iz)*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH)/vaa(rhoq_i(j))/vsol*(2*ir-1)
            
            if (eDmup==0d0) then

              bFmup2=bFmup2+0d0

            else

              bFmup2=bFmup2+rhop(ir,iz)*(-dlog(eDmup))*(2*ir-1)

            endif
            
          else
          
            bFmup=bFmup+avnaa(rhoq_i(j),ir,iz)*(1-faa(j,ir,iz))*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH)/vaa(rhoq_i(j))/vsol*(2*ir-1)
            
            if (eDmup==0d0) then

              bFmup2=bFmup2+0d0

            else

              bFmup2=bFmup2+rhop(ir,iz)*(-dlog(eDmup)+cnaa(j)*dlog((1-faa(j,ir,iz))/faa(j,ir,iz)))*(2*ir-1)

            endif
            
          endif
        enddo
     enddo
  enddo



  eDmui=xb_I/((1-xb_I)*a_l)**a_I*(1-fb_I)

do ir=iRpore+1,dimr
  bFms=bFms+xNs1(ir)/a_N/a_l*(dlog(xNs1(ir)/a_N)-1)*(2*ir-1)*pi*deltar**2
  bFms=bFms+xIs1(ir)/a_I/a_l*((dlog(xIs1(ir)/a_I))-1)*(2*ir-1)*pi*deltar**2
  
  bFms=bFms+xIs1(ir)/a_I/a_l*f_Is1(ir)*dlog(f_Is1(ir))*(2*ir-1)*pi*deltar**2
  bFms=bFms+xIs1(ir)/a_I/a_l*(1-f_Is1(ir))*dlog(1-f_Is1(ir))*(2*ir-1)*pi*deltar**2
  
  bFmul2=bFmul2+xIs1(ir)/a_I/a_l*(-dlog(eDmui*a_l)-bmuH-dlog(10**(-pK_I)*Na*vsol))*(2*ir-1)*pi*deltar**2 
  
enddo

do iz=1,ihm2
  bFms=bFms+xNs2(iz)/a_N/a_l*(dlog(xNs2(iz)/a_N)-1)*2*pi*iRpore*deltar*deltaz
  bFms=bFms+xIs2(iz)/a_I/a_l*((dlog(xIs2(iz)/a_I))-1)*2*pi*iRpore*deltar*deltaz
  
  bFms=bFms+xIs2(iz)/a_I/a_l*f_Is2(iz)*dlog(f_Is2(iz))*2*pi*iRpore*deltar*deltaz
  bFms=bFms+xIs2(iz)/a_I/a_l*(1-f_Is2(iz))*dlog(1-f_Is2(iz))*2*pi*iRpore*deltar*deltaz
 
  bFmul2=bFmul2+xIs2(iz)/a_I/a_l*(-dlog(eDmui*a_l)-bmuH-dlog(10**(-pK_I)*Na*vsol))*2*pi*iRpore*deltar*deltaz 
   
enddo






  do iz=1,ihm2 ! Inside Pore
     do ir=1,iRpore


	!!! Derivadas backward
	kz=iz-1
        if (iz==1) kz=iz

        dpsiz2=psi(ir,iz)-psi(ir,kz)
        dpsiz2=dpsiz2**2/deltaz**2
         
        if (ir==1) then
           dpsir2=psi(ir,iz)-psi(ir,iz)
           dpsir2=dpsir2**2/deltar**2           

        else
           dpsir2=psi(ir,iz)-psi(ir-1,iz)
           dpsir2=dpsir2**2/deltar**2
           
        endif

        
        do j=1,ntyi           
           rhoq(ir,iz)=rhoq(ir,iz)+faa(j,ir,iz)*zaa(j)*avnaa(rhoq_i(j),ir,iz)/vaa(rhoq_i(j))
        enddo
               
        Ue=Ue+(rhoq(ir,iz)*psi(ir,iz)/vsol-0.5/4d0/pi/lb_w*(dpsiz2+dpsir2))*(2*ir-1)
     enddo
  enddo
  
  
  
  
  do iz=ihm2+1,dimz ! Above pore
     do ir=1,dimr
        

	!!! Derivadas backward
        kz=iz
        if (iz==dimz) kz=dimz
        
        dpsiz2=psi(ir,kz)-psi(ir,iz-1)
        dpsiz2=dpsiz2**2/deltaz**2
        
        if (iz==ihm2+1.and.ir>iRpore) then           
           dpsiz2=(psi(ir,iz)-psis1(ir))/deltaz 
           dpsiz2=dpsiz2**2        
        endif
        
        
        
         
        if (ir==1) then
           dpsir2=psi(ir,iz)-psi(ir,iz)
           dpsir2=dpsir2**2/deltar**2            

        else
           dpsir2=psi(ir,iz)-psi(ir-1,iz)
           dpsir2=dpsir2**2/deltar**2
           
        endif

        
        do j=1,ntyi           
           rhoq(ir,iz)=rhoq(ir,iz)+faa(j,ir,iz)*zaa(j)*avnaa(rhoq_i(j),ir,iz)/vaa(rhoq_i(j))
        enddo
            
            

        Ue=Ue+(rhoq(ir,iz)*psi(ir,iz)/vsol-0.5/4d0/pi/lb_w*(dpsiz2+dpsir2))*(2*ir-1)
        
     enddo
  enddo
  

  
  do iz=1,ihm2              ! Inside membrane
     do ir=iRpore+1,dimr
                    


	!!! Derivadas backward
        if (iz==1) then
           dpsiz2=psi(ir,iz)-psi(ir,iz)
           dpsiz2=dpsiz2**2/deltaz**2
        else
          dpsiz2=psi(ir,iz)-psi(ir,iz-1)
          dpsiz2=dpsiz2**2/deltaz**2
        endif




        if (ir==iRpore+1) then
           dpsir2=psi(ir,iz)-psis2(iz)
           dpsir2=dpsir2**2/deltar**2
        else
           dpsir2=psi(ir,iz)-psi(ir-1,iz)
           dpsir2=dpsir2**2/deltar**2
        endif

 
        Ue=Ue-0.5/4d0/pi/lb_m*(dpsiz2+dpsir2)*(2*ir-1)                   

     enddo
  enddo

  
  

  
  do ir=iRpore+1,dimr ! Planar surface
     
     Ues=Ues+xIs1(ir)/a_I/a_l*f_Is1(ir)*z_I*psis1(ir)*(2*ir-1)*pi*deltar**2
     
  enddo
  
  
  do iz=1,ihm2 ! Pore surface
     
     Ues=Ues+xIs2(iz)/a_I/a_l*f_Is2(iz)*z_I*psis2(iz)*2*pi*iRpore*deltar*deltaz
      
  enddo


  !!!--------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(spH,"(f6.3)")pH
  if (pH<10d0.and.spH/="10.000") then
     write(spH,"(f5.3)")pH
     spH="0"//spH
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="umf-"//spH//".dat",unit=33)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        write(33,*)(dble(iz)-0.5)*deltaz,-dlog(rhop(ir,iz)/rhopbulk)
     enddo
     write(33,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        write(33,*)(dble(iz)-0.5)*deltaz,-dlog(rhop(ir,iz)/rhopbulk)
     enddo
     write(33,*)
  enddo
       
  close(33)
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="rhop-"//spH//".dat",unit=34)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        write(34,*)(dble(iz)-0.5)*deltaz,rhop(ir,iz)
     enddo
     write(34,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        write(34,*)(dble(iz)-0.5)*deltaz,rhop(ir,iz)
     enddo
     write(34,*)
  enddo
       
  close(34)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="rhomin-"//spH//".dat",unit=35)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        write(35,*)(dble(iz)-0.5)*deltaz,rhomin(ir,iz)/Na/vmin/vsol
     enddo
     write(35,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        write(35,*)(dble(iz)-0.5)*deltaz,rhomin(ir,iz)/Na/vmin/vsol
     enddo
     write(35,*)
  enddo
       
  close(35)
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="rhopls-"//spH//".dat",unit=36)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        write(36,*)(dble(iz)-0.5)*deltaz,rhopls(ir,iz)/Na/vpls/vsol
     enddo
     write(36,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        write(36,*)(dble(iz)-0.5)*deltaz,rhopls(ir,iz)/Na/vpls/vsol
     enddo
     write(36,*)
  enddo
       
  close(36)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="rhoH-"//spH//".dat",unit=30)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        write(30,*)(dble(iz)-0.5)*deltaz,rhoH(ir,iz)/Na/vH/vsol
     enddo
     write(30,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        write(30,*)(dble(iz)-0.5)*deltaz,rhoH(ir,iz)/Na/vH/vsol
     enddo
     write(30,*)
  enddo
       
  close(30)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="rhoOH-"//spH//".dat",unit=31)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        write(31,*)(dble(iz)-0.5)*deltaz,rhoOH(ir,iz)/Na/vOH/vsol
     enddo
     write(31,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        write(31,*)(dble(iz)-0.5)*deltaz,rhoOH(ir,iz)/Na/vOH/vsol
     enddo
     write(31,*)
  enddo
       
  close(31)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="rhoq-"//spH//".dat",unit=37)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        write(37,*)(dble(iz)-0.5)*deltaz,rhoq(ir,iz)
     enddo
     write(37,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        write(37,*)(dble(iz)-0.5)*deltaz,rhoq(ir,iz)
     enddo
     write(37,*)
  enddo
       
  close(37)
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="gamma-"//spH//".dat",unit=66) !isotermas de adsorción
  
  gamma=0d0
 
  do ir=1,iRpore
     do iz=1,dimz
        gamma=gamma+(rhop(ir,iz)-rhopbulk)*(2*ir-1)
     enddo
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
        gamma=gamma+(rhop(ir,iz)-rhopbulk)*(2*ir-1)
     enddo
  enddo
       
  gamma=gamma*pi*deltaz*deltar**2/vsol

  write(66,*) gamma
  close(66)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="gamma2-"//spH//".dat",unit=67) !isotermas de adsorción
  
  gamma=0d0
 
  do ir=1,dimr
     do iz=ihm2+1,dimz
        gamma=gamma+(rhop(ir,iz)-rhop(dimr,iz))*(2*ir-1)
     enddo
  enddo
  
  do ir=1,iRpore
     do iz=1,ihm2
        gamma=gamma+rhop(ir,iz)*(2*ir-1)
     enddo
  enddo
       
  gamma=gamma*pi*deltaz*deltar**2/vsol

  write(67,*) gamma
  close(67)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="q-"//spH//".dat",unit=3)
   
  
 
  do ir=1,iRpore
   volir=(2*ir-1)*deltar**2*deltaz*pi
    do iz=1,dimz
      q=0d0
      sumavnaa=0d0
      do j=1,ntyi
        q=q+cnaa(j)*faa(j,ir,iz)*zaa(j)
        sumavnaa=sumavnaa+avnaa(j,ir,iz)
      enddo
      write(3,*)(dble(iz)-0.5)*deltaz,q
    enddo
    write(3,*)
  enddo
  
  do ir=iRpore+1,dimr
    volir=(2*ir-1)*deltar**2*deltaz*pi
    do iz=ihm2+1,dimz
      q=0d0
      sumavnaa=0d0
      do j=1,ntyi
        q=q+cnaa(j)*faa(j,ir,iz)*zaa(j)
        sumavnaa=sumavnaa+avnaa(j,ir,iz)
      enddo
      write(3,*)(dble(iz)-0.5)*deltaz,q
    enddo
    write(3,*)
  enddo
  close(3)
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(file="qpore-"//spH//".dat",unit=4)
   
 
  qpore=0d0
  do ir=1,iRpore
    do iz=1,ihm2
      qpore=qpore+rhoq(ir,iz)*(2*ir-1)
    enddo
  enddo




  qlip=0d0
  do iz=1,ihm2
    qlip=qlip+xIs2(iz)/a_I/a_l*f_Is2(iz)*z_I*2*pi*iRpore*deltar*deltaz
  enddo



  
  write(*,*) qpore*pi*deltaz*deltar**2/vsol+qlip

  write(4,*)(dble(iRpore)-0.5)*deltar,qpore*pi*deltaz*deltar**2/vsol+qlip
  write(4,*)

  do iz=ihm2+1,dimz

    do ir=1,iRpore

      qpore=qpore+rhoq(ir,iz)*(2*ir-1)

    enddo

    qpore2=qpore
    qlip2=qlip
    write(4,*)(dble(iRpore)-0.5)*deltar,qpore*pi*deltaz*deltar**2/vsol+qlip
    do ir=iRpore+1,dimr

      do jz=ihm2+1,iz
        qpore2=qpore2+rhoq(ir,jz)*(2*ir-1)
      enddo
      qlip2=qlip2+xIs1(ir)/a_I/a_l*f_Is1(ir)*z_I*(2*ir-1)*pi*deltar**2
      write(4,*)(dble(ir)-0.5)*deltar,qpore2*pi*deltaz*deltar**2/vsol+qlip2
    enddo
    write(4,*)
  enddo


  close(4)



  q_sol=0d0

  do ir=1,iRpore
     do iz=1,ihm2
        q_sol=q_sol+rhoq(ir,iz)*(2*ir-1)
     enddo
  enddo

  do ir=1,dimr
    do iz=ihm2+1,dimz
       q_sol=q_sol+rhoq(ir,iz)*(2*ir-1)
    enddo
  enddo

  q_sol=q_sol*pi*deltaz*deltar**2/vsol



  q_solRef=0d0
  
  do iz=ihm2+1,dimz
     q_solRef=q_solRef+rhoq(dimr,iz)*(2*dimr-1)
  enddo

  q_solRef=q_solRef*pi*deltaz*deltar**2/vsol




  q_lip=0d0

  do iz=1,ihm2
     q_lip=q_lip+xIs2(iz)/a_I/a_l*f_Is2(iz)*z_I*dble(2*iRpore)*pi*deltar*deltaz
  enddo
  
  do ir=iRpore+1,dimr
    q_lip=q_lip+xIs1(ir)/a_I/a_l*f_Is1(ir)*z_I*dble(2*ir-1)*pi*deltar**2
  enddo



  q_lipRef=xIs1(dimr)/a_I/a_l*f_Is1(dimr)*z_I*dble(2*dimr-1)*pi*deltar**2









  write(*,*)
  write(*,*)"System:     ",real(q_sol),real(q_lip),real(q_sol+q_lip),&
       &real(abs((q_sol+q_lip)/q_sol))
  write(*,*)"Ref. System:",real(q_solRef),real(q_lipRef),real(q_solRef+q_lipRef),&
       &real(abs((q_solRef+q_lipRef)/q_solRef))

  

  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
 do j=1,ntyi 
 write(sj,"(i2.1)")j
 if (j<10) then

    write(sj,"('0'i1)")j
  endif 

 
  open(file="avnaa"//sj//"-"//spH//".dat",unit=80+j) ! write electrostatic potential for 3D plotting
 
  do ir=1,dimr
    do iz=1,dimz
      write(80+j,*)(dble(iz)-0.5)*deltaz,avnaa(j,ir,iz)
    enddo
    write(80+j,*)
  enddo
  close(80+j)
 enddo
 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
if (.false.) then
  open(file="psi3d-"//spH//".dat",unit=52) ! write electrostatic potential for 3D plotting
  
  do ir=1,iRpore 
     do iz=1,dimz
        write(52,*)(dble(ir)-0.5)*deltar,(dble(iz)-0.5)*deltaz,psi(ir,iz)
     enddo
  enddo


  do ir=iRpore+1,dimr

     do iz=1,ihm2
        write(52,*)(dble(ir)-0.5)*deltar,(dble(iz)-0.5)*deltaz,psi(ir,iz)
     enddo

     write(52,*)(dble(ir)-0.5)*deltar,dble(ihm2)*deltaz,psis1(ir)

     do iz=ihm2+1,dimz
        write(52,*)(dble(ir)-0.5)*deltar,(dble(iz)-0.5)*deltaz,psi(ir,iz)
     enddo
  enddo

  close(52)
  
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 


 open(file="psi_z-"//spH//".dat",unit=54)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        
	   if (iz==1) then
	     write(54,*)-(dble(iz)-0.5)*deltaz,psi(ir,iz)
	   endif

       write(54,*)(dble(iz)-0.5)*deltaz,psi(ir,iz)
          
     enddo
     write(54,*)
  enddo
  
  do ir=iRpore+1,dimr

    do iz=1,ihm2

	   if (iz==1) then
	     write(54,*)-(dble(iz)-0.5)*deltaz,psi(ir,iz)
	   endif

       write(54,*)(dble(iz)-0.5)*deltaz,psi(ir,iz)
    enddo
  
    write(54,*) ihm2*deltaz,psis1(ir)

    do iz=ihm2+1,dimz

      write(54,*)(dble(iz)-0.5)*deltaz,psi(ir,iz)
     enddo
     write(54,*)
  enddo
       
  close(54)





  open(file="psi_r-"//spH//".dat",unit=154)
  

  do iz=1,dimz
     
     do ir=1,iRpore
        
	if (ir==1) then
	  write(154,*)-(dble(ir)-0.5)*deltar,psi(ir,iz)
	endif

        write(154,*)(dble(ir)-0.5)*deltar,psi(ir,iz)
          
     enddo


     do ir=iRpore+1,dimr

        write(154,*)(dble(ir)-0.5)*deltar,psi(ir,iz)

     enddo

     write(154,*)

  enddo
       
  close(154)
  
  
  


  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  open(file="xwater-"//spH//".dat",unit=51)
  
 
  do ir=1,iRpore
     do iz=1,dimz
        
        write(51,*)(dble(iz)-0.5)*deltaz,phiw(ir,iz)
          
     enddo
     write(51,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz

        write(51,*)(dble(iz)-0.5)*deltaz,phiw(ir,iz)
             
     enddo
     write(51,*)
  enddo
       
  close(51)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  open(file="xlip-"//spH//".dat",unit=53) ! write surface density of lipids
  
  do iz=1,ihm2
     write(53,*)(iz-0.5)*deltaz,xNs2(iz),xIs2(iz),xIs2(iz)*(1d0-f_Is2(iz)),xIs2(iz)*f_Is2(iz)
  enddo

  write(53,*)

  do ir=iRpore+1,dimr  
     write(53,*)(ir-0.5)*deltar,xNs1(ir),xIs1(ir),xIs1(ir)*(1d0-f_Is1(ir)),xIs1(ir)*f_Is1(ir)
  enddo
  
  close(53)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ph0=log10(Na*vH*vsol)  

  open(file="ph-"//spH//".dat",unit=11)
  do ir=1,iRpore
     do iz=1,dimz
           
        epsi=dexp(-psi(ir,iz))   
        bpsi=psi(ir,iz)
        xsol=phiw(ir,iz)
        
        ysol=xsol/xsolbulk
        
        xH=xHbulk*ysol**vH*epsi**zH          ! volume fraction H3O+ locales
        xOH=xOHbulk*ysol**vOH*epsi**zOH      ! volume fraction OH-   locales
        xmin=xminbulk*ysol**vmin*epsi**zmin  ! volume fraction anion  
        xpls=xplsbulk*ysol**vpls*epsi**zpls
        
        
        write(11,*)(dble(iz)-0.5)*deltaz,-log10(xH)+ph0
          
     enddo
     write(11,*)
  enddo
  
  do ir=iRpore+1,dimr
     do iz=ihm2+1,dimz
     
        epsi=dexp(-psi(ir,iz))   
        bpsi=psi(ir,iz)
        xsol=phiw(ir,iz)
        
        ysol=xsol/xsolbulk
        
        xH=xHbulk*ysol**vH*epsi**zH          ! volume fraction H3O+ locales
        xOH=xOHbulk*ysol**vOH*epsi**zOH      ! volume fraction OH-   locales
        xmin=xminbulk*ysol**vmin*epsi**zmin  ! volume fraction anion  
        xpls=xplsbulk*ysol**vpls*epsi**zpls

        write(11,*)(dble(iz)-0.5)*deltaz,-log10(xH)+ph0
             
     enddo
     write(11,*)
  enddo
  
  close(11)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!! reference system
  rhop=0d0
  rhoq=0d0

  do iz=ihm2+1,iz_max  ! Above membrane
     do ir=1,ir_max

        if (iz>dimz) then           
            epsi=dexp(-psi(dimr,dimz))
            bpsi=psi(dimr,dimz)
            xsol=phiw(dimr,dimz)     
        else
            epsi=dexp(-psi(dimr,iz))
            bpsi=psi(dimr,iz)
            xsol=phiw(dimr,iz)    
        endif
        
        ysol=xsol/xsolbulk 
        
        xH=xHbulk*ysol**vH*epsi**zH          ! volume fraction H3O+ locales
        xOH=xOHbulk*ysol**vOH*epsi**zOH      ! volume fraction OH-   locales
        xmin=xminbulk*ysol**vmin*epsi**zmin  ! volume fraction anion  
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
                 
                 if (hh<=ihm2) cycle
                    

                ! endif
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
           

              prz=dexp(lnprz)

            
           if (iz<=dimz .and. ir<=dimr) then
            rhop(ir,iz)=rhop(ir,iz)+prz*eDmup
            
            if (eDmup==0d0) then
              Stpb=Stpb+0d0
            else
              Stpb=Stpb+prz*eDmup*(dlog(prz*eDmup*vsol)-1d0)*(2*ir-1)
            endif
           endif

        enddo
        
        if (iz>dimz.or.ir>dimr) cycle
        
                   
        xtot(ir,iz)=xH+xOH+xpls+xmin+xsol      ! volume fraction de todo menos el peptido
        
        ! Densidad de carga de todo menos la contribucion del peptido
        rhoq(ir,iz)=xH*zH/vH+xOH*zOH/vOH+xpls*zpls/vpls+xmin*zmin/vmin 




        Stfb=Stfb+xH/vH*(dlog(xH/vH)-1d0)*(2*ir-1)
        Stfb=Stfb+xOH/vOH*(dlog(xOH/vOH)-1d0)*(2*ir-1)
        Stfb=Stfb+xmin/vmin*(dlog(xmin/vmin)-1d0)*(2*ir-1)
        Stfb=Stfb+xpls/vpls*(dlog(xpls/vpls)-1d0)*(2*ir-1)
        Stfb=Stfb+xsol*(dlog(xsol)-1d0)*(2*ir-1)

        bFbchm=bFbchm+xH/vH*bmuH*(2*ir-1)
        bFbchm=bFbchm+xOH/vOH*bmuOH*(2*ir-1)
        bFbchm=bFbchm+xpls/vpls*bmupls*(2*ir-1)
        bFbchm=bFbchm+xmin/vmin*bmumin*(2*ir-1)

     enddo
  enddo

  do iz=ihm2+1,dimz ! Above membrane
     do ir=1,dimr
        do j=1,ntyi
          
          sumavnaabulk=sumavnaabulk+avnaa(rhoq_i(j),dimr,iz)*faa(j,dimr,iz)*dlog(faa(j,dimr,iz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          sumavnaabulk=sumavnaabulk+avnaa(rhoq_i(j),dimr,iz)*(1-faa(j,dimr,iz))*dlog(1-faa(j,dimr,iz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          
          
          if (sign(1d0,zaa(j))>0) then
          
            bFmupb=bFmupb+avnaa(rhoq_i(j),dimr,iz)*faa(j,dimr,iz)*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH)/vaa(rhoq_i(j))/vsol*(2*ir-1)
            
            if (eDmup==0d0) then

              bFmup2b=bFmup2b+0d0

            else

              bFmup2b=bFmup2b+rhop(dimr,iz)*(-dlog(eDmup))*(2*ir-1)

            endif
            
          else
          
            bFmupb=bFmupb+avnaa(rhoq_i(j),dimr,iz)*(1-faa(j,dimr,iz))*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH) &
		/vaa(rhoq_i(j))/vsol*(2*ir-1)

            if (eDmup==0d0) then

              bFmup2b=bFmup2b+0d0

            else

              bFmup2b=bFmup2b+rhop(dimr,iz)*(-dlog(eDmup)+cnaa(j)*dlog((1-faa(j,dimr,iz))/faa(j,dimr,iz)))*(2*ir-1)

            endif
          endif         
        enddo
     enddo
  enddo


  do iz=ihm2+1,dimz
    do j=1,ntyi           
      rhoq(dimr,iz)=rhoq(dimr,iz)+faa(j,dimr,iz)*zaa(j)*avnaa(rhoq_i(j),dimr,iz)/vaa(rhoq_i(j))
    enddo
  enddo



  do iz=ihm2+1,dimz ! Above membrane
     do ir=1,dimr


	!!! Derivadas backward
        kz=iz
        if (iz==dimz) kz=dimz
        
        dpsiz2=psi(dimr,kz)-psi(dimr,iz-1)
        dpsiz2=dpsiz2**2/deltaz**2
        
        if (iz==ihm2+1) then           
           dpsiz2=(psi(dimr,iz)-psis1(dimr))/deltaz 
           dpsiz2=dpsiz2**2        
        endif
        
    
        dpsir2=psi(dimr,iz)-psi(dimr,iz) !Just making the logic explicit here
        dpsir2=dpsir2**2/deltar**2


        
            

        Ueb=Ueb+(rhoq(dimr,iz)*psi(dimr,iz)/vsol-0.5/4d0/pi/lb_w*(dpsiz2+dpsir2))*(2*ir-1)
        
     enddo
  enddo


  do iz=1,ihm2              ! Inside membrane
     do ir=1,dimr
                    

	!!! Derivadas backward
        if (iz==1) then
           dpsiz2=psi(dimr,iz)-psi(dimr,iz)
           dpsiz2=dpsiz2**2/deltaz**2
        else
          dpsiz2=psi(dimr,iz)-psi(dimr,iz-1)
          dpsiz2=dpsiz2**2/deltaz**2
        endif



        dpsir2=psi(dimr,iz)-psi(dimr,iz) !Just making the logic explicit here
        dpsir2=dpsir2**2/deltar**2
 
 
        Ueb=Ueb-0.5/4d0/pi/lb_m*(dpsiz2+dpsir2)*(2*ir-1)                   

     enddo
  enddo




  eDmui=xb_I/((1-xb_I)*a_l)**a_I*(1-fb_I)


  do ir=1,dimr ! Planar surface
    bFmsb=bFmsb+xNs1(dimr)/a_N/a_l*(dlog(xNs1(dimr)/a_N)-1)*(2*ir-1)*pi*deltar**2
    bFmsb=bFmsb+xIs1(dimr)/a_I/a_l*((dlog(xIs1(dimr)/a_I))-1)*(2*ir-1)*pi*deltar**2
  
    bFmsb=bFmsb+xIs1(dimr)/a_I/a_l*f_Is1(dimr)*dlog(f_Is1(dimr))*(2*ir-1)*pi*deltar**2
    bFmsb=bFmsb+xIs1(dimr)/a_I/a_l*(1-f_Is1(dimr))*dlog(1-f_Is1(dimr))*(2*ir-1)*pi*deltar**2
  
    bFmul2b=bFmul2b+xIs1(dimr)/a_I/a_l*(-dlog(eDmui*a_l)-bmuH-dlog(10**(-pK_I)*Na*vsol))*(2*ir-1)*pi*deltar**2  

    Uesb=Uesb+xIs1(dimr)/a_I/a_l*f_Is1(dimr)*z_I*psis1(dimr)*(2*ir-1)*pi*deltar**2
  enddo

  !!!--------------------------------------------------------------




  bFcp=sumavnaa*pi*deltaz*deltar**2
  bFcpb=sumavnaabulk*pi*deltaz*deltar**2
  bFcp=bFcp-bFcpb
  
  
  bFmup=bFmup*pi*deltaz*deltar**2
  bFmupb=bFmupb*pi*deltaz*deltar**2
  bFmup=bFmup-bFmupb
  bFmup=-bFmup
  
  bFmup2=bFmup2*pi*deltaz*deltar**2
  bFmup2b=bFmup2b*pi*deltaz*deltar**2
  bFmup2=bFmup2-bFmup2b
  bFmup2=-bFmup2



  Ue=Ue*pi*deltaz*deltar**2
  Ueb=Ueb*pi*deltaz*deltar**2
  Ue=Ue-Ueb

  bFms=bFms-bFmsb

  
  bFmul2=bFmul2-bFmul2b
  bFmul2=-bFmul2
  
  Stf=Stf*pi*deltaz*deltar**2/vsol
  Stfb=Stfb*pi*deltaz*deltar**2/vsol
  Stf=Stf-Stfb
  
  bFchm=bFchm*pi*deltaz*deltar**2/vsol
  bFbchm=bFbchm*pi*deltaz*deltar**2/vsol
  bFchm=bFchm-bFbchm
  

  Stp=Stp*pi*deltaz*deltar**2/vsol
  Stpb=Stpb*pi*deltaz*deltar**2/vsol
  Stp=Stp-Stpb
  
  Ues=Ues-Uesb
  
  voltot=pi*(dimz*deltaz)*(dimr*deltar)**2


  FreeE=Stf+Stp+bFcp+bFms+Ue+Ues+bFchm+bFmup+bFmup2+bFmul2 !mitad de la energía libre


  open(file="freeE-"//spH//".dat",unit=39)
    write(39,*)"Stf:",Stf
    write(39,*)"Stp:",Stp
    write(39,*)"bFcp:",bFcp
    write(39,*)"bFms:",bFms
    write(39,*)"Ue:",Ue
    write(39,*)"Ues:",Ues
    write(39,*)"bFchm:",bFchm
    write(39,*)"bFmup:",bFmup
    write(39,*)"bFmup2:",bFmup2
    write(39,*)"bFmul2:",bFmul2    
    write(39,*)
    write(39,*)
    write(39,*)"FreeE:",FreeE  
  close(39)
  
  open(file="energy-"//spH//".dat",unit=38)
    write(38,*)"FreeE:",FreeE  
  close(38)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Free Energy of Water cilinder

  rhop=0d0
  umf=0d0
  
  Stf=0d0
  Stfb=0d0
  Stp=0d0
  Stpb=0d0
  bFms=0d0
  bFmsb=0d0
  bFchm=0d0
  bFbchm=0d0
  sumavnaa=0d0
  sumavnaabulk=0d0
  bFmup=0d0
  bFmupb=0d0
  bFmup2=0d0
  bFmup2b=0d0
  bFmul2=0d0
  bFmul2b=0d0
  Ue=0d0
  Ueb=0d0
  Ues=0d0
  Uesb=0d0
  avnaa=0d0

  ! Inside Pore
  do iz=1,ihm2-h_min
     do ir=1,iRpore-r_min,1

        epsi=dexp(-psi(ir,dimz))
        bpsi=psi(ir,dimz)
        xsol=phiw(ir,dimz)    
           
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

                 hh=dimz
                 
                 if (hh<=ihm2.and.rr>iRpore) cycle
                    
                 jz=iz
                 jr=ir
                 
                 if (iz>iz_mxx) jz=iz_mxx
                 if (ir>ir_mxx) jr=ir_mxx
                    
                 

                 do j=1,ntyi
                  lnprz=lnprz+naarz(i,j,jr,r,jz,h)*dlog(phiw(rr,hh))*vaa(j)/dble(ntht)
                    lnprz=lnprz-naarz(i,j,jr,r,jz,h)*dlog(1d0-faa(j,rr,hh))/dble(ntht)
                 enddo
                 
              enddo
           enddo
           
          
              prz=dexp(lnprz)
          
           if (iz<=dimz .and. ir<=dimr) then
            rhop(ir,iz)=rhop(ir,iz)+prz*eDmup
            
            if (eDmup==0d0) then
              Stp=Stp+0d0
            else
              Stp=Stp+prz*eDmup*(dlog(prz*eDmup*vsol)-1d0)*(2*ir-1)
            endif
           endif

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
        
        if (iz>ihm2.or.ir>iRpore) cycle
        
    
                   
        xtot(ir,iz)=xH+xOH+xpls+xmin+xsol
        
        rhoq(ir,iz)=xH*zH/vH+xOH*zOH/vOH+xpls*zpls/vpls+xmin*zmin/vmin 
        

        Stf=Stf+xH/vH*(dlog(xH/vH)-1d0)*(2*ir-1)
        Stf=Stf+xOH/vOH*(dlog(xOH/vOH)-1d0)*(2*ir-1)
        Stf=Stf+xmin/vmin*(dlog(xmin/vmin)-1d0)*(2*ir-1)
        Stf=Stf+xpls/vpls*(dlog(xpls/vpls)-1d0)*(2*ir-1)
        Stf=Stf+xsol*(dlog(xsol)-1d0)*(2*ir-1)
        
        
        
        bFchm=bFchm+xH/vH*bmuH*(2*ir-1)
        bFchm=bFchm+xOH/vOH*bmuOH*(2*ir-1)
        bFchm=bFchm+xpls/vpls*bmupls*(2*ir-1)
        bFchm=bFchm+xmin/vmin*bmumin*(2*ir-1)
        

     enddo
  enddo

  do j=1,ntyi
     avnaa(j,:,:)=avnaa(j,:,:)*vaa(j)*eDmup*vsol
  enddo


  do iz=1,ihm2 ! Inside Pore
     do ir=1,iRpore
        do j=1,ntyi
          sumavnaa=sumavnaa+avnaa(rhoq_i(j),ir,dimz)*faa(j,ir,dimz)*dlog(faa(j,ir,dimz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          sumavnaa=sumavnaa+avnaa(rhoq_i(j),ir,dimz)*(1-faa(j,ir,dimz))*dlog(1-faa(j,ir,dimz))/vaa(rhoq_i(j))/vsol*(2*ir-1)
          
          
          if (sign(1d0,zaa(j))>0) then
          
            bFmup=bFmup+avnaa(rhoq_i(j),ir,dimz)*faa(j,ir,dimz)*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH)/vaa(rhoq_i(j))/vsol*(2*ir-1)
            
            if (eDmup==0d0) then

              bFmup2=bFmup2+0d0

            else

              bFmup2=bFmup2+rhop(ir,dimz)*(-dlog(eDmup))*(2*ir-1)

            endif
            
          else
          
            bFmup=bFmup+avnaa(rhoq_i(j),ir,dimz)*(1-faa(j,ir,dimz))*(dlog(10**(-pKaa(j))*Na*vsol)+bmuH)/vaa(rhoq_i(j))/vsol*(2*ir-1)
            
            if (eDmup==0d0) then

              bFmup2=bFmup2+0d0

            else

              bFmup2=bFmup2+rhop(ir,dimz)*(-dlog(eDmup)+cnaa(j)*dlog((1-faa(j,ir,dimz))/faa(j,ir,dimz)))*(2*ir-1)

            endif
            
          endif
        enddo
     enddo
  enddo

  do iz=1,ihm2 ! Inside Pore
    do ir=1,iRpore
        
	

	!!! Derivadas backward  
        dpsiz2=psi(ir,dimz)-psi(ir,dimz)
        dpsiz2=dpsiz2**2/deltaz**2
        

        if (ir==1) then
           dpsir2=psi(ir,dimz)-psi(ir,dimz)
           dpsir2=dpsir2**2/deltar**2            

        else
           dpsir2=psi(ir,dimz)-psi(ir-1,dimz)
           dpsir2=dpsir2**2/deltar**2
           
        endif

        
        do j=1,ntyi           
           rhoq(ir,iz)=rhoq(ir,dimz)+faa(j,ir,dimz)*zaa(j)*avnaa(rhoq_i(j),ir,dimz)/vaa(rhoq_i(j))
        enddo
            
            

        Ue=Ue+(rhoq(ir,dimz)*psi(ir,dimz)/vsol-0.5/4d0/pi/lb_w*(dpsiz2+dpsir2))*(2*ir-1)
        
    enddo
  enddo


  !!!--------------------------------------------------------------




  bFcp=sumavnaa*pi*deltaz*deltar**2

  
  
  bFmup=bFmup*pi*deltaz*deltar**2
  bFmup=-bFmup
  
  bFmup2=bFmup2*pi*deltaz*deltar**2
  bFmup2=-bFmup2





  Ue=Ue*pi*deltaz*deltar**2


  
  bFmul2=-bFmul2
  
  Stf=Stf*pi*deltaz*deltar**2/vsol

  
  bFchm=bFchm*pi*deltaz*deltar**2/vsol
  

  Stp=Stp*pi*deltaz*deltar**2/vsol

  
  voltot=pi*(dimz*deltaz)*(dimr*deltar)**2


  FreeE=Stf+Stp+bFcp+bFms+Ue+Ues+bFchm+bFmup+bFmup2+bFmul2 !mitad de la energía libre


  open(file="wcfreeE-"//spH//".dat",unit=40)
    write(40,*)"Stf:",Stf
    write(40,*)"Stp:",Stp
    write(40,*)"bFcp:",bFcp
    write(40,*)"bFms:",bFms
    write(40,*)"Ue:",Ue
    write(40,*)"Ues:",Ues
    write(40,*)"bFchm:",bFchm
    write(40,*)"bFmup:",bFmup
    write(40,*)"bFmup2:",bFmup2
    write(40,*)"bFmul2:",bFmul2    
    write(40,*)
    write(40,*)
    write(40,*)"FreeE:",FreeE  
  close(40)
  
  open(file="wcenergy-"//spH//".dat",unit=41)
    write(41,*)"FreeE:",FreeE  
  close(41)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Free Energy of Membrane cilinder


  rhop=0d0
  umf=0d0
  
  Stf=0d0
  Stfb=0d0
  Stp=0d0
  Stpb=0d0
  bFms=0d0
  bFmsb=0d0
  bFchm=0d0
  bFbchm=0d0
  sumavnaa=0d0
  sumavnaabulk=0d0
  bFmup=0d0
  bFmupb=0d0
  bFmup2=0d0
  bFmup2b=0d0
  bFmul2=0d0
  bFmul2b=0d0
  Ue=0d0
  Ueb=0d0
  Ues=0d0
  Uesb=0d0



  do iz=1,ihm2              ! Inside membrane
     do ir=1,iRpore
                    


 	!!! Derivadas backward
        if (iz==1) then
          dpsiz2=psi(dimr,iz)-psi(dimr,iz)
          dpsiz2=dpsiz2**2/deltaz**2
        
        else
          dpsiz2=psi(dimr,iz)-psi(dimr,iz-1)
          dpsiz2=dpsiz2**2/deltaz**2
        endif



        dpsir2=psi(dimr,iz)-psi(dimr,iz)
        dpsir2=dpsir2**2/deltar**2
 
        Ue=Ue-0.5/4d0/pi/lb_m*(dpsiz2+dpsir2)*(2*ir-1)                   

     enddo
  enddo


  do ir=1,iRpore ! Planar surface
    bFms=bFms+xNs1(dimr)/a_N/a_l*(dlog(xNs1(dimr)/a_N)-1)*(2*ir-1)*pi*deltar**2
    bFms=bFms+xIs1(dimr)/a_I/a_l*((dlog(xIs1(dimr)/a_I))-1)*(2*ir-1)*pi*deltar**2
  
    bFms=bFms+xIs1(dimr)/a_I/a_l*f_Is1(dimr)*dlog(f_Is1(dimr))*(2*ir-1)*pi*deltar**2
    bFms=bFms+xIs1(dimr)/a_I/a_l*(1-f_Is1(dimr))*dlog(1-f_Is1(dimr))*(2*ir-1)*pi*deltar**2
  
    bFmul2=bFmul2+xIs1(dimr)/a_I/a_l*(-dlog(eDmui*a_l)-bmuH-dlog(10**(-pK_I)*Na*vsol))*(2*ir-1)*pi*deltar**2  

    Ues=Ues+xIs1(dimr)/a_I/a_l*f_Is1(dimr)*z_I*psis1(dimr)*(2*ir-1)*pi*deltar**2
  enddo
  !!!--------------------------------------------------------------




  bFcp=sumavnaa*pi*deltaz*deltar**2

  
  
  bFmup=bFmup*pi*deltaz*deltar**2
  bFmup=-bFmup
  
  bFmup2=bFmup2*pi*deltaz*deltar**2
  bFmup2=-bFmup2




  Ue=Ue*pi*deltaz*deltar**2


  
  bFmul2=-bFmul2
  
  Stf=Stf*pi*deltaz*deltar**2/vsol

  
  bFchm=bFchm*pi*deltaz*deltar**2/vsol
  

  Stp=Stp*pi*deltaz*deltar**2/vsol

  
  voltot=pi*(dimz*deltaz)*(dimr*deltar)**2


  FreeE=Stf+Stp+bFcp+bFms+Ue+Ues+bFchm+bFmup+bFmup2+bFmul2 !mitad de la energía libre


  open(file="mcfreeE-"//spH//".dat",unit=42)
    write(42,*)"Stf:",Stf
    write(42,*)"Stp:",Stp
    write(42,*)"bFcp:",bFcp
    write(42,*)"bFms:",bFms
    write(42,*)"Ue:",Ue
    write(42,*)"Ues:",Ues
    write(42,*)"bFchm:",bFchm
    write(42,*)"bFmup:",bFmup
    write(42,*)"bFmup2:",bFmup2
    write(42,*)"bFmul2:",bFmul2    
    write(42,*)
    write(42,*)
    write(42,*)"FreeE:",FreeE  
  close(42)
  
  open(file="mcenergy-"//spH//".dat",unit=43)
    write(43,*)"FreeE:",FreeE  
  close(43)


end subroutine freen
