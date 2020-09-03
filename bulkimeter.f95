
subroutine bulkimeter(vsol,vH,zH,vOH,zOH,vpls,zpls,vmin,zmin,pKw,pH,csalt,&
     &xsolbulk,xHbulk,xOHbulk,xplsbulk,xminbulk,rhopbulk,vaa,zaa,pKaa,naac,cnaa,ntyi)

  !#######################################################################
  !     This subroutine calculates bulk quantities
  !#######################################################################
  implicit none

  real(8),parameter:: Na=0.602214,pi=acos(-1.0d0)
  
  real(8),intent(in):: vsol,vH,zH,vOH,zOH,vpls,zpls,vmin,zmin
  real(8),intent(in):: pKw,pH,csalt
  real(8),intent(in):: rhopbulk,vaa(ntyi),zaa(ntyi),pKaa(ntyi)

  integer,intent(in):: naac,ntyi,cnaa(ntyi)

  real(8),intent(out):: xsolbulk,xHbulk,xOHbulk,xplsbulk,xminbulk

  integer:: i

  real(8) :: zxcs,pOH
  real(8) :: fbaa(ntyi)
  real(8) :: total_volume




  xHbulk=10**(-pH)*Na*vH  ! bulk H+ volume fraction
      
  xOHbulk=10**(pH-pKw)*Na*vOH ! bulk OH- volume fraction
      
  xplsbulk=csalt*Na*vpls
  xminbulk=csalt*Na*vmin


  do i=1,ntyi
     fbaa(i)=1d0/(1d0+10d0**(sign(1d0,zaa(i))*(pH-pKaa(i))))
  enddo


      
  zxcs=0d0
  total_volume=0d0

  ! add salt ions to preserve charge balance in the bulk
  
  zxcs=xHbulk*zH/vH+xOHbulk*zOH/vOH+xplsbulk*zpls/vpls+xminbulk*zmin/vmin ! Excess charge

  do i=1,ntyi
     zxcs=zxcs+rhopbulk*fbaa(i)*zaa(i)*dble(cnaa(i))
  enddo

  if (zxcs>0.d0) then     ! pH < 7 - Acidic
     xminbulk=xminbulk-zxcs/zmin*vmin ! add HCl
         
  elseif (zxcs<0.d0) then ! pH > 7 - Basic
     xplsbulk=xplsbulk-zxcs/zpls*vpls ! add NaOH
  endif


  xHbulk=xHbulk*vsol ! Volume fractions propiamente
  xOHbulk=xOHbulk*vsol
  xplsbulk=xplsbulk*vsol
  xminbulk=xminbulk*vsol

  do i=1,ntyi
     total_volume=total_volume+dble(cnaa(i))*vaa(i)
  enddo

  ! now apply packing constraint in the bulk to get xsolbulk  
  xsolbulk=1d0-xHbulk-xOHbulk-xplsbulk-xminbulk-rhopbulk*total_volume*vsol
   
   





  return
end subroutine bulkimeter
