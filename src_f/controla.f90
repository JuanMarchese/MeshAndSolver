subroutine control()
use def_solver
use def_variables
implicit none
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,qe
integer :: ns(nodpel),ns2d(2*nodpel)
integer  ::inbwt,nbwt,i,ij,j,nb,kk,ncase,inode,ipoin,jnode,ii,jj,j1,nicio,k,nno,npas,ne,inbwe,nbwe
double precision, allocatable :: ud(:),un(:)
integer, allocatable :: cont(:),ip(:),iu(:),iup(:),ju(:),iut(:),consim(:)
integer, allocatable :: cont2d(:),ip2d(:),iu2d(:),iup2d(:),ju2d(:),iut2d(:),consim2d(:)
integer:: unit_sist


! armo las estructuras logicas sparce
allocate(ad(nnodes),ia(nnodes+1),cont(nnodes+1),cx(nnodes+1),solucion(nnodes),rhs(nnodes),consim(nnodes))
allocate(ad2d(2*nnodes))
allocate(ia2d(2*nnodes+1),cont2d(2*nnodes+1),cx2d(2*nnodes+1),solucion2d(2*nnodes),rhs2d(2*nnodes),consim2d(2*nnodes))

rhs=0.0
ad=0.0
solucion=0.0
!ud=0.0
ia=0
!ip=0
cx=0
!iu=0
!iut=0
cont=0
consim=0
ne=nelements

      inbwt=0
      nbwt =0
      do kk=1,ne
        do i=1,nodpel
          ns(i)=conect(i,kk)
        enddo
        do  i=1,nodpel-1
		  ij=i+1
          do j=ij,nodpel
            nb=iabs(ns(i)-ns(j))
            if(nb.eq.0) then
			   write(*,*) 'elemento  ',kk,' tiene dos nodos identicos'
			   write(*,*) 'elemento  ',kk,' tiene dos nodos identicos'
			   stop
			endif
            if(nb.gt.nbwt) then
               inbwt=kk
               nbwt =nb
            endif
          enddo
        enddo
      enddo
      nbwt=nbwt+1

      !write(*,*) ' bandwidth: ',nbwt,'  en elemento  ',inbwt
   
   
! determino la forma simbolica del problema sparce. para eso ensamblo un sistema de unos y ceros



do kk=1,nnodes+1
  cx(kk)=0.0
  cont(kk)=0.0
  ia(kk) = 0
enddo


ncase=0
sigma_el=1.0
qe=1.0
ia(1)=1

do nno=1,nnodes
  consim=0
  do kk=1,nelements
     npas=0
     do i=1,nodpel
       ns(i)=conect(i,kk)
       if(nno==ns(i)) npas=1
     enddo
     if(npas==1) then   
       esm=1.0
 

       do i=1,nodpel
          ii=ns(i)
            do j=1,nodpel
              jj=ns(j)+1-ii
              if(jj.gt.0) then
	             if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
              endif  
	        enddo
        enddo

      endif
   enddo

   do i=1,nnodes
      if(consim(i)/=0) then
         cont(nno)=cont(nno)+1
      endif
   enddo
    ia(nno+1) = ia(nno) + cont(nno)
enddo




nonull = ia(nnodes+1)-1 !- nnodes


allocate (an(nonull),ja(nonull))
an=0.0
ja=0

!write(*,*) 'nonulos del sistema  ',nnodes,nonull,nonull*(1.0/real(nnodes))*(1/real(nnodes))

do nno=1,nnodes
  consim=0
  do kk=1,nelements
     npas=0
     do i=1,nodpel
       ns(i)=conect(i,kk)
       if(nno==ns(i)) npas=1
     enddo
     if(npas==1) then   
       esm=1.0
   
     do i=1,nodpel
       ii=ns(i)
       do j=1,nodpel
          jj=ns(j)+1-ii
          if(jj.gt.0) then
            if(nno==ns(i).and.nno/=ns(j)) consim(ns(j)) = consim(ns(j))+1          
          endif  
	   enddo
     enddo

   endif
  enddo
  do i=1,nnodes
      
      if(consim(i)/=0) then
          ja(ia(nno)+cx(nno)) =   i
          cx(nno)=cx(nno)+1
      endif
   enddo
enddo


! caso mecanico

rhs2d=0.0
ad2d=0.0
solucion2d=0.0
ia2d=0
cx2d=0
cont2d=0
consim2d=0
ne=nelements

      inbwe=0
      nbwe =0
      do kk=1,ne
        do i=1,nodpel
          ns2d(2*i-1)= 2*conect(i,kk)-1
          ns2d(2*i)  = 2*conect(i,kk)
        enddo
        do  i=1,2*nodpel-1
		  ij=i+1
          do j=ij,2*nodpel
            nb=iabs(ns2d(i)-ns2d(j))
            if(nb.eq.0) then 
			   write(*,*) 'elemento  ',kk,' tiene dos nodos identicos'
			   write(*,*) 'elemento  ',kk,' tiene dos nodos identicos'
			   stop
			endif    
            if(nb.gt.nbwe) then
               inbwe=kk
               nbwe =nb
            endif
          enddo
        enddo
      enddo
      nbwe=2*nbwe+1

      !write(*,*) ' bandwidth: ',nbwe,'  en elemento  ',inbwe
   
   
! determino la forma simbolica del problema sparce. para eso ensamblo un sistema de unos y ceros



do kk=1,2*nnodes+1
  cx2d(kk)=0.0
  cont2d(kk)=0.0
  ia2d(kk) = 0
enddo


ncase=0
sigma_el=1.0
qe=1.0
ia2d(1)=1

do nno=1,2*nnodes
  consim2d=0
  do kk=1,nelements
     npas=0
     do i=1,nodpel
       ns2d(2*i-1)=2*conect(i,kk)-1
       ns2d(2*i)=2*conect(i,kk)
       if(nno==ns2d(2*i-1) .or. nno==ns2d(2*i) ) npas=1
     enddo
     if(npas==1) then   
       do i=1,2*nodpel
            ii=ns2d(i)
            do j=1,2*nodpel
              jj=ns2d(j)+1-ii
              if(jj.gt.0) then
	             if(nno==ns2d(i) .and. nno/=ns2d(j)) consim2d(ns2d(j)) = consim2d(ns2d(j))+1          
              endif  
	        enddo
        enddo

      endif
   enddo

   do i=1,2*nnodes
      if(consim2d(i)/=0) then
         cont2d(nno)=cont2d(nno)+1
      endif
   enddo
    ia2d(nno+1) = ia2d(nno) + cont2d(nno)
enddo




nonull2d = ia2d(2*nnodes+1) !- nnodes


allocate (an2d(nonull2d),ja2d(nonull2d))
an2d=0.0
ja2d=0

!write(*,*) 'nonulos del sistema  ',2*nnodes,nonull2d,nonull2d*(1.0/real(2*nnodes))*(1/real(2*nnodes))

do nno=1,2*nnodes
  consim2d=0
  do kk=1,nelements
     npas=0
     do i=1,nodpel
       ns2d(2*i-1)=2*conect(i,kk)-1
       ns2d(2*i)=2*conect(i,kk)
       if(nno==ns2d(2*i-1) .or. nno==ns2d(2*i) ) npas=1
     enddo
     if(npas==1) then   
   
     do i=1,2*nodpel
       ii=ns2d(i)
       do j=1,2*nodpel
          jj=ns2d(j)+1-ii
          if(jj.gt.0) then
            if(nno==ns2d(i).and.nno/=ns2d(j)) consim2d(ns2d(j)) = consim2d(ns2d(j))+1          
          endif  
	   enddo
     enddo

   endif
  enddo
  
  do i=1,2*nnodes
     if(consim2d(i)/=0) then
          ja2d(ia2d(nno)+cx2d(nno)) = i
          cx2d(nno)=cx2d(nno)+1
     endif
  enddo

enddo

return
end

