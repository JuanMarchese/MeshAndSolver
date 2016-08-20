subroutine campo2d(nnodes,nelements,nodpel,solution,material,conect,coor_x,coor_y,sigmaext,sigmaint,sigmamem,grad_x,grad_y,gradxel_x,gradxel_y)


implicit none
integer :: nnodes,nelements,nodpel,conect(nelements,nodpel),material(nelements)
double precision :: grad_x(nnodes),grad_y(nnodes)
double precision :: gradxel_x(nelements),gradxel_y(nelements)

double precision :: solution(nnodes),coor_x(nnodes),coor_y(nnodes),sigmaext,sigmaint,sigmamem
! local
integer :: nope,ns(nodpel),jel,mat
double precision :: x(nodpel),y(nodpel),sol(nodpel),sigma_el,ex_el,ey_el,ez_el,ex(nodpel),ey(nodpel)
double precision:: b(3),c(3),deter,a,rmed
double precision:: pi=3.14159
double precision,allocatable :: gausspt(:),gausswt(:)
double precision,allocatable :: phi(:,:),dphi(:,:,:),gxcod(:,:),phidx(:,:,:),ctei(:) 

integer kk,jj,i,j,ii,k,nle,kgaus,i2,ngaus,ndimension,ndime,ndi
double precision:: t,s,sm,tm,sq,tp,ajaco(2,2),ajacoi(2,2),dphix(nodpel),dphiy(nodpel)


ngaus=2
ndimension=2

allocate(gausspt(ngaus),gausswt(ngaus),phi(2*ngaus,nodpel),dphi(ndimension,2*ngaus,nodpel),gxcod(ndimension,2*ngaus),phidx(ndimension,2*ngaus,nodpel),ctei(ngaus*2)) 


grad_x=0.0
grad_y=0.0
gradxel_x=0.0
gradxel_y=0.0

do jel=1,nelements

   mat=material(jel)        
   if(mat==3) then
       sigma_el = sigmaint
   elseif(mat==2) then
       sigma_el = sigmamem
   elseif(mat==1 ) then
       sigma_el = sigmaext
   endif
   
   
   do i=1,nodpel
        ns(i)=conect(jel,i)
	    j=ns(i)
        x(i)=coor_x(j)
        y(i)=coor_y(j)
        sol(i)=solution(j)
        ex(i)=0.0
        ey(i)=0.0
   enddo
        
   ex_el=0
   ey_el=0


  if(nodpel==3) then    
    !  evaluacion del gradiente en los elementos
       b(1)=y(2)-y(3)
       b(2)=y(3)-y(1)
       b(3)=y(1)-y(2)
       c(1)=x(3)-x(2)
       c(2)=x(1)-x(3)
       c(3)=x(2)-x(1)
   
       deter=x(2)*y(3)+x(3)*y(1)+x(1)*y(2)-x(2)*y(1)-x(3)*y(2)-x(1)*y(3)
       gradxel_x(jel)=-(b(1)*sol(1)+b(2)*sol(2)+b(3)*sol(3))/deter
       gradxel_y(jel)=-(c(1)*sol(1)+c(2)*sol(2)+c(3)*sol(3))/deter
   
    !   gradxel_x(jel)=sigma_el* gradxel_x(jel)*(-1.)
    !   gradxel_y(jel)=sigma_el* gradxel_y(jel)*(-1.)   

   
       do i=1,nodpel
            j=ns(i)
            grad_x(j) =  grad_x(j) + gradxel_x(i)/real(nodpel)
            grad_y(j) =  grad_y(j) + gradxel_y(i)/real(nodpel)
       enddo
    
 else


      gausspt(1)=-0.57735027
      gausspt(2)= 0.57735027
      gausswt(1)= 1.0
      gausswt(2)= 1.0 
 
       kgaus=0
   do kk=1,ngaus
     do jj=1,ngaus
         kgaus=kgaus+1

		 t = gausspt(kk)
         s = gausspt(jj)
      

             
              sm = 0.5*(1.0-s)
              tm = 0.5*(1.0-t)
              sq = 0.5*(1.0+s)
              tp = 0.5*(1.0+t)
              
              phi(kgaus,1)    = sm*tm
              dphi(1,kgaus,1) =-0.5*tm
              dphi(2,kgaus,1) =-0.5*sm
              phi(kgaus,2)     = sq*tm
              dphi(1,kgaus, 2) = 0.5*tm
              dphi(2,kgaus, 2) =-0.5*sq
              phi(kgaus,3)     = sq*tp
              dphi(1,kgaus, 3) = 0.5*tp
              dphi(2,kgaus, 3) = 0.5*sq
              phi(kgaus,4)     = sm*tp
              dphi(1,kgaus, 4) =-0.5*tp
              dphi(2,kgaus, 4) = 0.5*sm
              

!   jacobiano y su inversa
        ajaco=0.0

          do k=1,nodpel
              ajaco(1,1)=ajaco(1,1)+dphi(1,kgaus,k)*x(k)
              ajaco(1,2)=ajaco(1,2)+dphi(1,kgaus,k)*y(k)
              ajaco(2,1)=ajaco(2,1)+dphi(2,kgaus,k)*x(k)
              ajaco(2,2)=ajaco(2,2)+dphi(2,kgaus,k)*y(k)
          enddo

          deter= ajaco(1,1)*ajaco(2,2)-ajaco(1,2)*ajaco(2,1)


          write(*,*) 'Derminante es ',deter


          ajacoi(1,1)=ajaco(2,2)/deter
          ajacoi(2,2)=ajaco(1,1)/deter
          ajacoi(1,2)=-ajaco(1,2)/deter
          ajacoi(2,1)=-ajaco(2,1)/deter

          do i=1,nodpel
            dphix(i)=ajacoi(1,1)*dphi(1,kgaus,i) + ajacoi(1,2)*dphi(2,kgaus,i)
            dphiy(i)=ajacoi(2,1)*dphi(1,kgaus,i) + ajacoi(2,2)*dphi(2,kgaus,i) 
          enddo

         do i=1,nodpel
        
             ex(i)=ex(i)+ dphix(i)*sol(i)
             ey(i)=ey(i)+ dphiy(i)*sol(i)
             
             ex_el=ex_el+ dphix(i)*sol(i)
             ey_el=ey_el+ dphiy(i)*sol(i)
         enddo


      enddo

    enddo
      
    do i=1,nodpel
        j=ns(i)

        grad_x(j) =  grad_x(j) - ex(i)*0.25
        grad_y(j) =  grad_y(j) - ey(i)*0.25
    enddo
    gradxel_x(jel) = -ex_el*0.25
    gradxel_y(jel) = -ey_el*0.25
    



 endif  





enddo



end subroutine campo2d
