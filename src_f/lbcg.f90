
subroutine cg(np,ia,ja,a_spa,ad,b,x,tol,itmax,iter,err)
      implicit none 
	double precision a_spa(*),ad(*),x(*),b(*)
      integer ia(*),ja(*)
      double precision bnrm, znrm, bknum, bkden, xant, akden,bk,ak,zm1nrm,dxnrm,xnrm,alfa,b_norma,r_norma,aknum
	double precision tol,eps,err,funnorm
	integer itmax,iter,np,j,kk

      double precision, allocatable ::z(:),p(:),pp(:),r(:),zz(:),rr(:)

	allocate (p(np),pp(np),zz(np),z(np),r(np),rr(np)) 

      
	call matxvecsim(ia,ja,a_spa,ad,x,r,np)

      b_norma = 0.0
	do kk=1,np
        r(kk)  = b(kk) - r(kk)
        b_norma =b_norma  + b(kk)*b(kk) 
	enddo
      b_norma=sqrt(b_norma) 
      

	do kk=1,np
	  z(kk)=r(kk)/ad(kk)
        rr(kk)=z(kk)
	enddo

      znrm=funnorm(np,z)
      
      do while(err.gt.tol .and. iter.lt.itmax)

        iter=iter+1

	  bknum=0.d0
  	  do j=1,np
          bknum=bknum + z(j)*r(j)
        enddo
	  call matxvecsim(ia,ja,a_spa,ad,rr,z,np)
        
	  bkden=0.0
	  do j=1,np
          bkden=bkden + rr(j)*z(j)
        enddo
	   
        alfa= bknum/bkden    	  

        do kk=1,np
          x(kk)=x(kk)+alfa*rr(kk)          
	    r(kk)= r(kk) - alfa * z(kk)
	  enddo

     	  do kk=1,np
	    z(kk)=r(kk)/ad(kk)
        enddo

    	  akden=bknum
	  aknum=0.0
        do  j=1,np
          aknum = aknum+ z(j)*r(j)
        enddo

        ak=aknum/akden
          
        do j=1,np
          rr(j) = z(j)  + ak*rr(j)
        enddo

        r_norma = 0.0
	  do kk=1,np
          r_norma =r_norma  + r(kk)*r(kk)    
	  enddo
        r_norma=sqrt(r_norma) 

        err=r_norma/b_norma
        !write(2,*) iter,err,r_norma,ak,alfa

      enddo


	deallocate (p,pp,zz,z,r,rr) 
	
	end subroutine cg
      
   double precision function funnorm(n,y)
   implicit none
   double precision:: y(*)
   integer:: n
   !local
   integer :: isamax,i
   double precision:: sum2
   
        !isamax = 1
        !do  i=1,n
        !  if(abs(y(i)).gt.abs(y(isamax))) isamax=i
        !enddo


	  !funnorm=abs(y(isamax))

	  sum2=0.0

!$omp parallel do 
!$omp& shared(y), private(i), reduction(+: sum2) 
        do  i=1,n
          sum2=sum2 + y(i)*y(i)
        end do
!$omp end parallel do

	  funnorm=sqrt(sum2)


	return
	end

      

subroutine lin_bcg(np,ia,ja,a_spa,ad,b,x,unit_cont)
use def_constantes
integer :: np,unit_cont
double precision :: a_spa(*),ad(np),x(np),b(np)
integer :: ia(np+1),ja(*)
!local
double precision ::bnrm, znrm, bknum, bkden, xant, eps,akden,bk,ak,err,zm1nrm,dxnrm,xnrm,funnorm
double precision, allocatable ::z(:),p(:),pp(:),r(:),zz(:),rr(:)



allocate (p(np),pp(np),zz(np),z(np),r(np),rr(np)) 

iter=0
eps= 1.d-14
err = 1.0
      
call matxvecsim(ia,ja,a_spa,ad,x,r,np)

	do kk=1,np
      r(kk)  = b(kk) - r(kk)
	  rr(kk) = r(kk)   
	enddo

      call matxvecsim(ia,ja,a_spa,ad,r,rr,np)


 	do kk=1,np
	  z(kk)=b(kk)/ad(kk)
      if(abs(z(kk))>1.0e-9) then
         write(6,*) kk,zz(kk),b(kk)/ad(kk)
      endif
    enddo

      bnrm=funnorm(np,z)

	do kk=1,np
	  z(kk)=r(kk)/ad(kk)
      enddo

      znrm=funnorm(np,z)

      do while(err.gt.toler .and. iter.lt.itmax)

        iter=iter+1
      
 	  do kk=1,np
	    zz(kk)=rr(kk)/ad(kk)
        enddo
   	  
	  bknum=0.d0
  	  do j=1,np
          bknum=bknum + z(j)*rr(j)
        enddo

        if(iter.eq.1) then
          do  j=1,np
             p (j)=z(j)
             pp(j)=zz(j)
          enddo
        else
          bk=bknum/bkden
          do  j=1,np
            p(j) = bk*p(j)+z(j)
            pp(j)= bk*pp(j)+zz(j)
          enddo
        endif

        bkden=bknum

	  call matxvecsim(ia,ja,a_spa,ad,p,z,np)
      
    	  akden=0.d0
        do  j=1,np
          akden = akden+ z(j)*pp(j)
        enddo

        ak=bknum/akden
          
        call matxvecsim(ia,ja,a_spa,ad,pp,zz,np)

        do j=1,np
          x(j) = x(j)  + ak*p(j)
          r(j) = r(j)  - ak*z(j)
          rr(j)= rr(j) - ak*zz(j)
        enddo

     	  do kk=1,np
	    z(kk)=r(kk)/ad(kk)
        enddo

        zm1nrm=znrm
        znrm=funnorm(np,z)
        
	  if(abs(zm1nrm-znrm).gt.eps*znrm) then
          
		 dxnrm=abs(ak)*funnorm(np,p)
         err=znrm/abs(zm1nrm-znrm)*dxnrm
                  
      	 xnrm=funnorm(np,x)
          
	     if(err.le.0.5d0*xnrm) then
             err=err/xnrm
         else
             err=znrm/bnrm
         endif

	  else
        
	    err=znrm/bnrm
       
	  endif

    enddo

    write(unit_cont,*) 'iteraciones internas ', iter,err  

deallocate (p,pp,zz,z,r,rr) 
	
end subroutine lin_bcg

!double precision function funnorm(n,y) 
!integer :: n
!double precision :: y(n)
! local
!integer :: isamax,i

!  isamax=1
!  do  i=1,n
!      if(abs(y(i)).gt.abs(y(isamax))) isamax=i
!  enddo
        
!  funnorm=abs(y(isamax))

!end function funnorm
	
subroutine matxvecsim(ia,ja,an,ad,b,c,np)
implicit none
integer :: np
integer :: ia(np+1),ja(*)
double precision :: b(np),an(*),ad(np),c(np)
! local
integer :: k,i,iaf,iai,j

do k=1,np
     c(k)= ad(k)*b(k)
enddo
      
do i=1,np
     iai = ia(i)
     iaf = ia(i+1)- 1


     if(iaf.ge.iai) then
          do j=iai,iaf
            c(i) = c(i) + an(j)*b(ja(j))
            c(ja(j))=   c(ja(j)) + an(j)*b(i)
          enddo
      endif
enddo
     
end subroutine matxvecsim
      

