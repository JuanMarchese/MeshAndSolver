subroutine poisson()
use def_constantes
use def_solver
use def_variables
implicit none
double precision :: x(nodpel),y(nodpel),z(nodpel),esm(nodpel,nodpel),ef(nodpel),adiag,sigma_el,err,qe
!double precision, allocatable  :: esm_tot(:,:)

integer :: ns(nodpel)
integer  npaso,kk,jel,i,ii,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter
double precision, allocatable :: solucion_ant(:)
double precision :: error,epsil,denom,numer,sol(nodpel),funsigma1,campoxl
integer :: nconta,ncota


call control()

!allocate(solucion_ant(nnodes),esm_tot(nnodes,nnodes))
allocate(solucion_ant(nnodes))

!esm_tot=0.0

error=1.0
nconta=0
epsil = 1e-3
ncota=10
solucion_ant=0.0

do while(error>epsil .and. nconta< ncota )
    
    nconta=nconta+1
    an=0
    ad=0
    rhs=0
        
    do jel=1,nelements

        mat=material(jel)        
        if(mat==3 ) then
            sigma_el = sigmaext
        elseif(mat==2 ) then
            sigma_el = sigmamem
        elseif(mat==1) then
            sigma_el = sigmaint
        endif
   
        qe=0.0
        do i=1,nodpel
            ns(i)=conect(i,jel)
	        j=ns(i)
            x(i)=coor_x(j)
            y(i)=coor_y(j)

            sol(i)=solucion_ant(j)
        enddo
        
        gradxel_x(jel)=0.0
        gradxel_y(jel)=0.0


        call armado4(jel,x,y,ns,nodpel,esm,ef,sigma_el,qe,sol,gradxel_x(jel),gradxel_y(jel))

        do inode=1,nodpel
            ipoin=ns(inode)
            if(  vec_tierra(ipoin)/=-1 ) then
              adiag=esm(inode,inode)
              do jnode=1,nodpel
                 esm(inode,jnode)=0.0
                 ef (jnode)=ef(jnode)-esm(jnode,inode)*tierra
                 esm(jnode,inode)=0.0
              end do
              esm(inode,inode)= adiag
              ef(inode)       = adiag* tierra
            end if
            if(  vec_poten(ipoin)/=-1 ) then
              adiag=esm(inode,inode)
              do jnode=1,nodpel
                 esm(inode,jnode)=0.0
                 ef (jnode)=ef(jnode)-esm(jnode,inode)* potencial
                 esm(jnode,inode)=0.0
              end do
              esm(inode,inode)= adiag
              ef(inode)       = adiag* potencial
            end if
        end do


! ensamblo
        do ii=1,nodpel
            rhs(ns(ii))=rhs(ns(ii)) + ef(ii)
	        ad(ns(ii)) = ad(ns(ii)) + esm(ii,ii)
            	     
	        do jj=1,nodpel

                !esm_tot(ns(ii),ns(jj))=esm_tot(ns(ii),ns(jj)) + esm(ii,jj)

                jj2=ns(jj)+1-ns(ii)
	            if(jj2.gt.0) then
		 
	                do iaux = 1,cx(ns(ii))
         
		                keje = ia(ns(ii))+iaux-1  
	          	  	        
       		            if( ja(keje) .eq. ns(ii) + jj2 -1) then           
 		       
		                    an( keje ) = an(keje ) +  esm(ii,jj)
                
		                endif
          
		           enddo

	            endif
	
	        enddo
        enddo	  
    

    enddo

   
     
    solucion=rhs
    iter=0
    err=1.0
    call cg(nnodes,ia,ja,an,ad,rhs,solucion,toler,itermax,iter,err)
       do jj=1,nnodes
         !write(111,*) 'sol ',jj,solucion(jj)
        
         
        enddo

   error=epsil*0.5

    !write(unit_cont,*) 'iteraciones internas ', nconta, error, iter,err
    !write(6,*) 'iteraciones internas del loop poisson ', nconta, error, iter,err


enddo ! end while

call campo2d(nnodes,nelements,nodpel,solucion,material,conect,coor_x,coor_y,sigmaext,sigmaint,sigmamem,grad_x,grad_y,gradxel_x,gradxel_y)

!call salida_sol(solucion)


deallocate(solucion_ant)

end subroutine poisson	
      


double precision function funsigma1(mat,e)
implicit none
integer :: mat
double precision :: e

if(mat==4 .or. mat==2 .or. mat==3) then
   if(e>= 46.0 .and. e<=70.0) then
      funsigma1 = 1 + 2.5  
   else
      funsigma1 = 1
   endif
else
   funsigma1 = 1
endif

end function funsigma1

