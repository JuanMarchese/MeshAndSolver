subroutine deforma
use def_variables
use def_solver
use def_constantes
implicit none
double precision :: x(nodpel),y(nodpel),esm(nodpel2,nodpel2),ef(nodpel2),adiag,sigma_el,err,sigma0(ntension),cte0
double precision :: young,poiss,funyoung,funposs
double precision :: ecte,ecuad   !!! OJO
integer :: ns2d(nodpel2)
integer  NPASO,KK,JEL,I,II,ncase,inode,ipoin,jnode,jj,iaux,keje,jj2,mat,j,iter


!cte0=299.79*299.79
ecte = 8.8419E-7   !!! OJO
an2d=0
ad2d=0
rhs2d=0
        
DO JEL=1,nelements

    mat=material(jel)        
    
       
    young= funyoung(mat)
    poiss= funposs(mat)

    DO I=1,nodpel
	   j=conect(JEL,I)
       ns2d(2*I-1)=2*j-1
       ns2d(2*I)=2*j
       X(i)=coor_x(j)
       y(i)=coor_y(j)

   ENDDO
   sigma0=0.0  ! tensor de maxwell
   if (mat==2) then
   !!  
       ecuad = (gradxel_x(jel)*gradxel_x(jel) + gradxel_y(jel)*gradxel_y(jel))   !!! OJO
       
       sigma0(1) = ecte * (gradxel_x(jel)*gradxel_x(jel) - ecuad*0.5) ! gradxel_x(jel)*gradxel_x(jel)*0.5/cte0   !!! OJO
       sigma0(2) = ecte * (gradxel_y(jel)*gradxel_y(jel) - ecuad*0.5) ! gradxel_y(jel)*gradxel_y(jel)*0.5/cte0   !!! OJO
       sigma0(3) = ecte * (gradxel_x(jel)*gradxel_y(jel))             ! gradxel_x(jel)*gradxel_y(jel)/cte0       !!! OJO
       sigma0(4) = 0.0
   endif

   CALL ARMADO4_elas(JEL,X,Y,ns2d,nodpel2,ESM,EF,young,poiss,sigma0)
! INTRODUZCO EN LAS MATRICES LAS CONDICIONES DE CONTORNO

        do inode=1,nodpel2
            ipoin=NS2d(inode)
            if(  vec_quieto(ipoin)/=-1 ) then
              adiag=ESM(inode,inode)
              do jnode=1,nodpel2
                 ESM(inode,jnode)=0.0
                 EF (jnode)=EF(jnode)-ESM(jnode,inode)*0.0
                 ESM(jnode,inode)=0.0
              end do
              ESM(inode,inode)= adiag
              EF(inode)       = adiag* 0.0
            end if
        end do


! ENSAMBLO
        DO II=1,nodpel2
            RHS2d(NS2d(II))=RHS2d(NS2d(II)) + EF(II)
	        AD2d(NS2d(II)) = AD2d(NS2d(II)) + ESM(II,II)
            	     
	        DO JJ=1,nodpel2


                JJ2=NS2d(JJ)+1-NS2d(II)
	            IF(JJ2.GT.0) THEN
		 
	                DO IAUX = 1,CX2d(NS2d(II))
         
		                KEJE = IA2d(NS2d(II))+IAUX-1  
	          	  	        
       		            IF( JA2d(KEJE) .EQ. NS2d(II) + JJ2 -1) THEN           
 		       
		                    AN2d( KEJE ) = AN2d(KEJE ) +  ESM(II,JJ)
                
		                ENDIF
          
		           ENDDO

	            ENDIF
	
	        ENDDO
        ENDDO	  
    

    ENDDO

     
    solucion2d=RHS2d
    iter=0
    err=1.0
    call CG(nnodes*2,IA2d,JA2d,AN2d,AD2d,RHS2d,solucion2d,toler,itermax*10,ITER,ERR)

    call cal_tension(nnodes*2,solucion2d,solucion)

end subroutine deforma



subroutine cal_tension(np,solucion2d,solucion) 
use def_variables
implicit none
integer  :: np
double precision :: solucion2d(np),solucion(nnodes)
!local
double precision :: B(ntension,nodpel2),D(ntension,ntension),aux1(nodpel2,ntension)
double precision :: em,pr,funem,funpos,R
double precision :: X(nodpel),Y(nodpel),U(2*nodpel)
DOUBLE PRECISION,allocatable :: gauspt(:),gauswt(:),phi(:,:),dphi(:,:,:),gxcod(:,:),PHIdX(:,:,:),cteI(:) 
DOUBLE PREcIsION:: PI=3.14159
integer ::i,j,Ngaus,kgaus,jj,ii,kk,mat,ndime,k,ndi,nope,ns(nodpel2)
	 DOUBLE PREcIsION:: t,s,sm,tm,sq,tp,AJACO2d(2,2),AJACOI2d(2,2),deter,funyoung,funposs


ntension=4
Ngaus =2
ndimension=2
nope=nodpel2*0.5

allocate(gauspt(Ngaus),gauswt(Ngaus),phi(Ngaus*Ngaus,nope),dphi(ndimension,Ngaus*Ngaus,nope),gxcod(ndimension,Ngaus*Ngaus),PHIdX(ndimension,Ngaus*Ngaus,nope),cteI(Ngaus*Ngaus)) 

  GAUSPT(1)= -0.57735027
  GAUSPT(2)= 0.57735027 
  GAUSWT(1)=1.0
  GAUSWT(2)=1.0



open(unit=111,file='malla_final.dat')

write(111,* )'nodos membrana externa', nod_mem_ext+nod_mem_int
do kk=1,nod_mem_ext
   jj=nodos_mem(2,kk)
   write(111,'(i4,5E15.5)') jj,coor_x(jj),coor_x(jj) + solucion2d(2*jj-1),coor_y(jj),coor_y(jj) + solucion2d(2*jj),solucion(jj)
enddo
do kk=1,nod_mem_int
   jj=nodos_mem(1,kk)
   write(111,'(i4,5E15.5)') jj,coor_x(jj),coor_x(jj) + solucion2d(2*jj-1),coor_y(jj),coor_y(jj) + solucion2d(2*jj),solucion(jj)
enddo


write(111,* )'    '
write(111,* )'    '
write(111,* )'    '
DO I=1,nnodes
    write(111,'(i4,4E15.5)') i,coor_x(i),solucion2d(2*i-1),coor_y(i), solucion2d(2*i)
    coor_x(i)=coor_x(i) + solucion2d(2*i-1)
    coor_y(i)=coor_y(i) + solucion2d(2*i)
ENDDO



do kk=1,nelements
  
  mat = material(kk)

  EM = funyoung(mat)
  pr = funposs(mat)
  
    
  do i=1,nope
        
       j = conect(kk,i)
       ns(2*i-1)= 2*j-1
       ns(2*i  )= 2*j
       x(i)=coor_x(j)
       y(i)=coor_y(j)
       u(2*i-1) = solucion2d(ns(2*i-1))
       u(2*i  ) = solucion2d(ns(2*i  ))

  enddo
 
 
 R=EM/(1.+PR)

      
  D(1,1)=R*(1.-PR)/(1-2*PR)
  D(2,2)=R*(1.-PR)/(1-2*PR)
  D(4,4)=R*(1.-PR)/(1-2*PR) 
      
  D(3,3)=R/2.
  D(1,2)=PR*R/(1.-2*PR)
  D(2,1)=D(1,2)
    
  D(1,4)=PR*R/(1.-2*PR) 
  D(4,1)=D(1,4)
  D(2,4)=D(1,4)
  D(4,2)=D(1,4)
  D(1,3)=0.0
  D(2,3)=0.0
  D(4,3)=0.0
  D(3,1)=0.0
  D(3,2)=0.0
  D(3,4)=0.0

   kgaus=0 
        do jj=1,Ngaus
          do ii=1,Ngaus
              kgaus=kgaus+1
              t=GAUSPT(jj)
              s=GAUSPT(ii)
             
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
              
               do ndime=1,ndimension
                 gxcod(ndime,kgaus)=0.0
                 do i=1,nope
                     gxcod(ndime,kgaus) = gxcod(ndime,kgaus) + x(i)*phi(kgaus,i)
                 enddo
               enddo


               AJACO2d=0.0
               AJACOI2d=0.0
               DO K=1,nope
                    AJACO2d(1,1)=AJACO2d(1,1)+dphi(1,kgaus,K)*X(K)
                    AJACO2d(1,2)=AJACO2d(1,2)+dphi(1,kgaus,K)*Y(K)
                    AJACO2d(2,1)=AJACO2d(2,1)+dphi(2,kgaus,K)*X(K)
                    AJACO2d(2,2)=AJACO2d(2,2)+dphi(2,kgaus,K)*Y(K)
               ENDDO

               DETER= AJACO2d(1,1)*AJACO2d(2,2)-AJACO2d(1,2)*AJACO2d(2,1)
               
              
               AJACOI2d(1,1)=AJACO2d(2,2)/DETER
               AJACOI2d(2,2)=AJACO2d(1,1)/DETER
               AJACOI2d(1,2)=-AJACO2d(1,2)/DETER
               AJACOI2d(2,1)=-AJACO2d(2,1)/DETER

 
!   DERIVADAS DE LAS FUNCIONES DE INTERPOLACION EN FUNCION DE X-Y-Z

               DO i=1,nope
                    PHIdX(1,kgaus,i)=AJACOI2d(1,1)*dphi(1,kgaus,i) + AJACOI2d(1,2)*dphi(2,kgaus,i) 
                    PHIdx(2,kgaus,i)=AJACOI2d(2,1)*dphi(1,kgaus,i) + AJACOI2d(2,2)*dphi(2,kgaus,i) 
               ENDDO

               cteI(kgaus)= deter * GAUSWT(ii)* GAUSWT(jj)*2*PI*GXCOD(1,kgaus) 

           enddo
        enddo 


  DO Kgaus=1,Ngaus*Ngaus

        B=0.0
     
        
        DO i=1,nope
               B(1,2*i-1) = PHIdX(1,kgaus,i)
               B(2,2*i)   = PHIdx(2,kgaus,i)
               B(3,2*i-1) = B(2,2*i)
               B(3,2*i)   = B(1,2*i-1)
               B(4,2*i-1) = PHI(kgaus,i)/GXCOD(1,kgaus) 
        ENDDO   
      
      
    
       do k=1,ntension
           defor(kk,kgaus,k) = 0.0
           do i=1,ndimension*nope                             
                defor(kk,kgaus,k) = defor(kk,kgaus,k)+ B(K,I)*U(i)
           enddo
	   enddo

       DO k=1,ntension 
          tension(kk,kgaus,k)=0.0
           DO i=1,ntension
             tension(kk,kgaus,k)=tension(kk,kgaus,k) + D(k,i)*defor(KK,Kgaus,i)
             enddo
       enddo

     
  enddo     ! fin loop ngaus
        

enddo ! fin   loop elementos 
   


end subroutine cal_tension


double precision function funyoung(mat)
use def_variables
implicit none
integer mat

if(mat==1) then
   funyoung=funyoung_in
elseif(mat==2) then
   funyoung=funyoung_membrane  ! dyn/micron**2
elseif(mat==3) then
   funyoung=funyoung_out
endif

end function funyoung


double precision function funposs(mat)
use def_variables
implicit none
integer mat

if(mat==1) then
   funposs=funposs_in
elseif(mat==2) then
    funposs=funposs_membrane
elseif(mat==3) then
   funposs=funposs_out
endif

end function funposs

