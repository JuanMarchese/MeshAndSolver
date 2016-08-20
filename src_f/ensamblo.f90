!   calcula las matrices de cada elemento
SUBROUTINE ARMADO4(NLE,X,Y,ns,NOPE,ESM,EF,sigma_el,qe,sol,Ex_el,Ey_el)
implicit none
INTEGER :: nope,NS(NOPE),NLE
DOUBLE PRECISION :: X(NOPE),Y(NOPE),EF(NOPE),ESM(NOPE,NOPE),sigma_el,qe,sol(nope),Ex_el,Ey_el
DOUBLE PRECISION,allocatable :: gauspt(:),gauswt(:),phi(:,:),dphi(:,:,:),gxcod(:,:),PHIdX(:,:,:),cteI(:)
DOUBLE PREcIsION:: PI=3.14159
integer ::i,j,Ngaus,kgaus,jj,ii,ndimension,ndime,k,ndi
     DOUBLE PREcIsION:: t,s,sm,tm,sq,tp,AJACO2d(2,2),AJACOI2d(2,2),deter

 Ngaus =2
 ndimension=2

 esm=0.0
 ef=0.0

 allocate(gauspt(Ngaus),gauswt(Ngaus),phi(2*Ngaus,nope),dphi(ndimension,2*Ngaus,nope),gxcod(ndimension,2*Ngaus),PHIdX(ndimension,2*Ngaus,nope),cteI(Ngaus*2))

        GAUSPT(1)= -0.57735027
        GAUSPT(2)= 0.57735027
        GAUSWT(1)=1.0
        GAUSWT(2)=1.0

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


     do kgaus=1,Ngaus*Ngaus

       DO I=1,NOPE
          DO J=1,NOPE
             do ndi=1,ndimension
                ESM(I,J)=ESM(I,J)+ sigma_el*PHIdX(ndi,kgaus,I)*PHIdX(ndi,kgaus,J)*cteI(kgaus)
             enddo
         ENDDO
         EF(I)=EF(I) + cteI(kgaus)*PHI(kgaus,I)*QE
       ENDDO

     enddo


end subroutine armado4

subroutine ARMADO4_elas(NLE,X,Y,ns2d,nodpel2,ESM,EF,EM,pr,sigma0)

implicit none
integer :: NLE,nodpel2,ns2d(nodpel2)
double precision :: em,pr,r,sigma0(nodpel2),x(nodpel2/2),y(nodpel2/2),esm(nodpel2,nodpel2),ef(nodpel2)
!local
double precision :: B(4,nodpel2),D(4,4),aux1(nodpel2,4)

DOUBLE PRECISION,allocatable :: gauspt(:),gauswt(:),phi(:,:),dphi(:,:,:),gxcod(:,:),PHIdX(:,:,:),cteI(:)
DOUBLE PREcIsION:: PI=3.14159
integer ::i,j,Ngaus,kgaus,jj,ii,ndimension,ndime,k,ndi,ntension,nope
     DOUBLE PREcIsION:: t,s,sm,tm,sq,tp,AJACO2d(2,2),AJACOI2d(2,2),deter


ESM=0.0
EF=0.0
D=0.0

ntension=4
Ngaus =2
ndimension=2
nope=nodpel2*0.5

allocate(gauspt(Ngaus),gauswt(Ngaus),phi(2*Ngaus,nope),dphi(ndimension,2*Ngaus,nope),gxcod(ndimension,2*Ngaus),PHIdX(ndimension,2*Ngaus,nope),cteI(Ngaus*2))

  GAUSPT(1)= -0.57735027
  GAUSPT(2)= 0.57735027
  GAUSWT(1)=1.0
  GAUSWT(2)=1.0



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


     do kgaus=1,Ngaus*Ngaus

          B=0.0

         DO K=1,nope
               B(1,2*K-1) = PHIdX(1,kgaus,K)
               B(2,2*K)   = PHIdx(2,kgaus,K)
               B(3,2*K-1) = B(2,2*K)
               B(3,2*K)   = B(1,2*K-1)
               B(4,2*K-1) = PHI(kgaus,K)/GXCOD(1,kgaus)
         ENDDO


      do k=1,ndimension*NOPE
        do j=1,ntension
          aux1(k,j)=0.0
          do i=1,ntension
             aux1(k,j) = aux1(k,j) + B(i,k)*D(i,j)
          enddo
        enddo
      enddo


      DO I=1,ndimension*NOPE
         DO J=1,ndimension*NOPE
           do k=1,ntension
              ESM(I,J)=ESM(I,J)+ aux1(I,k)*b(k,J)*cteI(kgaus)
           enddo
         enddo
         DO j=1,ntension
           EF(I)=EF(I) -  B(j,i)*sigma0(j)*cteI(kgaus)
         ENDDO


      ENDDO


   enddo




end subroutine ARMADO4_elas

