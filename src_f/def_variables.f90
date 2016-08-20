module def_variables

integer :: nelements,nnodes,nod_tierra,nod_poten,mat_externo,mat_membrana,mat_interno
integer :: num_mat
integer :: ntension=4
integer :: nod_mem_ext,nod_mem_int

double precision,allocatable :: coor_x(:),coor_y(:),grad_x(:),grad_y(:),cer(:),gradxel_x(:),gradxel_y(:)
double precision,allocatable :: tension(:,:,:),defor(:,:,:)

integer, allocatable :: conect(:,:), material(:),vec_tierra(:),vec_poten(:),vec_quieto(:)
integer, allocatable :: nodos_mem(:,:)


double precision :: sigmaint,sigmaext,sigmamem,permit
double precision :: potencial,tierra

integer :: nodpel=4   ! numero de elementos por nodo!!
integer :: nodpel2    ! para deforma.f90
integer :: ndimension

!constantes

double precision :: funyoung_in=0.0001
double precision :: funyoung_membrane=0.1
double precision :: funyoung_out=0.0001

double precision :: funposs_in=0.49
double precision :: funposs_membrane=0.49
double precision :: funposs_out=0.49

   
end module def_variables
