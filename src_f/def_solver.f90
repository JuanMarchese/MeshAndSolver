module def_solver


  double precision, allocatable :: ad(:),an(:),rhs(:),solucion(:)
  integer,allocatable :: cx(:),ia(:),ja(:)

  integer :: nonull
  
  ! caso tensiones
  double precision, allocatable :: ad2d(:),an2d(:),rhs2d(:),solucion2d(:)
  integer,allocatable :: cx2d(:),ia2d(:),ja2d(:)

  integer :: nonull2d
  
    
end module def_solver