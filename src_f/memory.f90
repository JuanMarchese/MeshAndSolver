subroutine allocate_memory(number_nodes,nodes_per_elem,number_elements)
use def_solver
use def_variables
implicit none
integer :: number_nodes,nodes_per_elem,number_elements


    nodpel = nodes_per_elem
    nodpel2 = 2 * nodes_per_elem
    nelements = number_elements
    nnodes = number_nodes

    allocate(coor_x(nnodes))
    allocate(coor_y(nnodes))
    allocate(conect(nodpel,nelements))
    allocate(material(nelements))
    allocate(vec_tierra(nnodes))
    allocate(vec_poten(nnodes))
    allocate(vec_quieto(2*nnodes))
    allocate(grad_x(nnodes))
    allocate(grad_y(nnodes))
    allocate(gradxel_x(nelements))
    allocate(gradxel_y(nelements))

    allocate(tension(nelements,nodpel,ntension))
    allocate(defor(nelements,nodpel,ntension))


end subroutine allocate_memory




subroutine de_allocate_memory()
use def_solver
use def_variables
implicit none

    deallocate(coor_x)
    deallocate(coor_y)
    deallocate(conect)
    deallocate(material)
    deallocate(vec_tierra)
    deallocate(vec_poten)
    deallocate(vec_quieto)
    deallocate(grad_x)
    deallocate(grad_y)
    deallocate(gradxel_x)
    deallocate(gradxel_y)

    deallocate(tension)
    deallocate(defor)

    deallocate(solucion)

end subroutine de_allocate_memory

subroutine de_allocate_external_memory()
use def_solver
use def_variables
implicit none

    deallocate(ad)
    deallocate(ia)
    deallocate(cx)
    deallocate(rhs)
    deallocate(ad2d)
    deallocate(ia2d)
    deallocate(cx2d)
    deallocate(solucion2d)
    deallocate(rhs2d)
    deallocate(an)
    deallocate(ja)
    deallocate(an2d)
    deallocate(ja2d)

end subroutine de_allocate_external_memory


subroutine print_memory()
use def_solver
use def_variables
implicit none
integer :: i,j

    write(*,*) 'Coordenadas'
    do i=1,15
        write(*,*) i,' - ',coor_x(i),' ',coor_y(i),' ',solucion(i),' ',solucion2d(2*i-1),' ',solucion2d(2*i)
    enddo

    !write(*,*) 'Elementos'
    !do j=1,nelements
    !    write(*,*) j,' - ',conect(1,j),';',conect(2,j),';',conect(3,j),';',conect(4,j),' -- ',material(j)
    !enddo


end subroutine print_memory

