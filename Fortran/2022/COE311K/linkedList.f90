! Aaron Pandian (anp3238) Program LinkedList
implicit none

  type node
     integer :: value
     type(node),pointer :: next
  end type node

  type list
     type(node),pointer :: head
  end type list

  integer,parameter :: listsize=7
  type(list) :: the_list
  integer,dimension(listsize) :: inputs = &
       [ 62, 75, 51, 12, 14, 15, 16 ]
  integer :: input,input_value

  nullify(the_list%head)
  do input=1,listsize
     input_value = inputs(input)
     !call attach(the_list,input_value)
     call insert(the_list,input_value)
     !print *,"|list| = ",length(the_list)
     call print(the_list)
  end do

contains

  subroutine attach( the_list,new_value )
    implicit none
    ! parameters
    type(list),intent(inout) :: the_list
    integer,intent(in) :: new_value
    type(node),pointer :: current_node

    ! local

    ! if the list has no head node, attached the new node
    if (.not.associated(the_list%head)) then
       allocate( the_list%head )
       the_list%head%value = new_value
    else
       !call node_attach( the_list%head,new_value )
       current_node => the_list%head
       do while( associated(current_node%next) )
          current_node => current_node%next
       end do
       allocate(current_node%next)
       current_node%next%value = new_value
    end if

  end subroutine attach

  integer function length( the_list )
    implicit none
    type(list),intent(in),target :: the_list
    ! local
    type(node),pointer :: current

    length = 0
    current => the_list%head
    if ( associated(current) ) then
       do while( associated(current%next) )
          length = length+1
          !print *,current%value
          current => current%next
       end do
    end if

  end function length

  subroutine print(the_list)
    implicit none
    type(list),intent(in) :: the_list
    type(node),pointer :: current

    write(*,'("List: [ ")',advance="no")
    if (associated(the_list%head)) then
       current => the_list%head
       do while (associated(current))
          write(*,'(i0",")',advance="no") current%value
          if (.not.associated(current%next)) exit
          current => current%next
       end do
    end if
    write(*,'(x"]")')

  end subroutine print

 subroutine insert( the_list,new_value )
    implicit none
    ! parameters
    type(list),intent(inout) :: the_list
    integer,intent(in) :: new_value
    type(node),pointer :: current_node, temp_node
    integer :: logistic
    logistic = 0

    ! if the list has no head node, attached the new node
    if (.not.associated(the_list%head)) then
       allocate( the_list%head )
       the_list%head%value = new_value
    else if ( new_value .lt. the_list%head%value ) then
       allocate(temp_node)
       temp_node%next => the_list%head
       temp_node%value = new_value
       the_list%head => temp_node
    else
       ! call node_attach( the_list%head,new_value )
       current_node => the_list%head
       ! while loop to sort
       do while( associated(current_node%next) )
          if (new_value .lt. current_node%next%value) then
             allocate(temp_node)
             temp_node%value = new_value
             temp_node%next => current_node%next
             current_node%next => temp_node
             logistic = 1
             exit
          end if
          current_node => current_node%next
       end do
       if (logistic .eq. 0) then
          allocate(current_node%next)
          current_node%next%value = new_value
       end if
    end if
end subroutine insert

end Program
