module writer_module
    use kind_module
    implicit none
    
contains

subroutine binary_write_array(v, a, time, array_name)
    !! сохраняет массивы расределения и скорости
    implicit none
    real(wp), intent(in) :: v(:,:)    
    real(wp), intent(in) :: a(:,:,:)
    real(wp), intent(in) :: time
    character(len=*), intent(in) :: array_name
    real(wp), allocatable :: gv(:), ga(:)
    integer i, N, nshape(3)
    character(120) fname
    integer, parameter :: iu = 21

    !if (MOD(INT(time*100000), 10) /= 0) then 
    !    return
    !end if
    
    nshape=shape(a) 
    print *, 'write arr:', array_name, nshape
    N = nshape(2)
    print *, N
    write(fname,'("lhcd/", A,"/", f9.7,".bin")') array_name, time
    print *, fname

    open(iu, file=fname ,status='new',action='write',access='stream',form='unformatted')
        write (iu),  nshape(1), nshape(2)
        write (iu), v
        write (iu), a
    close(iu)      

end subroutine

subroutine write_v_array(v, a, time, array_name)
    !! сохраняет массивы расределения и скорости
    implicit none
    real(wp), intent(in) :: v(:,:)    
    real(wp), intent(in) :: a(:,:,:)
    real(wp), intent(in) :: time
    character(len=*), intent(in) :: array_name
    real(wp), allocatable :: gv(:), ga(:)
    integer i, N, nshape(3)
    character(120) fname
    integer, parameter :: iu = 21

    !if (MOD(INT(time*100000), 10) /= 0) then 
    !    return
    !end if
    
    nshape=shape(a) 
    print *, 'write arr:', array_name, nshape
    N = nshape(2)
    print *, N
    write(fname,'("lhcd/", A,"/", f9.7,".dat")') array_name, time
    print *, fname
    open(iu, file=fname,position="append")
        do i=1, N
            gv = glue_v_axis(v(:,i))
            write (iu, '(20012(ES22.14))') gv(:)
            deallocate(gv)
            ga = glue_arrays(a(:,i,2), a(:,i,1)) ! negative - positive
            write (iu, '(20012(ES22.14))') ga(:)
            deallocate(ga)
        end do
    close(iu)

    contains
        function glue_v_axis(a) result(g)
            implicit none
            real(wp), intent(in) :: a(:)
            real(wp), allocatable :: g(:)
            integer i, N
            N = size(a)
            allocate(g(-N:N))
            g(-N:-1) = -a(N:1:-1)
            g(1:N) = a(:)
            g(0) = 0
        end function
        function glue_arrays(a, b) result(g)
            implicit none
            real(wp), intent(in) :: a(:), b(:)
            real(wp), allocatable :: g(:)
            integer i, N
            N = size(a)
            allocate(g(-N:N))
            g(-N:-1) = a(N:1:-1)
            g(1:N) = b(:)
            g(0) = (a(1) + b(1)) /2d0
        end function
end subroutine


subroutine write_x_array(x, arr, time, array_name)
    implicit none
    real(wp), intent(in) :: x(:,:)    
    real(wp), intent(in) :: arr(:,:)
    real(wp), intent(in) :: time
    character(len=*), intent(in) :: array_name
    
    integer i, N, nshape(2)
    character(120) fname
    integer, parameter :: iu = 21

    nshape=shape(arr) 
    print *, 'write arr:', array_name, nshape
    N = nshape(2)
    print *, N
    write(fname,'("lhcd/", A,"/xar", f9.7,".dat")') array_name, time
    print *, fname
    open(iu, file=fname,position="append")
        do i=1, N
            write (iu, '(20012(ES22.14))') x(: ,i)
            write (iu, '(20012(ES22.14))') arr(:, i)
        end do
    close(iu)
end subroutine

subroutine write_matrix(arr, time, array_name)
    implicit none
    real(wp), intent(in) :: arr(:,:)
    real(wp), intent(in) :: time
    character(len=*), intent(in) :: array_name
    
    integer i, N, nshape(2)
    character(120) fname
    integer, parameter :: iu = 21

    nshape=shape(arr) 
    print *, 'write_matrix:', array_name, nshape
    N = nshape(1)
    print *, N
    write(fname,'("lhcd/", A,"/", f9.7,".dat")') array_name, time
    print *, fname
    open(iu, file=fname,position="append")
        do i=1, N
            write (iu,' ( I4.4, 100(ES21.14))') i, arr(i, :)
        end do
    close(iu)
end subroutine

subroutine write_array(arr, N, array_name)
    implicit none
    real(wp), intent(in) :: arr(*)
    integer, intent(in) :: N
    character(len=*), intent(in) :: array_name
    
    integer i
    integer, parameter :: iunit = 21
    
    character(80) fname
    print *, 'write_array:', array_name, N
    write(fname,'("lhcd/distribution/", A,".dat")') array_name
    print *, fname
    open(iunit,file=fname,position="append")
        do i=1,n
            write(iunit,*) i, arr(i)
        end do
    close(iunit)
end subroutine

subroutine write_distribution(arr,N,time)
    implicit none
    real(wp), intent(in) :: arr(*)
    integer, intent(in) :: N
    real(wp), intent(in) :: time

    integer i
    integer itime
    integer, parameter :: iunit = 20
    character(120) fname

    itime = INT(time*100000)
    !print *, N, time, itime, MOD(itime, 10)
    if (MOD(itime, 10) == 0) then
        write(fname,'("lhcd/distribution/", f9.7,".dat")') time
        !print *, fname
        open(iunit,file=fname,position="append")
        do i=1, N
            if (arr(i)>0) then
                write(iunit,*) i, arr(i)
            else
                exit
            end if
        end do
        close(iunit)
    end if
end subroutine


end module writer_module