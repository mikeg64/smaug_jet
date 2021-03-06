program main
    integer :: n        ! size of the vector
    real,dimension(:),allocatable :: a  ! the vector
    real,dimension(:),allocatable :: r  ! the results
    real,dimension(:),allocatable :: e  ! expected results
    integer :: i
    character(10) :: arg1
    if( iargc() .gt. 0 )then
        call getarg( 1, arg1 )
        read(arg1,'(i10)') n
    else
        n = 100000
    endif
    if( n .le. 0 ) n = 100000
    allocate(a(n))
    allocate(r(n))
    allocate(e(n))
    do i = 1,n
        a(i) = i*2.0
    enddo
    !$acc region
        do i = 1,n
            r(i) = a(i) * 2.0
        enddo
    !$acc end region
        do i = 1,n
            e(i) = a(i) * 2.0
        enddo
    ! check the results
    do i = 1,n
        if( r(i) .ne. e(i) )then
            print *, i, r(i), e(i)
            stop 'error found'
        endif
    enddo
    print *, n, 'iterations completed'
end program
