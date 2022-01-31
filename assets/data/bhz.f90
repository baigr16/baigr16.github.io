    module pub
    implicit none
    integer yn,kn,ne,N,enn,hn
    real eta,dk,de
    parameter(yn = 50,hn =  4,kn = 50, ne = 50,N = yn*4,eta = 0.01,dk = 0.01,de = dk)
    real,parameter::pi = 3.1415926535
    complex,parameter::im = (0.,1.0)  
    complex Ham(N,N),one(N,N) 
!=================================
    real m0,mu   
    real tx,ty
    real ax,ay,gamma
    complex g1(hn,hn) 
    complex g2(hn,hn) 
    complex g3(hn,hn) 
!================cheev===============
    integer::lda = N
    integer,parameter::lwmax = 2*N + N**2
    real,allocatable::w(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   
    integer lrwork    
    integer liwork   
    integer info
    end module pub
!================= PROGRAM START ============================
    program sol
    use pub
    integer i1
!====空间申请==================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1+5*N+2*N**2))
    allocate(iwork(3+5*N))
!======parameter value setting =====
    m0 = 1.5
    mu = 0
    tx = 1.0
    ty = 1.0
    ax = 1.0
    ay = 1.0
    do i1 = 1,N
        one(i1,i1) = 1  !单位矩阵
    end do
    call main()
    stop 
    end program sol
!=============================================================================================
    subroutine main()
    !  Calculate edge spectrum function
    use pub
    integer m1,m2,i1
    real kx,ky,omega,re2
    complex h1(N,N),h2(N,N)
    complex re1
    open(30,file="openy-bhz.dat")
    !-------------------------------------------------
    !   y-direction is open
    do omega = -3.0,3.0,de
        do kx = -pi,pi,dk
            re1 = 0
            call openy(kx)
            h1 = omega*one - ham + im*eta*one
            call inv(h1,h2)
            do i1 = 1,N
                re1 = re1 + h2(i1,i1)
            end do
            re2 = -2*aimag(re1)/pi/N
            write(30,999)kx/pi,omega,re2
        end do
        write(30,*)" "
    end do
    close(30)
    !---------------------------------------------------------------------
    !  x-direction is open
    open(31,file="openx-bhz.dat")
    do omega = -3.0,3.0,de
        do ky = -pi,pi,dk
            re1 = 0
            call openx(ky)
            h1 = omega*one - ham + im*eta*one
            call inv(h1,h2)
            do i1 = 1,N
                re1 = re1 + h2(i1,i1)
            end do
            re2 = -2*aimag(re1)/pi/N
            write(31,999)ky/pi,omega,re2
        end do
        write(31,*)" "
    end do
    close(31)
999 format(3f11.5)
    return
    end subroutine main
!======================== Pauli Matrix driect product============================
    subroutine Pauli()
    use pub
    !   TI
!=== Kinetic energy
    g1(1,1) = 1
    g1(2,2) = -1
    g1(3,3) = 1
    g1(4,4) = -1
!====== SOC-x
    g2(1,2) = 1
    g2(2,1) = 1
    g2(3,4) = -1
    g2(4,3) = -1
!====== SOC-y
    g3(1,2) = -im
    g3(2,1) = im
    g3(3,4) = -im
    g3(4,3) = im
    return
    end subroutine Pauli
!==========================================================
    subroutine openx(ky)
    use pub
    real ky
    integer m,l,k
    call Pauli()
    Ham = 0
!========== Positive energy ========
    do k = 0,yn-1
        if (k == 0) then ! Only right block in first line
            do m = 1,hn
                do l = 1,hn
                    Ham(m,l) = (m0-ty*cos(ky))*g1(m,l) + ay*sin(ky)*g3(m,l) 

                    Ham(m,l + hn) = (-tx*g1(m,l) - im*ax*g2(m,l))/2
                end do
            end do
        elseif ( k==yn-1 ) then ! Only left block in last line
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) = (m0-ty*cos(ky))*g1(m,l) + ay*sin(ky)*g3(m,l)

                    Ham(k*hn + m,k*hn + l - hn) = -tx*g1(m,l)/2 + im*ax*g2(m,l)/2 
                end do
            end do
        else
            do m = 1,hn
                do l = 1,hn ! k start from 1,matrix block from 2th row
                    Ham(k*hn + m,k*hn + l) = (m0 - ty*cos(ky))*g1(m,l) + ay*sin(ky)*g3(m,l)

                    Ham(k*hn + m,k*hn + l + hn) = (-tx*g1(m,l) - im*ax*g2(m,l))/2 
                    Ham(k*hn + m,k*hn + l - hn) = -tx*g1(m,l)/2 + im*ax*g2(m,l)/2 
                end do
            end do
        end if
    end do
    !------------------------
    call isHermitian()
    ! call eigsol()
    return
    end subroutine openx
!==========================================================
    subroutine openy(kx)
    use pub
    real kx
    integer m,l,k
    call Pauli()
    Ham = 0
!========== Positive energy ========
    do k = 0,yn-1
        if (k == 0) then ! Only right block in first line
            do m = 1,hn
                do l = 1,hn
                    Ham(m,l) = (m0-tx*cos(kx))*g1(m,l) + ax*sin(kx)*g2(m,l)

                    Ham(m,l + hn) = (-ty*g1(m,l) - im*ay*g3(m,l))/2
                end do
            end do
        elseif ( k==yn-1 ) then ! Only left block in last line
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) = (m0-tx*cos(kx))*g1(m,l) + ax*sin(kx)*g2(m,l)

                    Ham(k*hn + m,k*hn + l - hn) = -ty*g1(m,l)/2 + im*ay*g3(m,l)/2 
                end do
            end do
        else
            do m = 1,hn
                do l = 1,hn ! k start from 1,matrix block from 2th row
                    Ham(k*hn + m,k*hn + l) = (m0-tx*cos(kx))*g1(m,l) + ax*sin(kx)*g2(m,l)

                    Ham(k*hn + m,k*hn + l + hn) = (-ty*g1(m,l) - im*ay*g3(m,l) )/2 
                    Ham(k*hn + m,k*hn + l - hn) = -ty*g1(m,l)/2 + im*ay*g3(m,l)/2 
                end do
            end do
        end if
    end do
    !---------------------------------
    call isHermitian()
    ! call eigsol()
    return
    end subroutine openy
!============================================================
    subroutine isHermitian()
    use pub
    integer i,j
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                open(160,file = 'hermitian.dat')
                write(160,*)i,j
                write(160,*)Ham(i,j)
                write(160,*)Ham(j,i)
                write(160,*)"===================="
                write(*,*)"Hamiltonian is not hermitian"
                stop
            end if
        end do
    end do
    close(160)
    return
    end subroutine isHermitian
!================= 矩阵本征值求解 ==============
    subroutine eigSol()
    use pub
    integer m
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork &
        ,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork &
        ,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(110,file="mes.dat",status="unknown")
        write(110,*)'The algorithm failed to compute eigenvalues.'
        close(110)
    end if
    ! open(120,file="eigval.dat")
    ! do m = 1,N
    !     write(120,*)w(m)
    ! end do
    ! close(120)
    return
    end subroutine eigSol
!=================================================================
    subroutine inv(matin,matout)
    use pub
    complex,intent(in) :: matin(N,N)
    complex:: matout(size(matin,1),size(matin,2))
    real:: work2(size(matin,1))            ! work array for LAPACK
    integer::info2,ipiv(size(matin,1))     ! pivot indices
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    matout = matin
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call CGETRF(N,N,matout,N,ipiv,info2)
    if (info2.ne.0)then
      write(*,*)'Matrix is numerically singular!'
      stop
    end if
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call CGETRI(N,matout,N,ipiv,work2,N,info2)
    if (info2.ne.0)then
        write(*,*)'Matrix inversion failed!'
        stop
    end if
    return
    end subroutine inv
