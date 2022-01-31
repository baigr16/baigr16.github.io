    program prog
    real, allocatable :: e(:,:),eup(:,:),edn(:,:)
    real, allocatable :: k(:,:)
    real, dimension(3) ::k0,a
    real, dimension(6) ::xxxx
    character(len = 32):: xx, yy
    write(6,*) 'fermi level (eV)'
    read(5,*) ef
    open(10,file = 'EIGENVAL', status = 'old')
    open(11,file = 'band_fermi.dat')
    open(13,file = 'band.dat')
    read(10,*) iii, iii, iii, ispin  ! 第一行数据读取
    read(10,*) (xxxx(i),i = 1,5) ! 第二行数据读取
    read(10,*) xxxx(6) ! 第三行数据读取
    read(10,*) xx ! 第四行数据读取
    read(10,*) yy ! 第五行数据读取
    read(10,*) nn,nk,nbands  ! 第六行数据读取
    allocate(e(nk,nbands))
    allocate(k(nk,3))
    if(ispin.eq.2) then ! 2表示自旋极化计算
        do i = 1,nk
            read(10,*)
            read(10,*) (k(i,j),j = 1,3),wtk  ! 读取k点坐标及权重
            do j = 1,nbands  ! 读取确定k点下面的每条能带对应的本征值,这里分别是spin-up与spin-down
                read(10,*)jj,eup(i,j),edn(i,j)
            end do
            write(13,9030) (eup(i,j),j = 1,nbands)
            write(14,9030) (edn(i,j),j = 1,nbands)
        end do
        write(6,*)(k(i,j),j = 1,3)
        read(10,*) (e(i,n),n = 1,nbands)
        write(13,9030) (e(i,n),n = 1,nbands)
    else  ! 非自旋极化能带
        do i = 1,nk
            read(10,*)
            read(10,*) (k(i,j),j = 1,3),wtk
            do j = 1,nbands
                read(10,*) jj,e(i,j)  ! 读取能带k点和对应本征值
            end do
            write(13,9030) (e(i,j),j = 1,nbands)
        end do
    end if
9030 format (8f9.5)
    do j = 1,nbands
        dk = 0
        do i = 1,nk
            if (i.eq.1) then
                k0 = k(i,:)
            end if
            a = k(i,:) - k0
            dk = dk + sqrt(dot_product(a,a))
            if(ispin.eq.2) then
                write(11,*)dk,eup(i,j)-ef, edn(i,j)-ef
            else
                write(11,*)dk,e(i,j)-ef
            end if
            k0 = k(i,:)
        end do
        write(11,*)
    end do
    stop
    end program prog