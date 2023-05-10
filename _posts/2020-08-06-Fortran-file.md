---
title: Fortran中文件批量操作方法
tags: Study Fortran
layout: article
license: true
toc: true
pageview: true
key: a20200806
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
Fortran属于比较古老的语言，自然不如现在大火的python等语言那么灵活，但是fortran的计算速度一直是其的优势，有时候再使用的过程中又会遇到要对不同的数据分开输出，这个时候就要利用Fortran的文件批量处理了，不过也是通过数据类型的转换，从而实现的。
{:.success}
<!--more-->
# 代码展示
```fortran
    program main
    implicit none
    integer m1
    do m1 = 1,100
        call eigsol(m1)
    end do
    stop
    end program main
!==================================================
    subroutine eigSol(input)
    integer m,input
    character*20::str1,str2,str3,str4
    str1 = "file"
    str3 = ".dat"
    write(str2,"(I4.4)")input
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(input,file=str4)
    write(input,*)input
    close(input)
    return
    end subroutine eigSol
```

再fortran中，通过对一个字符串进行数据写入**write(str2,"(I4.4)")input**，这样就可以将整型变量input转变成长度为4的字符串(1--->0001)，通过这样的方式就可以通过循环来实现批量文献读写。
![png](/assets/images/Fortran/fortran-file.png)

在这里，同时将input当作了打开文件的标示号，其实这里也可以只用确定的标示号来进行同样的操作。利用循环变量做标示号可能会遇到这个文件标示号你在这个循环中已经用过了，那么再次使用可能会报错或者将之前的文件内容修改。说这么多都是虚的，上代码演示一下。

```fortran
    program main
    implicit none
    integer m1
    open(30,file="val.dat")
    do m1 = 1,100
        call eigsol(m1)
    end do
    write(30,*)123
    close(30)
    stop
    end program main
!==================================================
    subroutine eigSol(input)
    integer m,input
    character*20::str1,str2,str3,str4
    str1 = "te"
    str3 = ".dat"
    write(str2,"(I4.4)")input
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(input,file=str4)
    write(input,*)input
    close(input)
    return
    end subroutine eigSol

```

![png](/assets/images/Fortran/fortran-file2.png)
这里可以看到，因为在循环之前先打开了30这个文件，但是还没有向30中写入数据，之后就开始了循环操作，由于循环操作中会重新打开以此30号文件，完成读写后关闭，这个时候再执行**write(30,*)123**就会产生**fort.30**这个文件，它里面的内容是123，而val.dat中则什么都没有。

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)