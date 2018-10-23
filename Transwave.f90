!!$ thanks: http://www.andrew.cmu.edu/user/feenstra/wavetrans/
!!$ modified by: yang zhi long
!!$ mail: yangzhilong26@gmail.com
!!$ WAVECAR data structure 
!!$ #Record-length #spin components #RTAG(a value specifying the precision)
!!$ #k-points #bands ENCUT(maximum energy for plane waves)
!!$ LatVec-A
!!$ LatVec-B
!!$ LatVec-C
!!$ Loop over spin
!!$    Loop over k-points
!!$        plane waves number, k vector
!!$        Loop over bands
!!$            band energy, band occupation
!!$        End loop over bands
!!$        Loop over bands
!!$            Loop over plane waves
!!$                Plane-wave coefficient
!!$            End loop over plane waves
!!$        End loop over bands
!!$    End loop over k-points
!!$ End loop over spin
!!$ 
!!$ Plane waves  order see below

!!$   input the WAVECAR file in binary format from VASP, and write
!!$   selected real space wavefunction in a3 direction to standard output
!!$   Plane wave coefficients are written to GCOEFF.txt
!!$   Compile with gfortran or ifort. Flag "-assume byterecl" is required
!!$   for ifort.
!!$   version 2.0 - Jul 11, 2012 - R. M. Feenstra and M. Widom, 
!!$                                using updated 'c' value
!!$   version 2.1 - Sep 17, 2014 - changed estimator for max. no. of
!!$                                plane waves
!!$   version 2.2 - Nov 12, 2014 - output (n1,n2,n3) indices in both first half
!!$                                and second half of GCOEFF.txt output file

program  readwavecar
    implicit none
    !!$ 平面波展开系数
    !! complex*8 等价于complex(kind=4) 兼容性问题
    !! 等价于complex(8), allocatable :: coeff(:)
    complex(kind=4), allocatable :: coeff(:)
    !!$ cener 能量  复数形式,虚部为0
    !!等价于complex*16, allocatable :: cener(:)
    complex(8), allocatable :: cener(:)
    !!$ occ 占据态数目
    real(8), allocatable :: occ(:)
    !!$ k点处波函数 所有的平面波  n1G1+n2G2+n3G3
    integer, allocatable :: igall(:,:)
    !!$ 实空间晶格矢量ai，倒空间晶格矢量bi,sumkg
    !!$倒格式k+nG笛卡尔坐标形式,vtmp=b1xb2
    real(8) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3),vtmp(3)
    !!$ wk(3) k格点三个方向分量;xyz(3) 实空间坐标，以基矢a1，a2，a3为单位
    real(8)    :: wk(3),xyz(3),wkpg(3),ig(3)
    !!$ csum1,csum2 波函数 
    !!等价于 complex*16 ::csum1,csum2
    complex(8)::csum1,csum2
    integer ::kpoint,band,spin
    character(75) :: filename,soc
    real(8) :: c,pi,x,y,z
    integer :: iost
    real(8) :: xnrecl,xnspin,xnprec
    integer :: nrecl,nspin,nprec
    real(8) :: xnwk,xnband
    real(8)    :: ecut 
    integer :: nwk,nband
    integer :: npmax
    real(8) :: a1mag,a2mag,a3mag,b1mag,b2mag,b3mag
    real(8) :: Vcell,vmag,phi123,sinphi123
    real(8) :: phi12,phi13,phi23
    integer :: nb1maxA,nb2maxA,nb3maxA
    integer :: nb1maxB,nb2maxB,nb3maxB
    integer :: nb1maxC,nb2maxC,nb3maxC
    integer :: npmaxA,npmaxB,npmaxC
    integer :: nb1max,nb2max,nb3max
    integer :: i,j,iband,irec,iplane,iz
    real(8) :: xnplane
    integer :: nplane
    integer :: ncnt
    integer :: ig1,ig2,ig3,ig1p,ig2p,ig3p
    real(8) :: gtot,etot


    !!$*  物理量单位  ev  Ang
    !!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
    !!$*   adjusted in final decimal places to agree with VASP value; program
    !!$*   checks for discrepancy of any results between this and VASP values)

    data c/0.262465831d0/ 
    pi=4.*atan(1.)

    !!$ parse arguments
    !! 读取控制参数
    call parse(filename,soc,spin,kpoint,band,x,y)
    xyz(1)=x
    xyz(2)=y

    !!$*   input
    nrecl=24
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost            
    !!$ 文件记录单元长度，spin，精确度
    !!$模块单元的长度是由nprec来定义
    read(unit=10,rec=1) xnrecl,xnspin,xnprec
    close(unit=10)
    nrecl=nint(xnrecl)
    nspin=nint(xnspin)
    nprec=nint(xnprec)
    !!$ WAVECAR_double 变量是否定义
    !!$
    !!$该程序只支持WAVECAR_double没有定义的vasp版本。
    !!$如果定义了，则只需要修改coeff由8字节存储改为16字节
    if(nprec.eq.45210) then
       write(6,*) '*** error - WAVECAR_double requires complex*16'
       write(6,*) 'need to change the array coeff(:) from complex*8 to complex*16'
       stop
    endif
    !!$  该代码块为共线版本  ispin=1 
    !!if(nspin.eq.2) then
    !!   write(6,*) '*** error - Not a spinor WAVECAR. ISPIN =',nspin
    !!   stop
    !!endif

    !! nspin=2 必定是没有考虑soc
    !! nspin=1 可能是考虑soc或者没有考虑
    if (nspin .eq. 2) then 
        npmax=1
        if (soc .eq. 't') then 
            write(6,*) "input wrong!!! this case doesn't include soc"
            call help
            stop
        end if
    else if(nspin .eq. 1) then
        if(soc .eq. 't') then
            npmax=2
            if(spin .eq. 2) then
                write(6,*) 'current caseinclude soc, spin must be 1'
                stop
            end if
        else if(soc .eq. 'f') then
            npmax=1
        end if 
    end if 
    !!$ 输出基本信息
    write(6,*) ''
    write(6,*) 'record length  =',nrecl,' spins =',nspin ,'  prec flag ',nprec
    !!$ 重新打开WAVECAR文件，文件单元长度为nprec
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost
    !!$ 输出文件
    open(unit=11,file='GCOEFF.txt',position="append")
    !!$ 第二个模块单元
    !!$ xnw:k点数目,xnband:能带数目,ecut:截断能量
    read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
         (a3(j),j=1,3)
    nwk=nint(xnwk)
    nband=nint(xnband)
    !!$  所选的k点能带指标n不能超过范围,能带同理
    if (kpoint.gt.nwk) then
       write(0,*) '*** error - selected k=',kpoint,' > max k=',nwk
       stop
    endif
    if (band.gt.nband) then
       write(0,*) '*** error - selected band=',band,' > max band=',nband
       stop
    endif
    
    allocate(occ(nband))
    allocate(cener(nband))

    write(6,*)  'no. k points =',nwk
    write(6,*)  'no. bands =',nband
    write(6,*)  'real space lattice vectors:'
    write(6,*)  'a1 =',(sngl(a1(j)),j=1,3)
    write(6,*)  'a2 =',(sngl(a2(j)),j=1,3)
    write(6,*)  'a3 =',(sngl(a3(j)),j=1,3)
    write(6,*)  ' '
    write(6,*)  nspin
    write(6,*)  nwk
    write(6,*)  nband
    write(6,*)  ' '

    !!$*   compute reciprocal properties
    call vcross(a2xa3,a2,a3)
    Vcell=dot_product(a1,a2xa3)
    a3mag=dsqrt(dot_product(a3,a3))
    call vcross(b1,a2,a3)
    call vcross(b2,a3,a1)
    call vcross(b3,a1,a2)
       b1=2.*pi*b1/Vcell
       b2=2.*pi*b2/Vcell
       b3=2.*pi*b3/Vcell
    !!$ 倒格式 模长
    b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
    b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
    b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
    
    write(6,*) 'volume unit cell =',sngl(Vcell)
    write(6,*) 'reciprocal lattice vectors:'
    write(6,*) 'b1 =',(sngl(b1(j)),j=1,3)
    write(6,*) 'b2 =',(sngl(b2(j)),j=1,3)
    write(6,*) 'b3 =',(sngl(b3(j)),j=1,3)
    write(6,*) 'reciprocal lattice vector magnitudes:'
    write(6,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)
    write(6,*) ' '
    write(11,*) (sngl(a1(j)),j=1,3)
    write(11,*) (sngl(a2(j)),j=1,3)
    write(11,*) (sngl(a3(j)),j=1,3)
    write(11,*) (sngl(b1(j)),j=1,3)
    write(11,*) (sngl(b2(j)),j=1,3)
    write(11,*) (sngl(b3(j)),j=1,3)


    phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
    call vcross(vtmp,b1,b2)
    vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
    sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
    phi123=abs(asin(sinphi123))
    nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
    nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
    nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
    npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
          
    phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
    call vcross(vtmp,b1,b3)
    vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
    sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
    phi123=abs(asin(sinphi123))
    nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
    nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
    nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
    npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
          
    phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
    call vcross(vtmp,b2,b3)
    vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
    sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
    phi123=abs(asin(sinphi123))
    nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
    nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
    nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
    npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

    nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
    nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
    nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
    !! 2* to handle two component spinors 
    !! ispin=1,平面波分解成分中含有spin up 和spin  down 的成分
    !!写成通用形式，这里需要修改
    !! npmax 估计的最大的平面波数目
    npmax=npmax*min0(npmaxA,npmaxB,npmaxC)

    write(6,*) 'max. no. G values: 1,2,3 =',nb1max,nb2max,nb3max
    write(6,*) ' '


    allocate (igall(3,npmax))
    write(6,*) 'estimated max. no. plane waves =',npmax
    allocate (coeff(npmax))

    !!$*   Find desired wavefunction
    !!$     定位到查询的模块单元位置
    irec=3+(kpoint-1)*(nband+1)*(nspin)
    !! xnplane  k点处波函数展开的平面波数目
    !! 只读取关心的部分的信息 
    read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
         (cener(iband),occ(iband),iband=1,nband)
    !! test
    !!write(11,*) 'xnplane',xnplane
    nplane=nint(xnplane)

    write(11,*) "nplane",nplane
    write(11,*) kpoint,band
    write(11,*) (sngl(wk(j)),j=1,3)
    write(11,*) cener(band),occ(band)
          
    !!$*   Calculate plane waves
    !!!平面波矢量是[(-nb1max,nb1max),(-nb2max,nb2max),(-nb3max,nb3max)]的立方体内满足k+nG小于G_cut的圆内的点nG

    ncnt=0
    !!  第一个数先颠倒符号   注意顺序  ig3 ig2 ig1
    do ig3=0,2*nb3max
       ig3p=ig3
       ! 负的
       if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
       do ig2=0,2*nb2max
          ig2p=ig2
          if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
          do ig1=0,2*nb1max
             ig1p=ig1
             if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
             do j=1,3
                 !! k+n1b1+n2b2+n3b3 的笛卡尔坐标形式
                 !! h^2|k+G|^2/2m
                sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                     (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
             enddo
             !!k+n1b1+n2b2+n3b3 矢量长度
             gtot=sqrt(dot_product(sumkg,sumkg))
             !! 动能E,  限制了在第一布里渊区的k点
             etot=gtot**2/c
             !! test
             !write(11,*)  "ig1p",ig1p
             if (etot.lt.ecut) then
                ncnt=ncnt+1
                !! test
                !write(11,*) "ncnt",ncnt
                !! 所有的平面波k+nG的系数n1 n2 n3
                igall(1,ncnt)=ig1p
                igall(2,ncnt)=ig2p
                igall(3,ncnt)=ig3p
             end if
          enddo
       enddo
    enddo

    if (nspin .eq. 2) then
        if (ncnt .ne. nplane) then
            write(6,*) "*** error - computed ncnt=",ncnt,&
                    '  != input nplane=',nplane
            if(ncnt .eq. (2*nplane-1)) write(6,*) '*** suspect Gamma-only WAVECAR'
            stop
        endif
    else if (nspin .eq. 1) then
        if(soc .eq. 't') then
            if (2*ncnt .ne. nplane) then
               write(0,*) '*** error - computed 2*ncnt=',2*ncnt, &
                    ' != input nplane=',nplane
               stop
            endif
        else if (soc .eq. 'f') then
            if (ncnt .ne. nplane) then
                write(6,*) "*** error - computed ncnt=",ncnt,&
                        '  != input nplane=',nplane
                if(ncnt.eq.(2*nplane-1)) write(6,*) '*** suspect Gamma-only WAVECAR'
                stop
            endif
        endif
    endif



    !! 波矢k，能带n处平面波展开系数
    irec=irec+band
    write(6,*) 'irec',irec,'nplane',nplane
    read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)

    !!$*   output G value and coefficients
    !!$ G矢量的顺序
    do iplane=1,nplane
       write(11,570) (igall(j,mod(iplane-1,ncnt)+1),j=1,3), &
            coeff(iplane)
    570 format(3i6,'  ( ',g14.6,' , ',g14.6,' )')     
    enddo


    !! 实空间的波函数phi(x,y,z)
    do iz=0,2*nb3max
       z=dble(iz)/dble(1+2*nb3max)
       xyz(3)=z
       csum1=cmplx(0.,0.)
       csum2=cmplx(0.,0.)
       do iplane=1,ncnt
          ig=real(igall(:,iplane))
          wkpg=wk+ig
          !! 平面波  psi=sum C e^(ikr)
          csum1=csum1+coeff(iplane)* &
               cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
          if(nspin==1 .and. soc=='t') then
              csum2=csum2+coeff(ncnt+iplane)* &
                   cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
          endif
       enddo
       !! 除以体积，归一化
       csum1=csum1/dsqrt(Vcell)
       if((nspin.eq.1) .and. (soc.eq.'t')) csum2=csum2/dsqrt(Vcell)
    !!$ output z*a3mag for units of Angstroms
    !! phi(z)  实部，虚部
       if((nspin.eq.1) .and. (soc.eq.'t')) then
           write(6,*) sngl(z),sngl(real(csum1)),sngl(dimag(csum1)),&
                sngl(real(csum2)),sngl(dimag(csum2))
       else
            write(6,*) sngl(z),sngl(real(csum1)),sngl(dimag(csum1))
       endif 
    enddo

    deallocate(igall)
    deallocate(coeff)
    stop
end program readwavecar

!!$*   routine for computing vector cross-product
subroutine vcross(a,b,c)
  implicit none
  real(8):: a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross      
      
!!$   parse command line arguments
subroutine parse(filename,soc,spin,kpoint,band,x,y)
    implicit none
    character(75) :: filename
    character(75) :: soc
    integer :: spin,kpoint,band
    real(8) :: x,y
    !!$ 键-值 键:--option,值:--value
    character(20) :: option,value
    !! iarg 表示程序后面参数key value总个数，narg表示参数中value个数
    integer :: iarg,narg,ia,nargs
    !!$  输入的参数必须是偶数
    iarg=iargc()
    nargs=iarg/2
    !!default value
    filename="WAVECAR"
    soc="t"
    spin = 1
    kpoint = 1
    band = 1
    x = 0.
    y = 0.
    if(iarg.ne.2*nargs) then
       call getarg(1,option)
       if (option .eq. "-h") then
           call help
           stop
       endif    
       write(6,*) 'program error, key number mismatches value number'
       call help
    endif
    do ia=1,nargs
       call getarg(2*ia-1,option)
       call getarg(2*ia,value)
       if(option == "-f") then
          read(value,*) filename
       else if(option == "-soc") then
          read(value,*) soc
          if ((soc .ne. 't') .and. (soc .ne. 'f')) then
            write(6,*) "wrong parameter!!! soc value must be chosen from t or f"
            stop
          end if
       else if(option == "-spin") then
          read(value,*) spin
       else if(option == "-k") then
          read(value,*) kpoint
       else if(option == "-b") then
          read(value,*) band
       else if(option == "-x") then
          read(value,*) x
       else if(option == "-y") then
          read(value,*) y
       else
          call help
       endif
    enddo
    return
end subroutine parse

subroutine help
  write(6,*) "syntax: Wavetran -f file  -soc t/f -spin spin -k k-point -b band -x coord -y coord"
  write(6,*) "defaults parameter: -f WAVECAR -soc t -spin 1 -k 1 -b 1 -x 0.0 -y 0.0"
  write(6,*) "-soc t: means true (consider spin orbit coupling)"
  write(6,*) "inputs: x and y are direct coordinates on axes a1 and a2"
  write(6,*) "output: wave function phi(x,y,z) (chi(x,y,z) consider soc) at position r(x,y,z) with z direct coordinate on a3 axis"
  stop
end subroutine help
