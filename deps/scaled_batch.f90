! Compiled separately with each backend's param module. Transfer native
! mantissas and integer decimal exponents without reconstructing large values.
module scaled_batch
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use param, only: knd
  use prolate_swf, only: profcn
  use oblate_swf, only: oblfcn
  use complex_prolate_swf, only: cprofcn
  use complex_oblate_swf, only: coblfcn
  implicit none
contains
  subroutine spheroidal_scaled_text(geometry, complex_input, mode, option, m, n, count, width, &
                                    ctext, xtext, output, exponents, status) bind(C, name="spheroidal_scaled_text")
    integer(c_int), value :: geometry, complex_input, mode, option, m, n, count, width
    character(c_char), intent(in) :: ctext(*), xtext(*)
    character(c_char), intent(out) :: output(*)
    integer(c_int), intent(out) :: exponents(*), status
    call scaled_degrees(geometry,complex_input,mode,option,m,n,n,count,width,ctext,xtext,output,exponents,status)
  end subroutine

  subroutine spheroidal_degrees_scaled_text(geometry,mode,option,m,n_min,n,count,width, &
                  ctext,xtext,output,exponents,status) bind(C,name="spheroidal_degrees_scaled_text")
    integer(c_int), value :: geometry,mode,option,m,n_min,n,count,width
    character(c_char), intent(in) :: ctext(*),xtext(*)
    character(c_char), intent(out) :: output(*)
    integer(c_int), intent(out) :: exponents(*),status
    call scaled_degrees(geometry,0,mode,option,m,n_min,n,count,width,ctext,xtext,output,exponents,status)
  end subroutine

  subroutine spheroidal_degrees_scaled_double(geometry,mode,option,m,n_min,n,count,width, &
                  ctext,xtext,output,exponents,status) bind(C,name="spheroidal_degrees_scaled_double")
    integer(c_int), value :: geometry,mode,option,m,n_min,n,count,width
    character(c_char), intent(in) :: ctext(*),xtext(*)
    real(c_double), intent(out) :: output(*)
    character(c_char) :: unused(1)
    integer(c_int), intent(out) :: exponents(*),status
    call scaled_degrees(geometry,0,mode,option,m,n_min,n,count,width,ctext,xtext,unused,exponents,status,output)
  end subroutine

  subroutine scaled_degrees(geometry, complex_input, mode, option, m, n_min, n, count, width, &
                            ctext, xtext, output, exponents, status, double_output)
    integer(c_int), value :: geometry, complex_input, mode, option, m, n, count, width
    integer(c_int), value :: n_min
    character(c_char), intent(in) :: ctext(*), xtext(*)
    character(c_char), intent(out) :: output(*)
    integer(c_int), intent(out) :: exponents(*), status
    real(c_double), optional, intent(out) :: double_output(*)
    integer :: li, lnum, narg, irad, iang, norm, i, j, repeats, threshold, first_li, offset
    real(knd) :: cr, ci, x
    complex(knd) :: c
    logical :: ok
    real(knd), allocatable :: points(:), arg(:), r(:,:), s(:,:), ds(:,:), eig(:)
    complex(knd), allocatable :: z(:,:), a(:,:), da(:,:), zeig(:)
    integer, allocatable :: er(:,:), es(:,:), eds(:,:), acc(:), ac(:,:), adc(:,:)
    status = 0
    if (m < 0 .or. n_min < m .or. n < n_min .or. count < 1 .or. width < 70 .or. &
        geometry < 0 .or. geometry > 1 .or. complex_input < 0 .or. complex_input > 1 .or. &
        mode < 1 .or. mode > 2) then
      status = -4
      return
    end if
    if ((mode == 1 .and. (option < 0 .or. option > 1)) .or. &
        (mode == 2 .and. (option < 1 .or. option > 4))) then
      status = -4
      return
    end if
    call decode(ctext,1,cr,ok)
    if (.not. ok) return
    call decode(ctext,2,ci,ok)
    if (.not. ok) return
    c = cmplx(cr,ci,knd)
    if (abs(c) == 0 .or. (complex_input == 0 .and. cr <= 0)) then
      status = -5
      return
    end if
    allocate(points(count))
    do i=1,count
      call decode(xtext,i,points(i),ok)
      if (.not. ok) return
      if ((mode == 1 .and. abs(points(i)) > 1) .or. (mode == 2 .and. points(i) < 0)) then
        status = -3
        return
      end if
    end do
    ! Radial prolate coordinates arrive as x-1, preserving boundary distance.
    li = n-m+1
    first_li = n_min-m+1
    lnum = li
    if (geometry == 1) then
      threshold = int(2*abs(c)/acos(-1.0_knd))
      if (lnum < threshold .and. mod(lnum,2) /= 0) lnum=lnum+1
    end if
    narg = 1
    if (mode == 1) narg=count
    allocate(arg(narg),r(lnum,4),z(lnum,4),s(lnum,narg),ds(lnum,narg),eig(lnum))
    allocate(a(lnum,narg),da(lnum,narg),zeig(lnum),er(lnum,4),es(lnum,narg),eds(lnum,narg))
    allocate(acc(lnum),ac(lnum,narg),adc(lnum,narg))
    arg=0
    if (mode == 1) arg=points
    irad=0
    iang=0
    norm=0
    repeats=1
    if (mode == 1) then
      iang=2
      norm=option
    else
      irad=2
      if (option == 1) irad=1
      repeats=count
    end if
    do j=1,repeats
      r=0; z=0; er=0; acc=-1
      x=1
      if (geometry == 1) x=10
      if (mode == 2) x=points(j)
      if (complex_input == 0) then
        if (geometry == 0) then
          call profcn(cr,m,lnum,irad,x,iang,norm,narg,arg,r(:,1),er(:,1),r(:,2),er(:,2), &
                      r(:,3),er(:,3),r(:,4),er(:,4),acc,s,es,ds,eds,ac,eig)
        else
          call oblfcn(cr,m,lnum,irad,x,iang,norm,narg,arg,r(:,1),er(:,1),r(:,2),er(:,2), &
                      r(:,3),er(:,3),r(:,4),er(:,4),acc,s,es,ds,eds,ac,eig)
        end if
        z=cmplx(r,0.0_knd,knd)
        if (mode == 1) then
          a=cmplx(s,0.0_knd,knd)
          da=cmplx(ds,0.0_knd,knd)
        end if
      else
        if (geometry == 0) then
          call cprofcn(c,m,lnum,irad,x,iang,norm,narg,arg,z(:,1),er(:,1),z(:,2),er(:,2), &
                       z(:,3),er(:,3),z(:,4),er(:,4),acc,a,es,da,eds,ac,adc,zeig)
        else
          call coblfcn(c,m,lnum,irad,x,iang,norm,narg,arg,z(:,1),er(:,1),z(:,2),er(:,2), &
                       z(:,3),er(:,3),z(:,4),er(:,4),acc,a,es,da,eds,ac,adc,zeig)
        end if
      end if
      do li=first_li,n-m+1
        offset=8*count*(li-first_li)
        if (mode == 1) then
          do i=1,count
            call store(offset+8*(i-1)+1,a(li,i),es(li,i))
            call store(offset+8*(i-1)+3,da(li,i),eds(li,i))
            call store(offset+8*(i-1)+5,cmplx(0.0_knd,0.0_knd,knd),0)
            call store(offset+8*(i-1)+7,cmplx(0.0_knd,0.0_knd,knd),0)
          end do
        else
          do i=1,4
            call store(offset+8*(j-1)+2*i-1,z(li,i),er(li,i))
          end do
        end if
      end do
    end do
  contains
    subroutine decode(buffer,index,value,ok)
      character(c_char), intent(in) :: buffer(*)
      integer, intent(in) :: index
      real(knd), intent(out) :: value
      logical, intent(out) :: ok
      character(len=width) :: text
      integer :: k,ios
      do k=1,width
        text(k:k)=buffer((index-1)*width+k)
      end do
      read(text,*,iostat=ios) value
      ok=ios == 0
      if (ok) ok=ieee_is_finite(value)
      if (.not. ok) status=-5
    end subroutine
    subroutine store(index,value,exponent)
      integer, intent(in) :: index,exponent
      complex(knd), intent(in) :: value
      character(len=width) :: text
      integer :: k,channel
      do channel=0,1
        if (present(double_output)) then
          if (channel == 0) then
            double_output(index+channel)=real(value,c_double)
          else
            double_output(index+channel)=real(aimag(value),c_double)
          end if
          exponents(index+channel)=exponent
          cycle
        end if
        if (channel == 0) then
          write(text,'(ES70.60E4)') real(value,knd)
        else
          write(text,'(ES70.60E4)') aimag(value)
        end if
        do k=1,width
          output((index+channel-1)*width+k)=text(k:k)
        end do
        exponents(index+channel)=exponent
      end do
    end subroutine
  end subroutine
end module
