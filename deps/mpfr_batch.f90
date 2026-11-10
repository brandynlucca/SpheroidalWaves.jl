module spheroidal_mpfr_batch_api
  use, intrinsic :: iso_c_binding
  use spheroidal_mpfr_real
  use spheroidal_mpfr_complex
  use prolate_swf_mpfr, only: profcn, clear_prolate => mpfr_clear_caches
  use oblate_swf_mpfr, only: oblfcn, clear_oblate => mpfr_clear_caches
  use complex_prolate_swf_mpfr, only: cprofcn, clear_complex_prolate => mpfr_clear_caches
  use complex_oblate_swf_mpfr, only: coblfcn, clear_complex_oblate => mpfr_clear_caches
  implicit none
  private
  public :: spheroidal_mpfr_batch
  interface
    function parse_real(output,text,rounding) bind(C,name='sw_mpfr_set_string') result(status)
      import
      integer(c_int64_t), intent(inout) :: output(*)
      character(c_char), intent(in) :: text(*)
      integer(c_int), value :: rounding
      integer(c_int) :: status
    end function
  end interface
contains
  subroutine clear_caches()
    call clear_prolate()
    call clear_oblate()
    call clear_complex_prolate()
    call clear_complex_oblate()
  end subroutine

  subroutine write_real(value,output,status)
    type(mp_real), intent(in) :: value
    character(c_char), intent(out) :: output(:)
    integer(c_int), intent(inout) :: status
    character(len=:), allocatable :: text
    integer :: index
    text=mp_to_string(value)
    output=c_null_char
    if (len(text)>=size(output)) then
      status=3
      return
    end if
    do index=1,len(text)
      output(index)=text(index:index)
    end do
  end subroutine

  subroutine write_complex(value,output,status)
    type(mp_complex), intent(in) :: value
    character(c_char), intent(out) :: output(:,:)
    integer(c_int), intent(inout) :: status
    call write_real(value%re,output(:,1),status)
    call write_real(value%im,output(:,2),status)
  end subroutine

  subroutine spheroidal_mpfr_batch(geometry,complex_parameter,m,n,mode,normalize,bits,rounding, &
      npoints,text_length,parameter,coordinates,output,accuracy,status) bind(C)
    integer(c_int), value :: geometry,complex_parameter,m,n,mode,normalize,bits,rounding,npoints,text_length
    character(c_char), intent(in) :: parameter(text_length,2),coordinates(text_length,npoints)
    character(c_char), intent(out) :: output(text_length,10,npoints)
    integer(c_int), intent(out) :: accuracy(npoints),status
    type(mp_real) :: c_real,c_imag,x
    type(mp_complex) :: c
    type(mp_real), allocatable :: eta(:),r1(:),r1d(:),r2(:),r2d(:),s(:,:),sd(:,:),eigenvalues(:)
    type(mp_complex), allocatable :: cr1(:),cr1d(:),cr2(:),cr2d(:),cs(:,:),csd(:,:),ceigenvalues(:)
    integer, allocatable :: r1e(:),r1de(:),r2e(:),r2de(:),se(:,:),sde(:,:),radial_accuracy(:)
    integer, allocatable :: angular_accuracy(:,:),angular_derivative_accuracy(:,:)
    integer :: saved_bits,saved_rounding,saved_status,lnum,narg,point,ioprad,iopang,calls,index
    integer(c_int64_t) :: element_bytes
    status=0
    output=c_null_char
    accuracy=0
    if (geometry<0.or.geometry>1.or.complex_parameter<0.or.complex_parameter>1.or. &
        m<0.or.n<m.or.mode<0.or.mode>3.or.normalize<0.or.normalize>1.or. &
        bits<2.or.rounding<0.or.rounding>4.or.npoints<1.or.text_length<2) then
      status=1
      return
    end if
    if (any(parameter(text_length,:)/=c_null_char).or.any(coordinates(text_length,:)/=c_null_char)) then
      status=1
      return
    end if
    element_bytes=(storage_size(c_real,kind=c_int64_t)+7_c_int64_t)/8_c_int64_t &
      +88_c_int64_t+8_c_int64_t*((int(bits,c_int64_t)+63_c_int64_t)/64_c_int64_t)
    if (int(n,c_int64_t)-m+1>mpfr_workspace_limit_bytes/(32*element_bytes).or. &
        int(npoints,c_int64_t)>mpfr_workspace_limit_bytes/(32*element_bytes).or. &
        int(n,c_int64_t)+bits>1000000_c_int64_t) then
      status=4
      return
    end if
    lnum=n-m+1
    narg=1
    if (mode==1) narg=npoints
    if (int(lnum,c_int64_t)>mpfr_workspace_limit_bytes/(32*element_bytes)/int(narg,c_int64_t)) then
      status=4
      return
    end if
    saved_bits=mpfr_bits
    saved_rounding=mpfr_rounding
    saved_status=mpfr_status
    mpfr_bits=bits
    mpfr_rounding=rounding
    mpfr_status=0
    call clear_caches()
    call mp_init(c_real)
    call mp_init(c_imag)
    if (parse_real(c_real%data,parameter(:,1),rounding)/=0.or. &
        parse_real(c_imag%data,parameter(:,2),rounding)/=0) then
      status=1
      goto 900
    end if
    if (.not.mp_is_finite(c_real).or..not.mp_is_finite(c_imag)) then
      status=1
      goto 900
    end if
    if (c_real<0.or.c_real==0.or.(complex_parameter==0.and.c_imag/=0)) then
      status=2
      goto 900
    end if
    c=mpcmplx(c_real,c_imag)
    if (abs(c_real)+abs(c_imag)>1000000) then
      status=4
      goto 900
    end if
    allocate(eta(narg),r1e(lnum),r1de(lnum),r2e(lnum),r2de(lnum),se(lnum,narg),sde(lnum,narg))
    allocate(radial_accuracy(lnum),angular_accuracy(lnum,narg),angular_derivative_accuracy(lnum,narg))
    if (complex_parameter==1) then
      allocate(cr1(lnum),cr1d(lnum),cr2(lnum),cr2d(lnum),cs(lnum,narg),csd(lnum,narg),ceigenvalues(lnum))
    else
      allocate(r1(lnum),r1d(lnum),r2(lnum),r2d(lnum),s(lnum,narg),sd(lnum,narg),eigenvalues(lnum))
    end if
    eta=0
    if (mode==1) then
      do point=1,narg
        call mp_init(eta(point))
        if (parse_real(eta(point)%data,coordinates(:,point),rounding)/=0) then
          status=1
          goto 900
        end if
        if (.not.mp_is_finite(eta(point)).or.abs(eta(point))>1) then
          status=2
          goto 900
        end if
      end do
    end if
    iopang=0
    if (mode==1) iopang=2
    ioprad=max(0,mode-1)
    calls=1
    if (mode>=2) calls=npoints
    do point=1,calls
      x=1
      if (mode>=2) then
        call mp_init(x)
        if (parse_real(x%data,coordinates(:,point),rounding)/=0) then
          status=1
          goto 900
        end if
        if (.not.mp_is_finite(x).or.x<geometry.or.(geometry==1.and.x==1.and.ioprad==2)) then
          status=2
          goto 900
        end if
        if (geometry==1) x=x-1
      end if
      radial_accuracy=0
      angular_accuracy=0
      angular_derivative_accuracy=0
      if (complex_parameter==1) then
        if (geometry==1) then
          call cprofcn(c,m,lnum,ioprad,x,iopang,normalize,narg,eta,cr1,r1e,cr1d,r1de, &
              cr2,r2e,cr2d,r2de,radial_accuracy,cs,se,csd,sde,angular_accuracy,angular_derivative_accuracy,ceigenvalues)
        else
          call coblfcn(c,m,lnum,ioprad,x,iopang,normalize,narg,eta,cr1,r1e,cr1d,r1de, &
              cr2,r2e,cr2d,r2de,radial_accuracy,cs,se,csd,sde,angular_accuracy,angular_derivative_accuracy,ceigenvalues)
        end if
        if (mpfr_status/=0) then
          status=mpfr_status
          goto 900
        end if
        if (mode<=1) then
          do index=1,npoints
            call write_complex(ceigenvalues(lnum),output(:,1:2,index),status)
            if (mode==0) cycle
            call write_complex(cs(lnum,index)*mpreal(10)**se(lnum,index),output(:,3:4,index),status)
            call write_complex(csd(lnum,index)*mpreal(10)**sde(lnum,index),output(:,5:6,index),status)
            accuracy(index)=min(angular_accuracy(lnum,index),angular_derivative_accuracy(lnum,index))
          end do
        else
          call write_complex(ceigenvalues(lnum),output(:,1:2,point),status)
          call write_complex(cr1(lnum)*mpreal(10)**r1e(lnum),output(:,3:4,point),status)
          call write_complex(cr1d(lnum)*mpreal(10)**r1de(lnum),output(:,5:6,point),status)
          if (ioprad==2) then
            call write_complex(cr2(lnum)*mpreal(10)**r2e(lnum),output(:,7:8,point),status)
            call write_complex(cr2d(lnum)*mpreal(10)**r2de(lnum),output(:,9:10,point),status)
          end if
          accuracy(point)=radial_accuracy(lnum)
        end if
      else
        if (geometry==1) then
          call profcn(c_real,m,lnum,ioprad,x,iopang,normalize,narg,eta,r1,r1e,r1d,r1de, &
              r2,r2e,r2d,r2de,radial_accuracy,s,se,sd,sde,angular_accuracy,eigenvalues)
        else
          call oblfcn(c_real,m,lnum,ioprad,x,iopang,normalize,narg,eta,r1,r1e,r1d,r1de, &
              r2,r2e,r2d,r2de,radial_accuracy,s,se,sd,sde,angular_accuracy,eigenvalues)
        end if
        if (mpfr_status/=0) then
          status=mpfr_status
          goto 900
        end if
        if (mode<=1) then
          do index=1,npoints
            call write_real(eigenvalues(lnum),output(:,1,index),status)
            if (mode==0) cycle
            call write_real(s(lnum,index)*mpreal(10)**se(lnum,index),output(:,3,index),status)
            call write_real(sd(lnum,index)*mpreal(10)**sde(lnum,index),output(:,5,index),status)
            accuracy(index)=angular_accuracy(lnum,index)
          end do
        else
          call write_real(eigenvalues(lnum),output(:,1,point),status)
          call write_real(r1(lnum)*mpreal(10)**r1e(lnum),output(:,3,point),status)
          call write_real(r1d(lnum)*mpreal(10)**r1de(lnum),output(:,5,point),status)
          if (ioprad==2) then
            call write_real(r2(lnum)*mpreal(10)**r2e(lnum),output(:,7,point),status)
            call write_real(r2d(lnum)*mpreal(10)**r2de(lnum),output(:,9,point),status)
          end if
          accuracy(point)=radial_accuracy(lnum)
        end if
      end if
      if (status/=0) exit
    end do
900 continue
    call clear_caches()
    mpfr_bits=saved_bits
    mpfr_rounding=saved_rounding
    mpfr_status=saved_status
  end subroutine
end module
