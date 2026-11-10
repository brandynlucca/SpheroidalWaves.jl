program native_mpfr_arithmetic
  use, intrinsic :: iso_c_binding, only: c_char,c_null_char,c_int,c_int64_t
  use spheroidal_mpfr_complex
  use spheroidal_mpfr_batch_api, only: spheroidal_mpfr_batch
  use prolate_swf_mpfr, only: prolate_qleg_cache => qleg_cache
  use complex_prolate_swf_mpfr, only: complex_prolate_qleg_cache => qleg_cache
  use omp_lib
  implicit none
  type(mp_real) :: a,b,c,tolerance,values(4),copies(4)
  type(mp_complex) :: z,w,v,complex_values(3),complex_copies(3)
  integer :: bits,mode,i,thread,converted
  integer, parameter :: width=512
  character(kind=c_char) :: parameter(width,2),coordinates(width,1),output(width,10,1)
  integer(c_int) :: status,accuracy(1)
  do bits=128,512,128
    mpfr_bits=bits
    mpfr_rounding=0
    tolerance=mpreal(2)**(-bits+8)
    if (sqrt(mpreal(2))/=mpreal('1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273501384623091229702492483605585073721264412149709993583141322266592750559275579995050115278206057147011')) error stop 1
    if (exp(mpreal(1))/=mpreal('2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011574')) error stop 2
    a=mpreal('0.1')
    b=a
    a=100
    if (b/=mpreal('0.1')) error stop 3
    values=[mpreal(1),mpreal(2),mpreal(3),mpreal(4)]
    copies=values
    values=0
    if (sum(copies)/=10.or.maxval(copies)/=4.or.minval(copies)/=1) error stop 4
    if (max(1,2,mpreal(3),4)/=4.or.min(1d0,2,3,mpreal(4))/=1) error stop 5
    if (abs(log(exp(mpreal(1)))-1)>tolerance) error stop 6
    if (abs(sin(mpreal(1))**2+cos(mpreal(1))**2-1)>tolerance) error stop 7
    if (abs(cosh(mpreal(1))**2-sinh(mpreal(1))**2-1)>tolerance) error stop 8
    if (abs(atan(mpreal(1))*4-acos(mpreal(-1)))>tolerance) error stop 9
    if (nint(mpreal('2.5'))/=3.or.nint(mpreal('-2.5'))/=-3) error stop 10
    converted=mpreal('-2.9')
    if (converted/=-2.or.int(mpreal('-2.9'))/=-2) error stop 11
    if (sign(mpreal(-2),mpreal(3))/=2.or.mod(mpreal(7),mpreal(3))/=1) error stop 12
    if (.not.(huge(b)>b.and.tiny(b)>0.and.tiny(b)<b)) error stop 13
    if (mpreal(mp_to_string(b))/=b) error stop 14
    if (.not.mp_is_finite(b).or.mp_is_finite(mpreal('Inf')).or.mp_is_finite(mpreal('NaN'))) error stop 36
    do mode=0,4
      mpfr_rounding=mode
      a=mpreal(1)/3
      if (a/=mpreal('1/3')) error stop 15
      if (a%data(1)/=bits) error stop 16
    end do
    mpfr_rounding=0
    z=mpcmplx(mpreal('0.125'),mpreal('-0.375'))
    w=mpcmplx(mpreal('0.75'),mpreal('0.0625'))
    v=z*w
    if (real(v)/=mpreal('0.1171875').or.aimag(v)/=mpreal('-0.2734375')) error stop 17
    if (abs(z/w*w-z)>tolerance.or.abs(z**(-3)*z**3-1)>tolerance) error stop 18
    if (abs(sqrt(z)**2-z)>tolerance.or.abs(exp(log(z))-z)>tolerance) error stop 19
    if (abs(sin(z)**2+cos(z)**2-1)>tolerance) error stop 20
    if (abs(cosh(z)**2-sinh(z)**2-1)>tolerance) error stop 21
    if (abs(z**mpreal('0.5')-sqrt(z))>tolerance) error stop 22
    if (abs(abs(z)**2-mpreal('0.15625'))>tolerance) error stop 23
    if (conjg(z)/=mpcmplx(mpreal('0.125'),mpreal('0.375'))) error stop 24
    v=cmplx(0.125,-0.375)
    if (v/=z) error stop 25
    v=cmplx(0.125d0,-0.375d0)
    if (v/=z) error stop 26
    complex_values=z
    complex_copies=complex_values
    complex_values=complex_values*2+mpreal(1)
    if (any(abs(complex_copies-z)>tolerance)) error stop 27
    if (any(abs(complex_values-(z*2+1))>tolerance)) error stop 28
    complex_values=mpreal('0.5')
    if (any(complex_values/=mpreal('0.5'))) error stop 29
    c=z
    if (c/=mpreal('0.125')) error stop 30
    do i=1,10000
      v=(z+w)*(z-w)/(z*z-w*w)
      if (abs(v-1)>tolerance) error stop 31
      c=(mpreal(i)/7+sqrt(mpreal(2)))**2
      b=c
      c=0
      if (b<=0) error stop 32
    end do
  end do
!$omp parallel private(thread,a,b,c,z,w)
  thread=omp_get_thread_num()
  mpfr_bits=160+32*thread
  mpfr_rounding=mod(thread,5)
  a=mpreal('1/7')
  b=a
  c=sqrt(a)
  if (c%data(1)/=mpfr_bits.or.a%data(1)/=mpfr_bits) error stop 33
  a=0
  if (b/=mpreal('1/7')) error stop 34
  z=mpcmplx(b,mpreal('1/11'))
  w=z*z
  if (w%re%data(1)/=mpfr_bits.or.w%im%data(1)/=mpfr_bits) error stop 35
!$omp end parallel
  call check_batch()
  print *, 'MPFR/MPC arithmetic, ownership, rounding, ABI accuracy and bounded failures passed'
contains
  subroutine check_batch()
    integer :: fixture,requested_bits,geometry,complex_parameter,slot
    integer(c_int64_t) :: saved_limit
    type(mp_real) :: error,previous_error,reference_re,reference_im
    character(len=112), parameter :: references_re(4)=[character(len=112) :: &
      '0.319000055146892739783981985871826494103814376510741729114970059238053009763670379743178473870059809926688', &
      '-0.348602399470269089476372382847354292920767569504992053128206743341920042616829943858732111881131616819626', &
      '0.303235577304105060785313344346188701012396500405900347760112187638888468163625765174587690518159352874124', &
      '-0.321883129754014593679374588303615476083721858296693323713886315177274424476460712582066203735686938450904']
    character(len=112), parameter :: references_im(4)=[character(len=112) :: '0','0', &
      '0.153352965506012687954704507348706414785976931529398947646079078147847292101453723983951842226194416353958', &
      '-0.181098665410776759271204693665910000710095166563737449848924342434395006230367264607585795845087562112543']
    mpfr_bits=384
    mpfr_rounding=0
    call fill_text(parameter(:,1),'1')
    call fill_text(coordinates(:,1),'0')
    ! Independent mpmath eigenvalue references and refinement provenance are in arbitrary_precision.jl.
    do fixture=1,4
      geometry=mod(fixture,2)
      complex_parameter=merge(1,0,fixture>2)
      if (complex_parameter==1) then
        call fill_text(parameter(:,2),'1/4')
      else
        call fill_text(parameter(:,2),'0')
      end if
      reference_re=mpreal(trim(references_re(fixture)))
      reference_im=mpreal(trim(references_im(fixture)))
      do requested_bits=160,256,96
        call evaluate(geometry,complex_parameter,0,0,0,requested_bits,0)
        error=max(abs(component(1)-reference_re),abs(component(2)-reference_im))
        if (error>mpreal(2)**(-requested_bits+12)) error stop 40
        if (requested_bits==256) then
          if (error>=previous_error*mpreal(2)**(-70)) error stop 41
        end if
        previous_error=error
      end do
    end do
    call fill_text(parameter(:,2),'0')
    call fill_text(coordinates(:,1),'1/4')
    ! Native angular output omits the Condon-Shortley phase in accuracy_references.tsv.
    do requested_bits=160,256,96
      call evaluate(1,0,1,1,1,requested_bits,0)
      if (abs(component(3)+mpreal('-0.981030364528480149520801385335891232434186368'))>mpreal('1e-40')) error stop 42
      if (abs(component(5)+mpreal('0.309649688659690203591415049178354606857683011'))>mpreal('1e-40')) error stop 43
      call evaluate(0,0,1,1,1,requested_bits,0)
      if (abs(component(3)+mpreal('-0.954376874057940932495574682020773761538332828'))>mpreal('1e-40')) error stop 44
      if (abs(component(5)+mpreal('0.205752668742519211917950856913568615675196084'))>mpreal('1e-40')) error stop 45
    end do
    call evaluate(0,0,0,1000000,0,160,4)
    call evaluate(1,0,0,0,0,64000000,4)
    call evaluate(1,0,0,0,0,160,0)
    if (abs(component(1)-mpreal(trim(references_re(1))))>mpreal(2)**(-148)) error stop 48
    saved_limit=mpfr_workspace_limit_bytes
    mpfr_workspace_limit_bytes=1048576_c_int64_t
    call evaluate(1,0,0,0,0,160,4)
    mpfr_workspace_limit_bytes=saved_limit
    call fill_text(parameter(:,1),'NaN')
    call evaluate(1,0,0,0,0,160,1)
    call fill_text(parameter(:,1),'1')
    call evaluate(1,0,0,0,0,160,0)
    call fill_text(parameter(:,2),'1/4')
    call fill_text(coordinates(:,1),'2')
    call evaluate(1,1,0,1,3,256,0)
    z=mpcmplx(component(3),component(4))
    w=mpcmplx(component(5),component(6))
    v=z*mpcmplx(component(9),component(10))-w*mpcmplx(component(7),component(8))
    if (.not.mp_is_finite(real(v)).or..not.mp_is_finite(aimag(v))) error stop 51
    if (abs(v-1/(3*mpcmplx(mpreal(1),mpreal('1/4'))))>mpreal(2)**(-80)) error stop 52
    call fill_text(parameter(:,2),'0')
    call evaluate(0,0,0,1,3,1024,0)
    error=component(3)*component(9)-component(5)*component(7)
    if (.not.mp_is_finite(error).or.abs(error-mpreal('1/5'))>mpreal(2)**(-80)) error stop 53
    do slot=1,size(prolate_qleg_cache)
      if (allocated(prolate_qleg_cache(slot)%x1%data).or.allocated(prolate_qleg_cache(slot)%qdml%data).or. &
          allocated(prolate_qleg_cache(slot)%qml%data).or.allocated(prolate_qleg_cache(slot)%termpq%data)) error stop 49
    end do
    do slot=1,size(complex_prolate_qleg_cache)
      if (allocated(complex_prolate_qleg_cache(slot)%x1%data).or.allocated(complex_prolate_qleg_cache(slot)%qdml%data).or. &
          allocated(complex_prolate_qleg_cache(slot)%qml%data).or.allocated(complex_prolate_qleg_cache(slot)%termpq%data)) error stop 50
    end do
  end subroutine

  subroutine evaluate(geometry,complex_parameter,m,n,mode,bits,expected_status)
    integer, intent(in) :: geometry,complex_parameter,m,n,mode,bits,expected_status
    call spheroidal_mpfr_batch(geometry,complex_parameter,m,n,mode,0,bits,0,1,width, &
      parameter,coordinates,output,accuracy,status)
    if (status/=expected_status) then
      print *, 'MPFR ABI status mismatch:',geometry,complex_parameter,m,n,mode,bits,status,expected_status
      error stop 46
    end if
    if (mpfr_bits/=384.or.mpfr_rounding/=0) error stop 47
  end subroutine

  function component(channel) result(value)
    integer, intent(in) :: channel
    type(mp_real) :: value
    character(len=width) :: text
    integer :: index
    text=' '
    do index=1,width
      if (output(index,channel,1)==c_null_char) exit
      text(index:index)=output(index,channel,1)
    end do
    value=0
    if (len_trim(text)>0) value=mpreal(trim(text))
  end function

  subroutine fill_text(output,text)
    character(kind=c_char), intent(out) :: output(:)
    character(len=*), intent(in) :: text
    integer :: index
    output=c_null_char
    do index=1,len(text)
      output(index)=text(index:index)
    end do
  end subroutine
end program native_mpfr_arithmetic
