#include "precision.inc"

program prolate_cache_test
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use SW_PRECISION_MODULE
  use SW_PROLATE_MODULE, only: auxiliary_cache_limit_bytes, cache_real_equal, &
    pleg_cache, qleg_cache, gauss_cache, clear_pleg_cache_slot, clear_qleg_cache_slot, &
    find_pleg_cache, store_pleg_cache, load_pleg_cache, find_qleg_cache, store_qleg_cache, &
    qleg_cached, gauss_cached, reserve_cache_bytes, add_cache_bytes, profcn
#ifdef SWF_PRECISION_MPFR
  use SW_PROLATE_MODULE, only: mpfr_clear_caches
#endif
  implicit none
  integer(int64) :: saved_budget
#ifdef SWF_PRECISION_MPFR
  integer, parameter :: bits(3) = [160,256,160], modes(3) = [0,3,2]
  integer :: pass,saved_bits,saved_rounding,saved_status
  integer(int64) :: first_bytes
  saved_bits = mpfr_bits
  saved_rounding = mpfr_rounding
  saved_status = mpfr_status
#endif
  saved_budget = auxiliary_cache_limit_bytes
#ifdef SWF_PRECISION_MPFR
  do pass=1,size(bits)
    mpfr_bits = bits(pass)
    mpfr_rounding = modes(pass)
    mpfr_status = 0
    if(pass==1) first_bytes = sw_real_bytes()
    if(pass==2.and.sw_real_bytes()<=first_bytes) error stop "precision-dependent cache bytes"
    if(pass==3.and.sw_real_bytes()/=first_bytes) error stop "restored cache bytes"
    call check_caches()
    if(mpfr_status/=0) error stop "cache arithmetic status"
  end do
  mpfr_bits = saved_bits
  mpfr_rounding = saved_rounding
  mpfr_status = saved_status
#else
  call check_caches()
#endif
  auxiliary_cache_limit_bytes = saved_budget
  call check_radial_fallback()
  print *, "Prolate cache policies, workspace accounting, endpoint and radial fallback passed"
contains
  subroutine check_radial_fallback()
#ifdef SWF_PRECISION_QUAD
    integer, parameter :: parameter_value=500,order=100,degrees=300,selected_degree=264
#elif defined(SWF_PRECISION_MPFR)
    integer, parameter :: parameter_value=200,order=50,degrees=80,selected_degree=107
#else
    integer, parameter :: parameter_value=200,order=50,degrees=80,selected_degree=108
#endif
    integer, parameter :: selected=selected_degree-order+1
    SW_REAL :: c,x1,arguments(1),r1(degrees),r1d(degrees),r2(degrees),r2d(degrees)
    SW_REAL :: s1(degrees,1),s1d(degrees,1),eigenvalues(degrees),values(4),ten,tolerance,wronskian
    integer :: r1e(degrees),r1de(degrees),r2e(degrees),r2de(degrees),se(degrees,1),sde(degrees,1)
    integer :: radial_accuracy(degrees),angular_accuracy(degrees,1)
    integer(int64) :: cache_budget
#ifdef SWF_PRECISION_MPFR
    integer :: caller_bits,caller_rounding,caller_status
    integer(int64) :: workspace_budget
    caller_bits = mpfr_bits
    caller_rounding = mpfr_rounding
    caller_status = mpfr_status
    workspace_budget = mpfr_workspace_limit_bytes
#endif
    cache_budget = auxiliary_cache_limit_bytes
    call clear_fallback_caches()
#ifdef SWF_PRECISION_MPFR
    mpfr_bits = 113
    mpfr_rounding = 0
    mpfr_status = 0
    mpfr_workspace_limit_bytes = 536870912_int64
    tolerance = SW_LITERAL(1.0e-24_sw_kind,"1e-24")
#elif defined(SWF_PRECISION_QUAD)
    tolerance = SW_LITERAL(1.0e-13_sw_kind,"1e-13")
#else
    tolerance = SW_LITERAL(1.0e-7_sw_kind,"1e-7")
#endif
    c = sw_real(parameter_value)
    x1 = SW_LITERAL(0.1_sw_kind,"0.1")
    arguments = 0
    ten = SW_LITERAL(10.0_sw_kind,"10.0")
    call profcn(c,order,degrees,2,x1,0,0,1,arguments,r1,r1e,r1d,r1de,r2,r2e,r2d,r2de, &
      radial_accuracy,s1,se,s1d,sde,angular_accuracy,eigenvalues)
#ifdef SWF_PRECISION_MPFR
    if(mpfr_status/=0) error stop "radial fallback status"
#endif
    values = [r1(selected)*ten**r1e(selected),r1d(selected)*ten**r1de(selected), &
      r2(selected)*ten**r2e(selected),r2d(selected)*ten**r2de(selected)]
    wronskian = c*x1*(x1+2)*(values(1)*values(4)-values(2)*values(3))
#ifdef SWF_PRECISION_MPFR
    if(any(.not.mp_is_finite(values)).or..not.mp_is_finite(wronskian)) error stop "finite radial fallback"
#else
    if(any(.not.ieee_is_finite(values)).or..not.ieee_is_finite(wronskian)) error stop "finite radial fallback"
#endif
    if(abs(wronskian-1)>tolerance) error stop "radial fallback Wronskian"
    call clear_fallback_caches()
    auxiliary_cache_limit_bytes = cache_budget
#ifdef SWF_PRECISION_MPFR
    mpfr_bits = caller_bits
    mpfr_rounding = caller_rounding
    mpfr_status = caller_status
    mpfr_workspace_limit_bytes = workspace_budget
#endif
  end subroutine

  subroutine clear_fallback_caches()
#ifdef SWF_PRECISION_MPFR
    call mpfr_clear_caches()
#else
    integer(int64) :: budget
    logical :: fits
    budget = auxiliary_cache_limit_bytes
    auxiliary_cache_limit_bytes = 0_int64
    fits = reserve_cache_bytes(0_int64)
    auxiliary_cache_limit_bytes = budget
    if(.not.fits) error stop "radial fallback cache cleanup"
#endif
  end subroutine

  subroutine check_caches()
    SW_REAL :: one,near,copy
    integer(int64) :: bytes
    logical :: fits
    auxiliary_cache_limit_bytes = 1048576_int64
    one = SW_LITERAL(1.0_sw_kind,"1.0")
    near = one+epsilon(one)
    if(.not.sw_cache_equal(one,one).or..not.cache_real_equal(one,one)) error stop "equal cache keys"
#ifdef SWF_PRECISION_MPFR
    call mpfr_clear_caches()
    if(sw_cache_equal(one,near).or.cache_real_equal(one,near)) error stop "exact MPFR cache keys"
    if(sw_real_bytes()<8_int64*size(one%data,kind=int64)) error stop "MPFR payload accounting"
#else
    if(.not.sw_cache_equal(one,near).or..not.cache_real_equal(one,near)) error stop "native cache tolerance"
    if(sw_real_bytes()/=storage_size(one,kind=int64)/8) error stop "native element bytes"
#endif
    near = one+SW_LITERAL(0.01_sw_kind,"0.01")
    if(sw_cache_equal(one,near).or.cache_real_equal(one,near)) error stop "distinct cache keys"
    copy = one
    call sw_release(one)
#ifdef SWF_PRECISION_MPFR
    if(allocated(one%data)) error stop "released scalar storage"
#else
    if(one/=0) error stop "released native scalar"
#endif
    if(copy/=1) error stop "released scalar copy ownership"
    call check_workspace()
    call check_pleg()
    call check_qleg()
    call check_eviction()
    call check_gauss()
    bytes = 0_int64
    fits = add_cache_bytes(bytes,huge(0_int64),2)
    if(fits.or.bytes/=0) error stop "overflowing cache request"
    fits = add_cache_bytes(bytes,8_int64,-1)
    if(fits.or.bytes/=0) error stop "negative cache dimension"
    auxiliary_cache_limit_bytes = saved_budget
    call check_endpoint()
  end subroutine

  subroutine check_endpoint()
    SW_REAL :: c,x1,arguments(1),r1(3),r1d(3),r2(3),r2d(3),s1(3,1),s1d(3,1),eigenvalues(3)
    SW_REAL :: expected(3),ten,tolerance
    integer :: r1e(3),r1de(3),r2e(3),r2de(3),se(3,1),sde(3,1),radial_accuracy(3),angular_accuracy(3,1)
    c = SW_LITERAL(1.0e-40_sw_kind,"1.0e-40")
    x1 = 0
    arguments = 0
    ten = SW_LITERAL(10.0_sw_kind,"10.0")
    tolerance = 512*epsilon(c)
    call profcn(c,0,3,1,x1,0,0,1,arguments,r1,r1e,r1d,r1de,r2,r2e,r2d,r2de, &
      radial_accuracy,s1,se,s1d,sde,angular_accuracy,eigenvalues)
#ifdef SWF_PRECISION_MPFR
    if(mpfr_status/=0) error stop "radial endpoint status"
#endif
    r1 = r1*ten**r1e
    r1d = r1d*ten**r1de
#ifdef SWF_PRECISION_MPFR
    if(any(.not.mp_is_finite(r1)).or.any(.not.mp_is_finite(r1d)).or. &
      any(.not.mp_is_finite(eigenvalues))) error stop "finite radial endpoint"
#else
    if(any(.not.ieee_is_finite(r1)).or.any(.not.ieee_is_finite(r1d)).or. &
      any(.not.ieee_is_finite(eigenvalues))) error stop "finite radial endpoint"
#endif
    if(any(r1==0)) error stop "nonzero radial endpoint"
    expected = (eigenvalues-c*c)*r1/2
    if(any(abs(r1d-expected)>tolerance*max(abs(r1d),abs(expected)))) error stop "radial endpoint identity"
#ifdef SWF_PRECISION_MPFR
    call mpfr_clear_caches()
#endif
  end subroutine

  subroutine check_workspace()
    integer(int64) :: total
    logical :: fits
    total = 7_int64
    fits = sw_account_workspace(total,8_int64,2_int64,[2_int64,3_int64],200_int64)
    if(.not.fits.or.total/=103) error stop "workspace accounting"
    total = 7_int64
    fits = sw_account_workspace(total,8_int64,huge(0_int64),[1_int64],200_int64)
    if(fits.or.total/=7) error stop "workspace copy overflow"
    fits = sw_account_workspace(total,8_int64,2_int64,[huge(0_int64)],200_int64)
    if(fits.or.total/=7) error stop "workspace dimension overflow"
    fits = sw_account_workspace(total,8_int64,2_int64,[2_int64,-1_int64],200_int64)
    if(fits.or.total/=7) error stop "negative workspace dimension"
    fits = sw_account_workspace(total,8_int64,2_int64,[0_int64],200_int64)
    if(.not.fits.or.total/=7) error stop "empty workspace preservation"
    fits = sw_account_workspace(total,8_int64,25_int64,[0_int64],200_int64)
    if(fits.or.total/=7) error stop "empty workspace initial budget"
    total = 201_int64
    fits = sw_account_workspace(total,8_int64,2_int64,[1_int64],200_int64)
    if(fits.or.total/=201) error stop "invalid workspace total"
  end subroutine

  subroutine check_pleg()
    integer, parameter :: rows=2,capacity=8,retained=4
    SW_REAL :: arguments(rows),pr(rows,capacity),pdr(rows,capacity),coefficients(capacity,8)
    SW_REAL :: actual_arguments(rows),actual_pr(rows,capacity),actual_pdr(rows,capacity)
    SW_REAL :: actual_coefficients(capacity,8),pnorm(rows),pdnorm(rows),actual_pnorm(rows),actual_pdnorm(rows)
    integer :: ipnorm(rows),ipdnorm(rows),actual_ipnorm(rows),actual_ipdnorm(rows),i,j,slot,ndec,nex
    integer(int64) :: expected_bytes
    arguments = [SW_LITERAL(0.25_sw_kind,"0.25"),SW_LITERAL(0.5_sw_kind,"0.5")]
    ndec = precision(arguments(1))
    nex = range(arguments(1))-1
    do j=1,capacity
      do i=1,rows
        pr(i,j) = sw_real(i+2*j)/8
        pdr(i,j) = sw_real(i+3*j)/4
      end do
      do i=1,8
        coefficients(j,i) = sw_real(j+5*i)
      end do
    end do
    pnorm = 2
    pdnorm = 3
    ipnorm = 4
    ipdnorm = 5
    expected_bytes = int(rows+2*rows*retained+2*rows+8*retained,int64)*sw_real_bytes() &
      +2_int64*rows*(storage_size(1)/8)
    auxiliary_cache_limit_bytes = expected_bytes-1
    call store_pleg_cache(0,3,ndec,nex,arguments,rows,rows,retained,capacity,pr,pdr, &
      pdnorm,ipdnorm,pnorm,ipnorm,coefficients(:,1),coefficients(:,2),coefficients(:,3), &
      coefficients(:,4),coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
    if(any(pleg_cache%valid)) error stop "pleg over-budget insertion"
    auxiliary_cache_limit_bytes = expected_bytes
    call store_pleg_cache(0,3,ndec,nex,arguments,rows,rows,retained,capacity,pr,pdr, &
      pdnorm,ipdnorm,pnorm,ipnorm,coefficients(:,1),coefficients(:,2),coefficients(:,3), &
      coefficients(:,4),coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
    slot = find_pleg_cache(0,3,ndec,nex,arguments,rows,rows,retained,capacity)
    if(slot==0) error stop "pleg insertion"
    if(pleg_cache(slot)%bytes/=expected_bytes) error stop "pleg byte accounting"
    if(size(pleg_cache(slot)%pr,2)/=retained) error stop "pleg retained prefix"
    if(find_pleg_cache(0,3,ndec,nex,arguments,rows,rows,retained/2,capacity/2)/=slot) &
      error stop "pleg prefix lookup"
    if(find_pleg_cache(0,3,ndec,nex,arguments,rows,rows,retained+1,capacity)/=0) &
      error stop "pleg uncached extent"
    if(find_pleg_cache(0,2,ndec,nex,arguments,rows,rows,retained,capacity)/=0) error stop "pleg mode key"
    if(find_pleg_cache(0,3,ndec+1,nex,arguments,rows,rows,retained,capacity)/=0) error stop "pleg precision key"
    call load_pleg_cache(slot,rows,rows,retained,capacity,actual_arguments,actual_pr,actual_pdr, &
      actual_pdnorm,actual_ipdnorm,actual_pnorm,actual_ipnorm,actual_coefficients(:,1), &
      actual_coefficients(:,2),actual_coefficients(:,3),actual_coefficients(:,4), &
      actual_coefficients(:,5),actual_coefficients(:,6),actual_coefficients(:,7),actual_coefficients(:,8))
    if(any(actual_arguments/=arguments).or.any(actual_pr(:,1:retained)/=pr(:,1:retained)).or. &
      any(actual_pdr(:,1:retained)/=pdr(:,1:retained)).or. &
      any(actual_coefficients(1:retained,:)/=coefficients(1:retained,:))) error stop "pleg loaded values"
    if(any(actual_pnorm/=pnorm).or.any(actual_pdnorm/=pdnorm).or. &
      any(actual_ipnorm/=ipnorm).or.any(actual_ipdnorm/=ipdnorm)) error stop "pleg loaded normalization"
    pr = 0
    actual_pr = 0
    if(pleg_cache(slot)%pr(1,1)/=SW_LITERAL(0.375_sw_kind,"0.375")) error stop "pleg copy ownership"
    call clear_pleg_cache_slot(slot)
    if(pleg_cache(slot)%valid.or.pleg_cache(slot)%bytes/=0.or.allocated(pleg_cache(slot)%barg).or. &
      allocated(pleg_cache(slot)%pr).or.allocated(pleg_cache(slot)%pdr).or. &
      allocated(pleg_cache(slot)%pnorm).or.allocated(pleg_cache(slot)%ipnorm).or. &
      allocated(pleg_cache(slot)%alpha).or.allocated(pleg_cache(slot)%coefe)) error stop "pleg cleanup"
    auxiliary_cache_limit_bytes = 1048576_int64
  end subroutine

  subroutine check_qleg()
    integer, parameter :: degrees=2,capacity=4
    SW_REAL :: x,near,qdr(capacity),qr(capacity),qdl(degrees),ql(degrees),qdml,qml,termpq
    integer :: iqdl(degrees),iql(degrees),iqdml,iqml,itermpq,slot,i,ndec
    integer(int64) :: expected_bytes
    x = SW_LITERAL(1.0e-30_sw_kind,"1e-30")
    near = x*2
    ndec = precision(x)
    qdml = 7
    qml = 8
    termpq = 9
    do i=1,capacity
      qdr(i) = sw_real(i)/4
      qr(i) = sw_real(i)/8
    end do
    qdl = 3
    ql = 4
    iqdl = 5
    iql = 6
    expected_bytes = int(4+2*capacity+2*degrees,int64)*sw_real_bytes() &
      +2_int64*degrees*(storage_size(1)/8)
    auxiliary_cache_limit_bytes = expected_bytes
    call store_qleg_cache(0,degrees,capacity,capacity,x,ndec,qdr,qdml,10,qdl,iqdl,qr,qml,11,ql,iql,termpq,12)
    slot = find_qleg_cache(0,degrees,capacity,capacity,x,ndec)
    if(slot==0) error stop "qleg insertion"
    if(qleg_cache(slot)%bytes/=expected_bytes) error stop "qleg byte accounting"
    if(find_qleg_cache(0,degrees,capacity,capacity,near,ndec)/=0) error stop "qleg exact boundary key"
    if(find_qleg_cache(0,degrees,capacity,capacity,x,ndec+1)/=0) error stop "qleg precision key"
    qdr = 0
    qdml = 0
    call qleg_cached(0,degrees,capacity,capacity,x,ndec,qdr,qdml,iqdml,qdl,iqdl,qr,qml,iqml,ql,iql,termpq,itermpq)
    if(qdr(1)/=SW_LITERAL(0.25_sw_kind,"0.25").or.qdml/=7.or.qml/=8.or.termpq/=9) &
      error stop "qleg loaded scalars"
    if(any(qdl/=3).or.any(ql/=4).or.any(iqdl/=5).or.any(iql/=6).or. &
      iqdml/=10.or.iqml/=11.or.itermpq/=12) error stop "qleg loaded normalization"
    do i=1,capacity
      if(qdr(i)/=sw_real(i)/4.or.qr(i)/=sw_real(i)/8) error stop "qleg loaded arrays"
    end do
    qdml = 0
    if(qleg_cache(slot)%qdml/=7) error stop "qleg copy ownership"
    call clear_qleg_cache_slot(slot)
    if(qleg_cache(slot)%valid.or.qleg_cache(slot)%bytes/=0.or.allocated(qleg_cache(slot)%qdr).or. &
      allocated(qleg_cache(slot)%qr).or.allocated(qleg_cache(slot)%qdl).or.allocated(qleg_cache(slot)%ql).or. &
      allocated(qleg_cache(slot)%iqdl).or.allocated(qleg_cache(slot)%iql)) error stop "qleg cleanup"
#ifdef SWF_PRECISION_MPFR
    if(allocated(qleg_cache(slot)%x1%data).or.allocated(qleg_cache(slot)%qdml%data).or. &
      allocated(qleg_cache(slot)%qml%data).or.allocated(qleg_cache(slot)%termpq%data)) &
      error stop "qleg scalar storage cleanup"
#else
    if(qleg_cache(slot)%x1/=0.or.qleg_cache(slot)%qdml/=0.or. &
      qleg_cache(slot)%qml/=0.or.qleg_cache(slot)%termpq/=0) error stop "qleg scalar reset"
#endif
    auxiliary_cache_limit_bytes = 1048576_int64
  end subroutine

  subroutine check_eviction()
    SW_REAL :: x,values(1),scalar
    integer :: entry,slot,misses,ndec,exponents(1)
    integer(int64) :: entry_bytes
    values = 1
    scalar = 2
    exponents = 0
    ndec = precision(scalar)
    do slot=1,size(qleg_cache)
      call clear_qleg_cache_slot(slot)
    end do
    entry_bytes = 8_int64*sw_real_bytes()+2_int64*(storage_size(1)/8)
    auxiliary_cache_limit_bytes = (size(qleg_cache,kind=int64)+1)*entry_bytes
    do entry=1,size(qleg_cache)+1
      x = sw_real(entry)/16
      call store_qleg_cache(0,1,1,1,x,ndec,values,scalar,0,values,exponents, &
        values,scalar,0,values,exponents,scalar,0)
    end do
    if(count(qleg_cache%valid)/=size(qleg_cache)) error stop "full qleg cache"
    if(sum(qleg_cache%bytes)/=size(qleg_cache,kind=int64)*entry_bytes) error stop "qleg eviction accounting"
    if(find_qleg_cache(0,1,1,1,x,ndec)==0) error stop "newest qleg entry"
    misses = 0
    do entry=1,size(qleg_cache)
      x = sw_real(entry)/16
      if(find_qleg_cache(0,1,1,1,x,ndec)==0) misses = misses+1
    end do
    if(misses/=1) error stop "qleg single-slot eviction"
    do slot=1,size(qleg_cache)
      call clear_qleg_cache_slot(slot)
    end do
    if(any(qleg_cache%valid).or.sum(qleg_cache%bytes)/=0) error stop "qleg eviction cleanup"
    auxiliary_cache_limit_bytes = 1048576_int64
  end subroutine

  subroutine check_gauss()
    SW_REAL :: x(4),w(4),saved_x(4),saved_w(4),one,tolerance
    integer :: ndec
    logical :: fits
    one = SW_LITERAL(1.0_sw_kind,"1.0")
    tolerance = 256*epsilon(one)
    ndec = precision(one)
    call gauss_cached(ndec,4,x,w)
    if(.not.gauss_cache%valid.or.gauss_cache%bytes/=8*sw_real_bytes()) error stop "gauss insertion"
    if(abs(sum(w)-2)>tolerance.or.abs(sum(w*x**2)-2*one/3)>tolerance) error stop "gauss moments"
    saved_x = x
    saved_w = w
    call gauss_cached(ndec,4,x,w)
    if(any(x/=saved_x).or.any(w/=saved_w)) error stop "gauss cache reuse"
    call gauss_cached(ndec,2,x(1:2),w(1:2))
    if(gauss_cache%n/=2.or.size(gauss_cache%x)/=2.or.gauss_cache%bytes/=4*sw_real_bytes()) &
      error stop "gauss cache replacement"
    auxiliary_cache_limit_bytes = 0_int64
    fits = reserve_cache_bytes(0_int64)
    if(.not.fits.or.gauss_cache%valid.or.allocated(gauss_cache%x).or.allocated(gauss_cache%w).or. &
      any(pleg_cache%valid).or.any(qleg_cache%valid)) error stop "zero-budget cleanup"
    call gauss_cached(ndec,2,x(1:2),w(1:2))
    if(gauss_cache%valid.or.allocated(gauss_cache%x).or.allocated(gauss_cache%w)) &
      error stop "zero-budget cache suppression"
    if(abs(sum(w(1:2))-2)>tolerance) error stop "uncached gauss result"
    auxiliary_cache_limit_bytes = 1048576_int64
#ifdef SWF_PRECISION_MPFR
    call gauss_cached(ndec,2,x(1:2),w(1:2))
    call mpfr_clear_caches()
    if(gauss_cache%valid.or.allocated(gauss_cache%x).or.allocated(gauss_cache%w)) &
      error stop "batch cache cleanup"
#endif
  end subroutine
end program prolate_cache_test
