#include "precision.inc"

program complex_oblate_cache_test
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use SW_PRECISION_MODULE
  use SW_COBLATE_MODULE, only: auxiliary_cache_limit_bytes, cache_real_equal, &
    pleg_cache, gauss_cache, clear_pleg_cache_slot, &
    find_pleg_cache, store_pleg_cache, load_pleg_cache, &
    gauss_cached, reserve_cache_bytes, add_cache_bytes, coblfcn,pleg,pleg_cached,r1eta,r2eta
#ifdef SWF_PRECISION_MPFR
  use SW_COBLATE_MODULE, only: mpfr_clear_caches
#endif
  implicit none
  integer, parameter :: pleg_rows=2,pleg_capacity=16
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
  mpfr_bits = 113
  mpfr_rounding = 0
  call check_waves()
  mpfr_bits = saved_bits
  mpfr_rounding = saved_rounding
  mpfr_status = saved_status
#else
  call check_caches()
  call check_waves()
#endif
  call check_radial_fallback()
  auxiliary_cache_limit_bytes = saved_budget
  print *, "Complex oblate caches, workspace, angular and radial identities passed"
contains
  subroutine check_radial_fallback()
    SW_REAL :: x
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
    call clear_caches()
    auxiliary_cache_limit_bytes = 67108864_int64
#ifdef SWF_PRECISION_MPFR
    mpfr_bits = 113
    mpfr_rounding = 0
    mpfr_status = 0
    mpfr_workspace_limit_bytes = 1073741824_int64
#endif
    x = SW_LITERAL(0.2_sw_kind,"0.2")
    call evaluate_radial_fallback(500,50,2,x)
#ifdef SWF_PRECISION_MPFR
    call clear_caches()
    x = SW_LITERAL(0.25_sw_kind,"0.25")
    call evaluate_radial_fallback(500,0,310,x)
#else
#ifndef SWF_PRECISION_QUAD
    call clear_caches()
    call evaluate_radial_fallback(500,0,310,x)
#endif
#endif
    call check_radial_series(0)
    call check_radial_series(1)
    call clear_caches()
    auxiliary_cache_limit_bytes = cache_budget
#ifdef SWF_PRECISION_MPFR
    mpfr_bits = caller_bits
    mpfr_rounding = caller_rounding
    mpfr_status = caller_status
    mpfr_workspace_limit_bytes = workspace_budget
#endif
  end subroutine

  subroutine evaluate_radial_fallback(parameter_value,order,degrees,x)
    integer, intent(in) :: parameter_value,order,degrees
    SW_REAL, intent(in) :: x
    SW_REAL :: arguments(1),ten,tolerance
    SW_COMPLEX :: c,r1(degrees),r1d(degrees),r2(degrees),r2d(degrees)
    SW_COMPLEX :: s(degrees,1),sd(degrees,1),eigenvalues(degrees)
    SW_COMPLEX :: values(degrees,4),reference(degrees,4),wronskian(degrees)
    integer :: stage,r1e(degrees),r1de(degrees),r2e(degrees),r2de(degrees)
    integer :: se(degrees,1),sde(degrees,1),radial_accuracy(degrees),angular_accuracy(degrees,1),derivative_accuracy(degrees,1)
    c = sw_complex(sw_real(parameter_value),SW_LITERAL(0.25_sw_kind,"0.25"))
    arguments = 0
    ten = SW_LITERAL(10.0_sw_kind,"10.0")
    tolerance = SW_LITERAL(5.0e-8_sw_kind,"5e-8")
    do stage=1,2
      call coblfcn(c,order,degrees,2,x,0,0,1,arguments,r1,r1e,r1d,r1de,r2,r2e,r2d,r2de, &
        radial_accuracy,s,se,sd,sde,angular_accuracy,derivative_accuracy,eigenvalues)
#ifdef SWF_PRECISION_MPFR
      if(mpfr_status/=0) error stop "oblate radial fallback status"
#endif
      values(:,1) = r1*ten**r1e
      values(:,2) = r1d*ten**r1de
      values(:,3) = r2*ten**r2e
      values(:,4) = r2d*ten**r2de
      wronskian = c*(x*x+1)*(values(:,1)*values(:,4)-values(:,2)*values(:,3))
      if(.not.all(finite_complex(values)).or..not.all(finite_complex(wronskian))) &
        error stop "finite complex oblate radial fallback"
      if(any(abs(wronskian-1)>tolerance)) error stop "oblate radial fallback Wronskian"
      if(stage==1) then
        reference = values
      else
        if(any(abs(values-reference)>tolerance*max(sw_real(1),abs(reference)))) &
          error stop "oblate radial fallback cache consistency"
      end if
    end do
  end subroutine

  subroutine check_radial_series(terms)
    integer, intent(in) :: terms
    integer, parameter :: maxd=4,maxlp=2,maxn=8,maxp=10
    SW_REAL :: one,eta,wm,root7,pratb(maxp),pratt(maxp),pdratt(maxp)
    SW_COMPLEX :: cc,unit,enr(maxd),sneuf(maxn),sneudf(maxn),sneun(maxn),sneudr(maxlp)
    SW_COMPLEX :: r1,r2,r2d,actual(3),expected(3)
    integer :: ineue(maxn),r2e,r2de,naccnmax,nacceta,jeta,naccd,ndec,naccrpl,naccs1,naccs2
    one = 1
    cc = sw_complex(one,sw_real(0))
    unit = cc
    eta = one/2
    wm = 3*one/4
    root7 = sqrt(7*one)
    enr = 0
    sneuf = 1
    sneudf = 1
    sneun = 1
    sneudr = 1
    pratb = 1
    pratt = 1
    pdratt = 1
    ineue = 0
    ndec = precision(one)
    naccnmax = 0
    ! Synthetic unit ratios with a zero tail give R2=1 and R2'=2/sqrt(7).
    expected = [unit,2*unit/root7,unit]
    r1 = 3*root7/4
    if(terms==1) then
      ! One nonzero coefficient ratio gives R2=7/9 and R2'=95/(63*sqrt(7)).
      enr(1) = one/8
      expected(1) = 7*one/9
      expected(2) = 95*one/(63*root7)
      r1 = 161*root7/190
    end if
    call r2eta(0,0,cc,one,eta,1,8,ndec,maxd,maxlp,maxn,maxp,8,wm,enr,sneuf,sneun,ineue, &
      sneudf,sneudr,pdratt,pratb,pratt,one,0,one,0,r1,0,unit,0,ndec,naccrpl,naccnmax,ndec,r2,r2e,r2d,r2de,nacceta,jeta,naccd)
#ifdef SWF_PRECISION_MPFR
    if(mpfr_status/=0) error stop "r2eta direct arithmetic status"
#endif
    r2 = r2*sw_real(10)**r2e
    r2d = r2d*sw_real(10)**r2de
    ! The chosen R1 and R1'=1 make the normalized Wronskian equal one.
    actual = [r2,r2d,2*(r1*r2d-r2)]
    if(.not.all(finite_complex(actual))) error stop "nonfinite r2eta direct output"
    if(any(abs(actual-expected)>128*epsilon(one)*max(one,abs(expected)))) error stop "r2eta direct value"
    if(jeta/=terms) error stop "r2eta saved partial sum index"
    if(naccd/=ndec-2) error stop "r2eta denominator accuracy"
    if(nacceta<ndec-3.or.nacceta>ndec-1) error stop "r2eta Wronskian accuracy"
    call r1eta(0,0,cc,one,eta,1,8,ndec,range(one)-1,maxd,maxlp,maxn,maxp,8,wm, &
      enr,sneuf,sneun(1:maxlp),ineue(1:maxlp),sneudf,sneudr,pdratt,pratb,pratt, &
      one,0,one,0,0,r2,r2e,r2d,r2de,naccs1,naccs2,jeta)
    call require_status()
    actual(1) = r2*sw_real(10)**r2e
    actual(2) = r2d*sw_real(10)**r2de
    if(.not.all(finite_complex(actual(1:2)))) error stop "nonfinite r1eta direct output"
    if(any(abs(actual(1:2)-expected(1:2))>128*epsilon(one)*max(one,abs(expected(1:2))))) &
      error stop "r1eta direct value"
  end subroutine

  subroutine clear_caches()
#ifdef SWF_PRECISION_MPFR
    call mpfr_clear_caches()
#else
    integer(int64) :: budget
    logical :: fits
    budget = auxiliary_cache_limit_bytes
    auxiliary_cache_limit_bytes = 0_int64
    fits = reserve_cache_bytes(0_int64)
    auxiliary_cache_limit_bytes = budget
    if(.not.fits) error stop "complex oblate cache cleanup"
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
    call check_pleg_modes()
    call check_eviction()
    call check_gauss()
    bytes = 0_int64
    fits = add_cache_bytes(bytes,huge(0_int64),2)
    if(fits.or.bytes/=0) error stop "overflowing cache request"
    fits = add_cache_bytes(bytes,8_int64,-1)
    if(fits.or.bytes/=0) error stop "negative cache dimension"
    auxiliary_cache_limit_bytes = saved_budget
  end subroutine

  subroutine check_pleg_modes()
    SW_REAL :: arguments(pleg_rows),pr(pleg_rows,pleg_capacity),pdr(pleg_rows,pleg_capacity)
    SW_REAL :: coefficients(pleg_capacity,8),pnorm(pleg_rows),pdnorm(pleg_rows)
    SW_REAL :: saved_pr(pleg_rows,pleg_capacity),saved_pnorm(pleg_rows),saved_pdnorm(pleg_rows)
    integer :: ipnorm(pleg_rows),ipdnorm(pleg_rows),saved_ipnorm(pleg_rows),saved_ipdnorm(pleg_rows)
    integer :: mode,limcsav,limit,stage,slot
    do mode=0,2,2
      call clear_caches()
      arguments = [SW_LITERAL(0.25_sw_kind,"0.25"),SW_LITERAL(0.5_sw_kind,"0.5")]
      limcsav = 0
      limit = 8
      coefficients = 37
      pr = 37
      pdr = 37
      call evaluate_pleg(.true.,mode,limit,limcsav,arguments,pr,pdr,pnorm,pdnorm,ipnorm,ipdnorm,coefficients)
      if(any(pdr(:,1:limit+2)/=0)) error stop "omitted Legendre derivatives"
      if(any(pdnorm/=3).or.any(ipdnorm/=0)) error stop "Legendre derivative normalization"
      slot = find_pleg_cache(2,mode,precision(arguments(1)),range(arguments(1))-1, &
        arguments,pleg_rows,pleg_rows,limit,pleg_capacity)
      if(slot==0) error stop "Legendre mode insertion"
      if(any(pleg_cache(slot)%pdr/=0)) error stop "cached omitted Legendre derivatives"
      saved_pr(:,1:limit+2) = pr(:,1:limit+2)
      coefficients = 37
      pr = 37
      pdr = 37
      call evaluate_pleg(.true.,mode,limit,limcsav,arguments,pr,pdr,pnorm,pdnorm,ipnorm,ipdnorm,coefficients)
      if(any(pdr(:,1:limit)/=0).or.any(pr(:,1:limit)/=saved_pr(:,1:limit))) &
        error stop "Legendre mode cache hit"
      limit = 10
      arguments = [SW_LITERAL(0.375_sw_kind,"0.375"),SW_LITERAL(0.625_sw_kind,"0.625")]
      do stage=1,2
        call evaluate_pleg(.true.,mode,limit,limcsav,arguments,pr,pdr,pnorm,pdnorm,ipnorm,ipdnorm,coefficients)
        if(any(pdr(:,1:limit+2)/=0)) error stop "expanded omitted Legendre derivatives"
        if(stage==1) then
          saved_pr(:,1:limit+2) = pr(:,1:limit+2)
          saved_pnorm = pnorm
          saved_pdnorm = pdnorm
          saved_ipnorm = ipnorm
          saved_ipdnorm = ipdnorm
          call clear_caches()
          limcsav = 0
        else
          if(any(pr(:,1:limit+2)/=saved_pr(:,1:limit+2)).or.any(pnorm/=saved_pnorm).or. &
            any(pdnorm/=saved_pdnorm).or.any(ipnorm/=saved_ipnorm).or.any(ipdnorm/=saved_ipdnorm)) &
            error stop "Legendre cache lookahead and expanded reuse"
        end if
      end do
      arguments = [SW_LITERAL(0.4375_sw_kind,"0.4375"),SW_LITERAL(0.6875_sw_kind,"0.6875")]
      limcsav = limit
      call evaluate_pleg(.false.,mode,limit,limcsav,arguments,pr,pdr,pnorm,pdnorm,ipnorm,ipdnorm,coefficients)
      saved_pr(:,1:limit+2) = pr(:,1:limit+2)
      saved_pnorm = pnorm
      saved_pdnorm = pdnorm
      saved_ipnorm = ipnorm
      saved_ipdnorm = ipdnorm
      limcsav = 0
      call evaluate_pleg(.false.,mode,limit,limcsav,arguments,pr,pdr,pnorm,pdnorm,ipnorm,ipdnorm,coefficients)
      if(any(pr(:,1:limit+2)/=saved_pr(:,1:limit+2)).or.any(pnorm/=saved_pnorm).or. &
        any(pdnorm/=saved_pdnorm).or.any(ipnorm/=saved_ipnorm).or.any(ipdnorm/=saved_ipdnorm)) &
        error stop "Legendre reused-coefficient normalization"
    end do
    call clear_caches()
  end subroutine

  subroutine evaluate_pleg(cached,mode,limit,limcsav,arguments,pr,pdr,pnorm,pdnorm,ipnorm,ipdnorm,coefficients)
    logical, intent(in) :: cached
    integer, intent(in) :: mode,limit
    integer, intent(inout) :: limcsav,ipnorm(pleg_rows),ipdnorm(pleg_rows)
    SW_REAL, intent(inout) :: arguments(pleg_rows),pr(pleg_rows,pleg_capacity),pdr(pleg_rows,pleg_capacity)
    SW_REAL, intent(inout) :: pnorm(pleg_rows),pdnorm(pleg_rows),coefficients(pleg_capacity,8)
    integer :: ndec,nex
    ndec = precision(arguments(1))
    nex = range(arguments(1))-1
    if(cached) then
      call pleg_cached(2,limit,pleg_capacity,limcsav,mode,ndec,nex,arguments,pleg_rows,pleg_rows,pr,pdr, &
        pdnorm,ipdnorm,pnorm,ipnorm,coefficients(:,1),coefficients(:,2),coefficients(:,3), &
        coefficients(:,4),coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
    else
      call pleg(2,limit,pleg_capacity,limcsav,mode,ndec,nex,arguments,pleg_rows,pleg_rows,pr,pdr, &
        pdnorm,ipdnorm,pnorm,ipnorm,coefficients(:,1),coefficients(:,2),coefficients(:,3), &
        coefficients(:,4),coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
    end if
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

  subroutine check_eviction()
    SW_REAL :: arguments(1),values(1,1),normalization(1)
    integer :: entry,slot,misses,ndec,nex,exponents(1)
    integer(int64) :: entry_bytes
    call clear_caches()
    values = 1
    normalization = 2
    exponents = 0
    ndec = precision(normalization(1))
    nex = range(normalization(1))-1
    entry_bytes = 13_int64*sw_real_bytes()+2_int64*(storage_size(1)/8)
    auxiliary_cache_limit_bytes = (size(pleg_cache,kind=int64)+1)*entry_bytes
    do entry=1,size(pleg_cache)+1
      arguments = sw_real(entry)/64
      call insert_entry(arguments,values,normalization,exponents,ndec,nex)
    end do
    if(count(pleg_cache%valid)/=size(pleg_cache)) error stop "full pleg cache"
    if(sum(pleg_cache%bytes)/=size(pleg_cache,kind=int64)*entry_bytes) error stop "pleg eviction accounting"
    slot = find_pleg_cache(0,3,ndec,nex,arguments,1,1,1,1)
    if(slot==0) error stop "newest pleg entry"
    call check_minimal_load(slot)
    misses = 0
    do entry=1,size(pleg_cache)
      arguments = sw_real(entry)/64
      if(find_pleg_cache(0,3,ndec,nex,arguments,1,1,1,1)==0) misses = misses+1
    end do
    if(misses/=1) error stop "pleg single-slot eviction"
    auxiliary_cache_limit_bytes = entry_bytes
    arguments = SW_LITERAL(-0.5_sw_kind,"-0.5")
    call insert_entry(arguments,values,normalization,exponents,ndec,nex)
    if(count(pleg_cache%valid)/=1.or.sum(pleg_cache%bytes)/=entry_bytes) error stop "pleg budget eviction"
    if(find_pleg_cache(0,3,ndec,nex,arguments,1,1,1,1)==0) error stop "pleg budget replacement"
    do slot=1,size(pleg_cache)
      call clear_pleg_cache_slot(slot)
    end do
    if(any(pleg_cache%valid).or.sum(pleg_cache%bytes)/=0) error stop "pleg eviction cleanup"
    auxiliary_cache_limit_bytes = 1048576_int64
  end subroutine

  subroutine check_minimal_load(slot)
    integer, intent(in) :: slot
    SW_REAL :: arguments(1),pr(1,1),pdr(1,1),pnorm(1),pdnorm(1),coefficients(1,8)
    integer :: ipnorm(1),ipdnorm(1)
    call load_pleg_cache(slot,1,1,1,1,arguments,pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm, &
      coefficients(:,1),coefficients(:,2),coefficients(:,3),coefficients(:,4), &
      coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
    if(any(pr/=1).or.any(pdr/=1).or.any(coefficients/=1)) error stop "minimal cache capacity"
  end subroutine

  subroutine insert_entry(arguments,values,normalization,exponents,ndec,nex)
    SW_REAL, intent(in) :: arguments(1),values(1,1),normalization(1)
    integer, intent(in) :: exponents(1),ndec,nex
    call store_pleg_cache(0,3,ndec,nex,arguments,1,1,1,1,values,values,normalization,exponents, &
      normalization,exponents,values(1,:),values(1,:),values(1,:),values(1,:), &
      values(1,:),values(1,:),values(1,:),values(1,:))
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
      any(pleg_cache%valid)) error stop "zero-budget cleanup"
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
  subroutine check_waves()
    integer :: degrees
    do degrees=1,5
      call evaluate_waves(0,degrees,1,1,4,0)
      call evaluate_waves(1,degrees,5,2,1,2)
    end do
    call evaluate_waves(0,4,2,25,1,2)
    call evaluate_waves(0,3,10,1,4,2)
    call clear_caches()
  end subroutine

  subroutine evaluate_waves(order,degrees,cr,ci,ci_den,x_integer)
    integer, intent(in) :: order,degrees,cr,ci,ci_den,x_integer
    integer, parameter :: points=3
    SW_REAL :: x,arguments(points),ten,tolerance,one
    SW_COMPLEX :: c,r1(degrees),r1d(degrees),r2(degrees),r2d(degrees),eigenvalues(degrees)
    SW_COMPLEX :: s(degrees,points),sd(degrees,points),saved_eigenvalues(degrees)
    SW_COMPLEX :: saved_s(degrees,points),saved_sd(degrees,points),values(degrees,4),reference(degrees,4)
    integer :: r1e(degrees),r1de(degrees),r2e(degrees),r2de(degrees),se(degrees,points),sde(degrees,points)
    integer :: radial_accuracy(degrees),angular_accuracy(degrees,points),derivative_accuracy(degrees,points)
    integer :: angular,normalization,stage,degree,radial,saved_se(degrees,points),saved_sde(degrees,points)
    logical :: have_eigenvalues
    one = 1
    c = sw_complex(sw_real(cr),sw_real(ci)/ci_den)
    x = sw_real(x_integer)
    ten = 10
    arguments = [-one/4,sw_real(0),one/4]
    tolerance = 1024*epsilon(one)
    have_eigenvalues = .false.
    call clear_caches()
    do angular=0,2
      if(ci>20.and.angular/=2) cycle
      do normalization=0,min(1,angular)
        if(ci>20.and.normalization/=0) cycle
        do stage=1,2
          radial = 0
          call coblfcn(c,order,degrees,radial,x,angular,normalization,points,arguments, &
            r1,r1e,r1d,r1de,r2,r2e,r2d,r2de,radial_accuracy,s,se,sd,sde, &
            angular_accuracy,derivative_accuracy,eigenvalues)
          call require_status()
          if(.not.all(finite_complex(eigenvalues))) error stop "finite complex eigenvalues"
          if(have_eigenvalues) then
            if(any(abs(eigenvalues-saved_eigenvalues)>tolerance*max(one,abs(saved_eigenvalues)))) &
              error stop "complex eigenvalues independent of requested outputs"
          else
            saved_eigenvalues = eigenvalues
            have_eigenvalues = .true.
          end if
          if(angular==0) then
#ifdef SWF_PRECISION_MPFR
            if(allocated(s(1,1)%re%data).or.allocated(sd(1,1)%re%data)) error stop "eigenvalue-only angular fallthrough"
#endif
            cycle
          end if
          if(.not.all(finite_complex(s))) error stop "finite complex angular values"
          if(angular==1) then
            if(any(sd/=0).or.any(sde/=0).or.any(derivative_accuracy/=0)) &
              error stop "omitted complex angular derivatives"
          end if
          if(stage==1) then
            saved_s = s
            saved_sd = sd
            saved_se = se
            saved_sde = sde
          else
            if(any(s/=saved_s).or.any(sd/=saved_sd).or.any(se/=saved_se).or.any(sde/=saved_sde)) &
              error stop "complex angular cache consistency"
          end if
          s = s*ten**se
          sd = sd*ten**sde
          do degree=1,degrees
            if(abs(s(degree,1)-(-1)**(degree-1)*s(degree,3))>tolerance*max(one,abs(s(degree,3)))) &
              error stop "complex angular parity"
            if(angular==2) then
              if(.not.all(finite_complex(sd))) error stop "finite complex angular derivatives"
              if(abs(sd(degree,1)-(-1)**degree*sd(degree,3))>tolerance*max(one,abs(sd(degree,3)))) &
                error stop "complex angular derivative parity"
            end if
          end do
        end do
      end do
    end do
    radial = 2
    do stage=1,2
      radial_accuracy = -huge(0)
      call coblfcn(c,order,degrees,radial,x,0,0,points,arguments,r1,r1e,r1d,r1de,r2,r2e,r2d,r2de, &
        radial_accuracy,s,se,sd,sde,angular_accuracy,derivative_accuracy,eigenvalues)
      call require_status()
      if(any(radial_accuracy < -1).or.any(radial_accuracy>precision(one))) &
        error stop "complex radial accuracy initialization"
      values(:,1) = r1*ten**r1e
      values(:,2) = r1d*ten**r1de
      if(x_integer==0) then
        do degree=1,degrees
          if(mod(degree-1,2)==0) then
            if(values(degree,2)/=0) error stop "even complex origin derivative"
          else
            if(values(degree,1)/=0) error stop "odd complex origin value"
          end if
        end do
      end if
      values(:,3) = r2*ten**r2e
      values(:,4) = r2d*ten**r2de
      if(.not.all(finite_complex(values))) error stop "finite complex radial waves"
      ! The large-imaginary fixture loses more than 44 digits in the Wronskian subtraction.
      if(ci<=20) then
        if(any(abs(c*(x*x+1)*(values(:,1)*values(:,4)-values(:,2)*values(:,3))-1) &
          >SW_LITERAL(1.0e-8_sw_kind,"1e-8"))) error stop "complex radial Wronskian"
      end if
      if(stage==1) then
        reference = values
      else
        if(any(abs(values-reference)>tolerance*max(one,abs(reference)))) error stop "complex radial cache consistency"
      end if
    end do
  end subroutine

  impure elemental logical function finite_complex(value) result(finite)
    SW_COMPLEX, intent(in) :: value
#ifdef SWF_PRECISION_MPFR
    finite = mp_is_finite(value%re).and.mp_is_finite(value%im)
#else
    finite = ieee_is_finite(real(value,sw_kind)).and.ieee_is_finite(aimag(value))
#endif
  end function

  subroutine require_status()
#ifdef SWF_PRECISION_MPFR
    if(mpfr_status/=0) error stop "complex oblate arithmetic status"
#endif
  end subroutine

end program complex_oblate_cache_test
