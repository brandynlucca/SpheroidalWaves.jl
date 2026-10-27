! Exercise cache misses, reuse, and replacement on independent OpenMP threads.
! Guard cells catch a node-count/tolerance swap even for small quadrature rules.
program native_cache
  use oblate_swf, only: oblfcn, gauss_cached, gauss_cache, pleg_cache, store_pleg_cache, &
    load_pleg_cache, find_pleg_cache, auxiliary_cache_limit_bytes, reserve_cache_bytes, add_cache_bytes
  use prolate_swf, only: profcn
  use complex_oblate_swf, only: coblfcn
  use complex_prolate_swf, only: cprofcn
  use oblate_parameters, only: knd
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: iteration, failures

  failures = 0
!$omp parallel do reduction(+:failures)
  do iteration = 1, 256
    failures = failures + check_rule(iteration)
  end do
!$omp end parallel do
  failures = failures + check_pleg_storage()
  failures = failures + check_native_modes()
  if (failures /= 0) then
    print *, 'Native cache and mode failures:', failures
    stop 1
  end if
  print *, 'Passed quadrature, Legendre prefix, cache memory, and mode checks'

contains
  integer function check_rule(iteration) result(errors)
    integer, intent(in) :: iteration
    integer, parameter :: counts(4) = [2, 8, 32, 96]
    real(knd), parameter :: guard = -12345.0_knd, tolerance = 2.0e-13_knd
    real(knd) :: x(0:128), w(0:128), saved_x(96), saved_w(96), moment
    integer :: n, degree, repeat

    n = counts(mod(iteration-1, 4)+1)
    errors = 0
    do repeat = 1, 2
      x = guard
      w = guard
      call gauss_cached(15, n, x(1:n), w(1:n))
      if (x(0) /= guard .or. w(0) /= guard) errors = errors + 1
      if (any(x(n+1:) /= guard) .or. any(w(n+1:) /= guard)) errors = errors + 1
      if (.not. all(ieee_is_finite(x(1:n))) .or. .not. all(ieee_is_finite(w(1:n)))) errors = errors + 1
      if (any(abs(x(1:n)) >= 1) .or. any(w(1:n) <= 0)) errors = errors + 1
      do degree = 0, min(2*n-1, 10)
        moment = 0
        if (mod(degree,2) == 0) moment = 2.0_knd/(degree+1)
        if (abs(sum(w(1:n)*x(1:n)**degree)-moment) > tolerance) errors = errors + 1
      end do
      if (repeat == 1) then
        saved_x(1:n) = x(1:n)
        saved_w(1:n) = w(1:n)
      else
        if (any(x(1:n) /= saved_x(1:n)) .or. any(w(1:n) /= saved_w(1:n))) errors = errors + 1
      end if
    end do
  end function check_rule

  integer function check_pleg_storage() result(errors)
    integer, parameter :: capacity=512, retained=16
    real(knd) :: arguments(1), pr(1,capacity), pdr(1,capacity), coefficients(capacity,8)
    real(knd) :: actual_pr(1,capacity), actual_pdr(1,capacity), actual_coefficients(capacity,8)
    real(knd) :: pnorm(1), pdnorm(1), actual_pnorm(1), actual_pdnorm(1), actual_arguments(1)
    integer :: ipnorm(1), ipdnorm(1), actual_ipnorm(1), actual_ipdnorm(1), index, column, slot
    integer(int64) :: saved_budget, expected_bytes, bytes
    logical :: fits
    errors=0
    saved_budget=auxiliary_cache_limit_bytes
    auxiliary_cache_limit_bytes=0_int64
    fits=reserve_cache_bytes(0_int64)
    if(.not.fits) errors=errors+1
    auxiliary_cache_limit_bytes=saved_budget
    arguments=0.25_knd
    pnorm=2.0_knd
    pdnorm=3.0_knd
    ipnorm=4
    ipdnorm=5
    do index=1,capacity
      pr(1,index)=real(index,knd)/8
      pdr(1,index)=real(index,knd)/4
      do column=1,8
        coefficients(index,column)=real(index+column,knd)
      end do
    end do
    call store_pleg_cache(0,3,15,300,arguments,1,1,retained,capacity,pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm, &
      coefficients(:,1),coefficients(:,2),coefficients(:,3),coefficients(:,4), &
      coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
    slot=find_pleg_cache(0,3,15,300,arguments,1,1,retained,capacity)
    if(slot==0) then
      errors=errors+1
    else
      expected_bytes=int(10*retained+3,int64)*(storage_size(1.0_knd)/8)+2*(storage_size(1)/8)
      if(pleg_cache(slot)%bytes/=expected_bytes.or.pleg_cache(slot)%maxp/=capacity) errors=errors+1
      if(size(pleg_cache(slot)%pr,2)/=retained.or.size(pleg_cache(slot)%pdr,2)/=retained.or. &
        size(pleg_cache(slot)%alpha)/=retained.or.size(pleg_cache(slot)%coefe)/=retained) errors=errors+1
      call load_pleg_cache(slot,1,1,retained,capacity,actual_arguments,actual_pr,actual_pdr, &
        actual_pdnorm,actual_ipdnorm,actual_pnorm,actual_ipnorm,actual_coefficients(:,1), &
        actual_coefficients(:,2),actual_coefficients(:,3),actual_coefficients(:,4), &
        actual_coefficients(:,5),actual_coefficients(:,6),actual_coefficients(:,7),actual_coefficients(:,8))
      if(any(actual_arguments/=arguments).or.any(actual_pr(:,1:retained)/=pr(:,1:retained)).or. &
        any(actual_pdr(:,1:retained)/=pdr(:,1:retained)).or. &
        any(actual_coefficients(1:retained,:)/=coefficients(1:retained,:)).or. &
        any(actual_pnorm/=pnorm).or.any(actual_pdnorm/=pdnorm).or. &
        any(actual_ipnorm/=ipnorm).or.any(actual_ipdnorm/=ipdnorm)) errors=errors+1
      if(find_pleg_cache(0,3,15,300,arguments,1,1,retained/2,capacity/2)==0) errors=errors+1
      if(find_pleg_cache(0,3,15,300,arguments,1,1,retained+1,capacity)/=0) errors=errors+1
      if(find_pleg_cache(0,3,15,300,arguments,1,1,retained,capacity+1)/=0) errors=errors+1
      auxiliary_cache_limit_bytes=expected_bytes
      arguments=0.5_knd
      call store_pleg_cache(0,3,15,300,arguments,1,1,retained,capacity,pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm, &
      coefficients(:,1),coefficients(:,2),coefficients(:,3),coefficients(:,4), &
      coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
      if(sum(pleg_cache%bytes)+gauss_cache%bytes>auxiliary_cache_limit_bytes) errors=errors+1
      if(find_pleg_cache(0,3,15,300,arguments,1,1,retained,capacity)==0) errors=errors+1
      arguments=0.25_knd
      if(find_pleg_cache(0,3,15,300,arguments,1,1,retained,capacity)/=0) errors=errors+1
    end if
    auxiliary_cache_limit_bytes=0_int64
    fits=reserve_cache_bytes(0_int64)
    call store_pleg_cache(0,3,15,300,arguments,1,1,retained,capacity,pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm, &
      coefficients(:,1),coefficients(:,2),coefficients(:,3),coefficients(:,4), &
      coefficients(:,5),coefficients(:,6),coefficients(:,7),coefficients(:,8))
    if(.not.fits.or.any(pleg_cache%valid)) errors=errors+1
    errors=errors+check_rule(2)
    if(gauss_cache%valid.or.allocated(gauss_cache%x).or.allocated(gauss_cache%w)) errors=errors+1
    auxiliary_cache_limit_bytes=saved_budget
    bytes=0_int64
    if(add_cache_bytes(bytes,huge(0_int64),2)) errors=errors+1
    if(add_cache_bytes(bytes,8_int64,huge(0),huge(0))) errors=errors+1
  end function check_pleg_storage

  integer function check_native_modes() result(errors)
    integer, parameter :: degrees=2, radial_modes(5)=[2,0,0,1,2], angular_modes(5)=[2,0,2,2,2]
    real(knd) :: c,x,arguments(1),r1(degrees),r1d(degrees),r2(degrees),r2d(degrees)
    real(knd) :: s1(degrees,1),s1d(degrees,1),eigenvalues(degrees)
    complex(knd) :: cc,cr1(degrees),cr1d(degrees),cr2(degrees),cr2d(degrees)
    complex(knd) :: cs1(degrees,1),cs1d(degrees,1),ceigenvalues(degrees)
    complex(knd) :: values(degrees,7),reference(degrees,7)
    integer :: r1e(degrees),r1de(degrees),r2e(degrees),r2de(degrees),se(degrees,1),sde(degrees,1)
    integer :: radial_accuracy(degrees),angular_accuracy(degrees,1),angular_derivative_accuracy(degrees,1)
    integer :: geometry,stage,radial,angular,channel
    logical :: active(7)
    errors=0
    c=1.0_knd
    cc=cmplx(1.0_knd,0.25_knd,knd)
    arguments=0.25_knd
    do geometry=0,3
      x=1.0_knd
      if(mod(geometry,2)==1) x=2.0_knd
      do stage=1,5
        radial=radial_modes(stage)
        angular=angular_modes(stage)
        select case(geometry)
        case(0)
          call profcn(c,0,degrees,radial,x,angular,0,1,arguments,r1,r1e,r1d,r1de,r2,r2e,r2d,r2de, &
            radial_accuracy,s1,se,s1d,sde,angular_accuracy,eigenvalues)
        case(1)
          call oblfcn(c,0,degrees,radial,x,angular,0,1,arguments,r1,r1e,r1d,r1de,r2,r2e,r2d,r2de, &
            radial_accuracy,s1,se,s1d,sde,angular_accuracy,eigenvalues)
        case(2)
          call cprofcn(cc,0,degrees,radial,x,angular,0,1,arguments,cr1,r1e,cr1d,r1de,cr2,r2e,cr2d,r2de, &
            radial_accuracy,cs1,se,cs1d,sde,angular_accuracy,angular_derivative_accuracy,ceigenvalues)
        case(3)
          call coblfcn(cc,0,degrees,radial,x,angular,0,1,arguments,cr1,r1e,cr1d,r1de,cr2,r2e,cr2d,r2de, &
            radial_accuracy,cs1,se,cs1d,sde,angular_accuracy,angular_derivative_accuracy,ceigenvalues)
        end select
        active=[.true.,radial>0,radial>0,radial==2,radial==2,angular>0,angular>0]
        if(geometry<2) then
          values(:,1)=eigenvalues
          if(radial>0) then
            values(:,2)=r1*10.0_knd**r1e
            values(:,3)=r1d*10.0_knd**r1de
          end if
          if(radial==2) then
            values(:,4)=r2*10.0_knd**r2e
            values(:,5)=r2d*10.0_knd**r2de
          end if
          if(angular>0) then
            values(:,6)=s1(:,1)*10.0_knd**se(:,1)
            values(:,7)=s1d(:,1)*10.0_knd**sde(:,1)
          end if
        else
          values(:,1)=ceigenvalues
          if(radial>0) then
            values(:,2)=cr1*10.0_knd**r1e
            values(:,3)=cr1d*10.0_knd**r1de
          end if
          if(radial==2) then
            values(:,4)=cr2*10.0_knd**r2e
            values(:,5)=cr2d*10.0_knd**r2de
          end if
          if(angular>0) then
            values(:,6)=cs1(:,1)*10.0_knd**se(:,1)
            values(:,7)=cs1d(:,1)*10.0_knd**sde(:,1)
          end if
        end if
        if(stage==1) reference=values
        do channel=1,7
          if(.not.active(channel)) cycle
          if(.not.all(ieee_is_finite(real(values(:,channel)))).or. &
            .not.all(ieee_is_finite(aimag(values(:,channel))))) errors=errors+1
          if(any(abs(values(:,channel)-reference(:,channel))> &
            5.0e-12_knd*max(1.0_knd,abs(reference(:,channel))))) errors=errors+1
        end do
      end do
    end do
  end function check_native_modes
end program native_cache
