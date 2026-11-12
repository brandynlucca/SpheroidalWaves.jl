program mpfr_batch_context
  use, intrinsic :: iso_c_binding, only: c_char,c_null_char,c_int,c_int64_t
  use spheroidal_mpfr_real
  use spheroidal_mpfr_batch_api, only: spheroidal_mpfr_batch
  use prolate_swf_mpfr, only: p_pleg=>pleg_cache,p_qleg=>qleg_cache,p_gauss=>gauss_cache,p_clear=>mpfr_clear_caches
  use oblate_swf_mpfr, only: o_pleg=>pleg_cache,o_gauss=>gauss_cache,o_clear=>mpfr_clear_caches
  use complex_prolate_swf_mpfr, only: cp_pleg=>pleg_cache,cp_qleg=>qleg_cache,cp_gauss=>gauss_cache,cp_clear=>mpfr_clear_caches
  use complex_oblate_swf_mpfr, only: co_pleg=>pleg_cache,co_gauss=>gauss_cache,co_clear=>mpfr_clear_caches
  use omp_lib
  implicit none
  integer, parameter :: workers=4,width=192,points=2,cases=16,requested_bits(workers)=[113,160,192,256]
  integer(c_int64_t), parameter :: normal_budget=268435456_c_int64_t
  character(kind=c_char) :: reference(width,10,points,cases,workers)
  integer(c_int) :: reference_accuracy(points,cases,workers)
  integer :: worker,fixture,failures,actual_workers
  logical :: saved_dynamic

  failures=0
  do worker=1,workers
    do fixture=1,cases
      call evaluate(worker,fixture,0,reference(:,:,:,fixture,worker),reference_accuracy(:,fixture,worker),failures)
    end do
  end do
  if(failures/=0) error stop 'Serial mixed-geometry MPFR references failed'

  saved_dynamic=omp_get_dynamic()
  call omp_set_dynamic(.false.)
  actual_workers=0
!$omp parallel num_threads(workers) default(none) shared(reference,reference_accuracy,actual_workers) &
!$omp private(worker) reduction(+:failures)
!$omp single
  actual_workers=omp_get_num_threads()
!$omp end single
  worker=omp_get_thread_num()+1
  call check_worker(worker,reference(:,:,:,:,worker),reference_accuracy(:,:,worker),failures)
!$omp end parallel
  call omp_set_dynamic(saved_dynamic)
  if(actual_workers/=workers) error stop 'Mixed-geometry MPFR test requires four OpenMP workers'
  if(failures/=0) then
    print *, 'Mixed-geometry MPFR failures:',failures
    error stop 1
  end if
  print *, 'Mixed MPFR geometries, precisions, rounding modes, cache cleanup, and failure recovery passed'

contains
  recursive subroutine check_worker(worker,reference,reference_accuracy,failures)
    integer, intent(in) :: worker
    character(kind=c_char), intent(in) :: reference(width,10,points,cases)
    integer(c_int), intent(in) :: reference_accuracy(points,cases)
    integer, intent(inout) :: failures
    character(kind=c_char) :: output(width,10,points)
    integer(c_int) :: accuracy(points)
    integer :: iteration,step,fixture,failure

    do iteration=1,2
      do step=1,cases
        fixture=1+mod(step-1+4*(worker-1)+iteration,cases)
!$omp barrier
        call evaluate(worker,fixture,0,output,accuracy,failures)
        call compare(worker,fixture,output,accuracy,reference(:,:,:,fixture),reference_accuracy(:,fixture),failures)
      end do
    end do
    do failure=1,11
      fixture=1+mod(failure-1+4*(worker-1),cases)
!$omp barrier
      call evaluate(worker,fixture,failure,output,accuracy,failures)
      call evaluate(worker,fixture,0,output,accuracy,failures)
      call compare(worker,fixture,output,accuracy,reference(:,:,:,fixture),reference_accuracy(:,fixture),failures)
    end do
  end subroutine

  recursive subroutine compare(worker,fixture,output,accuracy,reference,reference_accuracy,failures)
    integer, intent(in) :: worker,fixture
    character(kind=c_char), intent(in) :: output(width,10,points),reference(width,10,points)
    integer(c_int), intent(in) :: accuracy(points),reference_accuracy(points)
    integer, intent(inout) :: failures
    integer :: point,field
    do point=1,points
      do field=1,10
        if(all(output(:,field,point)==reference(:,field,point))) cycle
!$omp critical
        print *, 'MPFR batch worker/fixture/point/field:',worker,fixture,point,field
        print *, 'Expected: ',pack(reference(:,field,point),reference(:,field,point)/=c_null_char)
        print *, 'Actual:   ',pack(output(:,field,point),output(:,field,point)/=c_null_char)
!$omp end critical
        failures=failures+1
      end do
    end do
    if(any(accuracy/=reference_accuracy)) then
      print *, 'MPFR batch worker/fixture/accuracy/expected:',worker,fixture,accuracy,reference_accuracy
      failures=failures+1
    end if
  end subroutine

  recursive subroutine evaluate(worker,fixture,failure,output,accuracy,failures)
    integer, intent(in) :: worker,fixture,failure
    character(kind=c_char), intent(out) :: output(width,10,points)
    integer(c_int), intent(out) :: accuracy(points)
    integer, intent(inout) :: failures
    character(kind=c_char) :: parameter(width,2),coordinates(width,points)
    character(kind=c_char) :: short_parameter(8,2),short_coordinates(8,points),short_output(8,10,points)
    integer :: family,geometry,complex_parameter,mode,degree,bits,rounding,normalize,expected_status,status
    integer :: saved_bits,saved_rounding,saved_status
    integer(c_int64_t) :: saved_budget,budget

    family=(fixture-1)/4
    geometry=mod(family,2)
    complex_parameter=family/2
    mode=mod(fixture-1,4)
    degree=mod(family,2)
    bits=requested_bits(1+mod(worker+fixture-2,workers))
    rounding=mod(worker+fixture,5)
    normalize=mod(fixture,2)
    expected_status=0
    budget=normal_budget
    parameter=c_null_char
    parameter(1,1)='1'
    parameter(1,2)='0'
    if(complex_parameter==1) parameter(1:3,2)=[character(kind=c_char) :: '1','/','4']
    coordinates=c_null_char
    if(mode<=1) then
      coordinates(1:3,1)=[character(kind=c_char) :: '1','/','4']
      coordinates(1:4,2)=[character(kind=c_char) :: '-','1','/','2']
    else
      coordinates(1,1)='2'
      coordinates(1,2)='3'
    end if
    select case(failure)
    case(1)
      parameter(1,1)='q'
      expected_status=1
    case(2)
      parameter(1:3,1)=[character(kind=c_char) :: 'N','a','N']
      expected_status=1
    case(3)
      parameter(1,1)='0'
      expected_status=2
    case(4)
      geometry=1
      mode=3
      coordinates=c_null_char
      coordinates(1,:)='0'
      expected_status=2
    case(5)
      budget=1048576_c_int64_t
      expected_status=4
    case(6)
      expected_status=3
    case(7)
      degree=huge(0_c_int)
      expected_status=4
    case(8)
      normalize=2
      expected_status=1
    case(9)
      parameter(width,1)='1'
      expected_status=1
    case(10)
      mode=3
      coordinates=c_null_char
      coordinates(1,1)='2'
      coordinates(1,2)='q'
      expected_status=1
    case(11)
      mode=1
      coordinates=c_null_char
      coordinates(1,1)='0'
      coordinates(1,2)='2'
      expected_status=2
    end select
    saved_bits=mpfr_bits
    saved_rounding=mpfr_rounding
    saved_status=mpfr_status
    saved_budget=mpfr_workspace_limit_bytes
    mpfr_bits=400+8*worker
    mpfr_rounding=mod(worker+2,5)
    mpfr_status=20+worker
    mpfr_workspace_limit_bytes=budget
    if(.not.caches_empty()) error stop 'MPFR batch retained a previous kernel cache'
    call seed_caches()
    if(failure==6) then
      short_parameter=parameter(:8,:)
      short_coordinates=coordinates(:8,:)
      call spheroidal_mpfr_batch(geometry,complex_parameter,0,degree,mode,normalize,bits,rounding, &
        points,8,short_parameter,short_coordinates,short_output,accuracy,status)
    else
      call spheroidal_mpfr_batch(geometry,complex_parameter,0,degree,mode,normalize,bits,rounding, &
        points,width,parameter,coordinates,output,accuracy,status)
    end if
    if(status/=expected_status) then
      print *, 'Unexpected MPFR status, worker/fixture/failure/status/expected:',worker,fixture,failure,status,expected_status
      failures=failures+1
    end if
    if(mpfr_bits/=400+8*worker.or.mpfr_rounding/=mod(worker+2,5).or.mpfr_status/=20+worker.or. &
      mpfr_workspace_limit_bytes/=budget) then
      print *, 'MPFR batch caller context changed for worker/fixture/failure:',worker,fixture,failure
      failures=failures+1
    end if
    if(failure>=7.and.failure<=9) then
      if(.not.seeded_caches_intact()) error stop 'Early rejection modified caller caches'
      call p_clear()
      call o_clear()
      call cp_clear()
      call co_clear()
    end if
    if(.not.caches_empty()) error stop 'MPFR batch did not release every kernel cache'
    if(status==0) then
      if(any(output(width,:,:)/=c_null_char).or.any(output(1,1,:)==c_null_char)) failures=failures+1
      if(mode>0.and.any(output(1,3,:)==c_null_char)) failures=failures+1
      if(mode==3.and.any(output(1,9,:)==c_null_char)) failures=failures+1
    end if
    mpfr_bits=saved_bits
    mpfr_rounding=saved_rounding
    mpfr_status=saved_status
    mpfr_workspace_limit_bytes=saved_budget
  end subroutine

  recursive subroutine seed_caches()
    ! Seed every kernel at the caller's precision to expose incomplete cross-kernel cleanup.
    allocate(p_pleg(1)%barg(1))
    allocate(o_pleg(1)%barg(1))
    allocate(cp_pleg(1)%barg(1))
    allocate(co_pleg(1)%barg(1))
    p_pleg(1)%barg=mpreal(1)/3
    o_pleg(1)%barg=mpreal(1)/3
    cp_pleg(1)%barg=mpreal(1)/3
    co_pleg(1)%barg=mpreal(1)/3
    p_pleg(1)%valid=.true.
    o_pleg(1)%valid=.true.
    cp_pleg(1)%valid=.true.
    co_pleg(1)%valid=.true.
    p_pleg(1)%bytes=128
    o_pleg(1)%bytes=128
    cp_pleg(1)%bytes=128
    co_pleg(1)%bytes=128
    p_qleg(1)%x1=mpreal(1)/3
    cp_qleg(1)%x1=mpreal(1)/3
    p_qleg(1)%valid=.true.
    cp_qleg(1)%valid=.true.
    p_qleg(1)%bytes=128
    cp_qleg(1)%bytes=128
    allocate(p_gauss%x(1))
    allocate(o_gauss%x(1))
    allocate(cp_gauss%x(1))
    allocate(co_gauss%x(1))
    p_gauss%x=mpreal(1)/3
    o_gauss%x=mpreal(1)/3
    cp_gauss%x=mpreal(1)/3
    co_gauss%x=mpreal(1)/3
    p_gauss%valid=.true.
    o_gauss%valid=.true.
    cp_gauss%valid=.true.
    co_gauss%valid=.true.
    p_gauss%bytes=128
    o_gauss%bytes=128
    cp_gauss%bytes=128
    co_gauss%bytes=128
  end subroutine

  recursive logical function seeded_caches_intact() result(intact)
    type(mp_real) :: expected
    intact=.false.
    if(.not.(p_pleg(1)%valid.and.o_pleg(1)%valid.and.cp_pleg(1)%valid.and.co_pleg(1)%valid)) return
    if(.not.(p_qleg(1)%valid.and.cp_qleg(1)%valid)) return
    if(.not.(p_gauss%valid.and.o_gauss%valid.and.cp_gauss%valid.and.co_gauss%valid)) return
    if(.not.(allocated(p_pleg(1)%barg).and.allocated(o_pleg(1)%barg).and. &
      allocated(cp_pleg(1)%barg).and.allocated(co_pleg(1)%barg))) return
    if(.not.(allocated(p_gauss%x).and.allocated(o_gauss%x).and.allocated(cp_gauss%x).and.allocated(co_gauss%x))) return
    expected=mpreal(1)/3
    if(any(p_pleg(1)%barg/=expected).or.any(o_pleg(1)%barg/=expected).or. &
      any(cp_pleg(1)%barg/=expected).or.any(co_pleg(1)%barg/=expected)) return
    if(p_qleg(1)%x1/=expected.or.cp_qleg(1)%x1/=expected) return
    if(any(p_gauss%x/=expected).or.any(o_gauss%x/=expected).or.any(cp_gauss%x/=expected).or.any(co_gauss%x/=expected)) return
    if(p_pleg(1)%bytes/=128.or.o_pleg(1)%bytes/=128.or.cp_pleg(1)%bytes/=128.or.co_pleg(1)%bytes/=128) return
    if(p_qleg(1)%bytes/=128.or.cp_qleg(1)%bytes/=128) return
    if(p_gauss%bytes/=128.or.o_gauss%bytes/=128.or.cp_gauss%bytes/=128.or.co_gauss%bytes/=128) return
    intact=.true.
  end function

  recursive logical function caches_empty() result(empty)
    integer :: slot
    empty=.false.
    if(any(p_pleg%valid).or.any(o_pleg%valid).or.any(cp_pleg%valid).or.any(co_pleg%valid)) return
    if(any(p_qleg%valid).or.any(cp_qleg%valid)) return
    if(sum(p_pleg%bytes)+sum(o_pleg%bytes)+sum(cp_pleg%bytes)+sum(co_pleg%bytes)/=0) return
    if(sum(p_qleg%bytes)+sum(cp_qleg%bytes)/=0) return
    do slot=1,size(p_pleg)
      if(allocated(p_pleg(slot)%barg).or.allocated(o_pleg(slot)%barg).or. &
        allocated(cp_pleg(slot)%barg).or.allocated(co_pleg(slot)%barg)) return
    end do
    do slot=1,size(p_qleg)
      if(allocated(p_qleg(slot)%x1%data).or.allocated(cp_qleg(slot)%x1%data)) return
    end do
    if(p_gauss%valid.or.o_gauss%valid.or.cp_gauss%valid.or.co_gauss%valid) return
    if(p_gauss%bytes+o_gauss%bytes+cp_gauss%bytes+co_gauss%bytes/=0) return
    if(allocated(p_gauss%x).or.allocated(o_gauss%x).or.allocated(cp_gauss%x).or.allocated(co_gauss%x)) return
    empty=.true.
  end function
end program
