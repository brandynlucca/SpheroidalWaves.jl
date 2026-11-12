program mpfr_oblate_threads
  use, intrinsic :: iso_c_binding, only: c_char,c_null_char,c_int,c_int64_t
  use spheroidal_mpfr_real, only: mpfr_bits,mpfr_rounding,mpfr_status,mpfr_workspace_limit_bytes
  use spheroidal_mpfr_batch_api, only: spheroidal_mpfr_batch
  use omp_lib
  implicit none
  integer, parameter :: workers=4,width=192,requested_bits(workers)=[159,160,192,256]
  integer(c_int64_t), parameter :: normal_budget=268435456_c_int64_t
  character(kind=c_char) :: reference_output(width,10,1,5,2,workers)
  integer(c_int) :: reference_accuracy(1,5,2,workers)
  integer :: worker,rounding,mode,failures,actual_workers
  logical :: saved_dynamic

  failures=0
  do worker=1,workers
    do rounding=0,4
      do mode=0,1
        call evaluate(worker,rounding,mode,normal_budget,0,.false., &
          reference_output(:,:,:,rounding+1,mode+1,worker), &
          reference_accuracy(:,rounding+1,mode+1,worker),failures)
      end do
    end do
  end do
  if (failures/=0) error stop 'Serial MPFR oblate references failed'

  saved_dynamic=omp_get_dynamic()
  call omp_set_dynamic(.false.)
  actual_workers=0
!$omp parallel num_threads(workers) default(none) &
!$omp shared(reference_output,reference_accuracy,actual_workers) private(worker) reduction(+:failures)
!$omp single
  actual_workers=omp_get_num_threads()
!$omp end single
  worker=omp_get_thread_num()+1
  call check_worker(worker,reference_output(:,:,:,:,:,worker),reference_accuracy(:,:,:,worker),failures)
!$omp end parallel
  call omp_set_dynamic(saved_dynamic)
  if (actual_workers/=workers) error stop 'MPFR oblate test requires four OpenMP workers'
  if (failures/=0) then
    print *, 'Concurrent MPFR oblate failures:',failures
    error stop 1
  end if
  print *, 'MPFR oblate kernels matched four-thread references; caller context and workspace recovery passed'

contains
  recursive subroutine check_worker(worker,reference,reference_accuracy,failures)
    integer, intent(in) :: worker
    character(kind=c_char), intent(in) :: reference(width,10,1,5,2)
    integer(c_int), intent(in) :: reference_accuracy(1,5,2)
    integer, intent(inout) :: failures
    character(kind=c_char) :: output(width,10,1)
    integer(c_int) :: accuracy(1)
    integer(c_int64_t) :: budget
    integer :: iteration,rounding,mode,expected_status

    do iteration=1,2
      do rounding=0,4
        do mode=0,1
          call evaluate(worker,rounding,mode,normal_budget,0,.false.,output,accuracy,failures)
          if (any(output/=reference(:,:,:,rounding+1,mode+1))) failures=failures+1
          if (any(accuracy/=reference_accuracy(:,rounding+1,mode+1))) failures=failures+1
        end do
      end do
    end do

    budget=normal_budget
    expected_status=0
    if (worker==1) then
      budget=1048576_c_int64_t
      expected_status=4
    end if
    call evaluate(worker,worker,1,budget,expected_status,.true.,output,accuracy,failures)
    if (expected_status==0) then
      if (any(output/=reference(:,:,:,worker+1,2))) failures=failures+1
      if (any(accuracy/=reference_accuracy(:,worker+1,2))) failures=failures+1
    end if
    call evaluate(worker,worker,1,normal_budget,0,.true.,output,accuracy,failures)
    if (any(output/=reference(:,:,:,worker+1,2))) failures=failures+1
    if (any(accuracy/=reference_accuracy(:,worker+1,2))) failures=failures+1
  end subroutine

  recursive subroutine evaluate(worker,rounding,mode,budget,expected_status,synchronize,output,accuracy,failures)
    integer, intent(in) :: worker,rounding,mode,expected_status
    integer(c_int64_t), intent(in) :: budget
    logical, intent(in) :: synchronize
    character(kind=c_char), intent(out) :: output(width,10,1)
    integer(c_int), intent(out) :: accuracy(1)
    integer, intent(inout) :: failures
    character(kind=c_char) :: parameter(width,2),coordinate(width,1)
    integer(c_int) :: status,saved_bits,saved_rounding,saved_status
    integer(c_int64_t) :: saved_budget

    parameter=c_null_char
    parameter(1,1)='1'
    parameter(1,2)='0'
    coordinate=c_null_char
    coordinate(1:3,1)=[character(kind=c_char) :: '1','/','4']
    saved_bits=mpfr_bits
    saved_rounding=mpfr_rounding
    saved_status=mpfr_status
    saved_budget=mpfr_workspace_limit_bytes
    mpfr_bits=384+8*worker
    mpfr_rounding=mod(worker+1,5)
    mpfr_status=20+worker
    mpfr_workspace_limit_bytes=budget
    if (synchronize) then
!$omp barrier
    end if
    call spheroidal_mpfr_batch(0,0,0,1,mode,0,requested_bits(worker),rounding,1,width, &
      parameter,coordinate,output,accuracy,status)
    if (synchronize) then
!$omp barrier
    end if
    if (status/=expected_status) failures=failures+1
    if (mpfr_bits/=384+8*worker.or.mpfr_rounding/=mod(worker+1,5).or. &
        mpfr_status/=20+worker.or.mpfr_workspace_limit_bytes/=budget) failures=failures+1
    if (status==0) then
      if (output(1,1,1)==c_null_char) failures=failures+1
      if (mode==1) then
        if (output(1,3,1)==c_null_char.or.output(1,5,1)==c_null_char) failures=failures+1
      end if
    end if
    mpfr_bits=saved_bits
    mpfr_rounding=saved_rounding
    mpfr_status=saved_status
    mpfr_workspace_limit_bytes=saved_budget
  end subroutine
end program
