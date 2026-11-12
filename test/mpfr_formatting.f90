program mpfr_formatting
  use spheroidal_mpfr_real
  use omp_lib
  implicit none
  integer, parameter :: workers=4,precisions(workers)=[159,160,192,256]
  character(len=128) :: reference(workers)
  integer :: worker,failures,actual_workers
  logical :: saved_dynamic
  do worker=1,workers
    mpfr_bits=precisions(worker)
    reference(worker)=mp_to_string(mpreal('1/7'))
  end do
  failures=0
  actual_workers=0
  saved_dynamic=omp_get_dynamic()
  call omp_set_dynamic(.false.)
!$omp parallel num_threads(workers) default(none) shared(reference,actual_workers) private(worker) reduction(+:failures)
!$omp single
  actual_workers=omp_get_num_threads()
!$omp end single
  worker=omp_get_thread_num()+1
  call check_worker(worker,reference(worker),failures)
!$omp end parallel
  call omp_set_dynamic(saved_dynamic)
  if (actual_workers/=workers) error stop 'Formatting test requires four workers'
  if (failures/=0) error stop 'Concurrent MPFR formatting differs from serial output'
  print *, 'Concurrent MPFR text lengths and round trips passed'
contains
  recursive subroutine check_worker(worker,reference,failures)
    integer, intent(in) :: worker
    character(len=*), intent(in) :: reference
    integer, intent(inout) :: failures
    type(mp_real) :: value
    character(len=:), allocatable :: text
    integer :: iteration
    mpfr_bits=precisions(worker)
    value=mpreal('1/7')
    do iteration=1,10000
      call mp_get_string(value,text)
      if (text/=trim(reference)) failures=failures+1
      if (mpreal(text)/=value) failures=failures+1
    end do
  end subroutine
end program
