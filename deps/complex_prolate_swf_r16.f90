! ---------------------------------------------------------------------------
! Modified Source Notice
! ---------------------------------------------------------------------------
! This file is a modified derivative of the original `complex_prolate_swf.f90`
! implementation developed by Arnie Lee Van Buren
! (Mathieu and Spheroidal Wave Functions project).
!
! Original upstream source:
!   https://github.com/MathieuandSpheroidalWaveFunctions/complex_prolate_swf
!
! Local modifications in this copy:
!   1) Added local `param` module defaults with debug/warn/output disabled.
!   2) Retained in-memory callable API (`cprofcn`) for integration use.
!   3) Added extensive caching subsystem for Legendre polynomials (pleg_cache),
!      Associated Legendre quotients (qleg_cache), and Gauss-Legendre quadrature
!      (gauss_cache) to avoid recomputation on repeated calls.
!   4) Integrated cached wrappers into all call sites within main cprofcn kernel.
!   5) Preserved original numerical kernels and attribution comments.
!
! Modified by: Brandyn M. Lucca; March 2026
! Note: This file is NOT a pristine upstream copy.
! ---------------------------------------------------------------------------

module param
    integer, parameter :: knd = selected_real_kind(16)
    integer, parameter :: knd1 = selected_real_kind(16)
    logical, parameter :: debug = .false.
    logical, parameter :: warn = .false.
    logical, parameter :: output = .false.
    logical, parameter :: suffix = .false.
end module param

module complex_prolate_swf
      use param

 integer, parameter :: pleg_cache_slots = 48
 integer, parameter :: qleg_cache_slots = 48

 type :: pleg_cache_entry
   logical :: valid = .false.
   integer :: m = -1
   integer :: iopd = -1
   integer :: narg = 0
   integer :: maxt = 0
   integer :: maxp = 0
   integer :: lim = 0
   integer :: ndec = 0
   integer :: nex = 0
   real(knd), allocatable :: barg(:)
   real(knd), allocatable :: pr(:,:), pdr(:,:)
   real(knd), allocatable :: pdnorm(:), pnorm(:)
   integer, allocatable :: ipdnorm(:), ipnorm(:)
   real(knd), allocatable :: alpha(:), beta(:), gamma(:)
   real(knd), allocatable :: coefa(:), coefb(:), coefc(:), coefd(:), coefe(:)
 end type

 type :: qleg_cache_entry
   logical :: valid = .false.
   integer :: m = -1
   integer :: lnum = 0
   integer :: limq = 0
   integer :: maxq = 0
   integer :: ndec = 0
   integer :: iqdml = 0
   integer :: iqml = 0
   integer :: itermpq = 0
   real(knd) :: x1 = 0.0_knd
   real(knd) :: qdml = 0.0_knd
   real(knd) :: qml = 0.0_knd
   real(knd) :: termpq = 0.0_knd
   real(knd), allocatable :: qdr(:), qr(:), qdl(:), ql(:)
   integer, allocatable :: iqdl(:), iql(:)
 end type

 type :: gauss_cache_entry
   logical :: valid = .false.
   integer :: ndec = 0
   integer :: n = 0
   real(knd), allocatable :: x(:), w(:)
 end type

 type(pleg_cache_entry), save :: pleg_cache(pleg_cache_slots)
 type(qleg_cache_entry), save :: qleg_cache(qleg_cache_slots)
 type(gauss_cache_entry), save :: gauss_cache
 integer, save :: pleg_cache_next = 1
 integer, save :: qleg_cache_next = 1

 contains

  logical function cache_real_equal(a, b)
    real(knd), intent(in) :: a, b
    real(knd) :: scale, tol
    scale = max(1.0_knd, max(abs(a), abs(b)))
    tol = 64.0_knd * epsilon(1.0_knd) * scale
    cache_real_equal = abs(a - b) <= tol
  end function

  logical function cache_real_vector_equal(a, b, n)
    real(knd), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    integer :: i
    cache_real_vector_equal = .false.
    if(size(a) < n .or. size(b) < n) return
    do i = 1, n
      if(.not. cache_real_equal(a(i), b(i))) return
    end do
    cache_real_vector_equal = .true.
  end function

  subroutine clear_pleg_cache_slot(idx)
    integer, intent(in) :: idx
    pleg_cache(idx)%valid = .false.
    pleg_cache(idx)%m = -1
    pleg_cache(idx)%iopd = -1
    pleg_cache(idx)%narg = 0
    pleg_cache(idx)%maxt = 0
    pleg_cache(idx)%maxp = 0
    pleg_cache(idx)%lim = 0
    pleg_cache(idx)%ndec = 0
    pleg_cache(idx)%nex = 0
    if(allocated(pleg_cache(idx)%barg)) deallocate(pleg_cache(idx)%barg)
    if(allocated(pleg_cache(idx)%pr)) deallocate(pleg_cache(idx)%pr)
    if(allocated(pleg_cache(idx)%pdr)) deallocate(pleg_cache(idx)%pdr)
    if(allocated(pleg_cache(idx)%pdnorm)) deallocate(pleg_cache(idx)%pdnorm)
    if(allocated(pleg_cache(idx)%pnorm)) deallocate(pleg_cache(idx)%pnorm)
    if(allocated(pleg_cache(idx)%ipdnorm)) deallocate(pleg_cache(idx)%ipdnorm)
    if(allocated(pleg_cache(idx)%ipnorm)) deallocate(pleg_cache(idx)%ipnorm)
    if(allocated(pleg_cache(idx)%alpha)) deallocate(pleg_cache(idx)%alpha)
    if(allocated(pleg_cache(idx)%beta)) deallocate(pleg_cache(idx)%beta)
    if(allocated(pleg_cache(idx)%gamma)) deallocate(pleg_cache(idx)%gamma)
    if(allocated(pleg_cache(idx)%coefa)) deallocate(pleg_cache(idx)%coefa)
    if(allocated(pleg_cache(idx)%coefb)) deallocate(pleg_cache(idx)%coefb)
    if(allocated(pleg_cache(idx)%coefc)) deallocate(pleg_cache(idx)%coefc)
    if(allocated(pleg_cache(idx)%coefd)) deallocate(pleg_cache(idx)%coefd)
    if(allocated(pleg_cache(idx)%coefe)) deallocate(pleg_cache(idx)%coefe)
  end subroutine

  integer function find_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp)
    integer, intent(in) :: m, iopd, ndec, nex, narg, maxt, lim, maxp
    real(knd), intent(in) :: barg(:)
    integer :: i
    find_pleg_cache = 0
    do i = 1, pleg_cache_slots
      if(.not. pleg_cache(i)%valid) cycle
      if(pleg_cache(i)%m /= m) cycle
      if(pleg_cache(i)%iopd /= iopd) cycle
      if(pleg_cache(i)%ndec /= ndec) cycle
      if(pleg_cache(i)%nex /= nex) cycle
      if(pleg_cache(i)%narg /= narg) cycle
      if(pleg_cache(i)%maxt /= maxt) cycle
      if(pleg_cache(i)%lim < lim) cycle
      if(pleg_cache(i)%maxp < maxp) cycle
      if(.not. cache_real_vector_equal(pleg_cache(i)%barg, barg, narg)) cycle
      find_pleg_cache = i
      return
    end do
  end function

  subroutine store_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    integer, intent(in) :: m, iopd, ndec, nex, narg, maxt, lim, maxp
    real(knd), intent(in) :: barg(maxt), pr(maxt, maxp), pdr(maxt, maxp), pdnorm(maxt), pnorm(maxt), alpha(maxp), beta(maxp), gamma(maxp), coefa(maxp), coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp)
    integer, intent(in) :: ipdnorm(maxt), ipnorm(maxt)
    integer :: idx
    idx = find_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp)
    if(idx == 0) then
      idx = 1
      do while(idx <= pleg_cache_slots)
        if(.not. pleg_cache(idx)%valid) exit
        idx = idx + 1
      end do
      if(idx > pleg_cache_slots) then
        idx = pleg_cache_next
        pleg_cache_next = pleg_cache_next + 1
        if(pleg_cache_next > pleg_cache_slots) pleg_cache_next = 1
      end if
    end if
    call clear_pleg_cache_slot(idx)
    allocate(pleg_cache(idx)%barg(narg))
    allocate(pleg_cache(idx)%pr(maxt, maxp))
    allocate(pleg_cache(idx)%pdr(maxt, maxp))
    allocate(pleg_cache(idx)%pdnorm(maxt))
    allocate(pleg_cache(idx)%pnorm(maxt))
    allocate(pleg_cache(idx)%ipdnorm(maxt))
    allocate(pleg_cache(idx)%ipnorm(maxt))
    allocate(pleg_cache(idx)%alpha(maxp))
    allocate(pleg_cache(idx)%beta(maxp))
    allocate(pleg_cache(idx)%gamma(maxp))
    allocate(pleg_cache(idx)%coefa(maxp))
    allocate(pleg_cache(idx)%coefb(maxp))
    allocate(pleg_cache(idx)%coefc(maxp))
    allocate(pleg_cache(idx)%coefd(maxp))
    allocate(pleg_cache(idx)%coefe(maxp))
    pleg_cache(idx)%valid = .true.
    pleg_cache(idx)%m = m
    pleg_cache(idx)%iopd = iopd
    pleg_cache(idx)%narg = narg
    pleg_cache(idx)%maxt = maxt
    pleg_cache(idx)%maxp = maxp
    pleg_cache(idx)%lim = lim
    pleg_cache(idx)%ndec = ndec
    pleg_cache(idx)%nex = nex
    pleg_cache(idx)%barg(1:narg) = barg(1:narg)
    pleg_cache(idx)%pr = 0.0_knd
    pleg_cache(idx)%pdr = 0.0_knd
    pleg_cache(idx)%alpha = 0.0_knd
    pleg_cache(idx)%beta = 0.0_knd
    pleg_cache(idx)%gamma = 0.0_knd
    pleg_cache(idx)%coefa = 0.0_knd
    pleg_cache(idx)%coefb = 0.0_knd
    pleg_cache(idx)%coefc = 0.0_knd
    pleg_cache(idx)%coefd = 0.0_knd
    pleg_cache(idx)%coefe = 0.0_knd
    pleg_cache(idx)%pr(1:maxt, 1:lim) = pr(1:maxt, 1:lim)
    pleg_cache(idx)%pdr(1:maxt, 1:lim) = pdr(1:maxt, 1:lim)
    pleg_cache(idx)%pdnorm(1:maxt) = pdnorm(1:maxt)
    pleg_cache(idx)%pnorm(1:maxt) = pnorm(1:maxt)
    pleg_cache(idx)%ipdnorm(1:maxt) = ipdnorm(1:maxt)
    pleg_cache(idx)%ipnorm(1:maxt) = ipnorm(1:maxt)
    pleg_cache(idx)%alpha(1:lim) = alpha(1:lim)
    pleg_cache(idx)%beta(1:lim) = beta(1:lim)
    pleg_cache(idx)%gamma(1:lim) = gamma(1:lim)
    pleg_cache(idx)%coefa(1:lim) = coefa(1:lim)
    pleg_cache(idx)%coefb(1:lim) = coefb(1:lim)
    pleg_cache(idx)%coefc(1:lim) = coefc(1:lim)
    pleg_cache(idx)%coefd(1:lim) = coefd(1:lim)
    pleg_cache(idx)%coefe(1:lim) = coefe(1:lim)
  end subroutine

  subroutine pleg_cached(m, lim, maxp, limcsav, iopd, ndec, nex, barg, narg, maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    integer, intent(in) :: m, lim, maxp, iopd, ndec, nex, narg, maxt
    integer, intent(inout) :: limcsav
    real(knd), intent(inout) :: barg(maxt)
    real(knd), intent(out) :: pr(maxt, maxp), pdr(maxt, maxp), pdnorm(maxt), pnorm(maxt), alpha(maxp), beta(maxp), gamma(maxp), coefa(maxp), coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp)
    integer, intent(out) :: ipdnorm(maxt), ipnorm(maxt)
    integer :: idx
    idx = find_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp)
    if(idx /= 0) then
      call load_pleg_cache(idx, maxt, narg, lim, maxp, barg, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
      limcsav = max(limcsav, lim)
      return
    end if
                call pleg(m, lim, maxp, ndec, nex, limcsav, iopd, barg, narg, maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    call store_pleg_cache(m, iopd, ndec, nex, barg, narg, maxt, lim, maxp, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
  end subroutine

  subroutine load_pleg_cache(idx, maxt, narg, lim, maxp, barg, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
    integer, intent(in) :: idx, maxt, narg, lim, maxp
    real(knd), intent(out) :: barg(maxt), pr(maxt, maxp), pdr(maxt, maxp), pdnorm(maxt), pnorm(maxt), alpha(maxp), beta(maxp), gamma(maxp), coefa(maxp), coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp)
    integer, intent(out) :: ipdnorm(maxt), ipnorm(maxt)
    barg(1:narg) = pleg_cache(idx)%barg(1:narg)
    pr(1:maxt, 1:lim) = pleg_cache(idx)%pr(1:maxt, 1:lim)
    pdr(1:maxt, 1:lim) = pleg_cache(idx)%pdr(1:maxt, 1:lim)
    pdnorm(1:maxt) = pleg_cache(idx)%pdnorm(1:maxt)
    ipdnorm(1:maxt) = pleg_cache(idx)%ipdnorm(1:maxt)
    pnorm(1:maxt) = pleg_cache(idx)%pnorm(1:maxt)
    ipnorm(1:maxt) = pleg_cache(idx)%ipnorm(1:maxt)
    alpha(1:lim) = pleg_cache(idx)%alpha(1:lim)
    beta(1:lim) = pleg_cache(idx)%beta(1:lim)
    gamma(1:lim) = pleg_cache(idx)%gamma(1:lim)
    coefa(1:lim) = pleg_cache(idx)%coefa(1:lim)
    coefb(1:lim) = pleg_cache(idx)%coefb(1:lim)
    coefc(1:lim) = pleg_cache(idx)%coefc(1:lim)
    coefd(1:lim) = pleg_cache(idx)%coefd(1:lim)
    coefe(1:lim) = pleg_cache(idx)%coefe(1:lim)
  end subroutine

  integer function find_qleg_cache(m, lnum, limq, maxq, x1, ndec)
    integer, intent(in) :: m, lnum, limq, maxq, ndec
    real(knd), intent(in) :: x1
    integer :: i
    find_qleg_cache = 0
    do i = 1, qleg_cache_slots
      if(.not. qleg_cache(i)%valid) cycle
      if(qleg_cache(i)%m /= m) cycle
      if(qleg_cache(i)%ndec /= ndec) cycle
      if(.not. cache_real_equal(qleg_cache(i)%x1, x1)) cycle
      if(qleg_cache(i)%lnum < lnum) cycle
      if(qleg_cache(i)%limq < limq) cycle
      if(qleg_cache(i)%maxq < maxq) cycle
      find_qleg_cache = i
      return
    end do
  end function

  subroutine clear_qleg_cache_slot(idx)
    integer, intent(in) :: idx
    qleg_cache(idx)%valid = .false.
    qleg_cache(idx)%m = -1
    qleg_cache(idx)%lnum = 0
    qleg_cache(idx)%limq = 0
    qleg_cache(idx)%maxq = 0
    qleg_cache(idx)%ndec = 0
    qleg_cache(idx)%iqdml = 0
    qleg_cache(idx)%iqml = 0
    qleg_cache(idx)%itermpq = 0
    qleg_cache(idx)%x1 = 0.0_knd
    qleg_cache(idx)%qdml = 0.0_knd
    qleg_cache(idx)%qml = 0.0_knd
    qleg_cache(idx)%termpq = 0.0_knd
    if(allocated(qleg_cache(idx)%qdr)) deallocate(qleg_cache(idx)%qdr)
    if(allocated(qleg_cache(idx)%qr)) deallocate(qleg_cache(idx)%qr)
    if(allocated(qleg_cache(idx)%qdl)) deallocate(qleg_cache(idx)%qdl)
    if(allocated(qleg_cache(idx)%ql)) deallocate(qleg_cache(idx)%ql)
    if(allocated(qleg_cache(idx)%iqdl)) deallocate(qleg_cache(idx)%iqdl)
    if(allocated(qleg_cache(idx)%iql)) deallocate(qleg_cache(idx)%iql)
  end subroutine

  subroutine store_qleg_cache(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
    integer, intent(in) :: m, lnum, limq, maxq, ndec, iqdml, iqml, itermpq
    real(knd), intent(in) :: x1, qdr(maxq), qdml, qdl(lnum), qr(maxq), qml, ql(lnum), termpq
    integer, intent(in) :: iqdl(lnum), iql(lnum)
    integer :: idx
    idx = find_qleg_cache(m, lnum, limq, maxq, x1, ndec)
    if(idx == 0) then
      idx = 1
      do while(idx <= qleg_cache_slots)
        if(.not. qleg_cache(idx)%valid) exit
        idx = idx + 1
      end do
      if(idx > qleg_cache_slots) then
        idx = qleg_cache_next
        qleg_cache_next = qleg_cache_next + 1
        if(qleg_cache_next > qleg_cache_slots) qleg_cache_next = 1
      end if
    end if
    call clear_qleg_cache_slot(idx)
    allocate(qleg_cache(idx)%qdr(maxq))
    allocate(qleg_cache(idx)%qr(maxq))
    allocate(qleg_cache(idx)%qdl(lnum))
    allocate(qleg_cache(idx)%ql(lnum))
    allocate(qleg_cache(idx)%iqdl(lnum))
    allocate(qleg_cache(idx)%iql(lnum))
    qleg_cache(idx)%valid = .true.
    qleg_cache(idx)%m = m
    qleg_cache(idx)%lnum = lnum
    qleg_cache(idx)%limq = limq
    qleg_cache(idx)%maxq = maxq
    qleg_cache(idx)%ndec = ndec
    qleg_cache(idx)%x1 = x1
    qleg_cache(idx)%qdml = qdml
    qleg_cache(idx)%iqdml = iqdml
    qleg_cache(idx)%qml = qml
    qleg_cache(idx)%iqml = iqml
    qleg_cache(idx)%termpq = termpq
    qleg_cache(idx)%itermpq = itermpq
    qleg_cache(idx)%qdr(1:maxq) = qdr(1:maxq)
    qleg_cache(idx)%qr(1:maxq) = qr(1:maxq)
    qleg_cache(idx)%qdl(1:lnum) = qdl(1:lnum)
    qleg_cache(idx)%ql(1:lnum) = ql(1:lnum)
    qleg_cache(idx)%iqdl(1:lnum) = iqdl(1:lnum)
    qleg_cache(idx)%iql(1:lnum) = iql(1:lnum)
  end subroutine

  subroutine qleg_cached(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
    integer, intent(in) :: m, lnum, limq, maxq, ndec
    real(knd), intent(in) :: x1
    real(knd), intent(out) :: qdr(maxq), qdml, qdl(lnum), qr(maxq), qml, ql(lnum), termpq
    integer, intent(out) :: iqdml, iqdl(lnum), iqml, iql(lnum), itermpq
    integer :: idx
    idx = find_qleg_cache(m, lnum, limq, maxq, x1, ndec)
    if(idx /= 0) then
      qdr(1:maxq) = qleg_cache(idx)%qdr(1:maxq)
      qr(1:maxq) = qleg_cache(idx)%qr(1:maxq)
      qdl(1:lnum) = qleg_cache(idx)%qdl(1:lnum)
      ql(1:lnum) = qleg_cache(idx)%ql(1:lnum)
      iqdl(1:lnum) = qleg_cache(idx)%iqdl(1:lnum)
      iql(1:lnum) = qleg_cache(idx)%iql(1:lnum)
      qdml = qleg_cache(idx)%qdml
      iqdml = qleg_cache(idx)%iqdml
      qml = qleg_cache(idx)%qml
      iqml = qleg_cache(idx)%iqml
      termpq = qleg_cache(idx)%termpq
      itermpq = qleg_cache(idx)%itermpq
      return
    end if
                call qleg(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
    call store_qleg_cache(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
  end subroutine

  subroutine gauss_cached(ndec, n, x, w)
    integer, intent(in) :: ndec, n
    real(knd), intent(out) :: x(n), w(n)
    if(gauss_cache%valid) then
      if(gauss_cache%ndec == ndec .and. gauss_cache%n == n) then
        x(1:n) = gauss_cache%x(1:n)
        w(1:n) = gauss_cache%w(1:n)
        return
      end if
      if(allocated(gauss_cache%x)) deallocate(gauss_cache%x)
      if(allocated(gauss_cache%w)) deallocate(gauss_cache%w)
    end if
        call gauss(n, ndec, x, w)
    allocate(gauss_cache%x(n))
    allocate(gauss_cache%w(n))
    gauss_cache%valid = .true.
    gauss_cache%ndec = ndec
    gauss_cache%n = n
    gauss_cache%x(1:n) = x(1:n)
    gauss_cache%w(1:n) = w(1:n)
  end subroutine

!
     
      subroutine cprofcn(cc, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                       r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                       s1c, is1e, s1dc, is1de, naccs, naccds)

!      version 1.06 August 2024
!
!  subroutine version of the fortran program cprofcn developed about 2005 by
!  arnie lee van buren with technical support from jeffrey boisvert. For more
!  information see the github repository:
!  GitHub.com/MathieuandSpheroidalWaveFunctions/complex_prolate_swf.
!  Especially see the readme file and example input and output files.
!
!
!  purpose:     To calculate the first and second kind prolate radial
!               functions r1 and r2 and their first derivatives with
!               respect to the radial coordinate x for a given value of
!               the complex size parameter c, for a given value of the
!               order m and for lnum degrees l = m, m+1, ...,m+lnum-1.
!               To calculate the first kind prolate angular functions
!               s1 and their first derivatives with respect to the
!               angle coordinate eta for a given value of the complex size
!               parameter c, for a given value of the order m and for lnum
!               degrees l = m, m+1, ..., m+lnum-1
!
!  Cprofcn can be run in double precision, quadruple precision or a hybrid
!  where the Bouwkamp procedure to refine the eigenvalues is run in quadruple
!  precision while the remainder of the calculations are performed in double
!  precision. In the latter case, cprofcn switches to quadruple precision for
!  the Bouwkamp procedure whenever double precision fails to provide highly
!  accurate eigenvalues. See the discussion below about when to choose which
!  arithmetic.The choice is set in the module param provided in the github
!  repository. If this is not available, then create param as follows:
!    module param
!    integer, parameter :: knd = selected_real_kind(16)
!    integer, parameter :: knd1 = selected_real_kind(16)
!    logical, parameter :: debug = .true.
!    logical, parameter :: warn = .true.
!    logical, parameter :: output = .true.
!    logical, parameter :: suffix = .true.
!    end module param
!  Set the value of knd in the parenthesis to either 8 for double precision
!  or 16 for quadruple precision arithmetic. Set the value of knd1 to that used
!  for the Bouwkamp procedure. Note that if knd = 16, knd1 should also be 16. The
!  logicals in param are described below in the discussion of the output files.
!
!  Cprofcn provides good results over reasonable ranges of input parameters. A
!  description of the expected accuracy of cprofcn is given in the readme file.
!  Cprofcn provides function values for c complex = real(c) + i aimag(c) = cr + i ci,
!  where the imaginary part ci often accounts for losses in wave propagation.
!  Ci is assumed positive in cprofcn. If the user has a negative value for ci,
!  just run cprofcn with ci positive instead and take the complex conjugate of
!  the results, including the function values, eigenvalues, expansion coefficients,
!  and normalization factors.
!
!     Input and Output
!
!    Input and output parameters from the subroutine call statement are
!    defined below:
!
!          cc     : desired complex value of the size parameter (= kd/2,
!                   where k is the complex wavenumber and d is the
!                   interfocal length) [complex(knd)]
!
!          m      : desired value for the order m (integer)
!
!          lnum   : number of values desired for the degree l (integer)
!                   if lnum is less than 2*(real(c)+aimag(c))/pi it
!                   should chosen to be an even integer.
!
!          ioprad : (integer)
!                 : =0 if radial functions are not computed
!                 : =1 if radial functions of only the first
!                      kind and their first derivatives are
!                      computed
!                 : =2 if radial functions of both kinds and
!                      their first derivatives are computed
!
!          x1     : x - 1, where x is the radial coordinate x (a nominal
!                   value of 1.0e0_knd can be entered for x1 if ioprad
!                   = 0) If x1 = 0.0e0_knd, i.e., x = 1.0e0_knd,
!                   only radial functions of the first kind and their
!                   first derivatives are calculated. They are equal to
!                   zero unless m = 0. Radial functions of the second kind
!                   and their first derivatives are infinite for all m.
!                   Thus ioprad must be set equal to 1 when x = 1.0.
!
!          iopang : (integer)
!                 : =0 if angular functions are not computed
!                 : =1 if angular functions of the first kind
!                      are computed
!                 : =2 if angular functions of the first kind and
!                      their first derivatives are computed
!
!          iopnorm: (integer)
!                 : =0 if not scaled. The angular functions have
!                      the same norm as the corresponding associated
!                      Legendre function [i.e., we use the Meixner
!                      -Schafke normalization scheme.] This norm
!                      becomes very large as m becomes large. The
!                      angular functions are computed below as
!                      a characteristic and an exponent to avoid
!                      overflow.
!                 : =1 if angular functions of the first kind
!                      (and their first derivatives if computed)
!                      are scaled by the square root of the
!                      normalization of the corresponding
!                      associated Legendre function. The resulting
!                      scaled angular functions have unity norm.
!                      This is very useful since it removes the
!                      need to calculate a normalization factor
!                      when using the angular function values given
!                      here. It also eliminates any chance for
!                      overflow when the characteristics and exponents
!                      are combined to form the angular functions.
!
!          narg   : number of values of the angular coordinate eta for
!                   which angular functions are calculated (integer)
!
!          arg:     vector containing the values of eta for which
!                   angular functions are desired [real(knd)]
!
!          r1c   :  vectors of length lnum containing the complex
!          r1dc     characteristics for the radial functions of the
!                   first kind r1 and their first derivatives [complex(knd)]
!
!          ir1e   : integer vectors of length lnum containing the
!          ir1de    exponents corresponding to r1c and r1dc
!
!          r2c    : vectors of length lnum containing the complex
!          r2dc     characteristics for the radial functions of the
!                   second kind r2 and their first derivatives [complex(knd)]
!
!          ir2e   : integer vectors of length lnum containing the
!          ir2de    exponents corresponding to r2c and r2dc
!
!          naccr  : integer vector of length lnum containing the estimated
!                   accuracy of the radial functions
!
!          s1c,   : two-dimensional arrays s1c(lnum,narg) and
!          s1dc     s1dc(lnum,narg) that contain narg calculated
!                   complex characteristics for the angular functions
!                   and their first derivatives for each of the lnum
!                   values of l [complex(knd)]
!                   For example, s1c(10,1) is the characteristic
!                   of the angular function for l = m + 10 -1 and
!                   the first value of eta given by arg(1)
!
!          is1e,  : integer arrays is1e(lnum,narg) and is1de(lnum,narg)
!          is1de    containing the exponents corresponding to s1c and
!                   s1dc
!
!          naccs  : two-dimensional array naccs(lnum,narg) containing
!                   narg estimated accuracy values for the angular functions
!                   for each of the lnum values of l
!
!          naccds : two-dimensional array naccds(lnum,narg) containing
!                   narg estimated accuracy values for the first derivatives
!                   of the angular functions for each of the lnum values of c
!
!     Output files
!
!  Cprofcn offers several several output files: Fort.20 and fort.30
!  list the calculated radial and angular functions. Fort.40 and
!  fort.50 are diagnostic files. Fort.60 provides warning whenever the
!  estimated accuracy falls below a specified minimum, currently set
!  equal to 6. Writing to these files is controlled by logicals specified
!  in the module param. False suppresses the file; true enables it.
!  Debug controls fort.30 and fort.40, warn controls fort.60 and output
!  controls fort.20 and fort.30. The logical suffix controls whether the
!  accuracy estimates given in fort.20 are followed by a letter designating
!  how the accuracy was determined. 'w' indicates it is based on the
!  Wronskian and 'e' designates it is based on subtraction errors involved
!  in the calculations. Setting suffix = false suppresses the letter.
!  Information about these files as well as a discussion about accuracy,
!  expansion d coefficients and eigenvalues is given in the readme file.
!
        complex(knd), intent(in)   ::  cc
        real(knd), intent (in)     ::  x1, arg(narg)
        integer, intent (in)       ::  m, lnum, ioprad, iopang, iopnorm, narg
        complex(knd), intent (out) ::  r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                                       s1c(lnum, narg), s1dc(lnum, narg)
        integer, intent (out)      ::  ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                                       is1e(lnum, narg), is1de(lnum, narg), &
                                       naccr(lnum), naccs(lnum, narg), naccds(lnum, narg)

        real(knd) c
!
!     ndec: number of decimal digits available
!     nex:  maximum exponent available
!
        ndec = precision(c)
        nex = range(c) - 1
!
!  open output files
if (output) then
        open(20, file='fort.20')
        open(30, file='fort.30')
end if
if (debug) then
        open(40, file='fort.40')
        open(50, file='fort.50')
end if
if (warn) then
        open(60, file='fort.60')
end if
!   
!  This is where the user sets kindd, the number of bytes available
!  in double precision data for the computer that cprofcn is run on.
!  Similarly kindq, the number of bytes available in quadruple
!  precision data is set here. This allows cprofcn to be run on
!  computers with values for kindd and kindq different than 8 and 16,
!  respectively.
!
5       kindd = 8
        kindq = 16
!
!  Here the minimum desired accuray minacc is set to 8 decimal
!  digits for double precision arithmetic (knd = kindd).
!  The minimum desired accuracy for quadruple precision arithmetic
!  (knd = kindq) is set to 15 decimal digits when aimag(c) <= 20
!  and to 8 decimal digits when aimag(c) > 20. These can be changed
!  if desired. See the readme file.
!
        if(knd == kindd) minacc = 8
        if(knd == kindq .and. aimag(cc) <= 20.0e0_knd) minacc = 15
        if(knd == kindq .and. aimag(cc) > 20.0e0_knd) minacc = 8
!
!  set array dimensions
        c = abs(cc)
        mnum = 1
        mmin = m
        minc = 0
        maxc = 2 * int((x1 + 1.0e0_knd) * (abs(real(cc)) + abs(aimag(cc)))) + 25
        maxm = mmin + minc * (mnum - 1)
        maxint = lnum + 3 * ndec + int(c) + 5
        maxj = maxint + maxm
        maxp = maxint
        maxn = maxp + maxm
        maxn = max(maxn, maxc)
        maxpdr = 4 * ndec + 5
        neta = 993
        ngau = 200
        if(ioprad /= 2) go to 10
        lnump = max(lnum + maxm, 1000)
        if(x1 >= 0.00065e0_knd) maxn = 2 * (lnump * (-18.5e0_knd - 20.0e0_knd* &
                              log10(x1)) + 5 * ndec + 4 * maxm + c + 05000) + maxm + 5
        if(x1 > 0.08e0_knd) maxn = 2 * (lnump * (0.5e0_knd - 3.0e0_knd* &
                              log10(x1)) + 5 * ndec + 4 * maxm + c + 01000) + maxm + 5
        if(x1 > 1.0e0_knd) maxn = 2 * (lnump * 0.5e0_knd + 5 * ndec + 4 * maxm + c+ &
                                 00500) + maxm + 5
        if(x1 <= 0.5e0_knd) maxpdr = maxpdr + int(2.e0_knd * c+ &
                                   100.0e0_knd * x1) + 400
        maxn = max(maxn, maxc)
        maxp = max(maxn, maxp)
        if(x1 < 1.0e-3_knd) ngau = 200 - 50 * int(log10(x1) - 1.0e-30_knd)
        if(x1 < 1.0e-10_knd) ngau = 250 - 50 * int(log10(x1) - 1.0e-30_knd)
        if(x1 < 1.0e-11_knd) ngau = 1200
        if(x1 < 1.0e-12_knd) ngau = 2400
10      maxq = maxint + maxm + maxm
        maxd = maxp / 2 + 1
        maxdr = maxpdr / 2 + 1
        maxp = max(maxp, maxpdr) + 5
        maxlp = lnum + maxm + 5
        maxmp = maxm + 5
        maxt = 1
        jnenmax = 10
        if(iopang /= 0) maxt = narg
!
        call main (mmin, minc, mnum, lnum, cc, ioprad, iopang, iopnorm, &
                   minacc, x1, ngau, narg, arg, maxd, maxdr, maxint, maxj, &
                   maxlp, maxm, maxmp, maxn, maxp, maxpdr, maxq, maxt, neta, &
                   jnenmax, ndec, nex, kindd, kindq, r1c, ir1e, r1dc, ir1de, r2c, &
                   ir2e, r2dc, ir2de, naccr, s1c, is1e, s1dc, is1de, naccs, naccds)
!
        end subroutine
!
!
        subroutine main (mmin, minc, mnum, lnum, cc, ioprad, iopang, iopnorm, &
                         minacc, x1, ngau, narg, barg, maxd, maxdr, maxint, &
                         maxj, maxlp, maxm, maxmp, maxn, maxp, maxpdr, maxq, &
                         maxt, neta, jnenmax, ndec, nex, kindd, kindq, &
                         r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, nar, &
                         s1, is1, s1d, is1d, nas, nads)
!
!  purpose:     to coordinate the calculation of both the prolate
!               spheroidal radial and angular functions and their
!               first derivatives using various algorithms.
!
!     input:    mmin   : minimum desired value of m
!               minc   : increment in m used to compute other values
!               mnum   : number of values of m that are desired
!               lnum   : desired number of values of l = m, m + 1, ...,
!                        m + lnum - 1
!               cc     : complex size parameter c
!               ioprad : equal to 0 if no radial functions are desired;
!                        equal to 1 if only radial functions of the
!                          first kind and their first derivatives are
!                          desired;
!                        equal to 2 if radial functions of both kinds
!                          and their first derivatives are desired
!               iopang : equal to 0 if no angular functions are desired;
!                        equal to 1 if only angular functions of the
!                          first kind are desired;
!                        equal to 2 if angular functions of the first
!                          kind and their first derivatives are desired
!               iopnorm: equal to 0 when the angular functions have
!                        the same norm as the corresponding associated
!                        Legendre functions;
!                        equal to 1 when the angular functions are
!                        scaled by the square root of the normalization
!                        of the corresponding Legendre function, giving
!                        them unity norm
!               minacc : desired minimum accuracy for the radial
!                        functions
!               x1     : x - 1, where x is the radial coordinate, also
!                        called the shape parameter
!               ngau   : order of the Gaussian quadrature to be used in
!                        computing integrals in subroutine pint for use
!                        in subroutine r2int where the integal method
!                        is used to calculate r2 and r2d
!               narg   : number of desired eta angle arguments
!               barg   : vector of desired eta values for which angular
!                        functions are to be computed
!               maxd   : dimension of enr array containing ratios of
!                        the expansion d coefficients
!               maxdr  : dimension of drhor array containing special d
!                        coefficient ratios used in subroutine r2leg
!                        when computing the sum of Legendre functions of
!                        the first kind that appear in the Legendre
!                        function expansion for r2 and r2d
!               maxint : maximum number of terms available for computing
!                        r2 and r2d in the subroutine r2int; dimension
!                        of the arrays of integrals computed in
!                        subroutine pint
!               maxj   : equal to the dimension of the array of ratios
!                        of spherical Bessel functions of the first kind
!                        and the array of ratios of the first derivative
!                        of this Bessel function to the corresponding
!                        Bessel function
!               maxlp  : maximum value desired for l
!               maxm   : maximum value desired for m
!               maxmp  : maxm + 5; dimension of the integer array norme
!                        used in scaling of the Neumann functions in
!                        the integrands in subroutine pint
!               maxn   : dimension of the arrays of Neumann function
!                        ratios used in computing r2 and r2d
!               maxp   : dimension of arrays of Legendre functions of
!                        the first kind used in computing angular
!                        functions, in computing integrands in
!                        subroutine pint and in computing r2 and r2d in
!                        subroutine r2eta
!               maxpdr : dimension of the arrays of ratios of both
!                        Legendre functions of the first kind and their
!                        first derivatives used in the sum of these
!                        functions that contribute to r2 and r2d in
!                        subroutine r2leg
!               maxq   : dimension of arrays of ratios of Legendre
!                        functions of the second kind and ratios of
!                        their first derivatives used in their sum in
!                        subroutine r2leg
!               maxt   : equal to narg if angular functions are
!                        computed where it is the maximum value of the
!                        first index in the arrays of Legendre functions
!                        used in subroutine s1leg;
!                        otherwise equal to 1 to specify the
!                        first index for the Legendre functions used
!                        in the variable eta method for computing r2
!                        and r2d in subroutine r2eta
!               neta   : number of values available for eta in the
!                        variable eta method for calculating r2 and r2d
!                        (subroutine r2eta) and the variable eta method
!                        for calculating r1 and r1d (subroutine r1eta);
!                        set equal to 30
!               jnenmax: number of arrays of ratios of Legendre and
!                        Neumann functions stored as eta is varied in
!                        subroutine r2eta; when x >= 0.05, it is set
!                        equal to 10 so that the previous 10 sets of
!                        ratios are available to use without having to
!                        recalculate them when one of the previous
!                        values for eta is used again for a later value
!                        of l. If x < 0.05 so that r2eta will never be
!                        called, jnenmax is set equal to unity.
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent for real(knd)
!               kindd  : kind value for double precision real data
!               kindq  : kind value for quadruple precision real data
!
!     output:   r1c    : array of lnum values for the characteristic
!                        of the radial function of the first kind
!               ir1e   : array of exponents corresponding to r1c
!               r1dc   : array of lnum values for the characteristic
!                        of the derivative of the radial function of
!                        the first kind)
!               ir1de  : array of exponents corresponding to r1dc
!               r2c    : array of lnum values for the characteristic
!                        of the radial function of the second kind
!               ir2e   : array of exponents corresponding to r2c
!               r2dc   : array of lnum values for the characteristic
!                        of the derivative of the radial function of
!                        the second kind
!               ir2de  : array of exponents corresponding to r2dc
!               nar    : integer vector of length lnum containing the
!                        estimated accuracy of the radial functions
!               s1     : two dimensional array of the characteristics
!                        of the angular functions of the first kind
!                        s1(i,j) is the characteristic for the jth
!                        value of eta and the degree = m + i -1
!                        (s1c in call to subroutine profcn)
!               is1    : array of exponents corresponding to s1 (is1e)
!               s1d    : two dimensional array of the characteristics
!                        of the first derivatives of the angular
!                        functions of the first kind; s1d(i,j) is the
!                        characteristic for the jth value of eta and
!                        the degree = m + i -1 (s1dc in call to
!                        subroutine profcn)
!               is1d   : array of exponents corresponding to s1d
!                       (is1de in call to subroutine profcn)
!               nas    : two-dimensional array nas(lnum,narg) containing
!                        narg estimated accuracy values for the angular functions
!                        for each of the lnum values of l
!               nads   : two-dimensional array nads(lnum,narg) containing
!                        narg estimated accuracy values for the first derivatives
!                        of the angular functions for each of the lnum values of
!
        use param
!
!  scalars
        real(knd) aj1, aj2, ang, apcoef, apcoefn, c, coefme, coefmo, &
                  coefn, coefr1e, coefr1o, dec, etaval, factor, pcoefe, pcoefet, &
                  pcoefn, pcoefo, pdcoefe, pdcoefet, pdcoefo, pi, qdml, qml, rl, &
                  rm, rm2, sgn, ten, termpq, t1, t2, t3, t4, t5, t6, t7, wm, x, xb, &
                  xbninp, x1
        real(knd1) t11, t21, t31, t41, t51, t61, t71
        complex(knd) cc, dmfnorm, dfnorm, dmlmf, dmlf, dmlms, dmlms1, &
                     c2, c4, d01, dneg, eigval, r1ca, r1cb, r1cin, r1dcin, &
                     r1dca, r1dcb, r2ec, r2dec, r2ic, r2dic, r2lc, r2dlc, r2nc, &
                     r2dnc, wronc, wronca, wroncb, wront
        complex(knd1) c21, c41
!
!  arrays with dimension lnum
        integer iqdl(lnum), iql(lnum), ifajo(lnum), ieig(lnum)
        real(knd) qdl(lnum), ql(lnum)
        complex(knd) fajo(lnum), eig(lnum), eigst(lnum)
        complex(knd) r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                     s1(lnum, narg), s1d(lnum, narg)
        integer ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                nar(lnum), is1(lnum, narg), is1d(lnum, narg), nas(lnum, narg), &
                nads(lnum, narg)
!
!  arrays with dimension maxd
        complex(knd) enr(maxd), bliste(maxd), gliste(maxd), blisto(maxd), &
                     glisto(maxd)
        complex(knd1) bliste1(maxd), blisto1(maxd), gliste1(maxd), &
                      glisto1(maxd)
!
!  arrays with dimension maxdr
        complex(knd) drhor(maxdr)
!
!  arrays with dimension maxint
        complex(knd) pint1(maxint), pint2(maxint), pint3(maxint), &
                     pint4(maxint), rpint1(maxint), rpint2(maxint)
!
!  arrays with dimension maxj
        complex(knd) sbesf(maxj), sbesf2(maxj), sbesdf(maxj), &
                     sbesdf2(maxj)
!
!  arrays with dimension maxlp
        integer ibese(maxlp), ibese2(maxlp), ipnormint(maxlp)
        real(knd) pnormint(maxlp)
        complex(knd) sbesdr(maxlp), sbesdr2(maxlp), sbesn(maxlp), &
                     sbesn2(maxlp)
        complex(knd) sneudre(maxlp), sneudrsv(jnenmax, maxlp)
        complex(knd) sneudr(maxlp)
!
!  arrays with dimension maxmp
        complex(knd) enrneg(maxmp)
!
!  arrays with dimension maxn
        integer ineue(maxn), ineuee(maxn), ineuesv(jnenmax, maxn)
        complex(knd) sneufe(maxn), sneudfe(maxn), &
                sneufsv(jnenmax, maxn), sneudfsv(jnenmax, maxn), &
                sneune(maxn), sneunsv(jnenmax, maxn)
        complex(knd) sneun(maxn), sneuf(maxn), sneudf(maxn)
!
!  arrays with dimension given by maxp
        real(knd) alpha(maxp), beta(maxp), coefa(maxp), coefb(maxp), &
                  coefc(maxp), coefd(maxp), coefe(maxp), gamma(maxp), &
                  pdr(maxt, maxp), pdrat(maxt, maxp), pdratt(maxp), &
                  pr(maxt, maxp), prat(maxt, maxp), pratb(maxp), pratt(maxp), &
                  prat1(maxp), pratbsv(jnenmax, maxp), &
                  prattsv(jnenmax, maxp), pdratsv(jnenmax, maxp)
!
!  arrays with dimension maxpdr
        real(knd) prx(maxpdr), pdrx(maxpdr)
!
!  arrays with dimension maxq
        real(knd) qr(maxq), qdr(maxq)
!
!  arrays with dimension maxt
        real(knd) barg(maxt), etainp(maxt), pdnorm(maxt), &
                  pdnorma(maxt), pnorm(maxt), pnorma(maxt), pdtempe(maxt), &
                  pdtempo(maxt), ptempe(maxt), ptempo(maxt), xin(maxt), &
                  xlninp(maxt)
        complex(knd) s1c(maxt), s1dc(maxt)
        integer ipdnorm(maxt), ipdnorma(maxt), ipnorm(maxt), &
                ipnorma(maxt), ipdtempe(maxt), ipdtempo(maxt), &
                iptempe(maxt), iptempo(maxt), is1e(maxt), is1de(maxt), &
                naccs(maxt), naccds(maxt)
!
!  arrays with dimension neta
        real(knd) eta(neta), wmeta2(neta), xbn(neta), xln(neta)
!
!  arrays with dimension ngau
        real(knd) wr(ngau), xr(ngau)
!
!  miscellaneous integer arrays
        integer nees(100), naccsav(100), neeb(jnenmax), limpsv(jnenmax), &
                limnsv(jnenmax), jelimsv(jnenmax)
!
        character chr_w, chr_e
        if (suffix) then
            chr_e = 'e'
            chr_w = 'w'
        else
            chr_e = ' '
            chr_w = ' '
        end if
!
        dec = 10.0e0_knd ** (-ndec - 1)
        ten = 10.0e0_knd
        if(ioprad /= 0) x = x1 + 1.0e0_knd
        jtest = ndec - minacc - 2
        pi = acos(-1.0_knd)
        c = abs(cc)
        c2 = cc * cc
        c4 = c2 * c2
        c21 = c2
        c41 = c4
        nbp = int(2.0e0_knd * (abs(real(cc)) + abs(aimag(cc))) / 3.14e0_knd)
        imax = max(50, nbp) + 5
        lical = imax + imax
!
!  begin loops
          igau = 0
if (debug) then
          if(knd == kindd .and. ioprad /= 0) write(40, 30) x, cc
30        format(1x,'x = ',e23.14,/,1x,'c = ',e23.14, e23.14)
          if(knd == kindq .and. ioprad /= 0) write(40, 35) x, cc
35        format(1x,'x = ',e39.30,/,1x,'c = ',e39.30, e39.30)
end if
          if(ioprad == 2) wront = 1.0e0_knd / (cc * x1 * (x1 + 2.0e0_knd))
40        continue
          maxe = max(100, nbp + nbp) + 10
          ibflag1 = 0
            do 900 mi = 1, mnum
            m = mmin + minc * (mi - 1)
            m2 = m + m
if (debug) then
            if(knd == kindd .and. iopang /= 0) write(50, 50) cc, m
50          format(1x,'c = ',e23.14, e23.14,'; m = ',i5)
            if(knd == kindq .and. iopang /= 0) write(50, 55) cc, m
55          format(1x,'c = ',e39.30, e39.30,'; m = ',i5)
            if(ioprad /= 0) write(40, 60) m
60          format(1x,'m = ',i5)
end if
if (output) then
            if(knd == kindd .and. iopang /= 0) write(30, 65) cc, m
65          format(1x,'c = ',e23.14, e23.14,'; m = ',i5)
            if(knd == kindq .and. iopang /= 0) write(30, 70) cc, m
70          format(1x,'c = ',e39.30, e39.30,'; m = ',i5)
end if
            rm = real(m, knd)
            rm2 = real(m + m)
            icounter = 0
            iopleg = 0
            iopneu = 0
            iopeta = 0
            iopint = 1
            jintm = 0
            iopd = 3
            limcsav = 0
            if(ioprad /= 2) go to 80
            if(x1 <= 0.4e0_knd .and. c <= 10.0e0_knd) iopleg = 1
            if(x1 > 0.4e0_knd .and. c <= 10.0e0_knd) iopneu = 1
            if(x1 <= 0.4e0_knd .and. c <= 15.0e0_knd .and. minacc <= 16) &
                  iopleg = 1
            if(x1 > 0.4e0_knd .and. c <= 20.0e0_knd .and. minacc <= 16) &
                  iopneu = 1
            if(iopleg == 1 .or. iopneu == 1) iopint = 0
            ioppsum = 1
            iopqnsum = 1
            if(m == 0) iopqnsum = 0
            if(knd == kindd) then
              neest = 897
              if(x1 > 0.01e0_knd) neest = 769
              if(x1 > 0.03e0_knd) neest = 705
              if(x1 > 0.04e0_knd) neest = 641
              if(x1 > 0.05e0_knd) neest = 577
              if(x1 > 0.06e0_knd) neest = 513
              if(x1 > 0.07e0_knd) neest = 449
              if(x1 > 0.08e0_knd) neest = 385
              if(x1 > 0.09e0_knd) neest = 1
              end if
              if(knd == kindq) then
              neest = 905
              if(x1 > 0.01e0_knd) neest = 1
              end if
            nee = neest
            jnen = 0
            incnee = 64
            if(knd == kindd .and. x1 < 0.2e0_knd) incnee = 32
            nacceta = 0
            msearch = 0
80          continue
            if(iopang == 0) go to 90
            limps1 = lnum + 3 * ndec + int(c)
            if(limps1 > maxp - 3) limps1 = maxp - 3
            iopd = 0
            if(iopang == 2) iopd = 1
            call pleg_cached(m, limps1, maxp, ndec, nex, limcsav, iopd, barg, narg, &
                      maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, &
                      beta, gamma, coefa, coefb, coefc, coefd, coefe)
            limcsav = limps1
            iopd = 3
90          if(ioprad == 0 .or. mi /= 1 .or. x1 == 0.0e0_knd) go to 100
            limj = lnum + 3 * ndec + int(c) + maxm
            xb = sqrt(x1 * (x1 + 2.0e0_knd))
            call sphbes(cc, xb, limj, maxj, maxlp, ndec, nex, sbesf, sbesdf, &
                        sbesn, ibese, sbesdr)
            if(aimag(cc) /= 0.0e0_knd) call sphbes(cc, x, limj, maxj, maxlp, &
                          ndec, nex, sbesf2, sbesdf2, sbesn2, ibese2, sbesdr2)
100         iflag = 0
            ibflag2 = 0
            legflag = 0
            jflagleg = 0
            naccleg = 0
            legstart = m
            nflag = 0
            lowacc = ndec
            lowtest = minacc
            nacctest = minacc
            naccintp = 0
            nacclegp = 0
            naccneup = 0
            naccr = minacc
            naccrp = minacc
            naccrplp = minacc
            ietacount = 0
            incnflag = 0
            iplflag = 0
            factor = 1.0e0_knd
            ijnet = 0
            iflagpc = 1
            iflagbesb = 0
            iopbesb = 0
            istartr2 = 1
            jeta = 0
            intlim = maxint - 4
            ipint = 0
            iflagnee = 0
            iflagq = 0
            iflagp = 0
            kflag = 0
            max1e = 0
            max1o = 0
            iflagc = 0
            kcor = 0
            jsubms = 0
            if(knd1 == knd) kflag = 2
110         continue
if (output) then
            if(knd == kindd .and. ioprad /= 0) write(20, 120) x, cc, m
120         format(1x,'x = 'e23.14,'; c = ',e23.14, e23.14,'; m = ',i5)
            if(knd == kindq .and. ioprad /= 0) write(20, 130) x, cc, m
130         format(1x,'x = 'e39.30,'; c = ',e39.30, e39.30,'; m = ',i5)
end if
              do 850 li = 1, lnum
              l = m + (li - 1)
if (output) then
              if(iopang /= 0) write(30, 140) l
140           format(1x, i6)
end if
if (debug) then
              if(iopang /= 0) write(50, 150) l
150           format(1x,'l = ',i6)
end if
              ix = l - m - 2 * ((l - m) / 2)
              nacciop = 0
              iopnee = 0
                if(iflagnee == 1) then
                incnee = 8
                  if(knd == kindd) then
                  if(x1 >= 0.05e0_knd) incnee = 16
                  if(x1 >= 0.2e0_knd) incnee = 32
                  end if
                  if(knd == kindq) then
                  if(x1 >= 0.05e0_knd) incnee = 16
                  if(x1 >= 0.1e0_knd) incnee = 32
                  end if
                iflagnee = 2
                end if
              naccetas = minacc
              naccr = -1
              nsav = ndec
              limdrad = 3 * ndec + int(c)
              if(ioprad /= 0 .and. li /= 1) limdrad = jbes + jbes + 20+ &
                                                  int(sqrt(c))
              if(iopint /= 0 .and. li /= 1 .and. jintm > jbes) &
                  limdrad = jintm + jintm + 20 + int(sqrt(c))
              limdbesb = 3 * ndec + int(c) + li
              if(iflagbesb == 1 .and. iopbesb == 0 .and.  &
                  limdbesb > limdrad) limdrad = limdbesb
              limdang = 3 * ndec + int(c)
              if(iopang /= 0 .and. li /= 1) limdang = jang + jang + 20 + &
                                                  int(sqrt(c))
              if(iopang == 0) limd = limdrad
              if(ioprad == 0) limd = limdang
              if(iopang /= 0 .and. ioprad /= 0) limd = max(limdang, limdrad)
              if(li == 1) limmf = limdang
              if(li > 1) limmf = jmf + jmf + 20 + c / 25
              limd = max(limd, limmf)
              if(ioprad /= 2) go to 155
              if(iopleg == 1) limdleg = l - m + 3 * ndec + int(c)
              if(iopleg == 2) limdleg = jleg + jleg + 20 + int(sqrt(c))
              if(iopleg /= 0) limd = max(limd, limdleg)
              limdneu = limd
              lplus = max(l, 1000)
              if(x1 >= 0.00065e0_knd) limdneu = 2 * ((lplus) * (-18.5e0_knd- &
                                        20.0e0_knd * log10(x1)) &
                                        +5 * ndec + 4 * m + c + 01000)
              if(x1 > 0.08e0_knd) limdneu = 2 * ((lplus) * (0.5e0_knd- &
                                        3.0e0_knd * log10(x1))+ &
                                        5 * ndec + 4 * m + c + 01000)
              if(x1 > 1.0e0_knd) limdneu = 2 * ((lplus) * 0.5e0_knd + 5 * ndec + 4 * m + c+ &
                                      00500)
              if(iopneu == 2 .and. naccneu > 0) &
                     limdneu = jneu + jneu + 20 + int(sqrt(c))
              if(iopneu /= 0) limd = max(limd, limdneu)
              limdeta = limd
              if(x1 >= 0.00065e0_knd) limdeta = 2 * ((lplus) * (-18.5e0_knd- &
                           20.0e0_knd * log10(x1)) + 5 * ndec + 4 * m + c + 05000)
              if(x1 > 0.08e0_knd) limdeta = 2 * ((lplus) * (0.5e0_knd - 3.0e0_knd* &
                                      log10(x1)) + 5 * ndec + 4 * m + c + 01000)
              if(x1 > 1.0e0_knd) limdeta = 2 * ((lplus) * 0.5e0_knd + 5 * ndec + 4 * m + c+ &
                                      00500)
              if(iopeta == 3 .and. naccrp > minacc) &
                              limdeta = jeta + jeta + 500 + c / 10
              if(iopeta == 3 .and. naccrp <= minacc) &
                              limdeta = jeta + jeta + 500 + c
              if(iopeta /= 0) limd = max(limd, limdeta)
155           continue
              if(limd > maxp) limd = maxp
              if(ix == 0) limd = max(limd, max1e)
              if(ix == 1) limd = max(limd, max1o)
              if(limd > maxn) limd = maxn - 4
              if(2 * (limd / 2) /= limd) limd = limd - 1
!
!        obtain estimates for the eigenvalues to be used as starting
!        values for the Bouwkamp procedure
              if(l /= m) go to 157
              call geteig(m, cc, nbp, ndec, lnum, maxe, eigst)
157           continue
!
!        use Bouwkamp procedure to obtain accurate eigenvalues
              ndec1 = precision(bliste1(1))
                if(li == 1) then
                jlowe = 1
                limdle = 2
                min1e = 2
                ienre = (3 * ndec + int(cc)) / 2
                if(kflag /= 2) ienre = (3 * ndec1 + int(cc)) / 2
                max1e = min(2 * ienre + 20 + 2 * int(aimag(cc)), limd)
                jlow1e = 1
                end if
                if(li == 2) then
                jlowo = 1
                limdlo = 3
                min1o = 3
                ienro = (3 * ndec + int(cc)) / 2
                if(kflag /= 2) ienro = (3 * ndec1 + int(cc)) / 2
                max1o = min(2 * ienro + 20 + 2 * int(aimag(cc)), limd)
                jlow1o = 1
                end if
                if(li > 2 .and. ix == 0) then
                max1e = min(limd, max1e)
                end if
                if(li > 2 .and. ix == 1) then
                max1o = min(limd, max1o)
                end if
              if(ix == 1) go to 159
!
!  beta coefficients (bliste) for l-m even
              if(limdle > limd) go to 158
              j = jlowe
                do i = limdle, limd, 2
                i2 = i + i
                t1 = i
                t2 = i - 1
                t3 = m2 + i
                t4 = m2 + i - 1
                t5 = m2 + i2 - 1
                t6 = m2 + i2 - 3
                t7 = m2 + i2 + 1
                bliste(j) = c4 * t1 * t2 * t3 * t4 / (t5 * t5 * t6 * t7)
                j = j + 1
                end do
                jsave = j
!
!  gamma coeficients (gliste) for l-m even
              j = jlowe
                do i = limdle - 1, limd + 1, 2
                i2 = i + i
                t1 = m + i - 1
                t2 = m + i
                t3 = m2 * m2 - 1
                t4 = m2 + i2 - 3
                t5 = m2 + i2 + 1
                gliste(j) = t1 * t2 + (0.5e0_knd) * c2 * ((1.0e0_knd) - t3/ &
                          (t4 * t5))
                j = j + 1
                end do
                limdle = limd + 2
                jlowe = jsave
158           continue
              go to 160
159           continue
!
!  beta coefficients (blisto) for l-m odd
              if(limdlo > limd) go to 160
              j = jlowo
                do i = limdlo, limd, 2
                i2 = i + i
                t1 = i
                t2 = i - 1
                t3 = m2 + i
                t4 = m2 + i - 1
                t5 = m2 + i2 - 1
                t6 = m2 + i2 - 3
                t7 = m2 + i2 + 1
                blisto(j) = c4 * t1 * t2 * t3 * t4 / (t5 * t5 * t6 * t7)
                j = j + 1
                end do
                jsave = j
!
!  gamma coeficient (glisto) for l-m odd
              j = jlowo
                do i = limdlo - 1, limd + 1, 2
                i2 = i + i
                t1 = m + i - 1
                t2 = m + i
                t3 = m2 * m2 - 1
                t4 = m2 + i2 - 3
                t5 = m2 + i2 + 1
                glisto(j) = t1 * t2 + (0.5e0_knd) * c2 * ((1.0e0_knd) - t3/ &
                          (t4 * t5))
                j = j + 1
                end do
                jlowo = jsave
                limdlo = limd + 1
160           continue
              if(kflag == 2) go to 163
              if(ix == 1) go to 162
!
!  beta coefficients (bliste1) for l-m even
              if(min1e > max1e) go to 161
              j = jlow1e
                do i = min1e, max1e, 2
                  if(knd1 /= knd) then
                  i2 = i + i
                  t11 = i
                  t21 = i - 1
                  t31 = m2 + i
                  t41 = m2 + i - 1
                  t51 = m2 + i2 - 1
                  t61 = m2 + i2 - 3
                  t71 = m2 + i2 + 1
                  bliste1(j) = c41 * t11 * t21 * t31 * t41 / (t51 * t51 * t61 * t71)
                  else
                  bliste1(j) = bliste(j)
                  end if
                j = j + 1
                end do
                jsave = j
!
!  gamma coeficients (gliste1) for l-m even
              j = jlow1e
              if(li > 2) j = jlow1e + 1
              imin = min1e - 1
              if(li > 2) imin = min1e + 1
                do i = imin, max1e + 1, 2
                  if(knd1 /= knd) then
                  i2 = i + i
                  t11 = m + i - 1
                  t21 = m + i
                  t31 = m2 * m2 - 1
                  t41 = m2 + i2 - 3
                  t51 = m2 + i2 + 1
                  gliste1(j) = t11 * t21 + (0.5e0_knd1) * c21 * ((1.0e0_knd1) - t31/ &
                             (t41 * t51))
                  else
                  gliste1(j) = gliste(j)
                  end if
                j = j + 1
                end do
                jlow1e = jsave
                min1e = max1e + 2
161           continue
              go to 163
162           continue
!
!  beta coefficients (blisto1) for l-m odd
              if(min1o > max1o) go to 163
              j = jlow1o
                do i = min1o, max1o, 2
                  if(knd1 /= knd) then
                  i2 = i + i
                  t11 = i
                  t21 = i - 1
                  t31 = m2 + i
                  t41 = m2 + i - 1
                  t51 = m2 + i2 - 1
                  t61 = m2 + i2 - 3
                  t71 = m2 + i2 + 1
                  blisto1(j) = c41 * t11 * t21 * t31 * t41 / (t51 * t51 * t61 * t71)
                  else
                  blisto1(j) = blisto(j)
                  end if
                j = j + 1
                end do
                jsave = j
!
!  gamma coeficient (glisto1) for l-m odd
              j = jlow1o
              if(li > 2) j = jlow1o + 1
              imin = min1o - 1
              if(li > 2) imin = min1o + 1
                do i = imin, max1o + 1, 2
                  if(knd1 /= knd) then
                  i2 = i + i
                  t11 = m + i - 1
                  t21 = m + i
                  t31 = m2 * m2 - 1
                  t41 = m2 + i2 - 3
                  t51 = m2 + i2 + 1
                  glisto1(j) = t11 * t21 + (0.5e0_knd1) * c21 * ((1.0e0_knd1) - t31/ &
                          (t41 * t51))
                  else
                  glisto1(j) = glisto(j)
                  end if
                j = j + 1
                end do
              min1o = max1o + 1
              jlow1o = jsave
163           continue
              if(li <= lical) eigval = eigst(li)
              if(li > lical) eigval = 4.0e0_knd * eig(li - 1)- &
                    6.0e0_knd * eig(li - 2) + 4.0e0_knd * eig(li - 3) - eig(li - 4)
              itestm = ndec
              idigc = ndec
              if(ix == 0) call conver (l, m, cc, limd, maxd, bliste, gliste, &
                                       bliste1, gliste1, ioprad, ienre, &
                                       kindd, kindq, ndec, eigval, enr, &
                                       idigc, itestm, kflag)
              if(ix == 1) call conver (l, m, cc, limd, maxd, blisto, glisto, &
                                       blisto1, glisto1, ioprad, ienro, &
                                       kindd, kindq, ndec, eigval, enr, &
                                       idigc, itestm, kflag)
              if(li == lnum + 1) write(40,*) eigval
              eig(li) = eigval
              ieig(li) = idigc
              if(ix == 0) max1e = 2 * ienre + 20 + 2 * int(aimag(cc))
              if(ix == 1) max1o = 2 * ienro + 20 + 2 * int(aimag(cc))
177           call dnorm (l, m, cc, limd, maxd, ndec, nex, ioprad, enr, sgn, d01, &
                          id01, dmfnorm, idmfe, dmlmf, idmlmfe, dfnorm, idfe, dmlf, &
                          idmlfe, jmf, nsubmf, jfla, nsubf)
              jmf = max(jmf, jfla)
              jsub = max(nsubmf, nsubf)
              if(ioprad == 0) go to 720
              if(l == m .and. nsubmf > jtest) iopint = 1
              if(l == m .and. nsubmf > jtest .and. iopneu /= 0) iopneu = 0
              if(l == m .and. nsubmf > jtest .and. iopleg /= 0) iopleg = 0
!
!  determine prolate radial functions of the first kind
!     calculation of r1 and r1d when x = 1.0
              if(x1 == 0.0e0_knd .and. m == 0) then
               if(l == 0) coefr1e = 1.0e0_knd
               if(l == 1) coefr1o = 1.0e0_knd
               rl = real(l, knd)
               if(ix == 0) then
                if(l > 0) coefr1e = coefr1e * rl / (rl - 1.0e0_knd)
                r1c(li) = coefr1e * d01 * dmlf
                iterm = int(log10(abs(r1c(li))))
                r1c(li) = r1c(li) * (10.0e0_knd ** (-iterm))
                ir1e(li) = id01 + idmlfe + iterm
                r1dc(li) = cc * cc * coefr1e * d01 * dmlf * ((enr(1) / 15.0e0_knd)- &
                           (1.0e0_knd / 3.0e0_knd))
                iterm = int(log10(abs(r1dc(li))))
                r1dc(li) = r1dc(li) * (10.0e0_knd ** (-iterm))
                ir1de(li) = id01 + idmlfe + iterm
               end if
               if(ix == 1) then
                if(l > 1) coefr1o = coefr1o * (rl - 1.0e0_knd) / rl
                r1c(li) = cc * coefr1o * d01 * dmlf / 3.0e0_knd
                iterm = int(log10(abs(r1c(li))))
                r1c(li) = r1c(li) * (10.0e0_knd ** (-iterm))
                ir1e(li) = id01 + idmlfe + iterm
                r1dc(li) = cc * cc * cc * coefr1o * d01 * dmlf * ((enr(1) / 35.0e0_knd)- &
                           (1.0e0_knd / 15.0e0_knd) + (1.0e0_knd / (3.0e0_knd * cc * cc)))
                iterm = int(log10(abs(r1dc(li))))
                r1dc(li) = r1dc(li) * (10.0e0_knd ** (-iterm))
                ir1de(li) = id01 + idmlfe + iterm
               end if
               if(abs(r1c(li)) < 1.0e0_knd) then
                r1c(li) = r1c(li) * 10.0e0_knd
                ir1e(li) = ir1e(li) - 1
               end if
               if(abs(r1dc(li)) < 1.0e0_knd) then
                r1dc(li) = r1dc(li) * 10.0e0_knd
                ir1de(li) = ir1de(li) - 1
               end if
               naccr1 = ndec - nsubf - 2
               if(naccr1 < 0) naccr1 = 0
              end if
              if(x1 == 0.0e0_knd .and. m /= 0) then
               r1c(li) = (0.0e0_knd, 0.0e0_knd)
               ir1e(li) = 0
               r1dc(li) = (0.0e0_knd, 0.0e0_knd)
               ir1de(li) = 0
               naccr1 = ndec
              end if
              if(x1 == 0.0e0_knd) then
if (debug) then
              write(40, 178) naccr1
              if(knd == kindd) write(40, 180) r1c(li), ir1e(li), r1dc(li), &
                     ir1de(li)
              if(knd == kindq) write(40, 181) r1c(li), ir1e(li), r1dc(li), &
                     ir1de(li)
end if
               go to 185
              end if
!     calculation of r1 and r1d when x /= 1.0
              if(li == 1) limr1 = 3 * ndec + int(c)
              if(li /= 1) limr1 = jbesa + jbesa + 20 + c / 25
              call r1besa(l, m, cc, x1, limr1, maxd, enr, maxj, maxlp, ndec, nex, &
                   iflag, sbesf, sbesdf, sbesn, ibese, sbesdr, d01, id01, &
                   dfnorm, idfe, r1ca, ir1ea, r1dca, ir1dea, jbesa, &
                   factor, nsubr1a)
              naccr1a = idigc - max(nsubr1a + 1, nsubf)
              if(naccr1a == ndec) naccr1a = ndec - 1
              if(naccr1a < 0) naccr1a = 0
if (debug) then
              write(40, 178) naccr1a
178           format(16x,'estimated accuracy of r1 and r1d is ',i3, &
                     ' digits.')
              if(knd == kindd) write(40, 180) r1ca, ir1ea, r1dca, ir1dea
              if(knd == kindq) write(40, 181) r1ca, ir1ea, r1dca, ir1dea
180           format(10x,'r1 = ',f17.14, 1x, f17.14, i6, 2x,'r1d = ', &
                      f17.14, 1x, f17.14, i6)
181           format(10x,'r1 = ',f33.30, 1x, f33.30, i6,/,10x,'r1d = ', &
                      f33.30, 1x, f33.30, i6)
end if
              ichoicer1 = 1
              iflagbesb = 0
              r1c(li) = r1ca
              ir1e(li) = ir1ea
              r1dc(li) = r1dca
              ir1de(li) = ir1dea
              naccr1 = naccr1a
              if(naccr1a >= minacc .or. min(idigc, itestm) - nsubmf <=  &
                  naccr1a) then
              iopbesb = 0
              jbesb = 0
              go to 185
              end if
              if(iflagpc == 0) go to 184
              prat1(1) = 1.0e0_knd
              prat1(2) = rm2 + 1.0e0_knd
                do jp = 3, limj - maxm
                aj1 = real(jp - 1)
                aj2 = real(jp - 2)
                prat1(jp) = (rm2 + aj1) * (rm2 + aj2) / (aj1 * aj2)
                end do
              pcoefn = x1 * (x1 + 2.0e0_knd) / (x * x)
              apcoefn = (rm / 2.0e0_knd) * log10(pcoefn)
              ipcoefn = int(apcoefn)
              pcoefn = 10.0e0_knd ** (apcoefn - ipcoefn)
              iflagpc = 0
184           continue
              if(li == 1 .or. iopbesb == 0) limr1 = l + 3 * ndec + int(c)
              if(li /= 1 .and. iopbesb == 1) limr1 = jbesb + jbesb + 20 + c / 25
              iopbesb = 1
              if(aimag(cc) /= 0.0e0_knd) call r1besb(l, m, cc, x1, limr1, &
                   maxd, maxlp, ndec, maxj, maxp, enr, sbesf2, sbesn2, ibese2, &
                   sbesdf2, sbesdr2, prat1, pcoefn, ipcoefn, dmfnorm, idmfe, &
                   r1cb, ir1eb, r1dcb, ir1deb, jbesb, nsubr1b)
if (debug) then
              if(knd == kindd) write(40, 180) r1cb, ir1eb, r1dcb, ir1deb
              if(knd == kindq) write(40, 181) r1cb, ir1eb, r1dcb, ir1deb
end if
              naccr1b = idigc - max(nsubr1b + 1, nsubmf)
              if(naccr1b == ndec) naccr1b = ndec - 1
              if(naccr1b < 0) naccr1b = 0
if (debug) then
              write(40, 178) naccr1b
end if
              if(naccr1a >= naccr1b) go to 185
              naccr1 = naccr1b
              r1c(li) = r1cb
              ir1e(li) = ir1eb
              r1dc(li) = r1dcb
              ir1de(li) = ir1deb
              ichoicer1 = 2
185           jbes = max(jbesa, jbesb)
if (output) then
              if(ioprad == 1) write(20, 186) l, r1c(li), ir1e(li), r1dc(li), &
                              ir1de(li), naccr1
186                           format(1x, i5, 2x, 2(f17.14, 1x, f17.14, i5, 2x), i4)
end if
                if(ioprad == 1 .and. naccr1 < 6) then
if (warn) then
                write(60, 715) naccr1, m, l, x, cc
end if
                end if
                if(ioprad == 1) then
                lisave = li
                go to 720
                end if
!
!  determine prolate radial functions of the second kind
!
!  check to see if r1 is large and can use r2 = i r1 to obtain
!  values for r2 accurate to naccrpl digits (or naccr1 digits
!  if naccr1 is less than naccrpl)
              naccrpl = ir1e(li) + ir1de(li) + int(log10(abs(r1c(li)* &
                      r1dc(li)) * c * (x * x1)))
              if(naccrpl < 0) naccrpl = 0
              if(naccrpl > ndec) naccrpl = ndec
                if(min(naccrpl, naccr1) > 0) then
                naccr = min(naccrpl, naccr1)
                r2c(li) = cmplx(-aimag(r1c(li)), real(r1c(li)), knd)
                ir2e(li) = ir1e(li)
                r2dc(li) = cmplx(-aimag(r1dc(li)), real(r1dc(li)), knd)
                ir2de(li) = ir1de(li)
                nacciop = 1
if (debug) then
                if(naccr > 1) write(40, 187) naccr
187             format(8x,'values for r2 and r2dc accurate to ',i2, &
                       ' digits are given by ir1 and ir1d')
end if
                end if
              ir2est = int(log10(abs(wront))) - ir1de(li) + 1
              if(naccrpl > 1) ir2est = ir1e(li)
                if(naccr >= minacc) then
                legstart = l
                if(iopint /= 0) iopint = 1
                go to 680
                else
                  if(naccr1 == naccr .and. naccr > 5) then
                  if(iopint /= 0) iopint = 1
                  go to 680
                  end if
                end if
!
!  calculation using integration technique
              if(ioprad /= 2) go to 680
              if(iopint == 0) go to 230
                if(li > intlim - 5) then
                iopint = 0
                go to 230
                end if
              if(iopint == 2) go to 190
              limint = lnum + 3 * ndec + int(c)
              if(igau == 0) call gauss_cached(ngau, ndec, xr, wr)
              igau = 1
              if(ipint == 1) go to 190
              ngqs = 10
              if(c > 2000.0e0_knd) ngqs = ngqs * (c / 2000.0e0_knd)* &
                                           (c / 2000.0e0_knd)
              lim1max = 4 * int(x * abs(real(cc)) + x * abs(aimag(cc))) + 25
              call pint(cc, m, lnum, x1, limint, maxint, maxlp, ndec, nex, &
                        lim1max, wr, xr, ngau, ngqs, intlim, rpint1, rpint2, &
                        pint1, pint2, pint3, pint4, norme, pnormint, &
                        ipnormint, coefme, coefmo)
              ipint = 1
190           continue
              if(iopint == 1) limint = 3 * ndec + int(c) + li
              if(iopint == 2) limint = jintm + jintm + 20 + int(sqrt(c))
              if(limint > maxint - 5) limint = maxint - 5
              call r2int(l, m, cc, x, limint, ndec, maxd, enr, d01, id01, &
                         maxint, maxmp, nex, maxlp, rpint1, rpint2, &
                         pint1, pint2, pint3, pint4, norme, pnormint, &
                         ipnormint, coefme, coefmo, r2ic, ir2ie, &
                         r2dic, ir2die, jint, coefn, icoefn, iflagc)
              iopint = 2
              if(jint > jintm) jintm = jint
              if(iopint == 3 .and. jint < jintm) jint = jintm
                if(ir2ie > ir2est + ndec) then
                iopint = 0
                naccint = 0
                naccintp = 0
                naccout = 0
                nacccor = 0
                r2ic = (0.0e0_knd, 0.0e0_knd)
                ir2ie = 0
                r2dic = (0.0e0_knd, 0.0e0_knd)
                ir2die = 0
                istartr2 = 1
                go to 201
                end if
              wronca = r1c(li) * r2dic * ten ** (ir1e(li) + ir2die)
              wroncb = r2ic * r1dc(li) * ten ** (ir2ie + ir1de(li))
              wronc = wronca - wroncb
              naccint = -int(log10(abs((wronc - wront) / wront) + dec)- &
                         0.5e0_knd)
              if(naccint < 0) naccint = 0
              if(naccint > ndec-1) naccint = ndec-1
              naccinto = naccint
              nacccor = -int(log10(abs((wronca - wroncb) / wronca) + dec))
              if(nacccor < 0) nacccor = 0
              nacccor = min(nacccor, naccrpl + 1)
              naccint = min(naccint + nacccor, naccr1)
              naccout = naccint
              if(naccout < naccr) go to 200
              naccr = naccout
              r2c(li) = r2ic
              ir2e(li) = ir2ie
              r2dc(li) = r2dic
              ir2de(li) = ir2die
              nacciop = 0
200           continue
              istartr2 = 1
              if(naccint > minacc .and. (naccintp > minacc .or. li == 1)) &
                  then
              iopleg = 0
              iopneu = 0
              iopeta = 0
              istartr2 = 0
              end if
              if(naccout >= minacc .and. naccintp >= minacc) then
              iopneu = 0
              iopeta = 0
              end if
              if(naccint >= minacc .and. ndec - jsub <= naccint .and.  &
                      iopleg /= 0) iopleg = 0
              if(naccint == 0 .and. naccintp == 0) iopint = 0
              naccintp = naccout
201           if(naccint < minacc .and. x1 <= 0.1e0_knd .and.  &
                 iopleg == 0 .and. l >= legstart .and.  &
                 jsub <= ndec - naccint) iopleg = 1
              nacccor = min(nacccor, naccint - naccinto, naccint)
              if(nacccor < 0) nacccor = 0
if (debug) then
                if (knd == kindd) then
                if(nacccor == 0) write(40, 205) naccout, r2ic, ir2ie, &
                                                r2dic, ir2die
205             format(15x,'accuracy =',i3, &
                       ' decimal digits.'/,10x,'r2 = ', f17.14, 1x, f17.14, &
                       i6, 2x,'r2d = ',f17.14, 1x, f17.14, i6)
                if(nacccor > 0) write(40, 210) naccout, nacccor, r2ic, &
                                                ir2ie, r2dic, ir2die
210             format(15x,'accuracy =',i3,' decimal digits, adjusted ', &
                       'for',i3,' subtraction digits in forming ', &
                       'Wronskian',/,10x,'r2 = ',f17.14, 1x, f17.14, i6, &
                       2x,'r2d = ',f17.14, 1x, f17.14, i6)
                end if
                if (knd == kindq) then
                if(nacccor == 0) write(40, 215) naccout, r2ic, ir2ie, &
                                                r2dic, ir2die
215             format(15x,'accuracy =',i3, &
                       ' decimal digits.'/,10x,'r2 = ', f33.30, 1x, f33.30, &
                       i6,/,10x,'r2d = ',f33.30, 1x, f33.30, i6)
                if(nacccor > 0) write(40, 220) naccout, nacccor, r2ic, &
                                                ir2ie, r2dic, ir2die
220             format(15x,'accuracy =',i3,' decimal digits, adjusted ', &
                       'for',i3,' subtraction digits in forming ', &
                       'Wronskian',/,10x,'r2 = ',f33.30, 1x, f33.30, i6,/, &
                       10x,'r2d = ',f33.30, 1x, f33.30, i6)
                end if
end if
230           continue
!
!  calculation using Legendre expansion and joining factor
              if(l >= legstart .and. naccr < minacc .and. naccrp < minacc &
                  .and. ndec - jsub >= naccrp .and. iopleg == 0) iopleg = 1
              if(l >= legstart .and. iopeta /= 0 .and. ndec - jsub >= naccrp &
                  .and. x1 < 0.1e0_knd .and. iopleg == 0) iopleg = 1
              if(iopleg == 0) go to 360
              if(iopleg == 2 .or. jflagleg == 1) go to 310
              jflagleg = 1
              limdr = c + ndec + 50.0e0_knd * x1 + 200
              if(limdr > maxdr - 2) limdr = maxdr - 2
              if(ioppsum == 0) go to 250
              xin(1) = x
              limpleg = limdr + limdr
              call pleg_cached(m, limpleg, maxp, ndec, nex, limcsav, iopd, xin, 1, maxt, &
                        prat, pdrat, pdnorma, ipdnorma, pnorma, ipnorma, &
                        alpha, beta, gamma, coefa, coefb, coefc, coefd, coefe)
              limcsav = max(limcsav, limpleg)
                do jj = 1, limpleg
                prx(jj) = prat(1, jj)
                pdrx(jj) = pdrat(1, jj)
                end do
250           limq = lnum + 3 * ndec + int(c)
              call qleg_cached(m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, &
                        iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
              fajo(1) = cc / (rm2 - 1.0e0_knd)
              ifajo(1) = 0
              if(m == 0) go to 280
                do im = 1, m
                fajo(1) = fajo(1) * real(im + im, knd) / cc
                if(abs(fajo(1)) < 1.0e+10_knd) go to 260
                fajo(1) = fajo(1) * (1.0e-10_knd)
                ifajo(1) = ifajo(1) + 10
260             continue
                if(abs(fajo(1)) > 1.0e-10_knd) go to 270
                fajo(1) = fajo(1) * (1.0e+10_knd)
                ifajo(1) = ifajo(1) - 10
270             continue
                end do
280           continue
              fajo(2) = -cc * fajo(1) / (rm2 - 3.0e0_knd)
              ifajo(2) = ifajo(1)
                do jl = 3, lnum - 1, 2
                fajo(jl) = fajo(jl - 2) * (real(jl + m + m - 1, knd) &
                         /real(jl - 2, knd))
                ifajo(jl) = ifajo(jl - 2)
                if(abs(fajo(jl)) < 1.0e10_knd) go to 290
                fajo(jl) = fajo(jl) * 1.0e-10_knd
                ifajo(jl) = ifajo(jl) + 10
290             fajo(jl + 1) = fajo(jl - 1) * (real(jl + m + m - 1, knd) &
                           /real(jl, knd))
                ifajo(jl + 1) = ifajo(jl - 1)
                if(abs(fajo(jl + 1)) < 1.0e10_knd) go to 300
                fajo(jl + 1) = fajo(jl + 1) * 1.0e-10_knd
                ifajo(jl + 1) = ifajo(jl + 1) + 10
300             end do
              if(2 * (lnum / 2) == lnum .or. lnum == 2) go to 310
              fajo(lnum) = fajo(lnum - 2) * real(lnum + m + m - 1, knd) / real(lnum - 2, knd)
              ifajo(lnum) = ifajo(lnum - 2)
310           continue
              limleg = l - m + 3 * ndec + int(c)
              limdr = c + ndec + 50.0e0_knd * x1 + 200
              if(iopleg == 2) limleg = jleg + jleg + 20 + int(sqrt(c))
              if(iopleg == 2) limdr = jlegp + 10 + int(0.5e0_knd * sqrt(c))
              if(limdr > maxdr) limdr = maxdr
              call dalt(l, m, cc, limdr, maxdr, maxmp, ndec, nex, ioppsum, &
                        eigval, enrneg, drhor, dneg, idneg, nsdneg, nsdrho)
              r1cin = r1c(li)
              ir1ein = ir1e(li)
              r1dcin = r1dc(li)
              ir1dein = ir1de(li)
              call r2leg(l, m, cc, x1, lnum, limleg, limdr, maxd, maxmp, ndec, &
                         nex, idigc, itestm, maxpdr, maxdr, maxq, enr, enrneg, &
                         drhor, nsdrho, d01, id01, dneg, idneg, nsdneg, dfnorm, &
                         idfe, nsubf, dmfnorm, idmfe, nsubmf, prx, pdrx, qdr, &
                         qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, iql, fajo, &
                         ifajo, termpq, itermpq, ioppsum, iopqnsum, r1cin, &
                         ir1ein, r1dcin, ir1dein, naccr1, wront, minacc, &
                         naccrpl, naccr, r2lc, ir2le, r2dlc, ir2dle, jleg, &
                         jlegp, naccleg, nacccor, jflagl, iflagq, iflagp)
              naccout = naccleg
              if(naccout <= naccr) go to 320
              nacciop = 0
              if(jflagl == 1) nacciop = 1
              naccr = naccout
              r2c(li) = r2lc
              ir2e(li) = ir2le
              r2dc(li) = r2dlc
              ir2de(li) = ir2dle
320           continue
                if(naccout >= naccrp) then
                iopleg = 2
                else
                  if(iopleg == 1 .and. l /= m .and. naccrp /=  &
                       min(naccrplp, naccr1p)) then
                  iopleg = 0
                  legstart = l + naccrp - naccout
                  end if
                end if
              if(naccout > minacc .and. nacclegp > minacc) then
              iopleg = 2
              iopneu = 0
              iopeta = 0
              iopint = 0
              end if
              nacclegp = naccleg
if (debug) then
                if(knd == kindd) then
                if(nacccor == 0) write(40, 330) naccout, r2lc, ir2le, &
                                                r2dlc, ir2dle
330             format(15x,'accuracy =',i3, &
                       ' decimal digits.'/,10x,'r2 = ', f17.14, 1x, f17.14, &
                       i6, 2x,'r2d = ',f17.14, 1x, f17.14, i6)
                if(nacccor > 0) write(40, 340) naccout, nacccor, r2lc, &
                                                ir2le, r2dlc, ir2dle
340             format(15x,'accuracy =',i3, &
                      ' decimal digits, adjusted for',i3,' subtraction', &
                      ' digits in forming Wronskian',/,10x,'r2 = ', &
                       f17.14, 1x, f17.14, i6, 2x,'r2d = ',f17.14, 1x, f17.14, &
                       i6)
                end if
                if(knd == kindq) then
                if(nacccor == 0) write(40, 350) naccout, r2lc, ir2le, &
                                                r2dlc, ir2dle
350             format(15x,'accuracy =',i3, &
                       ' decimal digits.'/,10x,'r2 = ', f33.30, 1x, f33.30, &
                       i6,/,10x,'r2d = ',f33.30, 1x, f33.30, i6)
                if(nacccor > 0) write(40, 355) naccout, nacccor, r2lc, &
                                                ir2le, r2dlc, ir2dle
355             format(15x,'accuracy =',i3, &
                      ' decimal digits, adjusted for',i3,' subtraction', &
                      ' digits in forming Wronskian',/,10x,'r2 = ', &
                       f33.30, 1x, f33.30, i6,/,10x,'r2d = ',f33.30, 1x, f33.30, &
                       i6)
                end if
end if
360           continue
!
!  calculation using conventional Neumann expansion (eta=1)
              if(iopneu == 0) go to 420
              if(iopneu == 2) go to 380
              if(ibflag1 == 1) go to 370
              ibflag1 = 1
              lnump = max(lnum + maxm, 1000)
              limn = 2 * (lnump * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
                    5 * ndec + 4 * m + c + 01000) + maxm
              if(x1 > 0.08e0_knd) limn = 2 * (lnump * (0.5e0_knd - 3.0e0_knd* &
                                     log10(x1)) + 5 * ndec + 4 * m + c + 01000) + maxm
              if(x1 > 1.0e0_knd) limn = 2 * (lnump * 0.5e0_knd + 5 * ndec + 4 * m + c+ &
                                       00500) + maxm
              limbesf = 4 * int(real(cc * x) + abs(aimag(cc * x))) + 25
              call sphneu(cc, x, limn, maxn, maxlp, ndec, nex, limbesf, sneuf, &
                          sneun, ineue, sneudf, sneudr)
370           if(ibflag2 == 1) go to 380
              ibflag2 = 1
              lp = max(lnum + m, 1000)
              limp = 2 * (lp * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
                   5 * ndec + 4 * m + c + 01000)
              if(x1 > 0.08e0_knd) limp = 2 * (lp * (0.5e0_knd - 3.0e0_knd * log10(x1))+ &
                                   5 * ndec + 4 * m + c + 01000)
              if(x1 > 1.0e0_knd) limp = 2 * (lp * 0.5e0_knd + 5 * ndec + 4 * m + c + 00500)
              if(limp > maxp) limp = maxp
              prat1(1) = 1.0e0_knd
              prat1(2) = rm2 + 1.0e0_knd
                do jp = 3, limp
                aj1 = real(jp - 1, knd)
                aj2 = real(jp - 2, knd)
                prat1(jp) = (rm2 + aj1) * (rm2 + aj2) / (aj1 * aj2)
                end do
              pcoefn = x1 * (x1 + 2.0e0_knd) / (x * x)
              apcoefn = (rm / 2.0e0_knd) * log10(pcoefn)
              ipcoefn = int(apcoefn)
              pcoefn = 10.0e0_knd ** (apcoefn - ipcoefn)
380           continue
              lplus = max(l, 1000)
              limneu = 2 * ((lplus) * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
                     5 * ndec + 4 * m + c + 00500)
              if(x1 >= 0.08e0_knd) limneu = 2 * ((lplus) * (0.5e0_knd- &
                                     log10(x1)) + 5 * ndec + 4 * m + c + 00500)
              if(x1 > 1.0e0_knd) limneu = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                                         4 * m + c + 00500)
              if(iopneu == 2 .and. naccneu > 0) limneu = jneu + jneu + 20+ &
                                                      int(sqrt(c))
              if(limneu > limp - 2) limneu = limp - 2
              r1dcin = r1dc(li)
              ir1dein = ir1de(li)
              call r2neu(l, m, cc, x1, limneu, maxd, maxlp, ndec, nex, maxn, &
                         maxp, minacc, enr, sneuf, sneun, ineue, sneudf, &
                         sneudr, prat1, pcoefn, ipcoefn, dmfnorm, idmfe, &
                         r1dcin, ir1dein, naccrpl, r2nc, ir2ne, r2dnc, &
                         ir2dne, jneu)
              wronca = r1c(li) * r2dnc * ten ** (ir1e(li) + ir2dne)
              wroncb = r2nc * r1dc(li) * ten ** (ir2ne + ir1de(li))
              wronc = wronca - wroncb
              naccneu = -int(log10(abs((wronc - wront) / wront) + dec)- &
                         0.5e0_knd)
              if(naccneu < 0) naccneu = 0
              if(naccneu > ndec-1) naccneu = ndec-1
              naccneuo = naccneu
              nacccor = -int(log10(abs((wronca - wroncb) / wronca) + dec))
              nacccor = min(nacccor, naccrpl + 1)
              if(nacccor < 0) nacccor = 0
              if(naccneu > 0) naccneu = min(naccneu + nacccor, naccr1, &
                                       ndec - nsubmf)
              if(naccneu == 0 .and. naccrpl > 0) then
              imat = -int(log10(abs((r1c(li) * (0.0e0_knd, 1.0e0_knd) - r2nc* &
              ten ** (ir2ne - ir1e(li))) / r1c(li)) + dec))
              imatd = -int(log10(abs((r1dc(li) * (0.0e0_knd, 1.0e0_knd)- &
                    r2dnc * ten ** (ir2dne - ir1de(li))) / r1dc(li)) + dec))
              naccneu = min(imat, imatd, naccrpl, naccr1, ndec - nsubmf)
              if(naccneu < 0) naccneu = 0
              end if
385           naccout = naccneu
              if(naccneup - naccneu > 8) naccneu = naccneup
              naccneup = naccneu
              if(naccout <= naccr) go to 390
              naccr = naccout
              r2c(li) = r2nc
              ir2e(li) = ir2ne
              r2dc(li) = r2dnc
              ir2de(li) = ir2dne
              nacciop = 0
390           continue
              if(naccneu >= minacc) then
              iopneu = 2
              iopeta = 0
              end if
              if(naccneu > minacc) then
              nflag = 1
              if(li > nbp) iopint = 0
              end if
              if(iopeta == 0 .and. naccr < minacc) then
              nflag = 0
              iopneu0 = 0
              iopeta = 1
              if(nacccor == 0) nee = max(nee, 993 - incnee * (minacc - naccneu &
                                   +1))
              end if
              nacccor = min(nacccor, naccneu - naccneuo, naccneu)
              if(nacccor < 0) nacccor = 0
if (debug) then
                if(knd == kindd) then
                if(nacccor == 0) write(40, 400) naccout, r2nc, ir2ne, &
                                                r2dnc, ir2dne
400             format(15x,'Wronskian accuracy =',i3, &
                      ' decimal digits.'/,10x,'r2 = ', f17.14, 1x, f17.14, &
                      i6, 2x,'r2d = ',f17.14, 1x, f17.14, i6)
                if(nacccor > 0) write(40, 405) naccout, nacccor, r2nc, &
                                                ir2ne, r2dnc, ir2dne
405             format(15x,'accuracy =',i3, &
                      ' decimal digits, adjusted for',i3,' subtraction', &
                      ' digits in forming Wronskian',/,10x,'r2 = ', &
                      f17.14, 1x, f17.14, 2x, i6, 2x,'r2d = ',f17.14, 1x, f17.14, &
                      2x, i6)
                end if
                if(knd == kindq) then
                if(nacccor == 0) write(40, 410) naccout, r2nc, ir2ne, &
                                                r2dnc, ir2dne
410             format(15x,'Wronskian accuracy =',i3, &
                      ' decimal digits.'/,10x,'r2 = ', f33.30, 1x, f33.30, &
                      i6,/,10x,'r2d = ',f33.30, 1x, f33.30, i6)
                if(nacccor > 0) write(40, 415) naccout, nacccor, r2nc, &
                                                ir2ne, r2dnc, ir2dne
415             format(15x,'accuracy =',i3, &
                      ' decimal digits, adjusted for',i3,' subtraction', &
                      ' digits in forming Wronskian',/,10x,'r2 = ', &
                      f33.30, 1x, f33.30, i6,/,10x,'r2d = ',f33.30, 1x, f33.30, &
                      i6)
                end if
end if
420           continue
!
!  calculation using the variable eta expansion
              if(iopeta == 0 .or. iopeta == 4) go to 670
                do 430 inn = 1, 100
                nees(inn) = 0
430             naccsav(inn) = 0
              inen = 0
              neemark = nee
              naccmax = 0
              neemax = nee
              naccnmax = 0
              netatry = 1
              naccdp = 0
              kounte = 0
              if(iopeta > 1) go to 440
              kounter = 0
                if(ijnet == 0) then
                  do jnet = 1, neta
                  ang = (neta + 1 - jnet) * pi * 0.5e0_knd / (neta + 1)
                  eta(jnet) = cos(ang)
                  wmeta2(jnet) = 2.0e0_knd * (1.0e0_knd + eta(jnet))* &
                               (sin(0.5e0_knd * ang) ** 2)
                  xbn(jnet) = sqrt(x1 * (x1 + 2.0e0_knd) + eta(jnet) * eta(jnet))
                  xln(jnet) = eta(jnet) * (x1 + 1.0e0_knd) / xbn(jnet)
                  end do
                ijnet = 1
                end if
              iopeta = 2
440           if(iopeta == 3) go to 540
              etaval = eta(nee)
450           xbninp = xbn(nee)
              netainp = 1
              etainp(1) = eta(nee)
              xlninp(1) = xln(nee)
              lplus = max(l, 1000)
              limn = 2 * ((lplus) * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
                   5 * ndec + 10 * incnee + 4 * m + c + 05000) + m
              if(x1 > 0.08e0_knd) limn = 2 * ((lplus) * (0.5e0_knd- &
                 3.0e0_knd * log10(x1)) + 10 * ndec + 10 * incnee+4 * m + c + 01000) + m
              if(x1 > 1.0e0_knd) limn = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                                       10 * incnee + 4 * m + c + 00500) + m
              if(limn > maxn - 2) limn = maxn - 2
              limp = limn - m
              if(jnen == 0) go to 510
              jnenlim = jnen
              if(jnen > jnenmax) jnenlim = jnenmax
              limplim = limp
              limnlim = limn
                do 500 jn = 1, jnenlim
                if(nee /= neeb(jn)) go to 500
                if(limplim > limpsv(jn)) limplim = limpsv(jn)
                if(limnlim > limnsv(jn)) limnlim = limnsv(jn)
                  do 460 je = 1, limplim
                  if(je <= limpd) pratb(je) = pratbsv(jn, je)
                  pratt(je) = prattsv(jn, je)
                  pdratt(je) = pdratsv(jn, je)
460               continue
                  do 470 je = 1, limnlim
                  sneufe(je) = sneufsv(jn, je)
                  sneudfe(je) = sneudfsv(jn, je)
470               continue
                  jelim = maxlp
                  if(maxlp > limn + 1) jelim = limn + 1
                  if(jelim > jelimsv(jn)) jelim = jelimsv(jn)
                  do 480 je = 1, jelim
                  sneune(je) = sneunsv(jn, je)
                  sneudre(je) = sneudrsv(jn, je)
                  ineuee(je) = ineuesv(jn, je)
480               continue
if (debug) then
                write(40, 490) etaval
490             format(8x,'r2eta: reused expansion functions for eta =' &
                       ,f13.9,'.')
end if
                go to 530
500             continue
510           continue
              jnen = jnen + 1
              jnencur = jnen - (jnenmax * int((jnen - 1) / jnenmax))
              neeb(jnencur) = nee
              limbesf = 4 * int(abs(real(cc * xbninp))+ &
                      abs(aimag(cc * xbninp))) + 25
              call sphneu(cc, xbninp, limn, maxn, maxlp, ndec, nex, limbesf, &
                          sneufe, sneune, ineuee, sneudfe, sneudre)
                do je = 1, limn
                sneufsv(jnencur, je) = sneufe(je)
                sneudfsv(jnencur, je) = sneudfe(je)
                limnsv(jnencur) = limn
                end do
              jelim = maxlp
              if(maxlp > limn + 1) jelim = limn + 1
                do 520 je = 1, jelim
                sneunsv(jnencur, je) = sneune(je)
                sneudrsv(jnencur, je) = sneudre(je)
                ineuesv(jnencur, je) = ineuee(je)
520             continue
              jelimsv(jnencur) = jelim
              iopd = 3
              if(limp > maxp - 2) limp = maxp - 2
              call pleg_cached(m, limp, maxp, ndec, nex, limcsav, iopd, xlninp, &
                        netainp, maxt, prat, pdrat, pdnorma, ipdnorma, pnorma, &
                        ipnorma, alpha, beta, gamma, coefa, coefb, coefc, &
                        coefd, coefe)
              limcsav = max(limcsav, limp)
                do je = 1, limp
                pratt(je) = prat(1, je)
                pdratt(je) = pdrat(1, je)
                prattsv(jnencur, je) = pratt(je)
                pdratsv(jnencur, je) = pdratt(je)
                limpsv(jnencur) = limp
                end do
              limpd = 2 * (lnum + int(c) + ndec)
              if(limpd > maxp - 2) limpd = maxp - 2
              iopd = 2
              call pleg_cached(m, limpd, maxp, ndec, nex, limcsav, iopd, etainp, &
                        netainp, maxt, prat, pdrat, pdnorma, ipdnorma, pnorma, &
                        ipnorma, alpha, beta, gamma, coefa, coefb, coefc, &
                        coefd, coefe)
              iopd = 3
                do je = 1, limpd
                pratb(je) = prat(1, je)
                pratbsv(jnencur, je) = pratb(je)
                end do
              pratb(limpd + 1) = 0.0e0_knd
              pratb(limpd + 2) = 0.0e0_knd
              pratbsv(jnencur, limpd + 1) = 0.0e0_knd
              pratbsv(jnencur, limpd + 2) = 0.0e0_knd
530           continue
              pcoefe = ((x1 * (x1 + 2.0e0_knd)) / (x1 * (x1 + 2.0e0_knd)+ &
                      eta(nee) ** 2))
              apcoef = (rm / 2.0e0_knd) * log10(pcoefe)
              ipcoefe = int(apcoef)
              pcoefe = 10.0e0_knd ** (apcoef - ipcoefe)
              pcoefo = pcoefe * pratt(2) / pratb(2)
              ipcoefo = ipcoefe
              pdcoefe = pcoefe
              if(m /= 0) pdcoefe = -pcoefe * rm * xln(nee) * xbn(nee) * xbn(nee)/ &
                                 (x1 * (x1 + 2.0e0_knd) * wmeta2(nee))
              ipdcoefe = ipcoefe
              pdcoefo = pdcoefe * pdratt(2) / pratb(2)
              ipdcoefo = ipdcoefe
              if(li < 3) go to 540
                do jl = 3, li + ix, 2
                pcoefe = pcoefe * pratt(jl) / pratb(jl)
                iterm = log10(abs(pcoefe))
                pcoefe = pcoefe * 10.0e0_knd ** (-iterm)
                ipcoefe = ipcoefe + iterm
                pdcoefe = pdcoefe * pdratt(jl) / pratb(jl)
                iterm = log10(abs(pdcoefe))
                pdcoefe = pdcoefe * 10.0e0_knd ** (-iterm)
                ipdcoefe = ipdcoefe + iterm
                end do
              continue
              if(li < 4) go to 540
                do jl = 4, li + 1 - ix, 2
                pcoefo = pcoefo * pratt(jl) / pratb(jl)
                iterm = log10(abs(pcoefo))
                pcoefo = pcoefo * 10.0e0_knd ** (-iterm)
                ipcoefo = ipcoefo + iterm
                pdcoefo = pdcoefo * pdratt(jl) / pratb(jl)
                iterm = log10(abs(pdcoefo))
                pdcoefo = pdcoefo * 10.0e0_knd ** (-iterm)
                ipdcoefo = ipdcoefo + iterm
                end do
540           continue
              if(ix == 0) go to 550
              pcoefet = pcoefo
              ipcoefet = ipcoefo
              pcoefo = pcoefo * pratt(li + 2) / pratb(li + 2)
              iterm = int(log10(abs(pcoefo)))
              pcoefo = pcoefo * 10.0e0_knd ** (-iterm)
              ipcoefo = ipcoefo + iterm
              pdcoefet = pdcoefo
              ipdcoefet = ipdcoefo
              pdcoefo = pdcoefo * pdratt(li + 2) / pratb(li + 2)
              iterm = int(log10(abs(pdcoefo)))
              pdcoefo = pdcoefo * 10.0e0_knd ** (-iterm)
              ipdcoefo = ipdcoefo + iterm
              go to 560
550           pcoefet = pcoefe
              ipcoefet = ipcoefe
              pcoefe = pcoefe * pratt(li + 2) / pratb(li + 2)
              iterm = int(log10(abs(pcoefe)))
              pcoefe = pcoefe * 10.0e0_knd ** (-iterm)
              ipcoefe = ipcoefe + iterm
              pdcoefet = pdcoefe
              ipdcoefet = ipdcoefe
              pdcoefe = pdcoefe * pdratt(li + 2) / pratb(li + 2)
              iterm = int(log10(abs(pdcoefe)))
              pdcoefe = pdcoefe * 10.0e0_knd ** (-iterm)
              ipdcoefe = ipdcoefe + iterm
560           continue
              lplus = max(l, 1000)
              limeta = 2 * ((lplus) * (-18.5e0_knd - 20.0e0_knd * log10(x1))+ &
                     5 * ndec + 4 * m + c + 05000)
              if(x1 > 0.08e0_knd) limeta = 2 * ((lplus) * (0.50e0_knd- &
                    3.0e0_knd * log10(x1)) + 5 * ndec + 4 * m + c + 01000)
              if(x1 > 1.0e0_knd) limeta = 2 * ((lplus) * 0.5e0_knd + 5 * ndec+ &
                                         4 * m + c + 00500)
              if(iopeta == 3 .and. naccrp > minacc) &
                              limeta = jeta + jeta + 500 + c
              if(iopeta == 3 .and. naccrp <= minacc) &
                              limeta = jeta + jeta + 500 + c
              if(iopeta == 2) limeta = max(limeta, jeta + jeta + 500 + int(c))
              if(limeta > limp - 2) limeta = limp - 2
              if(limeta > limd) limeta = limd
              r1cin = r1c(li)
              ir1ein = ir1e(li)
              r1dcin = r1dc(li)
              ir1dein = ir1de(li)
              wm = wmeta2(nee)
              call r2eta(l, m, cc, x1, etaval, nee, wm, limeta, maxd, maxlp, ndec, &
                         nex, maxn, maxp, minacc, lowtest, enr, sneufe, &
                         sneune, ineuee, sneudfe, sneudre, pdratt, pratb, &
                         pratt, pcoefet, ipcoefet, pdcoefet, ipdcoefet, &
                         nsubf, nsubmf, idigc, itestm, ichoicer1, naccr1, &
                         r1cin, ir1ein, r1dcin, ir1dein, naccmax, naccrpl, &
                         naccr, kcor, r2ec, ir2ee, r2dec, ir2dee, nacceta, &
                         nacciope, jeta, iopnee, neemark, naccd, naccn, &
                         naccnmax, nacccor, naccns)
              netatry = netatry + 1
              naccetas = nacceta
if (debug) then
580           if(nacciope == 0) write(40, 590)
590           format(15x,'r2eta accuracy is calculated using the', &
                     ' Wronskian.')
              if(nacciope == 1 .and. nacccor == 0) write(40, 600)
600           format(15x,'r2eta accuracy set equal to estimated' &
                     ' numerator accuracy.')
              if(nacciope == 1 .and. nacccor /= 0) write(40, 605)
605           format(15x,'r2eta accuracy = estimated numerator', &
                     ' accuracy minus sub. error in forming' &
                     ' Wronskian.')
end if
              iopeta = 3
              if(naccetas == 0 .and. naccmax == 0 .and. iopnee == 0) &
                  neemax = nee
                if(naccetas == naccmax .and. naccetas > 0) then
                kounte = kounte + 1
                if(kounte == 1) nee1 = nee
                if(kounte == 2) nee2 = nee
                  if(kounte > 2) then
                  neemax = nee1
                  nee1 = nee2
                  nee2 = nee
                  end if
                end if
                if(naccetas > naccmax) then
                naccmax = naccetas
                kounte = 0
                neemax = nee
                end if
                if(naccetas < naccmax) then
                kounte = 0
                end if
              if(naccetas < naccr) go to 610
              naccr = min(naccetas, naccr1)
              r2c(li) = r2ec
              ir2e(li) = ir2ee
              r2dc(li) = r2dec
              ir2de(li) = ir2dee
              if(nacciope /= 0) nacciop = 1
610           continue
if (debug) then
                if(knd == kindd) then
                if(nacciope /= 0 .or. nacccor == 0 .or. naccetas == 0) &
                   write(40, 620) naccetas, etaval, nee, r2ec, ir2ee, r2dec, &
                                 ir2dee
620             format(15x,'accuracy = ',i3,' decimal digits; eta', &
                       ' = ',f17.14,'; nee = ',i4,/,10x,'r2 = ', f17.14, &
                       f17.14, i6, 2x,'r2d = ',f17.14, f17.14, i6)
                if(nacciope == 0 .and. nacccor > 0 .and. naccetas /= 0) &
                   write(40, 625) naccetas, nacccor, etaval, nee, r2ec, ir2ee, &
                                 r2dec, ir2dee
625             format(15x,'accuracy =',i3,' decimal digits,' &
                       ' adjusted for',i3,' subtraction digits in' &
                       ' forming Wronskian;'/,15x,' eta = ',f17.14,';' &
                       ' nee = ',i4,/10x,'r2 = ',f17.14, f17.14, i5, &
                       2x,'r2d = ',f17.14, f17.14, 2x, i5)
                end if
                if(knd == kindq) then
                if(nacciope /= 0 .or. nacccor == 0 .or. naccetas == 0) &
                   write(40, 630) naccetas, etaval, nee, r2ec, ir2ee, r2dec, &
                                 ir2dee
630             format(15x,'accuracy = ',i3,' decimal digits; eta', &
                       ' = ',f17.14,'; nee = ',i4,/,10x,'r2 = ', f33.30, &
                       f33.30, i6,/,10x,'r2d = ',f33.30, f33.30, i6)
                if(nacciope == 0 .and. nacccor > 0 .and. naccetas /= 0) &
                   write(40, 635) naccetas, nacccor, etaval, nee, r2ec, ir2ee, &
                                 r2dec, ir2dee
635             format(15x,'accuracy =',i3,' decimal digits,' &
                       ' adjusted for',i3,' subtraction digits in' &
                       ' forming Wronskian;'/,15x,' eta = ',f17.14,';' &
                       ' nee = ',i4,/,10x,'r2 = ',f33.30, f33.30, i6, &
                       /,10x,'r2d = ',f33.30, f33.30, i6)
                end if
end if
              if(naccetas >= minacc) ietacount = ietacount + 1
              if(ietacount >= 5) incnflag = 1
                if(naccetas >= minacc .and. iplflag == 0) then
                nee = nee - incnee
                iopeta = 2
                if(nee < 1) nee = 1
                iplflag = 1
                end if
                if(naccetas >= minacc) then
                go to 660
                end if
              iopeta = 2
              if(iplflag == 1 .and. incnflag == 1 .and. netatry == 2) &
                     iopnee = 0
              ietacount = 0
              insflag = 0
              if(naccns == ndec .and. x <= 1.11e0_knd) insflag = 1
              if((naccd >= naccdp .or. naccd >= minacc .or. naccd >= naccrp &
                   .or. naccd >= naccr) .and. (naccd + naccdp > 0) .and.  &
                  nee /= neta .and. insflag == 0) iopnee = 0
              naccdp = naccd
              if(iopnee == 0) go to 650
              if(iopnee == 2) nee = neemax - incnee
              if(iopnee == 1) nee = max(neemark, neemax) - incnee
              if(nee < 1) nee = 1
              if(iopnee == 2 .and. naccetas < lowtest - 1) go to 640
              incnee = 8
                if(knd == kindd) then
                if(x1 >= 0.01e0_knd) incnee = 16
                end if
                if(knd == kindq) then
                if(x1 >= 0.05e0_knd) incnee = 16
                if(x1 >= 0.1e0_knd) incnee = 32
                end if
640           if(naccmax >= lowtest - 1) msearch = 1
              if(msearch == 0) iopeta = 0
              if(msearch == 0) lowtest = lowtest - 1
              go to 665
650           if(nee == neta) go to 660
              nee = nee + incnee
              if(nee > neta) nee = neta
              if(msearch /= 0) kounter = 0
              go to 440
660           continue
              if(naccetas < minacc .and. nee == neta) &
                     nee = nee - incnee
              if(naccetas < minacc) iopeta = 2
              if(nee /= neta) msearch = 1
              if(naccetas >= minacc) kounter = kounter + 1
              if(kounter >= (2 * incnee) .and. msearch /= 0) &
                     incnee = 2 * incnee
              if(incnee > 64) incnee = 64
                if(knd == kindd) then
                if(x1 <= 0.2e0_knd .and. incnee > 32) incnee = 32
                if(x1 <= 0.15e0_knd .and. incnee > 16) incnee = 16
                end if
                if(knd == kindq) then
                if(x1 <= 0.1e0_knd .and. incnee > 32) incnee = 32
                end if
              if(iopint /= 0 .and. naccetas < lowacc) iopeta = 0
              if(iopeta == 0) nacctest = naccetas
              if(naccetas < minacc) iplflag = 0
665           if(naccetas < minacc .and. iflagnee == 0) iflagnee = 1
              if(nee < neest) nee = neest
670           if(naccr > 0) go to 680
              naccr = 0
              r2c = 0.0e0_knd
              ir2e = 0
              r2dc = 0.0e0_knd
              ir2de = 0
680 continue
if (output) then
              if(ioprad == 2 .and. nacciop == 0) write(20, 690) l, r1c(li), &
                     ir1e(li), r1dc(li), ir1de(li), r2c(li), ir2e(li), &
                     r2dc(li), ir2de(li), naccr, chr_w
690           format(1x, i5, 2x, 2(f17.14, 1x, f17.14, i5, 2x),/,8x, &
                     2(f17.14, 1x, f17.14, i5, 2x), i2, ' ', a)
              if(ioprad == 2 .and. nacciop /= 0) write(20, 700) l, r1c(li), &
                     ir1e(li), r1dc(li), ir1de(li), r2c(li), ir2e(li), &
                     r2dc(li), ir2de(li), naccr, chr_e
700           format(1x, i5, 2x, 2(f17.14, 1x, f17.14, i5, 2x),/,8x, &
                     2(f17.14, 1x, f17.14, i5, 2x), i2,' ', a)
end if
              if(lowacc > naccr) lowacc = naccr
              if(istartr2 == 1 .and. naccr <= minacc .and. naccrp &
                   <= minacc) then
              if(idigc - jsub > naccr .and. x1 <= 0.4e0_knd .and.  &
                 iopleg == 0 .and. l >= legstart) iopleg = 1
              if(idigc - nsubmf > naccr .and. x1 >= (0.00065e0_knd) .and.  &
                 iopneu == 0 .and. iopleg /= 2) iopneu = 1
              if(iopeta == 0 .and. x1 >= 0.00065e0_knd .and. iopneu == 0 &
                  .and. iopleg /= 2) iopeta = 1
              end if
              if(nee == neta .and. iopeta /= 0 .and. idigc - nsubmf > naccr) &
                  iopneu = 1
              if(naccr == naccneu .and. iopleg == 0 .and. l >= legstart &
                  .and. x1 < 0.1e0_knd .and. idigc - jsub > naccr) iopleg = 1
710             if(ioprad == 2 .and. naccr < 6) then
if (warn) then
                write(60, 715) naccr, m, l, x, cc
715             format(1x,' est. acc. = ',i3,' digits for m = ',i5, ' l = ',i5, &
                       ' x = ',f20.14,' c = ',f25.14, 1x, f25.14)
end if
                end if
              naccrp = naccr
              naccrplp = naccrpl
              naccr1p = naccr1
              if(ioprad == 2) nar(li) = naccr
              if(ioprad == 1) nar(li) = naccr1
720           lisave = li
              if(iopang == 0) go to 850
!
!  determine first kind prolate angular function
              if(l == m) lims1 = 3 * ndec + int(c)
              if(l /= m) lims1 = jang + jang + 20 + int(sqrt(c))
              if(lims1 > maxp) lims1 = maxp
              call s1leg(l, m, cc, iopang, iopnorm, barg, narg, lims1, ndec, &
                         maxt, maxd, maxp, enr, sgn, pr, pdr, pdnorm, ipdnorm, &
                         pnorm, ipnorm, pdtempe, ipdtempe, pdtempo, ipdtempo, &
                         ptempe, iptempe, ptempo, iptempo, s1c, is1e, s1dc, &
                         is1de, naccs, naccds, jang, dmlms, idmlmse, dmlms1, &
                         idmlms1e, jsubms)
                do 810 jarg = 1, narg
                s1(li, jarg) = s1c(jarg)
                s1d(li, jarg) = s1dc(jarg)
                is1(li, jarg) = is1e(jarg)
                is1d(li, jarg) = is1de(jarg)
                nas(li, jarg) = naccs(jarg)
                nads(li, jarg) = 0
                if(iopang == 2) nads(li, jarg) = naccds(jarg)
if (debug) then
                if(iopang == 1) write(50, 740) barg(jarg), naccs(jarg)
                if(iopang == 2) write(50, 745) barg(jarg), naccs(jarg), naccds(jarg)
740             format(1x,'eta = ',f17.14,'  accuracy = ',i3, &
                       ' digits.')
745                    format(1x,'eta = ',f17.14,'  accuracies = ',i3, &
                       ' and',i3,' digits.')
end if
if (output) then
                if(iopang == 1) write(30, 750) &
                      barg(jarg), s1c(jarg), is1e(jarg), naccs(jarg)
                if(iopang == 2) write(30, 760) &
                      barg(jarg), s1c(jarg), is1e(jarg), &
                      s1dc(jarg), is1de(jarg), naccs(jarg), naccds(jarg)
750             format(1x, f19.14, 2x, f17.14, 1x, f17.14, 2x, i5, 2x, i2)
760             format(1x, f19.14, 2x, f17.14, 1x, f17.14, 2x, i5, 2x, f17.14, 1x, f17.14, 2x, &
                        i5, 2x, i2, ', ', i2)
end if
if (debug) then
                if(knd.eq.kindd.and.iopang.eq.1) write(50, 770) s1c(jarg),is1e(jarg)
                if(knd.eq.kindd.and.iopang.eq.2) write(50, 775) s1c(jarg),is1e(jarg),s1dc(jarg),is1de(jarg)
                if(knd.eq.kindq.and.iopang.eq.1) write(50, 780) s1c(jarg),is1e(jarg)
                if(knd.eq.kindq.and.iopang.eq.2) write(50, 785) s1c(jarg),is1e(jarg),s1dc(jarg),is1de(jarg)
770             format(12x,'s1 = ',f17.14,1x,f17.14,2x,i5)
775             format(12x,'s1 = ',f17.14,1x,f17.14,2x,i5,5x,'s1d = ',f17.14,1x,f17.14,2x,i5)
780             format(12x,'s1 = ',f33.30,1x,f33.30,2x,i5)
785             format(12x,'s1 = ',f33.30,1x,f33.30,2x,i5,/12x,'s1d = ',f33.30,1x,f33.30,2x,i5)
end if
810             continue
              if(ndec - jsubms - 1 < minacc) then
if (warn) then
              write(60,*) ' est. MS norm acc. = ',ndec - jsubms - 1, &
                          ' digits for m = ',m,' l = ', l,' c = ',cc
end if
              end if
850           continue
            ijmax = min(lisave, lical)
            if(ioprad == 1) ijmax = lisave
              do i = 1, ijmax - 2
                do j = i + 2, ijmax, 2
                iegch = -int(log10(abs((eig(i) - eig(j)) / eig(i)) + dec))
                  if(iegch >= min(8, ieig(i) - 1, ieig(j) - 1)) then
if (warn) then
                  write(60,*) m, cc, i + m - 1, j + m - 1,' duplicated eigenvalue'
end if
                  end if
                end do
              end do
900         continue
          return
          end subroutine
!
!
        subroutine s1leg (l, m, c, iopang, iopnorm, barg, narg, lims1, ndec, &
                          maxt, maxd, maxp, enr, sgn, pr, pdr, pdnorm, ipdnorm, &
                          pnorm, ipnorm, pdtempe, ipdtempe, pdtempo, &
                          ipdtempo, ptempe, iptempe, ptempo, iptempo, s1c, &
                          is1e, s1dc, is1de, naccs, naccds, jang, dmlms, &
                          idmlmse, dmlms1, idmlms1e, jsubms)
!
!  purpose:     to calculate the prolate angular functions of the first
!               kind and their first derivatives with respect to eta.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               iopang  : index = 1 when angular functions of the
!                         first kind are calculated; = 2 when the
!                         first derivatives with respect to eta are
!                         also calculated
!               iopnorm : = 1 when the angular functions (and
!                         first derivatives) are scaled by the
!                         square root of the normalization of the
!                         corresponding Legendre function, giving them
!                         unity norm; iopnorm = 0 otherwise
!               barg    : array of eta values for which angular
!                         functions are desired
!               narg    : number of eta values
!               lims1   : approximately twice the maximum number
!                         of terms available to be taken in the
!                         sums
!               ndec    : number of decimal digits for real(knd)
!               maxt    : dimension of barg, pdnorm, ipdnorm, pnorm,
!                         ipnorm, pdtempe, ipdtempe, pdtempo, ipdtempo,
!                         ptempe, iptempe, ptempo, iptempo, s1c, is1e,
!                         s1dc, is1de, and naccs arrays. first dimension
!                         of the doubly dimensioned arrays pr and pdr
!               maxd    : dimension of enr array
!               maxp    : second dimension of pr and pdr arrays
!               enr     : array of d coefficient ratios
!               sgn     : sign of d coefficient multiplying the
!                         associated Legendre of order m and
!                         degree l in the series for s1 and s1d
!               pr      : array of ratios of successive first kind
!                         associated Legendre functions of the same
!                         parity
!               pdr     : array of ratios of successive derivatives of
!                         first kind associated Legendre functions of
!                         the same parity
!               pdnorm  : array of characteristics of the first
!                         derivatives of associated Legendre functions
!                         of the first kind of order m and degree m
!               ipdnorm : array of exponents corresponding to pdnorm
!               pnorm   : array of characteristics of the associated
!                         Legendre functions of the first kind of order
!                         m and degree m
!               ipnorm  : array of exponents corresponding to pnorm
!               pdtempe : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is even
!                         or odd, to the first derivative of the
!                         function of order m and degree m
!               ipdtempe: array of exponents corresponding to pdtempe
!               pdtempo : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is odd
!                         or even, to the first derivtive of the
!                         function of order m and degree m
!               idptempo: array of exponents corresponding to pdtempo
!               ptempe  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is even or odd, to the function
!                         of order m and degree m
!               iptempe : array of exponents corresponding to ptempe
!               ptempo  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is odd or even, to the
!                         function of order m and degree m
!               iptempo : array of exponents corresponding to ptempo
!
!
!     output:   s1c    : array of characteristics of prolate
!                        angular functions of the first kind
!               is1e   : array of exponents of prolate angular
!                        functions of the first kind
!               s1dc   : array of characteristics of derivative with
!                        respect to eta of prolate angular functions
!                        of the first kind
!               is1de  : array of exponents of derivatives with respect
!                        to eta of prolate angular functions of first
!                        kind
!               naccs  : array of integer estimates of the number of
!                        accurate decimal digits in the values obtained
!                        for s1
!               naccds : array of integer estimates of the number of
!                        accurate decimal digits in the values obtained
!                        for s1d
!               jang   : maximum value of the index j in the forward
!                        sum for r1 and r1d, i.e., the highest enr(j)
!                        used
!               dmlms  : characteristic of the d coefficient with index
!                        l - m when using the Meixner and Schafke
!                        normalization for the angular functions
!               idmlmse: exponent associated with dmlms
!               dmlms1 : characteristic of the d coefficient with index
!                        l - m when using unit normalization for the
!                        angular functions
!               idmlms1e: exponent associated with dmlms1
!               jsubms  : effective number of decimal digits of
!                         subtraction error incurred in calculating
!                         dmsnorm
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) adec, aj, aj2, dcon, dec, factor, fterm, rm2, rm2m1, rm2m3, &
                  rm2p1, sgn
        complex(knd) c, coef, dnew, dnewd, dmsnorm, dmlms, dmlms1, dold, doldd, &
                     s1, s1d
        real(knd) barg(maxt), pdr(maxt, maxp), pdnorm(maxt), &
                  pnorm(maxt), pr(maxt, maxp), pdtemp(maxt), ptemp(maxt), &
                  pdtempe(maxt), ptempe(maxt), pdtempo(maxt), ptempo(maxt)
        complex(knd) enr(maxd), s1c(maxt), s1dc(maxt)
!
!  integer arrays
        integer ipdnorm(maxt), ipnorm(maxt), ipdtemp(maxt), iptemp(maxt), &
                ipdtempe(maxt), iptempe(maxt), ipdtempo(maxt), &
                iptempo(maxt), is1de(maxt), is1e(maxt), naccs(maxt), &
                naccds(maxt)
!
        dec = 10.0e0_knd ** (-ndec - 1)
        dcon = dec
        adec = 1000.0e0_knd * dec
        rm2 = real(m + m, knd)
        rm2m1 = real(m + m - 1, knd)
        rm2p1 = real(m + m + 1, knd)
        rm2m3 = real(m + m - 3, knd)
        if(l > (m + 1)) go to 30
          do 20 k = 1, narg
          if(pnorm(k) == 0.0e0_knd) go to 20
          if(l == (m + 1)) go to 10
          ptempe(k) = pr(k, 1)
          iptempe(k) = 0
          pdtempe(k) = pdr(k, 1)
          ipdtempe(k) = 0
          ptemp(k) = ptempe(k)
          pdtemp(k) = pdtempe(k)
          iptemp(k) = 0
          ipdtemp(k) = 0
          go to 20
10        ptempo(k) = pr(k, 2)
          iptempo(k) = 0
          pdtempo(k) = pdr(k, 2)
          ipdtempo(k) = 0
          ptemp(k) = ptempo(k)
          pdtemp(k) = pdtempo(k)
          iptemp(k) = 0
          ipdtemp(k) = 0
20        continue
30      continue
        lm2 = (l - m) / 2
        ix = l - m - 2 * lm2
        ixx = ix - 1
        ixx2 = ixx + 2
        if(l < (m + 2)) go to 110
          do 100 k = 1, narg
          if(pnorm(k) == 0.0e0_knd) go to 100
          if(ix /= 0) go to 60
          ptempe(k) = ptempe(k) * pr(k, l - m + 1)
          if(abs(ptempe(k)) < 1.0e+10_knd) go to 40
          ptempe(k) = ptempe(k) * (1.0e-10_knd)
          iptempe(k) = iptempe(k) + 10
40        ptemp(k) = ptempe(k)
          iptemp(k) = iptempe(k)
          if(abs(barg(k)) < adec) go to 100
          pdtempe(k) = pdtempe(k) * pdr(k, l - m + 1)
          if(abs(pdtempe(k)) < 1.0e+10_knd) go to 50
          pdtempe(k) = pdtempe(k) * (1.0e-10_knd)
          ipdtempe(k) = ipdtempe(k) + 10
50        pdtemp(k) = pdtempe(k)
          ipdtemp(k) = ipdtempe(k)
          go to 100
60        if(abs(barg(k)) < adec) go to 80
          ptempo(k) = ptempo(k) * pr(k, l - m + 1)
          if(abs(ptempo(k)) < 1.0e+10_knd) go to 70
          ptempo(k) = ptempo(k) * (1.0e-10_knd)
          iptempo(k) = iptempo(k) + 10
70        ptemp(k) = ptempo(k)
          iptemp(k) = iptempo(k)
80        pdtempo(k) = pdtempo(k) * pdr(k, l - m + 1)
          if(abs(pdtempo(k)) < 1.0e+10_knd) go to 90
          pdtempo(k) = pdtempo(k) * (1.0e-10_knd)
          ipdtempo(k) = ipdtempo(k) + 10
90        pdtemp(k) = pdtempo(k)
          ipdtemp(k) = ipdtempo(k)
100        continue
110     continue
        lim = lims1 / 2 - ix
!
!  compute the normalizing factor
        dmsnorm = (1.0e0_knd, 0.0e0_knd)
        fterm = 1.0e0_knd
        idmse = 0
        coef = (1.0e0_knd, 0.0e0_knd)
        jlow = l - m + 2
        jterm = lm2
          do 120 j = jlow, lims1, 2
          aj = real(j, knd)
          aj2 = aj + aj
          jterm = jterm + 1
          coef = coef * (aj + rm2) * enr(jterm) * (aj + rm2m1) * enr(jterm) &
               *(aj2 + rm2m3) / (aj * (aj - 1.0e0_knd) * (aj2 + rm2p1))
          dmsnorm = dmsnorm + coef
          if(abs(dmsnorm) > fterm) fterm = abs(dmsnorm)
          if(abs(coef / dmsnorm) < dcon) go to 130
120       continue
130     jlow = l - m
        jn = jterm
        if(jlow < 2) go to 150
        coef = (1.0e0_knd, 0.0e0_knd)
        jterm = lm2
        j = jlow
          do 140 jj = 2, jlow, 2
          aj = real(j, knd)
          aj2 = aj + aj
          coef = coef * aj * (aj - 1.0e0_knd) * (aj2 + rm2p1) / ((aj2 + rm2m3)* &
               enr(jterm) * enr(jterm) * (aj + rm2) * (aj + rm2m1))
          jterm = jterm - 1
          j = j - 2
          dmsnorm = dmsnorm + coef
          if(abs(dmsnorm) > fterm) fterm = abs(dmsnorm)
          if(abs(coef / dmsnorm) < dcon) go to 150
140       continue
150     jsubms = int(log10((fterm / abs(dmsnorm)) + dec))
        if(jsubms < 0) jsubms = 0
        if(jsubms > ndec) jsubms = ndec
        iterm = int(log10(abs(dmsnorm)))
        dmsnorm = dmsnorm * (10.0e0_knd ** (-iterm))
        idmse = idmse + iterm
          if(2 * (idmse / 2) /= idmse) then
          idmse = idmse - 1
          dmsnorm = 10.0_knd * dmsnorm
          end if
        dmlms = sgn / sqrt(dmsnorm)
        idmlmse = -idmse / 2
if (debug) then
        write(50, 155) jn, lim
155     format(5x,'Meixner-Schafke normalization series converged in ' &
                 ,i6,' terms; ',i6,' terms available.')
end if
        jlow = lm2 + 1
        jang = 0
!
!  compute the associated Legendre function normalization factor
        factor = 1.0e0_knd
        ifactor = 0
        if(iopnorm == 0) go to 210
        if(m == 0) go to 170
          do 160 j = 1, m
          aj = real(j, knd)
          factor = factor * (aj + aj) * (aj + aj - 1.0e0_knd)
          if(factor < 1.0e100_knd) go to 160
          factor = factor * 1.0e-100_knd
          ifactor = ifactor + 100
160       continue
170     if(l == m) go to 190
          do 180 j = 1, l - m
          aj = real(j, knd)
          factor = factor * (rm2 + aj) / (aj)
          if(factor < 1.0e100_knd) go to 180
          factor = factor * 1.0e-100_knd
          ifactor = ifactor + 100
180       continue
190     factor = factor * 2.0e0_knd / (l + l + 1.0e0_knd)
        factor = sqrt(factor)
        ifactor = ifactor / 2
        iterm = int(log10(factor))
        factor = factor * (10.0e0_knd ** (-iterm))
        ifactor = ifactor + iterm
        dmlms1 = dmlms / factor
        idmlms1e = idmlmse - ifactor
if (debug) then
        write(50, 200)
200     format(5x,'s1 is normalized to have unit norm.')
end if
210     continue
if (debug) then
        if(iopnorm == 0) write(50, 215)
215     format(5x,'s1 has the same normalization as the', &
               ' corresponding Legendre function.')
end if
!
!  compute the angular function s1
          do 380 k = 1, narg
          if(pnorm(k) == 0.0e0_knd) go to 220
          if((ix == 1) .and. (abs(barg(k)) < adec)) go to 220
          if(((abs(abs(barg(k)) - 1.0e0_knd)) < adec) &
               .and.(m /= 0)) go to 220
          go to 230
220       s1c(k) = (0.0e0_knd, 0.0e0_knd)
          is1e(k) = 0
          naccs(k) = ndec
          go to 300
230       dold = (1.0e0_knd, 0.0e0_knd)
          s1 = dold
          fterm = 1.0e0_knd
            do 240 j = jlow, lim
            dnew = dold * enr(j) * pr(k, j + j + ixx2)
            s1 = s1 + dnew
            if(abs(dnew) > fterm) fterm = abs(dnew)
            if(abs(dnew / s1) < dcon) go to 250
            dold = dnew
240         continue
250       if(j > jang) jang = j
if (debug) then
          write(50, 260) barg(k), j
260       format(8x,'s1 calculation for eta = ',f17.14,' converged in ', &
                 i6,' terms')
end if
          if(lm2 < 1) go to 280
          dold = (1.0e0_knd, 0.0e0_knd)
          j = lm2
            do 270 jj = 1, lm2
            dnew = dold / (pr(k, j + j + ixx2) * enr(j))
            s1 = s1 + dnew
            if(abs(dnew) > fterm) fterm = abs(dnew)
            if(abs(dnew / s1) < dcon) go to 280
            dold = dnew
            j = j - 1
270         continue
280       s1c(k) = s1 * dmlms * (ptemp(k) * pnorm(k) / factor)
          if(s1c(k) /= (0.0e0_knd, 0.0e0_knd)) &
                        iterm = int(log10(abs(s1c(k))))
          if(s1c(k) == (0.0e0_knd, 0.0e0_knd)) iterm = 0
          s1c(k) = s1c(k) * (10.0e0_knd ** (-iterm))
          is1e(k) = iptemp(k) + ipnorm(k) + iterm + idmlmse - ifactor
            if(abs(real(s1c(k))) > 10.0e0_knd .or. abs(aimag(s1c(k))) &
                > 10.0e0_knd) then
            s1c(k) = s1c(k) / 10.0e0_knd
            is1e(k) = is1e(k) + 1
            end if
            if(abs(real(s1c(k))) < 1.0e0_knd .and. abs(aimag(s1c(k))) &
                < 1.0e0_knd) then
            s1c(k) = s1c(k) * 10.0e0_knd
            is1e(k) = is1e(k) - 1
            end if
          if(s1 == (0.0e0_knd, 0.0e0_knd)) naccs(k) = 0
          if(s1 /= (0.0e0_knd, 0.0e0_knd)) naccs(k) = ndec - 2 - &
              log10(abs(fterm / s1))
          if(naccs(k) < 0) naccs(k) = 0
          if(naccs(k) > 0) go to 300
          naccs(k) = 0
          s1c(k) = (0.0e0_knd, 0.0e0_knd)
          is1e(k) = 0
          naccds(k) = 0
          s1dc(k) = (0.0e0_knd, 0.0e0_knd)
          is1de(k) = 0
          go to 380
!
!       compute the first derivative of the angular function when
!       iopang equals 2
300       if(iopang /= 2) go to 380
          if(pnorm(k) == 0.0e0_knd) go to 310
          if((ix == 0) .and. (abs(barg(k)) < adec)) go to 310
          if(((abs(abs(barg(k)) - 1.0e0_knd)) < adec) .and. (m /= 0) &
               .and. (m /= 2)) go to 310
          go to 320
310       s1dc(k) = (0.0e0_knd, 0.0e0_knd)
          is1de(k) = 0
          naccds(k) = ndec
          go to 380
320       doldd = (1.0e0_knd, 0.0e0_knd)
          s1d = doldd
          fterm = 1.0e0_knd
          if(l == 0) s1d = (0.0e0_knd, 0.0e0_knd)
            do 330 j = jlow, lim
            dnewd = doldd * enr(j) * pdr(k, j + j + ixx2)
            s1d = s1d + dnewd
            if(abs(dnewd) > fterm) fterm = abs(dnewd)
            if(abs(dnewd / s1d) < dcon) go to 340
            doldd = dnewd
330         continue
340       if(lm2 < 1) go to 360
          doldd = (1.0e0_knd, 0.0e0_knd)
          j = lm2
          ja = lm2
          if(m == 0 .and. ix == 0) ja = lm2 - 1
          if(ja == 0) go to 360
            do 350 jj = 1, ja
            dnewd = doldd / (pdr(k, j + j + ixx2) * enr(j))
            s1d = s1d + dnewd
            if(abs(dnewd) > fterm) fterm = abs(dnewd)
            if(abs(dnewd / s1d) < dcon) go to 360
            doldd = dnewd
            j = j - 1
350         continue
360       s1dc(k) = s1d * dmlms * (pdtemp(k) * pdnorm(k) / factor)
          if(s1d == (0.0e0_knd, 0.0e0_knd)) naccds(k) = 0
          if(s1d /= (0.0e0_knd, 0.0e0_knd)) naccds(k) = ndec - 2 - &
              log10(abs(fterm / s1d))
          if(naccds(k) < 0) naccds(k) = 0
          if(s1dc(k) /= (0.0e0_knd, 0.0e0_knd)) &
               iterm = int(log10(abs(s1dc(k))))
          if(s1dc(k) == (0.0e0_knd, 0.0e0_knd)) iterm = 0
          s1dc(k) = s1dc(k) * 10.0e0_knd ** (-iterm)
          is1de(k) = ipdtemp(k) + ipdnorm(k) + iterm + idmlmse - ifactor
            if(abs(real(s1dc(k))) > 10.0e0_knd .or. abs(aimag(s1dc(k))) &
                > 10.0e0_knd) then
            s1dc(k) = s1dc(k) / 10.0e0_knd
            is1de(k) = is1de(k) + 1
            end if
            if(abs(real(s1dc(k))) < 1.0e0_knd .and. abs(aimag(s1dc(k))) &
                < 1.0e0_knd) then
            s1dc(k) = s1dc(k) * 10.0e0_knd
            is1de(k) = is1de(k) - 1
            end if
380       continue
        return
        end subroutine
!
!
        subroutine r1besa (l, m, c, x1, limr1, maxd, enr, maxj, maxlp, ndec, nex, &
                           iflag, sbesf, sbesdf, sbesn, ibese, sbesdr, &
                           d01, id01, dfnorm, idfe, r1c, ir1e, r1dc, ir1de, &
                           jbesa, factor, nsubr1a)
!
!  purpose:     to calculate the prolate radial function of the
!               first kind and its first derivative with respect
!               to x, using an expansion of spherical Bessel
!               functions of the first kind with argument
!               c*sqrt(x*x-1).
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : complex c
!               x1     : x-1
!               limr1  : approximately twice the maximum number of
!                        terms available to be taken in the series
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               maxj   : dimension of sbesf and sbesdf arrays
!               maxlp  : maximum  l value desired; dimension
!                        of the sbesdr, sbesn, and ibese arrays
!               ndec   : number of decimal digits for kind = knd
!               nex    : maximum exponent for kind = knd
!               iflag  : integer = 1 if forward series not needed;
!                        =0 if the forward series is computed
!               sbesf  : array of ratios of successive first kind
!                        spherical Bessel functions of the same parity
!               sbesdf : array of ratios of successive derivatives of
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : array of characteristics for Bessel functions
!               ibese  : array of exponents corresponding to sbesn
!               sbesdr : value of ratio of first derivative of
!                        spherical Bessel function to the corresponding
!                        Bessel function
!               d01    : characteristic of the ratio of the first d
!                        coefficient, either d0 or d1, to the d
!                        coefficient for n = l - m
!               id01   : exponent for d01
!               dfnorm : characteristic of the Flammer normalization sum
!                        for the d coefficients
!               idfe   : exponent associated with dfnorm
!
!     output  : r1c    : characteristic of prolate radial function
!                        of the first kind
!               ir1e   : exponent of prolate radial function of the
!                        first kind
!               r1dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the first
!                        kind
!               ir1de  : exponent of derivative with respect to x of
!                        prolate radial function of the first kind
!               jbesa  : maximum value of the index j in the forward
!                        sum for r1 and r1d, i.e., the highest enr(j)
!                        used
!               factor : coefficient used in computing one of the two
!                        contributions to r1d for l odd and m equal to
!                        0. The value used for the previous odd l is
!                        input to the subroutine. The value for the
!                        current odd value of l is computed and output
!                        to be availabe for the next odd value of l.
!               nsubr1a: maximum number of digits of subtraction error
!                        in calculating r1 and r1d using this subroutine
!
        use param
!
!  real(knd) and complex*16 scalars and arrays
        real(knd) con, dec, factor, r1rp, r1ip, r1drp, r1dip, rj, teste, testeo, &
                  x1, fdterm, fterm
        complex(knd) c, dfnorm, dnew, dnewd, dold, doldd, d01, r1bot, r1c, r1d, &
                     r1dc, r1dstore, r1temp, r1top, r1dtop, r1topd, &
                     r1dtopd, term, termd
        complex(knd) enr(maxd), sbesdf(maxj), sbesdr(maxlp), sbesf(maxj), &
                     sbesn(maxlp)
!
!  integer array
        integer ibese(maxlp)
!
!  convergence ratio dec is set according to the requested accuracy
        dec = 10.0e0_knd ** (-ndec - 1)
        lm2 = (l - m) / 2
!
!  ix=0 for l-m even, ix=1 for l-m odd
        ix = l - m - 2 * lm2
        lim = limr1 / 2 - ix
        con = (x1 + 1.0e0_knd) / (sqrt(x1) * sqrt((x1 + 2.0e0_knd)))
        nfac = nex / 3
        teste = 10.0e0_knd ** nfac
        testeo = 1.0e0_knd / teste
        ir1tope = 0
        mml = m + m - 1 + ix
        iflagd = 0
        if(x1 < 0.1e0_knd .and. ix == 1 .and. m == 0) iflagd = 1
        if(iflagd == 1 .and. l /= 1) factor = factor * real(l, knd)/ &
            real(l - 1, knd)
!
!
!  compute radial function of the first kind r1 and its first
!  derivative r1d
!
!  forward summation of numerator series
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
        r1rp = 1.0e0_knd
        r1drp = 1.0e0_knd
        r1ip = 0.0e0_knd
        r1dip = 0.0e0_knd
        r1top = dold
        r1dtop = doldd
        jbesa = lm2
        fterm = 1.0e0_knd
        fdterm = 1.0e0_knd
        if(iflag == 1) go to 50
          do 20 j = lm2 + 1, lim
          jj = j + j + ix
          rj = real(jj)
          dnew = dold * enr(j) * sbesf(jj + m) * real((jj + mml), knd)/ &
                real((jj - ix), knd)
          dnewd = doldd * enr(j) * sbesdf(jj + m) * real((jj + mml), knd)/ &
                real((jj - ix), knd)
          if(abs(dnew) > fterm) fterm = abs(dnew)
          if(abs(dnewd) > fdterm) fdterm = abs(dnewd)
          r1top = r1top + dnew
          r1dtop = r1dtop + dnewd
          if(real(dnew) > 0.0e0_knd) r1rp = r1rp + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r1drp = r1drp + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r1ip = r1ip + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r1dip = r1dip + aimag(dnewd)
          if((abs(dnew / r1top) + abs(dnewd / r1dtop)) < dec) go to 30
          dold = dnew
          doldd = dnewd
20        continue
30      continue
        jbesa = min(j, lim)
        if(iflagd == 0 .or. l /= 1) go to 40
        r1topd = r1top - 1.0e0_knd
        r1dtopd = r1dtop - 1.0e0_knd
40      continue
50  continue
!
!  backward summation of numerator series
        if (lm2 < 1) go to 80
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
          do 70 j = lm2, 1,-1
          jj = j + j + ix
          rj = real(jj, knd)
          dnew = dold * (jj - ix) / (sbesf(jj + m) * real((jj + mml), knd) * enr(j))
          dnewd = doldd * (jj - ix) / (sbesdf(jj + m) * real((jj + mml), knd) * enr(j))
          if(abs(dnew) > fterm) fterm = abs(dnew)
          if(abs(dnewd) > fdterm) fdterm = abs(dnewd)
          r1top = r1top + dnew
          r1dtop = r1dtop + dnewd
          if(real(dnew) > 0.0e0_knd) r1rp = r1rp + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r1drp = r1drp + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r1ip = r1ip + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r1dip = r1dip + aimag(dnewd)
          if((abs(dnew / r1top) + abs(dnewd / r1dtop)) < dec) go to 80
          if(abs(r1top) < teste) go to 60
          r1top = r1top * testeo
          r1rp = max(r1rp * testeo, dec * real(r1top))
          r1ip = max(r1ip * testeo, dec * aimag(r1top))
          dnew = dnew * testeo
          ir1tope = ir1tope + nfac
          r1dtop = r1dtop * testeo
          r1drp = max(r1drp * testeo, dec * real(r1dtop))
          r1dip = max(r1dip * testeo, dec * aimag(r1dtop))
          dnewd = dnewd * testeo
60        dold = dnew
          doldd = dnewd
70        continue
          if(jj /= 3) iflagd = 0
          if(iflagd == 0) go to 80
          r1topd = r1top - dnew
          r1dtopd = r1dtop - dnewd
          go to 90
80        continue
          iflagd = 0
90        nsubr = 0
          if(r1rp /= 0.0e0_knd) nsubr = -log10(abs(real(r1top) / r1rp) + dec)
          if(nsubr < 0) nsubr = 0
          nsubi = 0
          if(r1ip /= 0.0e0_knd) nsubi=- &
                   log10(abs(aimag(r1top) / r1ip) + dec)
          if(nsubi < 0) nsubi = 0
          nsub = nsubr
          if(aimag(r1top) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(r1top) / real(r1top))))
          if(iterm < 0) nsub = max(nsubr, nsubi + iterm)
          if(iterm > 0) nsub = max(nsubi, nsubr - iterm)
          end if
          nsub1 = -int(log10(abs(r1top) / fterm + dec))
          nsubd1 = -int(log10(abs(r1dtop) / fdterm + dec))
          nsubdr = 0
          if(r1drp /= 0.0e0_knd) &
                nsubdr = -log10(abs(real(r1dtop) / r1drp) + dec)
          if(nsubdr < 0) nsubdr = 0
          nsubdi = 0
          if(aimag(r1dtop) /= 0.0e0_knd) nsubdi=- &
                              log10(abs(aimag(r1dtop) / r1dip) + dec)
          if(nsubdi < 0) nsubdi = 0
          nsubd = nsubdr
          if(aimag(r1dtop) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(r1dtop) / real(r1dtop))))
          if(iterm < 0) nsubd = max(nsubdr, nsubdi + iterm)
          if(iterm > 0) nsubd = max(nsubdi, nsubdr - iterm)
          end if
          r1bot = dfnorm
!
!  compute r1 and r1d
        r1temp = r1top * sbesn(l + 1) / r1bot
        if(ix == 1) r1temp = r1temp * con
        iterm = log10(abs(r1temp))
        ir1e = ir1tope + ibese(l + 1) + iterm - idfe
        r1c = r1temp * (10.0e0_knd ** (-iterm))
          if(abs(real(r1c)) > 10.0e0_knd .or. abs(aimag(r1c)) >  &
              10.0e0_knd) then
          r1c = r1c / 10.0e0_knd
          ir1e = ir1e + 1
          end if
          if(abs(real(r1c)) < 1.0e0_knd .and. abs(aimag(r1c)) <  &
              1.0e0_knd) then
          r1c = r1c * 10.0e0_knd
          ir1e = ir1e - 1
          end if
        if(iflagd == 1) r1temp = r1temp * r1topd / r1top
        if(iflagd == 1) r1dtop = r1dtopd
        r1d = r1dtop * sbesn(l + 1) * c * con * sbesdr(l + 1) / r1bot
        r1dstore = r1d * con
        term = r1temp / (x1 * (x1 + 1.0e0_knd) * (x1 + 2.0e0_knd))
          if(ix == 1) then
          r1d = r1dstore - term
          iterm = int(log10(abs(r1dstore / term) + dec))
          if(iterm > 0) nsubd = max(nsubd, nsub - iterm)
          if(iterm <= 0) nsubd = max(nsub, nsubd + iterm)
          iterm = int(log10(abs(r1dstore / r1d)))
          if(iterm < 0) iterm = 0
          nsubd = nsubd + iterm
          end if
        if(iflagd == 0) go to 110
        term = x1 * (x1 + 2.0e0_knd) * sbesdr(2) * sbesn(2)
        termd = term - sbesn(3) * (10.0e0_knd ** (ibese(3) - ibese(2)))
        termd = termd * (c * d01 / (x1 * (x1 + 2.0e0_knd) * r1bot * factor))* &
             (10.0e0_knd ** (id01 + ibese(2) - ibese(l + 1) - ir1tope))
        r1d = r1d + termd
          if(abs(termd) /= 0.0e0_knd) then
          iterm = int(log10(abs(termd / r1d)))
          if(iterm < 0) iterm = 0
          nsubd = nsubd + iterm
          end if
110     continue
        nsubr1a = max(nsub, nsubd)
        if(nsubr1a > ndec) nsubr1a = ndec
if (debug) then
        if(iflag == 0) write(40, 120) jbesa, lim, nsub, nsubd
120     format(8x,'r1besa: numerator conv. in ',i6,' terms with',i6, &
               ' terms available;',/,15x, i3,' and',i3,' digits of subt.' &
               ' error in r1 and r1d.')
        if(iflag == 1) write(40, 130) jbesa, lim, nsub, nsubd
130     format(8x,'r1besa: numerator conv. in ',i6,' terms with',i6, &
               ' terms available;',i3,' digits of'/,15x,' subt. error ' &
               'in r1,',i3,' digits in r1d. forward sum not required.')
end if
        iterm = log10(abs(r1d))
        ir1de = ir1tope + ibese(l + 1) + iterm - idfe
        r1dc = r1d * (10.0e0_knd ** (-iterm))
          if(abs(real(r1dc)) > 10.0e0_knd .or. abs(aimag(r1dc)) >  &
              10.0e0_knd) then
          r1dc = r1dc / 10.0e0_knd
          ir1de = ir1de + 1
          end if
          if(abs(real(r1dc)) < 1.0e0_knd .and. abs(aimag(r1dc)) <  &
              1.0e0_knd) then
          r1dc = r1dc * 10.0e0_knd
          ir1de = ir1de - 1
          end if
        mfac = ir1e - ibese(l + 1)
        if(mfac > (ndec + 5)) iflag = 1
        if(mfac <= (ndec + 5)) iflag = 0
        return
        end subroutine
!
!
        subroutine r1besb (l, m, c, x1, limbes, maxd, maxlp, ndec, maxj, &
                           maxp, enr, sbesf, sbesn, ibese, sbesdf, sbesdr, &
                           prat1, pcoefn, ipcoefn, dmfnorm, idmfe, r1c, ir1e, &
                           r1dc, ir1de, jbesb, nsubr1b)
!
!  purpose:     to calculate the prolate radial function of the
!               first kind and its first derivative with respect
!               to x, using the traditional expansions in terms of
!               spherical Bessel functions.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x1     : x-1
!               limbes : maximum number of terms to be taken in the
!                        series summations for r1 and r1d
!               maxd   : dimension of enr array
!               maxlp  : maximum  l value desired; dimension
!                        of the sbesn and ineue arrays
!               ndec   : number of decimal digits for kind = knd
!               maxj   : dimension of sbesf and sbesdr arrays
!               maxp   : dimension of prat1 array
!               enr    : array of ratios of successive d coefficients
!               sbesf  : array of ratios of successive spherical Bessel
!                        functions of the same parity
!               sbesn  : array of characteristics for Bessel functions
!               ibese  : array of exponents for Bessel functions
!               sbesdf : array of ratios of successive first derivatives
!                        of spherical Bessel functions of same parity
!               sbesdr : array of ratios of first derivatives of Bessel
!                        functions to the corresponding Bessel functions
!               prat1  : array of ratios of successive coefficients in
!                        r1 and r1d sum
!               pcoefn : characteristic of coefficient for term in both
!                        r1 and r1d sums that contains the Bessel function
!                        of order l + m
!               ipcoefn: exponent (to the base 10) corresponding to
!                        pcoefn
!               dmfnorm: characteristic of the Morse and Feshbach
!                        normalization factor of the d coefficients.
!                        equal to the reciprocal of the characteristic
!                        of the d coefficient d(n = l - m) using this
!                        normalization for the angular functions
!               idmfe  : expontent associated with dmfnorm
!
!     output:   r1c    : characteristic of prolate radial function
!                        of the first kind
!               ir1e   : exponent of prolate radial function of the
!                        first kind
!               r1dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the first
!                        kind
!               ir1de  : exponent of derivative with respect to x of
!                        prolate radial function of the first kind
!               jbesb  : index of term where convergence is
!                        achieved for r1 or for r1d, whichever index is
!                        larger
!               nsubr1b: maximum number of digits of subtraction error
!                        in calculating r1 and r1d using this subroutine
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) dcon, dec, dnewi, dnewr, dnewdi, dnewdr, pcoefn, rm, rm2, &
                  r1dcoef, r1dip, r1ip, r1drp, r1rp, x1
        complex(knd) c, dmfnorm, dnew, dnewd, dold, doldd, r1c, r1dc, &
                     r1dc1, r1dc2, r1dtemp, r1temp
        real(knd) prat1(maxp)
        complex(knd) enr(maxd), sbesdr(maxlp), sbesn(maxlp), &
                     sbesf(maxj), sbesdf(maxj)
!
!  integer arrays
        integer ibese(maxlp)
!
        rm = real(m, knd)
        dec = 10.0e0_knd ** (-ndec)
        dcon = dec
        r1dcoef = rm / ((x1 + 1.0e0_knd) * (x1 + 2.0e0_knd) * x1)
        rm2 = rm * 2.0e0_knd
        lm2 = (l - m) / 2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix = l - m - 2 * lm2
        lim = limbes / 2 - ix
!
!  compute radial function of the first kind
!
!  backward series
        r1temp = (1.0e0_knd, 0.0e0_knd)
        r1dtemp = (1.0e0_knd, 0.0e0_knd)
        r1rp = 1.0e0_knd
        r1ip = 0.0e0_knd
        r1drp = 1.0e0_knd
        r1dip = 0.0e0_knd
        if (lm2 < 1) go to 20
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
          do 10 j = lm2, 1,-1
          jj = j + j + ix
          dnew = -dold / (sbesf(jj + m) * prat1(jj + 1) * enr(j))
          dnewd = -doldd / (sbesdf(jj + m) * prat1(jj + 1) * enr(j))
          r1temp = r1temp + dnew
          r1dtemp = r1dtemp + dnewd
          dnewr = real(dnew)
          if(dnewr > 0.0e0_knd) r1rp = r1rp + dnewr
          dnewi = aimag(dnew)
          if(dnewi > 0.0e0_knd) r1ip = r1ip + dnewi
          dnewdr = real(dnewd)
          if(dnewdr > 0.0e0_knd) r1drp = r1drp + dnewdr
          dnewdi = aimag(dnewd)
          if(dnewdi > 0.0e0_knd) r1dip = r1dip + dnewdi
          dold = dnew
          doldd = dnewd
10      continue
20      continue
!
!  forward series
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
          do 30 j = lm2 + 1, lim
          jj = j + j + ix
          dnew = -dold * enr(j) * sbesf(jj + m) * prat1(jj + 1)
          dnewd = -doldd * enr(j) * sbesdf(jj + m) * prat1(jj + 1)
          r1temp = r1temp + dnew
          r1dtemp = r1dtemp + dnewd
          dnewr = real(dnew)
          if(dnewr > 0.0e0_knd) r1rp = r1rp + dnewr
          dnewi = aimag(dnew)
          if(dnewi > 0.0e0_knd) r1ip = r1ip + dnewi
          dnewdr = real(dnewd)
          if(dnewdr > 0.0e0_knd) r1drp = r1drp + dnewdr
          dnewdi = aimag(dnewd)
          if(dnewdi > 0.0e0_knd) r1dip = r1dip + dnewdi
          if((abs(dnew / r1temp) + abs(dnewd / r1dtemp)) < dcon) go to 40
          dold = dnew
          doldd = dnewd
30        continue
40      continue
        jbesb = min(j, lim)
!
!  combining results to form the radial function characteristics
!  r1c and r1dc and corresponding exponents ir1e and ir1de
        r1c = r1temp * sbesn(l + 1) * pcoefn / dmfnorm
        iterm = int(log10(abs(r1c)))
        ir1e = ibese(l + 1) + ipcoefn + iterm - idmfe
        r1c = r1c * 10.0e0_knd ** (-iterm)
          if(abs(real(r1c)) > 10.0e0_knd .or. abs(aimag(r1c)) >  &
              10.0e0_knd) then
          r1c = r1c / 10.0e0_knd
          ir1e = ir1e + 1
          end if
          if(abs(real(r1c)) < 1.0e0_knd .and. abs(aimag(r1c)) <  &
              1.0e0_knd) then
          r1c = r1c * 10.0e0_knd
          ir1e = ir1e - 1
          end if
50  continue
        isubr = -int(log10(abs(real(r1temp) / r1rp) + dec))
        if(isubr < 0) isubr = 0
        isubi = 0
        if(r1ip /= 0.0e0_knd) isubi = -int(log10(abs(aimag(r1temp) &
                                   /r1ip) + dec))
        if(isubi < 0) isubi = 0
        isub = isubr
          if(aimag(r1temp) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(r1temp) / real(r1temp))))
          if(iterm < 0) isub = max(isubr, isubi + iterm)
          if(iterm > 0) isub = max(isubi, isubr - iterm)
          end if
        isubdr = -int(log10(abs(real(r1dtemp) / r1drp) + dec))
        if(isubdr < 0) isubdr = 0
        isubdi = 0
        if(r1dip /= 0.0e0_knd) isubdi = -int(log10(abs(aimag(r1dtemp) &
                                     /r1dip) + dec))
        if(isubdi < 0) isubdi = 0
        isubd = isubdr
          if(aimag(r1dtemp) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(r1dtemp) / real(r1dtemp))))
          if(iterm < 0) isubd = max(isubdr, isubdi + iterm)
          if(iterm > 0) isubd = max(isubdi, isubdr - iterm)
          end if
        r1dc1 = r1dcoef * r1c
        r1dc2 = (c * r1dtemp * sbesn(l + 1) * sbesdr(l + 1) * pcoefn/ &
             dmfnorm) * 10.0e0_knd ** (ibese(l + 1) + ipcoefn - ir1e - idmfe)
        r1dc = r1dc1 + r1dc2
          if(m /= 0) then
          iterm = int(log10(abs(r1dc1 / r1dc2) + dec))
          if(iterm < 0) isubd = max(isubd, isub + iterm)
          if(iterm >= 0) isubd = max(isub, isubd - iterm)
          end if
        iterm = int(log10(abs(r1dc2 / r1dc)))
        if(iterm < 0) iterm = 0
        isubd = isubd + iterm
        nsubr1b = max(isub, isubd)
        if(nsubr1b > ndec) nsubr1b = ndec
if (debug) then
        write(40, 60) jbesb, lim, isub, isubd
60      format(8x,'r1besb: numerator conv. in ',i6,' terms with',i6, &
               ' available.',i3,' digits of sub. error',/,15x, &
               ' in r1 numerator. ',i3,' digits of sub. error in', &
               ' r1d numerator.')
end if
        iterm = int(log10(abs(r1dc)))
        ir1de = ir1e + iterm
    r1dc = r1dc * 10.0e0_knd ** (-iterm)
          if(abs(real(r1dc)) > 10.0e0_knd .or. abs(aimag(r1dc)) >  &
              10.0e0_knd) then
          r1dc = r1dc / 10.0e0_knd
          ir1de = ir1de + 1
          end if
          if(abs(real(r1dc)) < 1.0e0_knd .and. abs(aimag(r1dc)) <  &
              1.0e0_knd) then
          r1dc = r1dc * 10.0e0_knd
          ir1de = ir1de - 1
          end if
70  continue
        return
        end subroutine
!
!
        subroutine r2int (l, m, c, x, limint, ndec, maxd, enr, d01, id01, maxint, &
                          maxmp, nex, maxlp, rpint1, rpint2, pint1, pint2, &
                          pint3, pint4, norme, pnorm, ipnorm, coefme, coefmo, &
                          r2c, ir2e, r2dc, ir2de, jint, coefn, icoefn, iflagc)
!
!
!  purpose:     to calculate values of the radial function of the
!               second kind and its first derivative using an integral
!               representation of the radial functions in terms of the
!               angular function of the first kind together with a
!               Neumann function kernal. the angular function is
!               expanded in a series of associated Legendre functions.
!               gaussian quadrature is used (in subroutine pint) to
!               evaluate the resulting integrals involving associated
!               Legendre functions times the Neumann function kernel.
!               this subroutine r2int performs the summation of the
!               integrals times d coefficients to obtain r2 and r2d.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x      : x
!               limint : approximately twice the maximum number of
!                        terms available to be taken in the series
!               ndec   : number of decimal digits available in
!                        real(kind) arithmetic
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               d01    : characteristic of the first d coefficient,
!                        either d0 or d1, depending on whether l-m
!                        is even or odd
!               id01   : exponent (base 10) of the first d coefficient
!               maxint : dimension of pint and rpint arrays
!               maxmp  : dimension of norme array
!               nex    : maximum exponent available in real(knd)
!                        arithmetic
!               maxlp  : maximum  l value desired; dimension
!                        of the pnorm and ipnorm arrays
!               rpint1 : arrays of ratios of successive integrals of
!                        either the first or the third kind, depending
!                        on whether l-m is even or odd
!               rpint2 : array of ratios of successive integrals of
!                        either the second or the fourth kind,
!                        depending on whether l-m is even or odd
!               pint1  : array of scaled values for the integrals of
!                        the first kind
!               pint2  : array of scaled values for the integrals of
!                        the second kind
!               pint3  : array of scaled values for the integrals of
!                        the third kind
!               pint4  : array of scaled values for the integrals of
!                        the fourth kind
!               norme  : scaling exponent for the spherical Neumann
!                        function of order m
!               pnorm  : array of characteristics of the scaling factors
!                        used for the associated Legendre functions in
!                        the integrals to avoid overflow
!               ipnorm : array of exponents (base 10) corresponding to
!                        pnorm
!               coefme : coefficient used to multiply the r2 sum to
!                        obtain its contribution to r2d when l-m is even
!               coefmo : coefficient used to multiply the r2 sum to
!                        obtain its contribution to r2d when l-m is odd
!
!     output:   r2c    : characteristic of prolate radial function
!                        of the first kind
!               ir2e   : exponent of prolate radial function of the
!                        first kind
!               r2dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the first
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        prolate radial function of the first kind
!               jint   : maximum value of the index j in the forward
!                        sum for r2 and r2d, i.e., the highest enr(j)
!                        used
!               coefn  : characteristic of a coefficient calculated
!                        for l = m and used for higher values of l
!               icoefn : exponent for coefn
!               iflagc : equal to zero if have not yet calculated
!                        coefn because this is the first call to r2int;
!                        set equal to one after calculating coefn
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) coefa, coefme, coefmo, coefn, dec, ri, rm, rm2, r2dposi, &
                  r2dposr, r2posi, r2posr, x
        complex(knd) c, coefl, dnew, dnewd, dold, doldd, d01, r2c, r2dc, r2d1, &
                     r2d2, r2dtemp, r2temp
        complex(knd) enr(maxd), pint1(maxint), pint2(maxint), &
                     pint3(maxint), pint4(maxint), rpint1(maxint), &
                     rpint2(maxint)
        real(knd) pnorm(maxlp)
!
!  integer arrays
        integer ipnorm(maxlp)
!
        rm = real(m, knd)
        rm2 = rm + rm
        lm2 = (l - m) / 2
        ix = l - m - 2 * lm2
        ixx = ix - 1
        ixx2 = ixx + 2
        lim = limint / 2 - ix
!
!  compute the leading coefficient
        if(l > m .and. iflagc == 1) go to 20
        icoefn = norme
        coefn = 0.5e0_knd
        iflagc = 1
        if(m == 0) go to 20
          do 10 i = 1, m
          ri = real(i, knd)
          coefn = coefn / (ri + ri)
          iterm = int(log10(abs(coefn)))
          coefn = coefn * 10.0e0_knd ** (-iterm)
10        icoefn = icoefn + iterm
20      continue
        if(ix == 0) coefa = (rm2 + 1.0e0_knd) * coefn
        if(ix == 1) coefa = (rm2 + 3.0e0_knd) * coefn
        if((ix == 0) .and. (2 * (lm2 / 2) /= lm2)) coefa = -coefa
        if((ix == 1) .and. (2 * ((l - m - 1) / 4) /= (l - m - 1) / 2)) coefa = -coefa
        coefl = coefa / d01
        icoefl = -id01 + icoefn
        dec = 10.0e0_knd ** (-ndec - 1)
        jlow = lm2 + 1
!
!  compute the integrals of S1
!
!  forward summation for series for r2 and r2d
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
        r2dtemp = doldd
        r2temp = dold
        r2posr = 1.0e0_knd
        r2posi = 0.0e0_knd
        r2dposr = 1.0e0_knd
        r2dposi = 0.0e0_knd
          do j = jlow, lim
          dnew = dold * enr(j) * rpint1(j + j + ixx2)
          dnewd = doldd * enr(j) * rpint2(j + j + ixx2)
          r2temp = r2temp + dnew
          r2dtemp = r2dtemp + dnewd
          if(real(dnew) > 0.0e0_knd) r2posr = r2posr + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r2dposr = r2dposr + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r2posi = r2posi + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r2dposi = r2dposi + aimag(dnewd)
          if((abs(dnew / r2temp) + abs(dnewd / r2dtemp)) < dec) go to 40
          dold = dnew
          doldd = dnewd
          end do
!
!  backward summation for series for r2 and r2d
40      jint = j
        if(jint > lim) jint = lim
        if(lm2 == 0) go to 60
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
        j = lm2
          do 50 jj = 1, lm2
          dnew = dold / (rpint1(j + j + ixx2) * enr(j))
          dnewd = doldd / (rpint2(j + j + ixx2) * enr(j))
          r2temp = r2temp + dnew
          r2dtemp = r2dtemp + dnewd
          if(real(dnew) > 0.0e0_knd) r2posr = r2posr + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r2dposr = r2dposr + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r2posi = r2posi + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r2dposi = r2dposi + aimag(dnewd)
          if((abs(dnew / r2temp) + abs(dnewd / r2dtemp)) < dec) go to 60
          dold = dnew
          doldd = dnewd
          j = j - 1
50        continue
60      continue
        isubr = -int(log10(abs(real(r2temp) / r2posr) + dec))
        isubi = 0
        if(r2posi /= 0.0e0_knd) isubi = -int(log10(abs(aimag(r2temp) &
                                   /r2posi) + dec))
        isub = isubr
          if(aimag(r2temp) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(r2temp) / real(r2temp))))
          if(iterm < 0) isub = max(isubr, isubi + iterm)
          if(iterm > 0) isub = max(isubi, isubr - iterm)
          end if
        isubdr = -int(log10(abs(real(r2dtemp) / r2dposr) + dec))
        isubdi = 0
        if(r2dposi /= 0.0e0_knd) isubdi = -int(log10(abs(aimag(r2dtemp) &
                                     /r2dposi) + dec))
        isubd = isubdr
          if(aimag(r2dtemp) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(r2dtemp) / real(r2dtemp))))
          if(iterm < 0) isubd = max(isubdr, isubdi + iterm)
          if(iterm > 0) isubd = max(isubdi, isubdr - iterm)
          end if
        r2temp = r2temp * coefl * pnorm(l - m + 1)
        if(ix == 0) r2temp = r2temp * pint1(l - m + 1)
        if(ix == 1) r2temp = r2temp * pint3(l - m + 1) * x
        iterm = int(log10(abs(r2temp)))
        ir2e = iterm + ipnorm(l - m + 1) + icoefl
        r2c = r2temp
        r2c = r2temp * (10.0e0_knd ** (-iterm))
          if(abs(real(r2c)) >= 10.0e0_knd .or. abs(aimag(r2c)) >=  &
              10.0e0_knd) then
          r2c = r2c / 10.0e0_knd
          ir2e = ir2e + 1
          end if
          if(abs(real(r2c)) < 1.0e0_knd .and. abs(aimag(r2c)) <  &
              1.0e0_knd) then
          r2c = r2c * 10.0e0_knd
          ir2e = ir2e - 1
          end if
        r2dtemp = -r2dtemp * coefl * pnorm(l - m + 1) * c * x
          if(ix == 0) then
          r2d1 = r2dtemp * pint2(l - m + 1)
          r2d2 = r2temp * coefme
          end if
          if(ix == 1) then
          r2d1 = r2dtemp * pint4(l - m + 1) * x
          r2d2 = r2temp * coefmo
          end if
        r2dtemp = r2d1 + r2d2
        iterm = int(log10(abs(r2d1 / r2dtemp) + dec))
        if(iterm < 0) iterm = 0
        isubd = isubd + iterm
          if(abs(r2d2) /= 0.0e0_knd) then
          iterm = int(log10(abs(r2d1 / r2d2) + dec))
          if(iterm > 0) isubd = max(isubd, isub - iterm)
          if(iterm <= 0) isubd = max(isub, isubd + iterm)
          end if
        if(isubd > ndec) isubd = ndec
if (debug) then
        write(40, 80) jint, lim, isub, isubd
80      format(8x,'r2int: converged in ',i6,' terms; 'i6, &
               ' available; ',i3,' digits of subtraction error',/, &
               15x,'in r2, ',i3,' digits of subtraction error in r2d.')
end if
        jterm = int(log10(abs(r2dtemp)))
        ir2de = jterm + ipnorm(l - m + 1) + icoefl
        r2dc = r2dtemp * 10.0e0_knd ** (-jterm)
          if(abs(real(r2dc)) > 10.0e0_knd .or. abs(aimag(r2dc)) >  &
              10.0e0_knd) then
          r2dc = r2dc / 10.0e0_knd
          ir2de = ir2de + 1
          end if
          if(abs(real(r2dc)) < 1.0e0_knd .and. abs(aimag(r2dc)) <  &
              1.0e0_knd) then
          r2dc = r2dc * 10.0e0_knd
          ir2de = ir2de - 1
          end if
          if(abs(real(r2dc)) < 1.0e0_knd .and. abs(aimag(r2dc)) <  &
              1.0e0_knd) then
          r2dc = r2dc * 10.0e0_knd
          ir2de = ir2de - 1
          end if
        return
        end subroutine
!
!
        subroutine r2leg (l, m, c, x1, lnum, limleg, limdr, maxd, maxmp, ndec, &
                          nex, idigc, itestm, maxpdr, maxdr, maxq, enr, enrneg, &
                          drhor, nsdrho, d01, id01, dneg, idneg, nsdneg, &
                          dfnorm, idfe, nsubf, dmfnorm, idmfe, nsubmf, prx, &
                          pdrx, qdr, qdml, iqdml, qdl, iqdl, qr, qml, iqml, ql, &
                          iql, fajo, ifajo, termpq, itermpq, ioppsum, &
                          iopqnsum, r1c, ir1e, r1dc, ir1de, naccr1, wront, &
                          minacc, naccrpl, naccr, r2c, ir2e, r2dc, ir2de, jleg, &
                          jlegp, naccleg, nacccor, jflagl, iflagq, iflagp)
!
!  purpose:     to evaluate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x using the traditional expansion in associated
!               Legendre functions.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               x1      : x-1
!               lnum    : number of l values desired
!               limleg  : approximately twice the maximum number
!                         of terms available to be taken in qsum,
!                         (sum involving q's time d coefficients)
!               limdr   : maximum number of terms available to be
!                         taken in psum (sum involving p's time
!                         d rho coefficients)
!               maxd    : dimension of enr array
!               maxmp   : dimension of enrneg array
!               ndec    : number of decimal digits for real(knd)
!               nex     : maximum exponent available for real(knd)
!               idigc   : number of decimal digits of convergence
!                         of the eigenvalue
!               itestm  : match in decimal digits of the forward
!                         and backward recursions to calculate
!                         the d coefficient ratios
!               maxpdr  : dimension of prx and pdrx arrays
!               maxdr   : dimension of drhor array
!               maxq    : dimension of qr and qdr arrays
!               enr     : array of d coefficient ratios
!               enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               nsdrho  : maximum subtraction error in calculating
!                         the array drhor
!               d01     : characteristic of the ratio of the first d
!                         coefficient with nonnegative subscript, either
!                         d0 or d1 depending on whether l-m is even or
!                         odd, to the d coefficient with subscript l - m
!               id01    : exponent (base 10) corresponding to d01
!               dneg    : characteristic of the ratio of the first d
!                         coefficient with negative subscript (either
!                         -1 or -2, depending on whether l - m is odd
!                         or even) to the d coefficient with subscript
!                         l - m
!               idneg   : exponent corresponding to dneg
!               nsdneg  : subtraction error in calculating dneg
!               dfnorm  : characteristic of the Flammer normalization
!                         factor; equal to one over the characteristic
!                         of the d coefficient with subscript l - m
!                         using this normalization for the angular
!                         functions
!               idfe    : exponent associated with dfnorm
!               nsubf   : subtraction error in calculating dfnorm
!               dmfnorm : characteristic of the Morse and Feshbach
!                         normalization factor of the d coefficients.
!                         equal to one over the characteristic of
!                         the d coefficient d(n = l - m) using this
!                         this normalization for the angular functions
!               idmfe   : exponent associated with dmfnorm
!               nsubmf  : subtraction error in calculating dmfnorm
!               prx     : ratios of successive Legendre functions of
!                         the first kind of the same parity
!               pdrx    : ratios of successive first derivatives of
!                         Legendre functions of the first kind of the
!                         same parity
!               qdr     : ratios of first derivatives of successive
!                         Legendre functions of the second kind
!               qdml    : characteristic of the first derivative of
!                         the associated Legendre function of the second
!                         kind with order m and degree m-1
!               iqdml   : exponent corresponding to qdml
!               qdl     : characteristic of the first derivative of
!                         the associated Legendre function of the second
!                         kind with order m and degree m
!               iqdl    : exponent corresponding to qdl
!               qr      : array of ratios of successive associated
!                         Legendre functions of the second kind
!               qml     : characteristic of the associated Legendre
!                         function of the second kind with order m
!                         and degree m-1
!               iqml    : exponent corresponding to qml
!               ql      : characteristic of the associated Legendre
!                         function of the second kind with order m and
!                         degree m
!               iql     : exponent corresponding to ql
!               fajo    : characteristic of the joining factor of the
!                         second kind
!               ifajo   : exponent corresponding to fajo
!               termpq  : characteristic of the relative size of the
!                         maximum terms in the positive degree q series
!                         and the p series used to calculate r2 and r2d
!               itermpq : exponent corresponding to termpq
!               ioppsum : integer flag = 0 if psum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if psum is computed
!               iopqnsum: integer flag = 0 if qnsum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if qnsum is computed
!               r1c     : charcteristic of corresponding radial function
!                         of the first kind
!               ir1e    : exponent of corresponding radial function of
!                         the first kind
!               r1dc    : charcteristic of corresponding first
!                         derivative of the radial function of the first
!                         kind
!               ir1de   : exponent of corresponding first derivative of
!                         the radial function of the first kind
!               naccr1  : estimated accuracy of r1 and r1d
!               wront   : theoretical Wronskian
!               minacc  : minimum desired accuracy in decimal digits
!               naccrpl : estimated subtraction error in forming the
!                         Wronskian
!               naccr   : accuracy for this value of l obtained from
!                         either r2int or naccrpl
!
!     output:   r2c     : characteristic of prolate
!                         radial function of the second kind
!               ir2e    : exponent of prolate radial function of the
!                         second kind
!               r2dc    : characteristic of derivative with
!                         respect to x of prolate radial function
!                         of the second kind
!               ir2de   : exponent of derivative with respect to x of
!                         prolate radial function of second kind
!               jleg    : maximum number of terms taken in qsum
!               jlegp   : maximum number of terms taken in psum
!               naccleg : estimated accuracy of r2 and r2d
!               nacccor : subtraction error in decimal digits that
!                         occurs in forming Wronskian for calculated
!                         radial functions and their first derivatives
!               jflagl  : equal to unity if Wronskian is used to
!                         calculate more accurate values for leading
!                         coefficients and thus more accurate values
!                         for r2 and r2d; equal to zero otherwise
!               iflagq  : set equal to zero at l = m. Remains equal
!                         to zero if set function values and accuracy
!                         to zero because terms in backward qsum
!                         become too large to allow accurate results.
!                         Set equal to unity when this first does not
!                         occur.
!               iflagp  : set equal to zero at l = m. Remains equal
!                         to zero if set function values and accuracy
!                         to zero because terms in psum become too large
!                         to allow accurate results. Set equal to unity
!                         when this first does not occur.
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) dconp, dconq, dconqn, dec, dec1, qdml, qml, rm, ten, &
                  termpq, test, testd, teste, testeo, testm, testdm, &
                  tm, x1
        complex(knd) c, dfnorm, dmfnorm, dneg, dnegjf, dnew, dnewd, dold, doldd, &
                     d01, fac, psum, pdsum, qndsum, qdsum, qnsum, qsum, r1c, &
                     r1dc, r2c, r2dc, spsum, spdsum, wront, wronca, wroncb, &
                     wronc
        real(knd) prx(maxpdr), pdrx(maxpdr), qdl(lnum), qdr(maxq), ql(lnum), &
                  qr(maxq)
        real(knd) r2qposr, r2dqposr, r2qposi, r2dqposi, r2qnposr, r2dqnposr, &
                  r2qnposi, r2dqnposi, r2pposr, r2dpposr, r2pposi, r2dpposi, &
                  sr2pposr, sr2dpposr, sr2pposi, sr2dpposi
        complex(knd) drhor(maxdr), enr(maxd), enrneg(maxmp), fajo(lnum)
!
!  integer arrays
        integer ifajo(lnum), iqdl(lnum), iql(lnum)
!
        jsub = max(nsubf, nsubmf)
        jflagl = 0
        ten = 10.0e0_knd
        dec = ten ** (-ndec - 1)
        dec1 = ten ** (-ndec - 7)
        nfac = nex / 3
        if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
        teste = ten ** (nfac)
        testeo = 1.0e0_knd / teste
        lm2 = (l - m) / 2
        ix = l - m - 2 * lm2
        imxp = m + m + ix
        ixx = 1 - ix
        lim1 = limleg / 2 - ix
        lim2 = limdr - 1
        rm = real(m, knd)
        tm = rm + rm
        dconq = dec
        dconqn = dec
        dconp = dec
        dnegjf = dneg * d01
        if(m == 0) dnegjf = d01
        idnegjf = idneg + id01
        if(m == 0) idnegjf = id01
        fajo(l - m + 1) = fajo(l - m + 1) * dmfnorm * dnegjf / dfnorm
        iterm = int(log10(abs(fajo(l - m + 1))))
        fajo(l - m + 1) = fajo(l - m + 1) * (ten ** (-iterm))
        ifajo(l - m + 1) = ifajo(l - m + 1) + idnegjf + iterm + idmfe - idfe
        iterm = -int(log10(abs(c * x1 * (x1 + 2.0e0_knd) * r1dc)))
        itermq = int(log10(abs(fajo(l - m + 1) * termpq / ql(l - m + 1))))
        itestq = iterm + itermq - ir1de + itermpq + ifajo(l - m + 1) - iql(l - m + 1) + ndec + 3
        itermp = int(log10(abs(fajo(l - m + 1) / (dnegjf * termpq))))
        itestp = iterm + itermp - ir1de - idnegjf - itermpq + ifajo(l - m + 1) + ndec + 3
        nsubqr = 0
        nsubqi = 0
        nsubqdr = 0
        nsubqdi = 0
        nsubqnr = 0
        nsubqni = 0
        nsubqndr = 0
        nsubqndi = 0
        nsubpr = 0
        nsubpi = 0
        nsubpdr = 0
        nsubpdi = 0
!
!  begin calculation of series for r2
!
!  calculate d*q sum over positive n using pyramid summation
!
!  backward summation
        qsum = (1.0e0_knd, 0.0e0_knd)
        qdsum = (1.0e0_knd, 0.0e0_knd)
        r2qposr = 1.0e0_knd
        r2qposi = 0.0e0_knd
        r2dqposr = 1.0e0_knd
        r2dqposi = 0.0e0_knd
        iscale = 0
        if(lm2 < 1) go to 20
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
        j = lm2
          do 10 jj = 1, lm2
          dnew = dold / (qr(j + j + imxp) * qr(j + j + imxp - 1) * enr(j))
          qsum = qsum + dnew
          dnewd = doldd / (qdr(j + j + imxp) * qdr(j + j + imxp - 1) * enr(j))
          qdsum = qdsum + dnewd
          if(real(dnew) > 0.0e0_knd) r2qposr = r2qposr + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r2dqposr = r2dqposr + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r2qposi = r2qposi + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r2dqposi = r2dqposi + aimag(dnewd)
            if(int(log10(abs(qsum))) + iscale > itestq .and. iflagq == 0) &
                then
            r2c = 0.0e0_knd
            r2dc = 0.0e0_knd
            ir2e = 0
            ir2de = 0
            nsub = ndec
            nsubd = ndec
            jleg = j
            jlegp = 0
            naccleg = 0
            go to 180
            end if
            if(abs(qsum) > teste) then
            dnew = dnew * testeo
            dnewd = dnewd * testeo
            qsum = qsum * testeo
            qdsum = qdsum * testeo
            r2qposr = r2qposr * testeo
            r2dqposr = r2dqposr * testeo
            r2qposi = r2qposi * testeo
            r2dqposi = r2dqposi * testeo
            iscale = iscale + nfac
            end if
          if((abs(dnew / qsum) + abs(dnewd / qdsum)) < dconq) go to 20
          dold = dnew
          doldd = dnewd
          j = j - 1
10        continue
20      continue
!
!  forward summation
          iflagq = 1
          if(iscale /= 0) then
          jleg = lm2
          go to 45
          end if
        jlow = lm2 + 1
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
          do 30 j = jlow, lim1
          dnew = dold * enr(j) * qr(j + j + imxp) * qr(j + j + imxp - 1)
          qsum = qsum + dnew
          dnewd = doldd * enr(j) * qdr(j + j + imxp) * qdr(j + j + imxp - 1)
          qdsum = qdsum + dnewd
          if(real(dnew) > 0.0e0_knd) r2qposr = r2qposr + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r2dqposr = r2dqposr + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r2qposi = r2qposi + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r2dqposi = r2dqposi + aimag(dnewd)
          if((abs(dnew / qsum) + abs(dnewd / qdsum)) < dconq) go to 40
          dold = dnew
          doldd = dnewd
30        continue
40      continue
        jleg = j
        if(jleg > lim1) jleg = lim1
45      itestqsum = -int(log10(abs(dnew / qsum) + dec))
        if(real(qsum) /= 0.0e0_knd) &
           nsubqr = int(log10(abs(r2qposr / real(qsum) + dec)))
        if(nsubqr < 0) nsubqr = 0
        if(aimag(qsum) /= 0.0e0_knd) &
           nsubqi = int(log10(abs(r2qposi / aimag(qsum) + dec)))
        if(nsubqi < 0) nsubqi = 0
        nsubq = nsubqr
          if(aimag(qsum) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(qsum) / real(qsum))))
          if(iterm < 0) nsubq = max(nsubqr, nsubqi + iterm)
          if(iterm > 0) nsubq = max(nsubqi, nsubqr - iterm)
          end if
        if(real(qdsum) /= 0.0e0_knd) &
           nsubqdr = int(log10(abs(r2dqposr / real(qdsum) + dec)))
        if(nsubqdr < 0) nsubqdr = 0
        if(aimag(qdsum) /= 0.0e0_knd) &
           nsubqdi = int(log10(abs(r2dqposi / aimag(qdsum) + dec)))
        if(nsubqdi < 0) nsubqdi = 0
        nsubqd = nsubqdr
          if(aimag(qdsum) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(qdsum) / real(qdsum))))
          if(iterm < 0) nsubqd = max(nsubqdr, nsubqdi + iterm)
          if(iterm > 0) nsubqd = max(nsubqdi, nsubqdr - iterm)
          end if
        qsum = qsum * ql(l - m + 1) / (fajo(l - m + 1) * termpq)
        iterm = int(log10(abs(qsum)))
        qsum = qsum * (10.0e0_knd ** (-iterm))
        iqsum = iql(l - m + 1) - ifajo(l - m + 1) - itermpq + iterm + iscale
        qdsumsav = qdsum
        qdsum = qdsum * qdl(l - m + 1) / (fajo(l - m + 1) * termpq)
        iterm = int(log10(abs(qdsum)))
        qdsum = qdsum * (10.0e0_knd ** (-iterm))
        iqdsum = iqdl(l - m + 1) - ifajo(l - m + 1) - itermpq + iterm + iscale
        qdsum = qdsum * (10.0e0_knd ** (iqdsum - iqsum))
!
!  calculate d*q sum over negative n
        qnsum = (0.0e0_knd, 0.0e0_knd)
        qndsum = (0.0e0_knd, 0.0e0_knd)
        r2qnposr = 0.0e0_knd
        r2dqnposr = 0.0e0_knd
        r2qnposi = 0.0e0_knd
        r2dqnposi = 0.0e0_knd
        iqnsum = 0
        iqndsum = 0
        nsubqn = 0
        nsubqnd = 0
        nmterm = 0
        if(iopqnsum == 0) go to 90
        qnsum = enrneg(m)
        qndsum = enrneg(m)
        if(ix == 1) go to 50
        qnsum = qnsum * qr(m + m - 1)
        qndsum = qndsum * qdr(m + m - 1)
50      continue
        if(real(qnsum) > 0.0e0_knd) r2qnposr = real(qnsum)
        if(real(qndsum) > 0.0e0_knd) r2dqnposr = real(qndsum)
        if(aimag(qnsum) > 0.0e0_knd) r2qnposi = aimag(qnsum)
        if(aimag(qndsum) > 0.0e0_knd) r2dqnposi = aimag(qndsum)
        if(m == 1) go to 80
        dold = qnsum
        doldd = qndsum
          do 60 j = 2, m
          dnew = dold * enrneg(m - j + 1) * qr(imxp - j - j + 1) * qr(imxp - j - j + 2)
          qnsum = qnsum + dnew
          dnewd = doldd * enrneg(m - j + 1) * qdr(imxp - j - j + 1) * qdr(imxp - j - j + 2)
          qndsum = qndsum + dnewd
          if(real(dnew) > 0.0e0_knd) r2qnposr = r2qnposr + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r2dqnposr = r2dqnposr + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r2qnposi = r2qnposi + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r2dqnposi = r2dqnposi + aimag(dnewd)
          dold = dnew
60        doldd = dnewd
80      continue
        if(real(qnsum) /= 0.0e0_knd) &
           nsubqnr = int(log10(abs(r2qnposr / real(qnsum) + dec)))
        if(nsubqnr < 0) nsubqnr = 0
        if(aimag(qnsum) /= 0.0e0_knd) &
           nsubqni = int(log10(abs(r2qnposi / aimag(qnsum) + dec)))
        if(nsubqni < 0) nsubqni = 0
        nsubqn = nsubqnr
          if(aimag(qnsum) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(qnsum) / real(qnsum))))
          if(iterm < 0) nsubqn = max(nsubqnr, nsubqni + iterm)
          if(iterm > 0) nsubqn = max(nsubqni, nsubqnr - iterm)
          end if
        nsubqn = max(nsubqn, nsdneg)
        if(real(qndsum) /= 0.0e0_knd) &
           nsubqndr = int(log10(abs(r2dqnposr / real(qndsum) + dec)))
        if(nsubqndr < 0) nsubqndr = 0
        if(aimag(qndsum) /= 0.0e0_knd) &
           nsubqndi = int(log10(abs(r2dqnposi / aimag(qndsum) + dec)))
        if(nsubqndi < 0) nsubqndi = 0
        nsubqnd = nsubqndr
          if(aimag(qndsum) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(qndsum) / real(qndsum))))
          if(iterm < 0) nsubqnd = max(nsubqndr, nsubqndi + iterm)
          if(iterm > 0) nsubqnd = max(nsubqndi, nsubqndr - iterm)
          end if
        nsubqnd = max(nsubqnd, nsdneg)
        qnsum = qnsum * qml * d01 / (fajo(l - m + 1) * termpq)
        iterm = int(log10(abs(qnsum)))
        qnsum = qnsum * (10.0e0_knd ** (-iterm))
        iqnsum = iqml + id01 - ifajo(l - m + 1) - itermpq + iterm
        qnsum = qnsum * (10.0e0_knd ** (iqnsum - iqsum))
        qndsumsav = qndsum
        qndsum = qndsum * qdml * d01 / (fajo(l - m + 1) * termpq)
        iterm = int(log10(abs(qndsum)))
        qndsum = qndsum * (10.0e0_knd ** (-iterm))
        iqndsum = iqdml + id01 - ifajo(l - m + 1) - itermpq + iterm
        qndsum = qndsum * (10.0e0_knd ** (iqndsum - iqsum))
90      continue
!
!       calculate d(rho|n)*p summation
        psum = (0.0e0_knd, 0.0e0_knd)
        pdsum = (0.0e0_knd, 0.0e0_knd)
        iscale = 0
        ipsum = 0
        ipdsum = 0
        r2pposr = 0.0e0_knd
        r2dpposr = 0.0e0_knd
        r2pposi = 0.0e0_knd
        r2dpposi = 0.0e0_knd
        nsubp = 0
        nsubpd = 0
        jlegp = 0
        itestpsum = 0
        if(ioppsum == 0) go to 160
        psum = prx(ixx + 1) * drhor(1)
        pdsum = pdrx(ixx + 1) * drhor(1)
        dold = psum
        doldd = pdsum
        if(m /= 0 .or. ix /= 1) go to 100
        pdsum = (0.0e0_knd, 0.0e0_knd)
        doldd = drhor(1)
100     continue
        spsum = psum
        spdsum = pdsum
        testm = 1.0e0_knd
        testdm = 1.0e0_knd
        jlegpf = 1
        jlegpd = 1
        if(real(psum) > 0.0e0_knd) r2pposr = real(psum)
        if(real(pdsum) > 0.0e0_knd) r2dpposr = real(pdsum)
        if(aimag(psum) > 0.0e0_knd) r2pposi = aimag(psum)
        if(aimag(pdsum) > 0.0e0_knd) r2dpposi = aimag(pdsum)
        sr2pposr = r2pposr
        sr2pposi = r2pposi
        sr2dpposr = r2dpposr
        sr2dpposi = r2dpposi
          do 130 j = 2, lim2
          dnew = dold * drhor(j) * prx(j + j - ix)
          psum = psum + dnew
          dnewd = doldd * drhor(j) * pdrx(j + j - ix)
          pdsum = pdsum + dnewd
          test = abs(dnew / psum)
          testd = abs(dnewd / pdsum)
          if(real(dnew) > 0.0e0_knd) r2pposr = r2pposr + real(dnew)
          if(real(dnewd) > 0.0e0_knd) r2dpposr = r2dpposr + real(dnewd)
          if(aimag(dnew) > 0.0e0_knd) r2pposi = r2pposi + aimag(dnew)
          if(aimag(dnewd) > 0.0e0_knd) r2dpposi = r2dpposi + aimag(dnewd)
          if(test > testm .or. test == 0.0e0_knd) go to 110
          spsum = psum
          sr2pposr = r2pposr
          sr2pposi = r2pposi
          testm = test
          jlegpf = j
110       if(testd > testdm .or. testd == 0.0e0_knd) go to 120
          testdm = testd
          spdsum = pdsum
          sr2dpposr = r2dpposr
          sr2dpposi = r2dpposi
          jlegpd = j
120       continue
            if(int(log10(abs(qsum))) + iscale > itestp .and. iflagp == 0) &
                then
            r2c = 0.0e0_knd
            r2dc = 0.0e0_knd
            ir2e = 0
            ir2de = 0
            nsub = ndec
            nsubd = ndec
            jlegp = 0
            naccleg = 0
            go to 180
            end if
            if(abs(psum) > teste) then
            dnew = dnew * testeo
            dnewd = dnewd * testeo
            psum = psum * testeo
            pdsum = pdsum * testeo
            r2pposr = r2pposr * testeo
            r2dpposr = r2dpposr * testeo
            r2pposi = r2pposi * testeo
            r2dpposi = r2dpposi * testeo
            spsum = spsum * testeo
            sr2pposr = sr2pposr * testeo
            sr2pposi = sr2pposi * testeo
            spdsum = spdsum * testeo
            sr2dpposr = sr2dpposr * testeo
            sr2dpposi = sr2dpposi * testeo
            iscale = iscale + nfac
            end if
          if((test + testd) < dconp) go to 150
          dold = dnew
          doldd = dnewd
130       continue
150     iflagp = 1
        psum = spsum
        pdsum = spdsum
        r2pposr = sr2pposr
        r2pposi = sr2pposi
        r2dpposr = sr2dpposr
        r2dpposi = sr2dpposi
        jlegp = max(jlegpf, jlegpd)
        ntestm = -int(log10(testm))
        ntestdm = -int(log10(testdm))
        if(ntestm > ndec) ntestm = ndec
        if(ntestdm > ndec) ntestdm = ndec
        test = max(testm, testdm)
        itestpsum = -int(log10(test + dec))
        if(itestpsum > ndec) itestpsum = ndec
        if(itestpsum < 0) itestpsum = 0
        if(real(psum) /= 0.0e0_knd) &
           nsubpr = int(log10(abs(r2pposr / real(psum) + dec)))
        if(nsubpr < 0) nsubpr = 0
        if(aimag(psum) /= 0.0e0_knd) &
           nsubpi = int(log10(abs(r2pposi / aimag(psum) + dec)))
        if(nsubpi < 0) nsubpi = 0
        nsubp = nsubpr
          if(aimag(psum) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(psum) / real(psum))))
          if(iterm < 0) nsubp = max(nsubpr, nsubpi + iterm)
          if(iterm > 0) nsubp = max(nsubpi, nsubpr - iterm)
          end if
        nsubp = max(nsubp, nsdrho, ndec - ntestm)
        if(real(pdsum) /= 0.0e0_knd) &
           nsubpdr = int(log10(abs(r2dpposr / real(pdsum) + dec)))
        if(nsubpdr < 0) nsubpdr = 0
        if(aimag(pdsum) /= 0.0e0_knd) &
           nsubpdi = int(log10(abs(r2dpposi / aimag(pdsum) + dec)))
        if(nsubpdi < 0) nsubpdi = 0
        nsubpd = nsubpdr
          if(aimag(pdsum) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(pdsum) / real(pdsum))))
          if(iterm < 0) nsubpd = max(nsubpdr, nsubpdi + iterm)
          if(iterm > 0) nsubpd = max(nsubpdi, nsubpdr - iterm)
          end if
        nsubpd = max(nsubpd, nsdrho, ndec - ntestdm)
        psum = psum * dnegjf * termpq / fajo(l - m + 1)
        iterm = 0
        if(psum /= (0.0e0_knd, 0.0e0_knd)) iterm = int(log10(abs(psum)))
        psum = psum * (10.0e0_knd ** (-iterm))
        ipsum = idnegjf + itermpq - ifajo(l - m + 1) + iterm + iscale
        psum = psum * (10.0e0_knd ** (ipsum - iqsum))
        pdsum = pdsum * dnegjf * termpq / fajo(l - m + 1)
        if(m /= 0) then
        pdsum = pdsum * rm * (x1 + 1.0e0_knd) / (x1 * (x1 + 2.0e0_knd))
        end if
        iterm = 0
        if(pdsum /= 0.0e0_knd) iterm = int(log10(abs(pdsum)))
        pdsum = pdsum * (10.0e0_knd ** (-iterm))
        ipdsum = idnegjf + itermpq - ifajo(l - m + 1) + iterm + iscale
        pdsum = pdsum * (10.0e0_knd ** (ipdsum - iqsum))
!
160     r2c = qsum + qnsum + psum
        r2dc = qdsum + qndsum + pdsum
        nqs = 0
        if(qsum / r2c /= (0.0e0_knd, 0.0e0_knd)) &
                      nqs = int(log10(abs(qsum / r2c)))
        nqns = 0
        if(qnsum / r2c /= (0.0e0_knd, 0.0e0_knd)) &
                      nqns = int(log10(abs(qnsum / r2c)))
        nps = 0
        if(psum / r2c /= (0.0e0_knd, 0.0e0_knd)) &
                      nps = int(log10(abs(psum / r2c)))
        nsubq = nsubq + nqs
        if(nsubq < 0) nsubq = 0
        if(nsubq > ndec) nsubq = ndec
        if(qsum / r2c == (0.0e0_knd, 0.0e0_knd)) nsubq = 0
        nsubqn = nsubqn + nqns
        if(nsubqn < 0) nsubqn = 0
        if(nsubqn > ndec) nsubqn = ndec
        if(qnsum / r2c == (0.0e0_knd, 0.0e0_knd)) nsubqn = 0
        nsubp = nsubp + nps
        if(nsubp < 0) nsubp = 0
        if(nsubp > ndec) nsubp = ndec
        if(psum / r2c == (0.0e0_knd, 0.0e0_knd)) nsubp = 0
        nsub = max(nsubq, nsubqn, nsubp)
        nqds = 0
        if(qdsum / r2dc /= (0.0e0_knd, 0.0e0_knd)) &
                      nqds = int(log10(abs(qdsum / r2dc)))
        nqnds = 0
        if(qndsum / r2dc /= (0.0e0_knd, 0.0e0_knd)) &
                      nqnds = int(log10(abs(qndsum / r2dc)))
        npds = 0
        if(pdsum / r2dc /= (0.0e0_knd, 0.0e0_knd)) &
                      npds = int(log10(abs(pdsum / r2dc)))
        nsubqd = nsubqd + nqds
        if(nsubqd < 0) nsubqd = 0
        if(nsubqd > ndec) nsubqd = ndec
        if(qdsum / r2dc == (0.0e0_knd, 0.0e0_knd)) nsubqd = 0
        nsubqnd = nsubqnd + nqnds
        if(nsubqnd < 0) nsubqnd = 0
        if(nsubqnd > ndec) nsubqnd = ndec
        if(qndsum / r2dc == (0.0e0_knd, 0.0e0_knd)) nsubqnd = 0
        nsubpd = nsubpd + npds
        if(nsubpd < 0) nsubpd = 0
        if(nsubpd > ndec) nsubpd = ndec
        if(pdsum / r2dc == (0.0e0_knd, 0.0e0_knd)) nsubpd = 0
        nsubd = max(nsubqd, nsubqnd, nsubpd)
        if(qnsum == (0.0e0_knd, 0.0e0_knd) .and. qndsum ==  &
                    (0.0e0_knd, 0.0e0_knd)) iopqnsum = 0
        if(psum == (0.0e0_knd, 0.0e0_knd) .and. pdsum ==  &
                    (0.0e0_knd, 0.0e0_knd)) ioppsum = 0
        wronca = r1c * r2dc * 10.0e0_knd ** (ir1e+iqsum)
        wroncb = r2c * r1dc * 10.0e0_knd ** (iqsum + ir1de)
        wronc = wronca - wroncb
        naccleg = -int(log10(abs((wronc - wront) / wront) + dec))
        if(naccleg < 0) naccleg = 0
        if(naccleg > ndec-1) naccleg = ndec-1
        nacclego = naccleg
        nacccor = -int(log10(abs(wronc / wronca) + dec))
        if(nacccor < 0) nacccor = 0
        nacccor = min(nacccor, naccrpl + 1)
        if(naccleg > 1) naccleg = min(naccleg + nacccor, naccr1)
        itest = min(idigc, itestm) - 3 - nacccor - max(nsub, nsubd)
        if(itest < 0) itest = 0
          if(naccleg < minacc .and. naccleg < itest .and. itest >  &
               naccr .and. x1 <= 0.01e0_knd) then
          fac = wront / wronc
          qsum = qsum * fac
          qdsum = qdsum * fac
          qnsum = qnsum * fac
          qndsum = qndsum * fac
          psum = psum * fac
          pdsum = pdsum * fac
          r2c = qsum + qnsum + psum
          r2dc = qdsum + qndsum + pdsum
          jflagl = 1
          naccleg = itest
          end if
        if(naccleg > 3 .and. nps < (-ndec - 1) .and. npds < (-ndec - 1)) &
             ioppsum = 0
        if(naccleg > 3 .and. nqns < (-ndec - 1) .and. nqnds < (-ndec - 1)) &
            iopqnsum = 0
        if(naccleg < 0) naccleg = 0
          if(jflagl == 0) then
          nsub = max(nsub, jsub, nsdneg)
          nsubd = max(nsubd, jsub, nsdneg)
          end if
        iterm = int(log10(abs(r2c)))
        r2c = r2c * (10.0e0_knd ** (-iterm))
        ir2e = iqsum + iterm
          if(abs(real(r2c)) > 10.0e0_knd .or. abs(aimag(r2c)) >  &
              10.0e0_knd) then
          r2c = r2c / 10.0e0_knd
          ir2e = ir2e + 1
          end if
          if(abs(real(r2c)) < 1.0e0_knd .and. abs(aimag(r2c)) <  &
              1.0e0_knd) then
          r2c = r2c * 10.0e0_knd
          ir2e = ir2e - 1
          end if
        iterm = int(log10(abs(r2dc)))
        r2dc = r2dc * (10.0e0_knd ** (-iterm))
        ir2de = iqsum + iterm
          if(abs(real(r2dc)) > 10.0e0_knd .or. abs(aimag(r2dc)) >  &
              10.0e0_knd) then
          r2dc = r2dc / 10.0e0_knd
          ir2de = ir2de + 1
          end if
          if(abs(real(r2dc)) < 1.0e0_knd .and. abs(aimag(r2dc)) <  &
              1.0e0_knd) then
          r2dc = r2dc * 10.0e0_knd
          ir2de = ir2de - 1
          end if
180     continue
if (debug) then
        if(ioppsum == 1 .and. iopqnsum == 1) write(40, 190) jleg, jlegp, &
                                            m, lim1, lim2, m, nsub, nsubd
190     format(8x,'r2leg: qsum, psum and qnsum series converged in ',i6, &
              ',' i6,' and ',i4,' terms; ',i6,',' i6,' and ' i4, &
              ' terms avail.',/,15x, i2,' and ',i2,' digits of sub.', &
              ' error in r2 and r2d.')
        if(ioppsum == 1 .and. iopqnsum == 0) write(40, 200) jleg, jlegp, &
                                           lim1, lim2, nsub, nsubd
200     format(8x,'r2leg: qsum and psum series converged in ',i6, &
              ' and ',i6,' terms; ',i6,' and ',i6,' terms avail.',/, &
              15x, i2,' and ',i2,' digits of sub. error in r2 and r2d;', &
              ' qnsum is negligible.')
        if(ioppsum == 0 .and. iopqnsum == 1) write(40, 210) jleg, m, &
                                           lim1, m, nsub, nsubd
210     format(8x,'r2leg: qsum and qnsum series converged in ',i6, &
              ' and ',i4,' terms; ',i6,' and ',i4,' terms avail.',/, &
               15x, i2,' and ',i2,' digits of sub. error in r2 and r2d;' &
               ' psum is negligible.')
        if(ioppsum == 0 .and. iopqnsum == 0) write(40, 220) jleg, lim1, &
                                                  nsub, nsubd
220     format(8x,'r2leg: qsum series converged in ',i6,' terms with ', &
               i6,' terms avail.; 'i2,' and ',i2,' digits of',/,15x, &
               'sub. error in r2 and r2d; psum and qnsum are ', &
               'negligible.')
        if(jflagl == 1) write(40, 230)
230     format(15x,'Wronskian used to improve accuracy of the', &
               ' joining factor including dfnorm, dmfnorm and dneg.')
end if
        nacccor = min(nacccor, naccleg - nacclego, naccleg)
        if(nacccor < 0) nacccor = 0
        return
        end subroutine
!
!
        subroutine r2neu (l, m, c, x1, limneu, maxd, maxlp, ndec, nex, maxn, maxp, &
                          minacc, enr, sneuf, sneun, ineue, sneudf, &
                          sneudr, prat1, pcoefn, ipcoefn, dmfnorm, idmfe, &
                          r1dc, ir1de, naccrpl, r2c, ir2e, r2dc, ir2de, jneu)
!
!  purpose:     to calculate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x, using the traditional expansions in terms of
!               spherical Neumann functions.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : complex c
!               x1     : x-1
!               limneu : maximum number of terms to be taken in the
!                        series summations for r2 and r2d
!               maxd   : dimension of enr array
!               maxlp  : maximum  l value desired; dimension
!                        of the sneun, sneudn, ineue, and ineude arrays
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent available for real(knd)
!               maxn   : dimension of sneuf and sneudf arrays
!               maxp   : dimension of prat1 array
!               minacc : number of decimal digits of desired accuracy
!                        of the resulting radial functions
!               enr    : array of ratios of successive d coefficients
!               sneuf  : array of ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : array of characteristics for Neumann functions
!               ineue  : array of exponents for Neumann functions
!               sneudf : array of ratios of successive first derivatives
!                        of spherical Neumann functions of same parity
!               sneudr : array of ratios of first derivatives of Neumann
!                        functions to the corresponding functions
!               prat1  : array of ratios of successive coefficients in
!                        r2 and r2d sum
!               pcoefn : characteristic of coefficient for term in both
!                        r2 and r2d sums that contain the Neumann function
!                        of order l + m
!               ipcoefn: exponent (to the base 10) corresponding to
!                        pcoefn
!               dmfnorm: characteristic of the Morse and Feshbach
!                        normalization. equal to the reciprocal of the
!                        characteristic of the d coefficient with
!                        subscript l - m using this normalization for
!                        the angular functions
!               idmfe  : exponent associated with dmfnorm
!               r1dc   : charcteristic of corresponding first
!                        derivative of the radial function of the first
!                        kind
!               ir1de  : exponent of corresponding first derivative of
!                        the radial function of the first kind
!               naccrpl: expected subtraction error in forming the
!                        Wronskian
!
!     output:   r2c    : characteristic of prolate radial function
!                        of the second kind
!               ir2e   : exponent of prolate radial function of the
!                        second kind
!               r2dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the second
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        prolate radial function of the second kind
!               jneu   : index of term where best convergence is
!                        achieved for r2 or for r2d, whichever term is
!                        larger
!
        use param
!
!
!  real(knd) and complex(knd) scalars and arrays
        real(knd) dconb, dconf, dconi, dec, pcoefn, rm, rm2, r2dcoef, r2est, &
                  r2test, sumdpi, sumdpr, sumpi, sumpr, sumdpit, sumdprt, &
                  sumpit, sumprt, test, testa, testd, testdm, teste, &
                  testeo, testm, testr2, x1
        complex(knd) c, dmfnorm, dnew, dnewd, dold, doldd, r1dc, r2c, r2dc, &
                     r2dtemp, r2temp, sr2temp, sr2dtemp, sumcoef
        real(knd) prat1(maxp)
        complex(knd) enr(maxd), sneudr(maxlp), sneun(maxn), &
                     sneuf(maxn), sneudf(maxn)
!
!  integer arrays
        integer ineue(maxn)
!
        nfac = nex / 3
        if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
        teste = 10.0e0_knd ** (nfac)
        testeo = 1.0e0_knd / teste
        iscale = 0
        rm = real(m, knd)
        dconi = 10.0e0_knd ** (ndec + 2)
        dec = 10.0e0_knd ** (-ndec)
        dconf = 10.0e0_knd ** (-minacc - 4)
        sumcoef = (10.0e0_knd ** (-ir1de-ineue(l + 1) - ipcoefn + naccrpl))/ &
                (c * x1 * (x1 + 2.0e0_knd) * r1dc * sneun(l + 1) * pcoefn)
        r2est = abs(sumcoef * dmfnorm) * 10.0e0_knd ** (idmfe)
        dconb = r2est / dconi
        r2test = r2est * dconi
        testr2 = r2est * 0.001e0_knd
        r2dcoef = rm / ((x1 + 1.0e0_knd) * (x1 + 2.0e0_knd) * x1)
        rm2 = rm * 2.0e0_knd
        lm2 = (l - m) / 2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix = l - m - 2 * lm2
        lim = limneu / 2 - ix
!
!  compute radial function of the second kind
!
!  backward series
        r2temp = (1.0e0_knd, 0.0e0_knd)
        sumpr = 1.0e0_knd
        sumpi = 0.0e0_knd
        r2dtemp = (1.0e0_knd, 0.0e0_knd)
        sumdpr = 1.0e0_knd
        sumdpi = 0.0e0_knd
        if (lm2 < 1) go to 20
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
          do 10 j = lm2, 1,-1
          jj = j + j + ix
          dnew = -dold / (sneuf(jj + m) * prat1(jj + 1) * enr(j))
          dnewd = -doldd / (sneudf(jj + m) * prat1(jj + 1) * enr(j))
          r2temp = r2temp + dnew
          if(real(dnew) > 0.0e0_knd) sumpr = sumpr + real(dnew)
          if(aimag(dnew) > 0.0e0_knd) sumpi = sumpi + aimag(dnew)
          r2dtemp = r2dtemp + dnewd
          if(real(dnewd) > 0.0e0_knd) sumdpr = sumdpr + real(dnewd)
          if(aimag(dnewd) > 0.0e0_knd) sumdpi = sumdpi + aimag(dnewd)
          if(abs(dnew / r2temp) + abs(dnewd / r2dtemp) < dconb) go to 20
          dold = dnew
          doldd = dnewd
10      continue
20      continue
!
!  forward series
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd = (1.0e0_knd, 0.0e0_knd)
        testm = 1.0e0_knd
        testdm = 1.0e0_knd
        sr2temp = r2temp
        sumprt = sumpr
        sumpit = sumpi
        sr2dtemp = r2dtemp
        sumdprt = sumdpr
        sumdpit = sumdpi
        js = lim
        jds = lim
          do 70 j = lm2 + 1, lim
          jj = j + j + ix
          dnew = -dold * enr(j) * sneuf(jj + m) * prat1(jj + 1)
          dnewd = -doldd * enr(j) * sneudf(jj + m) * prat1(jj + 1)
          r2temp = r2temp + dnew
          if(real(dnew) > 0.0e0_knd) sumpr = sumpr + real(dnew)
          if(aimag(dnew) > 0.0e0_knd) sumpi = sumpi + aimag(dnew)
          r2dtemp = r2dtemp + dnewd
          if(real(dnewd) > 0.0e0_knd) sumdpr = sumdpr + real(dnewd)
          if(aimag(dnewd) > 0.0e0_knd) sumdpi = sumdpi + aimag(dnewd)
          test = abs(dnew / r2temp)
          testd = abs(dnewd / r2dtemp)
          testa = abs(r2temp)
          if(max(testa, sumpr, sumpi) < testr2) go to 65
          if(testa > r2test) go to 80
          if(test < testm) go to 30
          go to 40
30        testm = test
          sr2temp = r2temp
          sumprt = sumpr
          sumpit = sumpi
          js = j
40        continue
          if(testd < testdm) go to 50
          go to 60
50        testdm = testd
          sr2dtemp = r2dtemp
          sumdprt = sumdpr
          sumdpit = sumdpi
          jds = j
60        continue
          if((test + testd) < dconf) go to 90
          if(abs(r2temp) > teste) then
            r2temp = r2temp * testeo
            r2dtemp = r2dtemp * testeo
            sr2temp = sr2temp * testeo
            sr2dtemp = sr2dtemp * testeo
            sumpr = sumpr * testeo
            sumdpr = sumdpr * testeo
            sumpi = sumpi * testeo
            sumdpi = sumdpi * testeo
            dnew = dnew * testeo
            dnewd = dnewd * testeo
            iscale = iscale + nfac
            r2test = r2test * testeo
            testr2 = testr2 * testeo
            sumprt = sumprt * testeo
            sumdprt = sumdprt * testeo
            sumpit = sumpit * testeo
            sumdpit = sumdpit * testeo
            end if
65        dold = dnew
          doldd = dnewd
70        continue
80      r2temp = sr2temp
        r2dtemp = sr2dtemp
        sumpr = sumprt
        sumpi = sumpit
        sumdpr = sumdprt
        sumdpi = sumdpit
90      continue
        jneu = max(js, jds)
        if(j >= lim) jneu = jneu + int(abs(c))
        jtestm = -int(log10(testm + dec))
        if(jtestm < 0) jtestm = 0
        if(jtestm > ndec) jtestm = ndec
        jtestdm = -int(log10(testdm + dec))
        if(jtestdm < 0) jtestdm = 0
        if(jtestdm > ndec) jtestdm = ndec
        naccns1r = -int(log10(abs(real(r2temp) / sumpr) + dec))
        naccns1i = 0
        if(sumpi /= 0.0e0_knd) naccns1i = -int(log10(abs(aimag(r2temp) &
                                  /sumpi) + dec))
        if(naccns1r < 0) naccns1r = 0
        if(naccns1r > ndec) naccns1r = ndec
        if(naccns1i < 0) naccns1i = 0
        if(naccns1i > ndec) naccns1i = ndec
        naccns2r = -int(log10(abs(real(r2dtemp) / sumdpr) + dec))
        naccns2i = 0
        if(sumdpi /= 0.0e0_knd) naccns2i = -int(log10(abs(aimag(r2dtemp) &
                                  /sumdpi) + dec))
        if(naccns2r < 0) naccns2r = 0
        if(naccns2r > ndec) naccns2r = ndec
        if(naccns2i < 0) naccns2i = 0
        if(naccns2i > ndec) naccns2i = ndec
        naccns1 = max(naccns1r, naccns1i)
        naccns2 = max(naccns2r, naccns2i)
        naccns = max(naccns1, naccns2)
if (debug) then
        write(40, 100) j, lim, js, jtestm, naccns1, jds, jtestdm, naccns2
100     format(8x,'r2neu: numerator converged in ',i6,' terms; ',i6, &
               ' terms available. best r2 at ',i6,' terms with', &
               /15x,'convergence to',i3,' digits;',i3,' digits', &
               ' subtr. error. best r2d at ',i6,' terms with',/,15x, &
               'convergence to ',i2,' digits;',i3,' digits subtr.', &
               ' error.')
end if
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2c = r2temp * sneun(l + 1) * pcoefn / dmfnorm
        iterm = int(log10(abs(r2c)))
        ir2e = ineue(l + 1) + ipcoefn + iterm - idmfe + iscale
        r2c = r2c * 10.0e0_knd ** (-iterm)
          if(abs(real(r2c)) > 10.0e0_knd .or. abs(aimag(r2c)) >  &
              10.0e0_knd) then
          r2c = r2c / 10.0e0_knd
          ir2e = ir2e + 1
          end if
          if(abs(real(r2c)) < 1.0e0_knd .and. abs(aimag(r2c)) <  &
              1.0e0_knd) then
          r2c = r2c * 10.0e0_knd
          ir2e = ir2e - 1
          end if
        r2dc = r2dcoef * r2c + (c * r2dtemp * sneun(l + 1) * sneudr(l + 1) * pcoefn/ &
             dmfnorm) * 10.0e0_knd ** (ineue(l + 1) + ipcoefn - idmfe+iscale-ir2e)
        iterm = int(log10(abs(r2dc)))
        ir2de = ir2e + iterm
        r2dc = r2dc * 10.0e0_knd ** (-iterm)
          if(abs(real(r2dc)) > 10.0e0_knd .or. abs(aimag(r2dc)) >  &
              10.0e0_knd) then
          r2dc = r2dc / 10.0e0_knd
          ir2de = ir2de + 1
          end if
          if(abs(real(r2dc)) < 1.0e0_knd .and. abs(aimag(r2dc)) <  &
              1.0e0_knd) then
          r2dc = r2dc * 10.0e0_knd
          ir2de = ir2de - 1
          end if
        return
        end subroutine
!
!
        subroutine r2eta (l, m, c, x1, eta, nee, wm, limeta, maxd, maxlp, ndec, &
                          nex, maxn, maxp, minacc, lowtest, enr, sneuf, sneun, &
                          ineue, sneudf, sneudr, pdratt, pratb, pratt, &
                          pcoefn, ipcoefn, pdcoefn, ipdcoefn, nsubf, nsubmf, &
                          idigc, itestm, ichoicer1, naccr1, r1c, ir1e, r1dc, &
                          ir1de, naccmax, naccrpl, naccr, kcor, r2c, ir2e, &
                          r2dc, ir2de, nacceta, nacciop, jeta, iopnee, &
                          neemark, naccd, naccn, naccnmax, nacccor, naccns)
!
!  purpose:     to calculate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x, using an expansion of spherical Neumann
!               functions.
!
!  parameters:
!
!     input:   l        : l
!              m        : m
!              c        : complex c
!              x1       : x-1
!              eta      : value for eta used in calculation
!              nee      : index in the array of eta values in the main
!                         program that corresponds to the value of eta
!                         used in r2eta calculations
!              wm       : 1 - eta*eta = sin(nee)*sin(nee)
!              limeta   : maximum number of terms available in the sums
!                         for r2 and r2d
!              maxd     : dimension of enr array
!              maxlp    : maximum  l value desired; dimension
!                         of the sneudr arrays
!              ndec     : number of decimal digits for real(knd)
!              nex      : maximum exponent available for real(knd)
!              maxn     : dimension of sneuf, sneudf, sneun and ineue
!                         arrays
!              maxp     : dimension of pdratt, pratb, and pratt arrays
!              minacc   : minimum number of accurate decimal digits
!                         that are requested
!              lowtest  : minimum accuracy of calculated radial
!                         functions of the second kind for previous
!                         values of l
!              enr      : array of ratios of successive d coefficients
!              sneuf    : array of ratios of successive spherical
!                         Neumann functions of the same parity
!              sneun    : array of characteristics for Neumann functions
!              ineue    : array of exponents corresponding to sneun
!              sneudf   : array of ratios of successive first
!                         derivatives of spherical Neumann functions of
!                         the same parity
!              sneudr   : array of ratios of first derivatives of the
!                         spherical Neumann functions to the
!                         corresponding functions
!              pdratt   : array of ratios of successive first
!                         derivatives of the associated Legendre
!                         functions of the first kind of the same parity
!                         (used in numerator series)
!              pratb    : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in denominator series)
!              pratt    : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in numerator series)
!              pcoefn   : characteristic of the ratio of the numerator
!                         and denominator associated Legendre functions
!                         of the first kind of order m and degree l
!              ipcoefn  : exponent corresponding to pcoefn
!              pdcoefn  : characteristic of the ratio of the first
!                         derivative of the associated Legendre function
!                         of the first kind in the numerator and the
!                         associated Legendre function of the first kind
!                         in the denominator, both of order m and
!                         degree l
!              ipdcoefn : exponent corresponding to pdcoefn
!              nsubf    : number of decimal digits of subtraction error
!                         in computing the Flammer normalization factor
!              nsubmf   : number of decimal digits of subtraction error
!                         in computing the Morse-Feshbach normalization
!                         factor
!              idigc    : number of decimal digits of convergence of the
!                         eigenvalue using the Bouwkamp procedure
!              itestm   : agreement in decimal digits between the
!                         forward and backward recursions for the d
!                         coefficient ratios
!              ichoicer1: integer equal to 1 if the radial functions
!                         of the first kind are computed using the
!                         expansion with eta set equal to zero; equal
!                         to 2 if the traditional expansion with eta
!                         set equal to unity is used
!              naccr1   : estimate of number of decimal digits of
!                         accuracy of the radial function of the first
!                         kind (and its first derivative) calculated
!                         using the subtraction errors in their
!                         calculation
!              r1c      : first kind radial function (calculated in
!                         either r1besa or r1besb)
!              irie     : exponent corresponding to r1c
!              r1dc     : characteristic of the first derivative with
!                         respect to x of the radial function of the
!                         first kind (calculated in r1besa or r1besb)
!              ir1de    : exponent corresponding to r1dc
!              naccmax  : maximum accuracy (in decimal digits) obtained
!                         for the current value of l from previous
!                         r2eta calculations
!              naccrpl  : expected subtraction error in forming the
!                         Wronskian
!              naccr    : maximum accuracy obtained for the radial
!                         functions of the second kind for the current
!                         value of l using other methods as well as
!                         previous tries using r2eta
!
!     input/output: kcor: only applies if the denominator calculation
!                         is likely to suffer dec digits of subtraction
!                         error, but would have suffered ndec + kcor
!                         digits of subtraction if the precision were
!                         sufficiently large. Here the value for the
!                         denominator denom using ndec digits of
!                         precision will be larger that its correct
!                         value by a factor of about 10**kcor. Use of
!                         kcor provides a better estimate of the
!                         expected value for r2temp to aid in its
!                         calculation. The value of kcor obtained for
!                         l - 1 is used below for l. the value of kcor
!                         obtained below for this l is then available
!                         to be used for l + 1. Note that kcor is only
!                         calculated below if the Wronskian is used to
!                         obtain a useful value for the denominator when
!                         the numerator is sufficiently accurate so as
!                         to provide useful values for r2 and r2d.
!
!
!     output:  r2c      : characteristic of the prolate radial function
!                         of the second kind
!              ir2e     : exponent of the prolate radial function of the
!                         second kind
!              r2dc     : characteristic of the first derivative with
!                         respect to x of the prolate radial function
!                         of the second kind
!              ir2de    : exponent corresponding to r2dc
!              nacceta  : estimated number of accurate decimal digits in
!                         r2 and r2d. computed from the Wronskian if
!                         nacciop = 0; estimated from the series if
!                         nacciop = 1
!              nacciop  : integer flag = 1 if the denominator in the
!                         expressions for r2 and r2d is computed using
!                         the theoretical Wronskian and known values for
!                         r1 and r1d. nacciop = 0 otherwise
!              jeta     : maximum number of terms taken in the numerator
!                         sums for r2 and r2d
!              iopnee   : integer flag = 0 if none of the values used
!                         previously for eta in the variable eta method
!                         has led to estimated accuracy in the numerator
!                         of at least 3 digits (where iopnee is set to
!                         1) or has led to a subtraction error larger
!                         than either ndec-lowtest digits or ndec-naccr
!                         digits in either the r2 or the r2d numerator
!                         series (where iopnee is set to 2)
!              neemark  : index for the eta value in the array storing
!                         eta values in the main program that
!                         corresponds to the last eta value used for
!                         which iopnee = 0
!              naccd    : estimated accuracy of the denominator
!              naccn    : estimated accuracy of the numerator
!              naccnmax : maximum numerator accuracy (in decimal digits)
!                         obtained for the current value of l
!              nacccor  : minimum of naccrpl and the subtraction error
!                         during forming the Wronskian using r2 and r2d
!                         values computed using the variable eta method
!              naccns   : larger of the subtraction error in the
!                         calculation of the numerator for r2 and for
!                         r2d
!
!
        use param
!
!  real(knd) and complex(knd) scalars and arrays
        real(knd) ca, dcon, dconi, dec, eta, etas, pcoefn, pdcoefn, &
                  rm, rm2, r2dcoef1, r2est, r2test, sumdenpi, sumdenpr, &
                  sumdnp1r, sumdnp1i, sumdnp2r, sumdnp2i, sumnpi, sumnpr, &
                  ten, test, testd1, testd2, teste, testeo, testdm1, testdm2, &
                  testm, testr2, txi, txr, txd1i, txd2i, txd1r, txd2r, wm, xet, &
                  xets, x1
        complex(knd) c, denom, dnew, dnewd1, dnewd2, dnewsum, dnewdsum1, &
                     dnewdsum2, dold, doldd1, doldd2, r1c, r1dc, r2c, &
                     r2dc, r2dcoef2, r2dtemp1, r2dtemp2, r2dtemp, r2temp, &
                     reld12, sr2dtemp1, sr2dtemp2, sr2temp, sumcoef, wronc, &
                     wronca, wroncb, wront
        real(knd) pratb(maxp), pratt(maxp), pdratt(maxp)
        complex(knd) enr(maxd), sneudr(maxlp), sneun(maxn), sneuf(maxn), &
                     sneudf(maxn)
!
!  integer arrays
        integer ineue(maxn)
!
        ca = abs(c)
        wront = 1.0e0_knd / (c * x1 * (x1 + 2.0e0_knd))
        ten = 10.0e0_knd
        dec = ten ** (-ndec - 2)
        nfac = nex / 3
        if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
        teste = ten ** (nfac)
        testeo = 1.0e0_knd / teste
        iscale = 0
        rm = real(m, knd)
        dcon = ten ** (-ndec - 2)
        dconi = ten ** (ndec + 2)
        etas = eta * eta
        xet = sqrt(x1 * (x1 + 2.0e0_knd) + etas)
        ltest = 2 * int((abs(real(c)) + abs(aimag(c))) * xet)
        jtest = (ltest - m) / 2
        xets = xet * xet
        sumcoef = (ten ** (-ir1de - ineue(l + 1) - ipcoefn + naccrpl))/ &
                (c * x1 * (x1 + 2.0e0_knd) * r1dc * sneun(l + 1) * pcoefn)
        r2dcoef1 = -eta * wm / (xets * xet)
        r2dcoef2 = c * (x1 + 1.0e0_knd) / xet
        reld12 = (r2dcoef2 / r2dcoef1) * sneudr(l + 1) * (pcoefn/ &
               pdcoefn) * (ten ** (ipcoefn - ipdcoefn))
        rm2 = rm * 2.0e0_knd
        lm2 = (l - m) / 2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix = l - m - 2 * lm2
        lim = limeta / 2 - ix
!
!  compute radial function of the second kind and its first derivative
!
!  backward series for denominator
        denom = (1.0e0_knd, 0.0e0_knd)
        sumdenpr = 1.0e0_knd
        sumdenpi = 0.0e0_knd
        if (lm2 < 1) go to 20
        dold = (1.0e0_knd, 0.0e0_knd)
          do 10 j = lm2, 1, -1
          jj = j + j + ix
          dnew = dold / (pratb(jj + 1) * enr(j))
          denom = denom + dnew
          if(real(dnew) > 0.0e0_knd) sumdenpr = sumdenpr + real(dnew)
          if(aimag(dnew) > 0.0e0_knd) sumdenpi = sumdenpi + aimag(dnew)
          if(abs(dnew / denom) < dec) go to 20
          dold = dnew
10        continue
20      continue
!
!  forward series for denominator
        dold = (1.0e0_knd, 0.0e0_knd)
          do 30 j = lm2 + 1, lim
          jj = j + j + ix
          dnew = dold * enr(j) * pratb(jj + 1)
          denom = denom + dnew
          if(real(dnew) > 0.0e0_knd) sumdenpr = sumdenpr + real(dnew)
          if(aimag(dnew) > 0.0e0_knd) sumdenpi = sumdenpi + aimag(dnew)
          if(abs(dnew / denom) < dec) go to 40
          dold = dnew
30        continue
40      continue
        jden = j
        nsubdenr = -int(log10(abs(real(denom) / sumdenpr) + dec))
        nsubdeni = 0
        if(sumdenpi /= 0.0e0_knd) nsubdeni = -int(log10(abs(aimag(denom)/ &
                                     sumdenpi) + dec))
        if(nsubdenr < 0) nsubdenr = 0
        if(nsubdeni < 0) nsubdeni = 0
        if(nsubdenr > ndec) nsubdenr = ndec
        if(nsubdeni > ndec) nsubdeni = ndec
        nsubd = nsubdenr
          if(aimag(denom) /= 0.0e0_knd) then
          iterm = int(log10(abs(aimag(denom) / real(denom))))
          if(iterm < 0) nsubd = max(nsubdenr, nsubdeni + iterm)
          if(iterm > 0) nsubd = max(nsubdeni, nsubdenr - iterm)
          end if
         naccd = min(ndec - 2 - nsubd, idigc - 1 - nsubd)
         if(naccd < 0) naccd = 0
        r2est = abs(sumcoef * denom)
        r2test = r2est * dconi
          if(nsubd < ndec) then
          testr2 = r2est * 0.01e0_knd
          else
          testr2 = r2est * (10.0e0_knd ** (-kcor - 2))
          end if
!
!  backward series for numerator
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd1 = (1.0e0_knd, 0.0e0_knd)
        doldd2 = reld12
        r2temp = (1.0e0_knd, 0.0e0_knd)
        sumnpr = 1.0e0_knd
        sumnpi = 0.0e0_knd
        sumdnp1r = 0.0e0_knd
        sumdnp1i = 0.0e0_knd
        sumdnp2r = 0.0e0_knd
        sumdnp2i = 0.0e0_knd
        r2dtemp1 = (0.0e0_knd, 0.0e0_knd)
        r2dtemp2 = doldd2
        if(real(doldd2) > 0.0e0_knd) sumdnp2r = real(doldd2)
        if(aimag(doldd2) > 0.0e0_knd) sumdnp2i = aimag(doldd2)
        if(l /= 0) r2dtemp1 = (1.0e0_knd, 0.0e0_knd)
        if(l /= 0) sumdnp1r = 1.0e0_knd
        if(lm2 == 0) go to 60
          do 50 j = lm2, 1, -1
          jj = j + j + ix
          dnew = -dold / (sneuf(jj + m) * pratt(jj + 1) * enr(j))
          dnewd1 = -doldd1 / (sneuf(jj + m) * pdratt(jj + 1) * enr(j))
          dnewd2 = -doldd2 / (sneudf(jj + m) * pratt(jj + 1) * enr(j))
          r2temp = r2temp + dnew
          r2dtemp1 = r2dtemp1 + dnewd1
          r2dtemp2 = r2dtemp2 + dnewd2
          if(real(dnew) > 0.0e0_knd) sumnpr = sumnpr + real(dnew)
          if(aimag(dnew) > 0.0e0_knd) sumnpi = sumnpi + aimag(dnew)
          if(real(dnewd1) > 0.0e0_knd) sumdnp1r = sumdnp1r + real(dnewd1)
          if(aimag(dnewd1) > 0.0e0_knd) sumdnp1i = sumdnp1i + aimag(dnewd1)
          if(real(dnewd2) > 0.0e0_knd) sumdnp2r = sumdnp2r + real(dnewd2)
          if(aimag(dnewd2) > 0.0e0_knd) sumdnp2i = sumdnp2i + aimag(dnewd2)
          if(abs(dnew / r2temp) + abs(dnewd1 / r2dtemp1) + abs(dnewd2 / r2dtemp2) &
                < dcon) go to 60
          dold = dnew
          doldd1 = dnewd1
          doldd2 = dnewd2
50        continue
60      continue
        if(m == 0 .and. jj == 2) then
        r2dtemp1 = r2dtemp1 - dnewd1
        if(real(dnewd1) > 0.0e0_knd) sumdnp1r = sumdnp1r - real(dnewd1)
        if(aimag(dnewd1) > 0.0e0_knd) sumdnp1i = sumdnp1i - aimag(dnewd1)
        end if
!
!  forward series for numerator
        dold = (1.0e0_knd, 0.0e0_knd)
        doldd1 = (1.0e0_knd, 0.0e0_knd)
        doldd2 = reld12
        test = 1.0e0_knd
        testd1 = 1.0e0_knd
        testd2 = 1.0e0_knd
        testm = 1.0e0_knd
        testdm1 = 1.0e0_knd
        testdm2 = 1.0e0_knd
        js = lm2
        jds1 = lm2
        jds2 = lm2
        sr2temp = r2temp
        sr2dtemp1 = r2dtemp1
        sr2dtemp2 = r2dtemp2
        txr = sumnpr
        txi = sumnpi
        txd1r = sumdnp1r
        txd1i = sumdnp1i
        txd2r = sumdnp2r
        txd2i = sumdnp2i
        dnewsum = (0.0e0_knd, 0.0e0_knd)
        dnewdsum1 = (0.0e0_knd, 0.0e0_knd)
        dnewdsum2 = (0.0e0_knd, 0.0e0_knd)
        kount = 0
        kountd1 = 0
        kountd2 = 0
        lflag = 0
          do 110 j = lm2 + 1, lim - 1
          jj = j + j + ix
          if(j >= jtest .and. lflag == 0) lflag = 1
          kount = kount + 1
          kountd1 = kountd1 + 1
          kountd2 = kountd2 + 1
          dnew = -dold * enr(j) * sneuf(jj + m) * pratt(jj + 1)
          dnewd1 = -doldd1 * enr(j) * sneuf(jj + m) * pdratt(jj + 1)
          dnewd2 = -doldd2 * enr(j) * sneudf(jj + m) * pratt(jj + 1)
          if(((real(dnew) / real(dold)) <= 0.0e0_knd) .or. (kount == 100) &
                .or. j >= lim - 10 .or. lflag == 2) go to 70
          dnewsum = dnewsum + dnew
          go to 80
70        r2temp = r2temp + dnewsum
          kount = 0
          if(abs(r2temp) > r2test) go to 120
          if(real(dnewsum) > 0.0e0_knd) sumnpr = sumnpr + real(dnewsum)
          if(aimag(dnewsum) > 0.0e0_knd) sumnpi = sumnpi + aimag(dnewsum)
            if(lflag == 0) then
            if(max(abs(r2temp), sumnpr, sumnpi) > testr2) lflag = 1
            if(j >= lim - 10) lflag = 1
            end if
          dnewsum = dnew
          if(abs(dnewsum) /= 0.0e0_knd) test = abs(dnewsum / r2temp)
          if(test >= testm .or. lflag == 0) go to 80
          testm = test
          sr2temp = r2temp
          js = j
          txr = sumnpr
          txi = sumnpi
80        if(((real(dnewd1) / real(doldd1)) <= 0.0e0_knd) .or.  &
                (kountd1 == 100) .or. j >= lim - 10 .or. lflag == 2) go to 85
          dnewdsum1 = dnewdsum1 + dnewd1
          go to 90
85        r2dtemp1 = r2dtemp1 + dnewdsum1
          kountd1 = 0
          if(real(dnewdsum1) > 0.0e0_knd) sumdnp1r = sumdnp1r + &
                real(dnewdsum1)
          if(aimag(dnewdsum1) > 0.0e0_knd) sumdnp1i = sumdnp1i + &
                                                   aimag(dnewdsum1)
          dnewdsum1 = dnewd1
          if(abs(dnewdsum1) /= 0.0e0_knd) testd1 = abs(dnewdsum1 / r2dtemp1)
          if(testd1 >= testdm1 .or. lflag == 0) go to 90
          testdm1 = testd1
          sr2dtemp1 = r2dtemp1
          jds1 = j
          txd1r = sumdnp1r
          txd1i = sumdnp1i
90        if((real(dnewd2) / real(doldd2) <= 0.0e0_knd) .or.  &
             (kountd2 == 100) .or. j >= lim - 10 .or. lflag == 2) go to 95
          dnewdsum2 = dnewdsum2 + dnewd2
          go to 100
95        r2dtemp2 = r2dtemp2 + dnewdsum2
          kountd2 = 0
          if(real(dnewdsum2) > 0.0e0_knd) sumdnp2r = sumdnp2r + &
                real(dnewdsum2)
          if(aimag(dnewdsum2) > 0.0e0_knd) sumdnp2i = sumdnp2i + &
                aimag(dnewdsum2)
          if(abs(dnewdsum2) /= 0.0e0_knd) testd2 = abs(dnewdsum2 / r2dtemp2)
          dnewdsum2 = dnewd2
          if(testd2 >= testdm2 .or. lflag == 0) go to 100
          testdm2 = testd2
          sr2dtemp2 = r2dtemp2
          jds2 = j
          txd2r = sumdnp2r
          txd2i = sumdnp2i
100       if(test + testd1 + testd2 < dcon .and. lflag == 2) go to 130
          if(test + testd1 + testd2 < dcon .and. lflag == 1) lflag = 2
            if(abs(r2temp) > teste) then
            r2temp = r2temp * testeo
            r2dtemp1 = r2dtemp1 * testeo
            r2dtemp2 = r2dtemp2 * testeo
            sr2temp = sr2temp * testeo
            sr2dtemp1 = sr2dtemp1 * testeo
            sr2dtemp2 = sr2dtemp2 * testeo
            sumnpr = sumnpr * testeo
            sumdnp1r = sumdnp1r * testeo
            sumdnp2r = sumdnp2r * testeo
            sumnpi = sumnpi * testeo
            sumdnp1i = sumdnp1i * testeo
            sumdnp2i = sumdnp2i * testeo
            dnew = dnew * testeo
            dnewd1 = dnewd1 * testeo
            dnewd2 = dnewd2 * testeo
            dnewsum = dnewsum * testeo
            dnewdsum1 = dnewdsum1 * testeo
            dnewdsum2 = dnewdsum2 * testeo
            iscale = iscale + nfac
            r2test = r2test * testeo
            testr2 = testr2 * testeo
            txr = txr * testeo
            txi = txi * testeo
            txd1r = txd1r * testeo
            txd1i = txd1i * testeo
            txd2r = txd2r * testeo
            txd2i = txd2i * testeo
            end if
          dold = dnew
          doldd1 = dnewd1
          doldd2 = dnewd2
110       continue
120     r2temp = sr2temp
        r2dtemp1 = sr2dtemp1
        r2dtemp2 = sr2dtemp2
        go to 135
130     txr = sumnpr
        txi = sumnpi
        txd1r = sumdnp1r
        txd1i = sumdnp1i
        txd2r = sumdnp2r
        txd2i = sumdnp2i
        js = j
        jds1 = j
        jds2 = j
135     continue
        jmax = j
        jds = max(jds1, jds2)
        jeta = max(js, jds, jden)
        jtestm = ndec
        if(testm /= 0.0e0_knd) jtestm = -int(log10(abs(testm)))
        if(jtestm < 0) jtestm = 0
        if(jtestm > ndec) jtestm = ndec
        jtestdm1 = ndec
          if(testdm1 /= 0.0e0_knd) then
          jtestdm1 = -int(log10(abs(testdm1)))
          end if
        if(jtestdm1 < 0) jtestdm1 = 0
        if(jtestdm1 > ndec) jtestdm1 = ndec
        jtestdm2 = ndec
          if(testdm2 /= 0.0e0_knd) then
          jtestdm2 = -int(log10(abs(testdm2)))
          end if
        if(jtestdm2 < 0) jtestdm2 = 0
        if(jtestdm2 > ndec) jtestdm2 = ndec
        naccns1r = -int(log10(abs(real(r2temp) / txr) + dec))
        naccns1i = 0
        if(txi /= 0.0e0_knd) naccns1i = -int(log10(abs(aimag(r2temp) &
                                  /txi) + dec))
        if(naccns1r < 0) naccns1r = 0
        if(naccns1r > ndec) naccns1r = ndec
        if(naccns1i < 0) naccns1i = 0
        if(naccns1i > ndec) naccns1i = ndec
        r2dtemp = r2dtemp1 + r2dtemp2
        naccns2r = -int(log10(abs(real(r2dtemp) / (txd1r + txd2r) + dec)))
        naccns2i = 0
        if(txd1i + txd2i /= 0.0e0_knd) naccns2i= &
                   -int(log10(abs(aimag(r2dtemp) / (txd1i + txd2i) + dec)))
        if(naccns2r < 0) naccns2r = 0
        if(naccns2r > ndec) naccns2r = ndec
        if(naccns2i < 0) naccns2i = 0
        if(naccns2i > ndec) naccns2i = ndec
        naccns1 = max(naccns1r, naccns1i)
        naccns2 = max(naccns2r, naccns2i)
          if(abs(r2dtemp1 * r2dtemp2) /= 0.0e0_knd) then
          icord = int(log10(abs(r2dtemp1 / r2dtemp2)))
          if(icord > 0) jtestdm2 = jtestdm2 + icord
          if(icord < 0) jtestdm1 = jtestdm1 - icord
          jtestdm = min(jtestdm1, jtestdm2)
          end if
        if(abs(r2dtemp1) == 0.0e0_knd) jtestdm = jtestdm2
        if(abs(r2dtemp2) == 0.0e0_knd) jtestdm = jtestdm1
        ncorr = max(0,-int(log10(x1) - 0.001e0_knd)) + 1
        if(x1 >= 0.05e0_knd) ncorr = 1
          if(x1 <= 0.02e0_knd) then
          naccn1 = min(jtestm - 1, ndec - ncorr, ndec - naccns1 - 1, idigc - 1, &
                     itestm - 1)
          naccn2 = min(jtestdm - 1, ndec - ncorr, ndec - naccns2 - 1, idigc - 1, &
                     itestm - 1)
          else
          naccn1 = min(jtestm - 1, ndec - ncorr, ndec - naccns1 - 1, idigc - 1, &
                     itestm - 1)
          naccn2 = min(jtestdm - 1, ndec - ncorr, ndec - naccns2 - 1, idigc - 1, &
                     itestm - 1)
          end if
        naccn = min(naccn1, naccn2, ndec - 1 - max(0, int(log10(abs(c)))))
        naccns = max(naccns1, naccns2)
        if(naccn > ndec - 2) naccn = ndec - 2
        if(naccn < 0) naccn = 0
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2c = r2temp * sneun(l + 1) * pcoefn / denom
        iterm = 0
        if(r2c /= (0.0e0_knd, 0.0e0_knd)) iterm = int(log10(abs(r2c)))
        ir2e = ineue(l + 1) + ipcoefn + iscale + iterm
        r2c = r2c * 10.0e0_knd ** (-iterm)
        r2dc = (r2dcoef1 * r2dtemp * sneun(l + 1) * pdcoefn / denom)* &
             10.0e0_knd ** (ineue(l + 1) + ipdcoefn + iscale-ir2e)
        iterm = 0
        if(r2dc /= (0.0e0_knd, 0.0e0_knd)) iterm = int(log10(abs(r2dc)))
        ir2de = ir2e + iterm
        r2dc = r2dc * 10.0e0_knd ** (-iterm)
if (debug) then
        write(40, 140) jmax, jden, lim, js, jtestm, naccns1, jds, jtestdm, &
                      naccns2, naccn, naccd
140     format(8x,'r2eta: numerator, denominator converged in ', &
               i6,' ,',i6,' terms; ',i6,' terms available.',/, &
               15x,'best r2 at ',i6,' terms with convergence to',i3, &
               ' digits;',i3,' digits subtr. error.',/,15x, &
               'best r2d at ',i6,' terms with convergence to',i3, &
               ' digits;',i3,' digits subtr. error.',/,15x, &
               'estimated numerator and denominator accuracy is ',i4, &
               ' and',i4,' digits.')
end if
        wronca = r1c * r2dc * ten ** (ir1e + ir2de)
        wroncb = r2c * r1dc * ten ** (ir2e + ir1de)
        wronc = wronca - wroncb
        nacceta = -int(log10(abs((wronc - wront) / wront) + dec))
        if(nacceta < 0) nacceta = 0
        if(nacceta > ndec-1) nacceta = ndec-1
        naccetao = nacceta
        nacccor = -int(log10(abs((wronca - wroncb) / wronca) + dec))
        if(nacccor < 0) nacccor = 0
        nacccors = nacccor
        nacccor = min(nacccor, naccrpl + 1)
        if(nacceta > 0) nacceta = min(nacceta + nacccor, naccr1)
          if(nacceta == 0 .and. naccrpl > 0) then
          imat = -int(log10(abs((r1c * (0.0e0_knd, 1.0e0_knd) - r2c* &
          ten ** (ir2e - ir1e)) / r1c) + dec))
          imatd = -int(log10(abs((r1dc * (0.0e0_knd, 1.0e0_knd) - r2dc* &
          ten ** (ir2de - ir1de)) / r1dc) + dec))
          nacceta = min(imat, imatd, naccrpl, naccr1)
          if(nacceta < 0) nacceta = 0
          end if
145     continue
        nacciop = 0
        if(nacceta < minacc .and. naccn - nacccors > naccd .and.  &
            naccn > nacceta .and. (naccn >= min(naccr, 3) .or.  &
            (jtestm >= ndec - 1 .and. jtestdm >= ndec - 1))) nacciop = 1
        if(jtestm < 6 .and. jtestdm < 6) nacciop = 0
        if(naccn - nacccors <= 2) nacciop = 0
        if(nacciop /= 1) go to 160
if (debug) then
        write(40, 150)
150     format(15x,'denominator calculated using Wronskian.')
end if
        nacceta = naccn - nacccors
        r2c = r2c * wront / wronc
        kcor = int(log10(abs(wront / wronc)))
        iterm = 0
        if(r2c /= (0.0e0_knd, 0.0e0_knd)) iterm = int(log10(abs(r2c)))
        ir2e = ir2e + iterm
        r2c = r2c * 10.0e0_knd ** (-iterm)
        r2dc = r2dc * wront / wronc
        iterm = 0
        if(r2dc /= (0.0e0_knd, 0.0e0_knd)) iterm = int(log10(abs(r2dc)))
        ir2de = ir2de + iterm
        r2dc = r2dc * 10.0e0_knd ** (-iterm)
160     continue
          if(abs(real(r2c)) > 10.0e0_knd .or. abs(aimag(r2c)) >  &
              10.0e0_knd) then
          r2c = r2c / 10.0e0_knd
          ir2e = ir2e + 1
          end if
          if(abs(real(r2c)) < 1.0e0_knd .and. abs(aimag(r2c)) <  &
              1.0e0_knd) then
          r2c = r2c * 10.0e0_knd
          ir2e = ir2e - 1
          end if
          if(abs(real(r2dc)) > 10.0e0_knd .or. abs(aimag(r2dc)) >  &
              10.0e0_knd) then
          r2dc = r2dc / 10.0e0_knd
          ir2de = ir2de + 1
          end if
          if(abs(real(r2dc)) < 1.0e0_knd .and. abs(aimag(r2dc)) <  &
              1.0e0_knd) then
          r2dc = r2dc * 10.0e0_knd
          ir2de = ir2de - 1
          end if
        naccnmaxp = naccnmax
        if(naccnmax > naccn .and. naccnmax >= 3) iopnee = 1
        if(naccn > naccnmax) naccnmax = naccn
        if(ndec - naccns < naccr + 2) iopnee = 2
        if(nacciop == 0 .and. naccn > nacceta) naccnmax = max(naccnmaxp, &
            nacceta)
          if(nacceta < 3 .and. naccmax == 0 .and. ndec - naccns >=  &
              naccr) then
          iopnee = 0
          neemark = nee
          end if
        if(ndec - naccns >= naccr + 2) iopnee = 0
        if(jtestm == ndec .and. jtestdm == ndec .and. naccn <= max(naccr, 5)) &
            iopnee = 2
        if(naccns1 == ndec .or. naccns2 == ndec) iopnee = 2
        nacccor = min(nacccor, nacceta - naccetao, nacceta)
        if(nacccor < 0) nacccor = 0
        return
        end subroutine
!
!
        subroutine geteig (m, c, nbp, ndec, lnum, maxe, eigst)
!
!  purpose:     to determine estimates of the eigenvalues for
!               orders up to and beyond the break point for a given m
!               and complex c. These estimates will be starting values
!               for the Bouwkamp procedure.
!
!  parameters:
!
!     input:    m     : m
!               c     : complex c
!               nbp   : breakpoint value
!               ndec  : number of decimal digits for real(knd)
!               lnum  : number of values of l for which functions
!                       are desired
!               maxe  : dimension of matrix elements
!     output:   eigst : array of eigenvalue estimates
!
        use param
!
!  real(knd) scalars and complex(knd) scalars and arrays
        complex(knd) c, c2, c4, d(maxe), e(maxe), f(maxe), eigst(lnum), g(maxe)
        real(knd) cqr, cqi, c2m, xl, cm, em
!
        m2 = m + m
        em = m
        c2 = c * c
        c4 = c2 * c2
        cqi = abs(aimag(c))
        cqr = abs(real(c))
        cm = abs(c)
        c2m = abs(c2)
!  obtain starting eigenvalues eigst(i) for the Bouwkamp procedure
!
            lime = max(100, nbp + nbp)
            imax = max(50, nbp) + 5
              do 10 i = 1, lime
              xl = real(m + i + i - 2, knd)
              d(i) = xl * (xl + 1.0e0_knd) / c2 + (2.0e0_knd * xl * (xl + 1.0e0_knd)- &
                   2.0e0_knd * em * em - 1.0e0_knd) / ((2.0e0_knd * xl - 1.0e0_knd)* &
                  (2.0e0_knd * xl + 3.0e0_knd))
10            continue
            nm1 = lime - 1
              do 20 i = 1, nm1
              xl = real(m + i + i - 2, knd)
              e(i) = (1.0e0_knd / (2.0e0_knd * xl + 3.0e0_knd))* &
                   sqrt(((xl + 2.0e0_knd + em) * (xl + 1.0e0_knd + em) * (xl+ &
                   2.0e0_knd - em) * (xl + (1.0e0_knd) - em)) / ((2.0e0_knd * xl+ &
                   5.0e0_knd) * (2.0e0_knd * xl + 1.0e0_knd)))
20            continue
            call cmtql1(lime, maxe, ndec, d, e)
              do 30 i = 1, lime
              f(i) = d(i) * c2
30            continue
              do 40 i = 1, lime
              xl = real(m + i + i - 1, knd)
              d(i) = xl * (xl + 1.0e0_knd) / c2 + (2.0e0_knd * xl * (xl + 1.0e0_knd)- &
                   2.0e0_knd * em * em - 1.0e0_knd) / ((2.0e0_knd * xl - 1.0e0_knd)* &
                   (2.0e0_knd * xl + 3.0e0_knd))
40            continue
            nm1 = lime - 1
              do 50 i = 1, nm1
              xl = real(m + i + i - 1, knd)
              e(i) = (1.0e0_knd / (2.0e0_knd * xl + 3.0e0_knd))* &
                   sqrt(((xl + 2.0e0_knd + em) * (xl + 1.0e0_knd + em) * (xl+ &
                   2.0e0_knd - em) * (xl + (1.0e0_knd) - em)) / ((2.0e0_knd * xl+ &
                   5.0e0_knd) * (2.0e0_knd * xl + 1.0e0_knd)))
50            continue
            call cmtql1(lime, maxe, ndec, d, e)
              do 60 i = 1, lime
              g(i) = d(i) * c2
60            continue
            call eigorder(c, m, lime, maxe, nbp, f, g)
            lnumo2 = lnum / 2
            iupp = min(imax, lnumo2)
            if(iupp == 0) go to 80
              do 70 i = 1, iupp
              eigst(i + i - 1) = f(i)
              eigst(i + i) = g(i)
70            continue
            lnumo2 = lnum / 2
80          if(2 * lnumo2 /= lnum .and. imax > lnumo2) &
                 eigst(lnum) = f(lnumo2 + 1)
            return
            end subroutine
!
!
        subroutine cmtql1(n, maxe, ndec, d, e)
!
!  purpose:     to compute the eigenvalues of a complex symmetric
!               tridiagonal matrix
!
!  subroutine published by Jane Cullum and Ralph Willoughby in 'Lanczos
!  algorithsm for Large symmetric eigenvalue computations,' 2002;
!  I have removed the parameter ierr from the call statement and set
!  it initially equal to zero. I also relabeled the subdiagonal elements
!  e prior to input rather than after input.
!
!  parameters:
!
!     input:    n:    order of the matrix
!               maxe: dimension of both vectors d and e
!               ndec: number of decimal digits for real(knd)
!               d:    diagonal elements of the tridiagonal matrix
!               e:    subdiagonal elements of the matrix,
!                     beginning with e(1); e(n) is set equal to zero
!
!     output:   d  : eigenvalues
!
        use param
!
        real(knd) machep, eps, temp, t0, t1, zero, half, one, two
        complex(knd) d(maxe), e(maxe), b, c, f, g, p, r, s, w, czero, cone
!
        machep = 10.0e0_knd ** (-ndec)
        eps = 100.0e0_knd * machep
        zero = 0.0e0_knd
        half = 0.5e0_knd
        one = 1.0e0_knd
        two = 2.0e0_knd
        czero = cmplx(zero, zero, knd)
        cone = cmplx(one, zero, knd)
        ierr = 0
        e(n) = czero
          do 140 l = 1, n
          j = 0
20          do 30 m = l, n
            if(m == n) go to 40
            temp = abs(d(m)) + abs(d(m + 1))
            if(abs(e(m)) <= temp * machep) go to 40
30          continue
40        p = d(l)
          if(m == l) go to 100
          if(j == 100) go to 150
          j = j + 1
          g = (d(l + 1) - p) * half
          t0 = abs(g)
          t1 = abs(e(l))
          if(t0 > t1) go to 50
          w = g / e(l)
          r = sqrt(cone + w * w)
          t0 = abs(w + r)
          t1 = abs(w - r)
          temp = one
          if(t1 > t0) temp = -one
          g = d(m) - p + e(l) / (w + temp * r)
          go to 60
50        continue
          w = e(l) / g
          r = sqrt(cone + w * w)
          t0 = sqrt(cone + r)
          t1 = sqrt(cone - r)
          temp = one
          if(t1 > t0) temp = -one
          g = d(m) - p + w * e(l) / (cone + temp * r)
60        continue
          s = cone
          c = -cone
          p = czero
          mml = m - l
            do 90 i1 = 1, mml
            i = m - i1
            f = s * e(i)
            b = -c * e(i)
            t0 = abs(g)
            t1 = abs(f)
            if(t1 > t0) go to 70
            w = f / g
            r = sqrt(cone + w * w)
            e(i + 1) = g * r
            c = cone / r
            s = w * c
            go to 80
70          continue
            w = g / f
            r = sqrt(cone + w * w)
            e(i + 1) = f * r
            s = cone / r
            c = w * s
80          continue
            temp = abs(w) ** 2 + one
            t0 = sqrt(temp)
            t1 = abs(r)
            ierr = -l
            if(t1 <= eps * t0) go to 160
            ierr = 0
            g = d(i + 1) - p
            r = (d(i) - g) * s + two * c * b
            p = s * r
            d(i + 1) = g + p
            g = b - c * r
90          continue
          d(l) = d(l) - p
          e(l) = g
          e(m) = zero
          go to 20
100       if(l == 1) go to 120
            do 110 i1 = 2, l
            i = l + 2 - i1
            if(abs(p) >= abs(d(i - 1))) go to 130
            d(i) = d(i - 1)
110         continue
120       i = 1
130       d(i) = p
140       continue
        go to 160
150     ierr = l
160     return
        end subroutine
!
!
        subroutine eigorder(c, m, n, np, nbp, f, g)
!
!
!  purpose:     To order the eigenvalues estimates so that the
!               resulting radial function magnitudes have a similar
!               behavior with increasing order as the corresponding
!               spherical Bessel functions
!  parameters:
!
!     input:    c     : complex c
!               m     : m
!               n     : number of even eigenvalue estimates obtained
!                       from the tridiagonal matrix; also the number
!                       of odd eigenvalue estimates
!               np    : dimension of both the array of even estimates
!                       and the array of odd estimates
!               nbp   : breakpoint value
!
!     output:   f     : array of ordered even eigenvalue estimates
!               g     : array of ordered odd eigenvalue estimates
!
        use param
!
!  real(knd) and complex(knd) scalars and arrays
        complex(knd) c, c2, f(np), g(np), fm(np), gm(np), fs(np), gs(np), p, &
                     testp
        real(knd) q, testpr, testpi
!
        imax = 1.2 * max(25, nbp / 2) + 5
        c2 = c * c
!
!  order even eigenvalues in ascending real part
         do i = 1, n - 1
         k = i
         p = f(i)
           do j = i + 1, n
           if(real(f(j)) < real(p)) then
           k = j
           p = f(j)
           end if
           end do
         if(k /= i) then
         f(k) = f(i)
         f(i) = p
         end if
         end do
!
!  order odd eigenvalues in ascending real part
         do i = 1, n - 1
         k = i
         p = g(i)
           do j = i + 1, n
             if(real(g(j)) < real(p)) then
             k = j
             p = g(j)
             end if
           end do
           if(k /= i) then
           g(k) = g(i)
           g(i) = p
           end if
         end do
!
!  locate oblate-like eigenvalues
         limp = 0
         nmax = int(abs(aimag(c))) / 2
         if(nmax == 0) go to 50
         limp = 0
           do j = 1, nmax
           q = real(m + j + j - 1, knd)
           testp = c2 - 2 * q * c * (0.0e0_knd, 1.0e0_knd) - 0.5e0_knd* &
                 (q * q - m * m + 1) - (0.0e0_knd, 1.0e0_knd) * q* &
                 (q * q - m * m + 1) / (8.0e0_knd * c) + (5.0e0_knd * q * q * q * q+ &
                 10.0e0_knd * q * q + 1 - 2 * m * m * (3.0e0_knd * q * q + 1) + m * m * m * m)/ &
                 (64.0_knd * c2)
           testpr = real(testp)
           testpi = aimag(testp)
             do i = 1, imax
             if(abs((real(f(i)) - testpr) / real(f(i))) > 0.001e0_knd) &
                go to 10
             if(abs((aimag(f(i)) - testpi) / aimag(f(i))) > 0.001e0_knd) &
                go to 10
             limp = limp + 1
             fs(limp) = f(i)
             isav = i
               do kk = i, imax - limp
               f(kk) = f(kk + 1)
               end do
             go to 20
10           end do
           go to 50
20           il = max(isav - 1, 1)
             do i = il, isav
             if(abs((real(g(i)) - testpr) / real(g(i))) > 0.001e0_knd) &
                 go to 30
             if(abs((aimag(g(i)) - testpi) / aimag(g(i))) > 0.001e0_knd) &
                 go to 30
             gs(limp) = g(i)
               do kk = i, imax - limp
               g(kk) = g(kk + 1)
               end do
             go to 40
30           continue
             end do
           gs(limp) = g(isav)
           do kk = isav, imax - limp
           g(kk) = g(kk + 1)
           end do
           go to 50
40         continue
           end do
50       continue
!
!  order oblate-like even eigenvalues for insertion in eigenvalue
!  sequence
         if(limp == 0) go to 160
         limpc = limp / 2 + 1
         limpr = (limp + 1) / 2
           do i = 1, limpr
           fm(limpc + i - 1) = fs(i + i - 1)
           end do
         if(limp == 1) go to 60
           do i = 1, limp - limpr
           fm(limpc - i) = fs(i + i)
           end do
60       continue
!
!  order oblate-like odd eigenvalues for insertion in eigenvalue
!  sequence
         limpc = (limp + 1) / 2
         limpr = limp / 2
           do i = 1, limp - limpr
           gm(limpc - i + 1) = gs(i + i - 1)
           end do
           if(limp == 1) go to 70
           do i = 1, limpr
           gm(limpc + i) = gs(i + i)
           end do
70       continue
!
!  insert paired eigenvalues in proper location in sequence
         iflag = 1
         if(limp /= (2 * (limp / 2))) iflag = -iflag
           if(iflag == -1) then
             do j = 1, imax - limp - 1
             if(real(fm(1)) - real(f(j)) < 0.0e0_knd) exit
             end do
           jlowf = j
           jlowg = j
           end if
           if(iflag == 1) then
             do j = 1, imax - limp - 1
             if(real(gm(1)) - real(g(j)) < 0.0e0_knd) exit
             end do
           jlowg = j
           jlowf = j + 1
           end if
           do i =1, imax - jlowf - limp + 1
           f(imax - i + 1) = f(imax - limp - i + 1)
           end do
           do i =1, imax - jlowg - limp + 1
           g(imax - i + 1) = g(imax - limp - i + 1)
           end do
         k = jlowf - 1
           do i = 1, limp
           k = k + 1
           if(iflag == 1) f(k) = fm(i)
           if(iflag == -1) f(k) = fm(limp - i + 1)
           end do
         k = jlowg - 1
           do i = 1, limp
           k = k + 1
           if(iflag == -1) g(k) = gm(limp - i + 1)
           if(iflag == 1) g(k) = gm(i)
           end do
160      continue
         return
         end subroutine
!
!
        subroutine conver (l, m, c, limd, maxd, blist, glist, blist1, glist1, &
                           ioprad, ienr, kindd, kindq, ndec, eigval, enr, &
                           idigc, itestm, kflag)
!
!  purpose:     to determine the eigenvalue using the Boukwamp method.
!  parameters:
!
!     input:    l     : l
!               m     : m
!               c     : complex c
!               limd  : number of enr values computed
!               maxd  : dimension of enr,blist,glist arrays
!               blist : knd coefficients used in computing d coefficients
!               glist : knd coefficients used in computing d coefficients
!               blist1: knd1 coefficients used in computine d coefficients
!               glist1: knd1 coefficients used in computing d coefficients
!               eigval: estimated value of the eigenvalue
!               ioprad: integer equal to 0 if radial functions are not
!                       desired; equal to 1 if radial functions of the
!                       first kind only are required; equal to 2 if
!                       radial functions of both the first and second
!                       kinds are desired
!               ienr  : number of d coeficient ratios required for
!                       convergence of the sum involved in the
!                       Bouwkamp eigenvalue procedure for the previous
!                       value of l
!               kindd : kind value for double precision
!               kindq : kind value for quadruple precision
!               ndec  : number of decimal digits for kind = knd
!
!     output:   eigval: converged eigenvalue
!               enr   : array of limd values
!               idigc : number of decimal digits of convergence
!                       for the Boouwkamp procedure
!               itestm: maximum number of decimal digits of agreement
!                       for one of the d coefficient ratios calculated
!                       forward and for the same d coefficient ratio
!                       calculated backward using recursion
!               ienr  : number of d coeficient ratios required for
!                       convergence of the sum involved in the
!                       Bouwkamp eigenvalue procedure for the current
!                       value of l
!               kflag : flag used to control the arithmetic used in
!                       the Bouwkamp procedure.
!                       Initially set equal to 0 for l = m when running
!                       cprofcn in double precision arithmetic with knd1
!                       equal to the kind value for quadruple precision.
!                       Set equal to one when running cprofcn in double
!                       precision arithmetic but using quadruple
!                       precision for the Bouwkamp procedure. Here knd1
!                       has been set equal to the kind value for
!                       quadruple precision.
!                       set equal to 2 when running cprofcn where knd = knd1.
!
!
        use param
!
!  real(knd) scalars and complex(knd) scalars and arrays
        real(knd) dec, eigdec, eigtest, eigtests
        real(knd1) dec1, eigdec2, eigtest1, eigtest1s
        complex(knd) c, cora, corb, de, dl, eigval, eigvals, enrc, &
                     enrcs, enrs
        complex(knd) blist(maxd), eigst, enr(maxd), glist(maxd), enrf(maxd)
        complex(knd1) cora1, cora1s, corb1, corb1s, de1, dl1, eigval1, &
                      eigval1s, enrc1, enrc1s, blist1(maxd), enr1(maxd), &
                      glist1(maxd)
!
        eigst = eigval
        eigvals = eigval
        dec = 10.0e0_knd ** (-ndec - 1)
        nbp = int(2.0e0_knd * (abs(real(c)) + abs(aimag(c))) / 3.14e0_knd)
        eigdec = dec * 10.0e0_knd
        ndec1 = precision(eigval1)
        dec1 = 10.0e0_knd1 ** (-ndec1 - 1)
        eigdec2 = dec * 10.0e0_knd
        eigval1 = eigval
        if(l == m .or. l == m + 1) imax = 1
        ncsave = 0
        ienrs = ienr
        li = l - m + 1
        lm2 = (l - m) / 2
        limdb = 2 * ienr + 20 + 2 * int(aimag(c))
        if(l == m .or. l == m + 1) limdb = 2 * ienr
        if(limdb > limd - 2) limdb = limd - 2
!
!  begin Bouwkamp procedure
        ix = l - m - 2 * lm2
        ifc = 1
        ifcs = 1
        lim2 = limdb / 2 - ix
        iglim = lim2 + 1
        irio = lm2 + 1
        iw1 = lm2 + 2
        li = l - m + 1
        eigtests = 1.0e0_knd
        eigtest1s = 1.0e0_knd1
        if(kflag == 1) go to 110
10      enr(1) = eigval - glist(1)
        if(irio < 2) go to 30
!
!  evaluate the continued fraction
          do 20 i = 1, irio - 1
          enr(i + 1) = -blist(i) / enr(i) - glist(i + 1) + eigval
20        continue
30      enr(lim2) = -blist(lim2) / (glist(iglim) - eigval)
        iw15 = lim2 - 1
        ip = iw1 + iw15
!
!  evaluate the continued fraction
          do 40 i = iw1, iw15
          ipi = ip - i
          enr(ipi) = -blist(ipi) / (glist(ipi + 1) - eigval + enr(ipi + 1))
40        continue
50      enrc = -blist(irio) / (glist(irio + 1) - eigval + enr(irio + 1))
        de = enrc * enrc / blist(irio)
        corb = de
        if(lim2 < iw1) go to 70
!
!  compute denominator in the Bouwkamp
          do 60 i = iw1, lim2
          de = enr(i) * enr(i) / blist(i) * de
          corb = corb + de
          if(abs(de / corb) < dec .and. (l == m .or. l == m + 1 .or.  &
             i > ienr - 20)) go to 70
60        continue
70      if((l == m .or. l == m + 1) .and. i > imax) imax = i
        if(l /= m .and. l /= m + 1 .and. i > ienr) ienr = i
        de = (1.0e0_knd, 0.0e0_knd)
        cora = de
!
!  compute the denominator in the Bouwkamp
          do 80 i = 1, lm2
          de = blist(irio - i) / (enr(irio - i) * enr(irio - i)) * de
          cora = cora + de
          if(abs(de / cora) < dec) go to 90
80        continue
!
!  compute the correction to the eigenvalue
90      dl = (enrc - enr(irio)) / (cora + corb)
!
!  eigenvalue accurate enough?
        eigtest = abs(dl / eigval)
        if(eigtest <= eigdec) then
        eigval = dl + eigval
        enrcs = enrc
        enrs = enr(irio)
        ncsave = ndec
        go to 100
        end if
        if(eigtest > eigtests * 100.0e0_knd) go to 95
        eigval = eigval + dl
          if(eigtest < eigtests) then
          eigtests = eigtest
          enrcs = enrc
          enrs = enr(irio)
          eigvals = eigval
          ifcs = ifc
          end if
        ifc = ifc + 1
        if(ifc < 20) go to 10
95      ifc = ifcs
        eigval = eigvals
        ncsave = -int(log10(eigtests + dec))
100     continue
        icomp = -int(log10(abs((eigval - eigst) / eigst) + dec))
        if(l == m .or. l == m + 1) ienr = imax
        idigc = ncsave
        if(idigc > ndec) idigc = ndec
        if(kflag == 0 .and. idigc < ndec) kflag = 1
          if(kflag == 1) then
          eigval1 = eigst
          ifc = 1
          ncsave = 0
          ienr = ienrs
          if(l == m .or. l == m + 1) imax = 1
          go to 110
          end if
        itestm = -int(log10(abs((enrcs - enrs) / enrs)+ &
                 dec))
        if(itestm > ndec) itestm = ndec
        idcoef = 0
        if(itestm < ndec) idcoef = 1
        go to 150
110     enr1(1) = eigval1 - glist1(1)
        if(lm2 < 1) go to 115
!
!  evaluate the continued fraction
          do i = 1, lm2
          enr1(i + 1) = -blist1(i) / enr1(i) - glist1(i + 1) + eigval1
          end do
115     enr1(lim2) = -blist1(lim2) / (glist1(iglim) - eigval1)
        iw15 = lim2 - 1
        ip = iw1 + iw15
        if(iw15 < iw1) go to 120
!
!  evaluate the continued fraction
          do i = iw1, iw15
          ipi = ip - i
          enr1(ipi) = -blist1(ipi) / (glist1(ipi + 1) - eigval1 + enr1(ipi + 1))
          end do
120     enrc1 = -blist1(irio) / (glist1(irio + 1) - eigval1 + enr1(irio + 1))
        de1 = enrc1 * enrc1 / blist1(irio)
        corb1 = de1
        if(lim2 < iw1) go to 125
!
!  compute first sum in the denominator of the correction
          do i = iw1, lim2
          de1 = enr1(i) * enr1(i) / blist1(i) * de1
          corb1 = corb1 + de1
          if(abs(de1 / corb1) < dec1 .and. (l == m .or. l == m + 1 .or.  &
             i > ienr - 20)) go to 125
          end do
125     if((l == m .or. l == m + 1) .and. i > imax) imax = i
        if(l > m + 1 .and. i > ienr) ienr = i
        de1 = (1.0e0_knd1, 0.0e0_knd1)
        cora1 = de1
        if(lm2 == 0) go to 130
!
!  compute second term in the denominator of the correction
          do i = 1, lm2
          de1 = blist1(irio - i) / (enr1(irio - i) * enr1(irio - i)) * de1
          cora1 = cora1 + de1
          if(abs(de1 / cora1) < dec1) go to 130
          end do
!
!  compute the correction to the eigenvalue
130     dl1 = (enrc1 - enr1(irio) + dec1) / (cora1 + corb1 + dec1)
        eigtest1 = abs(dl1 / eigval1)
!
!  eigenvalue accurate enough?
        if(eigtest1 < eigdec2) then
        eigval1 = dl1 + eigval1
        ncsave = ndec
        go to 140
        end if
        if(eigtest1 > eigtest1s * 100.0e0_knd1) go to 135
        eigval1 = eigval1 + dl1
          if(eigtest1 < eigtest1s) then
          eigtest1s = eigtest1
          cora1s = cora1
          corb1s = corb1
          enrc1s = enrc1
          eigval1s = eigval1
          ifcs = ifc
          end if
        ifc = ifc + 1
        if(ifc < 20) go to 110
135     ifc = ifcs
        eigval1 = eigval1s
        cora1 = cora1s
        corb1 = corb1s
        enrc1 = enrc1s
        ncsave = -int(log10(eigtest1s + dec))
140     continue
        icomp = -int(log10(abs((eigval1 - eigst) / eigst) + dec))
        idigc = ncsave
        eigval = eigval1
        int5 = -int(log10(abs((cora1 + corb1) * eigval1 / enrc1) + dec))
        if(kflag == 1 .and. l - m + 1 > 4 * nbp / 3 .and. int5 < 1) kflag = 2
        idcoef = 1
150     continue
if (debug) then
        if(knd == kindd) then
          if(ioprad /= 0) write(40, 160) l, eigst, idigc, ifc
          if(ioprad == 0) write(50, 160) l, eigst, idigc, ifc
160       format(1x,'l = ',i4,' eigen. estimate',e23.14, e23.14, &
                 /,10x,'conv. to',i3,' dig. at ifc =',i3)
          end if
          if(knd == kindq) then
          if(ioprad /= 0) write(40, 170) l, eigst, idigc, ifc
          if(ioprad == 0) write(50, 170) l, eigst, idigc, ifc
170       format(1x,'l = ',i4,' eigen. estimate',e39.30, e39.30, &
                 /,10x,'conv. to',i3,' dig. at ifc =',i3)
          end if
end if
if (debug) then
              if(knd == kindd .and. ioprad /= 0) write(40, 180) eigval
              if(knd == kindd .and. ioprad == 0) write(50, 180) eigval
180           format(10x,'eigenvalue =',e23.14, e23.14)
              if(knd == kindq .and. ioprad /= 0) write(40, 190) eigval
              if(knd == kindq .and. ioprad == 0) write(50, 190) eigval
190           format(10x,'eigenvalue =',e39.30, e39.30)
end if
!
!  calculate the d coefficient ratios (enr)
        lim2 = limd / 2 - ix
        if(idcoef == 1) go to 210
        enr(1) = eigval - glist(1)
        if(lm2 < 2) go to 200
          do i = 1, lm2 - 1
          enr(i + 1) = -blist(i) / enr(i) - glist(i + 1) + eigval
          end do
200     continue
        enr(lim2) = -blist(lim2) / (glist(lim2 + 1) - eigval)
        ilow = lm2 + 1
        iupp = lim2 - 1
        ip = ilow + iupp
          do i = ilow, iupp
          ipi = ip - i
          enr(ipi) = -blist(ipi) / (glist(ipi + 1) - eigval + enr(ipi + 1))
          end do
if (debug) then
        if(ioprad /= 0) write(40, 250) lim2, itestm, lm2 + 1
        if(ioprad == 0) write(50, 250) lim2, itestm, lm2 + 1
end if
        go to 260
210     continue
        enr(lim2) = -blist(lim2) / (glist(lim2 + 1) - eigval)
          do i = lim2 - 1, 1,-1
          enr(i) = -blist(i) / (glist(i + 1) - eigval + enr(i + 1))
          end do
        enrf(1) = eigval - glist(1)
        nlim = 1
        itestm = -int(log10(abs((enrf(1) - enr(1)) / enrf(1))+ &
                 dec))
        if(itestm < 0) itestm = 0
        if(itestm >= ndec) go to 240
          do 230 n = 2, lim2
          enrf(n) = -blist(n - 1) / enrf(n - 1) - glist(n) + eigval
          itest = -int(log10(abs((enrf(n) - enr(n)) / enrf(n))+ &
                 dec))
          if(itest < 0) itest = 0
          if(itest < itestm - 4) go to 240
          if(itest <= itestm) go to 230
          itestm = itest
          nlim = n
          if(itestm >= ndec) go to 240
230       continue
240     nlimp = 2 * (nlim - 1) + ix
if (debug) then
        if(ioprad /= 0) write(40, 250) lim2, itestm, nlimp
        if(ioprad == 0) write(50, 250) lim2, itestm, nlimp
250     format(15x, i6,' "d" coefs.: forward and backward', &
               ' recursion match to ',i3,' digits at n = ',i6)
end if
          do n = 1, nlim
          enr(n) = enrf(n)
          end do
260     continue
270     continue
        if(icomp < 5) then
if (debug) then
        if(ioprad /= 0) write(40, 280) l, m, c, icomp
end if
if (warn) then
        write(60, 280) l, m, c, icomp
end if
280     format(1x,'possible error in eigenvalue at l = ',i5, 2x,'m = ',i5, 2x, &
               'c = ',e25.15, e25.15, /5x,'converged eigenvalue ', &
               'only agreed with its estimate to ', i3, ' digits')
        end if
        return
        end subroutine
!
!
        subroutine dnorm (l, m, c, limd, maxd, ndec, nex, ioprad, enr, sgn, d01, &
                          id01, dmfnorm, idmfe, dmlmf, idmlmfe, dfnorm, idfe, &
                          dmlf, idmlfe, jmf, nsubmf, jfla, nsubf)
!
!  purpose:     to compute unscaled d coefficient ratios using scaled
!               ratios, to compute the Morse-Feshbach normalization
!               factor for the d coefficients, to compute the Flammer
!               normalization factor for the d coefficients, to compute
!               the characteristic and exponent of the first d
!               coefficient, either d0 or d1, and to compute the sign
!               of the d coefficient for n=l-m.
!
!
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : complex c
!               limd    : approximately twice the maximum number
!                         of terms available to be taken in the sums
!               maxd    : dimension of enr array
!               ndec    : number of decimal digits for real(knd)
!               nex     : maximum exponent for real(knd)
!               ioprad  : = 0 if radial functions were not requested,
!                         equal to 1 or 2 otherwise
!               enr     : array of scaled d coefficient ratios
!
!     output:   enr     : array of unscaled d coefficient ratios
!                         enr(i) = ratio of the d coefficient with
!                         subscript 2*i+ix to the d coefficient with
!                         subscript 2*(i-1)+ix. Here ix =0 when l-m is
!                         even and ix=1 when l-m is odd.
!                         If the user needs the d coefficent ratios,
!                         they are available below right before
!                         statement 20.
!               sgn     : sign of the d coefficient for n=l-m
!               d01     : characteristic of ratio of first d
!                         coefficient (either d0 or d1, depending on
!                         whether l-m is even or odd) to the d
!                         coefficient of order equal to l-m
!               id01    : exponent corresponding to d01
!               dmfnorm : characteristic of the Morse-Feshbach
!                         normalization factor of the d coefficients.
!                         equal to the reciprocal of the characteristic
!                         of the d coefficient d(n = l - m) using
!                         this normalization for the angular functions
!               idmfe   : exponent associated with dmfnorm
!               dmlmf   : characteristic of the d coefficient with
!                         index l - m using the Morse and Feshbach
!                         normalization for the angular functions
!               idmmfe  : exponent associated with dmlmf
!               dfnorm  : characteristic of the Flammer normalization
!                         factor of the d coefficients. equal to the
!                         reciprocal of the characteristic of the d
!                         coefficient d(n = l - m) using this
!                         normalization for the angular functions
!               idfe    : exponent associated with dfnorm
!               dmlf    : characteristic of the d coefficient with
!                         index l - m using the Flammer normalization
!                         for the angular functions
!               idmfe   : exponent associated with dmlf
!               jmf     : maximum index of enr required for convergence
!                         of dmfnorm
!               nsubmf  : number of decimal digits of subtraction error
!                         incurred in calculating dmfnorm
!               jfla    : maximum index of enr required for convergence
!                         of dfnorm
!               nsubf   : number of decimal digits of subtraction error
!                         incurred in calculating dfnorm
!
        use param
!
!  real(knd) scalars and complex(knd) scalars and array
        real(knd) aj, arr, dec, ea, rm2, sgn, ten, teste, testeo, tmax
        complex(knd) c, csq, dfnorm, dmfnorm, dmlf, dmlmf, d01, term
        complex(knd) enr(maxd)
!
        rm2 = real(m + m, knd)
        ten = 10.0e0_knd
        dec = ten ** (-ndec - 1)
        nfac = nex / 3
        teste = 10.0e0_knd ** nfac
        testeo = 1.0e0_knd / teste
        csq = c * c
        lm2 = (l - m) / 2
        ix = l - m - 2 * lm2
        lim2 = limd / 2 - ix
        mml = m + m - 1 + ix
        sgn = 1.0e0_knd
          do 20 i = 1, lim2
          arr = real(ix + i + i, knd)
          ea = arr + arr + rm2
          enr(i) = (ea - 1.0e0_knd) * (ea + 1.0e0_knd) * enr(i) / ((arr + rm2)* &
                   (arr + rm2 - 1.0e0_knd) * csq)
          if(i > lm2) go to 10
          if(real(enr(i)) < 0.0e0_knd)  sgn = sgn * (-1.0e0_knd)
10        continue
20        continue
!
!  compute the Morse-Feshbach normalization factor
!    forward summation of series
        term = (1.0e0_knd, 0.0e0_knd)
        tmax = 1.0e0_knd
        dmfnorm = term
        jlow = l - m + 2
        jterm = lm2
        idmfe = 0
          do 30 j = jlow, limd, 2
          aj = real(j, knd)
          jterm = jterm + 1
          term = term * (aj + rm2) * enr(jterm) * (aj + rm2 - 1.0e0_knd)/ &
               (aj * (aj - 1.0e0_knd))
          if(abs(term) > tmax) tmax = abs(term)
          dmfnorm = dmfnorm + term
          if(abs(term) < dec) go to 40
30        continue
40      jlow = l - m
        jmf = jterm
!
!    backward summation of series
        if(jlow < 2) go to 60
        term = (1.0e0_knd, 0.0e0_knd)
        jterm = lm2
          do 50 j = jlow, 2,-2
          aj = real(j, knd)
          term = term * aj * (aj - 1.0e0_knd) / ((aj + rm2 &
               -1.0e0_knd) * (aj + rm2) * enr(jterm))
          if(abs(term) > tmax) tmax = abs(term)
          jterm = jterm - 1
          dmfnorm = dmfnorm + term
          if(abs(term) < dec) go to 60
50        continue
60      continue
        nsubmf = -int(log10(abs(dmfnorm) / tmax + dec))
        if(nsubmf < 0) nsubmf = 0
        if(nsubmf > ndec) nsubmf = ndec
        iterm = 0
        if(abs(dmfnorm) /= 0.0e0_knd) iterm = int(log10(abs(dmfnorm)))
        idmfe = idmfe + iterm
        dmfnorm = dmfnorm * (10.0e0_knd ** (-iterm))
        dmlmf = 1.0e0_knd / dmfnorm
        idmlmfe = -idmfe
if (debug) then
        if(ioprad /= 0) write(40, 70) jmf, nsubmf
70      format(28x'Morse-Feshbach norm. converged in ', &
               i6,' terms with ',i3,' digits subt. error.')
end if
!
!  compute the Flammer normalization factor
!    forward summation of series
        term = (1.0e0_knd, 0.0e0_knd)
        dfnorm = term
        tmax = 1.0e0_knd
          do 80 j = lm2 + 1, lim2
          jj = j + j + ix
          term = -term * enr(j) * real(jj + mml, knd)/ &
                real(jj - ix, knd)
          if(abs(term) > tmax) tmax = abs(term)
          dfnorm = dfnorm + term
          if(abs(term / dfnorm) < dec) go to 90
80        continue
90      continue
        jfla = min(j, lim2)
!
!    backward summation of series
        if(lm2 < 1) go to 110
        term = (1.0e0_knd, 0.0e0_knd)
          do 100 j = lm2, 1,-1
          jj = j + j + ix
          term = -term * (jj - ix) / (real(jj + mml, knd) &
               *enr(j))
          if(abs(term) > tmax) tmax = abs(term)
          dfnorm = dfnorm + term
          if(abs(term / dfnorm) < dec) go to 110
100       continue
110     continue
        nsubf = -int(log10(abs(dfnorm) / tmax + dec))
        if(nsubf < 0) nsubf = 0
        if(nsubf > ndec) nsubf = ndec
        iterm = 0
        if(abs(dfnorm) /= 0.0e0_knd) iterm = int(log10(abs(dfnorm)))
        idfe = iterm
        dfnorm = dfnorm * (10.0e0_knd ** (-iterm))
        dmlf = 1.0e0_knd / dfnorm
        idmlfe = -idfe
if (debug) then
        if(ioprad /= 0) write(40, 120) jfla, nsubf
120     format(28x,'Flammer norm. converged in ',i6,' terms; ', &
               'with ',i2,' digits subt. error.')
end if
!
!  compute the d0(c|ml) or d1(c|ml)
        id01 = 0
        d01 = (1.0e0_knd, 0.0e0_knd)
        if(lm2 == 0) go to 140
          do 130 kjl = 1, lm2
          kkjl = lm2 - kjl + 1
          d01 = d01 / enr(kkjl)
            if(abs(d01) > teste) then
            d01 = d01 * testeo
            id01 = id01 + nfac
            end if
            if(abs(d01) < testeo) then
            d01 = d01 * teste
            id01 = id01 - nfac
            end if
130       continue
        iterm = int(log10(abs(d01)))
        d01 = d01 * (10.0e0_knd ** (-iterm))
        id01 = id01 + iterm
140     continue
        return
        end subroutine
!
!
        subroutine dalt (l, m, cc, limdr, maxdr, maxmp, ndec, nex, ioppsum, &
                         eigval, enrneg, drhor, dneg, idneg, nsdneg, nsdrho)
!
!  purpose:     To calculate d ratios with negative subscripts
!               and d-rho ratios.
!  parameters:
!
!     input:    l       : l
!               m       : m
!               cc      : complex c
!               limdr   : number of ratios of successive d-rho
!                         coefficients calculated
!               maxdr   : dimension of drhor array
!               maxmp   : dimension of enrneg array
!               ndec    : number of decimal digits for real(knd)
!               nex     : maximum exponent available for real(knd)
!               ioppsum : integer index = 0 if no d rho coefficients
!                         are calculated (psum not needed for r2leg)
!               eigval  : eigenvalue
!
!     output:   enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               dneg    : characteristic of the ratio of the d
!                         coefficient with index -2m+ix to the
!                         d coefficient with index ix, where
!                         ix = 0 if l-m is even and ix = 1 if
!                         l-m is odd
!               idneg   : exponent (base 10) of dneg
!               nsdneg  : subtaction error in calculating dneg
!               nsdrho  : subtraction error in calculating drhor(1)
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) r, rm, rn, t, ten, teste, testeo, uterm, wterm
        real(knd) amnsdrho, ansdneg, ansdrho, asub, bsub
!
!  complex(knd) scalars and arrays
        complex(knd) cc, dneg, eigval
        complex(knd) enrneg(maxmp), drhor(maxdr)
        complex(knd) vterm
!
        ten = 10.0e0_knd
        nfac = nex / 3
        if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
        teste = ten ** (nfac)
        testeo = 1.0e0_knd / teste
!  if l-m is even, ix=0; if l-m is odd, ix=1
        ix = (l - m) - 2 * ((l - m) / 2)
        t = m + m - ix - ix
!
!  calculate ratios of d coefficients with negative subscripts
!
!       enrneg(k) = { d(-2m+2k-2)/d(-2m+2k), l-m even    }
!                   { d(-2m+1+2k-2)/d(-2m+1+2k), l-m odd }
!
        rm = m
          if(m == 0) then
          dneg = (1.0e0_knd, 0.0e0_knd)
          idneg = 0
          nsdneg = 0
          go to 30
          end if
          do 10 i = 1, m + 1
          enrneg(i) = (0.0e0_knd, 0.0e0_knd)
10        continue
!
!  first calculate enrneg(1)
        n = 2 - 2 * m + ix
        rn = n
        r = n + m + m
        uterm = r * (r - 1.0e0_knd) / ((rn + r + 1.0e0_knd) * (rn + r - 1.0e0_knd))
        r = n + m - 2
        vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm- &
              1.0e0_knd) / ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r- &
              1.0e0_knd)) + (r * (r + 1.0e0_knd) - eigval) / (cc * cc)
!
!       calculations continue up to and including
!       enrneg(k=m) = { d(-2)/d(0), l-m even }
!                     { d(-1)/d(1), l-m odd  }
!
        enrneg(1) = -uterm / vterm
        dneg = enrneg(1)
        idneg = 0
        ansdneg = 0.0e0_knd
        nsdneg = 0
!
!  backward recursion beginning with enrneg(1) and
!  ending with enrneg(m)
!
          do i = 2, 2 * m - 2, 2
          ii = i - 2 * m + ix
          n = ii + 2
          j = i / 2
          rn = n
          r = n + m + m
          uterm = r * (r - 1.0e0_knd) / ((rn + r + 1.0e0_knd) * (rn + r- &
                     1.0e0_knd))
          r = n + m - 2
          vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm- &
                1.0e0_knd) / ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r- &
                1.0e0_knd)) + (r * (r + 1.0e0_knd) - eigval) / (cc * cc)
          r = n - 4
          wterm = (r + 2.0e0_knd) * (r + 1.0e0_knd) / ((r + r + rm + rm+ &
                      3.0e0_knd) * (r + r + rm + rm + 1.0e0_knd))
          enrneg(j + 1) = -uterm / (wterm * enrneg(j) + vterm)
          dneg = dneg * enrneg(j + 1)
            if(wterm * enrneg(j) * vterm /= (0.0e0_knd, 0.0e0)) &
               then
            asub = log10(abs(vterm / (vterm+ &
                             wterm * enrneg(j))))
            if(asub > 0.0e0_knd) ansdneg = ansdneg + asub
            bsub = log10(abs(vterm / (wterm * enrneg(j))))
            if(bsub > 0.0e0_knd) ansdneg = max(0.0e0_knd, ansdneg - bsub)
            if(int(ansdneg) > nsdneg) nsdneg = ansdneg
            end if
            if(abs(dneg) > teste) then
            dneg = dneg * testeo
            idneg = idneg + nfac
            end if
            if(abs(dneg) < testeo) then
            dneg = dneg * teste
            idneg = idneg - nfac
            end if
          end do
          if(nsdneg > ndec) nsdneg = ndec
        iterm = int(log10(abs(dneg)))
        dneg = dneg * (ten ** (-iterm))
        idneg = idneg + iterm
!
!  calculate ratios of d rho coefficients
!
!       drhor(k-m) = { d(rho|2k)/d(rh0|2k-2), l-m even  }
!                    { d(rho|2k-1)/d(rho|2k-3), l-m odd }
!
30      if(ioppsum == 0) go to 60
          ansdrho = 0.0e0_knd
          amnsdrho = 0.0e0_knd
          nsdrho = 0
          mnsdrho = 0
          do 40 i = 1, limdr
          drhor(i) = (0.0e0_knd, 0.0e0_knd)
40        continue
          do 50 i = 2 * limdr, 6,-2
          n = 4 - i + ix - m - m
          ii = (i - 2) / 2
          rn = n
          r = n + m + m
          uterm = r * (r - 1.0e0_knd) / ((rn + r + 1.0e0_knd) * (rn + r - 1.0e0_knd))
          r = n + m - 2
          vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm- &
                 1.0e0_knd) / ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r- &
                 1.0e0_knd)) + (r * (r + 1.0e0_knd) - eigval) / (cc * cc)
          r = n - 4
          wterm = (r + 2.0e0_knd) * (r + 1.0e0_knd) / ((r + r + rm + rm + 3.0e0_knd)* &
                (r + r + rm + rm + 1.0e0_knd))
          drhor(ii) = -uterm / (wterm * drhor(ii + 1) + vterm)
            if(wterm * drhor(ii + 1) * vterm /= (0.0e0_knd, 0.0e0)) &
               then
            asub = log10(abs(vterm / (vterm + wterm * drhor(ii + 1))))
            if(asub > 0.0e0_knd) ansdrho = ansdrho + asub
            bsub = log10(abs(vterm / (wterm * drhor(ii + 1))))
            if(bsub > 0.0e0_knd) ansdrho = max(0.0e0_knd, ansdrho - bsub)
            if(ansdrho > amnsdrho) amnsdrho = ansdrho
            end if
50        continue
        n = -2 * m + ix
        r = n + m - 2
        vterm = (2.0e0_knd * r * (r + 1.0e0_knd) - 2.0e0_knd * rm * rm - 1.0e0_knd)/ &
              ((2.0e0_knd * r + 3.0e0_knd) * (2.0e0_knd * r - 1.0e0_knd))+ &
              (r * (r + 1.0e0_knd) - eigval) / (cc * cc)
        r = n - 4
        wterm = (r + 2.0e0_knd) * (r + 1.0e0_knd) / ((r + r + rm + rm + 3.0e0_knd)* &
              (r + r + rm + rm + 1.0e0_knd))
!
!       the final value of ii is 1;
!       drhor(1) has a special value:
!       drhor(1) = { d(rho|2m+2)/d(-2m), l-m even  }
!                  { d(rho|2m+1)/d(-2m+1), l-m odd }
!
        drhor(1) = 1.0e0_knd / ((t - 1.0e0_knd) * (t + 1.0e0_knd)* &
                 (wterm * drhor(2) + vterm))
          if(wterm * drhor(2) * vterm /= (0.0e0_knd, 0.0e0)) &
               then
          asub = log10(abs(vterm / (vterm + wterm * drhor(2))))
          if(asub > 0.0e0_knd) ansdrho = ansdrho + asub
          bsub = log10(abs(vterm / (wterm * drhor(2))))
          if(bsub > 0.0e0_knd) ansdrho = max(0.0e0_knd, ansdrho - bsub)
          if(ansdrho > amnsdrho) amnsdrho = ansdrho
          end if
        nsdrho = int(amnsdrho) + 1
        nsdrhor1 = asub
        if(nsdrhor1 < 0) nsdrhor1 = 0
        if(ix == 1) drhor(1) = -drhor(1)
60      continue
        return
        end subroutine
!
!
    subroutine gauss (n, ndec, x, w)
!
!  purpose:     To evaluate the coordinates and weighting factors
!               for an nth order Gaussian quadrature
!
!  parameters:
!
!               n   : order of quadrature
!               ndec: number of decimal digits for real(knd)
!
!     output:   x   : coordinate values for quadrature
!               w   : weighting factors
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) delta, der, pi, s, t, test, u, v, z
        real(knd) x(n), w(n)
!
        test = 10.0e0_knd ** (-ndec)
        imax = (n + 1) / 2
        pi = acos(-1.0_knd)
          do 40 i = 1, imax
          z = cos(pi * (i - 0.25e0_knd) / (n + 0.5e0_knd))
            do 20 j = 1, 30
            u = 0.0e0_knd
            v = 1.0e0_knd
              do 10 k = 1, n
              t = u
              u = v
              v = ((k + k - 1) * z * u - (k - 1) * t) / k
10            continue
            s = z * z - 1.0e0_knd
            der = n * (z * v - u) / s
            delta = -v / der - 0.5e0_knd * v * v * ((n * n * s - n * z * z - n) * v+ &
                    2.0e0_knd * n * z * u) / (der * der * der * s * s)
            z = z + delta
            if(abs(delta / z) < test) go to 30
20          continue
30        continue
          x(i) = -z
          x(n + 1 - i) = z
          w(i) = 2.0e0_knd / ((1.0e0_knd - z * z) * der * der)
          w(n + 1 - i) = w(i)
40        continue
        return
        end subroutine
!
!
        subroutine pleg (m, lim, maxp, ndec, nex, limcsav, iopd, barg, narg, &
                         maxt, pr, pdr, pdnorm, ipdnorm, pnorm, ipnorm, alpha, &
                         beta, gamma, coefa, coefb, coefc, coefd, coefe)
!
!  purpose:     To calculate ratios of successive associated Legendre
!               functions of the first kind for given arguments barg.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Legendre functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lim    : two less than the number of associated Legendre
!                        function ratios calculated for given arguments
!               maxp   : dimension of alpha, beta, gamma, coefa, coefb,
!                        coefc, coefd, and coefe arrays and second
!                        dimension of pr and pdr arrays
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent available for real(knd)
!               limcsav: integer equal to the number of coefficients in
!                        each of the arrays alpha, beta, gamma, coefa,
!                        coefb, coefc, coefd, and coefe arrays that
!                        have already been calculated in earlier calls
!                        to pleg for this value of m and will not be
!                        calculated again. [Note that the minimum
!                        array index for the coefficients is 3 and
!                        the maximum array index is limcsav+2]
!               iopd   : integer that is set = 0 if derivatives of
!                        Legendre functions (i.e., their ratios)
!                        are not required when iopang = 1 and the
!                        first derivatives of the angular functions
!                        are not requested.
!                        iopd is set = 1 when iopang = 2 and pleg is
!                        also being used to obtain ratios of first
!                        derivatives of Legendre functions for use in
!                        computing the first derivatives of the angular
!                        functions.
!                        iopd is set = 2 when pleg is being used to
!                        compute ratios of Legendre functions for use in
!                        the calculation of the denominator term used in
!                        calculating the radial functions of the second
!                        kind and their first derivatives in r2eta.
!                        iopd is set = 3 when pleg is being used to
!                        compute ratios of both the Legendre functions
!                        and their first derivatives for use in the
!                        calculation of the numerator terms used
!                        in r2eta and in r2leg to calculate the radial
!                        functions of the second kind and their first
!                        deriatives.
!               barg   : array of narg values of eta for which Legendre
!                        functions are to be calculated
!               narg   : number of specified values of eta in barg array
!               maxt   : dimension of barg array
!
!     output:   pr     : array of ratios of successive first kind
!                        associated Legendre functions of the same
!                        parity
!               pdr    : array of ratios of successive derivatives of
!                        first kind associated Legendre functions of
!                        the same parity
!               pdnorm : array of characteristics of the first
!                        derivatives of associated Legendre functions
!                        of the first kind of order m and degree m
!               ipdnorm: array of exponents of the first derivatives
!                        of associated Legendre functions of the first
!                        kind of order m and degree m
!               pnorm  : array of characteristics of the associated
!                        Legendre functions of the first kind of order m
!                        and degree m
!               ipnorm : array of exponents of the associated Legendre
!                        functions of the first kind of order m and
!                        degree m
!
!     input/output:
!               alpha  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               beta   : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               gamma  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               coefa  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefb  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefc  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefd  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefe  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) adec, ajterm, am2p1, anden1, anden2, an2tnp1, bargs, den, rm, &
                  rm2, temp1, temp2, temp3, ten, term, teste, testeo, ta, tb, tc, &
                  t1, t2
        real(knd) alpha(maxp), barg(maxt), beta(maxp), coefa(maxp), &
                  coefb(maxp), coefc(maxp), coefd(maxp), coefe(maxp), &
                  gamma(maxp), pdnorm(maxt), pdr(maxt, maxp), pdr1(maxp), &
                  pr(maxt, maxp), pnorm(maxt)
!
!  integer array
        dimension ipdnorm(maxt), ipnorm(maxt)
!
        ten = 10.0e0_knd
        adec = ten ** (-ndec + 2)
        nfac = nex / 2
        if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
        teste = ten ** nfac
        testeo = 1.0e0_knd / teste
        rm = m
        rm2 = m * m
        am2p1 = m + m + 1
        m2 = m + m
        m2m1 = m2 - 1
        mm1 = m - 1
        mm2 = m - 2
        msqp1 = 2 * m * m + 1
!
!  calculate the coefficients alpha(j), beta(j), and gamma(j) for
!  the three term recursion relating the Legendre function ratios
!
!              m                m
!   pr(k,j) = P    (barg(k)) / P  (barg(k))
!              m+j-1            m+j-3
!
!  and calculate the coefficients coefa(j), coefb(j), coefc(j),
!  coefd(j), and coefe(j) in the expression used to calculate
!  ratios of Legendre function derivatives
!
!               m                 m
!   pdr(k,j) = P'    (barg(k)) / P'  (barg(k))
!               m+j-1             m+j-3
!
!  Note that pr(k,1) and pr(k,2) are not ratios but actually equal to
!   m      m                                              m       m
!  P  and P   . Also, pdr(k,1) and pdr(k,2) are equal to P'  and P' .
!   m      m+1                                            m       m+1
!
          if(limcsav >= lim) go to 30
          do 10 j = limcsav + 3, lim + 2
          n = m + j - 3
          n2 = n + n
          n2p3 = n2 + 3
          n2p1 = n2 + 1
          n2m1 = n2 - 1
          nmmm2 = n - mm2
          nmmm1 = n - mm1
          npmm1 = n + mm1
          npm = n + m
          npmm1 = n + mm1
          npmp1 = n + m + 1
          npmp2 = npmp1 + 1
          an2tnp1 = real(2 * n, knd) * real((n + 1), knd)
          anden1 = real(nmmm2, knd) * real(nmmm1, knd)
          anden2 = real(n2m1, knd) * anden1
          alpha(j) = real(n2p3, knd) * real(n2p1, knd) / anden1
          beta(j) = real(n2p1, knd) * (real(msqp1, knd) - an2tnp1) / anden2
          gamma(j) = -real(n2p3, knd) * real(npm, knd) * real(npmm1, knd) / anden2
          coefa(j) = -real(npmp2, knd) / real(nmmm1, knd)
          coefb(j) = real(n2p3, knd) * real((n + 2), knd) / anden1
          coefc(j) = -real(npmp1, knd) * real(npmp2, knd) / anden1
          coefd(j) = real(npmp1, knd) / real(nmmm2, knd)
          coefe(j) = -real((n + 1), knd) * real(n2p3, knd) / anden1
10        continue
        gamma(3) = 0.0e0_knd
        gamma(4) = 0.0e0_knd
        term = 1.0e0_knd
        iterm = 0
        if(m < 2) go to 30
          do jm = 2, m
          term = (jm + jm - 1) * term
            if(term > teste) then
            term = term * testeo
            iterm = iterm + nfac
            end if
          end do
        jterm = int(log10(term))
        term = term * (ten ** (-jterm))
        iterm = iterm + jterm
30      continue
!
!   calculate the ratios of successive Legendre functions of the same
!   parity using the three term recursion relationship
!
!   pr(k,j) = alpha(j)*barg(k)*barg(k) + beta(j) + gamma(j)/pr(k,j-2)
!
          do 140 k = 1, narg
          pnorm(k) = term
          ipnorm(k) = iterm
          pdnorm(k) = term
          ipdnorm(k) = iterm
!
!   define the first two ratios equal to unity and (2m+1)*barg(k)
          pr(k, 1) = 1.0e0_knd
          pr(k, 2) = am2p1 * barg(k)
          jdelta = 1
          if(abs(barg(k)) < adec) jdelta = 2
          bargs = barg(k) * barg(k)
            do 40 j = 3, lim + 2, jdelta
            pr(k, j) = alpha(j) * bargs + beta(j) + (gamma(j) / pr(k, j - 2))
40          continue
!
!   calculate the corresponding ratios of first derviatives of
!   successive Legendre functions of the same parity using the
!   following relationship (except for (1) eta equal to zero or unity,
!   where special expressions are used and (2) when the magnitude of the
!   argument barg is <= 0.1 or (m+1)*barg*barg - 1 < 0.01, where recursion
!   on the ratios of successive first derivatives of the same parity is
!   used instead)
!
!              (coefa(j)+coefb(j)*barg(k)*barg(k))*pr(k,j)+coefc(j)
!   pdr(k,j) = ----------------------------------------------------
!                  pr(k,j)+coefd(j)+coefe(j)*barg(k)*barg(k)
!
          if(iopd == 0 .or. iopd == 2) go to 120
          pdr(k, 1) = 1.0e0_knd
          pdr(k, 2) = 1.0e0_knd
          if(abs(barg(k)) >= adec) go to 50
          pdr(k, 2) = am2p1
            do j = 4, lim + 2, 2
            pdr(k, j) = -real(m2m1 + j, knd) / real(j - 2, knd)
            end do
          go to 140
50        if(abs(abs(barg(k)) - 1.0e0_knd) >= adec) go to 70
          if(m == 0) go to 60
          if(m /= 2) go to 130
          pdr(k, 1) = -2.0e0_knd * barg(k)
          go to 80
60        temp1 = 1.0e0_knd
          temp2 = 3.0e0_knd
          pdr(k, 2) = 1.0e0_knd
          pdr(k, 3) = 3.0e0_knd * barg(k)
            do j = 4, lim + 2
            temp3 = temp2 + j - 1
            pdr(k, j) = temp3 / temp1
            temp1 = temp2
            temp2 = temp3
            end do
          go to 140
70        if(m /= 0) go to 80
          pdr(k, 1) = 1.0e0_knd
          pdr(k, 2) = 1.0e0_knd
          pdr(k, 3) = 3.0e0_knd * barg(k)
          jlow = 4
          go to 90
80        pdr(k, 2) = am2p1 * ((rm + 1.0e0_knd) * bargs - 1.0e0_knd) / (rm * barg(k))
          if(pdr(k, 2) == 0.0e0_knd) pdr(k, 2) = ten ** (-ndec)
          jlow = 3
          if(abs((rm + 1.0e0_knd) * bargs - 1.0e0_knd) < 0.01e0_knd) go to 110
90        continue
          if(abs(barg(k)) <= 0.1e0_knd) go to 110
            do 100 j = jlow, lim + 2
            den = (pr(k, j) + coefd(j) + coefe(j) * bargs)
            if(den == 0.0e0_knd) den = ten ** (- ndec)
            pdr(k, j) = ((coefa(j) + coefb(j) * bargs) * pr(k, j) + coefc(j)) / den
100         continue
          go to 120
110       continue
          if(m /= 0) pdr1(1) = pdr(k, 2)
          if(m == 0) pdr1(2) = pdr(k, 3)
            do j = jlow - 1, lim + 1
            n = j + m - 1
            t1 = bargs - 1.0e0_knd
            t2 = (n * n * t1 + rm2)
            ta = j * t2
            tb = (n + n + 1) * barg(k) * (t2 + n * t1)
            tc = -(n + m) * (t2 + (n + n + 1) * t1)
            pdr1(j) = (tb + tc / pdr1(j - 1)) / ta
            end do
            do j = jlow, lim + 2
            pdr(k, j) = pdr1(j - 2) * pdr1(j - 1)
            end do
120       if(m == 0 .or. iopd == 2 .or. iopd == 3) go to 140
          if(abs(abs(barg(k)) - 1.0e0_knd) < adec) go to 130
          ajterm = rm * log10(1.0e0_knd - bargs) / 2.0e0_knd
          jterm = int(ajterm)
          ipnorm(k) = ipnorm(k) + jterm
          pnorm(k) = pnorm(k) * (ten ** (ajterm - jterm))
          if(iopd == 0) go to 130
          ajterm = log10(rm * abs(barg(k))) + (rm - 2.0e0_knd)* &
                 log10(1.0e0_knd - bargs) / 2.0e0_knd
          jterm = int(ajterm)
          ipdnorm(k) = ipdnorm(k) + jterm
          pdnorm(k) = -pdnorm(k) * (ten ** (ajterm - jterm))
          if(barg(k) < 0.0e0_knd) pdnorm(k) = -pdnorm(k)
          go to 140
130       pnorm(k) = 0.0e0_knd
          ipnorm(k) = 0
          if(m /= 2) pdnorm(k) = 0.0e0_knd
          if(m /= 2) ipdnorm(k) = 0
140       continue
        return
        end subroutine
!
!
        subroutine qleg (m, lnum, limq, maxq, x1, ndec, qdr, qdml, iqdml, qdl, &
                         iqdl, qr, qml, iqml, ql, iql, termpq, itermpq)
!
!  purpose:     to calculate ratios of successive associated Legendre
!               functions of the second kind for given c,x, and m.
!               to calculate corresponding ratios of their first
!               derivatives. to calculate the characteristics and
!               exponents of both the Legendre functions of the second
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lnum   : number of l values desired (=lmax+1);
!                        also equal to the dimension of the arrays
!                        ql, iql, qdl, and iqdl
!               limq   : the number of associated Legendre function
!                        ratios calculated for given m,lnum,c,ndec,
!                        and x1
!               maxq   : dimension of qr and qdr arrays
!               x1     : x - 1.0e0_knd
!               ndec   : number of decimal digits in real(knd) arithmetic
!
!     output:   qdr    : ratios of derivatives of successive
!                        associated Legendre functions of the second
!                        kind
!               qdml   : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree m-1
!               iqdml  : exponent corresponding to qdml
!               qdl    : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree m
!               iqdl   : exponent corresponding to qdl
!               qr     : ratios of successive associated Legendre
!                        functions of the second kind
!               qml    : characteristic of the associated Legendre
!                        function of the second kind with order m
!                        and degree m-1
!               iqml   : exponent corresponding to qml
!               ql     : characteristic of the associated Legendre
!                        function of the second kind with order m and
!                                                           -m/2
!                        degree m, scaled by (2m-1)c!(x*x-1)
!               iql    : exponent corresponding to ql
!               termpq : characteristic of the relative size of the
!                        maximum terms in the positive degree q series
!                        and the p series used to calculate r2 and r2d
!                        in subroutine r2leg
!               itermpq: exponent corresponding to termpq
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) ajm, dec, qdml, qlow, qml, qupp, q00, q11, rin, rm, &
                term, termpq, tjm, tm, tmr, x, x1, x1d, xsqr
        real(knd) qdl(lnum), qdr(maxq), ql(lnum), qr(maxq)
!
!  integer arrays
        dimension iqdl(lnum), iql(lnum)
!
        dec = 10.0e0_knd ** (-ndec)
        rm = real(m, knd)
        tm = rm + rm
        tmr = tm / (tm + 1.0e0_knd)
        qflag = 0
        x = x1 + 1.0e0_knd
        x1d = (x + 1.0e0_knd) * x1
        xsqr = sqrt(x1d)
        mxqrest = limq + ndec * int((1.0e0_knd - 1.0e0_knd / log10(x - xsqr)))
        if(m == 0) mlimq = 50000 * ndec + limq
        if(m == 1) mlimq = 12000 * ndec + limq
        if(m == 2) mlimq = 5000 * ndec + limq
        if(m == 3) mlimq = 600 * ndec + limq
        if(m >= 4) mlimq = 100 * ndec + limq
        if(m == 1 .and. x1 < 1.0e-9_knd) mlimq = 50000 * ndec + limq
        mxqr = min(mxqrest, mlimq)
        mxqrpm = mxqr + m
if (debug) then
        write(40, 5) mxqrpm
5       format(15x,'used backward recursion to calculate ratios of q', &
               ' functions starting at order',i8)
end if
!
!                              m    m
!  calculate ratios qr(k+m) = q  / q    ,k=m+limq to k=m+1 using
!                              k    k-1
!
!                                         m       m               1/2
!  backward recursion from qr(maxm+2m) = q     / q      =x-(x*x-1)
!                                         mxqr+m  mxqr+m-1
!
        qupp = x - xsqr
          do jn = mxqr + m, limq + m + 1,-1
          rin = real(jn, knd)
          qlow = (rin + rm - 1.0e0_knd) / (x * (rin + rin - 1.0e0_knd) &
               -(rin - rm) * qupp)
          qupp = qlow
          end do
        qr(limq + m + m) = qupp
          do 10 jn = limq + m, m + 2,-1
          rin = real(jn, knd)
          qr(jn + m - 1) = (rin + rm - 1.0e0_knd) / (x * (rin + rin - 1.0e0_knd) &
                     -(rin - rm) * qr(jn + m))
10        continue
!
!                              m     m
!  calculate ratios qr(k+m) = q   / q   ,k=m-1 to k=-m+1 using
!                              k-1   k
!
!                                       m      m
!  backward recursion from qr(m-1+m) = q    / q     = x
!                                       m-2    m-1
!
20      if(m == 0) go to 100
        qr(m + m - 1) = x
        if(m == 1) go to 40
          do 30 jn = m - 1, 2 - m,-1
          rin = real(jn, knd)
          qr(jn + m - 1) = (x * (rin + rin - 1.0e0_knd) &
                     -((rin - rm) / qr(jn + m))) / (rin + rm - 1.0e0_knd)
30        continue
40      continue
!
!                  m
!  calculation of q    , m > 0 by forward division of qr ratios
!                  m-1
!
!                 m
!  starting with q  calculated from its closed form expression.
!                 0
!
        qml = rm * log10(x + 1.0e0_knd) - log10(2.0e0_knd)
        term = rm * log10(x1 / (x + 1.0e0_knd))
        if(term < -ndec) go to 50
        qml = qml + log10(1.0e0_knd - ((x1 / (x + 1.0e0_knd)) ** m))
50      continue
        term = 1.0e0_knd
        iterm = 0
        if(m < 2) go to 70
          do jm = 2, m
          ajm = real(jm, knd)
          term = term * (ajm - 1.0e0_knd) / (ajm + ajm - 1.0e0_knd)
          if(term > dec) go to 60
          term = term / dec
          iterm = iterm - ndec
60        continue
          end do
70      term = log10(term)
        qml = qml + term + iterm
        iqml = int(qml)
        qml = 10.0e0_knd ** (qml - iqml)
        if(2 * (m / 2) /= m) qml = -qml
        if(m < 2) go to 90
          do jm = 1, m - 1
          qml = qml / qr(jm + m)
          if(abs(qml) > dec) go to 80
          qml = qml / dec
          iqml = iqml - ndec
80        continue
          end do
90      continue
        iterm = int(log10(abs(qml)))
        qml = qml * 10.0e0_knd ** (-iterm)
        iqml = iqml + iterm
!
!                  m
!  calculation of q   by forward recursion in m starting with values
!                  m
!       0       1
!  for q   and q  obtained from their closed form expressions, scaled
!       0       1
!                    -m/2
!  by (2m-1)!!(x*x-1).
!
100     q00 = 0.5e0_knd * log((x + 1.0e0_knd) / x1)
        if(m /= 0) go to 110
        ql(1) = q00
        go to 130
110     q11 = x1d * q00 - x
        if(m /= 1) go to 120
        ql(1) = q11
        go to 130
120     qlow = q00
        qupp = q11
          do jm = 1, m - 1
          tjm = real(jm + jm, knd) / real(jm + jm + 1, knd)
          ql(1) = (x1d - tjm) * qupp + tjm * x1d * qlow
          qlow = qupp
          qupp = ql(1)
          end do
130     iql(1) = int(log10(abs(ql(1))))
        ql(1) = ql(1) * (10.0e0_knd ** (-iql(1)))
!
!  calculation of ratios of the first derivatives of q with respect
!  to x, using the relationships:
!
!                  m    m      [kx]qr(k+m)-(k+m)
!     qdr(k+m) = q'  / q'   =  ----------------- , k=m+lim to k=m+1
!                  k    k-1    [(k-m)]qr(k+m)-kx
!
!                  m      m    [(k-m)]qr(k+m)-kx
!     qdr(k+m) = q'   / q'  =  ----------------- , k=m-1 to k=-m+1
!                  k-1    k    [kx]qr(k+m)-(k+m)
!
          do jm = m + 1, m + limq
          ajm = real(jm, knd)
          qdr(jm + m) = (ajm * x * qr(jm + m) - (ajm + rm)) / ((ajm - rm) * qr(jm + m) - ajm * x)
          end do
!
        if(m == 0) go to 140
          do jm = 1 - m, m - 1
          ajm = real(jm, knd)
          qdr(jm + m) = (ajm * x * qr(jm + m) - (ajm - rm)) / ((ajm + rm) * qr(jm + m) - ajm * x)
          end do
140     continue
!
!                   m         m                      m        m
!  calculation of q'    and q'  from the values for q    and q  .
!                   m-1       m                      m-1      m
!
        if(m > 0) go to 150
        qdl(1) = -1.0e0_knd / x1d
        iqdl(1) = 0
        go to 160
150     qdml = -rm * x * qml / x1d
        iterm = int(log10(abs(qdml)))
        qdml = qdml * (10.0e0_knd ** (-iterm))
        iqdml = iqml + iterm
        qdl(1) = rm * (x * ql(1) - 2.0e0_knd * qml * (10.0e0_knd ** (iqml - iql(1))))/ &
                   x1d
        iqdl(1) = iql(1)
160     continue
        m2m1 = m + m - 1
          do jl = 2, lnum
          ql(jl) = ql(jl - 1) * qr(m2m1 + jl)
          iql(jl) = iql(jl - 1)
          if(abs(ql(jl)) > 1.0e-10_knd) go to 170
          ql(jl) = ql(jl) * 1.0e10_knd
          iql(jl) = iql(jl) - 10
170       qdl(jl) = qdl(jl - 1) * qdr(m2m1 + jl)
          iqdl(jl) = iqdl(jl - 1)
          if(abs(qdl(jl)) > 1.0e-10_knd) go to 180
          qdl(jl) = qdl(jl) * 1.0e10_knd
          iqdl(jl) = iqdl(jl) - 10
180       end do
        termpq = rm * log10(xsqr)
        itermpq = int(termpq)
        termpq = 10.0e0_knd ** (termpq - itermpq)
        return
        end subroutine
!
!
       subroutine pint (cc, m, lnum, x1, limint, maxint, maxlp, ndec, nex, &
                        lim1max, wg, xg, ngau, ngqs, limi, rpint1, rpint2, &
                        pint1, pint2, pint3, pint4, norme, pnorm, ipnorm, &
                        coefme, coefmo)
!
!  purpose:     to calculate integrals of the product of associated
!               Legendre functions and kernels containing spherical
!               Neumann functions and a window function. four
!               different kernel functions are involved leading to
!               integrals of four different types. the integrals are
!               calculated using gaussian quadrature.
!
!  parameters:
!
!     input:    cc     : complex c
!               m      : m
!               lnum   : number of l values desired
!               x1     : x - 1.0e0_knd
!               limint : number of integrals of each of the four types
!                        required
!               maxint : dimension of the integral arrays
!               maxlp  : dimension of characteristic and exponent
!                        arrays of integrals
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent available for real(knd)
!               lim1max: dimension of the spherical Bessel function
!                        array used in computing spherical Neumann
!                        functions
!               wg     : gaussian quadrature weighting factors
!               xg     : gaussian quadrature arguments
!               ngau   : order of gaussian quadrature
!               ngqs   : number of gaussian quadrature steps the
!                        integrand is divided into
!
!     output:   limi   : highest order of Legendre function for which
!                        integrals are calculated, since higher order
!                        integrals could overflow. series in subroutine
!                        r2int will be limited to order limi; initially
!                        set equal to maxint-4 or limint+1
!               rpint1 : array of ratios of successive integrals of
!                        the same parity of the first type (l-m even)
!                        or of the third type (l-m odd)
!               rpint2 : array of ratios of successive integrals of
!                        the same parity of the second type (l-m even)
!                        or of the fourth type (l-m odd)
!               pint1  : array of scaled values for the integrals of
!                        the first type
!               pint2  : array of scaled values for the integrals of
!                        the second type
!               pint3  : array of scaled values for the integrals of
!                        the third type
!               pint4  : array of scaled values for the integrals of
!                        the fourth type
!               norme  : scaling exponent for the spherical Neumann
!                        function of order m
!               pnorm  : array of characteristics for the scaling
!                        factors used for the associated Legendre
!                        functions
!               ipnorm : array of exponents for the scaling factors
!                        used for the associated Legendre functions
!               coefme : coefficient for the expression for r2 and r2d
!                        using the integration method (l-m even)
!               coefmo : coefficient for the expression for r2 and r2d
!                        using the integration method (l-m odd)
!
!  real(knd) scalars and arrays and complex(knd) scalars and arrays
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) ak, amo2, an, argb, arn, bn, coef, coefme, coefmo, coefo, dec, &
                  etai, etal, etau, etais, etaism1, etcoef1, etcoef2, factor, &
                  factor2, f2exp, rm, rn, sargb, step, step0, step1, step2, ten, &
                  term, tf2, term1, term2, test, teste, testeo, test1, twom, &
                  twomi, x, x1, x2, x2m1
        real(knd) alpha(maxint), beta(maxint), p(maxint), pnorm(maxlp), &
                  wg(ngau), xg(ngau)
!
!  complex(knd) scalars and arrays
        complex(knd) arg, cc, coef1, coef2, coef4, darg, darg2, sbes1, sbes2, &
                     sneuna, sneunb, sneunc, sneu1, sneu2, sneu3, ynorm1, &
                     ynorm2, ynorm3
        complex(knd) pint1(maxint), pint2(maxint), pint3(maxint), &
                     pint4(maxint), rpint1(maxint), rpint2(maxint), &
                     sbesf(lim1max)
!
!  integer array
        integer ipnorm(maxlp)
!
        ten = 10.0e0_knd
        nfac = nex / 3
        if(2 * (nfac / 2) /= nfac) nfac = nfac - 1
        teste = ten ** (nfac)
        testeo = 1.0e0_knd / teste
        test1 = ten ** (nex - 3)
        ifac = nex - 30
        if(aimag(cc) > 55.0e0_knd) ifac = nex - 40
        factor = ten ** (-ifac)
        test = 1.0e0_knd / factor
        x = x1 + 1.0e0_knd
        x2 = x * x
        x2m1 = x1 * (x1 + 2.0e0_knd)
        dec = 10.0e0_knd ** (-ndec - 1)
        rm = real(m, knd)
        amo2 = 0.5e0_knd * rm
        lim = limint
        step0 = 1.0e0_knd / ngqs
        coefme = (rm * (x1 + 1.0e0_knd)) / x2m1
        coefmo = (rm * (x1 + 1.0e0_knd) ** 2 + x2m1) / ((x1 + 1.0e0_knd) * x2m1)
        if(x1 >= 0.1e0_knd) go to 10
        step1 = 1.0e0_knd / (ngqs * (1.0e0_knd - 3.0e0_knd * log10(x1) * (1.0e0_knd &
                     +3.0e0_knd * log10(rm + 10.0e0_knd))))
        if(step1 > step0) step1 = step0
        go to 20
10      step1 = step0
20      step2 = (1.0e0_knd - step1) / (ngqs - 1)
!
!  calculation of scaling factors for the associated Legendre functions
        pnorm(1) = 1.0e0_knd
        pnorm(2) = 1.0e0_knd
        ipnorm(1) = 0
        ipnorm(2) = 0
        if(m == 0) go to 50
          do 40 n = 1, m
          an = real(n + n)
          bn = real(n + n + 1)
          pnorm(1) = pnorm(1) * an
          pnorm(2) = pnorm(2) * bn
          iterm1 = int(log10(pnorm(1)))
          iterm2 = int(log10(pnorm(2)))
          pnorm(1) = pnorm(1) * 10.0e0_knd ** (-iterm1)
          pnorm(2) = pnorm(2) * 10.0e0_knd ** (-iterm2)
          ipnorm(1) = ipnorm(1) + iterm1
          ipnorm(2) = ipnorm(2) + iterm2
40        continue
50      twom = real(m + m)
        pnorm(3) = pnorm(1) * (twom + real(2)) / real(2)
        iterm3 = int(log10(pnorm(3)))
        pnorm(3) = pnorm(3) * 10.0e0_knd ** (-iterm3)
        ipnorm(3) = iterm3 + ipnorm(1)
        if(lnum < 4) go to 70
          do 60 il = 4, lnum, 2
          pnorm(il) = pnorm(il - 2) * (twom + il - 1) / (il - 1)
          pnorm(il + 1) = pnorm(il - 1) * (twom + il) / (il)
          iterm1 = log10(pnorm(il))
          iterm2 = log10(pnorm(il + 1))
          ipnorm(il) = ipnorm(il - 2) + iterm1
          ipnorm(il + 1) = ipnorm(il - 1) + iterm2
          pnorm(il) = pnorm(il) * 10.0e0_knd ** (-iterm1)
          pnorm(il + 1) = pnorm(il + 1) * 10.0e0_knd ** (-iterm2)
60        continue
70       continue
!
!  calculation of the coefficients in the recursion relation used
!  for the scaled associated Legendre functions
        alpha(1) = (twom + 1.0e0_knd) * pnorm(1) / pnorm(2)
        alpha(1) = alpha(1) * 10.0e0_knd ** (ipnorm(1) - ipnorm(2))
        beta(1) = 0.0e0_knd
        alpha(2) = (twom + 3.0e0_knd) * pnorm(2) / (pnorm(3) * 2.0e0_knd)
        alpha(2) = alpha(2) * 10.0e0_knd ** (ipnorm(2) - ipnorm(3))
        beta(2) = -(twom + 1.0e0_knd) / (twom + 2.0e0_knd)
          do 80 il = 3, lim + 2
          alpha(il) = alpha(il - 2) * (twom + il - 1) * (twom + il + il - 1)* &
          (il - 2) / ((il - 1) * (twom + il) * (twom + il + il - 5))
          beta(il) = -(twom + il - 1) / (twom + il)
80        continue
!
          do 90 il = 1, lim + 2, 2
          pint1(il) = (0.0e0_knd, 0.0e0_knd)
          pint2(il) = (0.0e0_knd, 0.0e0_knd)
          pint3(il + 1) = (0.0e0_knd, 0.0e0_knd)
          pint4(il + 1) = (0.0e0_knd, 0.0e0_knd)
90        continue
!
!  calculation of the scaling exponents for the spherical Bessel
!  functions required for the four types of integrals
        twomi = 1.0e0_knd
        if(m == 0) go to 110
          do 100 n = 1, m
          twomi = twomi * real(n + n - 1, knd) / real(n + n, knd)
100       continue
110     continue
        arg = cc * sqrt(x2m1)
        darg = (1.0e0_knd, 0.0e0_knd) / arg
        ynorm1 = -cos(arg) * darg
        ynorm2 = ynorm1 * darg - sin(arg) * darg
        normy = 0
          do 120 n = 3, m + 3
          rn = n + n - 3
          arn = n + n - 1
          ynorm3 = -ynorm1 + darg * rn * ynorm2
            if(abs(ynorm3) > teste) then
            ynorm3 = ynorm3 * testeo
            ynorm2 = ynorm2 * testeo
            normy = normy + nfac
            end if
          ynorm1 = ynorm2
          ynorm2 = ynorm3
120       continue
        iterm = int(log10(abs(ynorm3)))
        normy = normy + iterm
        norme = normy
!
!  Gaussian quadrature integration loops. first dividing integrand
!  into ngqs steps
          do 180 k = 1, ngqs
          ak = real(k)
            if(k == 1) then
            etal = 0.0e0_knd
            etau = step1
            step = step1
            else
            etal = step1 + (ak - 2.0e0_knd) * step2
            etau = etal + step2
            step = step2
            end if
          etcoef1 = (etau + etal) / 2.0e0_knd
          etcoef2 = (etau - etal) / 2.0e0_knd
          coef = step
!
!  Gaussian quadrature integration over each step
            do 170 i = 1, ngau
            etai = etcoef1 + xg(i) * etcoef2
            etais = etai * etai
            etaism1 = 1.0e0_knd - etais
            argb = x2m1 + etais
            term2 = x2m1 * etaism1 * etaism1 / argb
            term1 = sqrt(term2)
            sargb = sqrt(argb)
            coefo = 1.0e0_knd / sargb
            arg = cc * sargb
            darg = (1.0e0_knd, 0.0e0_knd) / arg
            sneu1 = -cos(arg) * darg
            sneu2 = term1 * darg * (sneu1 - sin(arg))
            norm = ifac
            lown = 3
            limb = 2 * int(abs(real(arg)) + abs(aimag(arg)))
            lim1 = 2 * limb + 20
            if(abs(aimag(arg)) < 2.0e0_knd .or. m < 3) go to 130
            sbesf(lim1) = arg / real(lim1 + lim1 + 1, knd)
              do n = lim1 - 1, 1,-1
              rn = real(n + n + 1, knd)
              sbesf(n) = 1.0e0_knd / (rn * darg - sbesf(n + 1))
              end do
            darg2 = darg * darg
            sneuna = sneu1
            sneunb = sneu2 / term1
            sbes1 = sbesf(1) * sin(arg) * darg
            if(limb > m + 2) limb = m + 2
              do n = 2, limb
              sbes2 = sbesf(n) * sbes1
              sneunc = (sbes2 * sneunb - darg2) / sbes1
              sbes1 = sbes2
                if(n /= limb) then
                sneuna = sneunb
                sneunb = sneunc
                end if
              end do
              nf2 = int(real(limb - 2, knd) * abs(log10(term1))) / nfac + 1
              f2exp = real(limb - 2, knd) / real(nf2, knd)
              factor2 = 1.0e0_knd
              tf2 = term1 ** (f2exp)
                do ll = 1, nf2
                factor2 = factor2 * tf2
                  if(factor2 > teste) then
                  factor2 = factor2 * testeo
                  norm = norm + nfac
                  end if
                  if(factor2 < testeo) then
                  factor2 = factor2 * teste
                  norm = norm - nfac
                  end if
                end do
            sneu1 = sneuna * factor2
            sneu2 = term1 * sneunb * factor2
            sneu3 = term2 * sneunc * factor2
              if(limb == m + 2) then
              go to 150
              else
              sneu1 = sneu2
              sneu2 = sneu3
              lown = limb + 2
              end if
130         continue
              do n = lown, m + 3
              rn = real(n + n - 3, knd)
              sneu3 = -term2 * sneu1 + term1 * darg * rn * sneu2
              if(n == m + 3) go to 150
                if(abs(sneu3) > teste) then
                sneu3 = sneu3 * testeo
                sneu2 = sneu2 * testeo
                norm = norm + nfac
                end if
                if(abs(sneu3) < testeo) then
                sneu3 = sneu3 * teste
                sneu2 = sneu2 * teste
                norm = norm - nfac
                end if
              sneu1 = sneu2
              sneu2 = sneu3
140           end do
150         continue
            iterm = int(log10(abs(sneu3)))
            term = ten ** (-iterm)
            sneu3 = sneu3 * term
            sneu2 = sneu2 * term
            sneu1 = sneu1 * term
            norm = norm + iterm
            iexp = -norme + norm
            if(iexp <= -nex + 10 - ifac) go to 160
              if(iexp > -nex + 10) then
              term = ten ** iexp
              sneu3 = sneu3 * term
              sneu2 = sneu2 * term
              sneu1 = sneu1 * term
              coef1 = coef * wg(i) * sneu1
              coef2 = (coef / term1) * coefo * sneu2 * wg(i)
              coef4 = (coef / term2) * coefo * coefo * etai * sneu3 * wg(i)
              iflag = 0
              else
              term = ten ** (iexp + ifac)
              sneu3 = sneu3 * term
              sneu2 = sneu2 * term
              sneu1 = sneu1 * term
              iflag = 1
              end if
            p(1) = twomi * factor
            p(2) = alpha(1) * etai * p(1)
              if(iflag == 0) then
              pint1(1) = pint1(1) + coef1 * p(1)
              pint2(1) = pint2(1) + coef2 * p(1)
              pint4(2) = pint4(2) + coef4 * p(2)
              end if
              do il = 2, limi, 2
              p(il + 1) = alpha(il) * etai * p(il) + beta(il) * p(il - 1)
              p(il + 2) = alpha(il + 1) * etai * p(il + 1) + beta(il + 1) * p(il)
                if(iflag == 0) then
                pint1(il + 1) = pint1(il + 1) + coef1 * p(il + 1)
                pint2(il + 1) = pint2(il + 1) + coef2 * p(il + 1)
                pint4(il + 2) = pint4(il + 2) + coef4 * p(il + 2)
                  if(abs(pint1(il + 1)) > test1 .or. abs(pint2(il + 1)) >  &
                   test1 .or. abs(pint4(il + 2)) > test1) then
                  limi = il
                  go to 160
                  end if
                end if
                if(abs(p(il + 2)) > test) then
                p(il + 1) = p(il + 1) * factor
                p(il + 2) = p(il + 2) * factor
                  if(iflag == 0) then
                  coef1 = coef1 * test
                  coef2 = coef2 * test
                  coef4 = coef4 * test
                  else
                  coef1 = coef * sneu1 * wg(i)
                  coef2 = (coef / term1) * coefo * sneu2 * wg(i)
                  coef4 = (coef / term2) * coefo * coefo * etai * sneu3 * wg(i)
                  iflag = 0
                  end if
                end if
              end do
160         continue
170         continue
180       continue
190     continue
          do 200 il = 1, limi - 1, 2
          pint3(il + 1) = (pint2(il + 2) - beta(il + 1) * pint2(il)) &
                      /alpha(il + 1)
200       continue
210     continue
if (debug) then
        write(40, 220) ngau, step1, limi
220     format(10x,'order of gauss quadrature =',i4,'. first step', &
               ' size =',f10.6,/,10x,'integrals for ',i5, &
               ' lowest order Legendre functions will be used for r2.')
end if
!
!  calculation of ratios of integrals for ease in compution of r2 and
!  r2d in subroutine r2int
        rpint1(1) = (0.0e0_knd, 0.0e0_knd)
        rpint1(2) = (0.0e0_knd, 0.0e0_knd)
        rpint2(1) = (0.0e0_knd, 0.0e0_knd)
        rpint2(2) = (0.0e0_knd, 0.0e0_knd)
          do 240 il = 3, limi, 2
          rpint1(il) = (pint1(il) * (twom + il - 1)) / (pint1(il - 2) * (il - 1))
          rpint2(il) = (pint2(il) * (twom + il - 1)) / (pint2(il - 2) * (il - 1))
          rpint1(il + 1) = (pint3(il + 1) * (twom + il)) / (pint3(il - 1) * (il))
          rpint2(il + 1) = (pint4(il + 1) * (twom + il)) / (pint4(il - 1) * (il))
240       continue
        limi = limi - 1
        return
        end subroutine
!
!
        subroutine sphbes (cc, x, limj, maxj, maxlp, ndec, nex, sbesf, sbesdf, &
                           sbesn, ibese, sbesdr)
!
!  purpose:     To calculate ratios of successive first kind spherical
!               Bessel functions of the same parity for given c and x.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Bessel functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    cc     : complex c
!               x      : argument of spherical Bessel functions
!               limj   : the number of spherical Bessel function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxj   : dimension of sbesf vector
!               maxlp  : the number of scale factors that are
!                        calculated
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent available for real(knd)
!
!     output:   sbesf  : ratios of successive first kind spherical
!                        Bessel functions of the same parity
!               sbesdf : ratios of first derivatives of successive
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : characteristics for the spherical
!                        Bessel functions
!               ibese  : exponents for the spherical
!                        Bessel functions
!               sbesdr : ratios of first derivatives of spherical Bessel
!                        functions to the corresponding spherical functions
!
        use param
!
!  real(knd) scalars
        real(knd) adj, ci, cm, cr, rn, ten, x
!
!  complex(knd) scalars and arrays
        complex(knd) cc, cx, stemp0, stemp1
        complex(knd) sbesdf(maxj), sbesdr(maxlp), sbesf(maxj), &
                     sbesn(maxlp)
!
!  integer array
        integer ibese(maxlp)
!
        ten = 10.0e0_knd
        cx = cc * x
        ci = aimag(cc)
        cr = real(cc)
        cm = abs(cc)
        lim1 = 2 * cm * x + 20
!
!  compute first kind Bessel function ratios
!        sbesf(k)= j(n=k,c*x)/j(n=k-1,c*x)
!        sbesn(k)= (j(n=k-1),c*x))*10.0e0_knd**(-ibese(k))
!
        if(int(abs(cx)) < limj .or. aimag(cc) /= 0.0e0_knd) go to 20
!
!  for c*x >= limj and aimag(cc) = 0, use forward recursion to
!  get fcn. ratios:
!       j(n+1,c*x)/j(n,c*x)=(2*n+1)/(c*x)-1/(j(n,c*x)/j(n-1,c*x))
!
        stemp0 = sin(cx) / cx
        sbesf(1) = (stemp0 / cx - cos(cx) / cx) / stemp0
          do 10 n = 1, limj - 1
          rn = n + n + 1
          sbesf(n + 1) = (rn / cx) - (1.0e0_knd / sbesf(n))
10        continue
        go to 60
20      continue
!
!  for c*x < lim or aimag(cc) unequal to zero, use backward recursion
!  to get fcn. ratios:
!       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
!
        stemp0 = 0.0e0_knd
        if(lim1 <= limj) go to 40
          do 30 n = lim1, limj,-1
          rn = n + n + 1
          stemp1 = 1.0e0_knd / (rn / cx - stemp0)
          stemp0 = stemp1
30        continue
40      sbesf(limj) = stemp0
          do 50 n = limj - 1, 1,-1
          rn = n + n + 1
          sbesf(n) = 1.0e0_knd / (rn / cx - sbesf(n + 1))
50        continue
60      continue
!
!  for all c*x, calculate the amplitude and sign scale
!  factors by forward operation on the Bessel function
!  ratios.
          if(0.43e0_knd * x * ci > nex / 3) then
          adj = x * ci / log(10.0e0_knd)
          iadj = int(adj)
          sbesn(1) = (0.0e0_knd, 1.0e0_knd)* &
                     exp((0.0e0_knd,-1.0e0_knd) * cr * x) / (2.0e0_knd * cx)
          sbesn(1) = sbesn(1) * (10.0e0_knd ** (adj - iadj))
          iterm = int(log10(abs(sbesn(1))))
          sbesn(1) = sbesn(1) * (10.0e0_knd ** (-iterm))
          ibese(1) = iadj + iterm
          sbesn(2) = sbesn(1) / cx + (0.0e0_knd, 1.0e0_knd) * sbesn(1)
          iterm = log10(abs(sbesn(2)))
          sbesn(2) = sbesn(2) * (10.0e0_knd ** (-iterm))
          ibese(2) = ibese(1) + iterm
          else
          stemp0 = sin(cx) / cx
          stemp1 = stemp0 / cx - cos(cx) / cx
          ibese(1) = log10(abs(stemp0))
          sbesn(1) = stemp0 * ten ** (-ibese(1))
          if(abs(sin(cx)) < 0.5e0_knd .and. abs(cx) > 1.0e0_knd) &
             go to 70
          sbesn(2) = sbesn(1) * sbesf(1)
          ibese(2) = log10(abs(sbesn(2)))
          sbesn(2) = sbesn(2) * ten ** (-ibese(2))
          ibese(2) = ibese(2) + ibese(1)
          go to 80
70        ibese(2) = log10(abs(stemp1))
          sbesn(2) = stemp1 * ten ** (-ibese(2))
          sbesf(1) = stemp1 / stemp0
80        continue
          end if
           do 90 n = 3, maxlp
          sbesn(n) = sbesn(n - 1) * sbesf(n - 1)
          ibese(n) = log10(abs(sbesn(n)))
          sbesn(n) = sbesn(n) * ten ** (-ibese(n))
          ibese(n) = ibese(n) + ibese(n - 1)
90        continue
!
!  calculate the ratios of the first derivatives of successive
!  Bessel functions using corresponding function ratios
          do 100 n = 1, limj
          rn = n - 1
          sbesdf(n) = (cx - (rn + 2.0e0_knd) * sbesf(n)) / (rn - cx * sbesf(n))
100       continue
!
!  calculate the ratios of the first derivative to the corresponding
!  spherical Bessel function
          do 110 n = 1, maxlp
          rn = n - 1
          sbesdr(n) = (rn / cx) - sbesf(n)
110       continue
!
!  calculate the ratios of successive functions and derivatives
!  of the same parity
          do 120 n = limj, 2,-1
          sbesf(n) = sbesf(n - 1) * sbesf(n)
          sbesdf(n) = sbesdf(n - 1) * sbesdf(n)
120       continue
        return
        end subroutine
!
!
        subroutine sphneu (cc, x, limn, maxn, maxlp, ndec, nex, limbes, sneuf, &
                           sneun, ineue, sneudf, sneudr)
!
!  purpose:     to calculate ratios of spherical Neumann functions
!               and ratios of their first derivatives for given c and x.
!               to calculate the Neumann function characteristics
!               and exponents. to calculate ratios of the first
!               derivatives of the corresponding Neumann functions.
!
!  parameters:
!
!     input:    cc     : complex c
!               x      : argument of spherical Neumann functions
!               limn   : the number of spherical Neumann function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxn   : dimension of sneuf and sneudf arrays
!               maxlp  : the number of values of scale factors
!                        that are calculated
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent available for real(knd)
!               limbes : dimension of spherical Bessel function ratios
!                        calculated for use in calculating Neumann
!                        functions
!
!     output:   sneuf  : ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : characteristic for the spherical
!                        Neumann functions
!               ineue  : exponent for the spherical
!                        Neumann functions
!               sneudf : ratios of first derivatives of successive
!                        spherical Neumann functions of the same parity
!               sneudr : ratios of first derivatives of spherical
!                        Neumann functions to the corresponding
!                        function
!
        use param
!
!  real(knd) scalars
        real(knd) adj, ci, cr, rn, rnn, test, x
!
!  complex scalars and arrays
        complex(knd) cc, cx, cx2, sbes1, sbes2, stemp0, stemp1
        complex(knd) sneudf(maxn), sneudr(maxlp), sneuf(maxn), &
                     sneun(maxn), sbesf(limbes)
!
!  integer arrays
        dimension ineue(maxn)
!
!  compute first kind ratios of Neumann functions and ratios
!  of their first derivatives
!
!        sneuf(k)=y(n=k,c*x)/y(n=k-2,c*x)
!        sneun(k)=(y(n=k-1),c*x)*10.e0_knd**(-ineue(k))
!        sneudf(k)=y'(n=k,c*x)/y'(n=k-2,c*x)
!        sneudr(k)=(y'(n=k-1),c*x)/y(n=k-1),c*x))
!
!  calculate j ratios below turning point by backward recursion
!
!       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
!
        cx = cc * x
        cx2 = cx * cx
        cr = real(cc)
        ci = aimag(cc)
        test = 1.0e0_knd / abs(cx)
        itest = -int(log10(abs(cx)))
        limb = 2 * int(abs(real(cx)) + abs(aimag(cx)))
        lim1 = 2 * limb + 20
        if(limb < 2) limb = 2
        sbesf(lim1) = cx / (lim1 + lim1 + 1)
          do 10 n = lim1 - 1, 1,-1
          rn = real(n + n + 1, knd)
          sbesf(n) = 1.0e0_knd / (rn / cx - sbesf(n + 1))
10        continue
!
!  use relation with j's to compute y's from order zero
!  to the turning point. compute derivative ratios
!  at same time.
!
        if(0.43e0_knd * x * ci <= nex / 3) go to 30
        adj = x * ci / log(10.0e0_knd)
        iadj = int(adj)
        sbes1 = (0.0e0_knd, 1.0e0_knd)* &
               exp((0.0e0_knd,-1.0e0_knd) * cr * x) / (2.0e0_knd * cx)
        sbes1 = sbes1 * (10.0e0_knd ** (adj - iadj))
        iterm = int(log10(abs(sbes1)))
        ibes1 = iadj + iterm
        sbes1 = sbes1 * (10.0e0_knd ** (-iterm))
        sneun(1) = (0.0e0_knd, 1.0e0_knd) * sbes1
        ineue(1) = ibes1
        sbes2 = sbes1 / cx + (0.0e0_knd, 1.0e0_knd) * sbes1
        sneuf(1) = sbes2 / sbes1
        iterm = int(log10(abs(sbes2)))
        ibes2 = ibes1 + iterm
        sbes2 = sbes2 * (10.0e0_knd ** (-iterm))
        sneun(2) = (0.0e0_knd, 1.0e0_knd) * sbes2
        ineue(2) = ibes2
        sneudf(1) = -(cx - 2.0e0_knd * sneuf(1)) / (cx * sneuf(1))
        if(limb > limn) limb = limn
        j = 1
        sbes1 = sbes2
        ibes1 = ibes2
          do n = 2, limb
          j = j + 1
          rn = real(n, knd)
          rnn = real(n - 1, knd)
          sbes2 = sbesf(n) * sbes1
          iterm = int(log10(abs(sbes2)))
          sbes2 = sbes2 * (10.0e0 ** (-iterm))
          ibes2 = ibes1 + iterm
          sneun(n + 1) = sbesf(n) * sneun(n)
          if(ibes1 + ineue(n) < ndec + 10) &
          sneun(n + 1) = sneun(n + 1) - (10.0e0_knd ** (-ibes1 - ineue(n)))/ &
                     (cx2 * sbes1)
          sneuf(n) = sneun(n + 1) / sneun(n)
          iterm = int(log10(abs(sneun(n + 1))))
          sneun(n + 1) = sneun(n + 1) * (10.0e0_knd ** (-iterm))
          ineue(n + 1) = ineue(n) + iterm
          sneudf(n) = (cx - (rn + 1.0e0_knd) * sneuf(n)) / (rnn - cx * sneuf(n))
          if(ibes1 + ibes2 < itest) go to 20
          sbes1 = sbes2
          ibes1 = ibes2
          end do
20      limb = max(j, 2)
        go to 70
30      continue
        stemp0 = -cos(cx) / cx
        stemp1 = (stemp0 - sin(cx)) / cx
        sneuf(1) = stemp1 / stemp0
        sneun(1) = stemp0
        sneun(2) = stemp1
        sbes1 = sbesf(1) * sin(cx) / cx
        sneudf(1) = -(cx - 2.0e0_knd * sneuf(1)) / (cx * sneuf(1))
        if(limb > limn) limb = limn
        j = 1
          do 40 n = 2, limb
          j = j + 1
          rn = real(n, knd)
          rnn = real(n - 1, knd)
          sbes2 = sbesf(n) * sbes1
          sneun(n + 1) = sbesf(n) * sneun(n) - 1.0e0_knd / (cx2 * sbes1)
          sneuf(n) = sneun(n + 1) / sneun(n)
          sneudf(n) = (cx - (rn + 1.0e0_knd) * sneuf(n)) / (rnn - cx * sneuf(n))
          if(abs(sbes1 + sbes2) < test) go to 50
          sbes1 = sbes2
40        continue
50      limb = max(j, 2)
!
!  calculate characteristics and exponents for the Neumann functions
!  up to and including the turning point
        ineue(1) = int(log10(abs(stemp0)))
        sneun(1) = stemp0 * 10.0e0_knd ** (-ineue(1))
        ineue(2) = int(log10(abs(stemp1)))
        sneun(2) = stemp1 * 10.0e0_knd ** (-ineue(2))
          do 60 n = 3, limb
          ineue(n) = int(log10(abs(sneun(n))))
          sneun(n) = sneun(n) * 10.0e0_knd ** (-ineue(n))
60        continue
70      continue
!
!  use forward recursion from breakpoint to compute function ratios
!
!       y(n+1,c*x)/y(n,c*x)=(2*n+1)/(c*x)-1/(y(n,c*x)/y(n-1,c*x))
!
!  compute derivative ratios at same time using function ratios.
        if(limb == limn) go to 90
          do 80 n = limb + 1, limn
          rn = real(n - 1, knd)
          rnn = real(n + n - 1, knd)
          sneuf(n) = rnn / cx - 1.0e0_knd / sneuf(n - 1)
          sneudf(n) = (cx - (rn + 2.0e0_knd) * sneuf(n)) / (rn - cx * sneuf(n))
80        continue
90      continue
          sneuf(limn + 1) = (0.0e0_knd, 0.0e0_knd)
          sneuf(limn + 2) = (0.0e0_knd, 0.0e0_knd)
          sneudf(limn + 1) = (0.0e0_knd, 0.0e0_knd)
          sneudf(limn + 2) = (0.0e0_knd, 0.0e0_knd)
!
!  calculate the characteristics and exponents for Neumann
!  functions beyond the turning point by forward operation
!  on the Neumann function ratios:
        if(limb + 1 > maxlp) go to 110
        nlimit = maxlp
        if(maxlp > limn + 1) nlimit = limn + 1
          do 100 n = limb + 1, nlimit
          sneun(n) = sneun(n - 1) * sneuf(n - 1)
          ineue(n) = int(log10(abs(sneun(n))))
          sneun(n) = sneun(n) * 10.0e0_knd ** (-ineue(n))
          ineue(n) = ineue(n) + ineue(n - 1)
100       continue
110     continue
!
!  calculate the ratios of the first derivatives to the corresponding
!  spherical Neumann functions
          do 120 n = 1, maxlp
          rn = real(n - 1, knd)
          sneudr(n) = (rn / cx) - sneuf(n)
120       continue
!
!  calculate the ratios of successive functions and derivatives
!  of the same parity
          do 130 n = limn, 2,-1
          sneuf(n) = sneuf(n - 1) * sneuf(n)
          sneudf(n) = sneudf(n - 1) * sneudf(n)
130       continue
        return
        end subroutine
        end module complex_prolate_swf

