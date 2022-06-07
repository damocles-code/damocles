module random_routines

    use input

    implicit none

contains

 function normal(mean, sigma)

  ! mean  : mean of distribution
  ! sigma : number of standard deviations

  implicit none

  integer, parameter:: b8 = selected_real_kind(14)
  !real(b8), parameter :: pi = 3.141592653589793239_b8
  real(b8) normal
  real rand_num
  real(b8) tmp
  real(b8) mean
  real(b8) sigma
  real(b8) fac
  real(b8) gsave
  real(b8) rsq
  real(b8) r1
  real(b8) r2
  integer flag
  save flag
  save gsave
  data flag /0/

  if (flag.eq.0) then
    rsq=2.0_b8

    do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
      call random_number(rand_num)
      r1=2.0_b8*rand_num-1.0_b8
      call random_number(rand_num)
      r2=2.0_b8*rand_num-1.0_b8
      rsq=r1*r1+r2*r2
    enddo

    fac=sqrt(-2.0_b8*log(rsq)/rsq)
    gsave=r1*fac
    tmp=r2*fac
    flag=1
  else
    tmp=gsave
    flag=0
  endif

  normal=tmp*sigma+mean

  return
  endfunction normal

end module
