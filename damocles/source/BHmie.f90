! ***********************************************************************
! this subroutine is a modified version of the bh mie scattering routine
! subroutine taken from mocassin (barbara ercolano).
! <gsca> added from b. draine implementation by barbara ercolano.
! -----------------------------------------------------------------------
! __________________________________________________________________
!
!  subroutine bhmie calculates amplitude scattering matrix elements
!  & efficiencies for extinction, total scattering and bacscattering,
!  for a girven size parameter and relative refractive index
! __________________________________________________________________
      
subroutine bhmie (x,refrel,qext,qsca,ggsca)
    implicit none
        
    complex, intent(in) :: refrel
    real, intent(in)    :: x
        
    real, intent(out)   :: qext, qsca, ggsca
        
    real(kind = 8), dimension(100) ::  amu, theta,pii,tau,pi0,pi1
        
    double complex              :: d(100000),y,xi,xi0,xi1,an,bn,s1(200),s2(200),&
        &an1,bn1
    real(kind = 8)                 ::  psi0,psi1,psi,dn,dx,xstop,ymod,dang,chi0,&
        &chi1,apsi0,apsi1,apsi,fn,chi,p,t

    double complex :: dpcx

    real(kind = 8) ::  realpart
    real(kind = 8) :: imagpart
        
    integer              :: nang,j,nstop,nmx,i,n,nn,rn,jj


    realpart(dpcx)=real(dpcx)
    imagpart(dpcx)=imag(dpcx)


    nang = 2
    dx=x
    y=x*refrel

    ! series terminated after nstop terms
    ! ___________________________________________________________________
    xstop=x+4.*x**.3333 +2.0
    nstop=int(xstop)
    ymod=abs(y)
    nmx=int(max(xstop,ymod)) + 15
    dang=1.570796327/float(nang-1)
    do j = 1,nang
        theta(j)= (float(j)-1.)*dang
        amu(j)=cos(theta(j))

    end do


    ! __________________________________________________________________
    ! logarithmic derivative d(j) calculated by downward recurrence
    ! beginning with initial value 0.0 + i*0.0 at j = nmx
    !__________________________________________________________________
    d(nmx)=cmplx(0.0,0.0)
    nn=nmx-1
    do n=1,nn
        rn=nmx-n+1
        d(nmx-n)=(rn/y)-(1./(d(nmx-n+1)+rn/y))

    end do
    do j=1,nang
        pi0(j)=0.0
        pi1(j)=1.0
    end do
    nn=2*nang-1
      
    do j=1,nn
        s1(j)=cmplx(0.0,0.0)
        s2(j)=cmplx(0.0,0.0)
    end do
    ! __________________________________________________________________
    ! riccati bessel functions with real argument x calculated by upward
    ! recurrence
    ! __________________________________________________________________
    !
    psi0=cos(dx)
    psi1=sin(dx)
    chi0=-sin(x)
    chi1=cos(x)
    apsi0=psi0
    apsi1=psi1
    xi0=cmplx(apsi0,-chi0)
    xi1=cmplx(apsi1,-chi1)

    qsca=0.0
    ggsca = 0.0
    n=1
    do
        dn=n
        rn=n

        fn=(2.*rn+1.)/(rn*(rn+1.))
        psi=(2.*dn-1.)*psi1/dx-psi0

        apsi=psi
        chi=(2.*rn-1.)*chi1/x -  chi0
        xi = cmplx(apsi,-chi)

        if (n>1) then
            an1 = an
            bn1 = bn
        end if

        an=(d(n)/refrel+rn/x)*apsi - apsi1

        an=an/((d(n)/refrel+rn/x)*xi - xi1)

        bn=(refrel *d(n)+rn/x)*apsi - apsi1

        bn=bn/((refrel*d(n)+rn/x)*xi - xi1)

        qsca=qsca+(2.*rn+1.)*(abs(an)*abs(an)+abs(bn)*abs(bn))
         
        ggsca=ggsca+real(((2.*rn+1.)/(rn*(rn+1.))))*&
            &((realpart(an)*realpart(bn)+imagpart(an)*imagpart(bn)))

        if(n>1)then
            ggsca=ggsca+real(((rn-1.)*(rn+1.)/rn)*&
                &(realpart(an1)*realpart(an)+imagpart(an1)*imagpart(an)+&
                &realpart(bn1)*realpart(bn)+imagpart(bn1)*imagpart(bn)))
        endif

        do j=1,nang

            jj=2*nang-j
            pii(j)=pi1(j)
            tau(j)=rn*amu(j)*pii(j) - (rn+1.)*pi0(j)

            p=(-1.)**(n-1)
            s1(j)=s1(j)+fn*(an*pii(j)+bn*tau(j))
            t=(-1.)**n
            s2(j)=s2(j) + fn*(an*tau(j)+bn*pii(j))

            if (j == jj) exit
            s1(jj)=s1(jj) + fn*(an*pii(j)*p + bn*tau(j)*t)
            s2(jj)=s2(jj) + fn*(an*tau(j)*t + bn*pii(j)*p)

        end do
        psi0=psi1
        psi1=psi
        apsi1=psi1
        chi0=chi1
        chi1=chi
        xi1=cmplx(apsi1,-chi1)

        n=n+1
        rn=n
        do j=1,nang
            pi1(j)=((2.*rn-1.)/(rn-1.))*amu(j)*pii(j)
            pi1(j)=pi1(j) - rn*pi0(j)/(rn-1.)
            pi0(j) = pii(j)
        end do

        if (n-1-nstop>=0) exit
    end do

    ggsca=(2.*ggsca/qsca)
    qsca=(2./(x*x))*qsca
    qext=(4./(x*x))*realpart(s1(1))

end subroutine bhmie
