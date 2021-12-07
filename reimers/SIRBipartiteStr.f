c******************************************************************************
c program to find distributions of infection probabilities and in bipartite
c   structured problem of PIIs and PICs, version f. UNCORRELATED DEGREE DISTRns
c
c   04/06/20  edited from SIR-bipartite.f
c
c	--need populations to represent \pi^{(na)}(\tm) and \pi^{(an)}(\tv)
c	--need separate RNGs for all PC^{(na)} and all PD^{(an)}, which require
c --specification of degree ramges KCR,KDR, degrees KC, KD, and probabilities
c   PC and PD (to  avoid waste of storage, hold separate arrays for the sectors
c   with broad distributions?) E.g. KCR(1,n,a) < o => KCR(2,n,a) = ILL
c       where ILL is the index in the Long Degree Sequence List KL(*,ILL)
c       use NL to denote the number of long sequences
c
c   NPIC, NPII  number of PIC and PII types
c   KCR(2,NPII,NPIC)   with KCR(*,n,a) containing min degree and number 
c                                of degrees for a -> n statistics
c   KDR(2,NPIC,NPII)   with KDR(*,a,n) containing min degree and number 
c                                of degrees for n -> a statistics
c   KC(*,n,a),PC(*,n,a) list of degrees/probabilities for a -> n statistics
c   KD(*,a,n),PD(*,a,n) list of degrees/probabilities for a -> n statistics 
c   KL(*, l) , PL(*,l)  list of degrees/probabilities for l-th long sequence
c
c   KMax, KLMax   maximum number of degrees for short/long sequences
c   Np            populations size
c   ITW,ITE,ITM   number of warmup, equilibration and measurement sweeps
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c   current version go simpler taylored to the specific example
c   needs differend pmfs P(1/2),P(1),1+P(1),P(5/8),P(5/2),Pow(3,3,100), and
c   \delta_{k,1} as well as the corresponding Q(k)=k P(K)/c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c------------------------------------------------------------------------------
c******************************************************************************
        module StrSIRBPDefs
c******************************************************************************
        implicit none
        integer (kind=4):: Np,ITM,NPII,NPIC,NDat,NL,KMax,KLmax
        parameter (Np=5000,ITM=1000,NPII=3,NPIC=3,NDat=Np*ITM, 
     &     NL=1,KMax=20,KLMax=98) ! << initialize
c        integer (kind=4):: KCR(2,NPII,NPIC),KDR(2,NPIC,NPII)
c        integer (kind=4):: KC(KMax,NPII,NPIC),KD(KMax,NPIC,NPII)
c        real    (kind=4):: PC(KMax,NPII,NPIC),PD(KMax,NPIC,NPII)
        real    (kind=4):: tM(NP,NPII,NPIC),tV(NP,NPIC,NPII)   
        integer	(kind=4):: KL(KLMax),K1(KMax),K1P1(KMax)
        real	  (kind=4):: Pw(KLMax),P1(KMax),P58(KMax),P05(KMax)
     &        ,P52(KMax),P1P1(KMax)
        real	  (kind=4):: Qw(KLMax),Q1(KMax),Q58(KMax),Q05(KMax)
     &        ,Q52(KMax),Q1P1(KMax)
        real  	(kind=4):: rhv(NDat),rh1(NDat),rh2(NDat),rh3(Ndat)
        real  	(kind=4):: rn(3),sa(3),C,tb,mv,m1,m2,m3
        integer	(kind=4):: ioffs,iseed=-521245
        data rn /1.0,1.0,1.0/        ! fraction of susceptible PIIs by group
        data sa /1.0,1.0,1.0/        ! fraction of open PICs by class

        end module  StrSIRBPDefs
*******************************************************************************
        program SIRBStrU
c******************************************************************************
        use StrSIRBPDefs
        implicit none
        integer (kind=4):: nbins
        parameter (nbins=1001)
        integer (kind=4):: ITW,ITE,it,mode, i
        real	  (kind=4):: tb1,tb2,dtb

c-- open output file
        open(unit=10,file='../Data/SIRBStrPopDyn-P-tb011')
        open(unit=11,file='../Data/SIRBStrPopDyn-Av-P-tb011')
c
c some parameter and constant initializations
c
        ITW= 1000	! warmup
	      ITE= 1000	! equilibration
c
c ini value for p loops (transmission probability)
c
	      tb1=0.11
	      tb2=0.2
	      dtb=-0.005

	      tb= tb1

	      write(6,*) 'compute coordinations'

c---  define the various pmfs (not necessarily normalized)
        do i=1,KMax
          K1(i)=i-1
          K1P1(i) = i
        enddo

        P1(1)=exp(-1.)    !Poi(1)
        do i=1,KMax-1
         P1(i+1)=P1(i)/i
        enddo
        Q1=K1*P1
        
        P1P1=P1     ! same probabilities but K1P1 instead of K1
        Q1P1=K1P1*P1P1
        
        C=0.5
        P05(1) = exp(-C)
        do i=1,KMax-1
          P05(i+1)=P05(i)*C/i
        enddo
        Q05=K1*P05
  
        C=5./8.
        P58(1) = exp(-C)
        do i=1,KMax-1
          P58(i+1)=P58(i)*C/i
        enddo
        Q58=K1*P58

        C=5./2.
        P52(1) = exp(-C)
        do i=1,KMax-1
          P52(i+1)=P52(i)*C/i
        enddo
        Q52=K1*P52

c - power law        
c        do i=1,KLMax
c        KL(i) =i+2
c        Pw(i) =1./KL(i)**3
c        enddo
c        Pw = Pw/sum(Pw)
c        Qw = KL*Pw

c - change to P(5)
        C=5.0
        Pw(1) = exp(-C)
        do i=1,KLMax-1
          KL(i) = i-1
          Pw(i+1)=Pw(i)*C/i
        enddo
        Qw=KL*Pw

c-- actual initialization of RNGs 
        call iran_d(P1,KMax)
        call iran_d(Q1,KMax)
        call iran_d(P1P1,KMax)
        call iran_d(Q1P1,KMax)
        call iran_d(P05,KMax)
        call iran_d(Q05,KMax)
        call iran_d(P58,KMax)
        call iran_d(Q58,KMax)
        call iran_d(P52,KMax)
        call iran_d(Q52,KMax)
        call iran_d(Pw,KLMax)
        call iran_d(Qw,KLMax)

	      write(6,*) 'initialize population'
	      call ini_pop

c  warmup 
	      mode = 0
	      write(6,*) 'warmup '
	      do it=1,ITW
	        call sweep(mode)
	      enddo
	      write(6,*) 'warmup  finished   ',it-1,'  sweeps'

c        do while(tb .gt. tb2)
c        tb=max(tb+dtb,0.0)
c  equilibration 
	      mode = 0
	      do it=1,ITE
	        call sweep(mode)
	      enddo
	      write(6,*) 'equilibration  finished   ',it-1,'  sweeps'

c  measurement

	      mode = 1
	      ioffs = 0
	      do it=1,ITM
	        call sweep(mode)
	        ioffs =ioffs + Np
	      enddo
	      write(6,*) 'measurement  finished   ',it-1,'  sweeps'

         call data_to_hist(NDat,Nbins,tb,rh1,rh2,rh3,rhv)
         write(10,*)

         m1=sum(rh1)/NDat
         m2=sum(rh2)/NDat
         m3=sum(rh3)/NDat
         mv=sum(rhv)/NDat
	       write(6,100) tb,m1,m2,m3,mv
	       write(11,100) tb,m1,m2,m3,mv

c         enddo           ! do while(tb .gt. tb2)
	       close(unit=10)
	       close(unit=11)
100	     format(10e16.7)

         stop
c******************************************************************************
        end program SIRBStrU
c******************************************************************************
c	popdyn related subroutines
c******************************************************************************
c------------------------------------------------------------------------------
        subroutine ini_pop
c------------------------------------------------------------------------------
c initialization of population dynamics

        use StrSIRBPDefs
        implicit none
	      real	(kind=4):: ran0
	      integer	(kind=4):: ip,ia,in
c------------------------------------------------------------------------------
c- first population randomly generated

	      do ip=1,Np
          do ia=1,NPIC
            do in=1,NPII
              tm(ip,in,ia) = ran0(iseed)
              tv(ip,ia,in) = ran0(iseed)
            enddo
          enddo
        enddo	!ip
	      return
	      end
c------------------------------------------------------------------------------
	subroutine sweep(mode)
c------------------------------------------------------------------------------
c  	one sweep of updates: directed bipartite
c		no difference between marginal and cavity variabls
c------------------------------------------------------------------------------
	      use StrSIRBPDefs
        implicit none

        real	(kind=4):: fran_d,ran0,omg,omt,cm
        integer	(kind=4):: ip,ik,km,kn,ko,mode
c------------------------------------------------------------------------------
        do ip=1,Np

c--- update population of \tilde m 
c-- tm11
        km = fran_d(K1,Q58,KMax,iseed)
        kn = fran_d(K1,P58,KMax,iseed)
        ko = fran_d(K1,P52,KMax,iseed)
        omg = 1.0
         do ik=1,km-1
           omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,1))
         enddo
         do ik=1,kn
           omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,2))
         enddo
         do ik=1,ko
           omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,3))
         enddo
         tm(ip,1,1) = sa(1)*(1.0-omg)

c-- tm21
        km = fran_d(K1,P58,KMax,iseed)
        kn = fran_d(K1,Q58,KMax,iseed)
        ko = fran_d(K1,P52,KMax,iseed)
        omg = 1.0
         do ik=1,km
           omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,1))
         enddo
         do ik=1,kn-1
           omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,2))
         enddo
         do ik=1,ko
           omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,3))
         enddo
         tm(ip,2,1) =  sa(1)*(1.0-omg)

c-- tm31
        km = fran_d(K1,P58,KMax,iseed)
        kn = fran_d(K1,P58,KMax,iseed)
        ko = fran_d(K1,Q52,KMax,iseed)
        omg = 1.0
        do ik=1,km
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,1))
        enddo
        do ik=1,kn
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,2))
        enddo
        do ik=1,ko-1
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),1,3))
        enddo
        tm(ip,3,1) =  sa(1)*(1.0-omg)

c-- tm22
        kn = fran_d(K1,Q05,KMax,iseed)
        ko = fran_d(K1,P1,KMax,iseed)
        omg = 1.0
        do ik=1,kn-1
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),2,2))
        enddo
        do ik=1,ko
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),2,3))
        enddo
        tm(ip,2,2) = sa(2)*(1.0-omg)

c-- tm32
        kn = fran_d(K1,P05,KMax,iseed)
        ko = fran_d(K1,Q1,KMax,iseed)
        omg = 1.0
        do ik=1,kn
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),2,2))
        enddo
        do ik=1,ko-1
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),2,3))
        enddo
        tm(ip,3,2) = sa(2)*(1.0-omg)

c-- tm33
        ko = fran_d(KL,Qw,KLMax,iseed)
        omg = 1.0
        do ik=1,ko-1
          omg = omg*(1.0-tb*tv(1+int(Np*ran0(iseed)),3,3))
        enddo
        tm(ip,3,3) = sa(3)*(1.0-omg)

c -- tv11 NOTE d11=1 non-random and tv21=0 and tv31=0
        tv(ip,1,1) = 0.0

c -- tv12
        km = fran_d(K1,Q1,KMax,iseed)
        kn = fran_d(K1,P1,KMax,iseed)
        omg = 1.0
        do ik=1,km-1
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),2,1))
        enddo
        do ik=1,kn
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),2,2))
        enddo
        tv(ip,1,2) =rn(2)*(1.0-omg)
        
c -- tv22
        km = fran_d(K1,P1,KMax,iseed)
        kn = fran_d(K1,Q1,KMax,iseed)
        omg = 1.0
        do ik=1,km
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),2,1))
        enddo
        do ik=1,kn-1
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),2,2))
        enddo
        tv(ip,2,2) =rn(2)*(1.0-omg)
        
c -- tv13
        km = fran_d(K1P1,Q1P1,KMax,iseed)
        kn = fran_d(K1,P1,KMax,iseed)
        ko = fran_d(K1,P1,KMax,iseed)
        omg = 1.0
        do ik=1,km-1
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,1))
        enddo
        do ik=1,kn
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,2))
        enddo
        do ik=1,ko
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,3))
        enddo
        tv(ip,1,3) =rn(3)*(1.0-omg)

c -- tv23
        km = fran_d(K1P1,P1P1,KMax,iseed)
        kn = fran_d(K1,Q1,KMax,iseed)
        ko = fran_d(K1,P1,KMax,iseed)
        omg = 1.0
        do ik=1,km
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,1))
        enddo
        do ik=1,kn-1
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,2))
        enddo
        do ik=1,ko
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,3))
        enddo
        tv(ip,2,3) =rn(3)*(1.0-omg)

c -- tv33
        km = fran_d(K1P1,P1P1,KMax,iseed)
        kn = fran_d(K1,P1,KMax,iseed)
        ko = fran_d(K1,Q1,KMax,iseed)
        omg = 1.0
        do ik=1,km
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,1))
        enddo
        do ik=1,kn
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,2))
        enddo
        do ik=1,ko-1
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,3))
        enddo
        tv(ip,3,3) =rn(3)*(1.0-omg)
        
	      if(mode .eq. 1) then
	 
c -- rh1
        omg = tm(1+int(Np*ran0(iseed)),1,1)
        rh1(ioffs+ip) = rn(1)*omg
c -- rh2
        km = fran_d(K1,P1,KMax,iseed)
        kn = fran_d(K1,P1,KMax,iseed)
        omg = 1.0
        do ik=1,km
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),2,1))
        enddo
        do ik=1,kn
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),2,2))
        enddo
        rh2(ioffs+ip) = rn(2)*(1.0-omg)
c -- rh3
        km = fran_d(K1p1,P1P1,KMax,iseed)
        kn = fran_d(K1,P1,KMax,iseed)
        ko = fran_d(K1,P1,KMax,iseed)
        omg = 1.0
        do ik=1,km
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,1))
        enddo
        do ik=1,kn
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,2))
        enddo
        do ik=1,ko
          omg = omg*(1.0-tm(1+int(Np*ran0(iseed)),3,3))
        enddo
        rh3(ioffs+ip) = rn(3)*(1.0-omg)
c -- rhv	      
        if(mod(ip,4) .lt. 2) then
         rhv(ioffs+ip) = rh3(ioffs+ip)
        else if(mod(ip,4) .eq. 2) then
         rhv(ioffs+ip) = rh2(ioffs+ip)
        else
         rhv(ioffs+ip) = rh1(ioffs+ip)
        end if

	      endif                   ! mode
	 
	      enddo	!ip
	 
        return
	      end
c---------------------------------------------------------------------
	      subroutine data_to_hist(ND,Nb,p,D1,D2,D3,D4)
c---------------------------------------------------------------------
c Transforms Ddim data sets of different sizes containing in total ND 
c data points using Nb bins.
c
c Simplified version, exploiting the fact that we know xa=0, xb=1.
c
c Do normalization and output INSIDE the routine
c      nx(k)	# of events in k-th component of xdat
c
	      implicit none
	      integer (kind=4):: ND,id
	      integer	(kind=4):: Nb,ind,ik,i
	      real	(kind=4):: D1(ND),D2(ND),D3(ND),D4(ND),x,xa,xb,dx,dp,p
	      real	(kind=4):: PDF(Nb,4)

	      xa= -0.0005
	      xb= 1.0005
	      dx=(xb-xa)/Nb
	      dp = 1.0/(ND*dx)	! probility density
c	      dp = 1.0/ND   		! probabilities
	      PDF = 0.0

	      do id=1,ND
	        ind=1+int((D1(id)-xa)/dx) 
	        PDF(ind,1)=PDF(ind,1)+dp
	        ind=1+int((D2(id)-xa)/dx) 
	        PDF(ind,2)=PDF(ind,2)+dp
	        ind=1+int((D3(id)-xa)/dx) 
	        PDF(ind,3)=PDF(ind,3)+dp
	        ind=1+int((D4(id)-xa)/dx) 
	        PDF(ind,4)=PDF(ind,4)+dp
	      enddo

        do i=1,Nb
	        x=xa + (i-0.5)*dx
	        write(10,100) p,x,PDF(i,1),PDF(i,2),PDF(i,3),PDF(i,4)
	      enddo
	      return
100	    format(14e16.7)
	      end
c**********************************************************
      subroutine iran_d(p,n)
c----------------------------------------------------------
c  initializes discrete random 'number' generator
c    to prepare generating random elements from an
c    array of length n, with probabilities \propto p
c
c	initialization replaces possibly 
c	    UNNORMALIZED p with the cdf
c
      integer (kind=4):: n,i
      real	(kind=4):: p(n),s

      s=sum(p)
      p=p/s

      s=0.0
      do i=1,n
         s=s+p(i)
         p(i)=s
      end do
      return
      end
c----------------------------------------------------------
      function fran_d(k,p,n,ids)
c----------------------------------------------------------
c  generates an ineger random number from  (n values in) k, 
c  distributed according to p, 
c
      integer (kind=4):: n,ids,i,j,im,ir
      integer (kind=4):: k(n)
      real	(kind=4):: p(n),s,pm
      real	(kind=4):: ran0, fran_d

      s=ran0(ids)
      i=0
      j=n
      do while(i+1 .lt. j)
	      im=(i+j)/2
	      if(s .lt. p(im)) then
	        j=im
	      else
	        i=im
	      end if
	    end do
	    fran_d=k(j)
      return
      end
c******************************************************************************
c	includes
c******************************************************************************
       include 'Subs/ran0.f'                      
                                                            
