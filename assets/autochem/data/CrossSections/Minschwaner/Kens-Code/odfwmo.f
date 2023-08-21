c 
c	Same as makeodf.f, but uses WMO intervals 
c
c	Computes ODF for S-R bands
c
	parameter(nv=17000,nvmax=20000,nt=13,nbin=17)
	dimension v(nv),x2(nv,nt),tx2(nt),vg(nvmax),x2g(nvmax),xmid(6)
     *		 ,xx(nvmax),yy(nvmax),iwksp(nvmax)
c	Yoshino '88 Herz coeffs
	real a(5)/4.67604e3,-9.66666e1,0.7391704,-2.478075e-3,3.07613e-6/
c
c       WMO wavelength bins.  Note last bin extends to 48500 cm-1
        real vbin(18)/1754.,1770.,1786.,1802.,1818.,1835.,1852.,
     *  1869.,1887.,1905.,1923.,1942.,1961.,1980.,2000.,2020.,2041.,
     *  2062./
c
c       Estimated temperatures, based on closest Prather bin and/or interp
        real tbin(nbin)/205.,205.,210.,205.,235.,220.,225.,235.,
     *  250.,260.,255.,240.,230.,225.,225.,225.,225./
c
	real gint(7)/0.00,0.05,0.25,0.50,0.75,0.95,1.00/
c
	open(unit=1,name='spectrum.cold.plumb',
     *	form='formatted',status='old') 
	read(1,1)
1	format(/)
	read(1,5) (tx2(j),j=1,nt) 
5	format(3x,13f7.1) 
	do i=1,nv 
	     read(1,10) v(i),(x2(i,j),j=1,nt)
	     if(v(i).lt.52000) then
		wl=1.e7/v(i)
		herz=a(5)
		do j=4,1,-1
			herz=herz*wl+a(j)
		end do
		do j=1,nt
		   x2(i,j)=x2(i,j)+herz*1.e-24
		end do
	     end if	
	end do 
10	format(1x,f7.1,2x,13(e10.3))
	close(unit=1)
c
c
	do i=1,nbin
c		write(*,11) vbin(i), vbin(i+1), tbin(i)
11	format('Wavelength:',f5.0,'-',f5.0,10x,'Mean Temperature:',f4.0)
		call hunt(tx2,nt,tbin(i),klo) 
	   	if(klo.eq.0.or.klo.eq.nt) then
		   write(*,*) 'error in o2 xsec interp, vbin=',vbin(i)
		   pause
	   	end if
		khi=klo+1 
		k=0 
		do j=1,nv 
	   		wl=1.e8/v(j)
			start=vbin(i)
			end=vbin(i+1)
			if(wl.ge.start.and.wl.lt.end) then 
			   k=k+1
			   vg(k)=v(j) 
			   x2g(k)=x2(j,klo)*(tx2(khi)-tbin(i))+
     *			   x2(j,khi)*(tbin(i)-tx2(klo))
			   x2g(k)=x2g(k)/(tx2(khi)-tx2(klo))
			end if 
		end do 
		call gk(vg,x2g,k,yy,xx,nvmax,ng)
c		do j=1,ng
c			write(*,*) xx(j),yy(j)
c		end do
c
c		Get mid points for O2 ODF
		do j=1,6
			gmid=(gint(j)+gint(j+1))/2.
			do k=1,ng
			   if(gmid.ge.xx(k).and.gmid.lt.xx(k+1)) then
c			write(*,*) 'gmid,xx(k),xx(k+1)',gmid,xx(k),xx(k+1)
c			write(*,*) 'yy(k),yy(k+1)',yy(k),yy(k+1)
				den=xx(k+1)-xx(k)
				slope=(yy(k+1)-yy(k))/den
				yint=(yy(k)*xx(k+1)-yy(k+1)*xx(k))/den
				xmid(j)=slope*gmid+yint
			   end if
			end do
		end do
		write(*,12) vbin(i),vbin(i+1),tbin(i),(xmid(j),j=1,6)
	end do
12	format(f5.0,'-',f5.0,1x,f4.0,6(1pe11.3))
	end
c
c
	subroutine gk(wavnum,absorp,nv,y,x,nvmax,ng) 
	parameter(nmax=10000)
	dimension wavnum(nv),absorp(nv),y(nvmax),x(nvmax) 
	dimension iwksp(nmax)
	if(nv.gt.nmax) then
		write(*,*) 'error: nv > nmax in subr gk'
		pause
	end if
	call indexx(nv,absorp,iwksp)
	do i=1,nv
		y(i)=absorp(iwksp(i))
		x(i)=float(i)/float(nv)
	end do
	ng=nv
	return
	end
c
      SUBROUTINE HUNT(XX,N,X,JLO)
      DIMENSION XX(N)
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END
c
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

