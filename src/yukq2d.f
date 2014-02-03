cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c        This is the end of the debugging code and the beginning of the
c        Yukawa quadrature code proper.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Usage:
c	call yukq2d(ier,rky,ndigits,nlams,xs,zw,w,lw,lused)
c Inputs: 
c	rky---the frequency: rky \in (0.0,\infty); for frequency beyond
c		this interval, Yukawa is essentially meaningless.
c	ndigits---the accuracy,allowed value:1,2,3,4,5,6,7,8,9,10,11,12.
c	w,lw---workspace(real *8), set "lw" large first, use "lused" to
c		gauge the actual space used and to reset "lw".
c Outputs:
c	ier---tells you if there is some error, if not zero, bomb!
c	nlams,xs,zw---the returned the quadrature, the nodes xs(nlams)
c		are real, but I set the weights zw(nlams) to be complex,
c		HOWEVER, in all tests that I have conducted, their 
c		imaginary parts all come back zero. Therefore, if you
c		insist upon the weights zw to be real, you can simply
c		use another array and set it to be the real parts of zw.
c	lused---tell you how much work space is actually used.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	subroutine yukq2d(ier,rky,ndigits,nlams,xs,zw,w,lw,lused)
	subroutine yukq2d(ier,rky,nlams,xs,zw,w,lw,lused)
	implicit real *8 (a-h,o-z)
c
c Purpose:
c	This subroutine gives the nodes and weights for the quadrature
c	of H_0(k*||z||) through the integral representation:
c
c	   H_0(k*||z||)=1/pi*int_0^{\infty}\frac{(exp(eye*y*rk)+
c		exp(-eye*y*rk))*exp(-x*\sqrt{rk^2-k^2})}
c		{\sqrt{k^2-rk^2}} d(rk)
c	    		[note: z=x+eye*y and eye=\sqrt{-1}]
c
c	(the discrete version is:
c	   H_0(k*||z||)=sum_{i=1}^{nlams} zw(i)*\frac{
c		(exp(eye*y*zx(i))+exp(-eye*y*zx(i)))*exp(-x*
c		  \sqrt{zx(i)^2-k^2})}{\sqrt{k^2-zx(i)^2}} )
c
c	where k is the frequency value of the Yukawa equation.
c
c	The qudratures here are good for IMAGINARY k in the range of:
c		0 < rky <=8*pi.
c	Beyond this range, the qudratures here have no use.
c
c	The nodes and weights here will give "ndigits" accuracy for any 
c	z=(x,y) inside the following 2d FMM region (see lapquad.f):
c
c			     y
c			     |
c		      y=4__  |     ___
c		      y=3__  |   _|   |
c		             |  |     |
c		             |  |     |
c		             |__|_____|_____ x
c		            0   |     |  
c		                |     |
c		      y=-3__    |_    |
c		      y=-4__      |___|
c		
c		                  ^   ^
c		                ^ |   x=4
c		                | x=2
c		               x=1
c
c Inputs:
c	rky--- the frequency parameter.
c	ndigits--- the accuracy required.
c
c Outputs:
c	ier--- error return code
c		ier=0 successful return
c		ier=8 fatal error, rk not in the proper range
c		ier=16 fatal error, ndigits not proper
c	nlams--- integer indicating the number of nodes needed
c	xs(nlams)--- the real nodes
c	zw(nlams)--- the complex weights
c	lused--- it reports the amount of work space that is used
c
c Workspace:
c	(w,lw)
c
c
	dimension xs(*),w(*),ws(2 000)
	doublecomplex  rk7,zw(*),eye
c
	data eye/(0.0d0,1.0d0)/
c
        ndigits = 12
	ier=0
	rk7=rky*eye
c
c ... Calling for appropriate nodes (based on accuracy):
c	WARNING: Do not change the value of eps below!
c
	if (rky.le.0.0d0) then
	   ier=8
	   return
	endif
c
	if ((ndigits.lt.0).or.(ndigits.gt.12)) then
	   ier=16
	   return
	endif
c
	if ((ndigits.gt.0).and.(ndigits.le.3)) then
	   call yukq2d_tab3(jer,rky,nlams,xs)
	   eps=1.0d-16
	   call getwts(jer,eps,rk7,nlams,xs,zw,w,lw,lused)
	   return
	endif
c
	if ((ndigits.gt.3).and.(ndigits.le.6)) then
	   call yukq2d_tab6(jer,rky,nlams,xs)
	   eps=1.0d-16
	   call getwts(jer,eps,rk7,nlams,xs,zw,w,lw,lused)
	   return
	endif
c
	if ((ndigits.gt.6).and.(ndigits.le.9)) then
	   call yukq2d_tab9(jer,rky,nlams,xs)
	   eps=1.0d-19
	   if (rky.lt.0.0001d0) eps=1.0d-23
	   call getwts(jer,eps,rk7,nlams,xs,zw,w,lw,lused)
	   return
	endif
c
	if ((ndigits.gt.9).and.(ndigits.le.12)) then
	   call yukq2d_tab12(jer,rky,nlams,xs,ws)
	   do i=1,nlams
	      zw(i)=ws(i)
	   enddo
	   return
	endif
c
	end
c
c
c
c
c
	subroutine getwts(ier,eps,rk,nlams,xsc,wsc,w,lw,lused)
	implicit real *8 (a-h,o-z)
c
c Purpose:
c	This subroutine generates the weights for a generalized
c	Gaussian quadrature, given the nodes. The method is simply
c	the Least Square. This is the memory management routine,
c	actual code is in getwhts0 (see below)
c
c Inputs:
c	eps: the error control for LEAST SQUARE
c	rk: the frequency under consideration
c	nlams: the number of quad nodes
c	xsc(nlams): the nodes
c
c Outputs:
c	ier: error output code
c		=0 no error
c		=8 not enough memory in w(lw)
c	wsc(nlams): the weights
c
c Notes: xsc(nlams) the nodes are REAL, while wsc(nlams) are COMPLEX!
c
c
	real *8 xsc(*),w(*)
	doublecomplex  wsc(*),rk
c
c ... First creat the struture:
c
	n22=13
	ier=0
	iat=1
c
	ixat=iat
	lxat=n22*20+7
	iat=ixat+lxat
c
	iyat=iat
	lyat=n22*20+7
	iat=iyat+lyat
c
	iwhts=iat
c
	delta=0.01d0
	call creastr(n22,delta,npts,w(ixat),w(iyat),w(iwhts))
c
ccc	if (npts.ne.n22*20) STOP 'Error in getwts!'
c
	iamatr=iat
	lamatr=2*npts*nlams+7
	iat=iamatr+lamatr
c
	irhs=iat
	lrhs=2*npts+7
	iat=irhs+lrhs
c
	irnorms=iat
	lrnorms=nlams+1+7
	iat=irnorms+lrnorms
c
	iwhts=iat
	lwhts=nlams+7
	iat=iwhts+lwhts
c
	iu=iat
	lu=2*npts*nlams+7
	iat=iu+lu
c
	iv=iat
	lv=2*npts*nlams+7
	iat=iv+lv
c
	iw=iat
	lwww=2*npts*nlams+7
	iat=iw+lwww
c
	is=iat
	ls=2*npts*nlams+7
	iat=is+ls
c
	iwork=iat
	lwork=2*nlams+7
	iat=iwork+lwork
c
	lused=iat-1
	if (lused.gt.lw) goto 6008
c
	call getwts0(eps,rk,npts,w(ixat),w(iyat),nlams,xsc,wsc,
     1	   w(iamatr),w(irhs),w(irnorms),w(iwhts),w(iu),w(iv),w(iw),
     2	   w(is),w(iwork))
        return
c
 6008 ier=8
ccc	call prinf('iat=*',iat,1)
	return
c
        end
c
c
c
c
c
	subroutine getwts0(eps,rk,npts,xat,yat,nlams,xsc,wsc,amatr,rhs,
     1	   rnorms,whts,u,v,w,s,work)
	implicit real *8 (a-h,o-z)
c
	real *8 xsc(nlams)
	real *8 xat(npts),yat(npts)
	real *8 rnorms(1),whts(1)
	doublecomplex  wsc(nlams),rk,funinte
	doublecomplex  amatr(npts,nlams),rhs(npts)
	doublecomplex  u(1),v(1)
	doublecomplex  w(1),s(1),work(1)
c
	doublecomplex  zz,true,h1,zxs,eye
c
ccc	external funinte,funintin
c
	data eye/(0.0d0,1.0d0)/
c
c ... Create the matrix and rhs:
c
	do 1400 i=1,npts
	   xx=xat(i)
	   yy=yat(i)
c
c$$$	   yyy2=funintin(yy)
        do 1200 j=1,nlams
	   zxs=xsc(j)
c$$$	   amatr(i,j)=funinte(zxs,xx,rk)
	   amatr(i,j)=funinte(zxs,xx,yy,rk)
 1200 continue
	   zz=sqrt(xx*xx+yy*yy)*rk
	   call hank101(zz,true,h1)
c
	   rhs(i)=true
 1400 continue
ccc	call prin2('rhs=*',rhs,2*npts)
c
c ... Use least square to solve for the weights:
c	NOTE: the "eps" needs tobe changed for 9- and 12-digit cases!
c
ccc	eps=1.0d-16
	call cleastsq(amatr,u,w,s,npts,nlams,ncols,rnorms,eps,v)
	call cleasts2(u,w,s,npts,nlams,ncols,rhs,wsc,work,whts)
c
ccc	call prinf('after least square, ncols=*',ncols,1)
ccc	call prin2('and the rnorms=*',rnorms,ncols)
ccc	call prin2('the computed weights are, wsc=*',wsc,2*nlams)
c
        return
        end
c
c
c
c
c
	function funinte(zlam,x,y,rk)
	implicit real *8 (a-h,o-z)
	save
c
	doublecomplex  funinte,zlam,rk
	doublecomplex  eye,cd,zzz
	data eye/(0.0d0,1.0d0)/
c
c        this subroutine evaluates the function of rlam, x,y,rk
c        that, when integrated with respect to rlam, gives
c
c        H_0(rk * sqrt(x^2+y^2) ),
c
c        which is, obviously, the Green's function for the
c        Helmholtz equation in two dimensions. Please note
c        that for technical reasons, y is passed through the
c        auxiliary entry funintin (see below). The parameter
c        rk is assumed to be complex, with a positive imaginary
c        part - this is not the standard integral formula.
c
c
	eps2=1.0d-19
	cd=sqrt(zlam*zlam-rk*rk)
	if (abs(cd).le.eps2) then
	   funinte=0
	   return
	endif
c
	zzz=exp(eye*y*zlam)
	funinte=(zzz+1/zzz)*exp(-x*cd)/(cd*eye)
c
	return
	end
c
c
c
c
c

c$$$	function funinte(zlam,x,rk)
c$$$	implicit real *8 (a-h,o-z)
c$$$	save
c$$$c
c$$$	doublecomplex  funinte,zlam,rk
c$$$	doublecomplex  eye,cd,zzz
c$$$	data eye/(0.0d0,1.0d0)/
c$$$c
c$$$c        this subroutine evaluates the function of rlam, x,y,rk
c$$$c        that, when integrated with respect to rlam, gives
c$$$c
c$$$c        H_0(rk * sqrt(x^2+y^2) ),
c$$$c
c$$$c        which is, obviously, the Green's function for the
c$$$c        Helmholtz equation in two dimensions. Please note
c$$$c        that for technical reasons, y is passed through the
c$$$c        auxiliary entry funintin (see below). The parameter
c$$$c        rk is assumed to be complex, with a positive imaginary
c$$$c        part - this is not the standard integral formula.
c$$$c
c$$$c
c$$$	eps2=1.0d-19
c$$$	cd=sqrt(zlam*zlam-rk*rk)
c$$$	if (abs(cd).le.eps2) then
c$$$	   funinte=0
c$$$	   return
c$$$	endif
c$$$c
c$$$	zzz=exp(eye*y*zlam)
c$$$	funinte=(zzz+1/zzz)*exp(-x*cd)/(cd*eye)
c$$$c
c$$$	return
c$$$c
c$$$c
c$$$c
c$$$c
c$$$	entry funintin(y7)
c$$$	y=y7
c$$$	funintin=y
c$$$c
c$$$	return
c$$$	end
c
c
c
c
c


	subroutine creastr(n22,delta,npts,xat,yat,whts)
	implicit real *8 (a-h,o-z)
c
c	This subroutine creates a structure for FMM in two dismensions:
c	(see Rokhlin & Yarvin, YaleU/DCS/RR-1109, p24)
c
c				    (4,4)
c				____
c		      (1,3)  __|    |
c			    | (2,3) |
c			    |	    |
c			    |	    |
c			    |	    |
c			    | (2,-3)|
c			    |__     |
c		      (1,-3)   |____|
c				    (-4,4)
c
c	Note: unlike the Laplace case, the Helmhotz case does not have
c	the luxury of Maximum Principle, therefore, a "delta" coating
c	of the above structure is added to match not only the function
c	but also the the derivative.
c
	dimension xat(1),yat(1)
	dimension whts(1)
	dimension gpts(1 000),gwhts(1 000)
c
	zero=0.0d0
	half=1.0d0/2
	xmin=1.0d0
	xmid=2.0d0
	ymid=3.0d0
	rmax=4.0d0
c
c        put n22 Gaussian points on [-1,1]:
c
	itype=1
	call legeexps(itype,n22,gpts,uu,vv,gwhts)
c
	ncur=0
c
c	1) Put points on line segment (1,0)->(1,3):
c	 and its coating layer (1-delta,0)->(1-delta,3)
c
	aa=zero
	bb=ymid
	alpha=(bb-aa)/2
	beta=(bb+aa)/2
	do 1200 i=1,n22
	   iat=ncur+i
	   xat(iat)=xmin
	   yat(iat)=gpts(i)*alpha+beta
	   whts(iat)=gwhts(i)
 1200 continue
	ncur=ncur+n22
c
	do 1300 i=1,n22
	   iat=ncur+i
	   jat=iat-n22
	   xat(iat)=xat(jat)-delta
	   yat(iat)=yat(jat)
	   whts(iat)=whts(jat)
 1300 continue
	ncur=ncur+n22
c
c	2) Put points on line segment (1,3)->(2,3):
c	 and its coating layer (1,3+delta)->(2,3+delta)
c
	aa=xmin
	bb=xmid
	alpha=(bb-aa)/2
	beta=(bb+aa)/2
	do 1400 i=1,n22
	   iat=ncur+i
	   xat(iat)=gpts(i)*alpha+beta
	   yat(iat)=ymid
	   whts(iat)=gwhts(i)
 1400 continue
	ncur=ncur+n22
c
	do 1500 i=1,n22
	   iat=ncur+i
	   jat=iat-n22
	   xat(iat)=xat(jat)
	   yat(iat)=yat(jat)+delta
	   whts(iat)=whts(jat)
 1500 continue
	ncur=ncur+n22
c
c	3) Put points on line segment (2,3)->(2,4):
c	 and its coating layer (2-delta,3)->(2-delta,4)
c
	aa=ymid
	bb=rmax
	alpha=(bb-aa)/2
	beta=(bb+aa)/2
	do 2200 i=1,n22
	   iat=ncur+i
	   xat(iat)=xmid
	   yat(iat)=gpts(i)*alpha+beta
	   whts(iat)=gwhts(i)
 2200 continue
	ncur=ncur+n22
c
	do 2300 i=1,n22
	   iat=ncur+i
	   jat=iat-n22
	   xat(iat)=xat(jat)-delta
	   yat(iat)=yat(jat)
	   whts(iat)=whts(jat)
 2300 continue
	ncur=ncur+n22
c
c	4) Put points on line segment (2,4)->(4,4)
c	 and its coating layer (2,4+delta)->(4,4+delta)
c
	aa=xmid
	bb=rmax
	alpha=(bb-aa)/2
	beta=(bb+aa)/2
	do 2400 i=1,n22
	   iat=ncur+i
	   xat(iat)=gpts(i)*alpha+beta
	   yat(iat)=rmax
	   whts(iat)=gwhts(i)
 2400 continue
	ncur=ncur+n22
c
	do 2500 i=1,n22
	   iat=ncur+i
	   jat=iat-n22
	   xat(iat)=xat(jat)
	   yat(iat)=yat(jat)+delta
	   whts(iat)=whts(jat)
 2500 continue
	ncur=ncur+n22
c
c	5) Put points on line segment (4,4)->(4,0):
c	 and its coating layer (4+delta,4)->(4+delta,0)
c
	aa=zero
	bb=rmax
	alpha=(bb-aa)/2
	beta=(bb+aa)/2
	do 3200 i=1,n22
	   iat=ncur+i
	   xat(iat)=rmax
	   yat(iat)=gpts(n22-i+1)*alpha+beta
	   whts(iat)=gwhts(n22-i+1)
 3200 continue
	ncur=ncur+n22
c
	do 3300 i=1,n22
	   iat=ncur+i
	   jat=iat-n22
	   xat(iat)=xat(jat)+delta
	   yat(iat)=yat(jat)
	   whts(iat)=whts(jat)
 3300 continue
	ncur=ncur+n22
c
c	The bottom half can be obtained by flipping the up half
c
	do 4200 i=1,ncur
	   iat=ncur+i
	   xat(iat)=xat(i)
	   yat(iat)=-yat(i)
	   whts(iat)=whts(i)
 4200 continue
c
	npts=2*ncur
c
        return
        end
c
c
c
c
c
	subroutine yukq2d_tab3(ier,rky,nlam,xs)
	implicit real *8 (a-h,o-z)
c
c Purpose:
c	This is the quadrature table (nodes only) for yukawa 2D
c	with 3 digits of accuracy. It covers only the range
c	(0,8*pi]!
c
c Inputs:
c	rky: the frequency=rky*i
c
c Outputs:
c	ier: error output code
c		=0 no error
c		=8 rky out of range of (0,\infty)
c	nlam: number of nodes
c	xs(nlam): the quad nodes
c
	dimension xs(*)
c
	ier=0
	done=1
	pi=4*atan(done)
c
c ... If rky lies outside of (0,\infty), BOMB!
c
cccc	if ((rky.le.0.0d0).or.(rky.gt.8*pi)) then
	if (rky.lt.0.0d0) then
	   ier=8
	   return
	endif
c
c ... Choose the correct the nodes:
c
	if ((rky.ge.0.00d0).and.(rky.le.0.01d0))   goto 1001
	if ((rky.gt.0.01d0).and.(rky.le.0.10d0))   goto 1002
	if ((rky.gt.0.10d0).and.(rky.le.0.30d0))   goto 1003
	if ((rky.gt.0.30d0).and.(rky.le.0.70d0))   goto 1004
	if ((rky.gt.0.70d0).and.(rky.le.1.50d0))   goto 1005
	if ((rky.gt.1.50d0).and.(rky.le.2.50d0))   goto 1006
	if ((rky.gt.2.50d0).and.(rky.le.3.50d0))   goto 1007
	if (rky.gt.3.50d0)  goto 1008
cccc	if ((rky.gt.3.50d0).and.(rky.le.8d0*pi))   goto 1008
c
c ... If for some reason, the code gets here, bomb!
c
	ier=8
	return
c
c	Data for frequency k in region [ 0.00000, 0.01000]
c
 1001 continue
	nlam= 9
c
	xs( 1)= 0.2153145857164418D-06
	xs( 2)= 0.2415379757238156D-02
	xs( 3)= 0.1777099222612932D-01
	xs( 4)= 0.1737270157641291D+00
	xs( 5)= 0.6647465015718197D+00
	xs( 6)= 0.1467532694995345D+01
	xs( 7)= 0.2545555821530662D+01
	xs( 8)= 0.3876845321200079D+01
	xs( 9)= 0.5418927174712078D+01
c
	return
c
c	Data for frequency k in region [ 0.01000, 0.10000]
c
 1002 continue
	nlam= 8
c
	xs( 1)= 0.5083701804480256D-03
	xs( 2)= 0.6194121633175342D-01
	xs( 3)= 0.2827367271631154D+00
	xs( 4)= 0.7862437649170673D+00
	xs( 5)= 0.1585409055346098D+01
	xs( 6)= 0.2660620212775363D+01
	xs( 7)= 0.3998991025902821D+01
	xs( 8)= 0.5559536920492871D+01
c
	return
c
c	Data for frequency k in region [ 0.10000, 0.30000]
c
 1003 continue
	nlam= 7
c
	xs( 1)= 0.2571321690403749D-02
	xs( 2)= 0.1876479243831710D+00
	xs( 3)= 0.5755421898078836D+00
	xs( 4)= 0.1255058370685905D+01
	xs( 5)= 0.2224265217886611D+01
	xs( 6)= 0.3474018599608033D+01
	xs( 7)= 0.4986453908391322D+01
c
	return
c
c	Data for frequency k in region [ 0.30000, 0.70000]
c
 1004 continue
	nlam= 7
c
	xs( 1)= 0.5620131353737179D-02
	xs( 2)= 0.3587094028115381D+00
	xs( 3)= 0.8997156464662019D+00
	xs( 4)= 0.1691563879922619D+01
	xs( 5)= 0.2752988514647157D+01
	xs( 6)= 0.4083339196325625D+01
	xs( 7)= 0.5644176927058862D+01
c
	return
c
c	Data for frequency k in region [ 0.70000, 1.50000]
c
 1005 continue
	nlam= 6
c
	xs( 1)= 0.1016254186355070D-01
	xs( 2)= 0.6087354614430573D+00
	xs( 3)= 0.1387251540061051D+01
	xs( 4)= 0.2397677224410881D+01
	xs( 5)= 0.3654856147021022D+01
	xs( 6)= 0.5149802176621115D+01
c
	return
c
c	Data for frequency k in region [ 1.50000, 2.50000]
c
 1006 continue
	nlam= 5
c
	xs( 1)= 0.1599160736028438D-01
	xs( 2)= 0.9196686659673432D+00
	xs( 3)= 0.1990047287856541D+01
	xs( 4)= 0.3290168451588450D+01
	xs( 5)= 0.4809543603924954D+01
c
	return
c
c	Data for frequency k in region [ 2.50000, 3.50000]
c
 1007 continue
	nlam= 5
c
	xs( 1)= 0.1876033397081756D-01
	xs( 2)= 0.1064911899812323D+01
	xs( 3)= 0.2256316804999507D+01
	xs( 4)= 0.3625567807913974D+01
	xs( 5)= 0.5175314128520579D+01
c
	return
c
c	Data for frequency k in region [ 3.50000, 8*pi]
c
 1008 continue
	nlam= 4
c
	xs( 1)= 0.2551106321274332D-01
	xs( 2)= 0.1410093756179506D+01
	xs( 3)= 0.2892937104088659D+01
	xs( 4)= 0.4470619438480664D+01
c
	return
c
	end
c
c
c
c
c
	subroutine yukq2d_tab6(ier,rky,nlam,xs)
	implicit real *8 (a-h,o-z)
c
c Purpose:
c	This is the quadrature table (nodes only) for yukawa 2D
c	with 6 digits of accuracy. It covers only the range
c	(0,8*pi]!
c
c Inputs:
c	rky: the frequency=rky*i
c
c Outputs:
c	ier: error output code
c		=0 no error
c		=8 rky out of range of (0,\infty)
c	nlam: number of nodes
c	xs(nlam): the quad nodes
c
	dimension xs(*)
c
	ier=0
	done=1
	pi=4*atan(done)
c
c ... If rky lies outside of (0,\infty), BOMB!
c
cccc	if ((rky.le.0.0d0).or.(rky.gt.8*pi)) then
	if (rky.lt.0.0d0) then
	   ier=8
	   return
	endif
c
c ... Choose the correct the nodes:
c
	if ((rky.ge.0.000d0).and.(rky.le.0.001d0))   goto 1001
	if ((rky.gt.0.001d0).and.(rky.le.0.005d0))   goto 1002
	if ((rky.gt.0.005d0).and.(rky.le.0.020d0))   goto 1003
	if ((rky.gt.0.020d0).and.(rky.le.0.080d0))   goto 1004
	if ((rky.gt.0.080d0).and.(rky.le.0.300d0))   goto 1005
	if ((rky.gt.0.300d0).and.(rky.le.0.900d0))   goto 1006
	if ((rky.gt.0.900d0).and.(rky.le.1.800d0))   goto 1007
	if ((rky.gt.1.800d0).and.(rky.le.3.000d0))   goto 1008
	if ((rky.gt.3.000d0).and.(rky.le.3.d0*pi))   goto 1009
	if (rky.gt.3.d0*pi)   goto 1010
cccc	if ((rky.gt.3.d0*pi).and.(rky.le.8.d0*pi))   goto 1010
c
c ... If for some reason, the code gets here, bomb!
c
	ier=8
	return
c
c	Data for frequency k in region [ 0.00000, 0.00100]
c
 1001 continue
	nlam=16
c
	xs( 1)= 0.5855119859177194D-05
	xs( 2)= 0.1385730047285705D-02
	xs( 3)= 0.5799518172791807D-01
	xs( 4)= 0.2765056999495350D+00
	xs( 5)= 0.6576878692521575D+00
	xs( 6)= 0.1187529627477357D+01
	xs( 7)= 0.1852933920932802D+01
	xs( 8)= 0.2643678413701767D+01
	xs( 9)= 0.3552785652476516D+01
	xs(10)= 0.4576648332409285D+01
	xs(11)= 0.5715151510840393D+01
	xs(12)= 0.6970914622855577D+01
	xs(13)= 0.8346569624497910D+01
	xs(14)= 0.9850738591325090D+01
	xs(15)= 0.1148005593860510D+02
	xs(16)= 0.1326367977003878D+02
c
	return
c
c	Data for frequency k in region [ 0.00100, 0.00500]
c
 1002 continue
	nlam=16
c
	xs( 1)= 0.4230065586696696D-04
	xs( 2)= 0.5817329976284924D-02
	xs( 3)= 0.7090138428958781D-01
	xs( 4)= 0.2904521924384653D+00
	xs( 5)= 0.6717728805603791D+00
	xs( 6)= 0.1201605675115268D+01
	xs( 7)= 0.1866957214884323D+01
	xs( 8)= 0.2657628604520866D+01
	xs( 9)= 0.3566646046316531D+01
	xs(10)= 0.4590399825955533D+01
	xs(11)= 0.5728766611669425D+01
	xs(12)= 0.6984350956107754D+01
	xs(13)= 0.8359707617712228D+01
	xs(14)= 0.9863755292246605D+01
	xs(15)= 0.1149138490174764D+02
	xs(16)= 0.1327697161790279D+02
c
	return
c
c	Data for frequency k in region [ 0.00500, 0.02000]
c
 1003 continue
	nlam=16
c
	xs( 1)= 0.1743667945355298D-03
	xs( 2)= 0.1750281318077285D-01
	xs( 3)= 0.9882100989404563D-01
	xs( 4)= 0.3226194931807242D+00
	xs( 5)= 0.7047010999206904D+00
	xs( 6)= 0.1234523451997269D+01
	xs( 7)= 0.1899631899354574D+01
	xs( 8)= 0.2689959098881882D+01
	xs( 9)= 0.3598556604323441D+01
	xs(10)= 0.4621809459054997D+01
	xs(11)= 0.5759569213869600D+01
	xs(12)= 0.7014397107216631D+01
	xs(13)= 0.8388690071985138D+01
	xs(14)= 0.9891943241820773D+01
	xs(15)= 0.1151565281375130D+02
	xs(16)= 0.1330316257590545D+02
c
	return
c
c	Data for frequency k in region [ 0.02000, 0.08000]
c
 1004 continue
	nlam=16
c
	xs( 1)= 0.4653581501017356D-03
	xs( 2)= 0.3428455690866983D-01
	xs( 3)= 0.1140735459023823D+00
	xs( 4)= 0.2979898464223254D+00
	xs( 5)= 0.6286705525194165D+00
	xs( 6)= 0.1112065236093013D+01
	xs( 7)= 0.1738929755420060D+01
	xs( 8)= 0.2498562994691081D+01
	xs( 9)= 0.3382565963477695D+01
	xs(10)= 0.4385849080769798D+01
	xs(11)= 0.5507162200013983D+01
	xs(12)= 0.6748303323475835D+01
	xs(13)= 0.8112075339590554D+01
	xs(14)= 0.9600662439643365D+01
	xs(15)= 0.1123905119227398D+02
	xs(16)= 0.1292885540533359D+02
c
	return
c
c	Data for frequency k in region [ 0.08000, 0.30000]
c
 1005 continue
	nlam=16
c
	xs( 1)= 0.1440829094580209D-02
	xs( 2)= 0.9463734661340339D-01
	xs( 3)= 0.2487101400241301D+00
	xs( 4)= 0.5021410735587466D+00
	xs( 5)= 0.8812764350698430D+00
	xs( 6)= 0.1396135845926448D+01
	xs( 7)= 0.2044269782917588D+01
	xs( 8)= 0.2819030058437153D+01
	xs( 9)= 0.3714248071979089D+01
	xs(10)= 0.4726195598077243D+01
	xs(11)= 0.5854387001992340D+01
	xs(12)= 0.7101046123381430D+01
	xs(13)= 0.8467987099337385D+01
	xs(14)= 0.9965736549994418D+01
	xs(15)= 0.1157744061647766D+02
	xs(16)= 0.1336625658949028D+02
c
	return
c
c	Data for frequency k in region [ 0.30000, 0.90000]
c
 1006 continue
	nlam=16
c
	xs( 1)= 0.2940844763056560D-02
	xs( 2)= 0.1751744207125476D+00
	xs( 3)= 0.3873270735916883D+00
	xs( 4)= 0.6682619812182438D+00
	xs( 5)= 0.1045145989604332D+01
	xs( 6)= 0.1538687936944509D+01
	xs( 7)= 0.2159396023577582D+01
	xs( 8)= 0.2908658562714695D+01
	xs( 9)= 0.3783311305473152D+01
	xs(10)= 0.4779913730699093D+01
	xs(11)= 0.5897280067994870D+01
	xs(12)= 0.7136676069155296D+01
	xs(13)= 0.8498701928179784D+01
	xs(14)= 0.9994224520583327D+01
	xs(15)= 0.1160043026351689D+02
	xs(16)= 0.1339079875382269D+02
c
	return
c
c	Data for frequency k in region [ 0.90000, 1.80000]
c
 1007 continue
	nlam=15
c
	xs( 1)= 0.7191006118087984D-02
	xs( 2)= 0.4157316592269389D+00
	xs( 3)= 0.8794159843471263D+00
	xs( 4)= 0.1425302073029187D+01
	xs( 5)= 0.2072247943152981D+01
	xs( 6)= 0.2829264591583637D+01
	xs( 7)= 0.3699706756901708D+01
	xs( 8)= 0.4684469151067514D+01
	xs( 9)= 0.5784682825915322D+01
	xs(10)= 0.7002729241320779D+01
	xs(11)= 0.8341133376483807D+01
	xs(12)= 0.9801065691522098D+01
	xs(13)= 0.1138768233314175D+02
	xs(14)= 0.1309097523688911D+02
	xs(15)= 0.1497498341437167D+02
c
	return
c
c	Data for frequency k in region [ 1.80000, 3.00000]
c
 1008 continue
	nlam=13
c
	xs( 1)= 0.1106299115027909D-01
	xs( 2)= 0.6284342078440055D+00
	xs( 3)= 0.1298399101322833D+01
	xs( 4)= 0.2042381098776097D+01
	xs( 5)= 0.2880094457895659D+01
	xs( 6)= 0.3822354641883365D+01
	xs( 7)= 0.4876090232843193D+01
	xs( 8)= 0.6047551794269307D+01
	xs( 9)= 0.7341798692432253D+01
	xs(10)= 0.8760208608245911D+01
	xs(11)= 0.1031256594317095D+02
	xs(12)= 0.1197650825214881D+02
	xs(13)= 0.1379011250854491D+02
c
	return
c
c	Data for frequency k in region [ 3.00000, 9.42478]
c
 1009 continue
	nlam=11
c
	xs( 1)= 0.1511749431794840D-01
	xs( 2)= 0.8486795532140050D+00
	xs( 3)= 0.1732819001003794D+01
	xs( 4)= 0.2682887768100528D+01
	xs( 5)= 0.3721464681096727D+01
	xs( 6)= 0.4866230500274611D+01
	xs( 7)= 0.6130031933550999D+01
	xs( 8)= 0.7518409709507797D+01
	xs( 9)= 0.9034313653276023D+01
	xs(10)= 0.1066990121177804D+02
	xs(11)= 0.1244107808709449D+02
c
	return
c
c	Data for frequency k in region [ 9.42478,15.70796]
c
 1010 continue
	nlam=10
c
	xs( 1)= 0.2310735429121635D-01
	xs( 2)= 0.1280445527700792D+01
	xs( 3)= 0.2608421346890076D+01
	xs( 4)= 0.4002775470541412D+01
	xs( 5)= 0.5463876102045486D+01
	xs( 6)= 0.6992086265390597D+01
	xs( 7)= 0.8585642091493785D+01
	xs( 8)= 0.1027470294077428D+02
	xs( 9)= 0.1201654878189594D+02
	xs(10)= 0.1384643074122658D+02
c
	return
c
	end
c
c
c
c
c
	subroutine yukq2d_tab9(ier,rky,nlam,xs)
	implicit real *8 (a-h,o-z)
c
c Purpose:
c	This is the quadrature table (nodes only) for yukawa 2D
c	with 9 digits of accuracy. It covers only the range
c	(0,8*pi]!
c
c Inputs:
c	rky: the frequency=rky*i
c
c Outputs:
c	ier: error output code
c		=0 no error
c		=8 rky out of range of (0,\infty)
c	nlam: number of nodes
c	xs(nlam): the quad nodes
c
	dimension xs(*)
c
	ier=0
	done=1
	pi=4*atan(done)
c
c ... If rky lies outside of (0,\infty), BOMB!
c
cccc	if ((rky.le.0.0d0).or.(rky.gt.8*pi)) then
	if (rky.lt.0.0d0) then
	   ier=8
	   return
	endif
c
c ... Choose the correct the nodes:
c
	if ((rky.ge.0.0000d0).and.(rky.le.0.0001d0))   goto 1001
	if ((rky.gt.0.0001d0).and.(rky.le.0.0005d0))   goto 1002
	if ((rky.gt.0.0005d0).and.(rky.le.0.0020d0))   goto 1003
	if ((rky.gt.0.0020d0).and.(rky.le.0.0060d0))   goto 1004
	if ((rky.gt.0.0060d0).and.(rky.le.0.0150d0))   goto 1005
	if ((rky.gt.0.0150d0).and.(rky.le.0.0400d0))   goto 1006
	if ((rky.gt.0.0400d0).and.(rky.le.0.1000d0))   goto 1007
	if ((rky.gt.0.1000d0).and.(rky.le.0.3000d0))   goto 1008
	if ((rky.gt.0.3000d0).and.(rky.le.0.6000d0))   goto 1009
	if ((rky.gt.0.6000d0).and.(rky.le.1.0000d0))   goto 1010
	if ((rky.gt.1.0000d0).and.(rky.le.1.6000d0))   goto 1011
	if ((rky.gt.1.6000d0).and.(rky.le.2.4000d0))   goto 1012
	if ((rky.gt.2.4000d0).and.(rky.le.3.5000d0))   goto 1013
	if ((rky.gt.3.5000d0).and.(rky.le.5.0000d0))   goto 1014
	if ((rky.gt.5.0000d0).and.(rky.le.4.0d0*pi))   goto 1015
	if (rky.gt.4.0d0*pi)   goto 1016
cccc	if ((rky.gt.4.0d0*pi).and.(rky.le.8.0d0*pi))   goto 1016
c
c ... If for some reason, the code gets here, bomb!
c
	ier=8
	return
c
c	Data for frequency k in region [ 0.00000, 0.00010]
c
 1001 continue
	nlam=21
c
	xs( 1)= 0.1535117952755627d-05
	xs( 2)= 0.2622863114325469d-03
	xs( 3)= 0.4031193642404318d-01
	xs( 4)= 0.2048866890324526d+00
	xs( 5)= 0.4954951425841827d+00
	xs( 6)= 0.9046404661603660d+00
	xs( 7)= 0.1424151641748711d+01
	xs( 8)= 0.2046480915172257d+01
	xs( 9)= 0.2765284141800585d+01
	xs(10)= 0.3575577626029109d+01
	xs(11)= 0.4473737056250865d+01
	xs(12)= 0.5457498272316123d+01
	xs(13)= 0.6526001244221725d+01
	xs(14)= 0.7679780299012325d+01
	xs(15)= 0.8920466403627515d+01
	xs(16)= 0.1025010042590890d+02
	xs(17)= 0.1167039697375212d+02
	xs(18)= 0.1318254668780636d+02
	xs(19)= 0.1478707098181131d+02
	xs(20)= 0.1649310313440046d+02
	xs(21)= 0.1823354313957236d+02
c
	return
c
c	Data for frequency k in region [ 0.00010, 0.00050]
c
 1002 continue
	nlam=22
c
	xs( 1)= 0.4495141737947961d-05
	xs( 2)= 0.8271912007167259d-03
	xs( 3)= 0.4065494955817783d-01
	xs( 4)= 0.1973317751830841d+00
	xs( 5)= 0.4741592617361157d+00
	xs( 6)= 0.8644228955678983d+00
	xs( 7)= 0.1360712182311516d+01
	xs( 8)= 0.1956143185030555d+01
	xs( 9)= 0.2644899107198206d+01
	xs(10)= 0.3422372398735057d+01
	xs(11)= 0.4285151062989826d+01
	xs(12)= 0.5231013699226430d+01
	xs(13)= 0.6258997969810459d+01
	xs(14)= 0.7369494229891185d+01
	xs(15)= 0.8564173040036490d+01
	xs(16)= 0.9845530382915474d+01
	xs(17)= 0.1121618564606059d+02
	xs(18)= 0.1267841714850415d+02
	xs(19)= 0.1423511515892654d+02
	xs(20)= 0.1588974229162217d+02
	xs(21)= 0.1764243526549907d+02
	xs(22)= 0.1964596076586591d+02
c
	return
c
c	Data for frequency k in region [ 0.00050, 0.00200]
c
 1003 continue
	nlam=23
c
	xs( 1)= 0.1157629729675591d-04
	xs( 2)= 0.8841257579703665d-03
	xs( 3)= 0.4747890944948097d-02
	xs( 4)= 0.5105066805344549d-01
	xs( 5)= 0.2090220456892311d+00
	xs( 6)= 0.4870113541738093d+00
	xs( 7)= 0.8789190709063099d+00
	xs( 8)= 0.1377408183841528d+01
	xs( 9)= 0.1975530973526379d+01
	xs(10)= 0.2667337862741459d+01
	xs(11)= 0.3448070401884777d+01
	xs(12)= 0.4314185714687142d+01
	xs(13)= 0.5263363128605803d+01
	xs(14)= 0.6294557164602054d+01
	xs(15)= 0.7408065823168585d+01
	xs(16)= 0.8605445203561548d+01
	xs(17)= 0.9889060711835754d+01
	xs(18)= 0.1126137796194874d+02
	xs(19)= 0.1272455644859770d+02
	xs(20)= 0.1428111085556474d+02
	xs(21)= 0.1593432717161145d+02
	xs(22)= 0.1768821029389825d+02
	xs(23)= 0.1965301697989015d+02
c
	return
c
c	Data for frequency k in region [ 0.00200, 0.00600]
c
 1004 continue
	nlam=23
c
	xs( 1)= 0.3917914822793023d-04
	xs( 2)= 0.2833949658604951d-02
	xs( 3)= 0.1189598864391428d-01
	xs( 4)= 0.6584770675986817d-01
	xs( 5)= 0.2250649762565295d+00
	xs( 6)= 0.5022355276055066d+00
	xs( 7)= 0.8924860313249248d+00
	xs( 8)= 0.1388747337919430d+01
	xs( 9)= 0.1984210329945824d+01
	xs(10)= 0.2673046108909887d+01
	xs(11)= 0.3450605781214922d+01
	xs(12)= 0.4313435766968134d+01
	xs(13)= 0.5259285557232552d+01
	xs(14)= 0.6287175703697724d+01
	xs(15)= 0.7397485603926313d+01
	xs(16)= 0.8591881600222266d+01
	xs(17)= 0.9872865910812870d+01
	xs(18)= 0.1124305365462464d+02
	xs(19)= 0.1270475263837819d+02
	xs(20)= 0.1426069896160881d+02
	xs(21)= 0.1591420223421873d+02
	xs(22)= 0.1766909996193931d+02
	xs(23)= 0.1965759708019235d+02
c
	return
c
c	Data for frequency k in region [ 0.00600, 0.01500]
c
 1005 continue
	nlam=23
c
	xs( 1)= 0.1030476357577470d-03
	xs( 2)= 0.7140493181793772d-02
	xs( 3)= 0.2466911262313332d-01
	xs( 4)= 0.8933541075783680d-01
	xs( 5)= 0.2528786184798655d+00
	xs( 6)= 0.5317525804492682d+00
	xs( 7)= 0.9232221878652851d+00
	xs( 8)= 0.1420814305614044d+01
	xs( 9)= 0.2017932389874048d+01
	xs(10)= 0.2708820506465297d+01
	xs(11)= 0.3488820240276372d+01
	xs(12)= 0.4354403887106134d+01
	xs(13)= 0.5303182969952125d+01
	xs(14)= 0.6333965609024187d+01
	xs(15)= 0.7446833481138459d+01
	xs(16)= 0.8643085275134609d+01
	xs(17)= 0.9924861375752972d+01
	xs(18)= 0.1129453044290753d+02
	xs(19)= 0.1275433356609650d+02
	xs(20)= 0.1430703735824959d+02
	xs(21)= 0.1595575200560837d+02
	xs(22)= 0.1770707287710088d+02
	xs(23)= 0.1966764039304024d+02
c
	return
c
c	Data for frequency k in region [ 0.01500, 0.04000]
c
 1006 continue
	nlam=23
c
	xs( 1)= 0.2082948065549317d-03
	xs( 2)= 0.1349386971485700d-01
	xs( 3)= 0.3666176009605593d-01
	xs( 4)= 0.9194440160888462d-01
	xs( 5)= 0.2221003206373631d+00
	xs( 6)= 0.4602553726052285d+00
	xs( 7)= 0.8143191585786944d+00
	xs( 8)= 0.1280487489418999d+01
	xs( 9)= 0.1852130413933025d+01
	xs(10)= 0.2522748726217873d+01
	xs(11)= 0.3286866214209375d+01
	xs(12)= 0.4140263192394123d+01
	xs(13)= 0.5080038139999754d+01
	xs(14)= 0.6104670834932914d+01
	xs(15)= 0.7214090293585194d+01
	xs(16)= 0.8409588918288577d+01
	xs(17)= 0.9693372577557398d+01
	xs(18)= 0.1106785999899846d+02
	xs(19)= 0.1253525631269538d+02
	xs(20)= 0.1409825262054152d+02
	xs(21)= 0.1576156844030292d+02
	xs(22)= 0.1751929671551012d+02
	xs(23)= 0.1954563204712521d+02
c
	return
c
c	Data for frequency k in region [ 0.04000, 0.10000]
c
 1007 continue
	nlam=23
c
	xs( 1)= 0.4857380606431150d-03
	xs( 2)= 0.3050761140890756d-01
	xs( 3)= 0.7555715958671883d-01
	xs( 4)= 0.1572999181995804d+00
	xs( 5)= 0.3077081670720290d+00
	xs( 6)= 0.5540468760444206d+00
	xs( 7)= 0.9078116352334256d+00
	xs( 8)= 0.1369250553183225d+01
	xs( 9)= 0.1934105038627219d+01
	xs(10)= 0.2597130961459275d+01
	xs(11)= 0.3353509334142451d+01
	xs(12)= 0.4199357699693351d+01
	xs(13)= 0.5131934664462783d+01
	xs(14)= 0.6149782752980638d+01
	xs(15)= 0.7252851298602842d+01
	xs(16)= 0.8442454529240202d+01
	xs(17)= 0.9720838568126080d+01
	xs(18)= 0.1109045540360206d+02
	xs(19)= 0.1255349744437288d+02
	xs(20)= 0.1411252514716402d+02
	xs(21)= 0.1577203103335804d+02
	xs(22)= 0.1752744538374209d+02
	xs(23)= 0.1953781239370574d+02
c
	return
c
c	Data for frequency k in region [ 0.10000, 0.30000]
c
 1008 continue
	nlam=23
c
	xs( 1)= 0.1058264882679794d-02
	xs( 2)= 0.6438477024321543d-01
	xs( 3)= 0.1468043520804301d+00
	xs( 4)= 0.2658807812307487d+00
	xs( 5)= 0.4433539407813392d+00
	xs( 6)= 0.7016290095690536d+00
	xs( 7)= 0.1056645987868365d+01
	xs( 8)= 0.1514847620850635d+01
	xs( 9)= 0.2075730952914228d+01
	xs(10)= 0.2735556981228925d+01
	xs(11)= 0.3489803516430509d+01
	xs(12)= 0.4334380377645971d+01
	xs(13)= 0.5266184948045940d+01
	xs(14)= 0.6283394167136017d+01
	xs(15)= 0.7385633488090786d+01
	xs(16)= 0.8573932138937010d+01
	xs(17)= 0.9850286297227814d+01
	xs(18)= 0.1121695163563404d+02
	xs(19)= 0.1267601359058763d+02
	xs(20)= 0.1423004981508835d+02
	xs(21)= 0.1588281245101230d+02
	xs(22)= 0.1763594441160455d+02
	xs(23)= 0.1962129279152516d+02
c
	return
c
c	Data for frequency k in region [ 0.30000, 0.60000]
c
 1009 continue
	nlam=23
c
	xs( 1)= 0.2366336560417182d-02
	xs( 2)= 0.1398252914282345d+00
	xs( 3)= 0.3019178675621816d+00
	xs( 4)= 0.5077007380067045d+00
	xs( 5)= 0.7767810609572336d+00
	xs( 6)= 0.1125118595850669d+01
	xs( 7)= 0.1562854782446450d+01
	xs( 8)= 0.2094162755747121d+01
	xs( 9)= 0.2718978154147610d+01
	xs(10)= 0.3435048082134710d+01
	xs(11)= 0.4239423797941082d+01
	xs(12)= 0.5129359151267305d+01
	xs(13)= 0.6102839964916704d+01
	xs(14)= 0.7158920551244979d+01
	xs(15)= 0.8297897084033053d+01
	xs(16)= 0.9521189349463359d+01
	xs(17)= 0.1083084227088629d+02
	xs(18)= 0.1222893116201884d+02
	xs(19)= 0.1371746783854318d+02
	xs(20)= 0.1529804055573770d+02
	xs(21)= 0.1698256757655113d+02
	xs(22)= 0.1872704834894280d+02
	xs(23)= 0.2090472261725520d+02
c
	return
c
c	Data for frequency k in region [ 0.60000, 1.00000]
c
 1010 continue
	nlam=21
c
	xs( 1)= 0.3975152821883654d-02
	xs( 2)= 0.2313328560964649d+00
	xs( 3)= 0.4877615722850273d+00
	xs( 4)= 0.7928915445367970d+00
	xs( 5)= 0.1166660124900950d+01
	xs( 6)= 0.1623564833296637d+01
	xs( 7)= 0.2171784795984376d+01
	xs( 8)= 0.2814162378735475d+01
	xs( 9)= 0.3550180759429274d+01
	xs(10)= 0.4377810899547132d+01
	xs(11)= 0.5294819769797666d+01
	xs(12)= 0.6299617283708436d+01
	xs(13)= 0.7391776429061498d+01
	xs(14)= 0.8572188354450002d+01
	xs(15)= 0.9842696986917147d+01
	xs(16)= 0.1120536747996848d+02
	xs(17)= 0.1266202959103300d+02
	xs(18)= 0.1421489414826403d+02
	xs(19)= 0.1586767877416926d+02
	xs(20)= 0.1762191164017034d+02
	xs(21)= 0.1960260629081949d+02
c
	return
c
c	Data for frequency k in region [ 1.00000, 1.60000]
c
 1011 continue
	nlam=21
c
	xs( 1)= 0.5673727826685848d-02
	xs( 2)= 0.3263499838896848d+00
	xs( 3)= 0.6758475974764977d+00
	xs( 4)= 0.1069764003572374d+01
	xs( 5)= 0.1525436065189140d+01
	xs( 6)= 0.2055602573342085d+01
	xs( 7)= 0.2668182727746665d+01
	xs( 8)= 0.3366938305943169d+01
	xs( 9)= 0.4152731139136073d+01
	xs(10)= 0.5024873004998220d+01
	xs(11)= 0.5982255539112845d+01
	xs(12)= 0.7024179700712864d+01
	xs(13)= 0.8150863904981588d+01
	xs(14)= 0.9363519551917083d+01
	xs(15)= 0.1066391982256505d+02
	xs(16)= 0.1205377595556993d+02
	xs(17)= 0.1353454904631846d+02
	xs(18)= 0.1510805547742793d+02
	xs(19)= 0.1677583734170308d+02
	xs(20)= 0.1857065071659749d+02
	xs(21)= 0.2043047053701541d+02
c
	return
c
c	Data for frequency k in region [ 1.60000, 2.40000]
c
 1012 continue
	nlam=20
c
	xs( 1)= 0.7908470846459181d-02
	xs( 2)= 0.4507060505237845d+00
	xs( 3)= 0.9238376160177388d+00
	xs( 4)= 0.1439666618266866d+01
	xs( 5)= 0.2013912708684281d+01
	xs( 6)= 0.2657618717023279d+01
	xs( 7)= 0.3377589600470476d+01
	xs( 8)= 0.4177383366585843d+01
	xs( 9)= 0.5058489903850656d+01
	xs(10)= 0.6021442994505549d+01
	xs(11)= 0.7066732716700908d+01
	xs(12)= 0.8195421117175602d+01
	xs(13)= 0.9409289675326125d+01
	xs(14)= 0.1071043297432503d+02
	xs(15)= 0.1210065946198504d+02
	xs(16)= 0.1358136567970670d+02
	xs(17)= 0.1515413090154343d+02
	xs(18)= 0.1682033177615627d+02
	xs(19)= 0.1861245494464274d+02
	xs(20)= 0.2046084250477121d+02
c
	return
c
c	Data for frequency k in region [ 2.40000, 3.50000]
c
 1013 continue
	nlam=19
c
	xs( 1)= 0.1031362550703108d-01
	xs( 2)= 0.5832764675637385d+00
	xs( 3)= 0.1186948263684361d+01
	xs( 4)= 0.1829610856970760d+01
	xs( 5)= 0.2525568887751962d+01
	xs( 6)= 0.3284960933316992d+01
	xs( 7)= 0.4114331154066949d+01
	xs( 8)= 0.5017651539709050d+01
	xs( 9)= 0.5997415309644477d+01
	xs(10)= 0.7055597853803654d+01
	xs(11)= 0.8194329847158318d+01
	xs(12)= 0.9416077959283299d+01
	xs(13)= 0.1072327941954538d+02
	xs(14)= 0.1211783848652962d+02
	xs(15)= 0.1360111481769163d+02
	xs(16)= 0.1517457542725795d+02
	xs(17)= 0.1683895304758333d+02
	xs(18)= 0.1862937738025522d+02
	xs(19)= 0.2045445331903558d+02
c
	return
c
c	Data for frequency k in region [ 3.50000, 5.00000]
c
 1014 continue
	nlam=18
c
	xs( 1)= 0.1282200522414101d-01
	xs( 2)= 0.7202840270317452d+00
	xs( 3)= 0.1457728471942271d+01
	xs( 4)= 0.2228538893809137d+01
	xs( 5)= 0.3045592986105230d+01
	xs( 6)= 0.3918676012033906d+01
	xs( 7)= 0.4854866525377087d+01
	xs( 8)= 0.5859349115923688d+01
	xs( 9)= 0.6936307760576089d+01
	xs(10)= 0.8089590157093905d+01
	xs(11)= 0.9322839709665509d+01
	xs(12)= 0.1063909148600264d+02
	xs(13)= 0.1204034708698820d+02
	xs(14)= 0.1352771284705930d+02
	xs(15)= 0.1510233686657138d+02
	xs(16)= 0.1676210566488524d+02
	xs(17)= 0.1854952756943260d+02
	xs(18)= 0.2032525846463969d+02
c
	return
c
c	Data for frequency k in region [ 5.00000,4*pi]
c
 1015 continue
	nlam=16
c
	xs( 1)= 0.1606166986797319d-01
	xs( 2)= 0.8960773231289210d+00
	xs( 3)= 0.1807331516180948d+01
	xs( 4)= 0.2748749656275738d+01
	xs( 5)= 0.3732826729981642d+01
	xs( 6)= 0.4770103605157288d+01
	xs( 7)= 0.5869533363530103d+01
	xs( 8)= 0.7038994914131440d+01
	xs( 9)= 0.8285307637591547d+01
	xs(10)= 0.9613710275920965d+01
	xs(11)= 0.1102757246600709d+02
	xs(12)= 0.1252894438151515d+02
	xs(13)= 0.1412012173244783d+02
	xs(14)= 0.1580577030643560d+02
	xs(15)= 0.1758346459227530d+02
	xs(16)= 0.1961209382068556d+02
c
	return
c
c	Data for frequency k in region [4*pi,8*pi]
c
 1016 continue
	nlam=10
c
	xs( 1)= 0.2835582203611153d-01
	xs( 2)= 0.1547977664590988d+01
	xs( 3)= 0.3105649690558764d+01
	xs( 4)= 0.4683805342592347d+01
	xs( 5)= 0.6292155650485402d+01
	xs( 6)= 0.7939200773539802d+01
	xs( 7)= 0.9636934977159672d+01
	xs( 8)= 0.1139548968636039d+02
	xs( 9)= 0.1323099615974295d+02
	xs(10)= 0.1516035876055231d+02
c
	return
c
	end
c
c
c
c
c
c
c
c
	subroutine yukq2d_tab12(ier,rky,nlam,xs,ws)
	implicit real *8 (a-h,o-z)
c
c Purpose:
c	This routine returns a set of Gaussian nodes and weights for
c	integrating the functions:
c
c	      exp(-x*\sqrt(t**2 - rk**2))* cos(t*y)
c           -------------------------------------------
c		      \sqrt(t**2 - rk**2)
c
c	over the range t=0 to t=infinity. Here, rk=eye*rky
c	is pure imaginary, thus it is for the Yukawa kernal.
c
c 	They work only for values of x within the range [1,4] and 
c	y within the range [0,4].
c
c	These tables are for 12 digits of accuracy.
c
c Inputs:
c	rky: the frequency rk = rky*eye.
c
c Outputs:
c	ier: error output code
c		=0 no error
c		=8 rky out of range of [0.5d-6,\infty)
c	nlam: number of nodes
c	xs(nlam): the quad nodes
c	ws(nlam): the quad weights
c
	dimension xs(1),ws(1)
c
        dimension xs01( 3),ws01( 3),xs02( 4),ws02( 4),xs03( 6),ws03( 6),
     1            xs04( 7),ws04( 7),xs05( 8),ws05( 8),xs06( 8),ws06( 8),
     2            xs07( 9),ws07( 9),xs08(10),ws08(10),xs09(10),ws09(10),
     3            xs10(11),ws10(11),xs11(11),ws11(11),xs12(12),ws12(12),
     4            xs13(13),ws13(13),xs14(13),ws14(13),xs15(14),ws15(14),
     5            xs16(15),ws16(15),xs17(16),ws17(16),xs18(17),ws18(17),
     6            xs19(18),ws19(18),xs20(18),ws20(18),xs21(19),ws21(19),
     7            xs22(20),ws22(20),xs23(21),ws23(21),xs24(23),ws24(23),
     8            xs25(29),ws25(29),xs26(33),ws26(33),xs27(34),ws27(34),
     9            xs28(34),ws28(34),xs29(37),ws29(37),xs30(31),ws30(31),
     a            xs31(29),ws31(29)
c
        data xs01/
     1      0.4337271987788022D+00,0.1893878982463956D+01,
     2      0.3477834350051434D+01/
        data ws01/
     1      0.3915984315117753D+00,0.6339930334120912D+00,
     2      0.8589448414096248D+00/
c
        data xs02/
     1      0.9002078213524284D+00,0.2704207795480901D+01,
     2      0.4511567499416862D+01,0.6297736645624357D+01/
        data ws02/
     1      0.5948313840403550D+00,0.6024566175102140D+00,
     2      0.6243505788914511D+00,0.7535503823257744D+00/
c
        data xs03/
     1      0.2241315411848001D+00,0.1823147002626775D+01,
     2      0.3609840652245423D+01,0.5419679926600517D+01,
     3      0.7252439536267391D+01,0.9099696003425382D+01/
        data ws03/
     1      0.3025301833129270D+00,0.5671508570133170D+00,
     2      0.5790992539855463D+00,0.5898630639230499D+00,
     3      0.6114217697814978D+00,0.7056708152581415D+00/
c
        data xs04/
     1      0.8693795688834728D+00,0.2611236816302672D+01,
     2      0.4361656407886277D+01,0.6126385931777357D+01,
     3      0.7909868665727510D+01,0.9711710711389840D+01,
     4      0.1147256819664432D+02/
        data ws04/
     1      0.5545628616214894D+00,0.5566208655613182D+00,
     2      0.5605313998895761D+00,0.5670129296638948D+00,
     3      0.5770179162865163D+00,0.5936805673405133D+00,
     4      0.6191319688703032D+00/
c
        data xs05/
     1      0.8603913934970054D+00,0.2584171032143994D+01,
     2      0.4317424641985157D+01,0.6066193850217365D+01,
     3      0.7836798015616511D+01,0.9633857496848266D+01,
     4      0.1146311690522769D+02,0.1328317405378959D+02/
        data ws05/
     1      0.5480803170242237D+00,0.5501259376138710D+00,
     2      0.5542971131460643D+00,0.5605327974947940D+00,
     3      0.5689232745325713D+00,0.5801879371228671D+00,
     4      0.5976280020148437D+00,0.6291907593875732D+00/
c
        data xs06/
     1      0.8627840587172648D+00,0.2591851839954139D+01,
     2      0.4331864677170267D+01,0.6090080547435727D+01,
     3      0.7874020481892200D+01,0.9692021505225338D+01,
     4      0.1155241717579588D+02,0.1346272077688431D+02/
        data ws06/
     1      0.5495881083717208D+00,0.5518939478105152D+00,
     2      0.5565735379926828D+00,0.5637609842337140D+00,
     3      0.5736107834182939D+00,0.5872002147895956D+00,
     4      0.6083940342301878D+00,0.6805984629765324D+00/
c
        data xs07/
     1      0.8425041879324330D+00,0.2531027465168362D+01,
     2      0.4230099020586819D+01,0.5946745460929446D+01,
     3      0.7688141780758682D+01,0.9461745644338436D+01,
     4      0.1127632205297322D+02,0.1314151978952680D+02,
     5      0.1507163729454065D+02/
        data ws07/
     1      0.5365671959441241D+00,0.5388178378737477D+00,
     2      0.5432870546756012D+00,0.5500583709963521D+00,
     3      0.5591982449844027D+00,0.5709699417249471D+00,
     4      0.5866325911494517D+00,0.6093904787577754D+00,
     5      0.6795036879851589D+00/
c
        data xs08/
     1      0.7309691642975755D+00,0.2215031739367663D+01,
     2      0.3748687850079632D+01,0.5335402407307179D+01,
     3      0.6970890963350544D+01,0.8652781637607674D+01,
     4      0.1038226953071824D+02,0.1216364975176876D+02,
     5      0.1400363964669505D+02,0.1590917538784610D+02/
        data ws08/
     1      0.4666585537957219D+00,0.4796231841028948D+00,
     2      0.4968255008238749D+00,0.5130777132300673D+00,
     3      0.5280615353427175D+00,0.5428931043359737D+00,
     4      0.5587873531397004D+00,0.5771233462121232D+00,
     5      0.6011252112866533D+00,0.6621187269777797D+00/
c
        data xs09/
     1      0.8128755258384842D+00,0.2442452627066455D+01,
     2      0.4083459016310023D+01,0.5743381824795795D+01,
     3      0.7429616511825045D+01,0.9149595712895785D+01,
     4      0.1091126260137003D+02,0.1272360437483450D+02,
     5      0.1459866599010113D+02,0.1655131000029022D+02/
        data ws09/
     1      0.5176999997090170D+00,0.5201325201532916D+00,
     2      0.5249704408213474D+00,0.5321772694420557D+00,
     3      0.5417393396666967D+00,0.5537537933359328D+00,
     4      0.5684858546273517D+00,0.5869820157803626D+00,
     5      0.6128249101909814D+00,0.6783029216947134D+00/
c
        data xs10/
     1      0.4397365679118569D+00,0.1742790044037754D+01,
     2      0.3280045124290629D+01,0.4874425433049677D+01,
     3      0.6504524780895325D+01,0.8170532430888274D+01,
     4      0.9877215416682777D+01,0.1163114424681088D+02,
     5      0.1344067382286605D+02,0.1531711576557563D+02,
     6      0.1727765539299227D+02/
        data ws10/
     1      0.3234185402490417D+00,0.4715186488049498D+00,
     2      0.5008501027658292D+00,0.5134148804109376D+00,
     3      0.5244126613396737D+00,0.5364843136397529D+00,
     4      0.5504245623771843D+00,0.5667796541621446D+00,
     5      0.5865337155153481D+00,0.6132630281219112D+00,
     6      0.6776002133237198D+00/
c
        data xs11/
     1      0.7799234067553803D+00,0.2343973413038102D+01,
     2      0.3920512098201428D+01,0.5517505079708476D+01,
     3      0.7142497193657643D+01,0.8802660607900451D+01,
     4      0.1050504855828448D+02,0.1225711571473174D+02,
     5      0.1406769663459437D+02,0.1594906336457144D+02,
     6      0.1792023596576152D+02/
        data ws11/
     1      0.4967389451418022D+00,0.4994077846364660D+00,
     2      0.5046722577127967D+00,0.5124075497384656D+00,
     3      0.5224786087922808D+00,0.5347980872395960D+00,
     4      0.5494015952140132D+00,0.5665786777774623D+00,
     5      0.5872759203286858D+00,0.6150100314923930D+00,
     6      0.6801996788287310D+00/
c
        data xs12/
     1      0.7533365392404356D+00,0.2264374746918333D+01,
     2      0.3788349160953573D+01,0.5333393550370122D+01,
     3      0.6907038246729960D+01,0.8516200837897051D+01,
     4      0.1016735545192141D+02,0.1186692253621339D+02,
     5      0.1362192771952353D+02,0.1544116227148668D+02,
     6      0.1733767188194156D+02,0.1933497676241484D+02/
        data ws12/
     1      0.4798214612481737D+00,0.4825916980140330D+00,
     2      0.4880331951866413D+00,0.4959727875864585D+00,
     3      0.5062075183779131D+00,0.5185610088126115D+00,
     4      0.5329401194865709D+00,0.5494201216438677D+00,
     5      0.5683649673190375D+00,0.5908427567049604D+00,
     6      0.6204892882514176D+00,0.6870648741503630D+00/
c
        data xs13/
     1      0.6389094945657234D+00,0.1941985499119044D+01,
     2      0.3299790694857624D+01,0.4714272543993373D+01,
     3      0.6180517190602493D+01,0.7696149523662114D+01,
     4      0.9261709230901543D+01,0.1087954219673395D+02,
     5      0.1255321315721660D+02,0.1428758118289613D+02,
     6      0.1608947616316234D+02,0.1796919888909961D+02,
     7      0.1993785460281601D+02/
        data ws13/
     1      0.4082481082338956D+00,0.4228789219592761D+00,
     2      0.4414852330767394D+00,0.4587015092483641D+00,
     3      0.4746166317055216D+00,0.4903006480890068D+00,
     4      0.5064941213864286D+00,0.5236436218110074D+00,
     5      0.5421142588687967D+00,0.5624202747387799D+00,
     6      0.5856095142845739D+00,0.6145329191152086D+00,
     7      0.6668852408904583D+00/
c
        data xs14/
     1      0.6330332204631226D+00,0.1920277889109784D+01,
     2      0.3256016783517610D+01,0.4646072511819507D+01,
     3      0.6088977285299316D+01,0.7584049877258622D+01,
     4      0.9132498699010236D+01,0.1073687602156387D+02,
     5      0.1240078638482328D+02,0.1412915269113543D+02,
     6      0.1592919136825324D+02,0.1781257734213308D+02,
     7      0.1980152174493578D+02/
        data ws14/
     1      0.4042385413518549D+00,0.4167130935806414D+00,
     2      0.4338518223139906D+00,0.4509692861122436D+00,
     3      0.4675749929832744D+00,0.4842860972705150D+00,
     4      0.5016248772907705D+00,0.5199454861671128D+00,
     5      0.5395986335354940D+00,0.5611368549908327D+00,
     6      0.5857138488376861D+00,0.6165892609645446D+00,
     7      0.6811198030101050D+00/
c
        data xs15/
     1      0.6075221729325599D+00,0.1829209635891240D+01,
     2      0.3072228879243157D+01,0.4352438109762610D+01,
     3      0.5682881446911940D+01,0.7071369758096086D+01,
     4      0.8521596637102608D+01,0.1003536937966244D+02,
     5      0.1161420795769320D+02,0.1326029406691335D+02,
     6      0.1497728997211866D+02,0.1677151978248251D+02,
     7      0.1865441331514260D+02,0.2064929701574674D+02/
        data ws15/
     1      0.3871022807818025D+00,0.3914279399294173D+00,
     2      0.4007711595277486D+00,0.4149415787256892D+00,
     3      0.4324506210220707D+00,0.4516713514637193D+00,
     4      0.4716594908354853D+00,0.4921145750765137D+00,
     5      0.5131183524376622D+00,0.5350092229862505D+00,
     6      0.5584320524631859D+00,0.5846730276170546D+00,
     7      0.6169841689418908D+00,0.6821884133527737D+00/
c
        data xs16/
     1      0.5731456479358769D+00,0.1725110314259404D+01,
     2      0.2894235416237436D+01,0.4092580963018630D+01,
     3      0.5332711633160873D+01,0.6626147324837460D+01,
     4      0.7981390058042978D+01,0.9403399001490470D+01,
     5      0.1089461431609754D+02,0.1245644157595099D+02,
     6      0.1409054038089652D+02,0.1579995585129028D+02,
     7      0.1759053128643777D+02,0.1947347377891091D+02,
     8      0.2147243929646323D+02/
        data ws16/
     1      0.3651764295577102D+00,0.3687935724382133D+00,
     2      0.3761334804042009D+00,0.3874328466833262D+00,
     3      0.4026824060477422D+00,0.4211980498694760D+00,
     4      0.4418328670656090D+00,0.4635673128222583D+00,
     5      0.4858379564463574D+00,0.5085358313432928D+00,
     6      0.5319228943553758D+00,0.5566461304305450D+00,
     7      0.5840234909325783D+00,0.6173400418618044D+00,
     8      0.6822547548809375D+00/
c
        data xs17/
     1      0.5395806262399773D+00,0.1624400439980942D+01,
     2      0.2726080205914224D+01,0.3855645565341019D+01,
     3      0.5024075108494234D+01,0.6242050972226884D+01,
     4      0.7518936998885530D+01,0.8861628743272105D+01,
     5      0.1027426518842153D+02,0.1175893796495583D+02,
     6      0.1331684183580449D+02,0.1494943132498457D+02,
     7      0.1665954097819745D+02,0.1845291441870933D+02,
     8      0.2034072966012060D+02,0.2234914749643477D+02/
        data ws17/
     1      0.3438086827394166D+00,0.3474037024106396D+00,
     2      0.3545311213903358D+00,0.3651557969366697D+00,
     3      0.3792622955989329D+00,0.3966290658679965D+00,
     4      0.4166292240329886D+00,0.4383701640734485D+00,
     5      0.4610466847141394D+00,0.4841828567235468D+00,
     6      0.5076812491488425D+00,0.5318038975014301D+00,
     7      0.5572014918673597D+00,0.5851788052016255D+00,
     8      0.6192384825748790D+00,0.6865019896998276D+00/
c
        data xs18/
     1      0.2129807033618319D+00,0.1069405147330353D+01,
     2      0.2077202968452943D+01,0.3123171086676181D+01,
     3      0.4206837415010237D+01,0.5335857934572122D+01,
     4      0.6519110958960924D+01,0.7765090427697869D+01,
     5      0.9080568315490748D+01,0.1046982788423354D+02,
     6      0.1193496302597094D+02,0.1347691983908166D+02,
     7      0.1509673429917263D+02,0.1679678276996430D+02,
     8      0.1858234667714240D+02,0.2046442777248248D+02,
     9      0.2246706920004517D+02/
        data ws18/
     1      0.1866658662757055D+00,0.3117321044383928D+00,
     2      0.3274220116351086D+00,0.3386084809184605D+00,
     3      0.3517001439966007D+00,0.3675344199431456D+00,
     4      0.3862062249204605D+00,0.4073727659227647D+00,
     5      0.4303117423265594D+00,0.4542234714955994D+00,
     6      0.4785541574502829D+00,0.5031348310379862D+00,
     7      0.5281921234612573D+00,0.5543778372537952D+00,
     8      0.5830772552095165D+00,0.6178315428930617D+00,
     9      0.6843756811738616D+00/
c
        data xs19/
     1      0.4649050282328557D+00,0.1400708222445862D+01,
     2      0.2354163461009290D+01,0.3336049274472249D+01,
     3      0.4355968782092484D+01,0.5422366568163757D+01,
     4      0.6542793161356966D+01,0.7724160962940918D+01,
     5      0.8972617456937888D+01,0.1029300772837646D+02,
     6      0.1168844483952215D+02,0.1316047924767775D+02,
     7      0.1470984624011356D+02,0.1633748767236883D+02,
     8      0.1804571714535219D+02,0.1983982940707949D+02,
     9      0.2173052626225473D+02,0.2374082515291999D+02/
        data ws19/
     1      0.2962874999024707D+00,0.3000849011758663D+00,
     2      0.3074769772425589D+00,0.3181202825759740D+00,
     3      0.3316270498992511D+00,0.3476617120486491D+00,
     4      0.3659913273438296D+00,0.3864202585597354D+00,
     5      0.4086316052041948D+00,0.4321191395252142D+00,
     6      0.4563217579194220D+00,0.4808378785980597D+00,
     7      0.5055646282661518D+00,0.5307427450669253D+00,
     8      0.5570270655795977D+00,0.5858435936872881D+00,
     9      0.6205828274530917D+00,0.6847864761457598D+00/
c
        data xs20/
     1      0.4310131873815112D+00,0.1300075011973338D+01,
     2      0.2189744871866963D+01,0.3112260400484002D+01,
     3      0.4078087135733026D+01,0.5095938313442746D+01,
     4      0.6173146871386341D+01,0.7316133292561415D+01,
     5      0.8530567254761603D+01,0.9821025097820341D+01,
     6      0.1119052722986811D+02,0.1264054899300905D+02,
     7      0.1417166466334060D+02,0.1578458101515267D+02,
     8      0.1748137329018582D+02,0.1926712735816731D+02,
     9      0.2115279696035541D+02,0.2316433354812553D+02/
        data ws20/
     1      0.2747675356725984D+00,0.2792190892709974D+00,
     2      0.2878139166071109D+00,0.3000327336276316D+00,
     3      0.3152917621361939D+00,0.3330786104985851D+00,
     4      0.3530331845483589D+00,0.3749175301238786D+00,
     5      0.3984609451126882D+00,0.4232303228383628D+00,
     6      0.4486964295094368D+00,0.4744410090307526D+00,
     7      0.5003288287668549D+00,0.5265894704770495D+00,
     8      0.5538936701851433D+00,0.5836902532830273D+00,
     9      0.6193977262766783D+00,0.6847470482229011D+00/
c
        data xs21/
     1      0.3805049123127608D+00,0.1149438977324809D+01,
     2      0.1941452163441607D+01,0.2769840780469687D+01,
     3      0.3645424894404909D+01,0.4576544425712020D+01,
     4      0.5569519975253000D+01,0.6629343532589158D+01,
     5      0.7760336040470153D+01,0.8966453888304276D+01,
     6      0.1025106141730387D+02,0.1161645983890043D+02,
     7      0.1306374774279260D+02,0.1459328623169133D+02,
     8      0.1620562012060127D+02,0.1790268323721492D+02,
     9      0.1968940433183512D+02,0.2157705424339114D+02,
     a      0.2358887741742708D+02/
        data ws21/
     1      0.2426612108826140D+00,0.2476678564418039D+00,
     2      0.2572521076221768D+00,0.2706903750069641D+00,
     3      0.2871635978132027D+00,0.3059369755087225D+00,
     4      0.3264719757800226D+00,0.3484600540376705D+00,
     5      0.3717634859163877D+00,0.3962565112034149D+00,
     6      0.4216730853486884D+00,0.4476221979446817D+00,
     7      0.4737640345653161D+00,0.4999930972088836D+00,
     8      0.5265494021324794D+00,0.5540936855476652D+00,
     9      0.5839853810207363D+00,0.6199855205680492D+00,
     a      0.6826170050983182D+00/
c
        data xs22/
     1      0.3210090532410329D+00,0.9720613831953780D+00,
     2      0.1649325716975411D+01,0.2367595034032465D+01,
     3      0.3138481748848083D+01,0.3970389544133767D+01,
     4      0.4869055054458924D+01,0.5838372885264460D+01,
     5      0.6881252571962809D+01,0.8000297153416833D+01,
     6      0.9198078141055168D+01,0.1047689776201688D+02,
     7      0.1183831612626972D+02,0.1328296712041574D+02,
     8      0.1481093717523134D+02,0.1642260276766942D+02,
     9      0.1811977760468169D+02,0.1990735111892663D+02,
     a      0.2179615283034943D+02,0.2380908740958013D+02/
        data ws22/
     1      0.2048449590231995D+00,0.2105474000202533D+00,
     2      0.2214015173492093D+00,0.2364789625747589D+00,
     3      0.2547254841336221D+00,0.2751832364747371D+00,
     4      0.2971269212615204D+00,0.3201105469023526D+00,
     5      0.3439437797444267D+00,0.3686000399719085D+00,
     6      0.3940561525527044D+00,0.4201496732712054D+00,
     7      0.4465876683771589D+00,0.4731057446247343D+00,
     8      0.4996445178907279D+00,0.5264637696281906D+00,
     9      0.5542436948250807D+00,0.5844072994503561D+00,
     a      0.6203199026743286D+00,0.6838771294187156D+00/
c
        data xs23/
     1      0.2482661356626073D+00,0.7551014609533595D+00,
     2      0.1291786354389410D+01,0.1875022746416961D+01,
     3      0.2517683818203186D+01,0.3228700250355010D+01,
     4      0.4013606024637153D+01,0.4875444371737314D+01,
     5      0.5815766415262579D+01,0.6835533604788216D+01,
     6      0.7935786730588606D+01,0.9117919099489550D+01,
     7      0.1038344242490015D+02,0.1173346069590951D+02,
     8      0.1316836305942355D+02,0.1468807856058004D+02,
     9      0.1629284867488361D+02,0.1798437400688850D+02,
     a      0.1976748798205770D+02,0.2165319070946607D+02,
     1      0.2366559381644577D+02/
        data ws23/
     1      0.1586036321676004D+00,0.1651045818853231D+00,
     2      0.1774509483157463D+00,0.1945322307549624D+00,
     3      0.2150691648777504D+00,0.2378665265017659D+00,
     4      0.2619760184688583D+00,0.2867647275548428D+00,
     5      0.3119102085194389D+00,0.3373473757736840D+00,
     6      0.3631708377531386D+00,0.3894814646472840D+00,
     7      0.4162362003966941D+00,0.4432313387449853D+00,
     8      0.4702478747985498D+00,0.4972403564100715D+00,
     9      0.5244714930076714D+00,0.5526415743074017D+00,
     a      0.5831759919109494D+00,0.6197623824854233D+00,
     1      0.6835667156796343D+00/
c
        data xs24/
     1      0.1488189281465128D+00,0.4564693950412821D+00,
     2      0.7934408535746322D+00,0.1177033916867475D+01,
     3      0.1621566481950466D+01,0.2137788265836901D+01,
     4      0.2732786888252149D+01,0.3410415201999987D+01,
     5      0.4172089750185019D+01,0.5017722320503176D+01,
     6      0.5946596861339998D+01,0.6958090176870887D+01,
     7      0.8052164726650334D+01,0.9229506616390847D+01,
     8      0.1049119800562223D+02,0.1183811119838932D+02,
     9      0.1327051891935050D+02,0.1478827862697460D+02,
     a      0.1639156296064978D+02,0.1808198338130673D+02,
     1      0.1986425336409162D+02,0.2174924249749613D+02,
     2      0.2375929298550417D+02/
        data ws24/
     1      0.9527619504585659D-01,0.1016114863811062D+00,
     2      0.1138342309561022D+00,0.1311314926832140D+00,
     3      0.1524375585677698D+00,0.1765730468394757D+00,
     4      0.2024148031746038D+00,0.2290468423244383D+00,
     5      0.2558408336600186D+00,0.2824637307340618D+00,
     6      0.3088403932246697D+00,0.3350933906358408D+00,
     7      0.3614543683890729D+00,0.3881258386999940D+00,
     8      0.4151405030239033D+00,0.4423444951853907D+00,
     9      0.4695430154773929D+00,0.4966962147578915D+00,
     a      0.5240646739455831D+00,0.5523359509551321D+00,
     1      0.5829914238910520D+00,0.6194932106659727D+00,
     2      0.6863352963923220D+00/
c
        data xs25/
     1      0.2098412654706139D-01,0.6574622603503619D-01,
     2      0.1188313268527497D+00,0.1859411783598182D+00,
     3      0.2736094898889599D+00,0.3898729406403169D+00,
     4      0.5447092338494223D+00,0.7499792641675450D+00,
     5      0.1018489577219059D+01,0.1362095957451656D+01,
     6      0.1789635306017894D+01,0.2305930647637089D+01,
     7      0.2912273246995152D+01,0.3607683902909672D+01,
     8      0.4390137234925154D+01,0.5257431955983270D+01,
     9      0.6207750922721307D+01,0.7240026249478118D+01,
     a      0.8354139944055014D+01,0.9550844422595734D+01,
     1      0.1083128128257629D+02,0.1219629936645902D+02,
     2      0.1364607358235929D+02,0.1518034431472340D+02,
     3      0.1679920825003475D+02,0.1850430734713450D+02,
     4      0.2030042499477858D+02,0.2220135797246795D+02,
     5      0.2421541456343095D+02/
        data ws25/
     1      0.1350785957384881D-01,0.1528065344293252D-01,
     2      0.1881291873126457D-01,0.2425349133091318D-01,
     3      0.3198173843471552D-01,0.4256316708436612D-01,
     4      0.5664272363635930D-01,0.7472634340554724D-01,
     5      0.9685375533767977D-01,0.1223672927934196D+00,
     6      0.1500653067724500D+00,0.1786796308804676D+00,
     7      0.2072714118237611D+00,0.2353287699108961D+00,
     8      0.2626768419878028D+00,0.2893618209568780D+00,
     9      0.3155722054868349D+00,0.3415886866543625D+00,
     a      0.3677207409703431D+00,0.3941896439938748D+00,
     1      0.4210091085216632D+00,0.4479962878800188D+00,
     2      0.4749404147791777D+00,0.5018049527931516D+00,
     3      0.5288652996729649D+00,0.5568696625253119D+00,
     4      0.5869259625360876D+00,0.6256354188521704D+00,
     5      0.6827572362531625D+00/
c
        data xs26/
     1      0.2356588487387685D-02,0.7530956319001869D-02,
     2      0.1408173250831058D-01,0.2291954256833768D-01,
     3      0.3501933064436821D-01,0.5159464577824655D-01,
     4      0.7441672334121407D-01,0.1063894127033843D+00,
     5      0.1525607870750633D+00,0.2217315322860802D+00,
     6      0.3282971930989183D+00,0.4921255029477076D+00,
     7      0.7337273621046856D+00,0.1067516116079455D+01,
     8      0.1499384716851002D+01,0.2028947272557936D+01,
     9      0.2652762316534689D+01,0.3366534497941921D+01,
     a      0.4166239576845046D+01,0.5048625793902499D+01,
     1      0.6011437090396107D+01,0.7053527836564889D+01,
     2      0.8174890885350573D+01,0.9376478307751750D+01,
     3      0.1065969094191190D+02,0.1202570448150720D+02,
     4      0.1347507950042170D+02,0.1500795779336977D+02,
     5      0.1662479991199700D+02,0.1832753177111612D+02,
     6      0.2012072452335758D+02,0.2201728341051032D+02,
     7      0.2402857850560531D+02/
        data ws26/
     1      0.1524761510001431D-02,0.1817931941045159D-02,
     2      0.2400394254367102D-02,0.3276743721745344D-02,
     3      0.4488994527245189D-02,0.6154242357309768D-02,
     4      0.8523582667362566D-02,0.1209009095116789D-01,
     5      0.1775903636396073D-01,0.2703743243926098D-01,
     6      0.4188625234008935D-01,0.6354084442317633D-01,
     7      0.9104465800325363D-01,0.1217488727604254D+00,
     8      0.1531482880140613D+00,0.1837876272022165D+00,
     9      0.2131119230363022D+00,0.2410736686719048D+00,
     a      0.2678608232545373D+00,0.2937634140823849D+00,
     1      0.3191202635910268D+00,0.3442957341822134D+00,
     2      0.3696394952715237D+00,0.3953931309062216D+00,
     3      0.4215900890677036D+00,0.4480698163324660D+00,
     4      0.4746352733013586D+00,0.5012456910638461D+00,
     5      0.5281641223802437D+00,0.5561410942475963D+00,
     6      0.5860555831843230D+00,0.6258187810607486D+00,
     7      0.7106759721181988D+00/
c
        data xs27/
     1      0.8900059761929686D-03,0.2773500244549777D-02,
     2      0.4967067478571699D-02,0.7673519456239493D-02,
     3      0.1110105601463826D-01,0.1558016463638721D-01,
     4      0.2183882488451939D-01,0.3124268041089357D-01,
     5      0.4636534534156120D-01,0.7296337067468865D-01,
     6      0.1245253518762608D+00,0.2264844073925937D+00,
     7      0.4080168046365694D+00,0.6868036161971084D+00,
     8      0.1067356310926221D+01,0.1547374613019256D+01,
     9      0.2122333052858073D+01,0.2787562233662757D+01,
     a      0.3539052318545281D+01,0.4373718550020571D+01,
     1      0.5289475309695980D+01,0.6285258915359327D+01,
     2      0.7361002887971802D+01,0.8517451036489476D+01,
     3      0.9755678242776144D+01,0.1107641884011278D+02,
     4      0.1247957599663193D+02,0.1396420437528361D+02,
     5      0.1552889527641378D+02,0.1717227363795267D+02,
     6      0.1889319174335914D+02,0.2068819444943268D+02,
     7      0.2254639690481814D+02,0.2443918016091106D+02/
        data ws27/
     1      0.5720839241402956D-03,0.6379657563383231D-03,
     2      0.7692827204867902D-03,0.9642592398328269D-03,
     3      0.1233515720217732D-02,0.1655597299092846D-02,
     4      0.2399110366070954D-02,0.3714675887098084D-02,
     5      0.6193931478120942D-02,0.1141854210824250D-01,
     6      0.2280046828316033D-01,0.4376601013481069D-01,
     7      0.7273489523100148D-01,0.1049443721775837D+00,
     8      0.1371774025145449D+00,0.1681616633644519D+00,
     9      0.1976167829712175D+00,0.2256673700991595D+00,
     a      0.2525821307207455D+00,0.2786721328129627D+00,
     1      0.3042603334429509D+00,0.3296725228037677D+00,
     2      0.3552094226423027D+00,0.3810692681507067D+00,
     3      0.4072524784118721D+00,0.4335520642666400D+00,
     4      0.4596733502210447D+00,0.4853926658442933D+00,
     5      0.5106431714748835D+00,0.5355071956160309D+00,
     6      0.5599936603162547D+00,0.5823887819451633D+00,
     7      0.5876664762810564D+00,0.4809678617110528D+00/
c
        data xs28/
     1      0.3888409704980006D-04,0.1998592302697424D-03,
     2      0.4125944077019328D-03,0.6835670853213474D-03,
     3      0.1047305172038495D-02,0.1555260782451407D-02,
     4      0.2292947290855738D-02,0.3429360730359222D-02,
     5      0.5381101541500450D-02,0.9517469841729518D-02,
     6      0.2209434833321296D-01,0.6856704182172699D-01,
     7      0.1919492991961001D+00,0.4191552924010502D+00,
     8      0.7563127516631101D+00,0.1199969734864978D+01,
     9      0.1743900388144719D+01,0.2381785665981915D+01,
     a      0.3108162753170317D+01,0.3918712024446937D+01,
     1      0.4810324219690742D+01,0.5781140961629032D+01,
     2      0.6830612418131090D+01,0.7959489655853091D+01,
     3      0.9169561743115338D+01,0.1046300484127320D+02,
     4      0.1184156774757512D+02,0.1330612086749811D+02,
     5      0.1485691443345835D+02,0.1649441093720185D+02,
     6      0.1822046014637496D+02,0.2003642531824334D+02,
     7      0.2195698270221312D+02,0.2403800869315479D+02/
        data ws28/
     1      0.3418223415519108D-04,0.6090775660788185D-04,
     2      0.7550287164281223D-04,0.9882823014287491D-04,
     3      0.1353821466496369D-03,0.1923104021642913D-03,
     4      0.2856066195448874D-03,0.4573674950145027D-03,
     5      0.8445568176817392D-03,0.2040391401529968D-02,
     6      0.7252888960534763D-02,0.2491193978681773D-01,
     7      0.5504335678478973D-01,0.8985437172939705D-01,
     8      0.1245702112356292D+00,0.1575279033007831D+00,
     9      0.1884095190203754D+00,0.2173898262939608D+00,
     a      0.2448068077808498D+00,0.2710397489188958D+00,
     1      0.2964792474380519D+00,0.3215325965117131D+00,
     2      0.3466227680512424D+00,0.3721413373392838D+00,
     3      0.3983347459319631D+00,0.4251921334161726D+00,
     4      0.4524730891981678D+00,0.4799005094217993D+00,
     5      0.5073750924452041D+00,0.5351457760449594D+00,
     6      0.5640920094886253D+00,0.5906904462620981D+00,
     7      0.6058475603069529D+00,0.5895149890702301D+00/
c
        data xs29/
     1      0.1084451861691074D-04,0.3536954241355290D-04,
     2      0.9812941400469883D-04,0.9961972730029050D-04,
     3      0.1107412525476113D-03,0.1149525051207237D-03,
     4      0.1454894881871127D-03,0.1556300027690178D-03,
     5      0.1644067627631533D-03,0.1892078694732227D-03,
     6      0.1922658164754994D-03,0.6064701518630259D-03,
     7      0.6225254328946500D-03,0.6396617368353823D-03,
     8      0.3278405287697243D-02,0.4218104301734411D-01,
     9      0.1822184836857371D+00,0.4368456783825452D+00,
     a      0.8057356400665123D+00,0.1283534850692909D+01,
     1      0.1862604878491752D+01,0.2535073841129145D+01,
     2      0.3293986006707706D+01,0.4134033666955143D+01,
     3      0.5052084830019997D+01,0.6047308980306829D+01,
     4      0.7120545934582362D+01,0.8272752303423092D+01,
     5      0.9503035525659898D+01,0.1080781473236842D+02,
     6      0.1218285065899012D+02,0.1362828854988517D+02,
     7      0.1515324157132501D+02,0.1677475045094998D+02,
     8      0.1850684282021535D+02,0.2044947588480109D+02,
     9      0.2246904670606464D+02/
        data ws29/
     1      0.7038742644120121D-05,0.9038298575717543D-05,
     2      0.2493425609898636D-01,-.3462672196934129D-01,
     3      0.4136732021699359D-01,-.3640617985357597D-01,
     4      0.3563089059361323D-01,-.7286819400773069D-01,
     5      0.5105545379070175D-01,-.3804699345288361D-01,
     6      0.2905427689729759D-01,0.1393978339982059D-01,
     7      -.2836697801039164D-01,0.1471274445456523D-01,
     8      0.2449971284299176D-02,0.2688664848536731D-01,
     9      0.6265697369914192D-01,0.9941528902832976D-01,
     a      0.1351271058217804D+00,0.1686332957105248D+00,
     1      0.1995908988274685D+00,0.2281434128919084D+00,
     2      0.2547105046359893D+00,0.2799133097042467D+00,
     3      0.3044953213250279D+00,0.3291371383036638D+00,
     4      0.3541691997668450D+00,0.3793068298805426D+00,
     5      0.4037155270061492D+00,0.4266789965679773D+00,
     6      0.4486756657862525D+00,0.4720163812648909D+00,
     7      0.4998065290735255D+00,0.5339313979232739D+00,
     8      0.5719315769947966D+00,0.5483711732111221D+00,
     9      0.4142285670166452D+00/
c
        data xs30/
     1      0.9383692507469732D-06,0.2968598682429047D-05,
     2      0.5480135747859549D-05,0.8882309138869005D-05,
     3      0.1377190426410380D-04,0.2114938310171510D-04,
     4      0.3302260083515307D-04,0.5468129426589030D-04,
     5      0.1081243600147275D-03,0.5836914253691816D-03,
     6      0.4156465006196086D-01,0.2063886825058709D+00,
     7      0.4975273614327698D+00,0.9076387882262793D+00,
     8      0.1428666272140287D+01,0.2053190024713825D+01,
     9      0.2775060001145942D+01,0.3589570609289474D+01,
     a      0.4493439144917517D+01,0.5484765260961655D+01,
     1      0.6563007150093132D+01,0.7728820117640062D+01,
     2      0.8983516260886866D+01,0.1032823194557076D+02,
     3      0.1176336106931675D+02,0.1328873647843941D+02,
     4      0.1490454729377551D+02,0.1661250101533295D+02,
     5      0.1841850242192802D+02,0.2033750449232517D+02,
     6      0.2247265798978797D+02/
        data ws30/
     1      0.6054372141490441D-06,0.7041259772038848D-06,
     2      0.9159556239011105D-06,0.1280602891199475D-05,
     3      0.1882705379808585D-05,0.2911704842219567D-05,
     4      0.4887608940832018D-05,0.9746086943017523D-05,
     5      0.3005437383892944D-04,0.7547207395913029D-03,
     6      0.3179228166356317D-01,0.7290787551645030D-01,
     7      0.1120347716332354D+00,0.1486158734569309D+00,
     8      0.1826837362447185D+00,0.2145734484763226D+00,
     9      0.2447338660378923D+00,0.2736257890852822D+00,
     a      0.3016939836830113D+00,0.3293756807850941D+00,
     1      0.3570920857038262D+00,0.3851620571341180D+00,
     2      0.4136660604207481D+00,0.4424268641818212D+00,
     3      0.4711916838536064D+00,0.4999003472034489D+00,
     4      0.5288639675108805D+00,0.5587468699335812D+00,
     5      0.5914097394887634D+00,0.6259334005138946D+00,
     6      0.7296153285432272D+00/
c
        data xs31/
     1      0.9972760395271507D-07,0.3156571573015299D-06,
     2      0.5838252903345165D-06,0.9500795066230694D-06,
     3      0.1483668903290436D-05,0.2307753945274271D-05,
     4      0.3694502094475638D-05,0.6496082125551084D-05,
     5      0.1591466704553568D-04,0.4717869765901384D-02,
     6      0.1100737859075096D+00,0.3534385166986738D+00,
     7      0.7272474070150103D+00,0.1222830448929287D+01,
     8      0.1831750574320495D+01,0.2546675716168263D+01,
     9      0.3361677514594842D+01,0.4272295459680954D+01,
     a      0.5275580941527988D+01,0.6370121742138135D+01,
     1      0.7555835303076613D+01,0.8833300680315041D+01,
     2      0.1020283119112838D+02,0.1166393608636557D+02,
     3      0.1321556480646428D+02,0.1485694459473258D+02,
     4      0.1658995949638086D+02,0.1841450585647549D+02,
     5      0.2033091288214750D+02/
        data ws31/
     1      0.6434873100996694D-07,0.7497806362574902D-07,
     2      0.9811080557084429D-07,0.1385895301030744D-06,
     3      0.2072507269276423D-06,0.3303055991508012D-06,
     4      0.5891010078840259D-06,0.1359253857172051D-05,
     5      0.6699437420040374D-05,0.1076127530880279D-01,
     6      0.5586056497266390D-01,0.9867032522201409D-01,
     7      0.1388350707618756D+00,0.1762108519560343D+00,
     8      0.2110494721411024D+00,0.2437713284796949D+00,
     9      0.2748380614815450D+00,0.3047225964118315D+00,
     a      0.3339147000059104D+00,0.3628909245422213D+00,
     1      0.3919940383671215D+00,0.4212855555810948D+00,
     2      0.4505568001263470D+00,0.4795508270457126D+00,
     3      0.5082032926909292D+00,0.5366623189982458D+00,
     4      0.5661116256818458D+00,0.5889381420072928D+00,
     5      0.6834040961827335D+00/
c
	ier=0
	if (rky.lt.0.5d-6) then
	   ier=8
	   return
	endif
c
c     Quad for rky starts at:  0.24D+02
c
        if (rky.ge.0.24D+02) then
           nlam= 3
           do i=1,nlam
              xs(i)=xs01(i)
              ws(i)=ws01(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.23D+02
c
        if ((rky.ge.0.23D+02).and.(rky.lt.0.24D+02)) then
           nlam= 4
           do i=1,nlam
              xs(i)=xs02(i)
              ws(i)=ws02(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.22D+02
c
        if ((rky.ge.0.22D+02).and.(rky.lt.0.23D+02)) then
           nlam= 6
           do i=1,nlam
              xs(i)=xs03(i)
              ws(i)=ws03(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.21D+02
c
        if ((rky.ge.0.21D+02).and.(rky.lt.0.22D+02)) then
           nlam= 7
           do i=1,nlam
              xs(i)=xs04(i)
              ws(i)=ws04(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.20D+02
c
        if ((rky.ge.0.20D+02).and.(rky.lt.0.21D+02)) then
           nlam= 8
           do i=1,nlam
              xs(i)=xs05(i)
              ws(i)=ws05(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.19D+02
c
        if ((rky.ge.0.19D+02).and.(rky.lt.0.20D+02)) then
           nlam= 8
           do i=1,nlam
              xs(i)=xs06(i)
              ws(i)=ws06(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.18D+02
c
        if ((rky.ge.0.18D+02).and.(rky.lt.0.19D+02)) then
           nlam= 9
           do i=1,nlam
              xs(i)=xs07(i)
              ws(i)=ws07(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.17D+02
c
        if ((rky.ge.0.17D+02).and.(rky.lt.0.18D+02)) then
           nlam=10
           do i=1,nlam
              xs(i)=xs08(i)
              ws(i)=ws08(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.16D+02
c
        if ((rky.ge.0.16D+02).and.(rky.lt.0.17D+02)) then
           nlam=10
           do i=1,nlam
              xs(i)=xs09(i)
              ws(i)=ws09(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.15D+02
c
        if ((rky.ge.0.15D+02).and.(rky.lt.0.16D+02)) then
           nlam=11
           do i=1,nlam
              xs(i)=xs10(i)
              ws(i)=ws10(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.14D+02
c
        if ((rky.ge.0.14D+02).and.(rky.lt.0.15D+02)) then
           nlam=11
           do i=1,nlam
              xs(i)=xs11(i)
              ws(i)=ws11(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.13D+02
c
        if ((rky.ge.0.13D+02).and.(rky.lt.0.14D+02)) then
           nlam=12
           do i=1,nlam
              xs(i)=xs12(i)
              ws(i)=ws12(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.12D+02
c
        if ((rky.ge.0.12D+02).and.(rky.lt.0.13D+02)) then
           nlam=13
           do i=1,nlam
              xs(i)=xs13(i)
              ws(i)=ws13(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.11D+02
c
        if ((rky.ge.0.11D+02).and.(rky.lt.0.12D+02)) then
           nlam=13
           do i=1,nlam
              xs(i)=xs14(i)
              ws(i)=ws14(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.10D+02
c
        if ((rky.ge.0.10D+02).and.(rky.lt.0.11D+02)) then
           nlam=14
           do i=1,nlam
              xs(i)=xs15(i)
              ws(i)=ws15(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.90D+01
c
        if ((rky.ge.0.90D+01).and.(rky.lt.0.10D+02)) then
           nlam=15
           do i=1,nlam
              xs(i)=xs16(i)
              ws(i)=ws16(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.80D+01
c
        if ((rky.ge.0.80D+01).and.(rky.lt.0.90D+01)) then
           nlam=16
           do i=1,nlam
              xs(i)=xs17(i)
              ws(i)=ws17(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.70D+01
c
        if ((rky.ge.0.70D+01).and.(rky.lt.0.80D+01)) then
           nlam=17
           do i=1,nlam
              xs(i)=xs18(i)
              ws(i)=ws18(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.60D+01
c
        if ((rky.ge.0.60D+01).and.(rky.lt.0.70D+01)) then
           nlam=18
           do i=1,nlam
              xs(i)=xs19(i)
              ws(i)=ws19(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.50D+01
c
        if ((rky.ge.0.50D+01).and.(rky.lt.0.60D+01)) then
           nlam=18
           do i=1,nlam
              xs(i)=xs20(i)
              ws(i)=ws20(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.40D+01
c
        if ((rky.ge.0.40D+01).and.(rky.lt.0.50D+01)) then
           nlam=19
           do i=1,nlam
              xs(i)=xs21(i)
              ws(i)=ws21(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.30D+01
c
        if ((rky.ge.0.30D+01).and.(rky.lt.0.40D+01)) then
           nlam=20
           do i=1,nlam
              xs(i)=xs22(i)
              ws(i)=ws22(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.20D+01
c
        if ((rky.ge.0.20D+01).and.(rky.lt.0.30D+01)) then
           nlam=21
           do i=1,nlam
              xs(i)=xs23(i)
              ws(i)=ws23(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.10D+01
c
        if ((rky.ge.0.10D+01).and.(rky.lt.0.20D+01)) then
           nlam=23
           do i=1,nlam
              xs(i)=xs24(i)
              ws(i)=ws24(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.10D+00
c
        if ((rky.ge.0.20D+00).and.(rky.lt.0.10D+01)) then
           nlam=29
           do i=1,nlam
              xs(i)=xs25(i)
              ws(i)=ws25(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.50D-01
c
        if ((rky.ge.0.50D-01).and.(rky.lt.0.20D+00)) then
           nlam=33
           do i=1,nlam
              xs(i)=xs26(i)
              ws(i)=ws26(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.50D-02
c
        if ((rky.ge.0.50D-02).and.(rky.lt.0.50D-01)) then
           nlam=34
           do i=1,nlam
              xs(i)=xs27(i)
              ws(i)=ws27(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.50D-03
c
        if ((rky.ge.0.50D-03).and.(rky.lt.0.50D-02)) then
           nlam=34
           do i=1,nlam
              xs(i)=xs28(i)
              ws(i)=ws28(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.50D-04
c
        if ((rky.ge.0.50D-04).and.(rky.lt.0.50D-03)) then
           nlam=37
           do i=1,nlam
              xs(i)=xs29(i)
              ws(i)=ws29(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.50D-05 (11 digits only)
c
        if ((rky.ge.0.50D-05).and.(rky.lt.0.50D-04)) then
           nlam=31
           do i=1,nlam
              xs(i)=xs30(i)
              ws(i)=ws30(i)
           enddo
c
           return
        endif
c
c     Quad for rky starts at:  0.50D-06 (10 digits only)
c
        if ((rky.ge.0.50D-06).and.(rky.lt.0.50D-05)) then
           nlam=29
           do i=1,nlam
              xs(i)=xs31(i)
              ws(i)=ws31(i)
           enddo
c
           return
        endif
c
	return
	end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c        this is the end of the debugging code and the beginning of the
c        hankel function code proper.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine hank101(z,h0,h1)
        implicit real *8 (a-h,o-z)
        doublecomplex  z,h0,h1,h0u,h0r,h1u,h1r,
     1      fj0,fj1,y0,y1,com,zu,zr,ima,ser2,ser3,z2,
     2      cclog
        real *8 rea(2)
        equivalence (rea(1),com)
        data ima/(0.0d0,1.0d0)/,pi/0.31415926535897932D+01/,
     1      done/1.0d0/,two/2.0d0/
        data gamma/0.5772156649015328606d+00/
c
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for an arbitrary user-specified complex number z. its 
c        principal claim to fame is that it is valid on the whole
c        complex plane, and is reasonably accurate (14-digit 
c        relative accuracy) and reasonably fast.
c
c                      input parameters:
c
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c
c                      output parameters:
c
c  h0, h1 - the said Hankel functions
c        
c       
c        . . . if z in the upper half-plane - act accordingly
c
        com=z 
        if(rea(2) .lt. 0) goto 1400
        call hank101u(z,ier,h0,h1)
        return
 1400 continue
c
c       if z is in the right lower quadrant - act accordingly
c
        if(rea(1) .lt. 0) goto 2000
        call hank101r(z,ier,h0,h1)
        return
 2000 continue
c
c       z is in the left lower quadrant. compute 
c       h0, h1 at the points zu, zr obtained from z by reflection
c       in the x and y axis, respectively
c
        zu=dconjg(z)
        zr=-zu
c
        call hank101u(zu,ier,h0u,h1u)
        call hank101r(zr,ier,h0r,h1r)
cccc         call prin2('in hank101, h0u=*',h0u,2)
cccc         call prin2('in hank101, h1u=*',h1u,2)
cccc         call prin2('in hank101, h0r=*',h0r,2)
cccc         call prin2('in hank101, h1r=*',h1r,2)
c
c       compute the functions j0, j1, y0, y1
c       at the point zr
c
        half=1
        half=half/2
        y0=(h0u+h0r)*half/ima
        fj0=-(h0u-h0r)*half
c
        y1=-(h1u-h1r)*half/ima
        fj1=(h1u+h1r)*half
cccc        call prin2('fj0 as computed*',fj0,2)
cccc        call prin2('y0 as compoted*',y0,2)        
cccc        call prin2('fj1 as computed*',fj1,2)
cccc        call prin2('y1 as compoted*',y1,2)        
c
c        finally, compute h0, h1
c
cccc        two=2
cccc        done=1
cccc        pi=datan(done)*4
c
c       calculate ser2, ser3
c
         z2=-dconjg(z)
         cclog=cdlog(z2)
cccc         ser2=y0-fj0*2/pi*cdlog(z2)
cccc         ser3=y1-fj1*2/pi*cdlog(z2)
         ser2=y0-fj0*2/pi*cclog
         ser3=y1-fj1*2/pi*cclog
c
c       reflect all of these in the imaginary axis
c
        fj0=dconjg(fj0)
        fj1=-dconjg(fj1)
c
        ser2=dconjg(ser2)
        ser3=-dconjg(ser3)
c
c       reconstitute y0, y1
c
cccc        y0=ser2+fj0*2/pi*cdlog(z)
cccc        y1=ser3+fj1*2/pi*cdlog(z)
        cclog=cdlog(z)
        y0=ser2+fj0*2/pi*cclog
        y1=ser3+fj1*2/pi*cclog
c
cccc        call prin2('evaluating h0, h1, fj0=*',fj0,2)
cccc        call prin2('evaluating h0, h1, fj1=*',fj1,2)
cccc        call prin2('evaluating h0, h1, y0=*',y0,2)
cccc        call prin2('evaluating h0, h1, y1=*',y1,2)
c
        h0=fj0+ima*y0
        h1=fj1+ima*y1
        return
        end
c
c
c
c
c
        subroutine hank101u(z,ier,h0,h1)
        implicit real *8 (a-h,o-z)
        doublecomplex  z,com,ima,cd,h0,h1,ccex
        dimension rea(2)
c
        real *8 c0p1(34),c0p1b(36),buf01(2)
        equivalence (c0p1(34),buf01(1)),
     1      (c0p1b(1),buf01(2)),(rea(1),com)
c
        real *8 c1p1(34),c1p1b(36),buf11(2)
        equivalence (c1p1(34),buf11(1)),
     1      (c1p1b(1),buf11(2))
c
cccc        real *8 c0p2(34),c0p2b(36),buf02(2)
        real *8 c0p2(34),c0p2b(28),buf02(2)
        equivalence (c0p2(34),buf02(1)),
     1      (c0p2b(1),buf02(2))
c
        real *8 c1p2(34),c1p2b(28),buf12(2)
        equivalence (c1p2(34),buf12(1)),
     1      (c1p2b(1),buf12(2))
c
        data ima/(0.0d0,1.0d0)/
c
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for a user-specified complex number z in the upper half-plane.
c        it is reasonably accurate (14-digit relative accuracy) 
c        and reasonably fast.
c        
c
c                      input parameters:
c
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c
c                      output parameters:
c
c  ier - error return code. 
c         ier=0 means successful conclusion
c         ier=4 means that z is not in the upper half-plane
c  h0, h1 - the said Hankel functions
c        
        data c0p1/
     1     -.6619836118357782D-12,  -.6619836118612709D-12,
     2     -.7307514264754200D-21,  0.3928160926261892D-10,
     3     0.5712712520172854D-09,  -.5712712519967086D-09,
     4     -.1083820384008718D-07,  -.1894529309455499D-18,
     5     0.7528123700585197D-07,  0.7528123700841491D-07,
     6     0.1356544045548053D-16,  -.8147940452202855D-06,
     7     -.3568198575016769D-05,  0.3568198574899888D-05,
     8     0.2592083111345422D-04,  0.4209074870019400D-15,
     9     -.7935843289157352D-04,  -.7935843289415642D-04,
     a     -.6848330800445365D-14,  0.4136028298630129D-03,
     1     0.9210433149997867D-03,  -.9210433149680665D-03,
     2     -.3495306809056563D-02,  -.6469844672213905D-13,
     3     0.5573890502766937D-02,  0.5573890503000873D-02,
     4     0.3767341857978150D-12,  -.1439178509436339D-01,
     5     -.1342403524448708D-01,  0.1342403524340215D-01,
     6     0.8733016209933828D-02,  0.1400653553627576D-11,
     7     0.2987361261932706D-01,  0.2987361261607835D-01/
        data c0p1b/
     8     -.3388096836339433D-11,  -.1690673895793793D+00,
     9     0.2838366762606121D+00,  -.2838366762542546D+00,
     a     0.7045107746587499D+00,  -.5363893133864181D-11,
     1     -.7788044738211666D+00,  -.7788044738130360D+00,
     2     0.5524779104964783D-11,  0.1146003459721775D+01,
     3     0.6930697486173089D+00,  -.6930697486240221D+00,
     4     -.7218270272305891D+00,  0.3633022466839301D-11,
     5     0.3280924142354455D+00,  0.3280924142319602D+00,
     6     -.1472323059106612D-11,  -.2608421334424268D+00,
     7     -.9031397649230536D-01,  0.9031397649339185D-01,
     8     0.5401342784296321D-01,  -.3464095071668884D-12,
     9     -.1377057052946721D-01,  -.1377057052927901D-01,
     a     0.4273263742980154D-13,  0.5877224130705015D-02,
     1     0.1022508471962664D-02,  -.1022508471978459D-02,
     2     -.2789107903871137D-03,  0.2283984571396129D-14,
     3     0.2799719727019427D-04,  0.2799719726970900D-04,
     4     -.3371218242141487D-16,  -.3682310515545645D-05,
     5     -.1191412910090512D-06,  0.1191412910113518D-06/
c
        data c1p1/
     1     0.4428361927253983D-12,  -.4428361927153559D-12,
     2     -.2575693161635231D-10,  -.2878656317479645D-21,
     3     0.3658696304107867D-09,  0.3658696304188925D-09,
     4     0.7463138750413651D-19,  -.6748894854135266D-08,
     5     -.4530098210372099D-07,  0.4530098210271137D-07,
     6     0.4698787882823243D-06,  0.5343848349451927D-17,
     7     -.1948662942158171D-05,  -.1948662942204214D-05,
     8     -.1658085463182409D-15,  0.1316906100496570D-04,
     9     0.3645368564036497D-04,  -.3645368563934748D-04,
     a     -.1633458547818390D-03,  -.2697770638600506D-14,
     1     0.2816784976551660D-03,  0.2816784976676616D-03,
     2     0.2548673351180060D-13,  -.6106478245116582D-03,
     3     0.2054057459296899D-03,  -.2054057460218446D-03,
     4     -.6254962367291260D-02,  0.1484073406594994D-12,
     5     0.1952900562500057D-01,  0.1952900562457318D-01,
     6     -.5517611343746895D-12,  -.8528074392467523D-01,
     7     -.1495138141086974D+00,  0.1495138141099772D+00/
c
        data c1p1b/
     8     0.4394907314508377D+00,  -.1334677126491326D-11,
     9     -.1113740586940341D+01,  -.1113740586937837D+01,
     a     0.2113005088866033D-11,  0.1170212831401968D+01,
     1     0.1262152242318805D+01,  -.1262152242322008D+01,
     2     -.1557810619605511D+01,  0.2176383208521897D-11,
     3     0.8560741701626648D+00,  0.8560741701600203D+00,
     4     -.1431161194996653D-11,  -.8386735092525187D+00,
     5     -.3651819176599290D+00,  0.3651819176613019D+00,
     6     0.2811692367666517D+00,  -.5799941348040361D-12,
     7     -.9494630182937280D-01,  -.9494630182894480D-01,
     8     0.1364615527772751D-12,  0.5564896498129176D-01,
     9     0.1395239688792536D-01,  -.1395239688799950D-01,
     a     -.5871314703753967D-02,  0.1683372473682212D-13,
     1     0.1009157100083457D-02,  0.1009157100077235D-02,
     2     -.8997331160162008D-15,  -.2723724213360371D-03,
     3     -.2708696587599713D-04,  0.2708696587618830D-04,
     4     0.3533092798326666D-05,  -.1328028586935163D-16,
     5     -.1134616446885126D-06,  -.1134616446876064D-06/
c
        data c0p2/
     1     0.5641895835516786D+00,  -.5641895835516010D+00,
     2     -.3902447089770041D-09,  -.3334441074447365D-11,
     3     -.7052368835911731D-01,  -.7052368821797083D-01,
     4     0.1957299315085370D-08,  -.3126801711815631D-06,
     5     -.3967331737107949D-01,  0.3967327747706934D-01,
     6     0.6902866639752817D-04,  0.3178420816292497D-06,
     7     0.4080457166061280D-01,  0.4080045784614144D-01,
     8     -.2218731025620065D-04,  0.6518438331871517D-02,
     9     0.9798339748600499D-01,  -.9778028374972253D-01,
     a     -.3151825524811773D+00,  -.7995603166188139D-03,
     1     0.1111323666639636D+01,  0.1116791178994330D+01,
     2     0.1635711249533488D-01,  -.8527067497983841D+01,
     3     -.2595553689471247D+02,  0.2586942834408207D+02,
     4     0.1345583522428299D+03,  0.2002017907999571D+00,
     5     -.3086364384881525D+03,  -.3094609382885628D+03,
     6     -.1505974589617013D+01,  0.1250150715797207D+04,
     7     0.2205210257679573D+04,  -.2200328091885836D+04/
        data c0p2b/
     8     -.6724941072552172D+04,  -.7018887749450317D+01,
     9     0.8873498980910335D+04,  0.8891369384353965D+04,
     a     0.2008805099643591D+02,  -.2030681426035686D+05,
     1     -.2010017782384992D+05,  0.2006046282661137D+05,
     2     0.3427941581102808D+05,  0.3432892927181724D+02,
     3     -.2511417407338804D+05,  -.2516567363193558D+05,
     4     -.3318253740485142D+02,  0.3143940826027085D+05,
     5     0.1658466564673543D+05,  -.1654843151976437D+05,
     6     -.1446345041326510D+05,  -.1645433213663233D+02,
     7     0.5094709396573681D+04,  0.5106816671258367D+04,
     8     0.3470692471612145D+01,  -.2797902324245621D+04,
     9     -.5615581955514127D+03,  0.5601021281020627D+03,
     a     0.1463856702925587D+03,  0.1990076422327786D+00,
     1     -.9334741618922085D+01,  -.9361368967669095D+01/
c
        data c1p2/
     1     -.5641895835446003D+00,  -.5641895835437973D+00,
     2     0.3473016376419171D-10,  -.3710264617214559D-09,
     3     0.2115710836381847D+00,  -.2115710851180242D+00,
     4     0.3132928887334847D-06,  0.2064187785625558D-07,
     5     -.6611954881267806D-01,  -.6611997176900310D-01,
     6     -.3386004893181560D-05,  0.7146557892862998D-04,
     7     -.5728505088320786D-01,  0.5732906930408979D-01,
     8     -.6884187195973806D-02,  -.2383737409286457D-03,
     9     0.1170452203794729D+00,  0.1192356405185651D+00,
     a     0.8652871239920498D-02,  -.3366165876561572D+00,
     1     -.1203989383538728D+01,  0.1144625888281483D+01,
     2     0.9153684260534125D+01,  0.1781426600949249D+00,
     3     -.2740411284066946D+02,  -.2834461441294877D+02,
     4     -.2192611071606340D+01,  0.1445470231392735D+03,
     5     0.3361116314072906D+03,  -.3270584743216529D+03,
     6     -.1339254798224146D+04,  -.1657618537130453D+02,
     7     0.2327097844591252D+04,  0.2380960024514808D+04/
        data c1p2b/
     8     0.7760611776965994D+02,  -.7162513471480693D+04,
     9     -.9520608696419367D+04,  0.9322604506839242D+04,
     a     0.2144033447577134D+05,  0.2230232555182369D+03,
     1     -.2087584364240919D+05,  -.2131762020653283D+05,
     2     -.3825699231499171D+03,  0.3582976792594737D+05,
     3     0.2642632405857713D+05,  -.2585137938787267D+05,
     4     -.3251446505037506D+05,  -.3710875194432116D+03,
     5     0.1683805377643986D+05,  0.1724393921722052D+05,
     6     0.1846128226280221D+03,  -.1479735877145448D+05,
     7     -.5258288893282565D+04,  0.5122237462705988D+04,
     8     0.2831540486197358D+04,  0.3905972651440027D+02,
     9     -.5562781548969544D+03,  -.5726891190727206D+03,
     a     -.2246192560136119D+01,  0.1465347141877978D+03,
     1     0.9456733342595993D+01,  -.9155767836700837D+01/
c
c        if the user-specified z is in the lower half-plane
c        - bomb out
c
        ier=0
        com=z
        if(rea(2) .ge. 0) goto 1200
        ier=4
        return
 1200 continue
c
        done=1
        thresh1=1**2
        thresh2=3.7**2
        thresh3=20**2
c
c       check if if the user-specified z is in one of the 
c       intermediate regimes 
c
        d=z*dconjg(z)
        if( (d .lt. thresh1) .or. (d .gt. thresh3) ) goto 3000
c
c        the user-specified z is in one of the intermediate regimes.
c        act accordingly
c
c
        if(d .gt. thresh2) goto 2000
c
c       z is in the first intermediate regime: its absolute value is 
c       between 1 and 3.7. act accordingly
c
c       . . . evaluate the expansion
c
        cd=done/cdsqrt(z)
        ccex=cdexp(ima*z)*cd
        m=35
        call hank101p(c0p1,m,cd,h0)
cccc         call prin2('after hank101p, h0=*',h0,2)
cccc        h0=h0*cdexp(ima*z)*cd
        h0=h0*ccex
        h0=h0*z**9
c
        call hank101p(c1p1,m,cd,h1)
cccc         call prin2('after hank101p, h1=*',h1,2)
cccc        h1=h1*cdexp(ima*z)*cd
        h1=h1*ccex
        h1=h1*z**9
        return
 2000 continue
c
c       z is in the second intermediate regime: its absolute value is
c       between 3.7 and 20. act accordingly.
c
        cd=done/cdsqrt(z)
        ccex=cdexp(ima*z)*cd
        m=31
        call hank101p(c0p2,m,cd,h0)
cccc         call prin2('after hank101p, h0=*',h0,2)
cccc        h0=h0*cdexp(ima*z)*cd
        h0=h0*ccex
c
ccccc        cd=done/cdsqrt(z)
        m=31
        call hank101p(c1p2,m,cd,h1)
cccc         call prin2('after hank101p, h1=*',h1,2)
cccc        h1=h1*cdexp(ima*z)*cd
        h1=h1*ccex
        return
 3000 continue
c
c        z is either in the local regime or the asymptotic one.
c        if it is in the local regime - act accordingly.
c
        if(d .gt. 50.d0) goto 4000
        call hank101l(z,h0,h1)
        return
c
c        z is in the asymptotic regime. act accordingly.
c
 4000 continue
        call hank101a(z,h0,h1)
        return
        end
c
c
c
c
        subroutine hank101p(p,m,z,f)
        implicit real *8 (a-h,o-z)
        doublecomplex  p(1),z,f
c
c       evaluate a polynomial at a point
c
        f=p(m)
        do 1200 i=m-1,1,-1
        f=f*z+p(i)
 1200 continue
        return
        end




c
c
c
c
c
        subroutine hank101a(z,h0,h1)
        implicit real *8 (a-h,o-z)
        dimension p(18),q(18),p1(18),q1(18),rea(2)
        doublecomplex  z,zinv,pp,cd,qq,ima,h0,h1,pp1,qq1,
     1      com,z2,y0,y1,fj0,fj1,ser2,ser3
        equivalence (rea(1),com)
        data ima/(0.0d0,1.0d0)/,pi/0.31415926535897932D+01/,
     1      done/1.0d0/,two/2.0d0/
c
         data p/
     1     0.1000000000000000D+01,  -.7031250000000000D-01,
     2     0.1121520996093750D+00,  -.5725014209747314D+00,
     3     0.6074042001273483D+01,  -.1100171402692467D+03,
     4     0.3038090510922384D+04,  -.1188384262567833D+06,
     5     0.6252951493434797D+07,  -.4259392165047669D+09,
     6     0.3646840080706556D+11,  -.3833534661393944D+13,
     7     0.4854014686852901D+15,  -.7286857349377657D+17,
     8     0.1279721941975975D+20,  -.2599382102726235D+22,
     9     0.6046711487532401D+24,  -.1597065525294211D+27/
c
         data q/
     1     -.1250000000000000D+00,  0.7324218750000000D-01,
     2     -.2271080017089844D+00,  0.1727727502584457D+01,
     3     -.2438052969955606D+02,  0.5513358961220206D+03,
     4     -.1825775547429317D+05,  0.8328593040162893D+06,
     5     -.5006958953198893D+08,  0.3836255180230434D+10,
     6     -.3649010818849834D+12,  0.4218971570284096D+14,
     7     -.5827244631566907D+16,  0.9476288099260110D+18,
     8     -.1792162323051699D+21,  0.3900121292034000D+23,
     9     -.9677028801069847D+25,  0.2715581773544907D+28/

         data p1/
     1     0.1000000000000000D+01,  0.1171875000000000D+00,
     2     -.1441955566406250D+00,  0.6765925884246826D+00,
     3     -.6883914268109947D+01,  0.1215978918765359D+03,
     4     -.3302272294480852D+04,  0.1276412726461746D+06,
     5     -.6656367718817687D+07,  0.4502786003050393D+09,
     6     -.3833857520742789D+11,  0.4011838599133198D+13,
     7     -.5060568503314726D+15,  0.7572616461117957D+17,
     8     -.1326257285320556D+20,  0.2687496750276277D+22,
     9     -.6238670582374700D+24,  0.1644739123064188D+27/
c
         data q1/
     1     0.3750000000000000D+00,  -.1025390625000000D+00,
     2     0.2775764465332031D+00,  -.1993531733751297D+01,
     3     0.2724882731126854D+02,  -.6038440767050702D+03,
     4     0.1971837591223663D+05,  -.8902978767070679D+06,
     5     0.5310411010968522D+08,  -.4043620325107754D+10,
     6     0.3827011346598606D+12,  -.4406481417852279D+14,
     7     0.6065091351222699D+16,  -.9833883876590680D+18,
     8     0.1855045211579829D+21,  -.4027994121281017D+23,
     9     0.9974783533410457D+25,  -.2794294288720121D+28/
c
c        evaluate the asymptotic expansion for h0,h1 at
c        the user-supplied point z, provided it is not 
c        in the fourth quadrant
c
        com=z
cccc          if(2 .ne. 3) goto 2200
        if( (rea(1) .lt. 0) .and. (rea(2) .lt. 0) ) goto 2200

        m=10
cccc        done=1
        zero=0
        zinv=done/z
        pp=zero
        pp1=zero
        cd=done
cccc         pi=datan(done)*4
        do 1800 i=1,m
        pp=pp+p(i)*cd
c
        pp1=pp1+p1(i)*cd
c
        cd=cd*zinv**2
cccc         call prinf('i=*',i,1)
cccc         call prin2('and cd=*',cd,2)
cccc         call prin2('and cd times p(i)=*',cd*p(i),2)
 1800 continue
cccc         call prin2('and pp=*',pp,2)
c
        cd=zinv
        qq=zero
        qq1=zero
        do 2000 i=1,m
        qq=qq+q(i)*cd
c
        qq1=qq1+q1(i)*cd
c
        cd=cd*zinv**2
cccc         call prinf('i=*',i,1)
cccc         call prin2('and cd=*',cd,2)
 2000 continue
cccc         call prin2('and qq=*',qq,2)
c
        h0=pp+ima*qq
ccc        cdd=
        h0=cdsqrt(2/pi/z)*cdexp(ima*(z-pi/4)) * h0

        h1=pp1+ima*qq1

        h1=-cdsqrt(2/pi/z)*cdexp(ima*(z-pi/4)) * h1*ima
c
        return
c
 2200 continue
c
c       the point is in the third quadrant, and the asymptotic
c       expansions for h0, h1 converge slowly. evaluate functions
c       j0, j1, y0, y1 at the point z2=-dconjg(z)
c
        z2=-dconjg(z)
cccc        z2=z
c
        m=10
cccc        done=1
        zero=0
        zinv=done/z2
        pp=zero
        pp1=zero
        cd=done
cccc         pi=datan(done)*4
        do 2800 i=1,m
        pp=pp+p(i)*cd
c
        pp1=pp1+p1(i)*cd
c
        cd=cd*zinv**2
cccc         call prinf('i=*',i,1)
cccc         call prin2('and cd=*',cd,2)
cccc         call prin2('and cd times p(i)=*',cd*p(i),2)
 2800 continue
         call prin2('and pp=*',pp,2)
c
        cd=zinv
        qq=zero
        qq1=zero
        do 3000 i=1,m
        qq=qq+q(i)*cd
c
        qq1=qq1+q1(i)*cd
c
        cd=cd*zinv**2
cccc         call prinf('i=*',i,1)
cccc         call prin2('and cd=*',cd,2)
 3000 continue
         call prin2('and qq=*',qq,2)
c
        fj0=( pp*cdcos(z2-pi/4)-qq*cdsin(z2-pi/4) )*cdsqrt(2/pi/z2)
        y0= (pp*cdsin(z2-pi/4)+qq*cdcos(z2-pi/4) ) *cdsqrt(2/pi/z2)
c
        fj1=-( -pp1*cdsin(z2-pi/4)-qq1*cdcos(z2-pi/4) )*cdsqrt(2/pi/z2)
        y1= -(pp1*cdcos(z2-pi/4)-qq1*cdsin(z2-pi/4))*cdsqrt(2/pi/z2)

cccc        call prin2('fj0=*',fj0,2)
cccc        call prin2('y0=*',y0,2)
cccc        call prin2('fj1=*',fj1,2)
cccc        call prin2('y1=*',y1,2)
c
c       calculate ser2, ser3
c
         ser2=y0-fj0*2/pi*cdlog(z2)
         ser3=y1-fj1*2/pi*cdlog(z2)
c
c       reflect all of these in the imaginary axis
c
       fj0=dconjg(fj0)
        fj1=-dconjg(fj1)
c
        ser2=dconjg(ser2)
        ser3=-dconjg(ser3)
c
c       reconstitute y0, y1
c
        y0=ser2+fj0*2/pi*cdlog(z)
        y1=ser3+fj1*2/pi*cdlog(z)
c
                
        h0=fj0+ima*y0
        h1=fj1+ima*y1

 
        return
        end
c
c
c
c
c
        subroutine hank101l(z,h0,h1)
        implicit real *8 (a-h,o-z)
        dimension cj0(16),cj1(16),ser2(16),ser2der(16)
        doublecomplex  z,fj0,fj1,y0,y1,h0,h1,z2,cd,ima
c
        data gamma/0.5772156649015328606d+00/
        data ima/(0.0d0,1.0d0)/,pi/0.31415926535897932D+01/,
     1      done/1.0d0/,two/2.0d0/
c
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for a user-specified complex number z in the local regime,
c        i. e. for cdabs(z) < 1 in the upper half-plane, 
c        and for cdabs(z) < 4 in the lower half-plane, 
c        it is reasonably accurate (14-digit relative accuracy) and 
c        reasonably fast.
c
c                      input parameters:
c
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c
c                      output parameters:
c
c  h0, h1 - the said Hankel functions
c        
        data cj0/            
     1     0.1000000000000000D+01,  -.2500000000000000D+00,
     2     0.1562500000000000D-01,  -.4340277777777778D-03,
     3     0.6781684027777778D-05,  -.6781684027777778D-07,
     4     0.4709502797067901D-09,  -.2402807549524439D-11,
     5     0.9385966990329841D-14,  -.2896903392077112D-16,
     6     0.7242258480192779D-19,  -.1496334396734045D-21,
     7     0.2597802772107717D-24,  -.3842903509035085D-27,
     8     0.4901662639075363D-30,  -.5446291821194848D-33/
        data cj1/
     1     -.5000000000000000D+00,  0.6250000000000000D-01,
     2     -.2604166666666667D-02,  0.5425347222222222D-04,
     3     -.6781684027777778D-06,  0.5651403356481481D-08,
     4     -.3363930569334215D-10,  0.1501754718452775D-12,
     5     -.5214426105738801D-15,  0.1448451696038556D-17,
     6     -.3291935672814899D-20,  0.6234726653058522D-23,
     7     -.9991549123491221D-26,  0.1372465538941102D-28,
     8     -.1633887546358454D-31,  0.1701966194123390D-34/
        data ser2/
     1     0.2500000000000000D+00,  -.2343750000000000D-01,
     2     0.7957175925925926D-03,  -.1412850839120370D-04,
     3     0.1548484519675926D-06,  -.1153828185281636D-08,
     4     0.6230136717695511D-11,  -.2550971742728932D-13,
     5     0.8195247730999099D-16,  -.2121234517551702D-18,
     6     0.4518746345057852D-21,  -.8061529302289970D-24,
     7     0.1222094716680443D-26,  -.1593806157473552D-29,
     8     0.1807204342667468D-32,  -.1798089518115172D-35/
        data ser2der/
     1     0.5000000000000000D+00,  -.9375000000000000D-01,
     2     0.4774305555555556D-02,  -.1130280671296296D-03,
     3     0.1548484519675926D-05,  -.1384593822337963D-07,
     4     0.8722191404773715D-10,  -.4081554788366291D-12,
     5     0.1475144591579838D-14,  -.4242469035103405D-17,
     6     0.9941241959127275D-20,  -.1934767032549593D-22,
     7     0.3177446263369152D-25,  -.4462657240925946D-28,
     8     0.5421613028002404D-31,  -.5753886457968550D-34/
c
c        evaluate j0, j1
c
cccc        call prin2('in hank101l, z=*',z,2)
        m=16
        fj0=0
        fj1=0
        y0=0
        y1=0
        z2=z**2
        cd=1
c
cccc        done=1
cccc        two=2
cccc        pi=datan(done)*4
c        
        do 1800 i=1,m
        fj0=fj0+cj0(i)*cd
        fj1=fj1+cj1(i)*cd
        y1=y1+ser2der(i)*cd
        cd=cd*z2
        y0=y0+ser2(i)*cd
 1800 continue
        fj1=-fj1*z
cccc        call prin2('fj0=*',fj0,2)
cccc        call prin2('fj1=*',fj1,2)
c
        y0=(cdlog(z/two)+gamma)*fj0+y0
        y0=two/pi*y0
cccc        call prin2('and y0=*',y0,2)
c
        y1=y1*z
ccc        call prin2('and cd=*',cd,2)
c
        y1=-(cdlog(z/two)+gamma)*fj1+fj0/z+y1
        y1=-y1*two/pi
cccc        call prin2('and y1=*',y1,2)
c
        h0=fj0+ima*y0
        h1=fj1+ima*y1
        return
        end
c
c
c
c
c
        subroutine hank101r(z,ier,h0,h1)
        implicit real *8 (a-h,o-z)
        doublecomplex  z,com,ima,cd,h0,h1
        dimension rea(2)
        real *8 c0p1(34),c0p1b(36),buf01(2)
        equivalence (c0p1(34),buf01(1)),
     1      (c0p1b(1),buf01(2)),(rea(1),com)
c
        real *8 c1p1(34),c1p1b(36),buf11(2)
        equivalence (c1p1(34),buf11(1)),
     1      (c1p1b(1),buf11(2))
c
        real *8 c0p2(34),c0p2b(20),buf02(2)
        equivalence (c0p2(34),buf02(1)),
     1      (c0p2b(1),buf02(2))
c
        real *8 c1p2(34),c1p2b(28),buf12(2)
        equivalence (c1p2(34),buf12(1)),
     1      (c1p2b(1),buf12(2))
c
        data ima/(0.0d0,1.0d0)/
c
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for a user-specified complex number z in the right lower 
c        quadrant. it is reasonably accurate (14-digit relative 
c        accuracy) and reasonably fast.
c        
c
c                      input parameters:
c
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c
c                      output parameters:
c
c  ier - error return code. 
c         ier=0 means successful conclusion
c         ier=4 means that z is not in the right lower quadrant
c  h0, h1 - the said Hankel functions
c        
        data c0p1/
     1     -.4268441995428495D-23,  0.4374027848105921D-23,
     2     0.9876152216238049D-23,  -.1065264808278614D-20,
     3     0.6240598085551175D-19,  0.6658529985490110D-19,
     4     -.5107210870050163D-17,  -.2931746613593983D-18,
     5     0.1611018217758854D-15,  -.1359809022054077D-15,
     6     -.7718746693707326D-15,  0.6759496139812828D-14,
     7     -.1067620915195442D-12,  -.1434699000145826D-12,
     8     0.3868453040754264D-11,  0.7061853392585180D-12,
     9     -.6220133527871203D-10,  0.3957226744337817D-10,
     a     0.3080863675628417D-09,  -.1154618431281900D-08,
     1     0.7793319486868695D-08,  0.1502570745460228D-07,
     2     -.1978090852638430D-06,  -.7396691873499030D-07,
     3     0.2175857247417038D-05,  -.8473534855334919D-06,
     4     -.1053381327609720D-04,  0.2042555121261223D-04,
     5     -.4812568848956982D-04,  -.1961519090873697D-03,
     6     0.1291714391689374D-02,  0.9234422384950050D-03,
     7     -.1113890671502769D-01,  0.9053687375483149D-03/
        data c0p1b/
     8     0.5030666896877862D-01,  -.4923119348218356D-01,
     9     0.5202355973926321D+00,  -.1705244841954454D+00,
     a     -.1134990486611273D+01,  -.1747542851820576D+01,
     1     0.8308174484970718D+01,  0.2952358687641577D+01,
     2     -.3286074510100263D+02,  0.1126542966971545D+02,
     3     0.6576015458463394D+02,  -.1006116996293757D+03,
     4     0.3216834899377392D+02,  0.3614005342307463D+03,
     5     -.6653878500833375D+03,  -.6883582242804924D+03,
     6     0.2193362007156572D+04,  0.2423724600546293D+03,
     7     -.3665925878308203D+04,  0.2474933189642588D+04,
     8     0.1987663383445796D+04,  -.7382586600895061D+04,
     9     0.4991253411017503D+04,  0.1008505017740918D+05,
     a     -.1285284928905621D+05,  -.5153674821668470D+04,
     1     0.1301656757246985D+05,  -.4821250366504323D+04,
     2     -.4982112643422311D+04,  0.9694070195648748D+04,
     3     -.1685723189234701D+04,  -.6065143678129265D+04,
     4     0.2029510635584355D+04,  0.1244402339119502D+04,
     5     -.4336682903961364D+03,  0.8923209875101459D+02/
c
        data c1p1/
     1     -.4019450270734195D-23,  -.4819240943285824D-23,
     2     0.1087220822839791D-20,  0.1219058342725899D-21,
     3     -.7458149572694168D-19,  0.5677825613414602D-19,
     4     0.8351815799518541D-18,  -.5188585543982425D-17,
     5     0.1221075065755962D-15,  0.1789261470637227D-15,
     6     -.6829972121890858D-14,  -.1497462301804588D-14,
     7     0.1579028042950957D-12,  -.9414960303758800D-13,
     8     -.1127570848999746D-11,  0.3883137940932639D-11,
     9     -.3397569083776586D-10,  -.6779059427459179D-10,
     a     0.1149529442506273D-08,  0.4363087909873751D-09,
     1     -.1620182360840298D-07,  0.6404695607668289D-08,
     2     0.9651461037419628D-07,  -.1948572160668177D-06,
     3     0.6397881896749446D-06,  0.2318661930507743D-05,
     4     -.1983192412396578D-04,  -.1294811208715315D-04,
     5     0.2062663873080766D-03,  -.2867633324735777D-04,
     6     -.1084309075952914D-02,  0.1227880935969686D-02,
     7     0.2538406015667726D-03,  -.1153316815955356D-01/
c
        data c1p1b/       
     8     0.4520140008266983D-01,  0.5693944718258218D-01,
     9     -.9640790976658534D+00,  -.6517135574036008D+00,
     a     0.2051491829570049D+01,  -.1124151010077572D+01,
     1     -.3977380460328048D+01,  0.8200665483661009D+01,
     2     -.7950131652215817D+01,  -.3503037697046647D+02,
     3     0.9607320812492044D+02,  0.7894079689858070D+02,
     4     -.3749002890488298D+03,  -.8153831134140778D+01,
     5     0.7824282518763973D+03,  -.6035276543352174D+03,
     6     -.5004685759675768D+03,  0.2219009060854551D+04,
     7     -.2111301101664672D+04,  -.4035632271617418D+04,
     8     0.7319737262526823D+04,  0.2878734389521922D+04,
     9     -.1087404934318719D+05,  0.3945740567322783D+04,
     a     0.6727823761148537D+04,  -.1253555346597302D+05,
     1     0.3440468371829973D+04,  0.1383240926370073D+05,
     2     -.9324927373036743D+04,  -.6181580304530313D+04,
     3     0.6376198146666679D+04,  -.1033615527971958D+04,
     4     -.1497604891055181D+04,  0.1929025541588262D+04,
     5     -.4219760183545219D+02,  -.4521162915353207D+03/
c
        data c0p2/
     1     0.5641895835569398D+00,  -.5641895835321127D+00,
     2     -.7052370223565544D-01,  -.7052369923405479D-01,
     3     -.3966909368581382D-01,  0.3966934297088857D-01,
     4     0.4130698137268744D-01,  0.4136196771522681D-01,
     5     0.6240742346896508D-01,  -.6553556513852438D-01,
     6     -.3258849904760676D-01,  -.7998036854222177D-01,
     7     -.3988006311955270D+01,  0.1327373751674479D+01,
     8     0.6121789346915312D+02,  -.9251865216627577D+02,
     9     0.4247064992018806D+03,  0.2692553333489150D+04,
     a     -.4374691601489926D+05,  -.3625248208112831D+05,
     1     0.1010975818048476D+07,  -.2859360062580096D+05,
     2     -.1138970241206912D+08,  0.1051097979526042D+08,
     3     0.2284038899211195D+08,  -.2038012515235694D+09,
     4     0.1325194353842857D+10,  0.1937443530361381D+10,
     5     -.2245999018652171D+11,  -.5998903865344352D+10,
     6     0.1793237054876609D+12,  -.8625159882306147D+11,
     7     -.5887763042735203D+12,  0.1345331284205280D+13/
c
        data c0p2b/
     8     -.2743432269370813D+13,  -.8894942160272255D+13,
     9     0.4276463113794564D+14,  0.2665019886647781D+14,
     a     -.2280727423955498D+15,  0.3686908790553973D+14,
     1     0.5639846318168615D+15,  -.6841529051615703D+15,
     2     0.9901426799966038D+14,  0.2798406605978152D+16,
     3     -.4910062244008171D+16,  -.5126937967581805D+16,
     4     0.1387292951936756D+17,  0.1043295727224325D+16,
     5     -.1565204120687265D+17,  0.1215262806973577D+17,
     6     0.3133802397107054D+16,  -.1801394550807078D+17,
     7     0.4427598668012807D+16,  0.6923499968336864D+16/
c
c
        data c1p2/
     1     -.5641895835431980D+00,  -.5641895835508094D+00,
     2     0.2115710934750869D+00,  -.2115710923186134D+00,
     3     -.6611607335011594D-01,  -.6611615414079688D-01,
     4     -.5783289433408652D-01,  0.5785737744023628D-01,
     5     0.8018419623822896D-01,  0.8189816020440689D-01,
     6     0.1821045296781145D+00,  -.2179738973008740D+00,
     7     0.5544705668143094D+00,  0.2224466316444440D+01,
     8     -.8563271248520645D+02,  -.4394325758429441D+02,
     9     0.2720627547071340D+04,  -.6705390850875292D+03,
     a     -.3936221960600770D+05,  0.5791730432605451D+05,
     1     -.1976787738827811D+06,  -.1502498631245144D+07,
     2     0.2155317823990686D+08,  0.1870953796705298D+08,
     3     -.4703995711098311D+09,  0.3716595906453190D+07,
     4     0.5080557859012385D+10,  -.4534199223888966D+10,
     5     -.1064438211647413D+11,  0.8612243893745942D+11,
     6     -.5466017687785078D+12,  -.8070950386640701D+12,
     7     0.9337074941225827D+13,  0.2458379240643264D+13/
c
        data c1p2b/
     8     -.7548692171244579D+14,  0.3751093169954336D+14,
     9     0.2460677431350039D+15,  -.5991919372881911D+15,
     a     0.1425679408434606D+16,  0.4132221939781502D+16,
     1     -.2247506469468969D+17,  -.1269771078165026D+17,
     2     0.1297336292749026D+18,  -.2802626909791308D+17,
     3     -.3467137222813017D+18,  0.4773955215582192D+18,
     4     -.2347165776580206D+18,  -.2233638097535785D+19,
     5     0.5382350866778548D+19,  0.4820328886922998D+19,
     6     -.1928978948099345D+20,  0.1575498747750907D+18,
     7     0.3049162180215152D+20,  -.2837046201123502D+20,
     8     -.5429391644354291D+19,  0.6974653380104308D+20,
     9     -.5322120857794536D+20,  -.6739879079691706D+20,
     a     0.6780343087166473D+20,  0.1053455984204666D+20,
     1     -.2218784058435737D+20,  0.1505391868530062D+20/
c
c        if z is not in the right lower quadrant - bomb out
c
        ier=0
        com=z
        if( (rea(1) .ge. 0) .and. (rea(2) .le. 0) ) goto 1400
        ier=4
        return
 1400 continue
c
        done=1
        thresh1=4**2
        thresh2=8**2
        thresh3=20**2
c
c       check if if the user-specified z is in one of the 
c       intermediate regimes 
c
        d=z*dconjg(z)
        if( (d .lt. thresh1) .or. (d .gt. thresh3) ) goto 3000
c
c        if the user-specified z is in the first intermediate regime
c        (i.e. if its absolute value is between 4 and 8), act accordingly
c
        if(d .gt. thresh2) goto 2000
c
        cd=done/cdsqrt(z)
        cd=done/z
        m=35
        call hank101p(c0p1,m,cd,h0)
cccc         call prin2('after hank101p, h0=*',h0,2)
        h0=h0*cdexp(ima*z)/cdsqrt(z)
        h0=h0*z**18
c
        call hank101p(c1p1,m,cd,h1)
cccc         call prin2('after hank101p, h1=*',h1,2)
        h1=h1*cdexp(ima*z)/cdsqrt(z)
        h1=h1*z**18
        return
 2000 continue
c
c       z is in the second intermediate regime (i.e. its 
c       absolute value is between 8 and 20). act accordingly.
c
        cd=done/cdsqrt(z)
        cd=done/z
        m=27
c
        call hank101p(c0p2,m,cd,h0)
cccc         call prin2('after hank101p in second regime, h0=*',h0,2)
        h0=h0*cdexp(ima*z)/cdsqrt(z)
        m=31
        call hank101p(c1p2,m,cd,h1)
cccc         call prin2('after hank101p, h1=*',h1,2)
        h1=h1*cdexp(ima*z)/cdsqrt(z)
        return
 3000 continue
c
c
c        z is either in the local regime or the asymptotic one.
c        if it is in the local regime - act accordingly.
c
        if(d .gt. 50.d0) goto 4000
        call hank101l(z,h0,h1)
        return
c
c        z is in the asymptotic regime. act accordingly.
c
 4000 continue
        call hank101a(z,h0,h1)
        return
        end

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code and the beginning of the
c       least squares subroutines proper
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine cleastsq(a,u,w,t,n,m,ncols,rnorms,eps,v)
        implicit real *8 (a-h,o-z)
        doublecomplex  a(n,m),u(n,*),v(m,*),w(m,*),t(*)
        dimension rnorms(*)
c
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c
c   NOTE: this subroutine uses the subroutine leastsq1 (see) to perform 
c        almost all work. However, when m < n, it tends to be much 
c        more efficient than leastsq1, since it performs the two 
c        Gram-Schmidt processes in the optimal order (starting with 
c        the longer gram-schmidt involving shorter vectors)
c
c      
c                     input parameters: 
c
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c
c                     output parameters:
c
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first 
c        ncols columns of  u  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown, 
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1); On exit, it is structured as an
c        ncols * ncols matrix; however, on entry ncols is not
c        known, so it should be at least m*n complex *16 locations long.
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c
c    
c                      work arrays:
c
c  v - must be at least n*m real *8 locations long
c
c       . . . if m > n , construct the decomposition of the matrix as
c             specified by the user
c
        if(m .lt. n-2) goto 2000
c
        call cleast1(a,u,w,t,n,m,ncols,rnorms,eps,v)
cccc        call prin2('after first leastsq1, u=*',u,n*ncols)
c
        return
c
 2000 continue
c
c       n is greater than m. transpose the matrix, and decompose 
c       the transpose
c
        call cleastra(a,n,m,v)
        call cleascop(v,a,n*m)
c
        call cleast1(a,w,u,t,m,n,ncols,rnorms,eps,v)
c
cccc        call prin2('after leastsq1, u=*',u,n*ncols)
c       
c        transpose back everything that needs to be transposed back
c
        call cleastra(a,m,n,v)
        call cleascop(v,a,n*m)
c
        call cleastra(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c
        call cleasrer(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c
        call cleasrec(u,n,ncols,v)
        call cleascop(v,u,n*ncols)
c
c
c
        call cleasrec(w,m,ncols,v)
        call cleascop(v,w,m*ncols)
c
        call cleasrec(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c
        return
        end

c
c
c
c
c
        subroutine cleasts2(u,w,t,n,m,ncols,y,x,work,whts)
        implicit real *8 (a-h,o-z)
        doublecomplex  u(n,ncols),w(m,ncols),t(ncols,ncols)
        dimension whts(*)
        doublecomplex  x(m),y(n),work(ncols),d
c
c         this subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear 
c         system 
c
c                    A X = Y                                             (1)
c
c         The expansion used must have been produced by a prior call 
c         to the subroutine leastsq1 (see), and is of the form
c
c               A=U * T * W^*                                            (2)
c
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c      
c                     input parameters: 
c
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first 
c        ncols columns of  u  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown, 
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  n,m - the dimensionalities of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  y - the right-hand side in (1)
c
c                     output parameters:
c
c  x - the solution of the system (1) in the least squares sense
c    
c                      work arrays:
c
c  work - must be at least ncols real *8 locations long
c
c        . . . apply to the right-hand side the matrux  U^*
c
        do 1400 i=1,ncols
        d=0
        do 1200 j=1,n
cccc        d=d+u(j,i)*y(j)*whts(j)
cccccc        d=d+u(j,i)*y(j)
        d=d+dconjg(u(j,i))*y(j)
 1200 continue
        work(i)=d
 1400 continue
c
c       apply to the vector work the inverse of the matrix t
c
        x(1)=work(1)/t(1,1)
        do 2000 i=2,ncols
c
        d=0
        do 1600 j=1,i-1
        d=d+t(i,j)*x(j)
 1600 continue
c
        x(i)=(work(i)-d)/t(i,i)        
 2000 continue
c
        do 2200 i=1,ncols
        work(i)=x(i)
 2200 continue
c
c       multiply work by the matrix w
c
        do 2600 i=1,m
        d=0
        do 2400 j=1,ncols
        d=d+w(i,j)*work(j)
 2400 continue
        x(i)=d
 2600 continue
        return
        end
c
c
c
c
c
        subroutine cleastra(a,n,m,b)
        implicit real *8 (a-h,o-z)
        doublecomplex  a(n,m),b(m,n),x(*),y(*)
c
c       transpose a
c
        do 1400 i=1,n
        do 1200 j=1,m
cccc        b(j,i)=a(i,j)
        b(j,i)=dconjg(a(i,j))
 1200 continue
 1400 continue
        return
c
c
c
c
        entry cleascop(x,y,n)
c
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        end
c
c
c
c
c
        subroutine cleast1(a,u,w,t,n,m,ncols,rnorms,
     1    eps,v)
        implicit real *8 (a-h,o-z)
        doublecomplex  a(n,m),u(n,*),v(m,*),w(m,*),t(*),rnorms(*)
c
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c
c      
c                     input parameters: 
c
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c
c                     output parameters:
c
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first 
c        ncols columns of  u  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown, 
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c
c    
c                      work arrays:
c
c  v - must be at least n*m real *8 locations long
c
c        . . . using gram-schmidt process with pivoting, decompose 
c              the matrix a in the form
c
c          a=U  V^*,
c
c        with  u an orthogonal matrix of minimum rank, 
c        and v whatever it wants to be
c
        ifpivot=1
        call cleaspiv(a,n,m,u,v,ncols,rnorms,eps,ifpivot)
ccc         call prin2('after first leaspiv0, rnorms=*',rnorms,ncols)
c
c        using gram-schmidt process without pivoting, decompose 
c        the matrix v in the form
c
c          v=w t^*,
c
c        with  w an orthogonal matrix of minimum rank, 
c        and t a triangular matrix of dimensionality ncols * ncols
c
        ifpivot=0
        call cleaspiv(v,m,ncols,w,t,ncols2,rnorms,eps,ifpivot)
ccc         call prin2('after second leaspiv0, rnorms=*',rnorms,ncols2)


        if(2 .ne. 3) return
c
c       if the need be - restructure the matrix t
c
        if(ncols2 .eq. ncols) return
c
cccc        call leastres(t,v,t,ncols,ncols2)
        ncols=ncols2

        return
        end
c
c
c
c
c
c
        subroutine cmatve(a,x,y,n,m)
        implicit real *8 (a-h,o-z)
        doublecomplex  a(n,m),x(n),y(m),d
c
c       apply the matrix a to the vector x getting y
c
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c
c
c
c
c
        subroutine cmulback(u,w,t,n,ncols,m,c,work)
        implicit real *8 (a-h,o-z)
        doublecomplex  w(m,ncols),u(n,ncols),c(n,m),t(ncols,ncols),
     1      work(n,ncols),d
c
c       multiply u by t
c
        do 2600 i=1,n
        do 2400 j=1,ncols
        d=0
        do 2200 k=1,ncols
        d=d+u(i,k)*t(k,j)
 2200 continue
        work(i,j)=d
 2400 continue
 2600 continue
c
c        multiply work by w ^*
c
        do 3600 i=1,n
        do 3400 j=1,m
c
        d=0
        do 3200 k=1,ncols
        d=d+work(i,k)*dconjg(w(j,k))
 3200 continue
        c(i,j)=d
 3400 continue
 3600 continue
c
        return
        end
c
c
c
c
c
c
        subroutine cleaspiv(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        doublecomplex  a(n,m),b(n,m),v(m,*),prod
        dimension rnorms(*)
c
c       this matrix applies the compressing gram-schmidt 
c       process to the matrix a, obtaining its decomposition
c       in the form
c
c             a=b v^T,                                        (1)
c
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the 
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c      
c                     input parameters: 
c
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c
c                     output parameters:
c
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first 
c        ncols columns of  b  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown, 
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c
c        . . . copy the user-supplied matrix a into b
c
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c
c        apply the gram-schmidt proces (with pivoting) to b
c
         call cleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c
c        project the original matrix on the obtained orthogonal
c        columns
c
        do 2400 i=1,ncols
        do 2200 j=1,m
c
        call cleascap(a(1,j),b(1,i),n,prod)
cccc        v(j,i)=prod
        v(j,i)=dconjg(prod)
 2200 continue
 2400 continue
        return
        end
c
c
c
c
c
        subroutine cleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        doublecomplex  b(n,m),cd
        dimension rnorms(*)
c
c       this subroutine applies a pivoted double gram-schmidt 
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is 
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c
c                    input paramneters:
c
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy 
c  
c                     output parameters:
c
c  b - the matrix of gram-schmidt vectors of the matrix a. note 
c        that on exit from this subroutine only the first 
c        ncols columns of  b  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c
c       . . . prepare the array of values of norms 
c             of columns
c
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
cccc        d=d+b(j,i)**2
        d=d+b(j,i)*dconjg(b(j,i))
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c
c       . . . conduct gram-schmidt iterations
c
        thresh=dtot*eps**2
        do 4000 i=1,m
c
c       find the pivot  
c
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c
c       put the column number ipivot in the i-th place
c
        do 2600 j=1,n
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c
c       orthogonalize the i-th column to all preceeding ones
c
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c
        call cleascap(b(1,i),b(1,j),n,cd)
c
cccc        cd=dconjg(cd)
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c
c       normalize the i-th column
c
        call cleascap(b(1,i),b(1,i),n,cd)
c
c       if the norm of the remaining part of the matrix 
c       is sufficiently small - terminate the process
c
        d=cd
        if(d .lt. thresh) return
        ncols=i
c
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c
c        orthogonalize everything else to it 
c
        do 3200 j=i+1,m
c
        if(rnorms(j) .lt. thresh) goto 3200
c
        call cleascap(b(1,i),b(1,j),n,cd)
        cd=dconjg(cd)
c
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*dconjg(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue 
 3400 continue
 4000 continue

 4200 continue
c
        return
        end
c
c
c
c
c
        subroutine cleascap(x,y,n,prod)
        implicit doublecomplex  (a-h,o-z)
        dimension x(*),y(*)
c
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*dconjg(y(i))
 1200 continue
        return
        end
c
c
c
c
c
        subroutine cleasrec(a,n,m,b)
        implicit real *8 (a-h,o-z)
        doublecomplex  a(n,m),b(n,m)
c
c       reverse the columns of a 
c
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,m-i+1)
 1200 continue
 1400 continue
        return
c
c
c
c
        entry cleasrer(a,n,m,b)
c
c       reverse the rows of a 
c
        do 2400 i=1,m
        do 2200 j=1,n
        b(j,i)=a(n-j+1,i)
 2200 continue
 2400 continue
        return
        end

c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the legendre expansion routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c        This file contains a set of subroutines for the handling 
c        of Legendre expansions. It contains 19 subroutines that are 
c        user-callable. Following is a brief description of these 
c        subroutines.
c
c   legeexps - constructs Legendre nodes, and  corresponding Gaussian
c        weights. Also constructs the matrix v converting the 
c         coefficients of a legendre expansion into its values at 
c         the n Gaussian nodes, and its inverse u, converting the
c         values of a function at n Gaussian nodes into the
c         coefficients of the corresponding Legendre series.
c
c   legepol - evaluates a single Legendre polynomial (together
c         with its derivative) at the user-provided point
c
c   legepols - evaluates a bunch of Legendre polynomials
c         at the user-provided point
c
c   legeinmt - for the user-specified n, constructs the matrices of 
c        spectral indefinite integration differentiation on the n 
c        Gaussian nodes on the interval [-1,1]. 
c
c   legeinte - computes the indefinite integral of the legendre 
c        expansion polin getting the expansion polout
c
c   legediff -  differentiates the legendre expansion polin getting 
c        the expansion polout
c
c   legefder - computes the value and the derivative of a Legendre 
c        expansion at point X in interval [-1,1]; this subroutine 
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c
c   legefde2 - the same as legefder, except it is desigmed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c   
c   legeexev - computes the value of a Legendre expansion with 
c        at point X in interval [-1,1]; same as legefder, but does
c        not compute the derivative of the expansion
c
c   legeexe2 - the same as legeexev, except it is desigmed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c
c   lematrin - constructs the matrix interpolating functions from 
c        the n-point Gaussian grid on the interval [-1,1] to an 
c        arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c   
c   levecin - constructs the coefficients of the standard 
c        interpolation formula connecting the values of a 
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c
c   legeodev - evaluates at the point x a Legendre expansion
c        having only odd-numbered elements; this is a fairly 
c        efficient code, using external arrays that are 
c        precomputed
c
c   legeevev - evaluates at the point x a Legendre expansion
c        having only even-numbered elements; this is a fairly 
c        efficient code, using external arrays that are 
c        precomputed
c
c   legepeven - evaluates even-numbered Legendre polynomials 
c        of the argument x; this is a fairly efficient code, 
c        using external arrays that are precomputed
c
c   legepodd - evaluates odd-numbered Legendre polynomials 
c        of the argument x; this is a fairly efficient code, 
c        using external arrays that are precomputed
c
C   legefdeq - computes the value and the derivative of a
c        Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion; this subroutine 
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c
c   legeq - calculates the values and derivatives of a bunch 
c        of Legendre Q-functions at the user-specified point 
c        x on the interval (-1,1)
c
c   legeqs - calculates the value and the derivative of a single
c        Legendre Q-function at the user-specified point 
c        x on the interval (-1,1)
c
c
        subroutine legeexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        dimension x(*),whts(*),u(n,n),v(n,n)
c
c         this subroutine constructs the gaussiaqn nodes 
c         on the interval [-1,1], and the weights for the 
c         corresponding order n quadrature. it also constructs
c         the matrix v converting the coefficients
c         of a legendre expansion into its values at the n
c         gaussian nodes, and its inverse u, converting the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding legendre series.
c         no attempt has been made to make this code efficient, 
c         but its speed is normally sufficient, and it is 
c         mercifully short.
c
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c  
c                 output parameters:
c
c  x - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its 
c         legendre expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term legendre expansion into its values at
c         n legendre nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only 
c         if itype .ge. 1
c
c       . . . construct the nodes and the weights of the n-point gaussian 
c             quadrature
c
        ifwhts=0
        if(itype. gt. 0) ifwhts=1
        call legewhts(n,x,whts,ifwhts)
c
c       construct the matrix of values of the legendre polynomials
c       at these nodes        
c
        if(itype .ne. 2) return
        do 1400 i=1,n
c
        call legepols(x(i),n-1,u(1,i) )
 1400 continue
c
        do 1800 i=1,n
        do 1600 j=1,n
        v(i,j)=u(j,i)
 1600 continue
 1800 continue
c
c       now, v converts coefficients of a legendre expansion
c       into its values at the gaussian nodes. construct its 
c       inverse u, converting the values of a function at 
c       gaussian nodes into the coefficients of a legendre 
c       expansion of that function
c
        do 2800 i=1,n
        d=1
        d=d*(2*i-1)/2
        do 2600 j=1,n
        u(i,j)=v(j,i)*whts(j)*d
 2600 continue
 2800 continue
        return
        end
c
c
c
c
c
        subroutine legewhts(n,ts,whts,ifwhts)
        implicit real *8 (a-h,o-z)
        dimension ts(*),whts(*)
c
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on 
c        the interval [-1,1]
c
c                input parameters:
c
c  n - the number of nodes in the quadrature
c
c                output parameters:
c
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n) 
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c
c         use newton to find all roots of the legendre polynomial
c
        ts(n/2+1)=0
        do 2000 i=1,n/2
c
        xk=ts(i)
        deltold=1
        do 1400 k=1,10
        call legepol(xk,n,pol,der)
        delta=-pol/der
cccc         call prin2('delta=*',delta,1)
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=dabs(delta)
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c
c       now, use the explicit integral formulae 
c       to obtain the weights
c
        if(ifwhts .eq. 0) return
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call prodend(a,ts,n,i,fm)
        call prodend(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c
c
c
c
c
        subroutine legepol(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c
        pol=x
        der=1
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c
c
c
c
c
        subroutine prodend(x,xs,n,i,f)
        implicit real *8 (a-h,o-z)
        dimension xs(*)
c
c      evaluate the product
c
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c
        large=0
        do 2000 j=1,n
        dd=dabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
c
c
c
c
c
        subroutine legepols(x,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(*)
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=x
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w)
        implicit real *8 (a-h,o-z)
        dimension ainte(*),w(*),x(*),whts(*),adiff(*),endinter(*)
c
c
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes 
c        on the interval [-1,1]. Actually, this is omnly a 
c        memory management routine. All the actual work is done
c        by the subroutine legeinm0 (see)
c
c                           input parameters:
c
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION: 
c       itype=1 means that only the matrix ainte will 
c               be constructed
c       itype=2 means that only the matrix adiff will 
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c
c                           output paramaters:
c
c  ainte - the matrix of spectral indefinite integration on 
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on 
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c  endinter - the interpolation coefficients converting the 
c          values of a function at n Gaussian nodes into its
c          value at 1 (the right end of the interval)
c
c                           work arrays:
c
c  w - must be 3* n**2 + 2*n +50 *8 locations long
c
c        . . . allocate memory for the construction of the integrating
c              matrix
c
        ipolin=1
        lpolin=n+5
c
        ipolout=ipolin+lpolin
        lpolout=n+5
c
        iu=ipolout+lpolout
        lu=n**2+1
c
        iv=iu+lu
        lv=n**2+1
c
        iw=iv+lv
        lw=n**2+1
c
        ltot=iw+lw
c
c        construct the integrating matrix
c
        call legeinm0(n,ainte,adiff,w(ipolin),w(ipolout),
     1      x,whts,w(iu),w(iv),w(iw),itype,endinter)
c
        return
        end
c
c
c
c
c
        subroutine legeinm0(n,ainte,adiff,polin,polout,
     1      x,whts,u,v,w,itype,endinter)
        implicit real *8 (a-h,o-z)
        dimension ainte(n,n),u(n,n),v(n,n),w(n,n),
     1      endinter(*),x(n),whts(n),polin(n),polout(n),
     2      adiff(n,n)
c
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes 
c        on the interval [-1,1]
c
c                           input parameters:
c
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION: 
c       itype=1 means that only the matrix ainte will 
c               be constructed
c       itype=2 means that only the matrix adiff will 
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c
c                           output paramaters:
c
c  ainte - the matrix of spectral indefinite integration on 
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on 
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c
c                           work arrays:
c
c  polin, polout - must be n+3 real *8 locations each
c
c  u, v, w - must be n**2+1 real *8 locations each
c
c        . . . construct the matrices of the forward and inverse 
c              Legendre transforms
c
        itype2=2
        call legeexps(itype2,n,x,u,v,whts)
c
cccc         call prin2('after legeexps, u=*',u,n*n)
c
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the indefinite integral of that function
c
        if(itype. eq. 2) goto 2000
c
        do 1600 i=1,n
c
        do 1200 j=1,n+2
        polin(j)=0
 1200 continue
c
        polin(i)=1
c
        call legeinte(polin,n,polout)
c
        do 1400 j=1,n
        ainte(j,i)=polout(j)
 1400 continue
c
 1600 continue
c
cccc         call prin2('ainte initially is*',ainte,n*n)
c
c        multiply the three, obtaining the integrating matrix
c
        call matmul(ainte,u,w,n)
        call matmul(v,w,ainte,n)
c
 2000 continue
c
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the derivative of that function
c
        if(itype. eq. 1) goto 3000
c
        do 2600 i=1,n
c
        do 2200 j=1,n
        polin(j)=0
 2200 continue
c
        polin(i)=1
c
        call legediff(polin,n,polout)
c
        do 2400 j=1,n
        adiff(j,i)=polout(j)
cccc        ainte(i,j)=polout(j)
 2400 continue
c
 2600 continue
c
cccc         call prin2('adiff initially is*',adiff,n*n)
c
c        multiply the three, obtaining the integrating matrix
c
        call matmul(adiff,u,w,n)
        call matmul(v,w,adiff,n)
c
 3000 continue
c
c        construct the vector of interpolation coefficients
c        converting the values of a polynomial at the Gaussian
c        nodes into its value at the right end of the interval
c
        do 3400 i=1,n
c
        d=0
        do 3200 j=1,n
        d=d+u(j,i)
 3200 continue
        endinter(i)=d
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine legeinte(polin,n,polout)
        implicit real *8 (a-h,o-z)
        dimension polin(*),polout(*)
c
c       this subroutine computes the indefinite integral of the 
c       legendre expansion polin getting the expansion polout
c
c
c                       input parameters:
c
c  polin - the legendre expansion to be integrated
c  n - the order of the expansion polin 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c
c                       output parameters:
c
c  polout - the legendre expansion of the integral of the function 
c         represented by the expansion polin
c
        do 1200 i=1,n+2
        polout(i)=0
 1200 continue
c
        do 2000 k=2,n+1
        j=k-1
c
cccc        polout(k+1)=polin(k)/(2*j+1)+polout(k+1)
        polout(k+1)=polin(k)/(2*j+1)
        polout(k-1)=-polin(k)/(2*j+1)+polout(k-1)
c
 2000 continue
c
        polout(2)=polin(1)+polout(2)
c
        dd=0
        sss=-1
        do 2200 k=2,n+1
c
        dd=dd+polout(k)*sss
        sss=-sss
 2200 continue
c
        call prin2('dd=*',dd,1)
        polout(1)=-dd
c
        return
        end
c
c
c
c
c
        subroutine legediff(polin,n,polout)
        implicit real *8 (a-h,o-z)
        dimension polin(*),polout(*)
c
c       this subroutine differentiates the legendre 
c       expansion polin getting the expansion polout
c
c
c                       input parameters:
c
c  polin - the legendre expansion to be differentiated
c  n - the order of the expansion polin 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c
c                       output parameters:
c
c  polout - the legendre expansion of the derivative of the function 
c         represented by the expansion polin
c
        do 1200 k=1,n+1
        polout(i)=0
 1200 continue
c
        pk=polin(n+1)
        pkm1=polin(n)
        pkm2=0
        do 2000 k=n+1,2,-1
c
        j=k-1
c         
        polout(k-1)=pk*(2*j-1)
        if(k .ge. 3) pkm2=polin(k-2)+pk
c
        pk=pkm1
        pkm1=pkm2
c
 2000 continue
         return
         end
c
c
c
c
c
      SUBROUTINE legeFDER(X,VAL,der,pexp,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 pexp(*)
C
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients pexp
C     at point X in interval [-1,1].
c
c                input parameters:
c
C     X = evaluation point
C     pexp = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
C


        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c
        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c
      RETURN
      END


c
c
c
c
c
      SUBROUTINE legeFDE2(X,VAL,der,pexp,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 pexp(*),pjcoefs1(*),pjcoefs2(*)
c
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients pexp
C     at point X in interval [-1,1].
c
c                input parameters:
c
C  X - evaluation point
C  pexp - expansion coefficients
C  N  - order of expansion 
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call 
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below) 
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first 
c       call to this subroutine, ninit should be set to the maximum 
c       order n for which this subroutine might have to be called; 
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE 
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C  VAL - computed value
C  der - computed value of the derivative
C
C
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 1600 J = 2,N
c
cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2


        val=val+pexp(j+1)*pj
c
cccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
        derj=pjcoefs1(j)*(pjm1+x*derjm1)+pjcoefs2(j)*derjm2

ccc         call prin2('derj=*',derj,1)


cccc        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 1600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
      SUBROUTINE legeexev(X,VAL,pexp,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 pexp(*)
C
C     This subroutine computes the value o a Legendre
c     expansion with coefficients pexp at point X in interval [-1,1]
C
c                input parameters:
c
C     X = evaluation point
C     pexp = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C
        done=1
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c
        RETURN
        END
c
c
c
c
c
      SUBROUTINE legeexe2(X,VAL,pexp,N,
     1      pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 pexp(*),pjcoefs1(*),pjcoefs2(*)
c
C     This subroutine computes the value o a Legendre
c     expansion with coefficients pexp at point X in interval [-1,1]
C
c                input parameters:
c
C     X = evaluation point
C     pexp = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C
        done=1
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2

cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c
        RETURN
        END
c
c
c
c
c
        subroutine lematrin(n,m,xs,amatrint,ts,w)
        implicit real *8 (a-h,o-z)
        dimension amatrint(m,n),xs(*),w(*),ts(*)
c
c
c        This subroutine constructs the matrix interpolating
c        functions from the n-point Gaussian grid on the interval [-1,1]
c        to an arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c
c                 Input parameters:
c
c  n - the number of interpolation nodes
c  m - the number of nodes to which the functions will be interpolated
c  xs - the points at which the function is to be interpolated
c
c                  Output parameters:
c
c  amatrint - the m \times n matrix conerting the values of a function
c        at the n Legendre nodes into its values at m user-specified
c        (arbitrary) nodes
c  ts - the n Gaussian nodes on the interval [-1,1]
c
c                  Work arrays:
c
c  w - must be at least 2*n**2+n + 100 real *8 locations long 
c

        icoefs=1
        lcoefs=n+2
c
        iu=icoefs+lcoefs
        lu=n**2+10
c
        iv=iu+lu
c
        ifinit=1
        do 2000 i=1,m
c
        call levecin(n,xs(i),ts,w(iu),w(iv),w(icoefs),ifinit)
c
        do 1400 j=1,n
        amatrint(i,j)=w(j)
 1400 continue
c
        ifinit=0
 2000 continue
c
        return
        end

c
c
c
c
c
        subroutine levecin(n,x,ts,u,v,coefs,ifinit)
        implicit real *8 (a-h,o-z)
        dimension u(n,n),v(n,n),ts(*),coefs(*)
c
c        This subroutine constructs the coefficients of the 
c        standard interpolation formula connecting the values of a 
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c
c                 Input parameters:
c
c  n - the number of interpolation nodes
c  x - the points at which the function is to be interpolated
c  ts - the n Gaussian nodes on the interval [-1,1]; please note that
c        it is an input parameter only if the parameter ifinit (see 
c        below) has been set to 1; otherwise, it is an output parameter
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its 
c        legendre expansion; please note that
c        it is an input parameter only if the parameter ifinit (see 
c        below) has been set to 1; otherwise, it is an output parameter
c  ifinit - an integer parameter telling the subroutine whether it should
c        initialize the Legendre expander; 
c     ifinit=1 will cause the subroutine to perform the initialization
c     ifinit=0 will cause the subroutine to  skip the initialization
c
c                  Output parameters:
c
c  coefs - the interpolation coefficients
c
c                 Work arrays: 
c
c  v - must be at least n*n real *8 locations long
c
c       . . . construct the n Gausian nodes on the interval [-1,1]; 
c             also the corresponding Gaussian expansion-evaluation 
c             matrices
c
        itype=2
        if(ifinit .ne.0) call legeexps(itype,n,ts,u,v,coefs)
c
c       evaluate the n Legendre polynomials at the point where the
c       functions will have to be interpolated
c
        call legepols(x,n+1,v)
c
c       apply the interpolation matrix to the ector of values 
c       of polynomials from the right 
c
        call lematvec(u,v,coefs,n)
        return
        end
c
c
c
c
c
        subroutine lematvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        dimension a(n,n),x(n),y(n)
c
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(j,i)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c
c
c
c
c
        subroutine matmul(a,b,c,n)
        implicit real *8 (a-h,o-z)
        dimension a(n,n),b(n,n),c(n,n)
c
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1600 k=1,n
        d=d+a(i,k)*b(k,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
c
c
c
c
        entry matmua(a,b,c,n)
ccc          call prin2('in matmua, a=*',a,n**2)
ccc          call prin2('in matmua, b=*',b,n**2)
        do 3000 i=1,n
        do 2800 j=1,n
        d=0
        do 2600 k=1,n
        d=d+a(i,k)*b(j,k)
 2600 continue
        c(i,j)=d
 2800 continue
 3000 continue
ccc          call prin2('exiting, c=*',c,n**2)
        return
        end

c
c
c
c
c
        subroutine legeodev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension coepnm1(*),coepnp1(*),
     1            coexpnp1(*),coefs(*)
c
c
c       This subroutine evaluates at the point x a Legendre expansion
c       having only odd-numbered elements
c
c                  Input parameters:
c
c  x - point on the interval [-1,1] at which the Legendre expansion 
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - odd-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEPODD. IF these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  val - the value at the point x of the Legendre expansion with 
c       coefficients coefs (see above) 
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters
c 
c       
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c
        do 1200 nnn=2,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pi=x
        pip1=x*(2.5d0*x22-1.5d0)
c
        val=coefs(1)*pi+coefs(2)*pip1

        do 2000 i=1,nn/2-2
c
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pip1
c
        val=val+coefs(i+2)*pip2
c
        pi=pip1
        pip1=pip2


 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeevev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension coepnm1(*),coepnp1(*),
     1            coexpnp1(*),coefs(*)
c
c
c       This subroutine evaluates at the point x a Legendre expansion
c       having only even-numbered elements
c
c                  Input parameters:
c
c  x - point on the interval [-1,1] at which the Legendre expansion 
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - even-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEPEVEN. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  val - the value at the point x of the Legendre expansion with 
c       coefficients coefs (see above) 
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters
c       
        if(ninit .eq. 0) goto 1400
c
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pi=1
        pip1=1.5d0*x22-0.5d0
c
        val=coefs(1)+coefs(2)*pip1
c
c       n is greater than 2. conduct recursion
c
        do 2000 i=1,nn/2-2
c
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pip1
        val=val+coefs(i+2)*pip2
c
        pi=pip1
        pip1=pip2
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legepeven(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension pols(*),coepnm1(*),coepnp1(*),
     1            coexpnp1(*)
c
c       This subroutine evaluates even-numbered Legendre polynomials 
c       of the argument x, up to order nn+1
c
c                  Input parameters:
c
c  x - the argument for which the Legendre polynomials are 
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine ill initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c       PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEVEV. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c 
c                  Output parameters:
c
c  pols - even-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEVEV. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c 
c       
        if(ninit .eq. 0) goto 1400
c
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pols(1)=1
        pols(2)=1.5d0*x22-0.5d0
c
c       n is greater than 2. conduct recursion
c
        do 2000 i=1,nn/2
c
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pols(i+1)
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legepodd(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension pols(*),coepnm1(*),coepnp1(*),
     1            coexpnp1(*)
c
c       This subroutine evaluates odd-numbered Legendre polynomials 
c       of the argument x, up to order nn+1
c
c                  Input parameters:
c
c  x - the argument for which the Legendre polynomials are 
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  pols - the odd-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_1(x),
c       pols(2) = P_3(x), pols(3) = P_5 (x), etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS 
c       SUBROUTINE ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES 
c       SUSED BY THE UBROUTINE LEGEODEV. IF these arrays have been 
c       initialized by one of these two subroutines, they do not need 
c       to be initialized by the other one.
c       
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c
        do 1200 nnn=2,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pols(1)=x
        pols(2)=x*(2.5d0*x22-1.5d0)
c
        do 2000 i=1,nn/2
c
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pols(i+1)
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legefdeq(x,val,der,coefs,n)
        implicit real *8 (a-h,o-z)
        dimension coefs(*)
C
C     This subroutine computes the value and the derivative
c     of a Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion
c
c                input parameters:
c
C  X = evaluation point
C  coefs = expansion coefficients
C  N  = order of expansion 
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
c
        val=0
        der=0
c
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c
        pk=d
        pkp1=d*x-1

        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x 
c
        val=coefs(1)*pk+coefs(2)*pkp1
        der=coefs(1)*derk+coefs(2)*derkp1
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
c
        if(n .eq. 0) return
c
        return
 1200 continue
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
c
        derkm1=derk
        derk=derkp1
c
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
c
        val=val+coefs(k+2)*pkp1
        der=der+coefs(k+2)*derkp1
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeqs(x,n,pols,ders)
        implicit real *8 (a-h,o-z)
        dimension pols(*),ders(*)
c
c       This subroutine calculates the values and derivatives of 
c       a bunch of Legendre Q-functions at the user-specified point 
c       x on the interval (-1,1)
c
c                     Input parameters:
c
c  x - the point on the interval [-1,1] where the Q-functions and 
c       their derivatives are to be evaluated
c  n - the highest order for which the functions are to be evaluated
c  
c                     Output parameters:
c
c  pols - the values of the Q-functions (the evil twins of the 
c       Legeendre polynomials) at the point x (n+1 of them things)
c  ders - the derivatives of the Q-functions (the evil twins of the 
c       Legeendre polynomials) at the point x (n+1 of them things)
c  
c
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c
        pk=d
        pkp1=d*x-1

        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x 
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=pk
        ders(1)=derk
        if(n .eq. 0) return
c
        pols(2)=pkp1
        ders(2)=derkp1
        return
 1200 continue
c
        pols(1)=pk
        pols(2)=pkp1
c
c       n is greater than 2. conduct recursion
c
        ders(1)=derk
        ders(2)=derkp1
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
c
        derkm1=derk
        derk=derkp1
c
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
        ders(k+2)=derkp1
 2000 continue
c
        return
        end
c
c
c
c
c


        subroutine legeq(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the value and derivative of 
c       a Legendre Q-function at the user-specified point 
c       x on the interval (-1,1)
c
c
c                     Input parameters:
c
c  x - the point on the interval [-1,1] where the Q-functions and 
c       their derivatives are to be evaluated
c  n - the order for which the function is to be evaluated
c  
c                     Output parameters:
c
c  pol - the value of the n-th Q-function (the evil twin of the 
c       Legeendre polynomial) at the point x 
c  ders - the derivatives of the Q-function at the point x 
c  
c
        d= log( (1+x) /(1-x) ) /2
        pk=d
        pkp1=d*x-1
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pol=d

        der=(1/(1+x)+1/(1-x)) /2

        if(n .eq. 0) return
c
        pol=pkp1
        der=d + der *x 
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end

c
c
c
c
c
      SUBROUTINE legecFDE(x,val,der,pexp,n)
      implicit real *8 (a-h,o-z)
      complex *16  pexp(*),val,der
C
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with complex coefficients pexp
C     at point X in interval [-1,1].
c
c                input parameters:
c
C     X = evaluation point
C     pexp = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
C
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 j = 2,n
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c
        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
c
	SUBROUTINE PRINI(IP1,IQ1)
	CHARACTER *1 MES(1), AA(1)
	REAL *4 A(1)
	REAL *8 A2(1)
	INTEGER *4 IA(1)
	INTEGER *2 IA2(1)
	save
	IP=IP1
	IQ=IQ1
	RETURN

C
C
C
C
	ENTRY PRIN(MES,A,N)
	CALL  MESSPR(MES,IP,IQ)
	IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
	IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
	 RETURN
C
C
C
C
	ENTRY PRIN2(MES,A2,N)
	CALL MESSPR(MES,IP,IQ)
	IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
	IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
	RETURN
C
C
C
C
	ENTRY PRINF(MES,IA,N)
	CALL MESSPR(MES,IP,IQ)
	IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
	IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I8))
ccc 1600 FORMAT(10(1X,I7))
	RETURN
C
C
C
C
	ENTRY PRINF2(MES,IA2,N)
	CALL MESSPR(MES,IP,IQ)
	IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
	IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
	RETURN
C
C
C
C
	ENTRY PRINA(MES,AA,N)
	CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
	IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
	IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
	RETURN
	END
c
c
c
c
c
	SUBROUTINE MESSPR(MES,IP,IQ)
	CHARACTER *1 MES(1),AST
	DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
	I=0
	DO 1400 I=1,10000
	IF(MES(I).EQ.AST) GOTO 1600
	I1=I
 1400 CONTINUE
 1600 CONTINUE
	 IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
	 IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
	 RETURN
	 END
c
