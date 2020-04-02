
namespace funcs 
{
 

void ksol(int nright, int neq, int nsb, arma::vec& a, arma::vec& r, arma::vec&s, int& ising) {
		/*
		-----------------------------------------------------------------------

		Solution of a System of Linear Equations
		****************************************



		INPUT VARIABLES:

		nright,nsb       number of columns in right hand side matrix.
		for KB2D: nright=1, nsb=1
		neq              number of equations
		a()              upper triangular left hand side matrix (stored
		columnwise)
		r()              right hand side matrix (stored columnwise)
		for kb2d, one column per variable



		OUTPUT VARIABLES:

		s()              solution array, same dimension as  r  above.
		ising            singularity indicator
		0,  no singularity problem
		-1,  neq .le. 1
		k,  a null pivot appeared at the kth iteration



		PROGRAM NOTES:

		1. Requires the upper triangular left hand side matrix.
		2. Pivots are on the diagonal.
		3. Does not search for max. element for pivot.
		4. Several right hand side matrices possible.
		5. USE for ok and sk only, NOT for UK.


		c-----------------------------------------------------------------------
		*/

		float tol, piv, ap, ak;
		int i, j, k, kk, km1, lp, ll, ll1, llb, nm1, ii, ij, in, ijm, iv, m1, nn, nm;
	


		// If there is only one equation then set ising and return:


		if (neq <= 1)
		{
			ising = -1;
			return;
		}

		// Initialize:

		tol = 0.1e-06;
		ising = 0;
		nn = neq * (neq + 1) / 2;
		nm = nsb * neq;
		m1 = neq - 1;
		kk = 0;

		// Start triangulation:

		for (k = 1; k <= m1; ++k) {
			kk = kk + k;
			ak = a[kk];

			/*
			std::cout << "a = " << std::endl;
			a.print();
			std::cout << "ak = " << ak << std::endl;
			*/


			if (fabs(ak) < tol)
			{
				ising = k;
				return;
			}

			km1 = k - 1;
			for (iv = 1; iv <= nright; ++iv) 
			{
				nm1 = nm * (iv - 1);
				ii = kk + nn * (iv - 1);
				piv = 1/ a(ii);
				lp = 0;
				for (int i = k; i <= m1; ++i)
				{
					ll = ii;
					ii = ii + i;
					ap = a[ii] * piv;
					lp = lp + 1;
					ij = ii - km1;
					for (int j = i; j <= m1; ++j) 
					{
						ij = ij + j;
						ll = ll + j;
						a[ij] = a[ij] - ap * a[ll];
					}
					for (int llb = k; llb <= nm; llb = llb + neq) 
					{
						in = llb + lp + nm1;
						ll1 = llb + nm1;
						r[in] = r[in] - ap * r[ll1];
					}
				}
			}
		}



		// Error checking - singular matrix:

		ijm = ij - nn * (nright - 1);
		if (abs(a(ijm)) < tol)
		{

			ising = neq;
			return;
		}


		for (iv = 1; iv <= nright; ++iv) 
		{
			nm1 = nm * (iv - 1);
		 	 ij = ijm + nn * (iv - 1);
			piv = 1. / a[ij];
			for (llb = neq; llb <= nm; llb = llb + neq) {
				ll1 = llb + nm1;
			 s[ll1] = r[ll1] * piv;
			}
			i = neq;
			kk = ij;
			for (ii = 1; ii <= m1; ++ii) {
				kk = kk - i;
				piv = 1. / a[kk];
				i = i - 1;
				for (llb = i; llb <= nm; llb = llb + neq) {
					ll1 = llb + nm1;
					in = ll1;
					ap = r[in];
					ij = kk;
					for (j = i; j <= m1; ++j) {
						ij = ij + j;
						in = in + 1;
						ap = ap - a[ij] * s[in];
					}
					s[ll1] = ap * piv;
				}
			}
		}

		// Finished solving back, return:
	}

float sqdist(arma::vec point1, arma::vec point2, const arma::cube& rotmat, int irot) {
/*  This routine calculates the anisotropic distance between two points
    given the coordinates of each point and a definition of the
    anisotropy.

    This method only consider a single anisotropy senario.

    Parameters
    ----------
    point1 : tuple
        Coordinates of first point (x1,y1,z1)
    point2 : tuple
        Coordinates of second point (x2,y2,z2)
    rotmat : 3*3 ndarray
        matrix of rotation for this structure
	irot              The rotation matrix to use
    Returns
    -------
    sqdist : scalar
        The squared distance accounting for the anisotropy
        and the rotation of coordinates (if any).
    */

	

		float dx, dy, dz,cont;
		float  sqdist = 0.0;

		dx = point1[0] - point2[0];
		dy = point1[1] - point2[1];
		dz = point1[2] - point2[2];

		for (int i = 0; i < 3; i++) {
			cont = rotmat(i, 0,irot)  *  dx  +
				   rotmat(i,1,irot)   *  dy  +
				   rotmat(i,2,irot)   *  dz  ;
				sqdist += cont * cont;
								
		}
		return sqdist;

	};

float cova3(arma::vec point1, arma::vec point2,const Param& params,int ivarg, int irot) {
	/*	    Parameters
			 ----------
			point1, point2: tuple of 3
			coordinates of two points
	ivarg	number of nested structures (maximum of 4)

			Returns
			-------
			cova: scalar
			covariance between (x1,y1,z1) and (x2,y2,z2)
	*/
	float cova, ir, h, hr, pi, hsqd;
	int covtype;

		 pi = 3.14; 
		 hsqd = sqdist(point1, point2,params.rotmat,irot);
//Check for "zero" distance, return with cmax if so:
		if (hsqd < params.EPSLON) {
			cova = params.cmax;
			return cova;
		}
	
// looping over all structures
		cova = 0;
		for (int ist = 0; ist <= params.nst; ist++) {
			if (ist != 0) {
			//	ir = funcs::min((irot + ist - 1), params.MAXROT);
				hsqd = sqdist(point1, point2, params.rotmat,irot);				
			}

			// Compute the appropriate distance:
			h = std::sqrt(hsqd);

			covtype = params.it[ist]; // switch want integers specifically .. make u change param.it to integers

			switch (covtype) {
			case 1: // Spherical
				hr = h / params.aa_hmax[ist];   
				if (hr < 1) { cova = params.cc[ist] * (1. - hr * (1.5 - 0.5*hr*hr)); }
				break;
			case 2://gaussian
				cova = params.cc[ist] * exp(-3.0 * h / params.aa_hmax[ist]);  
				break;
			case 3://Exponential
				cova = params.cc[ist] *
					exp(-3.0 * (h / params.aa_hmax[ist]) * (h / params.aa_hmax[ist]));  
				break;
			case 4://Power
				cova = params.cmax - params.cc[ist] * (pow(h, (params.aa_hmax[ist]))); 
				break;
			case 5://Hole
				cova = params.cc[ist] - params.cc[ist] * cos(h / params.aa_hmax[ist] * pi); 
				break;
			}
			
		}
		return cova;
};

float gauinv(float p, float& ierr) {
	/*
	----------------------------------------------------------------------
	  Computes the inverse of the standard normal cumulative distribution
	  function with a numerical approximation from : Statistical Computing,
	  by W.J.Kennedy, Jr. and James E.Gentle, 1980, p. 95.

	 INPUT / OUTPUT :

		  p = float precision cumulative probability value : dble(psingle)
		 xp = G^-1 (p)in single precision
	   ierr = 1 - then error situation(p out of range), 0 - OK
	 ---------------------------------------------------------------------
	*/

	//  Coefficients of approximation :  
	float lim, p0, p1, p2, p3, p4, q0, q1, q2, q3, q4, xp, pp, y;
	lim = 1.0E-10;
	p0 = -0.322232431088;
	p1 = -1.0;
	p2 = -0.342242088547;
	p3 = -0.0204231210245;
	p4 = -0.0000453642210148;
	q2 = 0.531103462366;
	q0 = 0.0993484626060;
	q1 = 0.588581570495;
	q3 = 0.103537752850;
	q4 = 0.0038560700634;

	


	// Check for an error situation:
	ierr = 1;
	if (p < lim) 
	{
	xp = -1.0E10;
		return xp;
	}
	if (p >(1.0 - lim)) 
	{
		xp = 1.0E10;
		return xp;
	}
	
	ierr = 0;

	// Get k for an error situation:
	pp = p;
	if (p > 0.5) pp = 1 - pp;
	xp = float (0.0);
	if (p == double(0.5))
	{
		return xp;
	}

	// Approximate the function:
	y = sqrt(log(1.0 / (pp*pp)));
	xp = (float)(y + ((((y*p4 + p3)*y + p2)*y + p1)*y + p0) / ((((y*q4 + q3)*y + q2)*y + q1)*y + q0));
	if (float(p) == float(pp))
	{
		xp = -xp;
	}

	// Return with G^-1(p):
	return xp;
}

float gcum(float x) {	
		/*
		Evaluate the standard normal cdf given a normal deviate x.gcum is
		the area under a unit normal curve to the left of x.The results are
		accurate only to about 5 decimal places.
		*/
	float z,t,gcum,e2;
	if (x < 0) { z = -x; }
	else { z = x; }  // float check

		t = 1 / (1 + 0.2316419 * z);
		gcum = t*(0.31938153 + t*(-0.356563782 + t*(1.781477937 + 
			   t*(-1.821255978 + t*1.330274429))));
			
		
		if (z <= 6) { e2 = exp(-z*z / 2.)*0.3989422803;	}
		else {e2 = 0;}

		gcum = 1.0 - e2 * gcum;
		if (x >= 0) { return gcum; }	
		else { return 1.0 - gcum; }

}

float powint(float xlow, float xhigh, float ylow,  float yhigh, float xval, float potencia)
{
	float EPSLON =  1.0E-20;

	float valor;

	if ((xhigh - xlow) < EPSLON)
		valor = (yhigh + ylow) / float(2.0);
	else
		valor = ylow + (yhigh - ylow)*(float)pow(((xval - xlow) / (xhigh - xlow)), potencia);

	return valor;
}

float locate(arma::vec& xx, int n, int is, int ie, float x) {
	/*--------------------------------------------------------------------- -

		 Given an array "xx" of length "n", and given a value "x", this routine
		 returns a value "j" such that "x" is between xx(j) and xx(j + 1).
		 xx must be monotonic, either increasing or decreasing.
		 j = is - 1 or j = ie is returned to indicate that x is out of range.
		
	     Bisection Concept From "Numerical Recipes", Press et.al. 1986  pp 90.
		---------------------------------------------------------------------- -
		 Initialize lower and upper methods :
	*/


	int j, jl, ju, jm;

	if (is <= 0) { 
		is = 0 ;
	};
		jl = is-1;
	    ju = ie;
	
	if (xx(n) <= x) {
			j = ie; 
		return j;	
	}
	
		
		// If we are not done then compute a midpoint :

	L: if (ju - jl > 0)
		{
			jm = int((ju + jl) / 2);
		};

	//  Replace the lower or upper limit with the midpoint :

		if ((xx[ie] > xx[is])   ==   (x > xx[jm])  )
		{
			jl = jm;
		}		
		else 
		{
		ju = jm;
		}

		goto L;
	  j = jl; // Return with the array ivargex:
	return j;
}

float locate_c(arma::vec & xx, int  n, float x)
{
	float ju, jm, jl, j;
	int ascnd;

	jl = 0;      // lower
	ju = n + 1;  // upper
	ascnd = std::int16_t (xx[n] >= xx[1]);
	while (ju - jl > 1) {
		jm = int ( (ju + jl)/2 );
	
		if (x >= xx[jm] == ascnd)
		{
			jl = jm;
		}
		else
		{
			ju = jm;
		}
	}
	if (x == xx[1]) {
		j = 1;
	}
	else if (x == xx[n])
	{
		j = n - 1;
	}
	else { 
		j = jl;
	};

	
	return j;
}

void prints(std::string statement) {
	std::cout << statement <<std:: endl;
}

void sortem(arma::vec& vec1, arma::vec& vec2,int istart, int iend) {
	  arma::uvec  sort_ivargx;
	  int nsize = iend - istart + 1;
	  arma::vec tmp1(nsize), tmp2(nsize);

	  for (int i = istart; i <= iend; i++) {

		  tmp1[i-istart] = vec1[i];
		  tmp2[i-istart] = vec2[i];
	  }
	sort_ivargx = arma::sort_index(tmp1);
	tmp1 = tmp1(sort_ivargx);
	tmp2 = tmp2(sort_ivargx);

	for (int i = istart; i <= iend; i++) {

		 vec1[i] = tmp1[i - istart];
		 vec2[i] = tmp2[i - istart];
	}

	sort_ivargx = sort_ivargx + istart;
  }

double acorni(arma::vec& ixv) {
	  // I used the idum to know initialize ixv 

	  double MAXOP1, MAXINT, KORDEI, temp, Acorni;
	  KORDEI = 12;

	  MAXOP1 = KORDEI + 1;
	  MAXINT = 1073741824;


	  for (int i = 0; i < KORDEI; i++) {
		  ixv[i + 1] = ixv[i + 1] + ixv[i];
		  if (ixv[i + 1] >= MAXINT) {
			  ixv[i + 1] = ixv[i + 1] - MAXINT;
		  }
	  }

	  Acorni = double(ixv(KORDEI)) / MAXINT;
	  return Acorni;

	  // the func must be declared inside sgsim class
	  // cos the seeds must change all the time

  }

float max (float x, float y) {
	 return (x > y) ? x : y;
  }

float min(float x, float y) {
	 return (x < y) ? x : y;
 }

 template <typename T>
 T myMax(T x, T y) {
	 return (x > y) ? x : y;
 }

 template <typename T>
 T myMin(T x, T y) {
	 return (x < y) ? x : y;
 }
 
void sortem_c(float ib, float ie, arma::vec& a, int iperm,arma::vec& b) {

/*
-----------------------------------------------------------------------

                      Quickersort Subroutine
                      **********************

 This is a subroutine for sorting a real array in ascending order. This
 is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
 in collected algorithms of the ACM.

 The method used is that of continually splitting the array into parts
 such that all elements of one part are less than all elements of the
 other, with a third part in the middle consisting of one element.  An
 element with value t is chosen arbitrarily (here we choose the middle
 element). i and j give the lower and upper limits of the segment being
 split.  After the split a value q will have been found such that
 a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
 performs operations on the two segments (i,q-1) and (q+1,j) as follows
 The smaller segment is split and the position of the larger segment is
 stored in the lt and ut arrays.  If the segment to be split contains
 two or fewer elements, it is sorted and another segment is obtained
 from the lt and ut arrays.  When no more segments remain, the array
 is completely sorted.


 INPUT PARAMETERS:

   ib,ie        start and end index of the array to be sorteda
   a            array, a portion of which has to be sorted.
   iperm        0 no other array is permuted.
                1 array b is permuted according to array a
                2 arrays b,c are permuted.
                3 arrays b,c,d are permuted.
                4 arrays b,c,d,e are permuted.
                5 arrays b,c,d,e,f are permuted.
                6 arrays b,c,d,e,f,g are permuted.
                7 arrays b,c,d,e,f,g,h are permuted.
               >7 no other array is permuted.

   b,c,d,e,f,g,h  arrays to be permuted according to array a.

 OUTPUT PARAMETERS:

    a              = the array, a portion of which has been sorted.

    b,c,d,e,f,g,h  = arrays permuted according to array a (see iperm)

 NO EXTERNAL ROUTINES REQUIRED:

%-----------------------------------------------------------------------
*/

	arma::vec c, d, e, f, g, h;
	c.set_size(1); d.set_size(1);f.set_size(1); g.set_size(1), h.set_size(1);;

	// The dimensions for lt and ut have to be at least log (base 2) n
	int iring, lt[64], ut[64], i, j, k, m, p, q;
	float ta, tb, tc, td, te, tf, tg, th, xa, xb, xc, xd, xe, xf, xg, xh;

	// Initialize:
	j = ie;
	m = 1;
	i = ib;
	iring = iperm + 1;
	if (iperm > 7) iring = 1;

	// If this segment has more than two elements  we split it
A: if (j - i - 1 < 0)
	goto B;
   else if (j - i - 1 == 0)
	   goto D;
   else if (j - i - 1 > 0)
	   goto K;

   // p is the position of an arbitrary element in the segment we choose the
   // middle element. Under certain circumstances it may be advantageous
   // to choose p at random.
K: p = (j + i) / 2;
   ta = a[p];
   a[p] = a[i];
   if (iring >= 2) {
	   tb = b[p];
	   b[p] = b[i];
   }
   if (iring >= 3) {
	   tc = c[p];
	   c[p] = c[i];
   }
   if (iring >= 4) {
	   td = d[p];
	   d[p] = d[i];
   }
   if (iring >= 5) {
	   te = e[p];
	   e[p] = e[i];
   }
   if (iring >= 6) {
	   tf = f[p];
	   f[p] = f[i];
   }
   if (iring >= 7) {
	   tg = g[p];
	   g[p] = g[i];
   }
   if (iring >= 8) {
	   th = h[p];
	   h[p] = h[i];
   }

   //Start at the beginning of the segment, search for k such that a[k]>t
   q = j;
   k = i;
H: k = k + 1;
   if (k > q) goto F;
   if (a[k] <= ta) goto H;

   // Such an element has now been found now search for a q such that a[q]<t
   // starting at the end of the segment.
I:
   if (a[q] < ta) goto J;
   q = q - 1;
   if (q > k)     goto I;
   goto E;

   // a[q] has now been found. we interchange a[q] and a[k]
J: xa = a[k];
   a[k] = a[q];
   a[q] = xa;
   if (iring >= 2) {
	   xb = b[k];
	   b[k] = b[q];
	   b[q] = xb;
   }
   if (iring >= 3) {
	   xc = c[k];
	   c[k] = c[q];
	   c[q] = xc;
   }
   if (iring >= 4) {
	   xd = d[k];
	   d[k] = d[q];
	   d[q] = xd;
   }
   if (iring >= 5) {
	   xe = e[k];
	   e[k] = e[q];
	   e[q] = xe;
   }
   if (iring >= 6) {
	   xf = f[k];
	   f[k] = f[q];
	   f[q] = xf;
   }
   if (iring >= 7) {
	   xg = g[k];
	   g[k] = g[q];
	   g[q] = xg;
   }
   if (iring >= 8) {
	   xh = h[k];
	   h[k] = h[q];
	   h[q] = xh;
   }

   // Update q and search for another pair to interchange:
   q = q - 1;
   goto H;
E: q = k - 1;
F:

   // The upwards search has now met the downwards search:
   a[i] = a[q];
   a[q] = ta;
   if (iring >= 2) {
	   b[i] = b[q];
	   b[q] = tb;
   }
   if (iring >= 3) {
	   c[i] = c[q];
	   c[q] = tc;
   }
   if (iring >= 4) {
	   d[i] = d[q];
	   d[q] = td;
   }
   if (iring >= 5) {
	   e[i] = e[q];
	   e[q] = te;
   }
   if (iring >= 6) {
	   f[i] = f[q];
	   f[q] = tf;
   }
   if (iring >= 7) {
	   g[i] = g[q];
	   g[q] = tg;
   }
   if (iring >= 8) {
	   h[i] = h[q];
	   h[q] = th;
   }

   // The segment is now divided in three parts: (i,q-1),[q],(q+1,j)
   // store the position of the largest segment in lt and ut
   if (2 * q <= i + j) goto G;
   lt[m] = i;
   ut[m] = q - 1;
   i = q + 1;
   goto C;
G: lt[m] = q + 1;
   ut[m] = j;
   j = q - 1;

   //Update m and split the new smaller segment
C: m = m + 1;
   goto A;

   // We arrive here if the segment has  two elements we test to see if
   // the segment is properly ordered if not, we perform an interchange
D:
   if (a[i] <= a[j]) goto B;
   xa = a[i];
   a[i] = a[j];
   a[j] = xa;
   if (iring >= 2) {
	   xb = b[i];
	   b[i] = b[j];
	   b[j] = xb;
   }
   if (iring >= 3) {
	   xc = c[i];
	   c[i] = c[j];
	   c[j] = xc;
   }
   if (iring >= 4) {
	   xd = d[i];
	   d[i] = d[j];
	   d[j] = xd;
   }
   if (iring >= 5) {
	   xe = e[i];
	   e[i] = e[j];
	   e[j] = xe;
   }
   if (iring >= 6) {
	   xf = f[i];
	   f[i] = f[j];
	   f[j] = xf;
   }
   if (iring >= 7) {
	   xg = g[i];
	   g[i] = g[j];
	   g[j] = xg;
   }
   if (iring >= 8) {
	   xh = h[i];
	   h[i] = h[j];
	   h[j] = xh;
   }

   // If lt and ut contain more segments to be sorted repeat process:
B: m = m - 1;
   if (m <= 0) return;
   i = lt[m];
   j = ut[m];
   goto A;

}



}