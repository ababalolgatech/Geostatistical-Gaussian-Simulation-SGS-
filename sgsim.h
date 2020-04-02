#include"HelpFun.h"
#include"GEO.h"
#include<vector>


// In the next iteration , separate data and parameter class
// load parameters with Jason-like style
using namespace std;
class sgsim {

public:
	arma::Mat<float> simu;
	Param params;
	int ind; // for rotation matrix ...
	sgsim(Param& params_) : params(params_) {};
//	~sgsim();
	
	void check_params(Param&) {
		if (params.radius_hmax <= 0)		{
			std::cout << "Radius_hmax should be larger than zero" << std::endl;
		}
		if (params.nst <= 0) {
			std::cout << "nst must be at least 1" << std::endl;
		}
		if (params.sill != 1) {
			std::cout << "Sill is not equal to 1" << std::endl;
		}
	};

	void preprocess(Param&);
	void getindex(Param&, int ind,int&, int&, int&);
	void read_data(Param&);
	void set_rot(Param&, int);
	void ctable(Param&);
	void srchnd(Param&, float, float, float);
	inline void backtr(Param&, float , int , arma::vec& , arma::vec& , float& );
	void max_covariance(Param&,int);
	void simulate(Param&);
	void krigging(Param&,float,float,float,float, float,float);
	float acorni(Param&);
	void genrandpath(Param&,int,arma::vec&);

};

float sgsim::acorni(Param& params) {
	// I used the idum to know initialize ixv 

	double MAXOP1, MAXINT, KORDEI, temp, Acorni;
	KORDEI = 12;

	MAXOP1 = KORDEI + 1;
	MAXINT = 1073741824;


	for (int i = 0; i < KORDEI; i++)
	{
		params.ixv[i + 1] = params.ixv[i + 1] + params.ixv[i];
		if (params.ixv[i + 1] >= MAXINT)
		{
			params.ixv[i + 1] = params.ixv[i + 1] - MAXINT;
		}
	}

	Acorni = float(params.ixv(KORDEI)) / MAXINT;

	return Acorni;

	// the ixv parameters must change everytime
}

void sgsim::genrandpath(Param& params, int ns, arma::vec& order) {


	arma::vec sim;
	arma::vec tmp;

	order.zeros(ns);
	tmp.zeros(params.nxyz);
	sim.zeros(params.nxyz);

	for (int i = 1; i <= params.nxyz; i++)
	{
		sim[i] = sgsim::acorni(params);
		tmp[i] = i;
	}

	funcs::sortem_c(1, ns, sim, 1, tmp); // use order for now in debugging ... it should tmp

	for (int i = 0; i < ns; i++)
	{
		order[i] = tmp[i + 1];

	}
}
void sgsim::preprocess(Param&) {
	
	/* 
	   Calculate dimensional constants
	   sgsim_const = namedtuple('sgsim_const',
	   ['PMX', 'MAXNST', 'MAXDT', 'MAXSB',
		'MAXDIS', 'MAXSAM', 'UNEST'])
	*/

	params.maxbx = 1;
	if (params.nx > 1) { params.maxbx = int(params.nx / 2); }
	params.maxby = 1;
	if (params.ny > 1) { params.maxby = int(params.ny / 2); }
	params.maxbz = 1;
	if (params.nz > 1) { params.maxbz = int(params.nz / 2); }
	if (params.nz > 1) { params.maxbz = int(params.nz / 2); }

	params.PMX = 999,
	params.UNEST = -99.0;
	params.MAXNST = 4,
	params.MAXROT = 1 + params.MAXNST ;
	params.MAXDT = 9,
	params.MAXSB = (params.maxsbx, params.maxsby, params.maxsbz),
	params.MAXSAM = params.ndmax + 1;

	// GSLIB
	params.MAXCTX = params.mxctx;
	params.MAXCTY = params.mxcty;
	params.MAXCTZ = params.mxctz;
	params.MAXCXY = params.MAXCTX * params.MAXCTY;
	params.MAXXYZ = params.MAXCTX * params.MAXCTY * params.MAXCTZ;
	params.MAXX = params.nx;
	params.MAXY = params.ny;
	params.MAXZ = params.nz;
	params.MXYZ = params.MAXX * params.MAXY * params.MAXZ;
	if (params.MXYZ < 100) { params.MXYZ = 100; }
	params.MAXNOD = params.nodmax;
	params.MAXSAM = params.ndmax;
	params.MAXKR1 = params.MAXNOD + params.MAXSAM + 1;
	params.MAXKR2 = params.MAXKR1 * params.MAXKR1;
	params.MAXSBX = 1;
	if (params.nx > 1) { params.MAXSBX = int(params.nx / 2); }
	if (params.MAXSBX > 50) { params.MAXSBX = 50; }
	params.MAXSBY = 1;
	if (params.ny > 1) { params.MAXSBY = int(params.ny / 2); }
	if (params.MAXSBY > 50) { params.MAXSBY = 50; }
	params.MAXSBZ = 1;
	if (params.nz > 1) {params.MAXSBZ = int(params.nz / 2);	}
	if (params.MAXSBZ > 50) { 	params.MAXSBZ = 50; }
	params.MAXSB = params.MAXSBX * params.MAXSBY * params.MAXSBZ;
		arma::mat testData;
		testData.load(params.datafl);
		params.nd = testData.n_rows;
		int MAXDAT = params.nd;
		// pre-allocation
		std::cout << "...Pre-allocation..." << std::endl ;
		params.x.zeros(MAXDAT);
		params.y.zeros(MAXDAT);
		params.z.zeros(MAXDAT);
		params.vr.zeros(MAXDAT);
		params.wt.zeros(MAXDAT);
		params.vrtr.zeros(MAXDAT);
		params.vrgtr.zeros(MAXDAT);
		params.vrg.ones(MAXDAT);  // needs debugging  - used in backtr func
		params.close.zeros(MAXDAT);  // close samples from super search
		params.sec.zeros(MAXDAT);
		params.sim.zeros(params.MXYZ + 1);
		params.simdebug.zeros(params.MXYZ+1,params.nsynth + 1);
		params.synthdebug.zeros(params.MXYZ + 1,params.nsynth + 1);
		params.simu.zeros(params.MXYZ);   // add number of simulations later
		params.lvm.zeros(params.MXYZ);
		params.tmp.zeros(params.MAXXYZ);
		int MAXORD = params.MXYZ;

		if (MAXORD < params.MAXCXY) { MAXORD = params.MAXCXY; }
		params.order.zeros(MAXORD);
		params.covtab.zeros(params.MAXCTX+1, params.MAXCTY+1, params.MAXCTZ+1);
		params.cnodex.zeros(params.MAXNOD + 1);
		params.cnodey.zeros(params.MAXNOD + 1);
		params.cnodez.zeros(params.MAXNOD + 1);
		params.cnodev.zeros(params.MAXNOD + 1);
		params.vra.zeros(params.MAXKR1);
		params.vrea.zeros(params.MAXKR1);
		params.r.zeros(params.MAXKR1);
		params.s.zeros(params.MAXKR1);
		params.a.zeros(params.MAXKR2);
		params.nisb.zeros(params.MAXSB);
		params.icnode.zeros(params.MAXXYZ + 1);
		params.ixnode.zeros(params.MAXXYZ + 1);
		params.iynode.zeros(params.MAXXYZ + 1);
		params.iznode.zeros(params.MAXXYZ + 1);
		params.ixsbtosr.zeros(8*  params.MAXSB + 1);
		params.iysbtosr.zeros(8 * params.MAXSB + 1);
		params.izsbtosr.zeros(8 * params.MAXSB + 1);

		//params.rotmat.zeros(params.MAXROT, 3, 3);
		params.rotmat.zeros(3, 3, params.MAXROT); // MAXROT = MAXNST + 1, MAXNST = 4
		params.ixv.zeros(13);
		params.ixv[0] = params.seed;

		
	
}
void sgsim::getindex(Param&,int ind, int&x_index, int& y_index, int&z_index) {
	/*
-----------------------------------------------------------------------

     Gets the coordinate index location of a point within a grid
     ***********************************************************

 n       number of "nodes" or "cells" in this coordinate direction
 min     origin at the center of the first cell
 siz     size of the cells
 loc     location of the point being considered
 ind	 output index within [1,n]
 inflag  true if the location is actually in the grid (false otherwise
         e.g., if the location is outside then index will be set to
         nearest boundary

-----------------------------------------------------------------------
	*/	
	// I think the 1.5 is there because int rounds down 

	x_index = int( ((params.x[ind]-params.xmn)/params.xsiz) + 1.5);
	if (x_index < 1) 
	{	x_index = 1;		
		//std::cout << "X_index out of bounds" << std::endl;
	}
	else if (x_index > params.nx)
	{
		x_index = params.nx;
		//std::cout << "X_index out of bounds" << std::endl;
	};
	

	y_index = int( ((params.y[ind]-params.ymn)/params.ysiz) + 1.5);
	if (y_index < 1) 
	{	
	    y_index = 1;
		//std::cout << "Y_index out of bounds" << std::endl;	
	}
	else if (y_index > params.ny) {
		y_index = params.ny;
		//std::cout << "Y_index out of bounds" << std::endl;
	};

	z_index = int( ((params.z[ind]-params.zmn)/params.zsiz) + 1.5);
	if (z_index < 1)
	{		
		z_index = 1;
	//	std::cout << "Z_index out of bounds" << std::endl;
	}
	else if (z_index > params.nz) 
	{		
		z_index = params.nz;
		//std::cout << "Z_index out of bounds" << std::endl;
	};
};
void sgsim::read_data(Param&) {
	int ntr, nt, j, testfl, icolwt, icolvr;
	float vrr, twt, cp, oldcp, w, vrg, ierr;
	double p;
	std::cout << "NOTE: indexing starts from 0" << std::endl;

	// Perform some quick error checking:
	testfl = 0;
	if (params.ltail != 1 && params.ltail != 2) {
		std::cout << "ERROR invalid lower tail option %d \n" << params.ltail << std::endl;
		std::cout << "      only allow 1 or 2 " << std::endl;
		testfl = 1;
	}
	if (params.utail != 1 && params.utail != 2 && params.utail != 4) {
		std::cout << "ERROR invalid upper tail option %d \n" << params.utail << std::endl;
		std::cout << "      only allow 1,2 or 4 !" << std::endl;
		testfl = 1;
	}
	if (params.utail == 4 && params.utpar < 1.0) {
		std::cout << "ERROR invalid power for hyperbolic tail %f \n" << params.utpar << std::endl;
		std::cout << "      must be greater than 1.0!" << std::endl;
		testfl = 1;
	}
	if (params.ltail == 2 && params.ltpar < 0.0) {
		std::cout << "ERROR invalid power for power model %f \n" << params.ltpar << std::endl;
		std::cout << "      must be greater than 0.0!" << std::endl;
		testfl = 1;
	}
	if (params.utail == 2 && params.utpar < 0.0) {
		std::cout << "ERROR invalid power for power model %f \n" << params.utpar << std::endl;
		std::cout << "      must be greater than 0.0!" << std::endl;
		testfl = 1;
	}
	if (testfl == 1) exit(0);

	//Establish the reference histogram for the simulation (provided that
	// we have data, and we are transforming the data)
	if (params.itrans == 1)
	{
		std::cout << "Setting up transformation table..." << std::endl;



		if (params.ismooth == 1) {
			//decide which file to use for establishing the transformation table
			params.tmpfl = params.smthfl;
			icolvr = params.isvr;
			icolwt = params.iswt;
		}
		else { //data histogram, possibly with declustering weights is used for transformation.
			params.tmpfl = params.datafl;
			icolvr = params.ivrl;
			icolwt = params.iwt;
		}

		/***** read in file for transformation table*************/
		params.var.load(params.tmpfl);
		// Keep this data: Assign the data value and coordinate location:
		params.vrtr = params.var.col(icolvr - 1);
		int sz = params.vrtr.size();
		nt = 0;
		twt = 0.0;
		for (int ntr = 0; ntr < sz; ntr++) {  // dealing with the file

			if (icolwt <= 0) 
			{ params.vrgtr[ntr] = 1.0;
			}
			else { params.vrgtr[ntr] = params.var.col(icolwt - 1)[ntr]; }

			// Trim this data ?
			if (params.var[icolvr - 1] < params.tmin || params.var[icolvr - 1] > params.tmax)
			{
				nt = nt + 1;
				continue;
			}
			ntr = ntr + 1;

			if (params.vrgtr[ntr] <= 0.0) {
				ntr = ntr - 1;
				nt = nt + 1;
				continue;
			}
			twt = twt + params.vrgtr[ntr];
		}
		
	
	// what is vrtr & vrgtr ?
	//	funcs::sortem(params.vrtr, params.vrgtr, sort_ind);  
	// The data is not sorted in GSLIB - don't know why sortem func was used

	//Compute the cumulative probabilities and write transformation table
	twt = funcs::max(twt, params.EPSLON); 
	ntr = params.nd;
	std::cout << ntr << std::endl;
	oldcp = 0.0;
	cp = 0.0;
	float iend = -800000; 
	std::cout << "iend is a crazy number & I still don't know the purpose" << std::endl;
	for (j = 0; j <= iend; ++j) {
		cp = cp + float(params.vrgtr[j] / twt);
		w = (cp + oldcp)*0.5;
		vrg = funcs::gauinv(w, ierr);
		if (ierr == 1) { vrg = params.UNEST; }
		oldcp = cp;
		// Now, reset the weight to the normal scores value
		params.vrgtr[j] = vrg;
		//std::cout << "cp = " << cp << std::endl;
		//std::cout << "w = " << w << std::endl;
		//std::cout << "vrg= " << vrg << std::endl;
	}// vrtr - transformation table

}

	/***** read in hard data for simulation*************/
	std::cout << "Reading Data..." << std::endl;
	params.var.load(params.datafl);

	if (params.ixl < 0) { params.x.ones(ntr)*params.xmn; }
	else { params.x = params.var.col(params.ixl-1); }
	if (params.iyl < 0) { params.y.ones(ntr)*params.ymn; }
	else { params.y = params.var.col(params.iyl-1); }
	if (params.izl < 0) { params.z.ones(ntr)*params.zmn; }
	else{ params.z = params.var.col(params.izl-1);  }
	

	if (params.ivrl >= 0) {
		params.vr = params.var.col(params.ivrl-1);  
	}
	if (params.iwt >= 0) {
		params.wt = params.var.col(params.iwt);
	}
	if (params.isecvr >= 0) {
		params.sec = params.var.col(params.isecvr);
	}



	if (params.itrans ==1 ) {
		for (int nd = 0; nd<ntr; nd++) {
			if (params.itrans == 1);{
				vrr = params.vr[nd]; };
		
			//cout <<"vrtr = "<<params.vrtr[ntr - 1] << endl;
			//j = funcs::locate(params.vrtr, ntr-1,0, ntr-1, vrr);
			j = funcs::locate_c(params.vrtr,ntr-1, vrr);
		
			j = funcs::min(funcs::max(1, j), (ntr - 2));
		vrg = funcs::powint(params.vrtr[j],params.vrtr[j + 1], 
		params.vrgtr[j], params.vrgtr[j + 1], vrr, 1.0);
		params.vr[nd] = vrg;
		/*
		std::cout << "nd = " << nd << std::endl;
		std::cout << "j = " << j << std::endl;
		std::cout << "vrr = " << vrr << std::endl;
		std::cout << "vr[nd] = " << vrg<< std::endl;
		std::cout << "vrgtr[j] = " << params.vrgtr[j] << std::endl;
		std::cout << "vrgtr[j+1] = " << params.vrgtr[j+1] << std::endl;
		std::cout << "---------------------------------------------"<< std::endl;
		*/
		}
	}
	
	// random number generator (Acorni)
	for (int i = 0; i < 1000; i++) 
	{	
		p = sgsim::acorni(params);
	}



	// reading inversion data
	params.seis_grid.load(params.seisfl);
	params.wavlet.load(params.wavefl);

	};
void sgsim::set_rot(Param& params, int ind) {

	/*
	              Sets up an Anisotropic Rotation Matrix
              **************************************

 Sets up the matrix to transform cartesian coordinates to coordinates
 accounting for angles and anisotropy (see manual for a detailed
 definition):


 INPUT PARAMETERS:

   ang1             Azimuth angle for principal direction
   ang2             Dip angle for principal direction
   ang3             Third rotation angle
   anis1            First anisotropy ratio
   anis2            Second anisotropy ratio
   ind              matrix indicator to initialize
   MAXROT           maximum number of rotation matrices dimensioned
   rotmat           rotation matrices
	*/
	float DEG2RAD, alpha, beta, theta;
	float sina, sinb, cosb, cosa, sint, cost;
	float afac1, afac2, neq;

	/*
	Converts the input angles to three angles which make more
	 mathematical sense :
	
		         alpha   angle between the major axis of anisotropy and the
		         E - W axis.Note : Counter clockwise is positive.
		         beta    angle between major axis and the horizontal plane.
				(The dip of the ellipsoid measured positive down)
		         theta   Angle of rotation of minor axis about the major axis
		        of the ellipsoid.
		*/
	DEG2RAD = 3.141592654 / 180.0;
	if (ind > 3)  //  to compute rotmat for search ellipsoid
	{ 
		if (params.sang1 >= 0.0 & params.sang1 < 270) {
			alpha = (90 - params.sang1) * DEG2RAD;
		}
		else {
			alpha = (450 - params.sang1) * DEG2RAD;
		}
		beta = -1.0 * params.sang2 * DEG2RAD;
		theta = params.sang3* DEG2RAD;

		afac1 = 1.0 / float(funcs::max(params.sanis1, params.EPSLON));
		afac2 = 1.0 / float(funcs::max(params.sanis2, params.EPSLON));
	}
	else {
		if (params.ang1[ind] >= 0.0 & params.ang1[ind] < 270) {
			alpha = (90 - params.ang1[ind]) * DEG2RAD;
		}
		else {
			alpha = (450 - params.ang1[ind]) * DEG2RAD;
		}
		beta = -1.0 * params.ang2[ind] * DEG2RAD;
		theta = params.ang3[ind] * DEG2RAD;

		afac1 = 1.0 / float(funcs::max(params.anis1[ind], params.EPSLON));
		afac2 = 1.0 / float(funcs::max(params.anis2[ind], params.EPSLON));
	
	}



		// Get the required sines and cosines:
		sina = float(std::sin(alpha));
		sinb = float(std::sin(beta));
		sint = float(std::sin(theta));
		cosa = float(std::cos(alpha));
		cosb = float(std::cos(beta));
		cost = float(std::cos(theta));

		//  Construct the rotation matrix in the required memory:



		params.rotmat(0, 0, ind) = float(cosb * cosa);   // note its ijk 
		params.rotmat(0, 1, ind) = float(cosb * sina);
		params.rotmat(0, 2, ind) = float(-sinb);

		params.rotmat(1, 0, ind) = float(afac1 * (-cost * sina + sint * sinb*cosa));
		params.rotmat(1, 1, ind) = float(afac1 * (cost * cosa + sint * sinb * sina));
		params.rotmat(1, 2, ind) = float(afac1 * (sint * cosb));

		params.rotmat(2, 0, ind) = float(afac2 * (sint * sina + cost * sinb * cosa));
		params.rotmat(2, 1, ind) = float(afac2 * (-sint * cosa + cost * sinb * sina));
		params.rotmat(2, 2, ind) = float(afac2 * (cost * cosb));
		/*
		
		std::cout << "------------debugging-------------" << std::endl;
		
		cout << "sina = " << sina << endl;
		cout << "sinb = " << sinb << endl;
		cout << "sint = " << sint << endl;
		cout << "cosa = " << cosa << endl;
		cout << "cosb = " << cosb << endl;
		cout << "cost = " << cost << endl;	

		cout << "afac1 = " << afac1 << endl;
		cout << "afac2 = " << afac2 << endl;
	
		std::cout << "------------debugging-------------" << std::endl;
		cout << "params.rotmat(0, 0, ind) = " << params.rotmat(0, 0, ind) << endl;
		cout << "params.rotmat(0, 1, ind) =  " << params.rotmat(0, 1, ind) << endl;
		cout << "params.rotmat(0, 2, ind) =  " << params.rotmat(0, 2, ind) << endl;

		cout << "params.rotmat(1, 0, ind) = " << params.rotmat(1, 0, ind) << endl;
		cout << "params.rotmat(1, 1, ind) = " << params.rotmat(1, 1, ind) << endl;
		cout << "params.rotmat(1, 2, ind) = " << params.rotmat(1, 2, ind) << endl;

		cout << "params.rotmat(2, 0, ind) =  " << params.rotmat(2, 0, ind) << endl;
		cout << "params.rotmat(2, 1, ind) =  " << params.rotmat(2, 1, ind) << endl;
		cout << "params.rotmat(2, 2, ind) =  " << params.rotmat(2, 2, ind) << endl;


	     params.rotmat.print();
		 */
		
} 
void sgsim::ctable(Param&) {
	/*  Create covariance lookup table

		The idea is to establish a 3-D network that contains the covariance
		value for a range of grid node offsets that should be at as large
		as twice the search radius in each direction.  The reason it has to
		be twice as large as the search radius is because we want to use it
		to compute the data covariance matrix as well as the data-point
		covariance matrix.
		Secondly, we want to establish a search for nearby nodes that
		in order of closeness as defined by the variogram.

		OUTPUT: self.covtab						*/
	std::cout << "indexing starts @ 1 - just for my sanity for now " << std::endl;

	float tiny = 1.0e-10;
	float cmax, hsqd,temp,tmp1,tmp2;
	int xx, yy, zz, loc;
	int ix, iy, iz;
	int ic, jc, kc, ncts;
	arma::uvec sort_indx;
	arma::vec point1, point2, tmp, order;
	tmp.set_size(params.MAXXYZ);
	order.set_size(params.MAXXYZ) ;
	params.nctx = std::min(((params.mxctx - 1) / 2), (params.nx-1 ));
	params.ncty = std::min(((params.mxcty - 1) / 2), (params.ny-1 ));
	params.nctz = std::min(((params.mxctz - 1) / 2), (params.nz-1 ));
	params.nlooku = 0;
	point1.zeros(3, 1);
	point2.zeros(3, 1);

	std::cout << "Covariance Look up table and search for previously" << std::endl;
	std::cout << "simulated grid nodes.  The maximum range in each" << std::endl;
	std::cout << "coordinate direction for covariance look up is" << std::endl;
	std::cout << "X direction:" << params.nctx * params.xsiz << std::endl;
	std::cout << "Y direction:" << params.ncty * params.ysiz << std::endl;
	std::cout << "Z direction:" << params.nctz * params.zsiz << std::endl;
	std::cout << "Node Values are not searched beyond this distance!" << std::endl;


	cmax = params.cmax; // already computed in sgsim::maxcovariance function
	int irot = 0 ;
	int ivarg = 0;
	for (int i = -params.nctx; i <= params.nctx; i++) // investigate why <=
	{
		xx = i * params.xsiz;
		ic = params.nctx + i + 1;
		for (int j = -params.ncty; j <= params.ncty; j++)
		{
			yy = j * params.ysiz;
			jc = params.ncty + j + 1;
			for (int k = -params.nctz; k <= params.nctz; k++)
			{				
				zz = k * params.zsiz;
				kc = params.nctz + k + 1;
				point2[0] = xx; 
				point2[1] = yy;
				point2[2] = zz;
				temp = funcs::cova3(point1, point2,params,irot, ivarg);
			/*
				std::cout << "ic = " << ic << endl;
				std::cout << "jc = " << jc << endl;
				std::cout << "kc = " << kc << endl;
				std::cout << "-----------------------" << endl;
			*/
				params.covtab(ic, jc, kc) = temp;
				hsqd = funcs::sqdist(point1, point2, params.rotmat,params.isrot);
				if (hsqd <= params.radsq) {  
					params.nlooku = params.nlooku + 1;
				// We want to search by closest variogram distance (and use the
				// anisotropic Euclidean distance to break ties:
				//	std::cout << "hsqd = " << hsqd << endl;
					tmp1  = -(params.covtab(ic, jc, kc) - params.TINY*std::real(hsqd));
					tmp2 = std::real((kc - 1)*params.MAXCXY + (jc - 1)*params.MAXCTX + ic);

					/*
					std::cout << "tmp1 = " << tmp1 << endl;
					std::cout << "tmp2 = " << tmp2 << endl;
					std::cout << "-----------------------" << endl;
					*/
					 tmp[params.nlooku] = tmp1;
				   order[params.nlooku] = tmp2;			  
				}
			}
		}
	}

	    /* Finished setting up the look - up table, now order the nodes such
		   that the closest ones, according to variogram distance, are searched
		   first.Note: the "loc" array is used because I didn't want to make
		   special allowance for 2 byte integers in the sorting subroutine  
		 */

	funcs::sortem_c(1, params.nlooku, tmp, 1, order);

	for (int il = 1; il <= params.nlooku; il++) 
	{
		loc = int(order[il]);
		iz = int((loc-1)/params.MAXCXY) + 1;
		iy = int((loc-(iz-1)*params.MAXCXY-1) / params.MAXCTX) + 1;
		ix = loc - (iz-1)*params.MAXCXY - (iy-1)*params.MAXCTX;
		params.iznode(il) = int(iz);
		params.iynode(il) = int(iy);
		params.ixnode(il) = int(ix);	

	}

	if (params.nodmax > params.MAXNOD)
	{
		std::cout << "The maximum number of close nodes =" << params.nodmax << std::endl;
		std::cout << "Tthis was reset from your specification due to storage limitation" << std::endl;
		params.nodmax = params.MAXNOD;
	}

};
void sgsim::srchnd(Param&, float ix, float iy, float iz) {
	/*      
	Search for nearby Simulated Grid nodes in a Spiral fashion.

			The idea is to spiral away from the node being simulated and note all
			the nearby nodes that have been simulated.

			INPUT VARIABLES:
			ix,iy,iz:				 index of the point currently being simulated
			sim:					 the realization so far
			nodmax:					 the maximum number of nodes that we want
			nlookup:				 the number of nodes in the look up table
			ixnode, iynode, iznode:  the relative indices of those nodes.
			xmn, ymn, zmn:			 the origin of the global grid netwrok
			xsiz, ysiz, zsiz:        the spacing of the grid nodes.

			OUTPUT VARIABLES:
			  ncnode          the number of close nodes
			  icnode()        the number in the look up table
cnodex, cnodey, cnodez:       list the location of the nodes
			  cnodev()        the values at the nodes		
			  */
	int i, j, k, ind;
	int idx, idy, idz, iq;

	params.ninoct.zeros(8);
	params.ncnode = 0;
	 for (int il = 2; il <= params.nlooku; il++) {
		if (params.ncnode == params.nodmax)
		{
			return;
		} // exit fun
		 // calculate location of point in simulation grid
		i = ix + params.ixnode[il] - params.nctx - 1;
		j = iy + params.iynode[il] - params.ncty - 1;
		k = iz + params.iznode[il] - params.nctz - 1;
		// covariance lookup table is large so we need to make sure
		// the point is within the simulation grid.
		if (i < 1 || j < 1 || k < 1) { continue; }
		if (i > params.nx || j > params.ny || k > params.nz) { continue ; }
		ind = i + (j - 1)*params.nx + (k-1)*params.nxy;

		if (params.sim[ind] > params.UNEST)
		{
			// check the number of data already taken from this octant:
			if (params.noct > 0)
			{
				idx = ix - i;
				idy = iy - j;
				idz = iz - k;
				if (idz > 0) { iq = 4; }
				if (idx <= 0 && idy > 0) { iq = 1; }
				if (idx > 0 && idy >= 0) { iq = 2; }
				if (idx < 0 && idy <= 0) { iq = 3; }
				else {// params.ninoct[iq] += 1;
					iq = 8;
					if (idx <= 0 && idy >0 )  { iq = 5; }
					if (idx  > 0 && idy >= 0) { iq = 6; }
					if (idx  < 0 && idy <= 0) { iq = 7; }

				}
				params.ninoct[iq] = params.ninoct[iq] + 1;
				if (params.ninoct[iq] > params.noct) { continue; }
			}

			params.ncnode = params.ncnode + 1;
			params.icnode(params.ncnode) = il;
			params.cnodex(params.ncnode) = params.xmn + std::real(i - 1)*params.xsiz;
			params.cnodey(params.ncnode) = params.ymn + std::real(j - 1)*params.ysiz;
			params.cnodez(params.ncnode) = params.zmn + std::real(k - 1)*params.zsiz;
			params.cnodev(params.ncnode) = params.sim(ind);
			//params.nclose += 1;		 
		}
		continue;
	};
}
inline void sgsim::backtr(Param&, float simval, int nt, arma::vec& vr, arma::vec& vrg, float& Backtr) {

	/*	              
					  Back Transform Univariate Data from Normal Scores
					   *************************************************

			  This subroutine backtransforms a standard normal deviate from a
			  specified back transform table and option for the tails of the
			  distribution.Call once with "first" set to true then set to false
			 unless one of the options for the tail changes.

			  INPUT VARIABLES :

				simval             normal score value to be back transformed
				nt               number of values in the back transform table
				vr(nt)           original data values that were transformed
				vrg(nt)          the corresponding transformed values
				zmin, zmax        limits possibly used for linear or power model
				ltail            option to handle values less than vrg(1) :
				ltpar            parameter required for option ltail
				utail            option to handle values greater than vrg(nt) :
				utpar            parameter required for option utail

			 Parameters			
			Value in the lower tail ? 1 = linear, 2 = power, (3 and 4 are invalid) :  

			*/

	float  cdflo, cdfhi, cdfbt, cpow, lambda;
	int j;



	// values inthe lower tail
	if (simval <= vrg[0])
	{
		Backtr = vr[0];
		cdflo = funcs::gcum(vrg[0]);
		cdfbt = funcs::gcum(simval);
		if (params.ltail == 1)
		{
			Backtr = funcs::powint(0.0, cdflo, params.zmin,vr[0],cdfbt, 1.0);
		}
		else if (params.ltail == 2) {
			cpow = 1.0 / params.ltpar;
			Backtr = funcs::powint(0.0, cdflo, params.zmin, vr[0], cdfbt, cpow);
		}
	}

	// Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
	else if (simval >= vrg[nt])
	{
		Backtr = vr(nt);
		cdfhi = funcs::gcum(vrg[nt]);
		cdfbt = funcs::gcum(simval);
				if (params.utail == 1)
				{
					Backtr = funcs::powint(cdfhi,1.0, vr[nt], params.zmax, cdfbt, 1.0);
				}
				else if (params.utail == 2)
				{
					cpow = 1.0 / params.utpar;
				  Backtr = funcs::powint(cdfhi, 1.0, params.zmax, vr[nt] , cdfbt, cpow);
				}
				else if (params.utail == 4)
				{
					lambda = float(((pow(vr[nt], params.utpar))*(1.0 - funcs::gcum(vrg[nt]))));
					Backtr = float((pow((lambda /
						     (1.0 - funcs::gcum(simval))), (1.0 / params.utpar))));
				}

	}
	else 
	{
			//Value within the transformation table:
				j = funcs::locate_c(vrg, nt, simval);
				j = funcs::max(std::min((nt - 1), j), 1);
		   Backtr = funcs::powint(params.vrg(j), params.vrg(j + 1), vr(j), vr(j + 1), simval, 1.0);	
	}
	}
void sgsim::max_covariance(Param&, int ivarg) {
//  Calculate the maximum covariance value (used for zero distances and for power model covariance)	*/
//  ivarg			variogram number(set to 0 unless doing cokriging or indicator kriging) .. note c++ 1 is 0
//   cmax            maximum covariance
	int istart;
	params.cmax = params.c0[ivarg];
	for (int i = 0; i < params.nst; i++) {		
		if (params.it[i] == 4) {
			params.cmax += params.PMX;
		}
		else {
		//params.cmax += params.cc(params.ist);   // change later
		params.cmax += params.cc[i];   // change later
		}
		
	}
};
void sgsim::simulate(Param&) {
	/*
	Note : simulation indexing starts at 1 to conform to GSLIB
		   All Data starts at 0
	
	*/
	int id, id2, nd, index,test,test2, ind, cdp, ind2,idx2;
	float xx, yy, zz, ierr;
	float  p, xp, mask, temp, tmp_backtr;
	float  TINY = 0.0001;
	float simval, ne, av, ss, UNEST, EPSLON;
	int nx, nxy, nxyz, ix, iy, iz, nz;
	int ivarg; // variogram number (set to 0 unless doing cokriging or indicator kriging)
	arma::ivec new_ind, order;
	arma::mat tmpsynth, tmpwells, convmatrix;
	arma::vec trc, RC, corrL,tmpAI,synth;
	int maxcorr;
	int indexx, tmp_ind;


	   nz = params.nz ;
	   nx = params.nx ;
	  nxy = params.nxy;
	 nxyz = params.nxyz;
	   nd = params.nd  ;        // no of data .... change to ntr
	UNEST = params.UNEST;


	new_ind.zeros(params.nz + 1);  // juts to follow fortran convention for now 
	tmpsynth.zeros(nz, params.nsynth);
	tmpwells.zeros(nz, params.nsynth);
	RC.zeros(nz);
	corrL.zeros(params.nsynth);
	EPSLON = params.EPSLON;


	std::cout << "indexing starts from 1 for sim & nodes intentionally" << std::endl;
	std::cout << "I will fix later... for now , to get the code running without bugs" << std::endl;

	ivarg = 0;
	std::cout << "Preparing for Simulation" << std::endl;

	sgsim::preprocess(params);
	sgsim::read_data(params);



	// rotation matrix
	sgsim::max_covariance(params, ivarg);
	for (int i = 0; i < params.nst; i++)
	{
		sgsim::set_rot(params, i);
	}
	params.isrot = params.MAXNST;
	sgsim::set_rot(params, (params.isrot));
	// This is neccesary such that it'll be accessible by isrot , ind 4 cpp

	if (params.sstrat == 0) 
	{
		std::cout << "Setting up super block search..." << std::endl;
	}
	//sgsim::create_searcher(params);

	// Set up the covariance table and the spiral search:
	std::cout << "Setting up Covariance Lookup Table and Spiral Search..." << std::endl;
	sgsim::ctable(params);

	for (int isim = 1; isim <= params.nsim; isim++) 
	{
		// Read in the secondary data distribution for this realization:

	   // Work out a random path for this realization:		
		sgsim::genrandpath(params, nxy,params.order);



		//  The multiple grid search works with multiples of 4 (yes, that is somewhat arbitrary):
		//  deal with multigrid later


					// Initialize the simulation:
		for (int ind = 1; ind <= nxyz; ind++)
		{
			params.sim[ind] = params.UNEST;
		}

		std::cout << "working on realization = " << isim << std::endl;

		// Assign the data to the closest grid node

		for (int id = 1; id <= params.nd; id++)
		{
		   sgsim::getindex(params, id-1, ix, iy, iz);  // note id-1
	       ind = ix + (iy - 1)*nx + (iz - 1)*nxy    ;
			xx = params.xmn + std::real(ix - 1)*params.xsiz;
			yy = params.ymn + std::real(iy - 1)*params.ysiz;
			zz = params.zmn + std::real(iz - 1)*params.zsiz;
		  test = fabs(xx - params.x(id - 1)) +     // note id-1
			 	 fabs(yy - params.y(id - 1)) +     // note id-1
				 fabs(zz - params.z(id - 1));      // note id-1

	
   // Assign this data to the node (unless there is a closer data):

			if (params.sstrat == 1)
			{
 				if (params.sim[ind] >= 0.0)  // if true.. there node is assigned
				{
					  id2 = int(params.sim[ind] + 0.5);
					test2 = fabs(xx - params.x(id2-1)) +   // note id-1
						    fabs(yy - params.y(id2-1)) +   // data starts at 0
						    fabs(zz - params.z(id2-1));
					if (test <= test2) // this shows that there is a closer data
					{
						params.sim[ind] = std::real(id);
					}
				}
				else
				{
					params.sim[ind] = std::real(id);
				}

			}
			// Assign a flag so that this node does not get simulated:
			if (params.sstrat == 0 && test <= TINY)
			{
				params.sim[ind] = 10 * params.UNEST;
			}
		}
		

		// INPUTING THE DATA .... note the index conforms to the way the data is loaded.
		for (int ind = 1 ; ind <= nxyz; ind++)
		{
			id = int(params.sim[ind] + 0.5);
			if (id > 0)
			{
				params.sim[ind] = params.vr[id - 1]; // note id-1 for vr
			}
		}
		
		
		
		
	//----------------------------------  MAIN LOOP OVER ALL THE NODES:--------------------------------
	//	GP::convmtx2(params.wavlet, params.nz, convmatrix); // outputs convoution matrix

		for (int i = 0; i < nxy; i++)
		{
			cout << "working on nxy  = " << i << "  out of =" << nxy << endl;
			        corrL.zeros(params.nsynth); // reseting correlations
					for (int ii = 0; ii < params.nsynth; ii++) 
					{
						/* 
						No of synthetics generated
						CorrL must start at 0
						new_index ... 1 to nz
						*/
					
						cdp = params.order[i];
						new_ind.zeros(nz+1);
						
						for (int iii = 1; iii <= params.nz; iii++) {
							// looping over a well
						
							         index = cdp + (iii - 1)* nxy ;
							  new_ind[iii] = index                ;

							//  Figure out the location of this point and make sure it has not been assigned a value already :
							
							if (params.sim[index] > params.UNEST + params.EPSLON ||
								params.sim[index] < params.UNEST * 2)
							{
								continue;
							}

							// In this case , we are working backwards 
							// Inputing index & figuring out xx, yy , zz 
							iz = int((index - 1) / nxy) + 1                 ;
							iy = int((index - (iz - 1)*nxy - 1) / nx) + 1   ;
							ix = index - (iz - 1)*nxy - (iy - 1)*nx         ;
							xx = params.xmn + std::real(ix - 1)*params.xsiz ;    // xcoor
							yy = params.ymn + std::real(iy - 1)* params.ysiz;
							zz = params.zmn + std::real(iz - 1)*params.zsiz ;
	/*
	    Now, we'll simulate the point ix,iy,iz.
		1. First, get the close data and make sure that there are enough to actually simulate a value
		2. we'll only keep the closest "ndmax" data, and look for previously simulated grid nodes:
		*/

							if (params.sstrat == 0)
							{
								std::cout << "Super block search is not implimented" << std::endl;

								if (params.nclose < params.ndmin)
								{
									continue;
								}
								if (params.nclose > params.ndmax)
								{
									params.nclose = params.ndmax;
								}
							}
							else
							{
								params.nclose = 0;
							}

							/*	
							spiral search on the covariance lookup table
							calculate the conditional mean and standard deviation
							this will be done with kriging if there are data,
							otherwise, the glocal mean and standard deviation will be used.
							*/

							sgsim::srchnd(params, ix, iy, iz);  // output ixnode, iynode, iznode

							if (params.ktype == 2)
							{
								params.gmean = params.lvm[index-1];
							}
							else
							{
								params.gmean = 0;
							}


							if (params.nclose + params.ncnode < 1)
							{
								params.cmean = params.gmean;
								params.cstdev = 1;
							}
							else
							{
								/*
								Perform kriging
								Note that if there are fewer than four data then simple kriging is prefered so that the variance of the
								realization does not become artificially inflated :
								*/
								
								params.lktype = params.ktype;
								if (params.ktype == 1 && (params.nclose + params.ncnode) < 4)
								{
									params.lktype = 0;
								}
								sgsim::krigging(params, ix, iy, iz, xx, yy, zz);
							}

							 p = sgsim::acorni(params);
							xp = funcs::gauinv(p, ierr)   ;
							//params.sim[index] = xp *params.cstdev + params.cmean;
							params.sim[index] = xp * params.cstdev + params.cmean;
							
						} // end simulating a well


						  // do we need to re-assign data for al points........


						if (params.sstrat == 0)
						{
							for (int idx = 1; idx <= params.nz; idx++)
							{
								idx2 = new_ind[idx] -1 ; // data starts at 0
								sgsim::getindex(params, idx2, ix, iy, iz);
								ind = ix + (iy - 1)*nx + (iz - 1)*nxy;
								xx = params.xmn + std::real(ix - 1)*params.xsiz;
								yy = params.ymn + std::real(iy - 1)*params.ysiz;
								zz = params.zmn + std::real(iz - 1)*params.zsiz;
								test = fabs(xx - params.x[idx]) +
									   fabs(yy - params.y[idx]) +
									   fabs(zz - params.z[idx]);
								if (test <= TINY)
								{
									params.sim[idx2+1] = params.vr[idx2];// data starts at 0
								}

							}
						}



						//  back-transform each value and write results: 

							/*
							if (params.itrans == 1) {
							mask = params.sim > UNEST + EPSLON;
							temp = params.sim[mask];
							params.sim[mask] = params.zmax;
							*/

							
							int nt = params.nd - 1;
							for (int ind = 1; ind <= params.nz; ind++)		
							{
								ne = 0.0;
								av = 0.0;
								ss = 0.0;

								simval = params.sim(new_ind(ind));  // need debugging
								
								if (simval > -9.0 && simval < 9.0)
								{
									ne = ne + 1;
									av = av + simval;
									ss = ss + simval * simval;
								}
								if (params.itrans == 1 && simval >(params.UNEST + params.EPSLON))
								{
						    	  sgsim::backtr(params, simval, nt, params.vrtr, params.vrgtr, tmp_backtr);  // passing out Backtr				
								//  cout << "back-transformed data  = " << tmp_backtr << endl;
								//  cout << "simulated data  = " << simval << endl;
									if (tmp_backtr < params.zmin)
									{ 
										tmp_backtr = params.zmin;
									}
									if (tmp_backtr > params.zmax)
									{ 
										tmp_backtr = params.zmax;
									}
									params.sim(new_ind(ind)) = tmp_backtr;
								}

								/*
								cout << "new_ind(ind)) = " << new_ind(ind) << endl;
								cout << "simulation = " << simval << endl;
								cout << "backtransformed = " << tmp_backtr << endl;
								*/
							}

						
						//	tmpwells.col(ii) = params.sim(arma::span(new_ind[1], new_ind[nz]));
							GP::indV2M(tmpwells,ii, 0, nz, params.sim, new_ind,1,nz);
						
						// add debugging
						
		//-----------------------------------calculate synthetics-----------------------------			
						tmpAI = tmpwells.col(ii);
						GP::calc_RC(tmpAI, RC);  	
						synth = arma::conv(RC, params.wavlet, "same"); // this equals the size of the second argunment
						tmpsynth.col(ii) = arma::conv(RC, params.wavlet, "same");
						GP::indV2M(tmpwells, ii, 0, nz, params.sim, new_ind, 1, nz);
						trc.zeros(nz);
						//trc = params.seis_grid(arma::span(new_ind[1], new_ind[nz]));
						GP::indV2V(trc, 0, nz-1, params.seis_grid, new_ind, 1, nz);
						corrL[ii] = GP::corr(synth, trc);  // write a correlation func

		//---------------------------------------DEBUGGING----------------------------------------------
						
						/*
						float*pwav = GP::arma2array(params.wavlet);
						float*pAI = GP::arma2array(tmpAI);
						float*pRC = GP::arma2array(RC);
						float*psynth = GP::arma2array(synth);
						float*ptrc   = GP::arma2array(trc);
						*/		
						
		//---------------------------------------DEBUGGING----------------------------------------------
						for (int ir = 1; ir <= nz; ir++)  //  reseting
						{
							params.sim(new_ind[ir]) = params.UNEST;  
						}

					} // end for nsynth
					arma::vec test_inv;
					test_inv.zeros(nz);

					maxcorr = arma::index_max(corrL); // starts from 0
					for (int ir = 1; ir <= nz; ir++)   // sim & new_ind starts at 1 but tmpwells starts at 0
					{
			   		params.sim(new_ind[ir]) = tmpwells(ir-1, maxcorr);  // replacing with the max correlation co-efficint 
					test_inv[ir-1] = tmpwells(ir - 1, maxcorr);
					}
				//	float*pinv = GP::arma2array(test_inv);

					tmpsynth.zeros(nz, params.nsynth); // resetting to zeros ;
				} // end for nxy

		

				params.sim.save("sim.bin");  // you need to add using namespace arma to the main project.
		
				

			/*
			params.simu.save("simu.txt"); ;
			params.sim.save("sim.txt");
			cout << "sim = " << endl;
			params.simu.print();
			*/
	}
}
inline void sgsim::krigging(Param& params, float ix, float iy, float iz, float xx, float yy, float zz)
	{
		/*
		-----------------------------------------------------------------------

		Builds and Solves the SK or OK Kriging System
		*********************************************

		INPUT VARIABLES:

		ix,iy,iz        index of the point currently being simulated
		xx,yy,zz        location of the point currently being simulated


		OUTPUT VARIABLES:

		cmean           kriged estimate
		cstdev          kriged standard deviation


		EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
		-----------------------------------------------------------------------
		*/


		bool first;
		int na, neq;
		float cov;
		int in, ind, nx, nxy, nxyz, index, irot, ising;
		int x1, y1, z1, x2, y2, z2;
		int  ix1, iy1, iz1, ix2, iy2, iz2;
		int ii, jj, kk, nclose, ncnode;
		arma::vec point1, point2;
		
		arma::vec a, r, rr,s, vra;
		// pre-allocating 
		a.zeros(params.MAXKR2);
		vra.zeros(params.MAXKR1);
		r.zeros(params.MAXKR1);
		rr.zeros(params.MAXKR1);
		s.zeros(params.MAXKR1);
		point1.zeros(2); point2.zeros(2);

	




		nx = params.nx;
		nxy = params.nxy;
		nxyz = params.nxyz;
		nclose = params.nclose;
		ncnode = params.ncnode;

		irot = 0;

		first = false;
		na = params.nclose + params.ncnode;
		//L:   continue
		if (params.lktype == 0)  	{neq = na;}
		if (params.lktype == 1)		{neq = na + 1;	}
		if (params.lktype == 2)		{neq = na;		}
		if (params.lktype == 3)		{neq = na + 2;	}
		if (params.lktype == 4) 	{neq = na + 1;	}
		if (params.lktype >= 3) 	
		{
			ind = ix + (iy - 1)*nx + (iz - 1)*nxy;
		}

	

		// Set up kriging matrices:
		in = 0;
		for ( int j = 1; j <= na; ++j) 
		{

			// Sort out the actual location of point "j"
			if (j <= params.nclose) 
			{
				index = int(params.close[j]);
				x2 = params.x[index-1]; // check d indexing of x -  it starts from 0
				y2 = params.y[index-1];
				z2 = params.z[index-1];
				params.vra(j) = params.vr[index-1]; // check indexing
						params.vrea(j-1) = params.sec[index - 1];
						if (params.lktype == 2)
						{
							params.vra[j-1] = params.vra[j-1] - params.vrea(j-1);
						}
			}

			// It is a previously simulated node (keep index for table look-up):
			else
			{
				index = j - nclose;
				x1 = params.cnodex(index);
				y1 = params.cnodey(index);
				z1 = params.cnodez(index);
				params.vra(j-1) = params.cnodev(index-1);
				ind = params.icnode(index);
				ix1 = ix + (int(params.ixnode(ind)) - params.nctx - 1);
				iy1 = iy + (int(params.iynode(ind)) - params.ncty - 1);
				iz1 = iz + (int(params.iznode(ind)) - params.nctz - 1);
				index = ix1 + (iy1 - 1)*nx + (iz1 - 1)*nxy;
				params.vrea(j-1) = params.lvm[index -1];
						if (params.lktype == 2)
						{
							params.vra(j-1) = params.vra(j-1) - params.vrea(j-1);
						}
			
			}

			// Sort out the actual location of point "i"
			for (int i = 1; i <= j; ++i) {
							if (i <= nclose)
							{
								index = int(params.close(i));
								x2 = params.x[index-1];
								y2 = params.y[index-1];
								z2 = params.z[index-1];
							}
				// It is a previously simulated node (keep index for table look-up):
						else {
							index = i - nclose;
							x2 = params.cnodex(index);
							y2 = params.cnodey(index);
							z2 = params.cnodez(index);
							ind = params.icnode(index);
							ix2 = ix + (int(params.ixnode(ind)) - params.nctx - 1);
							iy2 = iy + (int(params.iynode(ind)) - params.ncty - 1);
							iz2 = iz + (int(params.iznode(ind)) - params.nctz - 1);
							
						}

				// Now, get the covariance value:
				in = in + 1;

				// Decide whether or not to use the covariance look-up table:
						if (j <= nclose || i <= nclose)
						{
							point1[0] = x1; point1[1] = y1; point1[2] = z1;
							point2[0] = x2; point2[1] = y2; point2[2] = z2;
							cov = funcs::cova3(point1, point2, params, irot, irot);
						  ta[in] = double(cov);
						}

				// Try to use the covariance look-up (if the distance is in range):
					else {
						ii = params.nctx + 1 + (ix1 - ix2);
						jj = params.ncty + 1 + (iy1 - iy2);
						kk = params.nctz + 1 + (iz1 - iz2);
						if (ii < 1 || ii > params.MAXCTX ||
							jj < 1 || jj > params.MAXCTY ||
							kk < 1 || kk > params.MAXCTZ) {
							point1[0] = x1; point1[1] = y1; point1[2] = z1;
							point2[0] = x2; point2[1] = y2; point2[2] = z2;
							cov = funcs::cova3(point1, point2, params, irot, irot);
						}
						else
						{
							cov = params.covtab(ii, jj, kk);
						}
					  	  a[in] = double(cov);
				}
			}

			// Get the RHS value (possibly with covariance look-up table):
			if (j <= params.nclose) {
				point1[0] = xx; point1[1] = yy; point1[2] = zz;
				point2[0] = x1; point2[1] = y1; point2[2] = z1;
				cov = funcs::cova3(point1, point2, params, irot, irot);	
				 r[j] = double(cov);
				//rr[j] = r[j];  ?
			}

			// Try to use the covariance look-up (if the distance is in range):
			else 
			{
				ii = params.nctx + 1 + (ix - ix1);
				jj = params.ncty + 1 + (iy - iy1);
				kk = params.nctz + 1 + (iz - iz1);

				if (ii < 1 || ii > params.MAXCTX ||
					jj < 1 || jj > params.MAXCTY ||
					kk < 1 || kk > params.MAXCTZ)
				{
					point1[0] = xx; point1[1] = yy; point1[2] = zz;
					point2[0] = x1; point2[1] = y1; point2[2] = z1;
					cov = funcs::cova3(point1, point2, params, irot, irot);
				}
				else
				{
					cov = params.covtab(ii, jj, kk);
				}
				r[j] = double(cov);
			}
			rr[j] = r[j];
	}

		// Addition of OK constraint:

	if (params.lktype == 1 || params.lktype == 3)
	{
		for (int i = 1; i <= na; i++) {
			in = in + 1;
			a[in] = 1.0;
		}
		in = in + 1;
		a[in] = 0.0;
		r[na + 1] = 1.0;
		rr[na + 1] = 1.0;
	}

	// Solve the Kriging System :

	if (neq == 1 && params.lktype != 3)
	{
		s[1] = r[1] / a[1];
		ising = 0;
	}
	else
	{
		funcs::ksol(1, neq, 1, a, r, s, ising);

	}

	//Write a warning if the matrix is singular

				if (ising != 0) 
				{
					std::cout << "SGSIM is singular for node  = " << ix << iy << iz << std::endl;
				}



				// Compute the estimate and kriging variance.Recall that kriging type
				//  0 = Simple Kriging :
				//  1 = Ordinary Kriging :
				//  2 = Locally Varying Mean :
				//  3 = External Drift :
				//  4 = Collocated Cosimulation :

				double sumwts      ;
				params.cmean  = 0.0;
				params.cstdev = params.cmax;
				sumwts = 0.0;
				for (int i = 1; i <= na; i++)
				{
					params.cmean = params.cmean + real(s(i))*params.vra(i);
					params.cstdev = params.cstdev - real(s(i)*rr(i));
					sumwts = sumwts + std::real(s(i));
				}
				if (params.lktype == 1)
				{
					params.cstdev = params.cstdev - std::real(s(na + 1));
				}
				if (params.lktype == 2) 
				{
					params.cmean = params.cmean + params.gmean;
				}
				if (params.lktype == 4) 
				{
					         ind = ix + (iy - 1)*nx + (iz - 1)*nxy;
					params.cmean = params.cmean + real(s(na + 1))*params.lvm(ind-1);
				   params.cstdev = params.cstdev - real(s(na + 1) *rr(na + 1))   ;
				}



				// Error message if negative variance:

				if (params.cstdev < 0.0) 
				{
					std::cout << "ERROR: Negative Variance:  = " << params.cstdev << std::endl;
					params.cstdev = 0;
				}

				params.cstdev = std::sqrt(funcs::max(params.cstdev, 0));
	}



