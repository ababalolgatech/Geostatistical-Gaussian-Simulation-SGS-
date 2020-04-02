
// I use a nested class to create a parameter file

class Param {
public:
	// -------------------VARIOGRAM-------------------
	int nst;   // The number of semivariogram structures
	arma::vec  it;   //  it:Type of nested structures (1-sph,2-exp, 3 - gau, 3-pow) */
	arma::vec c0;    //	c0 : the isotropic nugget constant.
	arma::vec cc;    //  multiplicative factor for each nested structure       
	arma::vec ang1, ang2, ang3;  // ang1`, `ang2`, `ang3`, the angles defining the geometric anisotropy;
	arma::vec aa_hmax;  // the maximum horizontal range;
	arma::vec aa_hmin;  // the minimum horizontal range
	arma::vec aa_vert;  // the vertical range.
	arma::vec  anis1; // Anisotropy for the dip angle
	arma::vec  anis2; // Anisotropy for the plunge angle
	float utail, utpar, varred;
	float sill;
	arma::cube covtab;
	float cmax; // maximum covariance
	float ltail, ltpar;
	// NOTE
	// aa = aa_hmax
	// aa1 = aa_hmin
	// aa2 = aa_vert

	//------------------- SEARCH----------------------
	float radius_hmax, radius_hmin, radius_vert;	// maximum search radius 
	float sang1, sang2, sang3; 	// the angle parameters that describe the orientation of the search ellipsoid.
	float sanis1, sanis2, radsq;
	float noct;   // Maximum number per octant if an octant search is desired (if <= 0, then no octant search)
	float nclose;
	float ndmin; // Minimum number of data required before sim
	float ndmax; // Maximum number of samples for simulation
	float nlookup;
	arma::cube rotmat;
	arma::vec cnodex, cnodey, cnodez, icnode, cnodev;
	arma::vec ixnode, iynode, iznode, ninoct;
	int nctx, ncty, nctz, nlooku, ncnode;

	// -------------------CONSTANTS-------------------
	float  EPSLON, UNEST, MXDT, PMX, MAXNST;
	float MAXSAM, maxsbx, maxsby, maxsbz, MAXDT, MAXSB;
	float MAXXYZ, MAXY, MAXZ, MAXX, MXYZ, MAXKR2;

	// -------------------DATA-------------------
	/*
	iwt - declustering weights
	isecvr - for external drift if used
	zmin& zmax allowable data values used for backtransform
	*/
	int icollvm, icolsec, iswt, isvr, ivrl, iwt;
	int isecvr, idbg;
	int icolvr, nd;
	int ivr, ixl, iyl, izl;
	std::string tmpfl, datafl, smthfl, dbgfl, secfl;
	std::string wavefl, seisfl;
	std::string transfl, outfl;
	int icolx, icoly, icolz;
	bool ismooth;
	float seed, tmin, tmax, zmin, zmax, rho;
	int sstrat;
	arma::mat var;
	arma::vec sec, wt, vrtr, vrgtr, vrg, vrgs, vra, vr, lvm;
	arma::vec tmp, vrea, r, s, a, nisb, ixsbtosr, iysbtosr, izsbtosr;
	arma::vec ixv;  // for arcorni 

	// SIMULATIONS
	//	noct : maximum number to retain from an octant	
	bool trans, multgrid;
	int mxctx, mxcty, mxctz;
	int nxctx, nxcty, nxctz, nmult;
	int nodmax;
	int ktype, lktype;
	float cmean, gmean, cstdev, Backtr;
	arma::mat left, right; /// for krig matrix
	arma::vec  close;
	arma::vec sim, order;
	arma::mat simu; // the backtransformed
	arma::uvec sort_index;
	// CONSTANTS
	int KORDEI, MAXINT, MAXOP1;
	float TINY;

	// -------------------GRID-------------------
	int  nx, ny, nz;
	int xsiz, ysiz, zsiz;
	int xmn, ymn, zmn;
	int nxy, nxyz;
	int maxbx, maxby, maxbz;
	int xloc, yloc, zloc;
	arma::vec x, y, z;
	arma::vec xx, yy, zz;
	float MAXNOD, MAXCXY, MAXCTX, MAXCTY, MAXCTZ;
	float MAXKR1, MAXSBX, MAXSBY, MAXSBZ;
	float MAXCUT, MXCUT, MAXROT,isrot;
	// GENERAL
	int nsim, itrans;
	int ntr, idb, lin;
	int lout, ldbg, llvm;
	int nvaril;
	int colorcorr, test;


	//------------- Geostat inversion
	arma::vec seis_grid, wavlet, trc,RC;
	arma::mat synth, simdebug,synthdebug;
	int nsynth;


};



