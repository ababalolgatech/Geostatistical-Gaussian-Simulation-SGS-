#include "ConFig.h"
#include"SGS_Param.h"



void read_param(Param& param, std::string fp) {
	ConfigFile cfg(fp);
	bool exist;

	std::string datafl = cfg.getValueOfKey<std::string>("datafl");
	exist = cfg.keyExists("datafl");
	if (exist == 0) { std::cout << " Check Variable Name in parameter file" << datafl << "\n"; }
	else { param.datafl = datafl; }

	param.wavefl = cfg.getValueOfKey<std::string>("wavefl");
	param.seisfl = cfg.getValueOfKey<std::string>("seisfl");
	param.ixl = cfg.getValueOfKey<int>("ixl");
	param.iyl = cfg.getValueOfKey<int>("iyl");
	param.izl = cfg.getValueOfKey<int>("izl");
	param.ivrl = cfg.getValueOfKey<int>("ivrl");
	param.icolsec = cfg.getValueOfKey<int>("icolsec");
	param.iwt = cfg.getValueOfKey<int>("iwt");
	param.tmin = cfg.getValueOfKey<float>("tmin");
	param.tmax = cfg.getValueOfKey<float>("tmax");
	param.itrans = cfg.getValueOfKey<bool>("itrans");
	param.transfl = cfg.getValueOfKey<std::string>("transfl");
	param.ismooth = cfg.getValueOfKey<bool>("ismooth");
	param.seed    = cfg.getValueOfKey<double>("seed");
	param.smthfl = cfg.getValueOfKey<std::string>("smthfl");
	param.isvr = cfg.getValueOfKey<int>("isvr");
	param.iswt = cfg.getValueOfKey<int>("iswt");
	// allowable data values used for backtransform
	param.zmin = cfg.getValueOfKey<float>("zmin");
	param.zmax = cfg.getValueOfKey<float>("zmax");
	//lower and upper tail model specification for backtransform
	param.ltail = cfg.getValueOfKey<float>("ltail");
	param.ltpar = cfg.getValueOfKey<float>("ltpar");
	param.utail = cfg.getValueOfKey<float>("utail");
	// debug and output data file
	param.idbg = cfg.getValueOfKey<int>("idbg");
	param.dbgfl = cfg.getValueOfKey<std::string>("dbgfl");
	param.outfl = cfg.getValueOfKey<std::string>("outfl");
	param.nsim = cfg.getValueOfKey<int>("nsim"); // nsim
	param.nsynth = cfg.getValueOfKey<int>("nsynth"); // nsim
	// Grid definition
	param.nx = cfg.getValueOfKey<int>("nx");
	param.ny = cfg.getValueOfKey<int>("ny");
	param.nz = cfg.getValueOfKey<int>("nz");
	param.xmn = cfg.getValueOfKey<int>("xmn");
	param.ymn = cfg.getValueOfKey<int>("ymn");
	param.zmn = cfg.getValueOfKey<int>("zmn");
	param.xsiz = cfg.getValueOfKey<int>("xsiz");
	param.ysiz = cfg.getValueOfKey<int>("ysiz");
	param.zsiz = cfg.getValueOfKey<int>("zsiz");

	param.ndmin = cfg.getValueOfKey<int>("ndmin"); //% minimum data points used in kriging
	param.ndmax = cfg.getValueOfKey<int>("ndmax"); //% maximum data points used in kriging
	if (param.sstrat == 1) { param.ndmax = 0; }; // line 234 GSLIB
	param.nodmax = cfg.getValueOfKey<int>("nodmax"); //  previously simulated nodes to use
	param.sstrat = cfg.getValueOfKey<int>("sstrat");  //% search strategy

	if (param.sstrat == 10) {
		param.ndmax = 0;
	}
	param.multgrid = cfg.getValueOfKey<int>("sstrat");
	// search radii
	param.radius_hmax = cfg.getValueOfKey<float>("radius_hmax");
	param.radius_hmin = cfg.getValueOfKey<float>("radius_hmin");
	param.radius_vert = cfg.getValueOfKey<float>("radius_vert");
	// search anisotropy angles
	param.sang1 = cfg.getValueOfKey<float>("sang1");
	param.sang2 = cfg.getValueOfKey<float>("sang2");
	param.sang3 = cfg.getValueOfKey<float>("sang3");
	//size of covariance lookup table size
	param.mxctx = cfg.getValueOfKey<int>("mxctx");
	param.mxcty = cfg.getValueOfKey<int>("mxcty");
	param.mxctz = cfg.getValueOfKey<int>("mxctz");
	param.ktype = cfg.getValueOfKey<int>("ktype");

	param.trans = true;
	if (param.ktype < 0) {
		param.trans = false;
		param.ktype = std::abs(param.ktype);
	}

	// self.skmean = param['skmean']
	param.rho = cfg.getValueOfKey<float>("rho"); //  correlation coefficient for COCOK
	param.varred = cfg.getValueOfKey<float>("varred"); //variance reduction factor for COCOK
	param.secfl = cfg.getValueOfKey<std::string>("secfl");
	param.icollvm = cfg.getValueOfKey<int>("icollvm");

	// variography definition
	// this is a bad code
	param.nst = cfg.getValueOfKey<int>("nst");

	float c01, it1, cc1, ang11, ang21, ang31, aa_hmax1, aa_hmin1, aa_vert1;
	float c02, it2, cc2, ang12, ang22, ang32, aa_hmax2, aa_hmin2, aa_vert2;
	float c03, it3, cc3, ang13, ang23, ang33, aa_hmax3, aa_hmin3, aa_vert3;
	c01 = cfg.getValueOfKey<float>("c01");
	it1 = cfg.getValueOfKey<int>("it1");
	cc1 = cfg.getValueOfKey<float>("cc1");
	ang11 = cfg.getValueOfKey<float>("ang11");
	ang21 = cfg.getValueOfKey<float>("ang21");
	ang31 = cfg.getValueOfKey<float>("ang31");
	aa_hmax1 = cfg.getValueOfKey<float>("aa_hmax1");
	aa_hmin1 = cfg.getValueOfKey<float>("aa_hmin1");
	aa_vert1 = cfg.getValueOfKey<float>("aa_vert1");

	c02 = cfg.getValueOfKey<float>("c02");
	it2 = cfg.getValueOfKey<int>("it2");
	cc2 = cfg.getValueOfKey<float>("cc2");
	ang12 = cfg.getValueOfKey<float>("ang12");
	ang22 = cfg.getValueOfKey<float>("ang22");
	ang32 = cfg.getValueOfKey<float>("ang32");
	aa_hmax2 = cfg.getValueOfKey<float>("aa_hmax2");
	aa_hmin2 = cfg.getValueOfKey<float>("aa_hmin2");
	aa_vert2 = cfg.getValueOfKey<float>("aa_vert2");

	c03 = cfg.getValueOfKey<float>("c03");
	it3 = cfg.getValueOfKey<int>("it3");
	cc3 = cfg.getValueOfKey<float>("cc3");
	ang13 = cfg.getValueOfKey<float>("ang13");
	ang23 = cfg.getValueOfKey<float>("ang23");
	ang33 = cfg.getValueOfKey<float>("ang33");
	aa_hmax3 = cfg.getValueOfKey<float>("aa_hmax3");
	aa_hmin3 = cfg.getValueOfKey<float>("aa_hmin3");
	aa_vert3 = cfg.getValueOfKey<float>("aa_vert3");


	// Pre-allocating the vectors
	param.c0.set_size(param.nst);
	param.it.set_size(param.nst);
	param.cc.set_size(param.nst);
	param.ang1.set_size(param.nst);
	param.ang2.set_size(param.nst);
	param.ang3.set_size(param.nst);
	param.aa_hmax.set_size(param.nst);
	param.aa_hmin.set_size(param.nst);
	param.aa_vert.set_size(param.nst);
	param.anis1.set_size(param.nst);
	param.anis2.set_size(param.nst);



	if (param.nst == 1) {
		param.c0[0] = c01;		param.it[0] = it1;
		param.cc[0] = cc1;		param.ang1[0] = ang11;
		param.ang1[0] = ang11;  param.ang2[0] = ang21;
		param.ang3[0] = ang31;	param.aa_hmax[0] = aa_hmax1;
		param.aa_hmin[0] = aa_hmin1; param.aa_vert[0] = aa_vert1;
	}



	else if (param.nst == 2) {
		param.c0[0] = c01;		param.it[0] = it1;
		param.cc[0] = cc1;		param.ang1[0] = ang11;
		param.ang1[0] = ang11;  param.ang2[0] = ang21;
		param.ang3[0] = ang31;	param.aa_hmax[0] = aa_hmax1;
		param.aa_hmin[0] = aa_hmin1; param.aa_vert[0] = aa_vert1;

		param.c0[1] = c02;		param.it[1] = it2;
		param.cc[1] = cc2;		param.ang1[1] = ang12;
		param.ang1[1] = ang12;  param.ang2[1] = ang22;
		param.ang3[1] = ang32;	param.aa_hmax[1] = aa_hmax2;
		param.aa_hmin[0] = aa_hmin2; param.aa_vert[1] = aa_vert2;
	}
	else if (param.nst == 3)
	{
		param.c0[0] = c01;		param.it[0] = it1;
		param.cc[0] = cc1;		param.ang1[0] = ang11;
		param.ang1[0] = ang11;  param.ang2[0] = ang21;
		param.ang3[0] = ang31;	param.aa_hmax[0] = aa_hmax1;
		param.aa_hmin[0] = aa_hmin1; param.aa_vert[0] = aa_vert1;

		param.c0[1] = c02;		param.it[1] = it2;
		param.cc[1] = cc2;		param.ang1[1] = ang12;
		param.ang1[1] = ang12;  param.ang2[1] = ang22;
		param.ang3[1] = ang32;	param.aa_hmax[1] = aa_hmax2;
		param.aa_hmin[1] = aa_hmin2; param.aa_vert[1] = aa_vert2;

		param.c0[2] = c03;		param.it[2] = it3;
		param.cc[2] = cc3;		param.ang1[2] = ang13;
		param.ang1[2] = ang13;  param.ang2[2] = ang23;
		param.ang3[2] = ang33;	param.aa_hmax[2] = aa_hmax3;
		param.aa_hmin[2] = aa_hmin3; param.aa_vert[2] = aa_vert3;
	}


	
	// send to sgs::preprocess eventually 
	param.nxy = param.nx * param.ny;
	param.nxyz = param.nx * param.ny* param.nz;
	param.MAXSAM = param.ndmax + 1;
	param.KORDEI = 12;
	param.MAXINT = 2 ^ 30;
	param.MAXOP1 = param.KORDEI + 1;
	param.TINY = 1.0e-10;
	param.tmin = -1.0e21;
	param.tmax = 1.0e21;
	param.EPSLON = 1.0e-20;



	for (int i = 0; i < param.nst; i++) {
		param.sill = param.c0[i] + param.cc[i];
		param.anis1[i] = param.aa_hmin[i] / param.aa_hmax[i];
		param.anis2[i] = param.aa_vert[i] / param.aa_hmax[i];
	}
	// aa = aa_hmax
	// aa1 = aa_hmin
	// aa2 = aa_vert

	if (param.radius_hmax < param.EPSLON) {
		std::cout << "Radius must be greater than zero" << std::endl;
		exit(0);
	}

	param.radsq = param.radius_hmax*param.radius_hmax;
	param.sanis1 = param.radius_hmin / param.radius_hmax;
	param.sanis2 = param.radius_vert / param.radius_hmax;

	/*
	std::cout << "radius_hmax = " << param.radius_hmax << std::endl;
	std::cout << "radius_hmin = " << param.radius_hmin << std::endl;
	std::cout << "sanis1 = " << param.sanis1 << std::endl;
	std::cout << "sanis2 = " << param.sanis2 << std::endl;

	*/
	


};




