
	//////////////////////////
// Copyright (c) 2015-2019 Julian Adamek/ 2019 -2023 Farbod Hassani
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//////////////////////////

//////////////////////////
// main.cpp
//////////////////////////
//
// main control sequence of k-evolution N-body code based on gevolution.
// Author (k-evolution): Farbod Hassani (Université de Genève & Universitetet i Oslo)
// Author (gevolution): Julian Adamek  (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: Sept 22 2023  by Farbod Hassani
//
//////////////////////////
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>    // For std::setw, std::left, std::right
#include <cstring>    // For snprintf
#include <filesystem>
#include <set>
#include <vector>
#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)
#include "class.h"
#undef MAX			// due to macro collision this has to be done BEFORE including LATfield2 headers!
#undef MIN
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#ifdef HAVE_CLASS
#include "class_tools.hpp"
#endif
#ifdef HAVE_HICLASS
#define HAVE_HICLASS_BG HAVE_HICLASS
#include "hiclass_tools.hpp"
#endif
#include "tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"
#include "gevolution.hpp"
#include "ic_basic.hpp"
#include "ic_read.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "radiation.hpp"
#include "parser.hpp"
#include "output.hpp"
#include "hibernation.hpp"
#ifdef VELOCITY
#include "velocity.hpp"
#endif

using namespace std;
using namespace LATfield2;

int main(int argc, char **argv)
{
	#ifdef BENCHMARK
		//benchmarking variables
		double ref_time, ref2_time, cycle_start_time;
		double initialization_time;
		double run_time;
		double kgb_update_time=0; 
		double cycle_time=0;
		double projection_time = 0;
		double snapshot_output_time = 0;
		double spectra_output_time = 0;
		double lightcone_output_time = 0;
		double gravity_solver_time = 0;
		double fft_time = 0;
		int fft_count = 0;
		double update_q_time = 0;
		int update_q_count = 0;
		double moveParts_time = 0;
		int  moveParts_count =0;
		// KGB
		double a_kgb;
		double Hc;
	#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;
	int i, j, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numsteps, numspecies, done_hij;
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dx, tau, a, fourpiG, tmp, start_time;
	double maxvel[MAX_PCL_SPECIES];
	FILE * outfile;
	FILE * check_file;
	char filename[2*PARAM_MAX_LENGTH+24];
	string h5filename;
	char * settingsfile = NULL;
	char * precisionfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	double T00hom;
	double T00KGBhom;
	#ifdef HAVE_HICLASS_BG
		gsl_interp_accel * acc = gsl_interp_accel_alloc();
		gsl_spline * H_spline = NULL;
		gsl_spline * H_prime_spline = NULL;
		gsl_spline * H_prime_prime_spline = NULL;
		gsl_spline * rho_cdm_spline = NULL;
		gsl_spline * rho_b_spline = NULL;
		gsl_spline * rho_g_spline = NULL;
		gsl_spline * rho_crit_spline = NULL;
		gsl_spline * rho_ur_spline = NULL;
		gsl_spline * cs2_spline = NULL;
		gsl_spline * cs2_prime_spline = NULL;
		gsl_spline * rho_smg_spline = NULL;
		gsl_spline * rho_smg_prime_spline = NULL;
		gsl_spline * p_smg_spline = NULL;
		gsl_spline * p_smg_prime_spline = NULL;
		gsl_spline * alpha_K_spline = NULL;
		gsl_spline * alpha_K_prime_spline = NULL;
  		gsl_spline * alpha_B_spline = NULL;
  		gsl_spline * alpha_B_prime_spline = NULL;
		gsl_spline * cs2num_spline = NULL;
		gsl_spline * kin_D_spline = NULL;
		gsl_spline * lambda_2_spline = NULL;
	#endif

	#ifndef H5_DEBUG
		H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
	#endif

	for (i=1 ; i < argc ; i++ )
	{
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) 
		{
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
				break;
			case 'p':
			#if !defined(HAVE_CLASS) && !defined(HAVE_HICLASS)
			cout << "HAVE_CLASS needs to be set at compilation to use CLASS precision files" << endl;
			exit(-100);
			#endif
			precisionfile = argv[++i];
			break;
			case 'i':
			#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
			#endif
			io_size =  atoi(argv[++i]);
			break;
			case 'g':
			#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
			#endif
			io_group_size = atoi(argv[++i]);
		}
	}

	#ifndef EXTERNAL_IO
		parallel.initialize(n,m);
	#else
	if (!io_size || !io_group_size)
	{
		cout << "invalid number of I/O tasks and group sizes for I/O server (-DEXTERNAL_IO)" << endl;
		exit(-1000);
	}
		parallel.initialize(n,m,io_size,io_group_size);
		if(parallel.isIO()) ioserver.start();
		else
	{
	#endif

	COUT << COLORTEXT_CYAN << endl;
	COUT << "                                                      "<<endl;
	COUT <<COLORTEXT_LIGHT_BROWN<<"KKKKKKKKK    KKKKKKK             GGGGGGGGGGGGG     BBBBBBBBBBBBBBBBB" <<COLORTEXT_CYAN  <<endl; 
	COUT <<"K:::::::K    K:::::K          GGG::::::::::::G     B::::::::::::::::B"  <<endl;
	COUT << "K:::::::K    K:::::K        GG:::::::::::::::G     B::::::BBBBBB:::::B"  <<endl;
	COUT <<"K:::::::K   K::::::K       G:::::GGGGGGGG::::G     BB:::::B     B:::::B"<<endl;
	COUT <<"KK::::::K  K:::::KKK      G:::::G       GGGGGG     B:::::B     B:::::B"<<endl;
	COUT <<"K:::::K K:::::K        G:::::G                     B::::B     B:::::B"  <<endl;
	COUT <<"K::::::K:::::K         G:::::G                     B::::BBBBBB:::::B"   <<endl;
	COUT <<"K:::::::::::K          G:::::G      GGGGGGGGGG     B:::::::::::::BB"  <<COLORTEXT_RESET <<" - evolution" << COLORTEXT_CYAN <<endl; 
	COUT <<"K:::::::::::K          G:::::G      G::::::::G     B::::BBBBBB:::::B"   <<endl;
	COUT <<"K::::::K:::::K         G:::::G      GGGGG::::G     B::::B     B:::::B"  <<endl;
	COUT <<"K:::::K K:::::K        G:::::G          G::::G     B::::B     B:::::B"  << COLORTEXT_CYAN << endl;
	COUT <<"KK::::::K  K:::::KKK      G:::::G       G::::G     B:::::B     B:::::B"<<endl;
	COUT <<"K:::::::K   K::::::K       G:::::GGGGGGGG::::G     BB:::::BBBBBB::::::B"<<endl;
	COUT <<"K:::::::K    K:::::K        GG:::::::::::::::G     B:::::::::::::::::B" <<endl;
	COUT <<"K:::::::K    K:::::K          GGG::::::GGG:::G     B::::::::::::::::B"  <<endl;
	COUT <<COLORTEXT_LIGHT_BROWN<<"KKKKKKKKK    KKKKKKK            GGGGGGGGGGGGGG     BBBBBBBBBBBBBBBBB"   <<endl;
	COUT <<COLORTEXT_RESET << endl;


	COUT << "running on " << n*m << " cores." << endl;

	if (settingsfile == NULL)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		parallel.abortForce();
	}

	COUT << " initializing..." << endl;

	start_time = MPI_Wtime();

	numparam = loadParameterFile(settingsfile, params);

	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);

	COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;

	sprintf(filename, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
	saveParameterFile(filename, params, numparam);
	sprintf(filename, "%s%s_background.dat", sim.output_path, sim.basename_generic);
	outfile = fopen(filename, "w");
	
	if (sim.check_bg_file == 1)
	{
	sprintf(filename, "%s%s_check_bg.dat", sim.output_path, sim.basename_generic);
	check_file = fopen(filename, "w");
	}

	free(params);

  	#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)
		background class_background;
		thermo class_thermo;
		perturbs class_perturbs;

		if (precisionfile != NULL)
		{
			numparam = loadParameterFile(precisionfile, params);
		if(parallel.isRoot())  COUT << " Precision parameters are found and being read!" <<" The number of precision parameters are:"<<numparam << endl;
		}
		else numparam = 0;
		#ifdef HAVE_HICLASS_BG
			initializeCLASSstructures(sim, ic, cosmo, class_background, class_thermo, class_perturbs, params, numparam);
			loadBGFunctions(class_background, H_spline, "H [1/Mpc]", sim.z_in);
			loadBGFunctions(class_background, H_prime_spline, "H_prime", sim.z_in);
  			loadBGFunctions(class_background, H_prime_prime_spline, "H_prime_prime", sim.z_in);
			loadBGFunctions(class_background, rho_cdm_spline, "(.)rho_cdm", sim.z_in);
			loadBGFunctions(class_background, rho_b_spline, "(.)rho_b", sim.z_in);
			loadBGFunctions(class_background, rho_g_spline, "(.)rho_g", sim.z_in);
			loadBGFunctions(class_background, rho_crit_spline, "(.)rho_crit", sim.z_in);
			loadBGFunctions(class_background, rho_ur_spline, "(.)rho_ur", sim.z_in);
			loadBGFunctions(class_background, cs2_spline, "c_s^2", sim.z_in);
			loadBGFunctions(class_background, cs2_prime_spline, "c_s^2_prime", sim.z_in);
			loadBGFunctions(class_background, rho_smg_spline, "(.)rho_smg", sim.z_in);
			loadBGFunctions(class_background, rho_smg_prime_spline, "(.)rho_smg_prime", sim.z_in);
			loadBGFunctions(class_background, p_smg_spline, "(.)p_smg", sim.z_in);
			loadBGFunctions(class_background, p_smg_prime_spline, "(.)p_smg_prime", sim.z_in);
			loadBGFunctions(class_background, alpha_K_spline, "kineticity_smg", sim.z_in);
			loadBGFunctions(class_background, alpha_K_prime_spline, "kineticity_prime_smg", sim.z_in);
			loadBGFunctions(class_background, alpha_B_spline, "braiding_smg", sim.z_in);
			loadBGFunctions(class_background, alpha_B_prime_spline, "braiding_prime_smg", sim.z_in);
			loadBGFunctions(class_background, cs2num_spline, "cs2num", sim.z_in);
			loadBGFunctions(class_background, kin_D_spline, "kin (D)", sim.z_in);
			loadBGFunctions(class_background, lambda_2_spline, "lambda_2", sim.z_in);
		#endif
	#endif

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.output_path);

	box[0] = sim.numpts;
	box[1] = sim.numpts;
	box[2] = sim.numpts;

	Lattice lat(3,box,1);
	Lattice latFT;
	latFT.initializeRealFFT(lat,0);

	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_cdm;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_b;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_ncdm[MAX_PCL_SPECIES-2];
	Field<Real> * update_cdm_fields[3];
	Field<Real> * update_b_fields[3];
	Field<Real> * update_ncdm_fields[3];
	double f_params[5];
	set<long> IDbacklog[MAX_PCL_SPECIES];

	Field<Real> phi;
	Field<Cplx> scalarFT;
	phi.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> plan_phi(&phi, &scalarFT);

	Field<Real> source;
	source.initialize(lat,1);
	PlanFFT<Cplx> plan_source(&source, &scalarFT);

	Field<Real> chi;
	chi.initialize(lat,1);
	PlanFFT<Cplx> plan_chi(&chi, &scalarFT);

	Field<Real> Sij;
	Field<Cplx> SijFT;
	Sij.initialize(lat,3,3,symmetric);
	SijFT.initialize(latFT,3,3,symmetric);
	PlanFFT<Cplx> plan_Sij(&Sij, &SijFT);

	Field<Real> Bi;
	Field<Cplx> BiFT;
	Bi.initialize(lat,3);
	BiFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_Bi(&Bi, &BiFT);


	// KGB
	Field<Real> phi_old;
	Field<Cplx> scalarFT_phi_old;
	phi_old.initialize(lat,1);
	scalarFT_phi_old.initialize(latFT,1);
	PlanFFT<Cplx> plan_phi_old(&phi_old, &scalarFT_phi_old);
	
	Field<Real> zeta_half;
	Field<Cplx> scalarFT_zeta_half;
	zeta_half.initialize(lat,1);
	scalarFT_zeta_half.initialize(latFT,1);
	PlanFFT<Cplx> plan_zeta_half(&zeta_half, &scalarFT_zeta_half);
  
	Field<Real> chi_old;
	Field<Cplx> scalarFT_chi_old;
	chi_old.initialize(lat,1);
	scalarFT_chi_old.initialize(latFT,1);
	PlanFFT<Cplx> plan_chi_old(&chi_old, &scalarFT_chi_old);


	Field<Real> phi_prime;
	Field<Cplx> phi_prime_scalarFT;
	phi_prime.initialize(lat,1);
	phi_prime_scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> phi_prime_plan(&phi_prime, &phi_prime_scalarFT);


	Field<Real> psi_prime;
	Field<Cplx> psi_prime_scalarFT;
	psi_prime.initialize(lat,1);
	psi_prime_scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> psi_prime_plan(&psi_prime, &psi_prime_scalarFT);

	Field<Real> psi_half;
	Field<Cplx> psi_half_scalarFT;
	psi_half.initialize(lat,1);
	psi_half_scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> psi_half_plan(&psi_half, &psi_half_scalarFT);


	Field<Real> pi_k;
	Field<Cplx> scalarFT_pi;
	pi_k.initialize(lat,1);
	scalarFT_pi.initialize(latFT,1);
	PlanFFT<Cplx> plan_pi_k(&pi_k, &scalarFT_pi);

	Field<Real> T00_kgb;
	Field<Cplx> T00_kgbFT;
	T00_kgb.initialize(lat,1);
	T00_kgbFT.initialize(latFT,1);
	PlanFFT<Cplx> plan_T00_kgb(&T00_kgb, &T00_kgbFT);

	Field<Real> T0i_kgb;
	Field<Cplx> T0i_kgbFT;
	T0i_kgb.initialize(lat,3);
	T0i_kgbFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_T0i_kgb(&T0i_kgb, &T0i_kgbFT);

	Field<Real> Tij_kgb;
	Field<Cplx> Tij_kgbFT;
	Tij_kgb.initialize(lat,3,3,symmetric);
	Tij_kgbFT.initialize(latFT,3,3,symmetric);
	PlanFFT<Cplx> plan_Tij_kgb(&Tij_kgb, &Tij_kgbFT);	


	Field<Real> deltaPm;
	Field<Cplx> scalarFT_deltaPm;
	deltaPm.initialize(lat,1);
	scalarFT_deltaPm.initialize(latFT,1);
	PlanFFT<Cplx> plan_deltaPm(&deltaPm, &scalarFT_deltaPm);
	
	#ifdef CHECK_B
		Field<Real> Bi_check;
		Field<Cplx> BiFT_check;
		Bi_check.initialize(lat,3);
		BiFT_check.initialize(latFT,3);
		PlanFFT<Cplx> plan_Bi_check(&Bi_check, &BiFT_check);
	#endif
	#ifdef VELOCITY
		Field<Real> vi;
		Field<Cplx> viFT;
		vi.initialize(lat,3);
		viFT.initialize(latFT,3);
		PlanFFT<Cplx> plan_vi(&vi, &viFT);
		double a_old;
	#endif
	
	update_cdm_fields[0] = &phi;
	update_cdm_fields[1] = &chi;
	update_cdm_fields[2] = &Bi;

	update_b_fields[0] = &phi;
	update_b_fields[1] = &chi;
	update_b_fields[2] = &Bi;

	update_ncdm_fields[0] = &phi;
	update_ncdm_fields[1] = &chi;
	update_ncdm_fields[2] = &Bi;

	Site x(lat);
	rKSite kFT(latFT);

	dx = 1.0 / (double) sim.numpts;
	numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;

	for (i = 0; i < 3; i++) // particles may never move farther than to the adjacent domain
	{
		if (lat.sizeLocal(i)-1 < sim.movelimit)
			sim.movelimit = lat.sizeLocal(i)-1;
	}
	parallel.min(sim.movelimit);
	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT; // Just a definition to make Friedmann equation simplified! and working with normal numbers
	a = 1. / (1. + sim.z_in);
  	tau = particleHorizon
	(a, fourpiG,
    #ifdef HAVE_HICLASS_BG
    gsl_spline_eval(H_spline, 1., acc), class_background
    #else
    cosmo
    #endif
  	);

  	if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG,
    #ifdef HAVE_HICLASS_BG
      H_spline, acc
    #else
      cosmo
    #endif
  	) )
    dtau = sim.Cf * dx;
  	else
    dtau = sim.steplimit / 	Hconf(a, fourpiG,
      #ifdef HAVE_HICLASS_BG
        H_spline, acc
      #else
        cosmo
      #endif
      );

	dtau_old = 0.;
	if (ic.generator == ICGEN_BASIC)
		generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &pi_k, &zeta_half, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_zeta_half, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, params, numparam);
	// generates ICs on the fly
	else if (ic.generator == ICGEN_READ_FROM_DISK)
  	readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount, IDbacklog, params, numparam);
	#ifdef ICGEN_PREVOLUTION
		else if (ic.generator == ICGEN_PREVOLUTION)
			generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, params, numparam);
	#endif
	#ifdef ICGEN_FALCONIC
		else if (ic.generator == ICGEN_FALCONIC)
			maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel+1, &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
	#endif
	else
	{
		COUT << " error: IC generator not implemented!" << endl;
		parallel.abortForce();
	}

	if (sim.baryon_flag > 1)
	{
		COUT << " error: baryon_flag > 1 after IC generation, something went wrong in IC generator!" << endl;
		parallel.abortForce();
	}

	numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;
	parallel.max<double>(maxvel, numspecies);

	if (sim.gr_flag > 0)
	{
		for (i = 0; i < numspecies; i++)
			maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
	}

	#ifdef CHECK_B
		if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
			for (kFT.first(); kFT.test(); kFT.next())
			{
				BiFT_check(kFT, 0) = BiFT(kFT, 0);
				BiFT_check(kFT, 1) = BiFT(kFT, 1);
				BiFT_check(kFT, 2) = BiFT(kFT, 2);
			}
		}
	#endif
	#ifdef VELOCITY
		a_old = a;
		projection_init(&vi);
	#endif

	#ifdef BENCHMARK
		initialization_time = MPI_Wtime() - start_time;
		parallel.sum(initialization_time);
		COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
	#else
		COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
	#endif

	#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)

		if(sim.fluid_flag > 0 )
		{
			if(parallel.isRoot())  cout << " \033[1;31merror:\033[0m"<< " \033[1;31merror: You are using k-evolution and asking for fluid KGB treatment at the same time! Don't know what to do!  \033[0m" << endl;
			parallel.abortForce();
		}
			if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
			{
				initializeCLASSstructures(sim, ic, cosmo, class_background, class_thermo, class_perturbs, params, numparam);
				if (sim.gr_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
				{
					prepareFTchiLinear(class_background, class_perturbs, scalarFT, sim, ic, cosmo, fourpiG, a);
					plan_source.execute(FFT_BACKWARD);
					for (x.first(); x.test(); x.next())
						chi(x) += source(x);
					chi.updateHalo();
				}
			}	
		if (numparam > 0) free(params);
	#endif

	while (true)    // main loop
	{
		for (x.first(); x.test(); x.next())
			{
				phi_old(x) = phi(x);
				chi_old(x) = chi(x);
			}
		#ifdef BENCHMARK
				cycle_start_time = MPI_Wtime();
		#endif
				// construct stress-energy tensor
				projection_init(&source);
		#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)
				if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
					projection_T00_project(class_background, class_perturbs, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
		#endif
		if (sim.gr_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			if (sim.baryon_flag)
				projection_T00_project(&pcls_b, &source, a, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_T00_project(pcls_ncdm+i, &source, a, &phi);
				else if (sim.radiation_flag == 0 || (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] == 0))
				{
					tmp = bg_ncdm(a, cosmo, i);
					for(x.first(); x.test(); x.next())
						source(x) += tmp;
				}
			}
		}
		else
		{
			scalarProjectionCIC_project(&pcls_cdm, &source);
			if (sim.baryon_flag)
				scalarProjectionCIC_project(&pcls_b, &source);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
			}
		}
		projection_T00_comm(&source);
		#ifdef VELOCITY
			if ((sim.out_pk & MASK_VEL) || (sim.out_snapshot & MASK_VEL))
			{
				projection_init(&Bi);
				projection_Ti0_project(&pcls_cdm, &Bi, &phi, &chi);
				vertexProjectionCIC_comm(&Bi);
				compute_vi_rescaled(cosmo, &vi, &source, &Bi, a, a_old
				#ifdef HAVE_HICLASS_BG
				, H_spline, acc
				#endif
				);
				a_old = a;
			}
		#endif

		if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
			projection_init(&Bi);
			projection_T0i_project(&pcls_cdm, &Bi, &phi);
			if (sim.baryon_flag)
				projection_T0i_project(&pcls_b, &Bi, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_Bncdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_T0i_project(pcls_ncdm+i, &Bi, &phi);
			}
			projection_T0i_comm(&Bi);
		}

		projection_init(&Sij);
		projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
		if (sim.baryon_flag)
			projection_Tij_project(&pcls_b, &Sij, a, &phi);
		if (a >= 1. / (sim.z_switch_linearchi + 1.))
		{
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
			}
		}
		projection_Tij_comm(&Sij);

		#ifdef BENCHMARK
			projection_time += MPI_Wtime() - cycle_start_time;
			ref_time = MPI_Wtime();
		#endif

	if (sim.kgb_source_gravity==1)
		{
			Hc = Hconf(a, fourpiG,
			#ifdef HAVE_HICLASS_BG
				H_spline, acc
			#else
				cosmo
			#endif
			);
		// KGB projection Tmunu
		// In the projection zeta_integer comes, since synched with particles..
		#ifdef HAVE_HICLASS_BG // hiclass used to provide quantities!
			// CHECK! the coeffs and etc!
			if (sim.vector_flag == VECTOR_ELLIPTIC)
			{
				projection_Tmunu_kgb(T00_kgb, T0i_kgb, Tij_kgb, dx, a, fourpiG, gsl_spline_eval(H_spline, 1., acc), phi, chi, phi_prime, pi_k, zeta_half, deltaPm, 
				Hconf(a, fourpiG, H_spline, acc), Hconf_prime(a, fourpiG, H_spline, acc), Hconf_prime_prime(a, fourpiG, H_spline, acc),
				gsl_spline_eval(rho_smg_spline, a, acc), gsl_spline_eval(p_smg_spline, a, acc), gsl_spline_eval(p_smg_prime_spline, a, acc), gsl_spline_eval(rho_crit_spline, 1., acc),
				gsl_spline_eval(alpha_K_spline, a, acc), gsl_spline_eval(alpha_B_spline, a, acc), gsl_spline_eval(alpha_K_prime_spline, a, acc), 
				gsl_spline_eval(alpha_B_prime_spline, a, acc));
			}
			else
			{
				projection_Tmunu_kgb(T00_kgb, T0i_kgb, Tij_kgb, dx, a, fourpiG, gsl_spline_eval(H_spline, 1., acc), phi, chi, phi_prime, pi_k, zeta_half, deltaPm, 
				Hconf(a, fourpiG, H_spline, acc), Hconf_prime(a, fourpiG, H_spline, acc), Hconf_prime_prime(a, fourpiG, H_spline, acc),
				gsl_spline_eval(rho_smg_spline, a, acc),  gsl_spline_eval(p_smg_spline, a, acc), gsl_spline_eval(p_smg_prime_spline, a, acc), gsl_spline_eval(rho_crit_spline, 1., acc),
				gsl_spline_eval(alpha_K_spline, a, acc), gsl_spline_eval(alpha_B_spline, a, acc), gsl_spline_eval(alpha_K_prime_spline, a, acc), 
				gsl_spline_eval(alpha_B_prime_spline, a, acc));
			}
		#else // default KGB-evolution or CLASS // No hiclass BG used
			if (sim.vector_flag == VECTOR_ELLIPTIC)
			{
				projection_Tmunu_kgb(T00_kgb, T0i_kgb, Tij_kgb, dx, a, phi, pi_k, zeta_half, cosmo.Omega_kgb * pow(a , -3. * cosmo.w_kgb) * (1. + cosmo.w_kgb) / (cosmo.cs2_kgb), cosmo.Omega_kgb * pow(a , -3. * cosmo.w_kgb) * (1. + cosmo.w_kgb), cosmo.w_kgb, cosmo.cs2_kgb, Hc, sim.NL_kgb, 1);
			}
			else
			{
				projection_Tmunu_kgb(T00_kgb, T0i_kgb, Tij_kgb, dx, a, phi, pi_k, zeta_half, cosmo.Omega_kgb * pow(a , -3. * cosmo.w_kgb) * (1. + cosmo.w_kgb) / (cosmo.cs2_kgb), cosmo.Omega_kgb * pow(a , -3. * cosmo.w_kgb) * (1. + cosmo.w_kgb), cosmo.w_kgb, cosmo.cs2_kgb, Hc, sim.NL_kgb, 0);
			}
		#endif

		for (x.first(); x.test(); x.next())
		{   // CHECK! the coeffs and etc!
			// The coefficient is because it wanted to to be sourced according to eq C.2 of gevolution paper
			// Note that it is multiplied to dx^2 and is divived by -a^3 because of definition of T00 which is scaled by a^3
			// We have T00 and Tij according to code's units, but source is important to calculate potentials and moving particles.
			// There is coefficient between Tij and Sij as source.
			source(x) += T00_kgb(x);
			if (sim.vector_flag == VECTOR_ELLIPTIC) for (int c=0;c<3;c++) Bi(x,c) += T0i_kgb(x,c);
			for(int c=0;c<6;c++) Sij(x,c) += (2.) * Tij_kgb(x,c);
		}
		}
		#ifdef BENCHMARK
			kgb_update_time += MPI_Wtime() - ref_time;
			ref_time = MPI_Wtime();
		#endif
		// KGB projection Tmunu end	

		if (sim.gr_flag > 0)
		{
			T00hom = 0.;
			for (x.first(); x.test(); x.next())
				T00hom += source(x);
			parallel.sum<double>(T00hom);
			T00hom /= (double) numpts3d;

			T00KGBhom= 0.;
			for (x.first(); x.test(); x.next())
				T00KGBhom+= T00_kgb(x);
			parallel.sum<double>(T00KGBhom);
			T00KGBhom /= (double) numpts3d;

			if (cycle % CYCLE_INFO_INTERVAL == 0)
			{
				COUT << " cycle " << cycle << ", background information: z = " << (1./a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
			}

			if (dtau_old > 0.)
			{
				Hc = Hconf(a, fourpiG,
				#ifdef HAVE_HICLASS_BG
				H_spline, acc
				#else
				cosmo
				#endif
				);
				prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hc * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hc * Hc * dx * dx);  // prepare nonlinear source for phi update
				#ifdef BENCHMARK
					ref2_time= MPI_Wtime();
				#endif
					plan_source.execute(FFT_FORWARD);  // go to k-space
				#ifdef BENCHMARK
					fft_time += MPI_Wtime() - ref2_time;
					fft_count++;
				#endif

				solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hc / dtau_old);  // phi update (k-space)

				#ifdef BENCHMARK
					ref2_time= MPI_Wtime();
				#endif
					plan_phi.execute(FFT_BACKWARD);	 // go back to position space
				#ifdef BENCHMARK
					fft_time += MPI_Wtime() - ref2_time;
					fft_count++;
				#endif
			}
		}
		else
		{
			#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
			#endif
				plan_source.execute(FFT_FORWARD);  // Newton: directly go to k-space
			#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
			#endif

			solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Newton: phi update (k-space)

			#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
			#endif
				plan_phi.execute(FFT_BACKWARD);	 // go back to position space
			#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
			#endif
		}

		phi.updateHalo();  // communicate halo values

		// // Remove dark energy contribution from source after solving the Poisson equation
		// for (x.first(); x.test(); x.next())
		// {   
		// 	source(x) -= T00_kgb(x);
		// 	for(int c=0;c<6;c++) Sij(x,c) -= (2.) * Tij_kgb(x,c);
		// }

		// record some background data
		if (kFT.setCoord(0, 0, 0))
		    {
        // Use snprintf for safer string formatting
        char filename[512]; // Adjust the buffer size as needed
        snprintf(filename, sizeof(filename), "%s%s_background.dat", sim.output_path, sim.basename_generic);
        FILE* outfile = fopen(filename, "a");
        
        if (outfile == NULL)
        {
            std::cerr << "Error opening file for background output!" << std::endl;
            // Handle the error appropriately (e.g., exit or clean up)
        }
        else
        {
    #ifdef HAVE_HICLASS_BG
            if (cycle == 0)
            {
                // Write multi-line header with fixed-width fields
                fprintf(outfile, "# background statistics\n");
                fprintf(outfile, "# All the values except 'a' and 'T00' are computed in the hiclass.\n");
                fprintf(outfile, "# The quantities rho_x, rho_x_prime are in unit of hiclass, T00 is in the k-evolution code's unit\n");		
				fprintf(outfile, "# Constant values at end of file:\n");
				fprintf(outfile, "# fourpiG   = %24e\n", fourpiG);
				fprintf(outfile, "# H0[1/Mpc] = %24e\n", gsl_spline_eval(H_spline, 1, acc));
                
                // Header line with fixed-width fields (25 characters each)
                fprintf(outfile, "\n# %-12s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s\n",
                    "0:cycle",
                    "1:tau/boxsize",
                    "2:a",
					"3:z",
					"4:H [1/Mpc]",
                    "5:Hconf",
                    "6:Hconf_prime",
                    "7:Hconf_prime_prime",
					"8:rho_cdm",
					"9:rho_b",
					"10:rho_g",
					"11:rho_ur",
					"12:rho_crit",
					"13:rho_smg",
                    "14:p_smg",
                    "15:rho_smg_prime",
                    "16:p_smg_prime",
                    "17:alpha_K",
                    "18:alpha_B",
                    "19:alpha_K_prime",
                    "20:alpha_B_prime",
                    "21:cs2",
					"22:cs2num",
					"23:kin (D)",
					"24:phi(k=0)",
					"25:T00_hom",
					"26:T00_KGB_hom"
                );
            }
    
            // Define a format string with fixed-width fields for alignment (25 characters each)
            // Left-align each field using the '-' flag
            const char* format = " %-15d %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e\n";
    
            // Write the data with alignment
            fprintf(outfile, format,
                    cycle,                                                // 0:cycle (int)
                    tau,                                                  // 1:tau/boxsize (double)
                    a,                                                    // 2:a (double)
					1./a-1.,                                              // 3:z
                    gsl_spline_eval(H_spline, a, acc),                    // 4:H [1/Mpc]
                    Hconf(a, fourpiG, H_spline, acc),                     // 5:Hconf
                    Hconf_prime(a, fourpiG, H_spline, acc),               // 6:Hconf_prime
					Hconf_prime_prime(a, fourpiG, H_spline, acc),         // 7:Hconf_prime_prime
                    gsl_spline_eval(rho_cdm_spline, a, acc),              // 8:rho_cdm
                    gsl_spline_eval(rho_b_spline, a, acc),                // 9:rho_b
                    gsl_spline_eval(rho_g_spline, a, acc),                // 10:rho_g
                    gsl_spline_eval(rho_ur_spline, a, acc),               // 11:rho_ur
                    gsl_spline_eval(rho_crit_spline, a, acc),             // 12:rho_crit
                    gsl_spline_eval(rho_smg_spline, a, acc),              // 13:rho_smg
                    gsl_spline_eval(p_smg_spline, a, acc),                // 14:p_smg
                    gsl_spline_eval(rho_smg_prime_spline, a, acc),        // 15:rho_smg_prime
                    gsl_spline_eval(p_smg_prime_spline, a, acc),          // 16:p_smg_prime
                    gsl_spline_eval(alpha_K_spline, a, acc),              // 17:alpha_K
                    gsl_spline_eval(alpha_B_spline, a, acc),              // 18:alpha_B
                    gsl_spline_eval(alpha_K_prime_spline, a, acc),        // 19:alpha_K_prime
                    gsl_spline_eval(alpha_B_prime_spline, a, acc),        // 20:alpha_B_prime
					gsl_spline_eval(cs2_spline, a, acc),                  // 21:cs2
                    gsl_spline_eval(cs2num_spline, a, acc),               // 22:cs2num
                    gsl_spline_eval(kin_D_spline, a, acc),                // 23:kin (D)
					scalarFT(kFT).real(),                                 // 24:phi(k=0)
					T00hom,                                               // 25:T00hom
					T00KGBhom                                             // 26:T00KGBhom
				); 
				

            fclose(outfile);
    #else
            if (cycle == 0)
            {
                // Write multi-line header with fixed-width fields (25 characters each)
                fprintf(outfile, "# background statistics\n");
                fprintf(outfile, "# %-25s %-25s %-25s %-25s %-25s %-25s %-25s\n",
                    "cycle",
                    "tau/boxsize",
                    "a",
                    "conformal H/H0",
                    "Hconf'/H0^2",
                    "phi(k=0)",
                    "T00(k=0)"
                );
            }

            // Define a format string with fixed-width fields for alignment (25 characters each)
            const char* format = " %-25d %-25e %-25e %-25e %-25e %-25e %-25e\n";

            // Write the data with alignment
            fprintf(outfile, format,
                    cycle, 
                    tau, 
                    a, 
                    Hconf(a, fourpiG, cosmo) / Hconf(1.0, fourpiG, cosmo), 
                    Hconf_prime(a, fourpiG, cosmo) / (Hconf(1.0, fourpiG, cosmo) * Hconf(1.0, fourpiG, cosmo)), 
                    scalarFT(kFT).real(), 
                    T00hom);
            fclose(outfile);
    #endif
        }
			}
		// done recording background data

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Consistency Check //////////////////////////////////////////////////////////////////

		// record some background data for consistency check
		if (sim.check_bg_file == 1)
		{
			if (kFT.setCoord(0, 0, 0))
  {
        // Use snprintf for safer string formatting
        snprintf(filename, sizeof(filename), "%s%s_check_bg.dat", sim.output_path, sim.basename_generic);
        check_file = fopen(filename, "a");
        
        if (check_file == NULL)
        {
            std::cerr << "Error opening file for consistency check for background output!" << std::endl;
            // Handle the error appropriately (e.g., exit or clean up)
        }
        else
        {
    #ifdef HAVE_HICLASS_BG
            if (cycle == 0)
            {
                // Write multi-line header with fixed-width fields
                fprintf(check_file, "# background check for consistency check regarding the implementation\n");
                fprintf(check_file, "# All the values except 'a' and 'T00' are computed in the hiclass.\n");
                fprintf(check_file, "# The quantities rho_x, rho_x_prime are in unit of hiclass, T00 is in the k-evolution code's unit\n");
                
                // Header line with fixed-width fields (25 characters each)
                fprintf(check_file, "# %-12s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s %-24s\n",
                    "0:cycle",
                    "1:tau/boxsize",
                    "2:a",
                    "3:Hconf_F",
                    "4:Hconf",
                    "5:Hconf_prime_F",
                    "6:Hconf_prime",
                    "7:Hconf_prime_prime_F",
                    "8:Hconf_prime_prime",
                    "9:fourpiG",
                    "10:H0[1/Mpc]",
                    "11:alpha_K",
                    "12:alpha_B",
                    "13:alpha_K_prime",
                    "14:alpha_B_prime",
                    "15:cs2",
                    "16:rho_smg",
                    "17:p_smg",
                    "18:rho_smg_prime",
                    "19:p_smg_prime",
                    "20:rho_cdm ",
                    "21:rho_b",
                    "22:rho_crit",
					"23:cs2num",
					"24:kin (D)",
					"25:H [1/Mpc]",
					"26:phi(k=0)",
					"27:T00hom",
					"28:T00KGBhom"
                );
            }
    
            // Define a format string with fixed-width fields for alignment (25 characters each)
            // Right-align numeric data for better readability
            const char* format = " %-15d %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e %-24e\n";
    
            // Write the data with alignment
            fprintf(check_file, format,
                    cycle,                                                // 0:cycle (int)
                    tau,                                                  // 1:tau/boxsize (double)
                    a,                                                    // 2:a (double)
                    Hconf_F(a, fourpiG, cosmo),                           // 3:Hconf_F (double)
                    Hconf(a, fourpiG, H_spline, acc),                     // 4:Hconf (double)
                    Hconf_prime_F(a, fourpiG, cosmo),                     // 5:Hconf_prime_F (double)
                    Hconf_prime(a, fourpiG, H_spline, acc),               // 6:Hconf_prime (double)
                    Hconf_prime_prime_F(a, fourpiG, cosmo),               // 7:Hconf_prime_prime_F (double)
                    Hconf_prime_prime(a, fourpiG, H_spline, acc),         // 8:Hconf_prime_prime (double)
                    fourpiG,                                              // 9:fourpiG (double)
                    gsl_spline_eval(H_spline, 1, acc),                    // 10:H0[1/Mpc] (double)
                    gsl_spline_eval(alpha_K_spline, a, acc),              // 11:alpha_K (double)
                    gsl_spline_eval(alpha_B_spline, a, acc),              // 12:alpha_B (double)
                    gsl_spline_eval(alpha_K_prime_spline, a, acc),        // 13:alpha_K_prime (double)
                    gsl_spline_eval(alpha_B_prime_spline, a, acc),        // 14:alpha_B_prime (double)
                    gsl_spline_eval(cs2_spline, a, acc),                  // 15:cs2 (double)
                    gsl_spline_eval(rho_smg_spline, a, acc),              // 16:rho_smg (double)
                    gsl_spline_eval(p_smg_spline, a, acc),                // 17:p_smg (double)
                    gsl_spline_eval(rho_smg_prime_spline, a, acc),        // 18:rho_smg_prime (double)
                    gsl_spline_eval(p_smg_prime_spline, a, acc),          // 19:p_smg_prime (double)
                    gsl_spline_eval(rho_cdm_spline, a, acc),              // 20:rho_cdm (double)
                    gsl_spline_eval(rho_b_spline, a, acc),                // 21:rho_b (double)
                    gsl_spline_eval(rho_crit_spline, a, acc),             // 22:rho_crit (double)
					gsl_spline_eval(cs2num_spline, a, acc),               // 23:cs2num
					gsl_spline_eval(kin_D_spline, a, acc),                 // 23:kin D
					gsl_spline_eval(H_spline, a, acc),                     // 25:H [1/Mpc]
					scalarFT(kFT).real(),                                 // phi(k=0)
					T00hom,                                               // T00hom
					T00KGBhom                                             //T00KGBhom
					);                                              
            fclose(check_file);
    #else
            if (cycle == 0)
            {
                // Write multi-line header with fixed-width fields
                fprintf(check_file, "# background statistics\n");
                
                // Header line with fixed-width fields (25 characters each)
                fprintf(check_file, "# %-25s %-25s %-25s %-25s %-25s %-25s %-25s\n",
                    "cycle",
                    "tau/boxsize",
                    "a",
                    "conformal H/H0",
                    "Hconf'/H0^2",
                    "phi(k=0)",
                    "T00(k=0)"
                );
            }
    
            // Define a format string with fixed-width fields for alignment (25 characters each)
            const char* format = " %25d %25e %25e %25e %25e %25e %25e\n";
    
            // Write the data with alignment
            fprintf(check_file, format,
                    cycle, 
                    tau, 
                    a, 
                    Hconf(a, fourpiG, cosmo) / Hconf(1.0, fourpiG, cosmo), 
                    Hconf_prime(a, fourpiG, cosmo) / (Hconf(1.0, fourpiG, cosmo) * Hconf(1.0, fourpiG, cosmo)), 
                    scalarFT(kFT).real(), 
                    T00hom);
            fclose(check_file);
    #endif
        }
		}
		}
		// done recording background data for consistency check

		prepareFTsource<Real>(phi, Sij, Sij, deltaPm, 2. * fourpiG * dx * dx / a);  // prepare nonlinear source for additional equations

		#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
		#endif
			plan_Sij.execute(FFT_FORWARD);  // go to k-space
		#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count += 6;
		#endif

		#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)
			if (sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.))
			{
				prepareFTchiLinear(class_background, class_perturbs, scalarFT, sim, ic, cosmo, fourpiG, a);
				projectFTscalar(SijFT, scalarFT, 1);
			}
			else
		#endif
		projectFTscalar(SijFT, scalarFT);  // construct chi by scalar projection (k-space)

		#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
		#endif
			plan_chi.execute(FFT_BACKWARD);	 // go back to position space
		#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
		#endif
			chi.updateHalo();  // communicate halo values

				if (sim.vector_flag == VECTOR_ELLIPTIC)
				{
		#ifdef BENCHMARK
					ref2_time= MPI_Wtime();
		#endif
					plan_Bi.execute(FFT_FORWARD);
		#ifdef BENCHMARK
					fft_time += MPI_Wtime() - ref2_time;
					fft_count++;
		#endif
					projectFTvector(BiFT, BiFT, fourpiG * dx * dx); // solve B using elliptic constraint (k-space)
		#ifdef CHECK_B
					evolveFTvector(SijFT, BiFT_check, a * a * dtau_old);
		#endif
				}
				else
					evolveFTvector(SijFT, BiFT, a * a * dtau_old);  // evolve B using vector projection (k-space)

				if (sim.gr_flag > 0)
				{
		#ifdef BENCHMARK
					ref2_time= MPI_Wtime();
		#endif
					plan_Bi.execute(FFT_BACKWARD);  // go back to position space
		#ifdef BENCHMARK
					fft_time += MPI_Wtime() - ref2_time;
					fft_count += 3;
		#endif
					Bi.updateHalo();  // communicate halo values
				}

		#ifdef BENCHMARK
				gravity_solver_time += MPI_Wtime() - ref_time;
				ref_time = MPI_Wtime();
		#endif


		// lightcone output
		if (sim.num_lightcone > 0)
			writeLightcones(sim, cosmo, fourpiG, a, tau, dtau, dtau_old, maxvel[0], cycle, h5filename + sim.basename_lightcone,
			#ifdef HAVE_HICLASS_BG
			class_background, H_spline, acc,
			#endif
		&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &Sij, &BiFT, &SijFT, &plan_Bi, &plan_Sij, done_hij, IDbacklog);
		else done_hij = 0;

		#ifdef BENCHMARK
		lightcone_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
		#endif

		// snapshot output
		if (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

			writeSnapshots(sim, cosmo, fourpiG, a, dtau_old, done_hij, snapcount, h5filename + sim.basename_snapshot,
			#ifdef HAVE_HICLASS_BG
			H_spline, acc,
			#endif
			&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_kgb, &T0i_kgb, &Tij_kgb, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
				#ifdef CHECK_B
						, &Bi_check, &BiFT_check, &plan_Bi_check
				#endif
				#ifdef VELOCITY
						, &vi
				#endif
			);

			snapcount++;
		}

		#ifdef BENCHMARK
				snapshot_output_time += MPI_Wtime() - ref_time;
				ref_time = MPI_Wtime();
		#endif

		// power spectra
		if (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
			writeSpectra(sim, cosmo, fourpiG, a, pkcount,
			#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)
							class_background, class_perturbs, ic,
			#ifdef HAVE_HICLASS_BG
					H_spline, acc, gsl_spline_eval(rho_smg_spline, a, acc), gsl_spline_eval(rho_crit_spline, 1., acc),
			#endif
			#endif
							&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_kgb, &T0i_kgb, &Tij_kgb ,&source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_zeta_half, &BiFT, &T00_kgbFT, &T0i_kgbFT, &Tij_kgbFT, &SijFT, &plan_phi, &plan_pi_k , &plan_zeta_half, &plan_chi, &plan_Bi, &plan_T00_kgb, &plan_T0i_kgb, &plan_Tij_kgb, &plan_source, &plan_Sij
			#ifdef CHECK_B
							, &Bi_check, &BiFT_check, &plan_Bi_check
			#endif
			#ifdef VELOCITY
							, &vi, &viFT, &plan_vi
			#endif
				);

			if (sim.out_pk & MASK_PHI_PRIME)
			{
				writeSpectra_phi_prime(sim, cosmo, fourpiG,  a, pkcount,
					#ifdef HAVE_HICLASS_BG
					 H_spline, acc,
					#endif
					&phi_prime ,&phi_prime_scalarFT ,  &phi_prime_plan);

			}

			pkcount++;
		}

		#ifdef EXACT_OUTPUT_REDSHIFTS
				tmp = a;
			rungekutta4bg(tmp, fourpiG,
			#ifdef HAVE_HICLASS_BG
				H_spline, acc,
			#else
				cosmo,
			#endif
			0.5 * dtau);
			rungekutta4bg(tmp, fourpiG,
			#ifdef HAVE_HICLASS_BG
				H_spline, acc,
			#else
				cosmo,
			#endif
			0.5 * dtau);

			if (pkcount < sim.num_pk && 1. / tmp < sim.z_pk[pkcount] + 1.)
			{
					writeSpectra(sim, cosmo, fourpiG, a, pkcount,
		#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)
								class_background, class_perturbs, ic,
			#ifdef HAVE_HICLASS_BG
				H_spline, acc, gsl_spline_eval(rho_smg_spline, a, acc), gsl_spline_eval(rho_crit_spline, 1., acc),
			#endif
		#endif
								&pcls_cdm, &pcls_b, pcls_ncdm, &phi,&pi_k, &zeta_half, &chi, &Bi, &T00_kgb, &T0i_kgb, &Tij_kgb, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_zeta_half, &BiFT, &T00_kgbFT, &T0i_kgbFT, &Tij_kgbFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_T00_kgb, &plan_T0i_kgb, &plan_Tij_kgb, &plan_source, &plan_Sij
			#ifdef CHECK_B
								, &Bi_check, &BiFT_check, &plan_Bi_check
			#endif
			#ifdef VELOCITY
							, &vi, &viFT, &plan_vi
			#endif
					);
			
			
				if (sim.out_pk & MASK_PHI_PRIME)
				{
					writeSpectra_phi_prime(sim, cosmo, fourpiG,  a, pkcount,
						#ifdef HAVE_HICLASS_BG
						 H_spline, acc,
						#endif
						&phi_prime ,&phi_prime_scalarFT ,  &phi_prime_plan);
	
				}		
			}

  		#endif // EXACT_OUTPUT_REDSHIFTS


		#ifdef BENCHMARK
				spectra_output_time += MPI_Wtime() - ref_time;
		#endif

		if (pkcount >= sim.num_pk && snapcount >= sim.num_snapshot)
		{
			for (i = 0; i < sim.num_lightcone; i++)
			{
				if (sim.lightcone[i].z + 1. < 1. / a)
					i = sim.num_lightcone + 1;
			}
			if (i == sim.num_lightcone) break; // simulation complete
		}

		// compute number of step subdivisions for ncdm particle updates
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (dtau * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
				numsteps_ncdm[i] = (int) ceil(dtau * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
			else numsteps_ncdm[i] = 1;
		}

		if (cycle % CYCLE_INFO_INTERVAL == 0)
		{
			COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
			if (sim.baryon_flag)
			{
				COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx;
			}

      COUT << "), time step / Hubble time = " << Hconf(a, fourpiG,
      #ifdef HAVE_HICLASS_BG
        H_spline, acc
      #else
        cosmo
      #endif
      ) * dtau;

			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (i == 0)
				{
					COUT << endl << " time step subdivision for ncdm species: ";
				}
				COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel[i+1+sim.baryon_flag] << ")";
				if (i < cosmo.num_ncdm-1)
				{
					COUT << ", ";
				}
			}

			COUT << endl;
		}

		// KGB loop start!
		#ifdef BENCHMARK
			ref_time = MPI_Wtime();
		#endif
		#ifdef HAVE_HICLASS_BG // If we have BG vlaues from hicalss/CLASS!
			derivatives_update(dtau_old, cycle, phi, phi_old, chi, chi_old, phi_prime, psi_prime); // The derivatives of phi and psi computed at step n! At cycle 0 they are 0! We should use dtau not dtau_old to be the derivative at the requested time similar to the way we update the background a_n -> a_n+1 where we use dtau!

			a_kgb = a;
			if(cycle==0)
			{
				update_zeta(-dtau/ (2.0 * sim.n_kgb_numsteps), dx, a_kgb, fourpiG, gsl_spline_eval(H_spline, 1., acc), phi, chi, phi_prime, pi_k, zeta_half, deltaPm, 
				Hconf(a_kgb, fourpiG, H_spline, acc), Hconf_prime(a_kgb, fourpiG, H_spline, acc), Hconf_prime_prime(a_kgb, fourpiG, H_spline, acc),
				 gsl_spline_eval(rho_smg_spline, a_kgb, acc), gsl_spline_eval(p_smg_spline, a_kgb, acc), gsl_spline_eval(p_smg_prime_spline, a_kgb, acc), gsl_spline_eval(rho_crit_spline, 1., acc),
				 gsl_spline_eval(alpha_K_spline, a_kgb, acc), gsl_spline_eval(alpha_B_spline, a_kgb, acc), gsl_spline_eval(alpha_K_prime_spline, a_kgb, acc), 
				 gsl_spline_eval(alpha_B_prime_spline, a_kgb, acc), sim.NL_kgb);
				zeta_half.updateHalo();
			}
			for (i=0;i<sim.n_kgb_numsteps;i++)
			{
				update_zeta(dtau/ sim.n_kgb_numsteps, dx, a_kgb, fourpiG, gsl_spline_eval(H_spline, 1., acc), phi, chi, phi_prime, pi_k, zeta_half, deltaPm, 
				Hconf(a_kgb, fourpiG, H_spline, acc), Hconf_prime(a_kgb, fourpiG, H_spline, acc), Hconf_prime_prime(a_kgb, fourpiG, H_spline, acc),
				 gsl_spline_eval(rho_smg_spline, a_kgb, acc), gsl_spline_eval(p_smg_spline, a_kgb, acc), gsl_spline_eval(p_smg_prime_spline, a_kgb, acc), gsl_spline_eval(rho_crit_spline, 1., acc),
				 gsl_spline_eval(alpha_K_spline, a_kgb, acc), gsl_spline_eval(alpha_B_spline, a_kgb, acc), gsl_spline_eval(alpha_K_prime_spline, a_kgb, acc), 
				 gsl_spline_eval(alpha_B_prime_spline, a_kgb, acc), sim.NL_kgb);
				zeta_half.updateHalo();

				rungekutta4bg(a_kgb, fourpiG, H_spline, acc, dtau  / sim.n_kgb_numsteps / 2.0);
				update_pi(dtau/ sim.n_kgb_numsteps, dtau,  phi, chi, psi_prime, pi_k, zeta_half, Hconf(a_kgb, fourpiG, H_spline, acc)); // H_old is updated here in the function
				pi_k.updateHalo();
				rungekutta4bg(a_kgb, fourpiG, H_spline, acc, dtau  / sim.n_kgb_numsteps / 2.0);
			}
		#else // If not HAVE_HICLASS_BG We use  KGB-evolution with w, c_s^2 constants.
			derivatives_update(dtau_old, cycle, phi, phi_old, chi, chi_old, phi_prime, psi_prime); // The derivatives of phi and psi computed at step n! At cycle 0 they are 0! We should use dtau not dtau_old to be the derivative at the requested time similar to the way we update the background a_n -> a_n+1 where we use dtau!
			a_kgb = a;
			//First we update zeta_half to have it at -1/2 just in the first loop
			if(cycle==0)
			{
				update_zeta(-dtau/(2.0 * sim.n_kgb_numsteps) , dx, a_kgb, phi, chi, phi_prime, pi_k, zeta_half, cosmo.cs2_kgb, 0., cosmo.w_kgb, Hconf(a_kgb, fourpiG, cosmo), Hconf_prime(a_kgb, fourpiG, cosmo), sim.NL_kgb);
				zeta_half.updateHalo();
			}
			for (i=0;i<sim.n_kgb_numsteps;i++)
			{
				update_zeta(dtau/ sim.n_kgb_numsteps, dx, a_kgb, phi, chi, phi_prime, pi_k, zeta_half, cosmo.cs2_kgb, 0., cosmo.w_kgb, Hconf(a_kgb, fourpiG, cosmo), Hconf_prime(a_kgb, fourpiG, cosmo), sim.NL_kgb);
				zeta_half.updateHalo();
				rungekutta4bg(a_kgb, fourpiG, cosmo, dtau/sim.n_kgb_numsteps/2.0);
				update_pi(dtau/ sim.n_kgb_numsteps, dtau, phi, chi, psi_prime, pi_k, zeta_half, Hconf(a_kgb, fourpiG, cosmo)); 
				pi_k.updateHalo();
				rungekutta4bg(a_kgb, fourpiG, cosmo, dtau/sim.n_kgb_numsteps/2.0 );
			}
		#endif // KGB - LeapFrog: End


		#ifdef BENCHMARK
			kgb_update_time += MPI_Wtime() - ref_time;
			ref_time = MPI_Wtime();
		#endif
		
		#ifdef BENCHMARK
				ref2_time = MPI_Wtime();
		#endif
		for (i = 0; i < cosmo.num_ncdm; i++) // non-cold DM particle update
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;

			tmp = a;

			for (j = 0; j < numsteps_ncdm[i]; j++)
			{
				f_params[0] = tmp;
				f_params[1] = tmp * tmp * sim.numpts;
				if (sim.gr_flag > 0)
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				else
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

		#ifdef BENCHMARK
			update_q_count++;
			update_q_time += MPI_Wtime() - ref2_time;
			ref2_time = MPI_Wtime();
		#endif

		rungekutta4bg(tmp, fourpiG,
        #ifdef HAVE_HICLASS_BG
          H_spline, acc,
        #else
          cosmo,
        #endif
        0.5 * dtau / numsteps_ncdm[i]);
				f_params[0] = tmp;
				f_params[1] = tmp * tmp * sim.numpts;
		if (sim.gr_flag > 0)
			pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
		else
			pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
		#ifdef BENCHMARK
			moveParts_count++;
			moveParts_time += MPI_Wtime() - ref2_time;
			ref2_time = MPI_Wtime();
		#endif
				rungekutta4bg(tmp, fourpiG,
          #ifdef HAVE_HICLASS_BG
            H_spline, acc,
          #else
            cosmo,
          #endif
          0.5 * dtau / numsteps_ncdm[i]);
			}
		}

		// cdm and baryon particle update
		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;
		if (sim.gr_flag > 0)
		{
			maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
			if (sim.baryon_flag)
				maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
		}
		else
		{
			maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
			if (sim.baryon_flag)
				maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
		}

		#ifdef BENCHMARK
			update_q_count++;
			update_q_time += MPI_Wtime() - ref2_time;
			ref2_time = MPI_Wtime();
		#endif

		rungekutta4bg(a, fourpiG,
		#ifdef HAVE_HICLASS_BG
		H_spline, acc,
		#else
		cosmo,
		#endif
		0.5 * dtau);  // evolve background by half a time step

		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;
		if (sim.gr_flag > 0)
		{
			pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
			if (sim.baryon_flag)
				pcls_b.moveParticles(update_pos, dtau, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
		}
		else
		{
			pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
			if (sim.baryon_flag)
				pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
		}

		#ifdef BENCHMARK
			moveParts_count++;
			moveParts_time += MPI_Wtime() - ref2_time;
		#endif

		rungekutta4bg(a, fourpiG,
		#ifdef HAVE_HICLASS_BG
			H_spline, acc,
		#else
			cosmo,
		#endif
		0.5 * dtau);  // evolve background by half a time step

		parallel.max<double>(maxvel, numspecies);

		if (sim.gr_flag > 0)
		{
			for (i = 0; i < numspecies; i++)
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
		}
		// done particle update

		tau += dtau;

		if (sim.wallclocklimit > 0.)   // check for wallclock time limit
		{
			tmp = MPI_Wtime() - start_time;
			parallel.max(tmp);
			if (tmp > sim.wallclocklimit)   // hibernate
			{
				COUT << COLORTEXT_YELLOW << " reaching hibernation wallclock limit, hibernating..." << COLORTEXT_RESET << endl;
				COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
				if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
					plan_Bi.execute(FFT_BACKWARD);
		#ifdef CHECK_B
				if (sim.vector_flag == VECTOR_ELLIPTIC)
				{
					plan_Bi_check.execute(FFT_BACKWARD);
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi_check, a, tau, dtau, cycle);
				}
				else
		#endif
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi, a, tau, dtau, cycle);
				break;
			}
		}

		if (restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
			if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
				plan_Bi.execute(FFT_BACKWARD);
		#ifdef CHECK_B
			if (sim.vector_flag == VECTOR_ELLIPTIC)
			{
				plan_Bi_check.execute(FFT_BACKWARD);
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi, a, tau, dtau, cycle, restartcount);
			}
			else
		#endif
			hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi, a, tau, dtau, cycle, restartcount);
			restartcount++;
		}

		dtau_old = dtau;

    if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG,
    #ifdef HAVE_HICLASS_BG
      H_spline, acc
    #else
      cosmo
    #endif
    ))
		dtau = sim.Cf * dx;
	else
      dtau = sim.steplimit / Hconf(a, fourpiG,
      #ifdef HAVE_HICLASS_BG
        H_spline, acc
      #else
        cosmo
      #endif
      );

	cycle++;

	#ifdef BENCHMARK
		cycle_time += MPI_Wtime()-cycle_start_time;
	#endif
	}

		COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

	#ifdef BENCHMARK
		ref_time = MPI_Wtime();
	#endif

	#if defined(HAVE_CLASS) || defined(HAVE_HICLASS)
		if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
			freeCLASSstructures(class_background, class_thermo, class_perturbs);
	#endif

	#ifdef BENCHMARK
		lightcone_output_time += MPI_Wtime() - ref_time;
		run_time = MPI_Wtime() - start_time;

		parallel.sum(run_time);
		parallel.sum(cycle_time);
		parallel.sum(projection_time);
		parallel.sum(snapshot_output_time);
		parallel.sum(spectra_output_time);
		parallel.sum(lightcone_output_time);
		parallel.sum(gravity_solver_time);
		parallel.sum(kgb_update_time);
		parallel.sum(fft_time);
		parallel.sum(update_q_time);
		parallel.sum(moveParts_time);

		COUT << endl << "BENCHMARK" << endl;
		COUT << "total execution time  : "<<hourMinSec(run_time) << endl;
		COUT << "total number of cycles: "<< cycle << endl;
		COUT << "time consumption breakdown:" << endl;
		COUT << "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time <<"%."<<endl;
		COUT << "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time <<"%."<<endl;

		COUT << "----------- main loop: components -----------"<<endl;

		COUT << "projections                : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time <<"%."<<endl;
	// KGB update
	COUT << "kgb_update                : "<< hourMinSec(kgb_update_time) << " ; " << 100. * kgb_update_time/cycle_time <<"%."<<endl;
		COUT << "snapshot outputs           : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time <<"%."<<endl;
		COUT << "lightcone outputs          : "<< hourMinSec(lightcone_output_time) << " ; " << 100. * lightcone_output_time/cycle_time <<"%."<<endl;
		COUT << "power spectra outputs      : "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time <<"%."<<endl;
		COUT << "update momenta (count: "<<update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/cycle_time <<"%."<<endl;
		COUT << "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/cycle_time <<"%."<<endl;
		COUT << "gravity solver             : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time <<"%."<<endl;
		COUT << "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time <<"%."<<endl;
	#endif

	#ifdef EXTERNAL_IO
		ioserver.stop();
		}
	#endif
	return 0;
}
