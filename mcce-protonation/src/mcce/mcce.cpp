#include <stdio.h>
#include <cstring>
#include <cmath>
#include "mcce.h"

void welcome();
int set_env(char *pdb, char *param);

int main(int argc, char *argv[])
{
	/* Welcome */
	welcome();

	/* Check user input */
	if (argc != 3) {
		fprintf(stderr, "Error: exactly two arguments have to be provided\n");
		// tpl_direcory is the path the upper level of tpl folder, exclude "param04",
		//mcce appends it automatically.
		fprintf(stderr, "    Usage: mcce-protonation input.pdb tpl_directory\n");
		return USERERR;
	}

	set_env(argv[1], argv[2]);


	/* Do step 0, initialization */
	db_open();
//	printf("Step 0. Initialize enviroment\n"); fflush(stdout);
	if (init()) {
		db_close();
		printf("Help message: double check file \"run.prm\" in current directory.\n");
		return USERERR;
	}
//	else printf("Step 0 Done.\n\n");


	/* Do step 1, premcce */
	if (env.do_premcce) {
//		printf("Step 1. Test and format structral file\n"); fflush(stdout);
		if (premcce()) {db_close(); return USERERR;}
//		else printf("Step 1 Done.\n\n");
	}
//	else printf("Not doing \"Step 1. Test and format structral file\"\n\n");


	db_close();
	return 0;
}

void welcome()
{
	printf("===========================================================\n");
	printf("<<< MCCE Multi-Conformation Continuum Electrostatics >>>   \n");
	printf(" Marilyn Gunner's Lab at City College of New York, 2005    \n");
	printf("-----------------------------------------------------------\n");
	printf("Version:        2.5                                        \n");
	printf("MCCE Home Page: http://www.sci.ccny.cuny.edu/~gunner/mcce  \n");
	printf("Support:        gunner@sci.ccny.cuny.edu                   \n");
	printf("Developed by:   Junjun Mao, Yifan Song, Marilyn Gunner     \n");
	printf("Reference MCCE: If you publish data calculated with MCCE,  \n");
	printf("                you need to cite papers suggested in MCCE  \n");
	printf("                Home Page.                                 \n");
	printf("===========================================================\n\n");
//	printf("Last Updates:                                              \n");
	//   printf("   05/10, Step2, use relative energy to decide vdw clash   \n");
//	printf("   04/10, Added Pascal Comte's GA/SA algorithm for sidechain\n");
//	printf("          packing and sampling.     			      \n");
//	printf("   04/10, Added APBS electrostatic calculation option      \n");
//	printf("   08/08, Total pairwise is set to be 999.0 if vdw is 999.0\n");
//	printf("   08/02, step 3, zero interaction bug fixed               \n");
//	printf("   07/23, step 1 changes the policy on AltLoc backbone atoms\n");
//	printf("          now it is warning instead of fatal error.        \n");
//	printf("   06/22, step 1 recognizes FME and ACE as NTR caps        \n");
//	printf("   05/27, step 2 reduces rotation steps for surface residues\n");
//	printf("          The run.prm.default uses 0.25 SAS to decide suface\n");
//	printf("          residues.                                        \n");
//	printf("   05/04, Program doesn't stop at duplicate parameters,    \n");
//	printf("          it issues a warning message. This fixes the deadlock\n");
//	printf("          caused by multiple new cofactors                 \n");
//	printf("   04/24, Energy lookup table compressed, size is reduced  \n");
//	printf("          to about 1/12                                    \n");
//	printf("   04/24, Opp files are compressed. Temporary files are    \n");
//	printf("          stored under local directory.                    \n");
//	printf("   04/14, Step 2 rotamer pruning by pairwise corrected.    \n");
//	printf("   04/08, Step 2 repacking now ses better H bnd correction.\n");
//	printf("   01/26, Step 3 fixed incorrect errno message.            \n");
//	printf("   01/14, Step 3 dielectric boundary condition improved.   \n");
//	printf("   12/07, Step 3 dielectric boundary condition is corrected.\n");
//	printf("   11/17, pKa fitting writes \"pKa titration curve too sharp\"\n");
//	printf("          when occ jumps from 0.015 to 0.985.              \n");
//	printf("   11/17, Exposed conformer is marked as E in history.     \n");
//	printf("   11/02, Optional quick step 3 included.                  \n");
//	printf("   11/02, Pruning fuction is added at the end of step 2.   \n");
//	printf("          Parameter lines were added to run.prm.           \n");
//	printf("   10/29, Temporary gdbm files have file names.            \n");
//	printf("   10/11, Use biggest rxn energy from the last 3 delphi    \n");
//	printf("          focusing runs.                                   \n");
//	printf("===========================================================\n\n");
	fflush(stdout);

	return;
}


int set_env(char *pdb, char *param)
{
	memset(&env, 0, sizeof(ENV));

	/* Default values */
	env.minimize_size = 0;
	env.PI                = 4.*atan(1.);
	env.d2r               = env.PI/180.;

	strcpy(env.debug_log,    "mcce_debug.log");
	strcpy(env.new_tpl,      "new.tpl");
	strcpy(env.progress_log, "progress.log");
	strcpy(env.extra, "extra.tpl");
	env.reassign     = 0;
	env.pbe_start = 0;
	env.pbe_end   = 999999;
	strcpy(env.pbe_solver, "delphi");
	strcpy(env.rxn_method, "self");		/*surface or self energies*/
	env.rot_specif   = 0;
	env.prune_thr = 0.01;
	env.ngh_vdw_thr  = 0.1;
	env.repack_e_thr_exposed  = 0.5;
	env.repack_e_thr_buried  = 4.;
	env.repack_fav_vdw_off   = 0;
	env.nconf_limit       =    0;
	env.n_hv_conf_limit   =   20;
	env.relax_wat         =    1;
	env.trans_dist        =  0.5;

	env.hdirected         =    0;
	env.hdirdiff          =  1.0;
	env.hdirlimt          =   36;

	env.water_relax_thr   =  2.4;

	env.n_initial_relax     =  0;
	//env.initial_relax_rebuild = 0;

	env.hv_relax_ncycle     =  0;
	env.hv_relax_niter      = 50;
	env.hv_relax_vdw_thr    =  5;
	env.hv_relax_hv_vdw_thr =  5;
	env.hv_relax_elec_thr   =  -2.0;
	env.hv_relax_elec_crg_thr =  0.1;
	env.hv_relax_elec_dist_thr = 2.4;
	env.hv_relax_dt         =  1;
	env.hv_tors_scale       =  1;
	env.hv_relax_constraint =  1.;
	env.hv_relax_constraint_frc = 10.;
	env.hv_relax_n_shake    =  3000;
	env.hv_relax_shake_tol =  1e-4;  /* Ratio to constraint distance */
	env.hv_relax_include_ngh    =  0;
	env.hv_relax_ngh_thr    =  4.;
	env.prune_rmsd        = 2.0;
	env.prune_ele         = 2.0;
	env.prune_vdw         = 2.0;


	env.relax_n_hyd       =    6;
	env.relax_clash_thr   =  10.;

	env.recalc_tors     = 0;

	env.default_radius = 1.7;
	env.factor_14lj = 0.5;
	env.epsilon_coulomb = 6.;

	env.sas2vdw = -0.06;

	env.warn_pairwise     = 20.0;
	env.big_pairwise      = 5.0;

	env.monte_adv_opt     =    0;
	env.anneal_temp_start = ROOMT;
	env.anneal_nstep      =    1;
	env.monte_tsx         =    0;
	env.anneal_niter_step =   30;
	env.monte_niter_max   =   -1;
	env.adding_conf       =    0;
	env.monte_old_input   =    0;
	env.monte_niter_chk   =  100;
	env.monte_converge    = 1e-4;
	env.monte_do_energy   =    0;
	env.monte_print_nonzero =  1;
	strcpy(env.pbe_folder, "/tmp");
	env.delphi_clean      =  1;
	env.ionrad = 0.0;
	env.salt =0.00;


	// Only run the first step of MCCE.
	env.do_premcce = 1;

	// path of the pdb file
	strcpy(env.inpdb, pdb);

	// path of the param folder, which should be in the same folder as delphi and mcce.
	strcpy(env.param, param);
	if (env.param[strlen(env.param)-1] == '/') env.param[strlen(env.param)-1] = '\0';
	strcat(env.param, "/param04");

	// Rule to rename pdb file
	strcpy(env.rename_rules, param);
	if (env.rename_rules[strlen(env.rename_rules)-1] == '/') env.rename_rules[strlen(env.rename_rules)-1] = '\0';
	strcat(env.rename_rules, "/name.txt");

	// Don't label terminal
	env.terminals = 0;

	return 0;
}
