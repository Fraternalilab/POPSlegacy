/*==============================================================================
arg.c : parse command line arguments
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "arg.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
	
	fprintf(stdout, "\n"
					" __   __   __   __   *\n"
					"|__) /  \\ |__) /__` \n"
					"|    \\__/ |    .__/\n"); 

	fprintf(stdout, "\n%s %s\n", PACKAGE, VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header()
{
    fprintf(stdout, "\nPOPS* : Parameter OPtimised Surface of proteins and nucleic acids\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stdout, "\nCopyright (C) 2002-2012 Franca Fraternali (program author)\n"
			"Copyright (C) 2008-2012 Jens Kleinjung (modular C code)\n"
			"Copyright (C) 2002 Kuang Lin and Valerie Hindie (translation to C)\n"
			"Copyright (C) 2002 Luigi Cavallo (parametrisation)\n"
			"POPS* is free software and comes with ABSOLUTELY NO WARRANTY.\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the COPYING file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
	fprintf(stdout, "\nPOPS* METHOD\n");
    fprintf(stdout, "Fraternali, F. and van Gunsteren, W.F.\n"
			"An efficient mean solvation force model for use in molecular"
			" dynamics simulations of proteins in aqueous solution.\n"
			"Journal of Molecular Biology 256 (1996) 939-948.\n");
    fprintf(stdout, "\nFraternali, F. and Cavallo, L.\n"
			"Parameter optimized surfaces (POPS): analysis of key interactions "
			"and conformational changes in the ribosome.\n"
			"Nucleic Acids Research 30 (2002) 2950-2960.\n");
	fprintf(stdout, "\nPOPS SERVER\n");
	fprintf(stdout, "Cavallo, L., Kleinjung, J. and Fraternali, F.\n"
			"POPS: A fast algorithm for solvent accessible surface areas "
			"at atomic and residue level.\n"
			"Nucleic Acids Research 31 (2003) 3364-3366.\n");
	fprintf(stdout, "\nPOPSCOMP SERVER\n");
	fprintf(stdout, "Kleinjung, J. and Fraternali, F.\n"
			"POPSCOMP: an automated interaction analysis "
			"of biomolecular complexes.\n"
			"Nucleic Acids Research 33 (2005) W342-W346.\n\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg, Argpdb *argpdb)
{
    arg->pdbInFileName = "";
	arg->trajInFileName = 0; /* trajectory input file */ 
	argpdb->coarse = 0; /* Calpha-only computation [0,1] */
	argpdb->hydrogens = 0; /* hydrogens [0,1] */
	arg->multiModel = 0; /* input with multiple models */
	arg->rProbe = 1.4; /* probe radius (in Angstrom) */
	arg->silent = 0; /* suppress stdout */
    arg->sasaOutFileName = "pops.out";
    arg->sasatrajOutFileName = "popstraj";
	arg->compositionOut = 0; /* output of molecule composition */
    arg->sigmaOutFileName = "sigma.out";
    arg->sigmatrajOutFileName = "sigma.out";
	arg->typeOut = 0; /* output of atom/residue types */
	arg->topologyOut = 0; /* output of molecule topology */
	arg->atomOut = 0; /* output of atom areas */
	arg->residueOut = 0; /* output of residue areas */
	arg->chainOut = 0; /* output of chain areas */
    arg->neighbourOutFileName = "neli.out";
	arg->neighbourOut = 0; /* output of neighbour list */
    arg->parameterOutFileName = "para.out";
	arg->parameterOut = 0; /* output of atom parameters (for benchmarking) */
	arg->noTotalOut = 0; /* suppress output of total area (for benchmarking) */
	arg->noHeaderOut = 0; /* suppress output headers (for benchmarking) */
	arg->padding = 0; /* add lines to pad the missing hydrogen atom lines 
						(for benchmarking) */
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg, Argpdb *argpdb)
{
	if (strlen(arg->pdbInFileName) == 0)
		Error("Invalid PDB file name");
	assert(argpdb->coarse == 0 || argpdb->coarse == 1);
    assert(argpdb->hydrogens == 0 || argpdb->hydrogens == 1);
	assert(arg->multiModel == 0 || arg->multiModel == 1);
	assert(arg->rProbe > 0);
	assert(strlen(arg->sasaOutFileName) > 0);
	assert(strlen(arg->sasatrajOutFileName) > 0);
	assert(strlen(arg->sigmaOutFileName) > 0);
	assert(strlen(arg->sigmatrajOutFileName) > 0);
	assert(arg->compositionOut == 0 || arg->compositionOut == 1);
	assert(arg->typeOut == 0 || arg->typeOut == 1);
	assert(arg->topologyOut == 0 || arg->topologyOut == 1);
	assert(arg->atomOut == 0 || arg->atomOut == 1);
	assert(arg->residueOut == 0 || arg->residueOut == 1);
	assert(arg->chainOut == 0 || arg->chainOut == 1);
	assert(arg->neighbourOut == 0 || arg->neighbourOut == 1);
	assert(arg->parameterOut == 0 || arg->parameterOut == 1);
	assert(arg->noTotalOut == 0 || arg->noTotalOut == 1);
	assert(arg->silent == 0 || arg->silent == 1);
	assert(arg->noHeaderOut == 0 || arg->noHeaderOut == 1);
	assert(arg->padding == 0 || arg->padding == 1);
}

/*____________________________________________________________________________*/
/** print command line arguments */
static void print_args(Arg *arg, Argpdb *argpdb)
{
    time_t now;
    time(&now);

    fprintf(stdout, "\ndate: %s", ctime(&now));
    fprintf(stdout, \
                    "pdb: %s\n"
                    "traj: %s\n"
                    "coarse: %d\n"
                    "multiModel: %d\n"
                    "rProbe: %f\n",
        arg->pdbInFileName, arg->trajInFileName, argpdb->coarse,
        0, arg->rProbe);
    fflush(stdout);
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg, Argpdb *argpdb)
{
	int c;
	const char usage[] = "\npops [--pdb ...] [OPTIONS ...]\n\
	 INPUT OPTIONS\n\
	   --pdb <PDB input>\t\t(mode: mandatory, type: char  , default: void)\n\
	   --traj <trajectory input>\t(mode: optional, type: char  , default: void)\n\
	 MODE OPTIONS\n\
	   --coarse\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --multiModel (disabled)\t(mode: optional , type: no_arg, default: off)\n\
	   --rProbe <probe radius [A]>\t(mode: optional , type: float , default: 1.4)\n\
	   --silent\t\t\t(mode: optional , type: no_arg, default: off)\n\
	 OUTPUT OPTIONS\n\
	   --popsOut <POPS output>\t(mode: optional , type: char  , default: pops.out)\n\
	   --popstrajOut <POPS output>\t(mode: optional , type: char  , default: popstraj.out)\n\
	   --sigmaOut <SFE output>\t(mode: optional , type: char  , default: sigma.out)\n\
	   --sigmatrajOut <SFE output>\t(mode: optional , type: char  , default: sigmatraj.out)\n\
	   --compositionOut\t\t(mode: optional , type: no_arg, default: off)\n\
	   --typeOut\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --topologyOut\t\t(mode: optional , type: no_arg, default: off)\n\
	   --atomOut\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --residueOut\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --chainOut\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --neighbourOut\t\t(mode: optional , type: no_arg, default: off)\n\
	   --parameterOut\t\t(mode: optional , type: no_arg, default: off)\n\
	   --noTotalOut\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --noHeaderOut\t\t(mode: optional , type: no_arg, default: off)\n\
	   --padding\t\t\t(mode: optional , type: no_arg, default: off)\n\
	 INFO OPTIONS\n\
	   --cite\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --version\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --help\n";

    if (argc < 2) {
		print_header();
        fprintf(stderr, "%s", usage);
		print_license();
        exit(0);
    }

    set_defaults(arg, argpdb);

    /** long option definition */
    static struct option long_options[] = {
        {"pdb", required_argument, 0, 1},
        {"traj", required_argument, 0, 2},
        {"coarse", no_argument, 0, 3},
        {"multiModel", no_argument, 0, 4},
        {"rProbe", required_argument, 0, 5},
        {"popsOut", required_argument, 0, 6},
        {"popstrajOut", required_argument, 0, 7},
        {"sigmaOut", required_argument, 0, 8},
        {"sigmatrajOut", required_argument, 0, 9},
        {"silent", no_argument, 0, 10},
        {"compositionOut", no_argument, 0, 11},
        {"typeOut", no_argument, 0, 12},
        {"topologyOut", no_argument, 0, 13},
        {"atomOut", no_argument, 0, 14},
        {"residueOut", no_argument, 0, 15},
        {"chainOut", no_argument, 0, 16},
        {"neighbourOut", no_argument, 0, 17},
        {"parameterOut", no_argument, 0, 18},
        {"noTotalOut", no_argument, 0, 19},
        {"noHeaderOut", no_argument, 0, 20},
        {"padding", no_argument, 0, 21},
        {"cite", no_argument, 0, 22},
        {"version", no_argument, 0, 23},
        {"help", no_argument, 0, 24},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3 4 5:6:7:8:9:10 11 12 13 14 15 16 17 18 19 20 21 22 23 24", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->pdbInFileName = optarg;
                break;
            case 2:
                arg->trajInFileName = optarg;
                break;
            case 3:
                argpdb->coarse = 1;
                break;
            case 4:
                arg->multiModel = 1;
                break;
            case 5:
                arg->rProbe = atof(optarg);
                break;
            case 6:
                arg->sasaOutFileName = optarg;
                break;
            case 7:
                arg->sasatrajOutFileName = optarg;
                break;
            case 8:
                arg->sigmaOutFileName = optarg;
                break;
            case 9:
                arg->sigmatrajOutFileName = optarg;
                break;
            case 10:
                arg->silent = 1;
                break;
            case 11:
                arg->compositionOut = 1;
                break;
            case 12:
                arg->typeOut = 1;
                break;
            case 13:
                arg->topologyOut = 1;
                break;
            case 14:
                arg->atomOut = 1;
                break;
            case 15:
                arg->residueOut = 1;
                break;
            case 16:
                arg->chainOut = 1;
                break;
            case 17:
                arg->neighbourOut = 1;
                break;
            case 18:
                arg->parameterOut = 1;
                break;
            case 19:
                arg->noTotalOut = 1;
                break;
            case 20:
                arg->noHeaderOut = 1;
                break;
            case 21:
                arg->padding = 1;
                break;
            case 22:
                print_citation();
                exit(0);
            case 23:
				print_version();
				print_license();
                exit(0);
            case 24:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(1);
        }
    }

	check_input(arg, argpdb);
    print_header();
    print_args(arg, argpdb);

    return 0;
}

