/*==============================================================================
arg.h : parse command line arguments
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "argpdb.h"
#include "error.h"

/*____________________________________________________________________________*/
/* structures */

/* variables for commmand line arguments */
typedef struct  
{
    FILE *pdbInFile;
	char *pdbInFileName;
    FILE *trajInFile;
	char *trajInFileName;
    float rProbe;
	int silent;
	int multiModel;
    FILE *sasaOutFile;
    char *sasaOutFileName;
    FILE *sasatrajOutFile;
    char *sasatrajOutFileName;
    FILE *sigmaOutFile;
    char *sigmaOutFileName;
    FILE *sigmatrajOutFile;
    char *sigmatrajOutFileName;
	int compositionOut;
	int topologyOut;
	int typeOut;
	int atomOut;
	int residueOut;
	int chainOut;
    FILE *neighbourOutFile;
    char *neighbourOutFileName;
	int neighbourOut;
    FILE *parameterOutFile;
    char *parameterOutFileName;
	int parameterOut;
	int noTotalOut;
	int noHeaderOut;
	int padding;
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg, Argpdb *argpdb);

#endif
