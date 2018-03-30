/*===============================================================================
getmmCIF.h : read mmCIF structures 
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETMMCIF_H
#define GETMMCIF_H

#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arg.h"
#include "argpdb.h"
#include "error.h"
#include "modstring.h"
#include "pdb_structure.h"
#include "pattern.h"
#include "safe.h"
#include "seq.h"
#include "vector.h"

/*____________________________________________________________________________*/
/* prototypes */
int read_mmCIF(FILE *pdbfile, Str *str, int coarse, int hydrogens);
void read_structure(Arg *arg, Argpdb *argpdb, Str *pdb);

#endif
