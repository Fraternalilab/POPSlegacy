/*==============================================================================
getmmCIF.c : routines for reading mmCIF structures
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "getpdb.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/* match PDB residue name against constant residue name array */
__inline__ static char scan_array(char *code3, char *residue_array[], int shift)
{
	unsigned int i;
	char residue = ' ';

	for (i = 0; i < 26; ++ i)
		if (strncmp(code3, residue_array[i], 3) == 0) {
			residue = i + shift; /* shift=65 for UPPER, shift=97 for lower */
			break;
		}

	return residue;
}

/*____________________________________________________________________________*/
/** amino acid 3-letter to 1-letter code conversion */
__inline__ static char aacode(char *code3)
{
	char residue = ' '; /* 1-letter residue name */

	/* three-letter code of amino acid residues, exception HET -> X */
	char *aa3[] = {"ALA","---","CYS","ASP","GLU","PHE","GLY","HIS","ILE","---","LYS","LEU","MET","ASN","---","PRO","GLN","ARG","SER","THR","UNL","VAL","TRP","HET","TYR","UNK"};
	/* nucleotide residues */
	char *nuc[] = {"  A"," DA","  C"," DC","---","---","  G"," DG","  I"," DI","---","---","---","  N"," DN","---","---","---"," DT","  T","  U"," DU","---","---","---","---"};

	/* match against amino acid residues */
	residue = scan_array(code3, aa3, 65);

	/* match against nucleotide residues */
	if (residue == ' ')
		residue = scan_array(code3, nuc, 97);

	/* residue not found */
	if (residue == ' ') {
		Warning("Check the format of the PDB file: might be non-standard. See examples in the README file.");
		ErrorSpec("Unknown standard residue type in protein structure:", code3);
	}
	/*else
		fprintf(stderr, "%s:%d: %c ", __FILE__, __LINE__, residue);*/
	
	return residue;
}

/*____________________________________________________________________________*/
/** standardise non-standard atom names */
__inline__ static int standardise_name(char *residueName, char *atomName)
{
	/* GRO 'ILE CD' to PDB 'ILE CD1' */
	if ((strcmp(residueName, "ILE") == 0) && (strcmp(atomName, " CD ") == 0))
		strcpy(atomName, " CD1");

	return 0;
}

/*____________________________________________________________________________*/
/** process HET residues */
/* Differences between standard residues and HET residues:
 * 1. HET residues may be single atoms/molecules or part of a polymer;
 * 2. HET residues may not have a CA or P atom. */
__inline__ static int process_het(Str *str, char *line, regex_t *regexPattern, char (*hetAtomNewname)[32], int nHetAtom)
{
	int hetAtomNr = -1;

	/* residue name */
	strncpy(str->atom[str->nAtom].residueName, "HET", 3);

	/* atom name: assign only allowed atom elements, otherwise atom is skipped */
	if ((hetAtomNr = match_patterns(regexPattern, nHetAtom, &(str->atom[str->nAtom].atomName[0]))) >= 0) {
			sprintf(str->atom[str->nAtom].atomName, "%s", &(hetAtomNewname[hetAtomNr][0]));
	} else {
		WarningSpec("Skipping HEATM", str->atom[str->nAtom].atomName);
		return 1;
	}

	return 0;
}

/*____________________________________________________________________________*/
/** read mmCIF file */
int read_mmCIF(FILE *pdbInFile, Str *str, int coarse, int hydrogens)
{
}

/*____________________________________________________________________________*/
/** read CONECT record */
int read_conect(FILE *pdbInFile)
{
	char line[80];

	if (fseek(pdbInFile, 0L, SEEK_SET) != 0) {
		fprintf(stderr, "File pointer repositioning error\n");
		exit(1);
	}
	
    while(fgets(line, 80, pdbInFile) != 0) {
        if (strncmp(line, "CONECT", 6) == 0) {
			;
		}
	}

	return 0;
}

/*_____________________________________________________________________________*/
/** read PDB structure */
void read_structure(Arg *arg, Argpdb *argpdb, Str *pdb)
{
	FILE *pdbInFile = 0;

    pdbInFile = safe_open(arg->pdbInFileName, "r");
    pdb->sequence.name = safe_malloc((strlen(basename(arg->pdbInFileName)) + 1) * sizeof(char));
    strcpy(pdb->sequence.name, basename(arg->pdbInFileName));
    read_mmCIF(pdbInFile, pdb, argpdb->coarse, argpdb->hydrogens);
    fclose(pdbInFile);

    /* check for empty pdb structure and exit */
    if (pdb->nAtom == 0)
    {
        ErrorSpecNoexit("Invalid PDB file", arg->pdbInFileName);
        free(pdb->atom);
        free(pdb->sequence.res);
        free(pdb->sequence.name);
        exit(1);
    }

    if (! arg->silent) fprintf(stdout, "\tPDB file: %s\n"
										"\tPDB file content:\n"
										"\tnAtom = %d (excluding hydrogens)\n"
										"\tnAllAtom = %d (all atoms to match trajectory entries)\n",
							arg->pdbInFileName, pdb->nAtom, pdb->nAllAtom);
}

