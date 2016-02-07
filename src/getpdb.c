/*==============================================================================
getpdb.c : routines for reading PDB structures
Copyright (C) 2004 Jens Kleinjung
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
/** read PDB file */

/**
http://www.wwpdb.org/docs.html
The ATOM record:

COLUMNS        DATA TYPE       FIELD         DEFINITION
________________________________________________________________________________
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X
                                               in Angstroms.
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
                                               in Angstroms.
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
                                               in Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier,
                                               left-justified.
77 - 78        LString(2)      element       Element symbol,
                                               right-justified.
79 - 80        LString(2)      charge        Charge on the atom.
*/

int read_pdb(FILE *pdbInFile, Str *str, int coarse, int hydrogens)
{
	unsigned int i, j;
	unsigned int k = 0;
	char line[80];
	char stopline[80] = "";
    int stopflag = 0;
	unsigned int allocated_atom = 64;
	unsigned int allocated_residue = 64;
	char atomName[] = "    ";
	char resbuf;
	int ca_p = 0;
	/* for HETATM entries */
	regex_t *regexPattern; /* regular atom patterns */
	/* allowed HETATM atom types (standard N,CA,C,O) and elements (any N,C,O,P,S) */
	const int nHetAtom = 9;
	char hetAtomPattern[9][32] = {{" N  "},{" CA "},{" C  "},{" O  "},{".{1}C[[:print:]]{1,3}"},{".{1}N[[:print:]]{1,3}"},{".{1}O[[:print:]]{1,3}"},{".{1}P[[:print:]]{1,3}"},{".{1}S[[:print:]]{1,3}"}};
	char hetAtomNewname[9][32] = {{" N  "},{" CA "},{" C  "},{" O  "},{" C_ "},{" N_ "},{" O_ "},{" P_ "},{" S_ "}};

	/*____________________________________________________________________________*/
	/* initialise/allocate memory for set of (64) selected (CA) atom entries */
	str->nAtom = 0;
	str->nAllAtom = 0;
	str->nResidue = 0;
	str->nAllResidue = 0;
	str->nChain = 0;
	
	str->atom = safe_malloc(allocated_atom * sizeof(Atom));
	str->atomMap = safe_malloc(allocated_atom * sizeof(int));

	/* allocate memory for sequence residues */
	str->sequence.res = safe_malloc(allocated_residue * sizeof(char));

	/* compile allowed HETATM element patterns */
	regexPattern = safe_malloc(nHetAtom * sizeof(regex_t));
	compile_patterns(regexPattern, &(hetAtomPattern[0]), nHetAtom);

	/*____________________________________________________________________________*/
    /* count the number of models */
    while(fgets(line, 80, pdbInFile) != 0) {
        if (strncmp(line, "MODEL ", 6) == 0) {
            if (stopflag == 0) {
                stopflag = 1;
                continue;
            } else {
                strcpy(stopline, line);
                break;
            }
        }
    }

    /* rewind the file handle to the start */
	if (fseek(pdbInFile, 0L, SEEK_SET) != 0) {
		/* handle repositioning error */
	}

	/*____________________________________________________________________________*/
	/* not all PDB data types are used in this program to save resources */
    while(fgets(line, 80, pdbInFile) != 0) {
		ca_p = 0; /* CA or P flag */

		/*____________________________________________________________________________*/
		/* check conditions to start assigning this entry */
		/* skip other models */
		if((strcmp(line, stopline) == 0) && (stopflag == 1))
			break;

		/* read only ATOM/HETATM records */
		if((strncmp(line, "ATOM  ", 6) != 0) && (strncmp(line, "HETATM", 6) != 0))
			continue;

        /* skip alternative locations except for location 'A' */ 
		if (line[16] != 32 && line[16] != 65) {
			/*fprintf(stderr, "Warning: Skipping atom %d in alternative location %c\n",
				atoi(&line[6]), line[16]);*/
			continue;
		}

		/*____________________________________________________________________________*/
		/* read this entry */
		/* atom number */
		str->atom[str->nAtom].atomNumber = atoi(&line[6]);

		/* atom name */
		for (i = 12, j = 0; i < 16; )
			str->atom[str->nAtom].atomName[j++] = line[i++];
		str->atom[str->nAtom].atomName[j] = '\0';

		/* alternative location */
		/*str->atom[str->nAtom].alternativeLocation[0] = line[16];	
		str->atom[str->nAtom].alternativeLocation[1] = '\0';*/

		/* residue name */
		for (i = 17, j = 0; i < 20; )
			str->atom[str->nAtom].residueName[j++] = line[i++];
		str->atom[str->nAtom].residueName[j] = '\0';

		/* chain identifier */
		str->atom[str->nAtom].chainIdentifier[0] = line[21];
		str->atom[str->nAtom].chainIdentifier[1] = '\0';

		/* residue number */
		str->atom[str->nAtom].residueNumber = atoi(&line[22]);

		/* code for insertion of residues */
		/*str->atom[str->nAtom].icode[0] = line[26];
		str->atom[str->nAtom].icode[1] = '\0';*/

		/* coordinates */
		str->atom[str->nAtom].pos.x = atof(&line[30]);
		str->atom[str->nAtom].pos.y = atof(&line[38]);
		str->atom[str->nAtom].pos.z = atof(&line[46]);

		/*printf("x %6.4f, y %6.4f, z %6.4f\n", str->atom[str->nAtom].x,
			str->atom[str->nAtom].y, str->atom[str->nAtom].z);*/

		/* occupancy */
		/*str->atom[str->nAtom].occupancy = atof(&line[54]);*/

		/* temperature factor */
		/*str->atom[str->nAtom].temp_f = atof(&line[60]);*/

		/* segment identifier */
		/*for (i = 72, j = 0; i < 76; )
			str->atom[str->nAtom].segmentIdentifier[j++] = line[i++];
		str->atom[str->nAtom].segmentIdentifier[j] = '\0';*/

		/* element */
		for (i = 76, j = 0; i < 78; )
			str->atom[str->nAtom].element[j++] = line[i++];
		str->atom[str->nAtom].element[j] = '\0';

		/* charge */
		/*for (i = 78, j = 0; i < 80; )
			str->atom[str->nAtom].charge[j++] = line[i++];
		str->atom[str->nAtom].charge[j] = '\0';*/

		/* description: everything before coordinates */
		for (i = 0, j = 0; i < 30; )
			str->atom[str->nAtom].description[j++] = line[i++];
		str->atom[str->nAtom].description[j] = '\0';

		/*____________________________________________________________________________*/
		/* check conditions to record this entry */
		/* if no hydrogens set, skip hydrogen lines */
		if (! hydrogens) {
			strip_char(str->atom[str->nAtom].atomName, &(atomName[0]));
			/* skip patterns 'H...' and '?H..', where '?' is a digit */
			if ((atomName[0] == 'H') || \
				((atomName[0] >= 48) && (atomName[0] <= 57) && (atomName[1] == 'H'))) {
				++ str->nAllAtom;
				continue;
			}
		}

		/* check whether ATOM residue name is standard */
		if (strncmp(line, "ATOM  ", 6) == 0)
			resbuf = aacode(str->atom[str->nAtom].residueName);
			
		/* detect CA and P atoms of standard residues for residue allocation */
		if ((strncmp(line, "ATOM  ", 6) == 0) &&
			((strncmp(str->atom[str->nAtom].atomName, " CA ", 4) == 0) ||
			(strncmp(str->atom[str->nAtom].atomName, " P  ", 4) == 0))) {
			str->sequence.res[k++] = aacode(str->atom[str->nAtom].residueName);
			++ ca_p;
			if (k == allocated_residue)
				str->sequence.res = safe_realloc(str->sequence.res, (allocated_residue += 64) * sizeof(char));
		}

		/* standardise non-standard atom names */
		standardise_name(str->atom[str->nAtom].residueName, str->atom[str->nAtom].atomName);

		/* process HETATM entries */
		if (strncmp(line, "HETATM", 6) == 0)
			if (process_het(str, &(line[0]), regexPattern, &(hetAtomNewname[0]), nHetAtom) != 0)
				continue;

		/* in coarse mode record only CA and P entries */
		if (!ca_p && coarse)
			continue;

		/*____________________________________________________________________________*/
		/* count number of allResidues (including HETATM residues) */
        if (str->nAtom == 0 || str->atom[str->nAtom].residueNumber != str->atom[str->nAtom - 1].residueNumber)
			++ str->nAllResidue;

		/*____________________________________________________________________________*/
		/* count number of chains */
        if (str->nAtom == 0 || str->atom[str->nAtom].chainIdentifier[0] != str->atom[str->nAtom - 1].chainIdentifier[0])
			++ str->nChain;

		/*____________________________________________________________________________*/
		/* recors original atom order (count) */
		str->atomMap[str->nAtom] = str->nAllAtom;
		/* increment to next atom entry */
		++ str->nAtom;
		++ str->nAllAtom;

		/*____________________________________________________________________________*/
		/* allocate more memory if needed */
		if (str->nAtom == allocated_atom) {
			allocated_atom += 64;
			str->atom = safe_realloc(str->atom, allocated_atom * sizeof(Atom));
			str->atomMap = safe_realloc(str->atomMap, allocated_atom * sizeof(int));
		}
	}
	str->sequence.res[k] = '\0';
	str->nResidue = k;

	/*____________________________________________________________________________*/
	/* free the compiled regular expressions */
	free_patterns(regexPattern, nHetAtom);
	free(regexPattern);

	/*____________________________________________________________________________*/
	return 0;
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
    read_pdb(pdbInFile, pdb, argpdb->coarse, argpdb->hydrogens);
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

