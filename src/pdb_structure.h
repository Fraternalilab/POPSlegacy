/*===============================================================================
pdb_structure.h : PDB structure definition
Copyright (C) 2004-2008 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef PDB_STRUCTURE_H
#define PDB_STRUCTURE_H

#include "seq.h"
#include "vector.h"

/*____________________________________________________________________________*/
/* structures */

/* atom : definition of PDB atom format, numbers indicate columns */
typedef struct
{
	char recordName[8]; /* Record type; 1 -  6*/
	int atomNumber; /* Atom serial number;  7 - 11 */
	char atomName[8]; /* Atom name; 13 - 16 */
	/*char alternativeLocation[2];*/ /* Alternate location indicator; 17 */
	char residueName[4]; /* Residue name; 18 - 20 */
	char chainIdentifier[2]; /* Chain identifier; 22 */
	int residueNumber; /* Residue sequence number; 23 - 26 */
	char icode[2]; /* Code for insertion of residues; 27 */
	Vec pos; /* position vector (x, y, z) */
	/*float occupancy;*/ /* Occupancy; 55 - 60 */
	/*float temp_f;*/ /* Temperature factor; 61 - 66 */
	/*char segmentIdentifier[5];*/ /* Segment identifier; 73 - 76 */
	char element[3]; /* Element symbol; 77 - 78 */
	/*char charge[3];*/ /* Charge on the atom; 79 - 80 */
	char description[32]; /* everything before coordinates */
	Vec tpos; /* transformed position vector */
	int atomType; /* GROMOS atom type */
	int groupID; /* atom group ID */
} Atom;

/* molecular structure */
typedef struct
{
	Atom *atom; /* array of selected (CA) atoms constituting structure */
	int *atomMap; /* map of the selected atom count to the original atom count */
	int nAtom; /* number of selected (CA) atoms */
	int nAllAtom; /* number of all atoms */
	int nResidue; /* number of residues (CA and P atoms) */
	int nAllResidue; /* number of all residues (including HETATM) */
	int nChain; /* number of chains */
	Seq sequence; /* amino acid sequence of structure */
	Seq strSequence; /* sequence of string-encoded structure */
} Str;

#endif
