/*==============================================================================
sigma_data.h : sigma parameters for Solvation Free Energy computation
Copyright (C) 2011 Franca Fraternali
Copyright (C) 2011 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef SIGMADATA_H
#define SIGMADATA_H

#include "sigma_const.h"

/** The array 'constant_sigma_data' defines the constant SIGMA values. 
 * The first array element (0) defines the SIGMA values of POPS at atomic resolution.
 * The second array element (1) defines the SIGMA values of POPS at residuic (CA) resolution
 * (or coarse grained).
 * Note that the number of atoms per residue needs to be modified if atoms are
 * added/removed to/from the arrays. */
ConstantSigma constant_sigma_data[] = {
	/*____________________________________________________________________________*/
	/* array element 0 : atomic resolution */
	/*____________________________________________________________________________*/
	{
		/*____________________________________________________________________________*/
		39, /* number of residue types */

		/*____________________________________________________________________________*/
		/* number of atoms per residue */
		{
			/* protein */
			/*ALA*/ 5,
			/*ARG*/ 11, 
			/*ASN*/ 8,
			/*ASP*/ 8,
			/*CYS*/ 6,
			/*GLN*/ 9,
			/*GLU*/ 9,
			/*GLY*/ 4,
			/*HIS*/ 10,
			/*ILE*/ 8,
			/*LEU*/ 8,
			/*LYS*/ 9,
			/*MET*/ 8,
			/*PHE*/ 11,
			/*PRO*/ 7,
			/*SER*/ 6,
			/*THR*/ 7,
			/*TRP*/ 14,
			/*TYR*/ 12,
			/*VAL*/ 7,
			/*ANY*/ 3,
			/*SOL*/ 1,
			/* RNA */
			/* A */ 23,
			/* C */ 21,
			/* G */ 24,
			/* I */ 23,
			/* N */ 22, 
			/* T */ 22, 
			/* U */ 21,
			/* DNA */
			/*DA */ 22,
			/*DC */ 20,
			/*DG */ 23,
			/*DI */ 22,
			/*DN */ 21, 
			/*DT */ 21, 
			/*DU */ 20,
			/* unknown polymer residue */
			/*UNK*/ 5,
			/* ligand */
			/*HET*/ 6,
			/* unknown ligand */
			/*UNL*/ 5
		},

		/*____________________________________________________________________________*/
		/* atom-specific data:
			1: residue name, 2: atom name, 3: sigma atom parameter, 4. sigma group parameter.

			residue name: PDB standard residue.
			atom name: PDB standard atom.
            SGIMA atom parameter: free energy per atom area (kJ / (mol * nm^2)).
            SGIMA group parameter: free energy per atom  area (kJ / (mol * nm^2)). */
		{
			{	/*N   ALA*/   {"ALA", "N",    0. ,   0. }, 
				/*CA  ALA*/   {"ALA", "CA",   3.8,   4.1},
				/*C   ALA*/   {"ALA", "C",    0. ,   0. },
				/*O   ALA*/   {"ALA", "O",   -7.2,  -7.3},
				/*CB  ALA*/   {"ALA", "CB",   3.3,   4.1}},
			{	/*N   ARG*/   {"ARG", "N",    0. ,   0. },
				/*CA  ARG*/   {"ARG", "CA",   3.8,   4.1},
				/*C   ARG*/   {"ARG", "C",    0. ,   0. },
				/*O   ARG*/   {"ARG", "O",   -7.2,  -7.3},
				/*CB  ARG*/   {"ARG", "CB",   5.0,   4.1},
				/*CG  ARG*/   {"ARG", "CG",   5.0,   4.1},
				/*CD  ARG*/   {"ARG", "CD",   5.0,   4.1},
				/*NE  ARG*/   {"ARG", "NE",   0. ,   0.},
				/*CZ  ARG*/   {"ARG", "CZ",   0. ,   0.},
				/*NH1 ARG*/   {"ARG", "NH1",-13.3, -23.3},
				/*NH2 ARG*/   {"ARG", "NH2",-13.3,  -7.3}},
			{	/*N   ASN*/   {"ASN", "N",    0. ,   0. },
				/*CA  ASN*/   {"ASN", "CA",   3.8,   4.1},
				/*C   ASN*/   {"ASN", "C",    0. ,   0. },
				/*O   ASN*/   {"ASN", "O",   -7.2,  -7.3},
				/*CB  ASN*/   {"ASN", "CB",   5.0,   4.1},
				/*CG  ASN*/   {"ASN", "CG",   0. ,   0. },
				/*OD1 ASN*/   {"ASN", "OD1", -7.2,  -7.3},
				/*ND2 ASN*/   {"ASN", "ND2", -4.0 , -7.3}},
			{	/*N   ASP*/   {"ASP", "N",    0. ,   0. },
				/*CA  ASP*/   {"ASP", "CA",   3.8,   4.1},
				/*C   ASP*/   {"ASP", "C",    0. ,   0. },
				/*O   ASP*/   {"ASP", "O",   -7.2,  -7.3},
				/*CB  ASP*/   {"ASP", "CB",   5.0,   4.1},
				/*CG  ASP*/   {"ASP", "CG",   0. ,   0. },
				/*OD1 ASP*/   {"ASP", "OD1",-21.7, -23.3},
				/*OD2 ASP*/   {"ASP", "OD2",-21.7, -23.3}},
			{	/*N   CYS*/   {"CYS", "N",    0. ,   0. },
				/*CA  CYS*/   {"CYS", "CA",   3.8,   4.1},
				/*C   CYS*/   {"CYS", "C",    0. ,   0. },
				/*O   CYS*/   {"CYS", "O",   -7.2,  -7.3},
				/*CB  CYS*/   {"CYS", "CB",   5.0,   4.1},
				/*SG  CYS*/   {"CYS", "SG",   0. ,   0. }},
			{	/*N   GLN*/   {"GLN", "N",    0. ,   0. },
				/*CA  GLN*/   {"GLN", "CA",   3.8,   4.1},
				/*C   GLN*/   {"GLN", "C",    0. ,   0. },
				/*O   GLN*/   {"GLN", "O",   -7.2,  -7.3},
				/*CB  GLN*/   {"GLN", "CB",   5.0,   4.1},
				/*CG  GLN*/   {"GLN", "CG",   5.0,   4.1},
				/*CD  GLN*/   {"GLN", "CD",   0. ,   0. },
				/*OE1 GLN*/   {"GLN", "OE1", -7.2,  -7.3},
				/*NE2 GLN*/   {"GLN", "NE2", -4.0,  -7.3}},
			{	/*N   GLU*/   {"GLU", "N",    0. ,   0. },
				/*CA  GLU*/   {"GLU", "CA",   3.8,   4.1},
				/*C   GLU*/   {"GLU", "C",    0. ,   0. },
				/*O   GLU*/   {"GLU", "O",   -7.2,  -7.3},
				/*CB  GLU*/   {"GLU", "CB",   5.0,   4.1},
				/*CG  GLU*/   {"GLU", "CG",   5.0,   4.1},
				/*CD  GLU*/   {"GLU", "CD",   0. ,   0. },
				/*OE1 GLU*/   {"GLU", "OE1",-21.7, -23.3},
				/*OE2 GLU*/   {"GLU", "OE2",-21.7, -23.3}},
			{	/*N   GLY*/   {"GLY", "N",    0. ,   0. },
				/*CA  GLY*/   {"GLY", "CA",   3.8,   4.1},
				/*C   GLY*/   {"GLY", "C",    0. ,   0. },
				/*O   GLY*/   {"GLY", "O",   -7.2,  -7.3}},
			{	/*N   HIS*/   {"HIS", "N",    0. ,   0. },
				/*CA  HIS*/   {"HIS", "CA",   3.8,   4.1},
				/*C   HIS*/   {"HIS", "C",    0. ,   0. },
				/*O   HIS*/   {"HIS", "O",   -7.2,  -7.3},
				/*CB  HIS*/   {"HIS", "CB",   5.0,   4.1},
				/*CG  HIS*/   {"HIS", "CG",   4.5,   4.1},
				/*ND1 HIS*/   {"HIS", "ND1", -4.5,  -7.3},
				/*CD2 HIS*/   {"HIS", "CD2",  4.5,   4.1},
				/*CE1 HIS*/   {"HIS", "CE1",  4.5,   4.1},
				/*NE2 HIS*/   {"HIS", "NE2", -4.5,  -7.3}},
			{	/*N   ILE*/   {"ILE", "N",    0. ,   0. },
				/*CA  ILE*/   {"ILE", "CA",   3.8,   4.1},
				/*C   ILE*/   {"ILE", "C",    0. ,   0. },
				/*O   ILE*/   {"ILE", "O",   -7.2,  -7.3},
				/*CB  ILE*/   {"ILE", "CB",   5.0,   4.1},
				/*CG1 ILE*/   {"ILE", "CG1",  5.0,   4.1},
				/*CG2 ILE*/   {"ILE", "CG2",  5.0,   4.1},
				/*CD1 ILE*/   {"ILE", "CD1",  3.3,   4.1}},
			{	/*N   LEU*/   {"LEU", "N",    0. ,   0. },
				/*CA  LEU*/   {"LEU", "CA",   3.8,   4.1},
				/*C   LEU*/   {"LEU", "C",    0. ,   0. },
				/*O   LEU*/   {"LEU", "O",   -7.2,  -7.3},
				/*CB  LEU*/   {"LEU", "CB",   5.0,   4.1},
				/*CG  LEU*/   {"LEU", "CG",   5.0,   4.1},
				/*CD1 LEU*/   {"LEU", "CD1",  3.3,   4.1},
				/*CD2 LEU*/   {"LEU", "CD2",  3.3,   4.1}},
			{	/*N   LYS*/   {"LYS", "N",    0. ,   0. },
				/*CA  LYS*/   {"LYS", "CA",   3.8,   4.1},
				/*C   LYS*/   {"LYS", "C",    0. ,   0. },
				/*O   LYS*/   {"LYS", "O",   -7.2,  -7.3},
				/*CB  LYS*/   {"LYS", "CB",   5.0,   4.1},
				/*CG  LYS*/   {"LYS", "CG",   5.0,   4.1},
				/*CD  LYS*/   {"LYS", "CD",   5.0,   4.1},
				/*CE  LYS*/   {"LYS", "CE",   5.0,   4.1},
				/*NZ  LYS*/   {"LYS", "NZ", -26.1, -23.3}},
			{	/*N   MET*/   {"MET", "N",    0. ,   0. },
				/*CA  MET*/   {"MET", "CA",   3.8,   4.1},
				/*C   MET*/   {"MET", "C",    0. ,   0. },
				/*O   MET*/   {"MET", "O",   -7.2,  -7.3},
				/*CB  MET*/   {"MET", "CB",   5.0,   4.1},
				/*CG  MET*/   {"MET", "CG",   5.0,   4.1},
				/*SD  MET*/   {"MET", "SD",   0. ,   0. },
				/*CE  MET*/   {"MET", "CE",   3.3,   4.1}},
			{	/*N   PHE*/   {"PHE", "N",    0. ,   0. },
				/*CA  PHE*/   {"PHE", "CA",   3.8,   4.1},
				/*C   PHE*/   {"PHE", "C",    0. ,   0. },
				/*O   PHE*/   {"PHE", "O",   -7.2,  -7.3},
				/*CB  PHE*/   {"PHE", "CB",   5.0,   4.1},
				/*CG  PHE*/   {"PHE", "CG",   4.5,   4.1},
				/*CD1 PHE*/   {"PHE", "CD1",  4.5,   4.1},
				/*CD2 PHE*/   {"PHE", "CD2",  4.5,   4.1},
				/*CE1 PHE*/   {"PHE", "CE1",  4.5,   4.1},
				/*CE2 PHE*/   {"PHE", "CE2",  4.5,   4.1},
				/*CZ  PHE*/   {"PHE", "CZ",   4.5,   4.1}},
			{	/*N   PRO*/   {"PRO", "N",    0. ,   0. },
				/*CA  PRO*/   {"PRO", "CA",   3.8,   4.1},
				/*C   PRO*/   {"PRO", "C",    0. ,   0. },
				/*O   PRO*/   {"PRO", "O",   -7.2,  -7.3},
				/*CB  PRO*/   {"PRO", "CB",   5.0,   4.1},
				/*CG  PRO*/   {"PRO", "CG",   5.0,   4.1},
				/*CD  PRO*/   {"PRO", "CD",   5.0,   4.1}},
			{	/*N   SER*/   {"SER", "N",    0. ,   0. },
				/*CA  SER*/   {"SER", "CA",   3.8,   4.1},
				/*C   SER*/   {"SER", "C",    0. ,   0. },
				/*O   SER*/   {"SER", "O",   -7.2,  -7.3},
				/*CB  SER*/   {"SER", "CB",   5.0,   4.1},
				/*OG  SER*/   {"SER", "OG",  -7.0,  -7.3}},
			  {	/*N   THR*/   {"THR", "N",    0. ,   0. },
				/*CA  THR*/   {"THR", "CA",   3.8,   4.1},
				/*C   THR*/   {"THR", "C",    0. ,   0. },
				/*O   THR*/   {"THR", "O",   -7.2,  -7.3},
				/*CB  THR*/   {"THR", "CB",   5.0,   4.1},
				/*OG1 THR*/   {"THR", "OG1", -7.0,  -7.3},
				/*CG2 THR*/   {"THR", "CG2",  5.0,   4.1}},
			{	/*N   TRP*/   {"TRP", "N",    0. ,   0. },
				/*CA  TRP*/   {"TRP", "CA",   3.8,   4.1},
				/*C   TRP*/   {"TRP", "C",    0. ,   0. },
				/*O   TRP*/   {"TRP", "O",   -7.2,  -7.3},
				/*CB  TRP*/   {"TRP", "CB",   5.0,   4.1},
				/*CG  TRP*/   {"TRP", "CG",   4.5,   4.1},
				/*CD1 TRP*/   {"TRP", "CD1",  4.5,   4.1},
				/*CD2 TRP*/   {"TRP", "CD2",  4.5,   4.1},
				/*NE1 TRP*/   {"TRP", "NE1", -4.5,  -7.3},
				/*CE2 TRP*/   {"TRP", "CE2",  4.5,   4.1},
				/*CE3 TRP*/   {"TRP", "CE3",  4.5,   4.1},
				/*CH2 TRP*/   {"TRP", "CH2",  4.5,   4.1},
				/*CZ2 TRP*/   {"TRP", "CZ2",  4.5,   4.1},
				/*CZ3 TRP*/   {"TRP", "CZ3",  4.5,   4.1}},
			{	/*N   TYR*/   {"TYR", "N",    0. ,   0. },
				/*CA  TYR*/   {"TYR", "CA",   3.8,   4.1},
				/*C   TYR*/   {"TYR", "C",    0. ,   0. },
				/*O   TYR*/   {"TYR", "O",   -7.2,  -7.3},
				/*CB  TYR*/   {"TYR", "CB",   5.0,   4.1},
				/*CG  TYR*/   {"TYR", "CG",   4.5,   4.1},
				/*CD1 TYR*/   {"TYR", "CD1",  4.5,   4.1},
				/*CD2 TYR*/   {"TYR", "CD2",  4.5,   4.1},
				/*CE1 TYR*/   {"TYR", "CE1",  4.5,   4.1},
				/*CE2 TYR*/   {"TYR", "CE2",  4.5,   4.1},
				/*CZ  TYR*/   {"TYR", "CZ",   4.5,   4.1},
				/*OH  TYR*/   {"TYR", "OH",  -7.0,  -7.3}},
			{	/*N   VAL*/   {"VAL", "N",    0. ,   0. },
				/*CA  VAL*/   {"VAL", "CA",   3.8,   4.1},
				/*C   VAL*/   {"VAL", "C",    0. ,   0. },
				/*O   VAL*/   {"VAL", "O",   -7.2,  -7.3},
				/*CB  VAL*/   {"VAL", "CB",   5.0,   4.1},
				/*CG1 VAL*/   {"VAL", "CG1",  3.3,   4.1},
				/*CG2 VAL*/   {"VAL", "CG2",  3.3,   4.1}},
			{   /*OXT ANY*/   {"ANY", "OXT",-21.7, -23.3},
				/*O1  ANY*/   {"ANY", "O1" ,-21.7, -23.3},
				/*O2  ANY*/   {"ANY", "O2", -21.7, -23.3}},
			{	/*O2  SOL*/   {"SOL", "O2",   0.,    0. }},
			/* sigma values of nucleotides have not been parametrised yet */
			{	/*P     A*/   {"A",   "P",    0., 0.},
				/*OP1   A*/   {"A",   "OP1",  0., 0.},
				/*OP2   A*/   {"A",   "OP2",  0., 0.},
				/*OP3   A*/   {"A",   "OP3",  0., 0.},
				/*O5'   A*/   {"A",   "O5'",  0., 0.},
				/*C5'   A*/   {"A",   "C5'",  0., 0.},
				/*C4'   A*/   {"A",   "C4'",  0., 0.},
				/*O4'   A*/   {"A",   "O4'",  0., 0.},
				/*C3'   A*/   {"A",   "C3'",  0., 0.},
				/*O3'   A*/   {"A",   "O3'",  0., 0.},
				/*C2'   A*/   {"A",   "C2'",  0., 0.},
				/*O2'   A*/   {"A",   "O2'",  0., 0.},
				/*C1'   A*/   {"A",   "C1'",  0., 0.},
				/*N9    A*/   {"A",   "N9",   0., 0.},
				/*C8    A*/   {"A",   "C8",   0., 0.},
				/*N7    A*/   {"A",   "N7",   0., 0.},
				/*C5    A*/   {"A",   "C5",   0., 0.},
				/*C6    A*/   {"A",   "C6",   0., 0.},
				/*N6    A*/   {"A",   "N6",   0., 0.},
				/*N1    A*/   {"A",   "N1",   0., 0.},
				/*C2    A*/   {"A",   "C2",   0., 0.},
				/*N3    A*/   {"A",   "N3",   0., 0.},
				/*C4    A*/   {"A",   "C4",   0., 0.}},
			{	/*P     C*/   {"C",   "P",    0., 0.},
				/*OP1   C*/   {"C",   "OP1",  0., 0.},
				/*OP2   C*/   {"C",   "OP2",  0., 0.},
				/*OP3   C*/   {"C",   "OP3",  0., 0.},
				/*O5'   C*/   {"C",   "O5'",  0., 0.},
				/*C5'   C*/   {"C",   "C5'",  0., 0.},
				/*C4'   C*/   {"C",   "C4'",  0., 0.},
				/*O4'   C*/   {"C",   "O4'",  0., 0.},
				/*C3'   C*/   {"C",   "C3'",  0., 0.},
				/*O3'   C*/   {"C",   "O3'",  0., 0.},
				/*C2'   C*/   {"C",   "C2'",  0., 0.},
				/*O2'   C*/   {"C",   "O2'",  0., 0.},
				/*C1'   C*/   {"C",   "C1'",  0., 0.},
				/*N1    C*/   {"C",   "N1",   0., 0.},
				/*C2    C*/   {"C",   "C2",   0., 0.},
				/*O2    C*/   {"C",   "O2",   0., 0.},
				/*N3    C*/   {"C",   "N3",   0., 0.},
				/*C4    C*/   {"C",   "C4",   0., 0.},
				/*N4    C*/   {"C",   "N4",   0., 0.},
				/*C5    C*/   {"C",   "C5",   0., 0.},
				/*C6    C*/   {"C",   "C6",   0., 0.}},
			{	/*P     G*/   {"G",   "P",    0., 0.},
				/*OP1   G*/   {"G",   "OP1",  0., 0.},
				/*OP2   G*/   {"G",   "OP2",  0., 0.},
				/*OP3   G*/   {"G",   "OP3",  0., 0.},
				/*O5'   G*/   {"G",   "O5'",  0., 0.},
				/*C5'   G*/   {"G",   "C5'",  0., 0.},
				/*C4'   G*/   {"G",   "C4'",  0., 0.},
				/*O4'   G*/   {"G",   "O4'",  0., 0.},
				/*C3'   G*/   {"G",   "C3'",  0., 0.},
				/*O3'   G*/   {"G",   "O3'",  0., 0.},
				/*C2'   G*/   {"G",   "C2'",  0., 0.},
				/*O2'   G*/   {"G",   "O2'",  0., 0.},
				/*C1'   G*/   {"G",   "C1'",  0., 0.},
				/*N9    G*/   {"G",   "N9",   0., 0.},
				/*C8    G*/   {"G",   "C8",   0., 0.},
				/*N7    G*/   {"G",   "N7",   0., 0.},
				/*C5    G*/   {"G",   "C5",   0., 0.},
				/*C6    G*/   {"G",   "C6",   0., 0.},
				/*O6    G*/   {"G",   "O6",   0., 0.},
				/*N1    G*/   {"G",   "N1",   0., 0.},
				/*C2    G*/   {"G",   "C2",   0., 0.},
				/*N2    G*/   {"G",   "N2",   0., 0.},
				/*N3    G*/   {"G",   "N3",   0., 0.},
				/*C4    G*/   {"G",   "C4",   0., 0.}},
			{	/*P     I*/   {"I",   "P",    0., 0.},
				/*OP1   I*/   {"I",   "OP1",  0., 0.},
				/*OP2   I*/   {"I",   "OP2",  0., 0.},
				/*OP3   I*/   {"I",   "OP3",  0., 0.},
				/*O5'   I*/   {"I",   "O5'",  0., 0.},
				/*C5'   I*/   {"I",   "C5'",  0., 0.},
				/*C4'   I*/   {"I",   "C4'",  0., 0.},
				/*O4'   I*/   {"I",   "O4'",  0., 0.},
				/*C3'   I*/   {"I",   "C3'",  0., 0.},
				/*O3'   I*/   {"I",   "O3'",  0., 0.},
				/*C2'   I*/   {"I",   "C2'",  0., 0.},
				/*O2'   I*/   {"I",   "O2'",  0., 0.},
				/*C1'   I*/   {"I",   "C1'",  0., 0.},
				/*N9    I*/   {"I",   "N9",   0., 0.},
				/*C8    I*/   {"I",   "C8",   0., 0.},
				/*N7    I*/   {"I",   "N7",   0., 0.},
				/*C5    I*/   {"I",   "C5",   0., 0.},
				/*C6    I*/   {"I",   "C6",   0., 0.},
				/*O6    I*/   {"I",   "O6",   0., 0.},
				/*N1    I*/   {"I",   "N1",   0., 0.},
				/*C2    I*/   {"I",   "C2",   0., 0.},
				/*N3    I*/   {"I",   "N3",   0., 0.},
				/*C4    I*/   {"I",   "C4",   0., 0.}},
			{	/*P     N*/	  {"N",   "P",    0., 0.},
				/*OP1   N*/	  {"N",   "OP1",  0., 0.},
				/*OP2	N*/	  {"N",   "OP2",  0., 0.},
				/*OP3	N*/	  {"N",   "OP3",  0., 0.},
				/*O5'	N*/	  {"N",   "O5'",  0., 0.},
				/*C5'	N*/	  {"N",   "C5'",  0., 0.},
				/*C4'	N*/	  {"N",   "C4'",  0., 0.},
				/*O4'	N*/	  {"N",   "O4'",  0., 0.},
				/*C3'	N*/	  {"N",   "C3'",  0., 0.},
				/*O3'	N*/	  {"N",   "O3'",  0., 0.},
				/*C2'	N*/	  {"N",   "C2'",  0., 0.},
				/*O2'	N*/	  {"N",   "O2'",  0., 0.},
				/*C1'	N*/	  {"N",   "C1'",  0., 0.},
				/*N1 	N*/	  {"N",   "N1",   0., 0.},
				/*C2 	N*/	  {"N",   "C2",   0., 0.},
				/*O2 	N*/	  {"N",   "O2",   0., 0.},
				/*N3 	N*/	  {"N",   "N3",   0., 0.},
				/*C4 	N*/	  {"N",   "C4",   0., 0.},
				/*N4 	N*/	  {"N",   "N4",   0., 0.},
				/*O5 	N*/	  {"N",   "O4",   0., 0.},
				/*C5 	N*/	  {"N",   "C5",   0., 0.},
				/*C6 	N*/	  {"N",   "C6",   0., 0.}},
			{	/*P     T*/   {"T",   "P",    0., 0.},
				/*OP1   T*/   {"T",   "OP1",  0., 0.},
				/*OP2   T*/   {"T",   "OP2",  0., 0.},
				/*OP3   T*/   {"T",   "OP3",  0., 0.},
				/*O5'   T*/   {"T",   "O5'",  0., 0.},
				/*C5'   T*/   {"T",   "C5'",  0., 0.},
				/*C4'   T*/   {"T",   "C4'",  0., 0.},
				/*O4'   T*/   {"T",   "O4'",  0., 0.},
				/*C3'   T*/   {"T",   "C3'",  0., 0.},
				/*O3'   T*/   {"T",   "O3'",  0., 0.},
				/*C2'   T*/   {"T",   "C2'",  0., 0.},
				/*C1'   T*/   {"T",   "C1'",  0., 0.},
				/*N1    T*/   {"T",   "N1",   0., 0.},
				/*C2    T*/   {"T",   "C2",   0., 0.},
				/*O2    T*/   {"T",   "O2",   0., 0.},
				/*N3    T*/   {"T",   "N3",   0., 0.},
				/*C4    T*/   {"T",   "C4",   0., 0.},
				/*O4    T*/   {"T",   "O4",   0., 0.},
				/*C5    T*/   {"T",   "C5",   0., 0.},
				/*C7    T*/   {"T",   "C7",   0., 0.},
				/*C6    T*/   {"T",   "C6",   0., 0.}},
			{	/*P, 	U*/	  {"U",   "P",    0., 0.},
				/*OP1	U*/	  {"U",   "OP1",  0., 0.},
				/*OP2	U*/	  {"U",   "OP2",  0., 0.},
				/*OP3	U*/	  {"U",   "OP3",  0., 0.},
				/*O5'	U*/	  {"U",   "O5'",  0., 0.},
				/*C5'	U*/	  {"U",   "C5'",  0., 0.},
				/*C4'	U*/	  {"U",   "C4'",  0., 0.},
				/*O4'	U*/	  {"U",   "O4'",  0., 0.},
				/*C3'	U*/	  {"U",   "C3'",  0., 0.},
				/*O3'	U*/	  {"U",   "O3'",  0., 0.},
				/*C2'	U*/	  {"U",   "C2'",  0., 0.},
				/*O2'	U*/	  {"U",   "O2'",  0., 0.},
				/*C1'	U*/	  {"U",   "C1'",  0., 0.},
				/*N1 	U*/	  {"U",   "N1",   0., 0.},
				/*C2 	U*/	  {"U",   "C2",   0., 0.},
				/*O2 	U*/	  {"U",   "O2",   0., 0.},
				/*N3 	U*/	  {"U",   "N3",   0., 0.},
				/*C4 	U*/	  {"U",   "C4",   0., 0.},
				/*O4 	U*/	  {"U",   "O4",   0., 0.},
				/*C5 	U*/	  {"U",   "C5",   0., 0.},
				/*C6 	U*/	  {"U",   "C6",   0., 0.}},
			{	/*P    DA*/   {"DA",  "P",    0., 0.},
				/*OP1  DA*/   {"DA",  "OP1",  0., 0.},
				/*OP2  DA*/   {"DA",  "OP2",  0., 0.},
				/*OP3  DA*/   {"DA",  "OP3",  0., 0.},
				/*O5'  DA*/   {"DA",  "O5'",  0., 0.},
				/*C5'  DA*/   {"DA",  "C5'",  0., 0.},
				/*C4'  DA*/   {"DA",  "C4'",  0., 0.},
				/*O4'  DA*/   {"DA",  "O4'",  0., 0.},
				/*C3'  DA*/   {"DA",  "C3'",  0., 0.},
				/*O3'  DA*/   {"DA",  "O3'",  0., 0.},
				/*C2'  DA*/   {"DA",  "C2'",  0., 0.},
				/*C1'  DA*/   {"DA",  "C1'",  0., 0.},
				/*N9   DA*/   {"DA",  "N9",   0., 0.},
				/*C8   DA*/   {"DA",  "C8",   0., 0.},
				/*N7   DA*/   {"DA",  "N7",   0., 0.},
				/*C5   DA*/   {"DA",  "C5",   0., 0.},
				/*C6   DA*/   {"DA",  "C6",   0., 0.},
				/*N6   DA*/   {"DA",  "N6",   0., 0.},
				/*N1   DA*/   {"DA",  "N1",   0., 0.},
				/*C2   DA*/   {"DA",  "C2",   0., 0.},
				/*N3   DA*/   {"DA",  "N3",   0., 0.},
				/*C4   DA*/   {"DA",  "C4",   0., 0.}},
			{	/*P    DC*/   {"DC",  "P",    0., 0.},
				/*OP1  DC*/   {"DC",  "OP1",  0., 0.},
				/*OP2  DC*/   {"DC",  "OP2",  0., 0.},
				/*OP3  DC*/   {"DC",  "OP3",  0., 0.},
				/*O5'  DC*/   {"DC",  "O5'",  0., 0.},
				/*C5'  DC*/   {"DC",  "C5'",  0., 0.},
				/*C4'  DC*/   {"DC",  "C4'",  0., 0.},
				/*O4'  DC*/   {"DC",  "O4'",  0., 0.},
				/*C3'  DC*/   {"DC",  "C3'",  0., 0.},
				/*O3'  DC*/   {"DC",  "O3'",  0., 0.},
				/*C2'  DC*/   {"DC",  "C2'",  0., 0.},
				/*C1'  DC*/   {"DC",  "C1'",  0., 0.},
				/*N1   DC*/   {"DC",  "N1",   0., 0.},
				/*C2   DC*/   {"DC",  "C2",   0., 0.},
				/*O2   DC*/   {"DC",  "O2",   0., 0.},
				/*N3   DC*/   {"DC",  "N3",   0., 0.},
				/*C4   DC*/   {"DC",  "C4",   0., 0.},
				/*N4   DC*/   {"DC",  "N4",   0., 0.},
				/*C5   DC*/   {"DC",  "C5",   0., 0.},
				/*C6   DC*/   {"DC",  "C6",   0., 0.}},
			{	/*P    DG*/   {"DG",  "P",    0., 0.},
				/*OP1  DG*/   {"DG",  "OP1",  0., 0.},
				/*OP2  DG*/   {"DG",  "OP2",  0., 0.},
				/*OP3  DG*/   {"DG",  "OP3",  0., 0.},
				/*O5'  DG*/   {"DG",  "O5'",  0., 0.},
				/*C5'  DG*/   {"DG",  "C5'",  0., 0.},
				/*C4'  DG*/   {"DG",  "C4'",  0., 0.},
				/*O4'  DG*/   {"DG",  "O4'",  0., 0.},
				/*C3'  DG*/   {"DG",  "C3'",  0., 0.},
				/*O3'  DG*/   {"DG",  "O3'",  0., 0.},
				/*C2'  DG*/   {"DG",  "C2'",  0., 0.},
				/*C1'  DG*/   {"DG",  "C1'",  0., 0.},
				/*N9   DG*/   {"DG",  "N9",   0., 0.},
				/*C8   DG*/   {"DG",  "C8",   0., 0.},
				/*N7   DG*/   {"DG",  "N7",   0., 0.},
				/*C5   DG*/   {"DG",  "C5",   0., 0.},
				/*C6   DG*/   {"DG",  "C6",   0., 0.},
				/*O6   DG*/   {"DG",  "O6",   0., 0.},
				/*N1   DG*/   {"DG",  "N1",   0., 0.},
				/*C2   DG*/   {"DG",  "C2",   0., 0.},
				/*N2   DG*/   {"DG",  "N2",   0., 0.},
				/*N3   DG*/   {"DG",  "N3",   0., 0.},
				/*C4   DG*/   {"DG",  "C4",   0., 0.}},
			{	/*P    DI*/   {"DI",  "P",    0., 0.},
				/*OP1  DI*/   {"DI",  "OP1",  0., 0.},
				/*OP2  DI*/   {"DI",  "OP2",  0., 0.},
				/*OP3  DI*/   {"DI",  "OP3",  0., 0.},
				/*O5'  DI*/   {"DI",  "O5'",  0., 0.},
				/*C5'  DI*/   {"DI",  "C5'",  0., 0.},
				/*C4'  DI*/   {"DI",  "C4'",  0., 0.},
				/*O4'  DI*/   {"DI",  "O4'",  0., 0.},
				/*C3'  DI*/   {"DI",  "C3'",  0., 0.},
				/*O3'  DI*/   {"DI",  "O3'",  0., 0.},
				/*C2'  DI*/   {"DI",  "C2'",  0., 0.},
				/*C1'  DI*/   {"DI",  "C1'",  0., 0.},
				/*N9   DI*/   {"DI",  "N9",   0., 0.},
				/*C8   DI*/   {"DI",  "C8",   0., 0.},
				/*N7   DI*/   {"DI",  "N7",   0., 0.},
				/*C5   DI*/   {"DI",  "C5",   0., 0.},
				/*C6   DI*/   {"DI",  "C6",   0., 0.},
				/*O6   DI*/   {"DI",  "O6",   0., 0.},
				/*N1   DI*/   {"DI",  "N1",   0., 0.},
				/*C2   DI*/   {"DI",  "C2",   0., 0.},
				/*N3   DI*/   {"DI",  "N3",   0., 0.},
				/*C4   DI*/   {"DI",  "C4",   0., 0.}},
			{	/*P    DN*/	  {"DN",  "P",    0., 0.},
				/*OP1  DN*/	  {"DN",  "OP1",  0., 0.},
				/*OP2  DN*/	  {"DN",  "OP2",  0., 0.},
				/*OP3  DN*/	  {"DN",  "OP3",  0., 0.},
				/*O5'  DN*/	  {"DN",  "O5'",  0., 0.},
				/*C5'  DN*/	  {"DN",  "C5'",  0., 0.},
				/*C4'  DN*/	  {"DN",  "C4'",  0., 0.},
				/*O4'  DN*/	  {"DN",  "O4'",  0., 0.},
				/*C3'  DN*/	  {"DN",  "C3'",  0., 0.},
				/*O3'  DN*/	  {"DN",  "O3'",  0., 0.},
				/*C2'  DN*/	  {"DN",  "C2'",  0., 0.},
				/*C1'  DN*/	  {"DN",  "C1'",  0., 0.},
				/*N1   DN*/	  {"DN",  "N1",   0., 0.},
				/*C2   DN*/	  {"DN",  "C2",   0., 0.},
				/*O2   DN*/	  {"DN",  "O2",   0., 0.},
				/*N3   DN*/	  {"DN",  "N3",   0., 0.},
				/*C4   DN*/	  {"DN",  "C4",   0., 0.},
				/*N4   DN*/	  {"DN",  "N4",   0., 0.},
				/*O5   DN*/	  {"DN",  "O4",   0., 0.},
				/*C5   DN*/	  {"DN",  "C5",   0., 0.},
				/*C6   DN*/	  {"DN",  "C6",   0., 0.}},
			{	/*P    DT*/   {"DT",  "P",    0., 0.},
				/*OP1  DT*/   {"DT",  "OP1",  0., 0.},
				/*OP2  DT*/   {"DT",  "OP2",  0., 0.},
				/*OP3  DT*/   {"DT",  "OP3",  0., 0.},
				/*O5'  DT*/   {"DT",  "O5'",  0., 0.},
				/*C5'  DT*/   {"DT",  "C5'",  0., 0.},
				/*C4'  DT*/   {"DT",  "C4'",  0., 0.},
				/*O4'  DT*/   {"DT",  "O4'",  0., 0.},
				/*C3'  DT*/   {"DT",  "C3'",  0., 0.},
				/*O3'  DT*/   {"DT",  "O3'",  0., 0.},
				/*C2'  DT*/   {"DT",  "C2'",  0., 0.},
				/*C1'  DT*/   {"DT",  "C1'",  0., 0.},
				/*N1   DT*/   {"DT",  "N1",   0., 0.},
				/*C2   DT*/   {"DT",  "C2",   0., 0.},
				/*O2   DT*/   {"DT",  "O2",   0., 0.},
				/*N3   DT*/   {"DT",  "N3",   0., 0.},
				/*C4   DT*/   {"DT",  "C4",   0., 0.},
				/*O4   DT*/   {"DT",  "O4",   0., 0.},
				/*C5   DT*/   {"DT",  "C5",   0., 0.},
				/*C7   DT*/   {"DT",  "C7",   0., 0.},
				/*C6   DT*/   {"DT",  "C6",   0., 0.}},
			{	/*P,   DU*/	  {"DU",  "P",    0., 0.},
				/*OP1  DU*/	  {"DU",  "OP1",  0., 0.},
				/*OP2  DU*/	  {"DU",  "OP2",  0., 0.},
				/*OP3  DU*/	  {"DU",  "OP3",  0., 0.},
				/*O5'  DU*/	  {"DU",  "O5'",  0., 0.},
				/*C5'  DU*/	  {"DU",  "C5'",  0., 0.},
				/*C4'  DU*/	  {"DU",  "C4'",  0., 0.},
				/*O4'  DU*/	  {"DU",  "O4'",  0., 0.},
				/*C3'  DU*/	  {"DU",  "C3'",  0., 0.},
				/*O3'  DU*/	  {"DU",  "O3'",  0., 0.},
				/*C2'  DU*/	  {"DU",  "C2'",  0., 0.},
				/*C1'  DU*/	  {"DU",  "C1'",  0., 0.},
				/*N1   DU*/	  {"DU",  "N1",   0., 0.},
				/*C2   DU*/	  {"DU",  "C2",   0., 0.},
				/*O2   DU*/	  {"DU",  "O2",   0., 0.},
				/*N3   DU*/	  {"DU",  "N3",   0., 0.},
				/*C4   DU*/	  {"DU",  "C4",   0., 0.},
				/*O4   DU*/	  {"DU",  "O4",   0., 0.},
				/*C5   DU*/	  {"DU",  "C5",   0., 0.},
				/*C6   DU*/	  {"DU",  "C6",   0., 0.}},
			{	/*C   UNK*/   {"UNK", "C",    0., 0.},
				/*N   UNK*/   {"UNK", "N",    0., 0.},
				/*O   UNK*/   {"UNK", "O",    0., 0.},
				/*P   UNK*/   {"UNK", "P",    0., 0.},
				/*S   UNK*/   {"UNK", "S",    0., 0.}},
			{	/*N_  HET*/   {"HET", "N_",   0., 0.}, 
				/*CA  HET*/   {"HET", "CA",   0., 0.},
				/*C_  HET*/   {"HET", "C_",   0., 0.},
				/*O_  HET*/   {"HET", "O_",   0., 0.},
				/*P_  HET*/   {"HET", "P_",   0., 0.},
				/*S_  HET*/   {"HET", "S_",   0., 0.}},
			{	/*C   UNL*/   {"UNL", "C",    0., 0.},
				/*N   UNL*/   {"UNL", "N",    0., 0.},
				/*O   UNL*/   {"UNL", "O",    0., 0.},
				/*P   UNL*/   {"UNL", "P",    0., 0.},
				/*S   UNL*/   {"UNL", "S",    0., 0.}}
		},
	},

	/*____________________________________________________________________________*/
	/* array element 1 : residue resolution (coarse grained) */
	/* Note that the name 'atom' in the coarse grained parameters means 'CA' or 'P',
		which represents the sphere of a residue. That means, CA atom and residue
		are synonymous in this context. Residue AND atom specification are kept here
		to use the same data structure as for the atomic parametrisation above. */
	/*____________________________________________________________________________*/
	{
		/*____________________________________________________________________________*/
		39, /* number of residue types */

		/*____________________________________________________________________________*/
		/* number of atoms per residue */
		{
			/* protein */
			/*ALA*/ 1,
			/*ARG*/ 1,
			/*ASN*/ 1,
			/*ASP*/ 1,
			/*CYS*/ 1,
			/*GLN*/ 1,
			/*GLU*/ 1,
			/*GLY*/ 1,
			/*HIS*/ 1,
			/*ILE*/ 1,
			/*LEU*/ 1,
			/*LYS*/ 1,
			/*MET*/ 1,
			/*PHE*/ 1,
			/*PRO*/ 1,
			/*SER*/ 1,
			/*THR*/ 1,
			/*TRP*/ 1,
			/*TYR*/ 1,
			/*VAL*/ 1,
			/*ANY*/ 1,
			/*SOL*/ 1,
			/* RNA */
			/* A */ 1,
			/* G */ 1,
			/* I */ 1,
			/* N */ 1,
			/* T */ 1,
			/* U */ 1,
			/* DNA */
			/*DA */ 1,
			/*DC */ 1,
			/*DG */ 1,
			/*DI */ 1,
			/*DN */ 1, 
			/*DT */ 1, 
			/*DU */ 1,
			/* unknown polymer residue */
			/*UNK*/ 1,
			/* ligand */
			/*HET*/ 1,
			/* unknown ligand */
			/*UNL*/ 1
		},

		/*____________________________________________________________________________*/
		/* 1: residue name, 2: atom name, 3: sigma atom parameter, 4. sigma group parameter. */
		{
			{	/*CA	ALA*/ {"ALA", "CA",  0., 0.}},
			{	/*CA	ARG*/ {"ARG", "CA",  0., 0.}}, 
			{	/*CA	ASN*/ {"ASN", "CA",  0., 0.}}, 
			{	/*CA	ASP*/ {"ASP", "CA",  0., 0.}}, 
			{	/*CA	CYS*/ {"CYS", "CA",  0., 0.}}, 
			{	/*CA	GLN*/ {"GLN", "CA",  0., 0.}}, 
			{	/*CA	GLU*/ {"GLU", "CA",  0., 0.}}, 
			{	/*CA	GLY*/ {"GLY", "CA",  0., 0.}}, 
			{	/*CA	HIS*/ {"HIS", "CA",  0., 0.}}, 
			{	/*CA	ILE*/ {"ILE", "CA",  0., 0.}}, 
			{	/*CA	LEU*/ {"LEU", "CA",  0., 0.}}, 
			{	/*CA	LYS*/ {"LYS", "CA",  0., 0.}}, 
			{	/*CA	MET*/ {"MET", "CA",  0., 0.}}, 
			{	/*CA	PHE*/ {"PHE", "CA",  0., 0.}}, 
			{	/*CA	PRO*/ {"PRO", "CA",  0., 0.}}, 
			{	/*CA	SER*/ {"SER", "CA",  0., 0.}}, 
			{	/*CA	THR*/ {"THR", "CA",  0., 0.}}, 
			{	/*CA	TRP*/ {"TRP", "CA",  0., 0.}}, 
			{	/*CA	TYR*/ {"TYR", "CA",  0., 0.}}, 
			{	/*CA	VAL*/ {"VAL", "CA",  0., 0.}}, 
			{	/*OXT	ANY*/ {"ANY", "OXT", 0., 0.}},
			{	/*O2	SOL*/ {"SOL", "O2",  0., 0.}},
			{	/*P		  A*/ {"A",   "P",   0., 0.}},
			{	/*P		  C*/ {"C",   "P",   0., 0.}},
			{	/*P		  G*/ {"G",   "P",   0., 0.}},
			{	/*P		  I*/ {"I",   "P",   0., 0.}},
			{	/*P		  N*/ {"N",   "P",   0., 0.}},
			{	/*P		  T*/ {"T",   "P",   0., 0.}},
			{	/*P		  U*/ {"U",   "P",   0., 0.}},
			{	/*P		 DA*/ {"DA",  "P",   0., 0.}},
			{	/*P		 DC*/ {"DC",  "P",   0., 0.}},
			{	/*P		 DG*/ {"DG",  "P",   0., 0.}},
			{	/*P		 DI*/ {"DI",  "P",   0., 0.}},
			{	/*P		 DN*/ {"DN",  "P",   0., 0.}},
			{	/*P		 DT*/ {"DT",  "P",   0., 0.}},
			{	/*P		 DU*/ {"DU",  "P",   0., 0.}},
			{	/*CA	UNK*/ {"UNK", "CA",  0., 0.}},
			{	/*CA	HET*/ {"HET", "CA",  0., 0.}},
			{	/*CA	UNL*/ {"UNL", "CA",  0., 0.}}
		}
	}
};

#endif

