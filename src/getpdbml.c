/*==============================================================================
getpdbml.c : routines for reading PDBML structures
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/
#include "getpdbml.h"
#include "pdb_structure.h"

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
int parseXML(const char *filename, Str *pdb) {
    xmlDoc *doc; /* the resulting document tree */
    xmlNode *root_node = 0;
	xmlNode *cur_node = 0;
	xmlNode *site_node = 0;
	xmlNode *atom_node = 0;
	unsigned int allocated_atom = 64;
	unsigned int allocated_residue = 64;
	xmlChar *content = 0;
	unsigned int k = 0;
	int ca_p = 0;

	/*____________________________________________________________________________*/
    /*parse the file and get the document (DOM) */
	if ((doc = xmlReadFile(filename, NULL, 0)) == NULL) {
        fprintf(stderr, "XML Parser: Failed to read %s\n", filename);
		exit(-1);
	}

	/*____________________________________________________________________________*/
	/* parse document tree */
	/* set root node */
	root_node = xmlDocGetRootElement(doc);

	/* atom sites are set as child nodes */
	for (cur_node = root_node->children; cur_node; cur_node = cur_node->next) {
		if(strcmp("atom_siteCategory", (char *)cur_node->name) == 0) {
			site_node = cur_node;
			break;
		}
	}

	/*____________________________________________________________________________*/
	/* allocate PDB structure */
	pdb->atom = safe_malloc(allocated_atom * sizeof(Atom));
	/* array of residue-centric atom indices */
	pdb->resAtom = safe_malloc(allocated_residue * sizeof(int));
	/* allocate memory for sequence residues */
	pdb->sequence.res = safe_malloc(allocated_residue * sizeof(char));

	pdb->nAtom = 0;
	pdb->nAllAtom = 0;

	/*____________________________________________________________________________*/
	/* traverse XML tree (atom sites) */
	for (atom_node = site_node->children; atom_node; atom_node = atom_node->next) {
		if (strcmp("atom_site", (char *)atom_node->name) == 0) {
			/* children of atom this atom site */
			for (cur_node = atom_node->children; cur_node; cur_node = cur_node->next) {
				/* assign node content to string */
				content = xmlNodeGetContent(cur_node);

				/* copy string content to PDB data structure */
				/* temperature factor */
				if (strcmp((char *)cur_node->name, "B_iso_or_equiv") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].temperatureFactor));
				}
				/* x coordinate */
				if (strcmp((char *)cur_node->name, "Cartn_x") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].pos.x));
				}
				/* y coordinate */
				if (strcmp((char *)cur_node->name, "Cartn_y") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].pos.y));
				}
				/* z coordinate */
				if (strcmp((char *)cur_node->name, "Cartn_z") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].pos.z));
				}
				/* chain identifier */
				if (strcmp((char *)cur_node->name, "auth_asym_id") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].chainIdentifier);
				}
				/* atom name */
				if (strcmp((char *)cur_node->name, "auth_atom_id") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].atomName);
				}
				/* residue name */
				if (strcmp((char *)cur_node->name, "auth_comp_id") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].residueName);
				}
				/* occupancy */
				if (strcmp((char *)cur_node->name, "occupancy") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].occupancy));
				}
				/* model number */
				if (strcmp((char *)cur_node->name, "pdbx_PDB_model_num") == 0) {
					sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].modelNumber));
				}
				/* atom element */
				if (strcmp((char *)cur_node->name, "type_symbol") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].element);
				}
				/* charge */
				if (strcmp((char *)cur_node->name, "pdbx_formal_charge") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].charge);
					sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].formalCharge));
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].partialCharge));
				}
			}

			/*____________________________________________________________________________*/
			/* select entries to record */
			/* only first MODEL if several are persent in PDB entry */
			if (pdb->atom[pdb->nAtom].modelNumber > 1) {
				goto ENDPARSE;
			}

			/* skip hydrogen atoms */
			if (strcmp(pdb->atom[pdb->nAtom].element, "H") == 0) {
				continue;
			}

			/* detect CA and P atoms of standard residues for residue allocation */
			if ((strcmp(pdb->atom[pdb->nAtom].atomName, "CA") == 0) ||
				(strcmp(pdb->atom[pdb->nAtom].atomName, "P") == 0)) {
				pdb->resAtom[k] = pdb->nAtom;
				pdb->sequence.res[k ++] = aacode(pdb->atom[pdb->nAtom].residueName);
				if (k == allocated_residue) {
					allocated_residue += 64;
					pdb->resAtom = safe_realloc(pdb->resAtom, allocated_residue * sizeof(int));
					pdb->sequence.res = safe_realloc(pdb->sequence.res, allocated_residue * sizeof(char));
				}
				++ ca_p;
			}

			/* standardise non-standard atom names (here GRO ILE_CD) */
			standardise_name(pdb->atom[pdb->nAtom].residueName, pdb->atom[pdb->nAtom].atomName);

			/* increment atom number */
			++ pdb->nAtom;
			++ pdb->nAllAtom;

			/* allocate more memory if needed */
			if (pdb->nAtom == allocated_atom) {
				allocated_atom += 64;
				pdb->atom = safe_realloc(pdb->atom, allocated_atom * sizeof(Atom));
			}
		}
	}

	ENDPARSE:
		fprintf(stderr, "\tUsing only the first model of the PDB entry\n");

	/*____________________________________________________________________________*/
	/* free global variables */
    xmlFreeDoc(doc);

	return 0;
}

/*____________________________________________________________________________*/
/* parse PDB file in XML format */
void read_structure_xml(Arg *arg, Argpdb *argpdb, Str *pdb)
{
	const char filename[] = "5lff.xml";

    LIBXML_TEST_VERSION;

	fprintf(stderr, "Parsing XML input file\n");
    parseXML(filename, pdb);

    xmlCleanupParser();
    xmlMemoryDump();
}

