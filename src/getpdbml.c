/*==============================================================================
getpdbml.c : routines for reading PDBML structures
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/
#include "getpdbml.h"
#include "pdb_structure.h"

/*____________________________________________________________________________*/
int parseXML(const char *filename, Str *pdb) {
    xmlDoc *doc; /* the resulting document tree */
    xmlNode *root_node = 0;
	xmlNode *cur_node = 0;
	xmlNode *site_node = 0;
	xmlNode *atom_node = 0;
	unsigned int allocated_atom = 64;
	xmlChar *content = 0;

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
	pdb->nAtom = 0;

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

			}

			/* increment atom number */
			++ pdb->nAtom;
			/* allocate more memory if needed */
			if (pdb->nAtom == allocated_atom) {
				allocated_atom += 64;
				pdb->atom = safe_realloc(pdb->atom, allocated_atom * sizeof(Atom));
			}
		}
	}

	/*____________________________________*/
	/* write document to file */
    //xmlSaveFormatFileEnc("outfile.xml", doc, "UTF-8", 1);

	/*____________________________________*/
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

