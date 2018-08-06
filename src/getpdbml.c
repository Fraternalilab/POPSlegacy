/*==============================================================================
getpdbml.c : routines for reading PDBML structures
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/
#include "getpdbml.h"
#include "pdb_structure.h"

/*____________________________________________________________________________*/
int parseXML(const char *filename) {
    xmlDoc *doc; /* the resulting document tree */
    xmlNode *root_node = NULL;
	xmlNode *cur_node = NULL;
	xmlNode *site_node = NULL;
	xmlNode *atom_node = NULL;

	//char buff[256];

	/*____________________________________*/
    /*parse the file and get the document (DOM) */
	if ((doc = xmlReadFile(filename, NULL, 0)) == NULL) {
        fprintf(stderr, "XML Parser: Failed to read %s\n", filename);
		exit(-1);
	}

	/*____________________________________*/
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

	/*____________________________________*/
	/* traverse XML tree (atom sites) and copy content to PDB data structure */
	for (atom_node = site_node->children; atom_node; atom_node = atom_node->next) {
		if (strcmp("atom_site", (char *)atom_node->name) == 0) {
			;
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
    parseXML(filename);

    xmlCleanupParser();
    xmlMemoryDump();
}

