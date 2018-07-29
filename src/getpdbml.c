/*==============================================================================
getpdbml.c : routines for reading PDBML structures
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/
#include "getpdbml.h"
#include "pdb_structure.h"

/*____________________________________________________________________________*/
int parseXML(const char *filename) {
    xmlDocPtr doc; /* the resulting document tree */
    xmlNodePtr root_node = NULL;
	xmlNodePtr node = NULL;
	xmlNodePtr node1 = NULL; /* node pointers */
    char buff[256];
    int i, j;

	/*____________________________________*/
	/* read document from file */
    doc = xmlReadFile(filename, NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Failed to parse %s\n", filename);
		return 0;
    }

	/*____________________________________*/
	/* create tree */
    root_node = xmlNewNode(NULL, BAD_CAST "root");
    xmlDocSetRootElement(doc, root_node);

	/* DTD declaration, not mandatory */
	xmlCreateIntSubset(doc, BAD_CAST "root", NULL, BAD_CAST "tree2.dtd");

     /* xmlNewChild() creates a new node,
		"attached" as child node of root_node node */
    xmlNewChild(root_node, NULL, BAD_CAST "node1", BAD_CAST "content of node 1");
	/* same as above, but the new child node doesn't have a content */
    xmlNewChild(root_node, NULL, BAD_CAST "node2", NULL);

	/* loop to repeat node creation */
    for (i = 5; i < 7; ++ i) {
        sprintf(buff, "node %d", i);
        node = xmlNewChild(root_node, NULL, BAD_CAST buff, NULL);
        for (j = 1; j < 4; ++ j) {
            sprintf(buff, "node %d %d", i, j);
            node1 = xmlNewChild(node, NULL, BAD_CAST buff, NULL);
            xmlNewProp(node1, BAD_CAST "odd", BAD_CAST((j % 2) ? "no" : "yes"));
        }
    }

	/*____________________________________*/
	/* write document to file */
    xmlSaveFormatFileEnc("outfile.xml", doc, "UTF-8", 1);

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

    parseXML(filename);

    xmlCleanupParser();
    xmlMemoryDump();
}

