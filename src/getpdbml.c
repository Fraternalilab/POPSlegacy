/*==============================================================================
getpdbml.c : routines for reading PDBML structures
	Using to a large extent the ReadPDB.c routine of Bioplib as template,
		partially verbatim and partially with style and syntax adaptations,
		because POPS has its own core data structure ('Str') for PDB data.
	For a detailed code comparison use a tool like 'vimdiff'.
    Copyright (C) UCL / Dr. Andrew C. R. Martin 1988-2015
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "getpdbml.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/* parse PDB file in XML format */
int readpdbml(FILE *pdbInFile, Str *str, int coarse, int hydrogens)
{
	int ret = 0;
	int bufferSize = 1024;
	char xml_buffer[bufferSize];
	xmlParserCtxtPtr ctxt;
	xmlDoc *document;
	xmlNode *root_node = 0; 
	xmlNode *sites_node = 0;
    xmlNode *atom_node = 0; 
    xmlNode *n = 0;

	/* Generate Document From Filehandle */
	ret = fread(xml_buffer, 1, bufferSize, pdbInFile);
	ctxt = xmlCreatePushParserCtxt(NULL, NULL, xml_buffer, bufferSize, "file");
	while ((ret = fread(xml_buffer, 1, bufferSize, pdbInFile)) > 0) {
		xmlParseChunk(ctxt, xml_buffer, bufferSize, 0);
	}
	xmlParseChunk(ctxt, xml_buffer, 0, 1);
	document = ctxt->myDoc;
	xmlFreeParserCtxt(ctxt);

	if(document == 0) {
		/* Error: Failed to parse file */
		xmlFreeDoc(document); /* free document */
		xmlCleanupParser(); /* clean up xml parser */
		str->nAtom = -1; /* indicate error */
		return(str); /* return structure */
	}

	/* Parse Document Tree */
	root_node = xmlDocGetRootElement(document);
	if(root_node == 0) {
		/* Error: Failed to set root node */
		xmlFreeDoc(document); /* free document */
		xmlCleanupParser(); /* clean up xml parser */
		str->nAtom = -1; /* indicate error */
		return(str);  /* return structure */
	}

	return 0;
}

