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
#include "pdb_structure.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

static void init_atom(Atom *atom) {
	strcpy(atom->recordName, "      ");
	atom->entityId = 0;
	atom->atomNumber = 0;
	strcpy(atom->atomName, "    ");
	strcpy(atom->atomName_raw, "  ");
	strcpy(atom->residueName, "   ");
	atom->residueNumber = 0;
	strcpy(atom->icode, " ");
	strcpy(atom->chainIdentifier, " ");
	atom->pos.x = .0;
	atom->pos.y = .0;
	atom->pos.z = .0;
	atom->altpos = ' ';
	atom->occupancy = .0;
	atom->temperatureFactor = .0;
	atom->next = 0;
	strcpy(atom->element, "  ");
	strcpy(atom->segmentIdentifier, "    ");
	atom->secstr = ' ';
}

/*____________________________________________________________________________*/
/* parse PDB file in XML format */
int readpdbml(FILE *pdbInFile, Str *str, int coarse, int hydrogens)
{
	/* variables */
	int ret = 0;
	int bufferSize = 1024;
	char xml_buffer[bufferSize];
	xmlParserCtxtPtr ctxt;
	xmlDoc *document;
	xmlNode *root_node = 0; 
	xmlNode *sites_node = 0;
    xmlNode *atom_node = 0; 
    xmlNode *n = 0;
	xmlChar *content;
	double content_lf;
	Atom *atom = 0;

	str->multiModel = 0;
	str->nAtom = 0;
	atom = &(str->atom[str->nAtom]);

	/* Generate Document From Filehandle */
	ret = fread(xml_buffer, 1, bufferSize, pdbInFile);
	ctxt = xmlCreatePushParserCtxt(NULL, NULL, xml_buffer, bufferSize, "file");
	while ((ret = fread(xml_buffer, 1, bufferSize, pdbInFile)) > 0) {
		xmlParseChunk(ctxt, xml_buffer, bufferSize, 0);
	}
	xmlParseChunk(ctxt, xml_buffer, 0, 1);
	document = ctxt->myDoc;
	xmlFreeParserCtxt(ctxt);

	if (document == 0) {
		/* Error: Failed to parse file */
		xmlFreeDoc(document); /* free document */
		xmlCleanupParser(); /* clean up xml parser */
		str->nAtom = -1; /* indicate error */
		return 1; /* return error code */
	}

	/* Parse Document Tree */
	root_node = xmlDocGetRootElement(document);
	if (root_node == 0) {
		/* Error: Failed to set root node */
		xmlFreeDoc(document); /* free document */
		xmlCleanupParser(); /* clean up xml parser */
		str->nAtom = -1; /* indicate error */
		return 1;  /* return error code */
	}

   /* Parse Atom Coordinate Nodes */
   for (n = root_node->children; n != 0; n = n->next) {
      /* Find Atom Sites Node */
      if(!strcmp("atom_siteCategory", (char *)n->name)) {
         /* Found Atom Sites */
         sites_node = n;
         break;
      }
   }
   
   if (sites_node == 0) {
      /* Error: Failed to find atom sites */
      xmlFreeDoc(document); /* free document */
      xmlCleanupParser(); /* clean up xml parser */
      str->nAtom = -1; /* indicate error */
      return 1; /* return error code */
   }

   /* Scan through atom nodes and populate PDB list. */
	for (atom_node = sites_node->children; atom_node; atom_node = atom_node->next) {
		if (!strcmp("atom_site", (char *)atom_node->name)) {
			/* Set default values */
			init_atom(atom);

			/* Scan atom node children */
			for (n = atom_node->children; n != NULL; n = n->next) {
				if (n->type != XML_ELEMENT_NODE) {
					continue;
				}
				content = xmlNodeGetContent(n);
				if (content == 0) {
					/* Error: Failed to set node content */
					/*FREELIST(wpdb->pdb,PDB);*/ /* free pdb list */
					xmlFreeDoc(document); /* free document */
					xmlCleanupParser(); /* clean up xml parser */
					str->nAtom = -1; /* indicate error */
					return 1; /* return erroe code */
				}
            
				/* Set PDB values */
				if (!strcmp((char *)n->name, "B_iso_or_equiv")) {
				   sscanf((char *)content, "%lf", &content_lf);
				   atom->temperatureFactor = (float)content_lf;
				} else if (!strcmp((char *)n->name, "Cartn_x")) {
				   sscanf((char *)content, "%lf", &content_lf);
				   atom->pos.x = (float)content_lf;
				} else if (!strcmp((char *)n->name, "Cartn_y")) {
				   sscanf((char *)content,"%lf",&content_lf);
				   atom->pos.y = (float)content_lf;
				} else if (!strcmp((char *)n->name, "Cartn_z")) {
				   sscanf((char *)content, "%lf", &content_lf);
				   atom->pos.z = (float)content_lf;
				} else if (!strcmp((char *)n->name, "auth_asym_id")) {
				   strncpy(atom->chainIdentifier, (char *)content, 2);
				} else if (!strcmp((char *)n->name, "auth_atom_id")) {
				   strncpy(atom->atomName, (char *)content, 5);
				} else if (!strcmp((char *)n->name, "auth_comp_id")) {
				   strncpy(atom->residueName, (char *)content, 3);
				} else if (!strcmp((char *)n->name, "auth_seq_id")) {
				   sscanf((char *)content, "%lf", &content_lf);
				   atom->residueNumber = (float)content_lf;
				} else if (!strcmp((char *)n->name, "pdbx_PDB_ins_code")) {
				   strncpy(atom->icode, (char *)content, 1);
				} else if (!strcmp((char *)n->name, "group_PDB")) {
				   strncpy(atom->recordName, (char *)content, 6);
				} else if (!strcmp((char *)n->name, "occupancy")) {
				   sscanf((char *)content, "%lf", &content_lf);
				   atom->occupancy = (float)content_lf;
				} else if (!strcmp((char *)n->name, "label_alt_id")) {
				   /* Use strlen as test for alt position */
				   atom->altpos = strlen((char *)content) ? content[0]:' ';
				} else if (!strcmp((char *)n->name, "pdbx_PDB_model_num")) {
				   sscanf((char *)content, "%lf", &content_lf);
				   str->modelNumber = (int)content_lf;
				} else if (!strcmp((char *)n->name, "type_symbol")) {
				   strncpy(atom->element, (char *)content, 2);
				} else if (!strcmp((char *)n->name, "label_asym_id")) {
				   if (strlen(atom->chainIdentifier) == 0) {
					  strncpy(atom->chainIdentifier, (char *)content, 1);
				   }
				} else if (!strcmp((char *)n->name, "label_atom_id")) {
				   if(strlen(atom->atomName) == 0) {
					  strncpy(atom->atomName, (char *)content, 4);
				   }
				} else if (!strcmp((char *)n->name, "label_comp_id")) {
				   if(strlen(atom->residueName) == 0) {
					  strncpy(atom->residueName, (char *)content, 3);
				   }
				} else if (!strcmp((char *)n->name, "label_entity_id")) {
				   if((atom->entityId == 0) && (strlen((char *)content) > 0)) {
					  sscanf((char *)content, "%lf", &content_lf);
					  atom->entityId = (float)content_lf;
				   }
				} else if (!strcmp((char *)n->name, "label_seq_id")) {
				   if((strlen((char *)content) > 0)) {
					  sscanf((char *)content, "%lf", &content_lf);
					  atom->residueNumber = (float)content_lf;
				   }
				} else if (!strcmp((char *)n->name, "pdbx_formal_charge")) {
				   sscanf((char *)content, "%lf", &content_lf);
				   atom->formalCharge = (int)content_lf;
				   atom->partialCharge = (float)content_lf;
				} else if (!strcmp((char *)n->name, "seg_id")) {
				   if(strlen(atom->segmentIdentifier) == 0) {
					  strncpy(atom->segmentIdentifier, (char *)content, 4);
				   }
				}

				xmlFree(content);           
			 }
			 
			 /* Set raw atom name
				Note: The text pdb format uses columns 13-16 to store the atom
					  name. By convention, columns 13-14 contain the 
					  right-justified element symbol for the atom.
					  The raw atom name is equivalent to colums 13-16 of a
					  pdb-formatted text file .                             
			 */

			if (strlen(atom->atomName) == 1) {
			   /* copy 1-letter name atnam_raw */
			   strcpy((atom->atomName_raw), " ");
			   strncpy((atom->atomName_raw) + 1, atom->atomName, 7);
			}
			if (strlen(atom->atomName) == 4) {
			   /* copy 4-letter name atnam_raw */
			   strncpy(atom->atomName_raw, atom->atomName, 8);
			} else if (strlen(atom->element) == 1) {
			   strcpy((atom->atomName_raw), " ");
			   strncpy((atom->atomName_raw) + 1, atom->atomName, 7);
			} else {
			   strncpy(atom->atomName_raw, atom->atomName, 4);
			}

			/* Set chain to " " if not already set */
			if (strlen(atom->chainIdentifier) == 0) {
			   strcpy(atom->chainIdentifier, " ");
			}

			/* multi-model assignment */
			str->multiModel = str->modelNumber > str->multiModel ? str->modelNumber : str->multiModel;

			/* TODO: model numbers */
			/*
			if (str->modelNumber != ModelNum) {
				FREELIST(str, PDB);
				str = 0;
				   
				if (model_number > ModelNum) {
					break;
				} else {
					continue;
				}
			}
			*/

			/* Filter: All Atoms */
			if (strncmp(atom->recordName, "ATOM  ", 6)) {
				/* Free str and skip atom  */
				/*FREELIST(str,PDB);*/
				str = 0;
				continue; /* filter */
			}

			/* Add partial occ atom from temp storage to output PDB list */
			/* TODO: not sure about this */
			/*
			if ((NPartial != 0) && strcmp(atom->atomName, store_atomName)) {
				if (StoreOccRankAtom(OccRank, multi, NPartial, &str,
									 &end_pdb, &str->nAtom )) {
					LAST(end_pdb);
					NPartial = 0;
				} else {
					xmlFreeDoc(document);
					xmlCleanupParser();
					wpdb->nAtom = -1;
					return(str);
				}
			}
			*/

			/* Set atom number
			   Note: Cannot use atom site id for atom number
				so base atomNumber on number of atoms stored 
			*/
			++ str->nAtom;

			/* Add partial occupancy atom to temp storage */
			/* TODO: not sure about this */
			/*
			if ((str->altpos != ' ') && (NPartial < MAXPARTIAL)) {
			   blCopyPDB(&multi[NPartial], str);

			   gPDBPartialOcc = TRUE;
				   
			   strncpy(store_atomName, str->atnam, 8);
			   ++ NPartial;

			   str = NULL;
			   continue;
			}
			*/

			/* Store Atom */
			/* TODO: adjust this */
			/*
			if (str == 0) {
			   wpdb->pdb = str;
			   end_pdb = str;
			   str = NULL;
			   str->nAtom = 1;
			} else {
			   end_pdb->next = str;
			   end_pdb = str;
			   str = NULL;
			   ++ str->nAtom;
			}
			*/
		}
	}

	return 0;
}

