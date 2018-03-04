/*==============================================================================
json.c : JSON routines, using the cJSON library
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "json.h"
#include "cJSON.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/*
{
  "data_resource": "popscomp_asymmetric", # this could be the name of the "asymmetric" calculated data
  "resource_version": "", # you don't have a database, so there is no version number, if I'm correct
  "software_version": "1.0.0", # this would be the version of your software
  "resource_entry_url": "", # this would probably not apply to you - for databases, this is where the entry data would be found
  "release_date": "21/02/2018",
  "pdb_id": "0x00",
  "chains": [
    {
      "chain_label": "A",
      "additional_chain_annotations": {},
      "residues": [
        {
          "pdb_res_label": "1",
          "aa_type": "ALA",
          "additional_residue_annotations": {},
          "site_data": [
            {
              "raw_score": 99.3, # this is the SASA value
              "confidence_score": 1, # this is a confidence score that you calculate for your value; if cannot be applied, delete this item
              "confidence_classification": "high" # this is a textual confidence classification (low, medium, high, null)
            }
          ]
        },
        {
          "pdb_res_label": "2",
          "aa_type": "HIS",
          "additional_residue_annotations": {},
          "site_data": [
            {
              "site_id_ref": 1, # this can show that this particular residue belongs to a site
              "raw_score": 134.5, # this is the SASA value
              "confidence_score": 1,
              "confidence_classification": "high"
            }
          ]
        },
        {
          "pdb_res_label": "3",
          "aa_type": "GLU",
          "additional_residue_annotations": {},
          "site_data": [
            {
              "raw_score": 85.1, # this is the SASA value
              "confidence_score": 1,
              "confidence_classification": "high"
            }
          ]
        }
      ]
    }
  ],
  "sites": [
    {
      "site_id": 1,
      "label": "interface", # this can collect all residues on a specific interface
      "source_database": "pdb",
      "source_accession": "0x00",
      "source_release_date": "01/01/2017"
    }
  ],
  "additional_entry_annotations": {},
  "evidence_code_ontology": [
    {
      "eco_term": "computational combinatorial evidence used in automatic assertion",
      "eco_code": "ECO:0000246"
    }
  ]
}
*/

/*____________________________________________________________________________*/
void print_json(Arg *arg, cJSON *json)
{
	/* print JSON object to string */
	char *popsOutJson = cJSON_Print(json);

	/* print string to file */
	arg->jsonOutFile = safe_open(arg->jsonOutFileName, "w");
	fprintf(arg->jsonOutFile, "%s", popsOutJson);
	fclose(arg->jsonOutFile);

	free(popsOutJson);
}


/*____________________________________________________________________________*/
/* The natural way to make residues dependent on chains would be a
	double loop iterating over chains and residues. The current data
	is such that it is easier (but less natural) to iterate over residues
	and to create chains on the fly, to which new residue array are attached. */
void make_resSasaJson(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *json)
{
	/* 'json' is the root object to which everything else will be attached */

	/* indices to iterate through arrays defined below */
	unsigned int r = 0; /* residue index */
	char isChainLabel[2] = "@"; /* dummy chain name, never in structure */

	/* header: attached to 'json' */
	cJSON_AddStringToObject(json, "data_resource", "popscomp_asymmetric");
	cJSON_AddStringToObject(json, "resource_version", "");
	cJSON_AddStringToObject(json, "software_version", "3.0.0");
	cJSON_AddStringToObject(json, "resource_entry_url", "");
	cJSON_AddStringToObject(json, "release_date", "21/02/2018");
	cJSON_AddStringToObject(json, "pdb_id", "0x00");

	/* add Chain array */
	cJSON *chains = cJSON_AddArrayToObject(json, "chains");

	/* add Residue array for new Chain */
	cJSON *residues = cJSON_AddArrayToObject(chains, "residues");

	/* iterate over all Residues */
	for (r = 0; r < pdb->nResidue; ++ r) { 
		/*
		if (strcmp(pdb->atom[pdb->resAtom[r]].chainIdentifier, isChainLabel) != 0) {
			strcpy(isChainLabel, pdb->atom[pdb->resAtom[r]].chainIdentifier);
			cJSON *chain = cJSON_CreateObject();
			cJSON_AddItemToArray(chains, chain);
			cJSON_AddStringToObject(chain, "chain_label", isChainLabel); 
		}
		*/

		/* add Residue */
		cJSON *residue = cJSON_CreateObject();
		cJSON_AddItemToArray(residues, residue);
		cJSON_AddNumberToObject(residue, "pdb_res_label", r);
		cJSON_AddStringToObject(residue, "aa_type", \
								pdb->atom[pdb->resAtom[r]].residueName);
		cJSON_AddStringToObject(residue, "additional_chain_annotations", "");

		/* add Site_Data array */
		cJSON *site_data = cJSON_AddArrayToObject(residue, "site_data");

		/* add Sites */
		cJSON *phil = cJSON_CreateObject();
		cJSON_AddItemToArray(site_data, phil);
		cJSON_AddNumberToObject(phil, "raw_score", resSasa[r].philicSasa);
		cJSON_AddNumberToObject(phil, "confidence_score", 90.0);
		cJSON_AddStringToObject(phil, "confidence_classification", "high");

		cJSON *phob = cJSON_CreateObject();
		cJSON_AddItemToArray(site_data, phob);
		cJSON_AddNumberToObject(phob, "raw_score", resSasa[r].phobicSasa);
		cJSON_AddNumberToObject(phob, "confidence_score", 90.0);
		cJSON_AddStringToObject(phob, "confidence_classification", "high");

		cJSON *total = cJSON_CreateObject();
		cJSON_AddItemToArray(site_data, total);
		cJSON_AddNumberToObject(total, "raw_score", resSasa[r].sasa);
		cJSON_AddNumberToObject(total, "confidence_score", 90.0);
		cJSON_AddStringToObject(total, "confidence_classification", "high");
	}

	/* add Sites array */
	cJSON *sites = cJSON_AddArrayToObject(json, "sites");

	/* add Site information */
	cJSON_AddNumberToObject(sites, "site_id", 1);
	cJSON_AddStringToObject(sites, "label", "interface");
	cJSON_AddStringToObject(sites, "source_database", "pdb");
	cJSON_AddStringToObject(sites, "source_accession", "0x00");
	cJSON_AddStringToObject(sites, "source_release_data", "01/01/2017");
}

