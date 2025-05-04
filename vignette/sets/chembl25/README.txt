The ChEMBL20 Postgres database was loaded loaded into the zinc database by John Erwin Jan 2015


Preparation of ChEMBL20 for chemoinformatic analysis requires
selecting targets, compounds, and activities between the them. Each
target may be componsed of one or more components that each correspond
to an entry in UniProt. We propagate activities to each component of
the target. 

target criteria
  target type is 'SINGLE PROTEIN' or 'PROTEIN COMPLEX'
  propagate activities to each component of target

compound
  nchar(canonical smiles) < 1024
  
activity
  activity <= 10 uM and > 0 uM
  Ki, Kd, KD, AC, AD, CC, EC, ED, GI, IA, IC, ID, LD, MIC


1_postgres_export.R

  filter 

   run.sh:

   1) 1_mysql_export.sh
     - this exports activies, smiles, and targets from the MySQL database into data/my_sql exports.
     - See the top mysql_export.sh for more details

   2) 2_standardize_affinities.py
     - Converts affinity values to nM and applies an affinity cutoff

   3) 3_prepare_sea_files.R

     - targets [accession, perf_name, uniprot_entry]
        * Use both functional and binding data (but not ADMET)
	* Some targets have the same accession but slightly different perf_name. FOr each accession, take the perf_name given.
	* Check that there is only one uniprot_entry per accession

     - smiles [canonical_smiles, chembl_id]
        * 342 chembl_ids have the same canonical_smiles
	* filter out compounds wiht canonical_smiles >= 1024 characters

     - affinities [accession, chembl_id, affinity]
        * filter for affinity > 0
        * target and compound must be in the targets and smiles tables
	* use relations = <  ~  <=
        * take the lowest affinity for target compound pair reported

     - Make
        * target_10uM.set
	   * target: uniprot_entry (e.g. 1433G_HUMAN)
           * name: pref_name (e.g. 14-3-3 protein gamma)
           * compound: chembl_id (e.g. CHEMBL1306960)
	       binding or functional activity at <= 10uM:
	       EC50, AC50, GI50, LD50
	       Ki, Kd, IC50, IC60, IC70, IC80, IC90, IC95, IC99
	       and different log forms of these types.
        * compounds_10uM.smi
           * compound: chembl_id (e.g. CHEMBL1306960)
           * smiles: canonical_smiles (e.g. Clc1ccc(cc1)c2ccccc2C(=O)NCC3CCNCC3)
               neutralized
               less than 1024 characters

        * target.set
	   * target: <uniprot_entry>_<activity_threshold> (e.g. 1433G_HUMAN_100)
	       * activity thresholds are .1, 1, 10, 100, 1000, 10000
	   * name: perf_name (e.g. 14-3-3 protein gamma)
	   * compound: chembl_id (same as above)
	* compounds.smi
	   (same as above)
