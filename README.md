# MAF Workflow

This repository contains scripts for generating ESET (expression set) compatible files from MAFs (Mutation Annotation Format).


* `download_GDC_Functions.py` contains Request - Download - Control - Extraction of MAFs from GDC


* `Maf_to_esetfiles_Functions.py` loads mafs into a pd.DataFrame, writes eset files and merges TCGA phenodata together with case_ids. TCGA phenodata information has to be provided by the user.

* `Workflow_MAF.py` contains an example workflow how to use the scripts

Tested with python 3.9 
