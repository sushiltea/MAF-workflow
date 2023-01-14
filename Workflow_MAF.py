###################################################################################################
# Workflow step 1 & 2                                                                             #
# Step1: Request - Download - Control - Extract mafs                                              #
# Step2: load mafs in df - write eset files - merge TCGA pheno together with case_ids             #
# 09/2019                                                                                         #
###################################################################################################

from python_MA_project.step1_step2 import download_GDC_Functions, Maf_to_esetfiles_Functions
import glob

# --------------------------- STEP 1--------------------------- #
# request TCGA-GBM MAFs
request = download_GDC_Functions.request_data("TCGA-GBM", "MAF")
# download the data from request
file_uuid_list, file_path = download_GDC_Functions.download_data(request, download_path="GDC_Downloads")
# check if correct files are downloaded
download_GDC_Functions.check_tar_files(file_uuid_list, file_path)
# extract downloaded folder
folder_path = download_GDC_Functions.extract_files(file_path)


# --------------------------- STEP 2--------------------------- #
# set path to download folder, containing the mafs
files = glob.glob(folder_path + "/*.gz", recursive=True)

# select colnames see https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/ for available columnnames
colnames_selection = ["NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Reference_Allele",
                    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Transcript_ID", "case_id", "COSMIC",
                    "Existing_variation", "t_ref_count", "t_alt_count", "Variant_Type"]


# write all mafs into one dataframe add counts (t_ref_count und t_alt_count),
# add featurename (Chromosome.Start_Position.End_Position.Strand.Reference_Allele.Tumor_Seq_Allele1.Tumor_Seq_Allele2)
maf_df = Maf_to_esetfiles_Functions.prepare_maf_df(files, colnames=colnames_selection)

# create files for eset
df_data, df_featuredata, df_pheno= Maf_to_esetfiles_Functions.write_eset_files(maf_df, folder_path + "/eset/",
                                                                               use_binary_matrix=True)

# now merge TCGA pheno information to case_ids in pheno.csv
pathTCGA = "/pathtofile/TCGA-GBMs_cases.tsv"
phenodata = Maf_to_esetfiles_Functions.merge_pheno_data(folder_path + "/eset/pheno.csv", pathTCGA, folder_path + "/eset/phenodata.csv")

