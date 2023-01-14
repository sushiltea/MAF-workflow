###################################################################################################
# Functions from maf to eset compatible files                                                     #
# load mafs in df - write eset files - merge TCGA pheno together with case_ids                    #
#                                                                                                 #
# 09/2019                                                                                         #
###################################################################################################

import pandas as pd
import os

# Functions to create eset compatible files
# © Sarah Meitz 11.2019


def load_maf(filename, colnames):
    """
    load one maf with specified columns into a dataframe

    :param filename: filename of maf
    :param colnames: specifiy list of colnames
    :return: dataframe containing one maf

    """
    df = pd.read_csv(filename, comment="#", sep="\t", low_memory=False,
        skip_blank_lines=True, header=0, compression = "gzip", usecols=colnames)
    return df


def combine_mafs(files, colnames):
    """
    combine multiple mafs into one dataframe, concatinated by headers

    :param files: file names of all mafs in one folder
    :param colnames: list of colnames, default colnames list
    :return: dataframe with all mafs
    """

    mafs = []
    for f in files:
        df = load_maf(f, colnames)
        mafs.append(df)
        print("read:", f, "\n", "maf shape: ", df.shape)

    df = pd.concat(mafs, axis=0, ignore_index=True)
    print("combined maf dataframe shape: ",df.shape)
    return df


def add_counts(df):
    """
    add counts to dataframe as combination of "t_ref_count":"t_alt_count"
    :param df: combined maf dataframe
    :return: dataframe incl counts
    """
    counts = df["t_ref_count"].map(str) + ":" + df["t_alt_count"].map(str)
    df.insert(loc=0, column="counts", value=counts)
    return df


def add_featurename(df, drop=True):
    """
    add featurename at first column to dataframe as unique identifier and set as index
    :param df: combined maf dataframe
    :return: df incl featurename
    """
    featurename = df["Chromosome"].map(str) + "." + df["Start_Position"].map(str) + "." + df["End_Position"].map(str)\
                  + "." + df["Strand"].map(str) + "." + df["Reference_Allele"].map(str) + "."\
                  + df["Tumor_Seq_Allele1"].map(str) + "." + df["Tumor_Seq_Allele2"].map(str)
    df.insert(loc=0, column="featurename", value=featurename)
    df.set_index("featurename", drop=drop, inplace=True)
    return df


def prepare_maf_df(files, colnames):
    """
    read all maf files within a folder add counts and add featurename column
    :param files: list of mafs
    :return: dataframe of all mafs with additional counts and featurename column
    """
    df_mafs = combine_mafs(files, colnames)
    df_counts = add_counts(df_mafs)
    df = add_featurename(df_counts, drop=True)
    return df


def set_data_file(df, path, use_binary_matrix, filename="data.csv"):
    """
    write data.csv file containing unique featurename as index,
    case_id as header filled with counts (summed by case_id)
    :param df: dataframe containing at least counts, featurename and case_id
    :param path: path where file is stored
    :param use_binary_matrix: if set to True, will create a (0,1)-matrix
    :param filename: name of outputfile, default = data.csv
    :return: data dataframe
    """
    print("extracting assay data..")
    filepath = os.path.join(path + filename)
    df_small = df.drop(df.columns.difference(["counts","case_id"]), 1)
    df_small_unq = col2unq(df_small, groupby=["featurename", "case_id"]) # unique featurename/case_id
    df_sumcounts = df_small_unq.drop(columns="case_id").reset_index()
    df_sumcounts["counts_sum"] = df_sumcounts["counts"].apply(lambda x: counts_addition(x)) # add counts up
    df_pivot = df_sumcounts.pivot(index="featurename", columns="case_id", values="counts_sum")
    #df_pivot.set_index("featurename")
    if use_binary_matrix:
        df_pivot = df_pivot.notnull().astype("int")  # create (0,1)-matrix (nan = 0, else 1)
    else:
        df_pivot.fillna(value="0:0", inplace=True)  # fill nan with 0:0
    df_pivot.to_csv(filepath, sep="\t")
    print("written assay data to ", filepath)
    return df_pivot


def set_featuredata_file(df, df_data, path, filename ="featuredata.csv"):
    """
    write featuredata.csv file containing unique featurename as index and all previous selected columns
    :param df: dataframe with unique (df_uniq) entries containing featurename
    :param df_data: final dataframe of assayData
    :param path: folder path to save featuredata
    :param filename: name of outputfile, default = featuredata.csv
    :return: featuredata dataframe
    """
    print("extracting feature data..")
    filepath = os.path.join(path + filename)
    featuredata = df.copy()
    # drop all irrelevant columns
    featuredata.drop(columns=["counts", "case_id", "t_ref_count", "t_alt_count"], inplace=True)
    df_featuredata = col2unq(featuredata, groupby=["featurename"])

    #print("df_featuredata has shape: ", df_featuredata.shape, ", df_pivot has shape: ", df_data.shape)
    #df_featuredata.set_index("featurename", drop=False, inplace=True)

    # ensure df_featuredata has the same dimensions as df_data (=assayData)
    df_featuredata_joined = df_data.reset_index().join(df_featuredata, on="featurename", how="inner")
    df_featuredata_joined = df_featuredata_joined.set_index("featurename")
    df_featuredata_joined = df_featuredata_joined.drop(df_featuredata_joined.columns.difference(df_featuredata.columns), 1)
    print("shape df_featuredata: ", df_featuredata_joined.shape)
    df_featuredata.to_csv(filepath, sep="\t")
    df_featuredata_joined.to_csv(path+"featuredata.csv", sep="\t")

    print("written featuredata to ", filepath, path+"featuredata.csv")
    return df_featuredata_joined


def set_phenodata_file(df_data, path, filename="pheno.csv"):
    """
    write pheno.csv file case_id as column and headers in same order as in data file
    :param df_data: df output from get_data_file function
    :param filename: name of outputfile, default = pheno.csv
    :return: return pheno dataframe
    """
    print("extracting pheno data..")
    filepath = path+filename
    case = df_data.columns
    pheno = pd.DataFrame(index=case)
    pheno.to_csv(filepath, sep="\t")
    print("written pheno data to ", filepath)
    return pheno


def merge_pheno_data(pheno_data, TCGApheno, phenodata_path):
    """
    merge TCGApheno data to case_ids in pheno.csv file
    :param pheno_data: phenodata file containing only indices (case_ids)
    :param TCGApheno: list of all TCGA case_id and corresponding informations
    :param phenodata_path: path to final phenodata file for eset
    :return: phenodata case_id + added information from TCGA data
    """
    pheno = pd.read_csv(pheno_data, sep="\t")
    phenoTCGA = pd.read_csv(TCGApheno, sep="\t")
    # normal merge or join functions create duplicates
    # this workaround overcomes the problem
    pheno["g"] = pheno.groupby("case_id").cumcount()
    phenoTCGA["g"] = phenoTCGA.groupby("case_id").cumcount()
    phenodata = pheno.merge(phenoTCGA, how="inner").drop("g", 1)
    phenodata.set_index("case_id", inplace=True)
    phenodata.to_csv(phenodata_path, sep="\t")
    os.remove(pheno_data)
    return phenodata


def write_eset_files(df, path, use_binary_matrix=True):
    """
    get all files necessary for files2eset function (R function): data, featuredata, phenodata
    :param df: dataframe with unique featurename and all selected columns
    :return: three dataframes: df_data, df_featuredata, df_pheno
    """
    if not os.path.exists(path):
        os.makedirs(path)
        print("new directory for eset files: "+path)
    print("write data.csv")
    df_data = set_data_file(df, path, use_binary_matrix)
    print("write featuredata.csv")
    df_featuredata = set_featuredata_file(df, df_data, path)
    print("write phenodata.csv")
    df_pheno = set_phenodata_file(df_data, path)
    return df_data, df_featuredata, df_pheno


def col2unq(df, groupby):
    """
    makes columns unique in respect to selected groupby column/s ("featurename")
    :param df: dataframe which columns should be unique
    :param groupby: columns which should become unique in df
    :return: unique dataframe
    """
    grouped = df.drop_duplicates().groupby(groupby)

    def combine(entry):
        try:
            return ";".join(set(entry)) # join
        except: # Exception raised if x is not of type string
            return ";".join(str(x) for x in set(entry))
    dict_for_df = dict()
    for col in df.columns:
        dict_for_df[col] = grouped[col].apply(lambda x: combine(x))
    new_df = pd.DataFrame(dict_for_df)
    return new_df


def counts_addition(counts):
    """
    sum counts up when more than one count exists per featurename/case_id
    :param counts: column containing the counts separated by ";"
    :return: summed up counts column
    """
    counts = counts.split(";")
    result_first = 0
    result_second = 0
    for count in counts: # count looks like e.g. 10:15
        count = count.split(":") # now ["10", "15"]
        first = int(count[0])
        second = int(count[1])
        result_first += first
        result_second += second
    count_sum = str(result_first) + ":" + str(result_second)
    return count_sum




