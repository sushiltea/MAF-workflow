###################################################################################################
# Functions for downloading open access files form GDC data portal https://portal.gdc.cancer.gov/ #
# Request - Download - Control - Extract                                                          #
#                                                                                                 #
# 09/2019                                                                                         #
###################################################################################################

import requests
import json
import re
import os
import tarfile


def request_data(project, filetype):
    """
    request data from GDC data portal https://portal.gdc.cancer.gov/
    :param project: GDC projectname e.g.: "TCGA-GBM"
    :param filetype: Specify file type e.g.: MAF, TXT, SVS...
    :return: post-request response
    """
    print("Start request")
    fields = [
        "Id"
        # "file_name",
        # "cases.project.project_id"
        # "cases.submitter_id",
        # "cases.samples.sample_type"
        # "cases.disease_type",
 ]
    fields = ",".join(fields)
    files_endpt = "https://api.gdc.cancer.gov/files"

    # This set of filters is nested under an 'and' operator.
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [project]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.access",
                    "value": ["open"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_format",
                    "value": [filetype]
                }
            }
        ]
    }

    # A POST is used, so the filter parameters can be passed directly as a Dict object.
    params = {
        "filters": filters,
        "fields": fields,
        "format": "JSON",
        "size": "2000"
    }

    # The parameters are passed to 'json' rather than 'params' in this case
    response = requests.post(files_endpt, headers={"Content-Type": "application/json"}, json=params)

    #print(response.content.decode("utf-8"))
    print("Request complete")
    return response


def download_data(request, download_path="downloads"):
    """
    download files from post-request response
    :param request: post-request response
    :param download_path: define download path, default = Downloads
    :return: uuid list (file_uuid_list), download folder name (folder_name)
    """
    print("Download started")
    file_uuid_list = []
    # This step populates the download list with the file_ids from the previous query
    for file_entry in json.loads(request.content.decode("utf-8"))["data"]["hits"]:
        file_uuid_list.append(file_entry["file_id"])
    data_endpt = "https://api.gdc.cancer.gov/data"
    params = {"ids": file_uuid_list}
    response = requests.post(data_endpt, data=json.dumps(params), headers={"Content-Type": "application/json"})
    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]
    # write to disk
    data = response.content
    if not os.path.isdir(download_path):
        os.mkdir(download_path)
    path = os.path.join(download_path, file_name)
    print(path)
    with open(path, 'wb') as s:
        s.write(data)

    print("Download complete")
    return file_uuid_list, os.path.abspath(path)


def check_tar_files(file_uuid_list, file_name):
    """
    Quality control step: check if downloaded files exist and are not empty
    :param file_uuid_list: list of uuids
    :param file_name: downloaded file
    :return: print statement containing a list of all requestet files and their status
    """
    if os.path.exists(file_name):
        tar = tarfile.open(file_name, "r:gz")
        tar_dict = dict()
        for member in tar.getmembers():
            mem_split = member.name.split("/")
            size = member.size / 1000000  # get size in MB
            try:
                if not mem_split[1]:  # skip index errors (single files without uuid - filename pair)
                    pass
                uuid = mem_split[0]
                tar_dict[uuid] = {}
                tar_dict[uuid]["filename"] = mem_split[1]
                tar_dict[uuid]["filesize"] = size

            except:
                pass
        for uuid in file_uuid_list:
            # list(tar_dict.keys()) == file_uuid_list
            if uuid not in list(tar_dict.keys()):
                print("UUID", uuid, "is missing!")
                tar_dict[uuid] = {}
                tar_dict[uuid]["status"] = "missing"
                tar_dict[uuid]["filesize"] = "0       "
                tar_dict[uuid]["filename"] = "unkown"

            else:
                tar_dict[uuid]["status"] = "OK"
        print("Download Protocol:\n" + "UUID" + " ".ljust(len(list(tar_dict.keys())[0])) + " Size (MB)\t" + " Status" + "  " + "Filename")
        for key in tar_dict:
            print(key, "\t", tar_dict[key]["filesize"], "\t", tar_dict[key]["status"], "\t", tar_dict[key]["filename"])
    else:
        print("File: " + file_name + " Download folder does not exist! Please check download.")


def extract_files(file_path):
    """
    extract the downloaded folder (except MANIFEST.txt) for further analysis, but not the content within the folder
    :param file_path: download path containing the files
    :return:
    """
    # create directory with same name like the downloaded file if it doesn't exist
    folder_name = os.path.basename(file_path)
    directory = os.path.abspath(os.path.join(os.path.dirname(file_path), folder_name.split(".tar.gz")[0]))
    if not os.path.exists(directory):
        print("Creating directory: ", directory)
        os.makedirs(directory)

    tar = tarfile.open(file_path, "r:gz")
    for item in tar:
        item.name = os.path.basename(item.name)
        if item.name != "MANIFEST.txt":
            tar.extract(item, directory)
    print("Files Extracted to ", directory)
    print("Removing file: ", file_path)
    os.remove(file_path)
    return directory



