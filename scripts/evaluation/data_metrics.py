
import pandas
import re

dataset_name="summary"
####################################################
data = pandas.read_csv("./master_summary.txt", index_col=False)
####################################################

def tool_name(list):
    tool_list = []
    for item in list:
        print (item)
        item = item.replace("run.", '')
        item = item.replace(".sh", '')
        tool_list.append(item.capitalize())
    return tool_list


def length(list):
    cover = []
    for item in list:
        if "50" in item:
            cover.append("50L")
        elif "75" in item:
            cover.append("75L")
        elif "100" in item:
            cover.append("100L")
        else:
            cover.append("100")
    return cover


def coverage(list):
    coverage = []
    for item in list:
        try:
            m = re.search("cov_(\d+)", item)
            coverage.append(int(str(m.groups()[0])))
        except:
            coverage.append("NA")

    return coverage



print (data.head())
data["Tool"] = tool_name(data["Wrapper Name"])
data["Coverage"] = coverage(data["EC Filename"])
data["Length"] = length(data["EC Filename"])
data["Base Sensitivity"] = data["Base - TP"] / (data["Base - TP"] + data["Base - FN"])
data["Base Precision"] = data["Base - TP"] / (data["Base - TP"] + data["Base - FP"] + data["Base - FP INDEL"])
data["Base Gain"] = (data["Base - TP"] - (data["Base - FP"] + data["Base - FP INDEL"])) / \
                    (data["Base - TP"] + data["Base - FN"])

data["Base Accuracy"] = (data["Base - TP"] + data["Base - TN"])/(data["Base - TP"] + data["Base - TN"] + data["Base - FP"] + data["Base - FN"] +
                                                                data["Base - FP INDEL"] + data["Base - FP TRIM"] + data["Base - FN WRONG"])

data["Trim Percent"] = (data["Base - TP TRIM"] + data["Base - FP TRIM"])/data["Total Bases"]

data["Trim Efficiency"] = (data["Base - TP TRIM"]/data["Base - FP TRIM"])


####################################################

data.to_csv("./summary.csv")
####################################################
