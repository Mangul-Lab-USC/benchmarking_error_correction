
import pandas
import re

dataset_name="summary"
####################################################
data = pandas.read_csv("/Users/jaquejbrito/Dropbox/PosDoc/ErrorCorrection/"+dataset_name+"master_summary2.txt", index_col=False)
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


def dataset(list):
    dataset_list = []
    for item in list:
        if "TRA" in item:
            dataset_list.append("TRA")
        elif "IGH" in item:
            dataset_list.append("IGH")
        elif "t1" in item:
            dataset_list.append("T1")
        elif "t3" in item:
            dataset_list.append("T3")
        elif "SRR" in item:
            dataset_list.append("RSR")
        else:
            dataset_list.append("IGH")

    return dataset_list

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
data["Dataset"] = dataset(data["EC Filename"])

data["Trim Percent"] = (data["Base - TP TRIM"] + data["Base - FP TRIM"])/data["Total Bases"]


data["Trim Effeciency"] = (data["Base - TP TRIM"]/data["Base - FP TRIM"])



# z = data.groupby(["Tool", "Dataset", "Coverage"])["Base Sensitiviy"].mean()
# y = data.groupby(["Tool", "Dataset", "Coverage"])["Base Precision"].mean()
# x = data.groupby(["Tool", "Dataset", "Coverage"])["Base Gain"].mean()
# y = y.fillna(1)
# n = z.index

# tra = data[(data["Dataset"]=="TRA")]
# rsr = data[(data["Dataset"]=="RSR")]
# t1 = data[(data["Dataset"]=="T1")]
igh = data[(data["Dataset"]=="IGH")]


####################################################
# tra.to_csv("data/tra_from_master.csv")
# rsr.to_csv("data/rsr_from_master.csv")
# t1.to_csv("data/t1_from_master.csv")
igh.to_csv("/Users/jaquejbrito/Dropbox/PosDoc/ErrorCorrection/"+dataset_name+"_summary.csv")
####################################################
