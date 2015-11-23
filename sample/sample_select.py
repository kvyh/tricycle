import pandas as pd

def save_list(kiclist, filename):
    """
    Save a list of KIC numbers to a file.

    Parameters
    ----------
    kiclist : array_like
        KIC numbers.
    filename: str
        Name of output file.

    """
    text_file = open(filename, "w")

    for kic in kiclist:
        text_file.write("##############\n")
        text_file.write("KIC {0:d}\n".format(kic))
        text_file.write("##############\n\n")

    text_file.close()

df = pd.read_csv("../villanova-db.csv", comment="#")
df["pdepth/sdepth"] = df.pdepth / df.sdepth

df_new = df[(df.pdepth > 0.1) & (df.sdepth > 0.1) & (df.period > 3.0)]
kics = df_new.sort("pdepth/sdepth").KIC.values

list_a = kics[::2]
list_b = kics[1::2]

save_list(list_a, "list_a.txt")
save_list(list_b, "list_b.txt")