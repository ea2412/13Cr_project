import pandas as pd
import os

file_dir = os.path.dirname(__file__)
data_dir = os.path.join(file_dir, "data")
DP_data_folder = "depth_profiles_hex"
DP_file_name = "a23928_1.txt"


data = pd.read_csv(
    os.path.join(data_dir, DP_data_folder, DP_file_name),
    header=1,
    skiprows=[2, 3],
    sep="\t",
)

mass_list = list(data.columns.values)

with open("test_list.py", "w") as f:
    f.write("mass_list = [")
    f.write("\n")
    for item in mass_list[1:-1]:
        f.write("    ")
        f.write('"')
        f.write(item)
        f.write('"')
        f.write(",")
        f.write("\n")
    f.write("]")

# print(mass_list)
