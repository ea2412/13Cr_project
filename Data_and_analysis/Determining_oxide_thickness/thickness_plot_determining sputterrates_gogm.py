import matplotlib
import os
import re
import string
import math

matplotlib.use("TkAgg")  # Fixes some weird error
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from numpy.polynomial.polynomial import polyfit
import pandas as pd

import sys

sys.path.append(os.path.abspath(os.path.join(__file__, os.path.pardir, os.path.pardir)))

from mass_list_py_13Cr import mass_list


new_mass_list = ["O-"]


file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder = "depth_profiles_temp"

output_file_name = "tests.pdf"

# measured_depths = [100, 200, 300, 400] # no RT
measured_depths = [
    0.279085583763063,
    0.197906929552957,
    0.192951594104085,
    0.142693411048343,
]


class Files:
    untreated_RT = "a14704_4.txt"
    untreated_100 = "a14705_1.txt"
    untreated_200 = "a14707_1.txt"
    untreated_300 = "a14708_3.txt"
    untreated_400 = "a14710_2.txt"


class FilesNames:
    untreated_RT = "RT untreated"
    untreated_100 = "100$^o$C untreated"
    untreated_200 = "200$^o$C untreated"
    untreated_300 = "300$^o$C untreated"
    untreated_400 = "400$^o$C untreated"


class Temp:
    untreated_RT = 25
    untreated_100 = 100
    untreated_200 = 200
    untreated_300 = 300
    untreated_400 = 400


DP_file_name = [
    Files.untreated_100,
    Files.untreated_200,
    Files.untreated_300,
    Files.untreated_400,
]


sample_names = [
    FilesNames.untreated_100,
    FilesNames.untreated_200,
    FilesNames.untreated_300,
    FilesNames.untreated_400,
]

temperature = [
    Temp.untreated_100,
    Temp.untreated_200,
    Temp.untreated_300,
    Temp.untreated_400,
]


class SputterRates:
    T100 = 0.11
    T200 = 0.08
    T300 = 0.08
    T400 = 0.07


rate = SputterRates.T300

sputter_rates = [
    SputterRates.T100,
    SputterRates.T200,
    SputterRates.T300,
    SputterRates.T400,
]

number_of_plots = 1
# number_of_plots = len(DP_file_name)

fig_width = 7
mass_spectra_ar = 7 / 3
padding = [[0.7, 0.3], [0.5, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.03

spectra_width = fig_width - sum(padding[0])
spectra_height = spectra_width / mass_spectra_ar
fig_height = number_of_plots * (spectra_height + sum(padding[1]))

spectra_width_rel = spectra_width / fig_width
spectra_height_rel = spectra_height / fig_height
vertical_gap_rel = vertical_gap / fig_height
horizontal_gap_rel = horizontal_gap / fig_width
horizontal_padding_rel = [x / fig_width for x in padding[0]]
vertical_padding_rel = [x / fig_height for x in padding[1]]


fig = plt.figure(figsize=(fig_width, fig_height))

axs = np.array(
    [
        fig.add_axes(
            [
                horizontal_padding_rel[0],
                1
                - (
                    vertical_padding_rel[0] / 2
                    + spectra_height_rel
                    + i * (spectra_height_rel + vertical_padding_rel[0])
                ),
                spectra_width_rel,
                spectra_height_rel,
            ]
        )
        for i in range(number_of_plots)
    ]
)


# Plot data
panda_list = []

for file_name in DP_file_name:
    dataframe = pd.read_csv(
        os.path.join(data_dir, DP_data_folder, file_name),
        header=1,
        skiprows=[2, 3],
        sep="\t",
    )
    panda_list.append(dataframe)

cmap = plt.get_cmap("jet")
colors = cmap(np.linspace(0, 1.0, len(sample_names)))

ax = axs[0]
oxide_sputtertime = []
total_sputtertime = []

for l, data in enumerate(panda_list):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]][3:].values
        y_values = data[species][3:].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        finding_05 = norm_y_values - 0.5
        y2 = max(finding_05[finding_05 < 0])
        y1 = min(finding_05[finding_05 > 0])
        index_y2 = np.where(finding_05 == y2)[0][0]
        index_y1 = np.where(finding_05 == y1)[0][0]
        x1 = x_values[index_y1]
        x2 = x_values[index_y2]
        sputtertime = ((x2 - x1) / (y2 - y1)) * (-y1) + x1
        oxide_sputtertime.append(sputtertime)
        total_sputtertime.append(x_values[-1])
    # label = re.sub(r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species))

metal_sputtertime = [
    total - oxide for total, oxide in zip(total_sputtertime, oxide_sputtertime)
]

print(total_sputtertime)
print(oxide_sputtertime)
print(metal_sputtertime)
print(measured_depths)

some_metal_gradients = []

for i in range(3):
    gradient_metal = (
        (oxide_sputtertime[i] * measured_depths[i + 1])
        - (oxide_sputtertime[i + 1] * measured_depths[i])
    ) / (
        (metal_sputtertime[i + 1] * oxide_sputtertime[i])
        - (metal_sputtertime[i] * oxide_sputtertime[i + 1])
    )
    some_metal_gradients.append(gradient_metal)

x_list = []
y_list = []

for i in range(4):
    x = metal_sputtertime[i] / oxide_sputtertime[i]
    x_list.append(x)
    y = measured_depths[i] / oxide_sputtertime[i]
    y_list.append(y)

print(x_list, y_list)

array_x_list = np.array(x_list)
array_y_list = np.array(y_list)
b, m = np.polyfit(array_x_list, array_y_list, 1, full=True)
print("intercept: " + str(b))
print("gradient: " + str(m))
plt.plot(
    array_x_list, b + m * array_x_list, "-", linewidth=1,
)


ax.scatter(x_list, y_list)


print("Metal gradients =")
print(some_metal_gradients)


gradient_oxide_100 = (
    measured_depths[0] - (gradient_metal * metal_sputtertime[0])
) / oxide_sputtertime[0]
print(gradient_oxide_100)

gradient_oxide_200 = (
    measured_depths[1] - (gradient_metal * metal_sputtertime[1])
) / oxide_sputtertime[1]
print(gradient_oxide_200)

# for k in range(len(sample_names)):
#     print(oxide_thickness[k] * sputter_rates[k])

# colours = ["red", "blue", "pink", "black", "purple", "orange"]

# # for k in range(len(sample_names)):
# #     ax.scatter(
# #         temperature[k],
# #         oxide_thickness[k],
# #         color=colours[k],
# #         s=4,
# #         label=sample_names[k],
# #     )
# # ax.set_ylabel("Sputter time (s)")

# for k in range(len(sample_names)):
#     ax.scatter(
#         temperature[k],
#         oxide_thickness[k] * sputter_rates[k],
#         color=colours[k],
#         s=6,
#         label=sample_names[k],
#     )
ax.set_ylabel("Total Depth / Oxide sputtter time")

# ax.legend(
#     loc="upper right",
#     bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.5),
#     fancybox=False,
#     prop={"size": 4},
#     frameon=True,
#     ncol=1,
#     markerscale=2,
# )

ax.set_xlabel("Metal sputter time / Oxide sputter time")


fig.savefig(os.path.join(file_dir, output_file_name))
