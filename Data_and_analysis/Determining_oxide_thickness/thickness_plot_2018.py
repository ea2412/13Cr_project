import matplotlib
import os
import re
import string
import math

matplotlib.use("TkAgg")  # Fixes some weird error
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd

import sys

sys.path.append(os.path.abspath(os.path.join(__file__, os.path.pardir, os.path.pardir)))

from mass_list_py_13Cr import mass_list


new_mass_list = ["O-"]


file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder = "comparable_scans_neg"

output_file_name = "thicknesses_O-_untreated_treated_time.pdf"


class Files:
    untreated_RT = "a23914_1.txt"
    treated_RT_1 = "a23916_1.txt"
    treated_RT_2 = "a23924_1.txt"
    untreated_200 = "a23910_1.txt"
    treated_200_1 = "a23911_1.txt"
    treated_200_2 = "a23925_1.txt"
    untreated_300 = "a23907_1.txt"
    treated_300_1 = "a23906_1.txt"
    treated_300_2 = "a23928_1.txt"
    treated_300_3 = "a23929_1.txt"
    treated_300_4 = "a23930_1.txt"
    treated_300_5 = "a23931_1.txt"


class FilesNames:
    untreated_RT = "RT untreated"
    treated_RT_1 = "RT treated 1"
    treated_RT_2 = "RT treated 2"
    untreated_200 = "200$^o$C untreated"
    treated_200_1 = "200$^o$C treated 1"
    treated_200_2 = "200$^o$C treated 2"
    untreated_300 = "300$^o$C untreated"
    treated_300_1 = "300$^o$C treated 1 (23906)"
    treated_300_2 = "300$^o$C treated 2 (23928)"
    treated_300_3 = "300$^o$C treated 3 (23929)"
    treated_300_4 = "300$^o$C treated 4 (23930)"
    treated_300_5 = "300$^o$C treated 5 (23931)"


class Temp:
    untreated_RT = 25
    treated_RT_1 = 25
    treated_RT_2 = 25
    untreated_200 = 200
    treated_200_1 = 200
    treated_200_2 = 200
    untreated_300 = 300
    treated_300_1 = 300
    treated_300_2 = 300
    treated_300_3 = 300
    treated_300_4 = 300
    treated_300_5 = 300


DP_file_name = [
    Files.untreated_RT,
    Files.treated_RT_1,
    Files.untreated_200,
    Files.treated_200_2,
    Files.untreated_300,
    Files.treated_300_2,
]


sample_names = [
    FilesNames.untreated_RT,
    FilesNames.treated_RT_1,
    FilesNames.untreated_200,
    FilesNames.treated_200_2,
    FilesNames.untreated_300,
    FilesNames.treated_300_2,
]

temperature = [
    Temp.untreated_RT,
    Temp.treated_RT_1,
    Temp.untreated_200,
    Temp.treated_200_2,
    Temp.untreated_300,
    Temp.treated_300_2,
]


class SputterRates:
    RT = 0.11
    T100 = 0.11
    T200 = 0.08
    T300 = 0.08
    T400 = 0.07


rate = SputterRates.T300

sputter_rate = [
    SputterRates.RT,
    SputterRates.RT,
    SputterRates.T200,
    SputterRates.T200,
    SputterRates.T300,
    SputterRates.T300,
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
sputtertime = []

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
        sputtertime_05 = ((x2 - x1) / (y2 - y1)) * (-y1) + x1
        sputtertime.append(sputtertime_05)
    # label = re.sub(r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species))


colours = ["red", "blue", "pink", "black", "purple", "orange"]

for k in range(len(sample_names)):
    ax.scatter(
        temperature[k], sputtertime[k], color=colours[k], s=4, label=sample_names[k],
    )

# for k in range(len(sample_names)):
#     ax.scatter(
#         temperature[k],
#         oxide_thickness[k],
#         color=colours[k],
#         s=4,
#         label=sample_names[k],
#     )

ax.legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.5),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)

ax.set_xlabel("Temperature")
ax.set_ylabel("Sputter time")

fig.savefig(os.path.join(file_dir, output_file_name))
