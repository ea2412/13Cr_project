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
DP_data_folder = "comparable_scans_neg"

output_file_name = "effect_of_crater_size.pdf"

# used RT_1, 200_2, 300_2


class Files:
    treated_300_1 = "a23906_1.txt"
    treated_300_2 = "a23928_1.txt"
    treated_300_3 = "a23929_1.txt"
    treated_300_4 = "a23930_1.txt"
    treated_300_5 = "a23931_1.txt"


class FilesNames:
    treated_300_1 = "300$^o$C treated 1"
    treated_300_2 = "300$^o$C treated 300 1"
    treated_300_3 = "300$^o$C treated 300 2"
    treated_300_4 = "300$^o$C treated 400"
    treated_300_5 = "300$^o$C treated 500"


class Temp:
    t_RT = 25
    t_100 = 100
    t_200 = 200
    t_300 = 300
    t_400 = 400


DP_file_names = [
    Files.treated_300_2,
    Files.treated_300_3,
    Files.treated_300_4,
    Files.treated_300_5,
]


sample_names = [
    FilesNames.treated_300_2,
    FilesNames.treated_300_3,
    FilesNames.treated_300_4,
    FilesNames.treated_300_5,
]

temperature = 300
crater_width = [300, 300, 400, 500]
crater_area = [width ** 2 for width in crater_width]


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

for file_name in DP_file_names:
    dataframe = pd.read_csv(
        os.path.join(data_dir, DP_data_folder, file_name),
        header=1,
        skiprows=[2, 3],
        sep="\t",
    )
    panda_list.append(dataframe)


ax = axs[0]
oxide_sputtertime = []

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


# Adding a zero point
sample_names.append("0")
crater_area.append(0)
oxide_sputtertime.append(0)


cmap = plt.get_cmap("jet")
colors = cmap(np.linspace(0, 1.0, len(sample_names)))

crater_size = crater_area

for k in range(len(sample_names)):
    ax.scatter(
        crater_size[k],
        oxide_sputtertime[k],
        color=colors[k],
        s=7,
        marker="x",
        linewidths=1,
        label=sample_names[k],
    )
ax.set_ylabel("Sputter time (s)")
ax.set_xlabel("Crater area ($\mu$m$^2$)")

plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

array_crater_size = np.array(crater_size)
array_oxide_sputtertime = np.array(oxide_sputtertime)
b, m = polyfit(array_crater_size, array_oxide_sputtertime, 1)
print(b, m)
plt.plot(
    array_crater_size, b + m * array_crater_size, "-", linewidth=1,
)


ax.legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.5),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)


fig.savefig(os.path.join(file_dir, output_file_name))
