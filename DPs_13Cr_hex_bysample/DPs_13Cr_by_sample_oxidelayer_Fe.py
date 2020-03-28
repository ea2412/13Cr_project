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
from test_list import mass_list

mass_list = [
    "Fe_3O-",
    "Fe_2O-",
    "Fe_3O_2-",
    "FeO-",
    "Fe_2O_2-",
    "Fe_3O_3-",
    "FeO_2-",
    "FeO_3-",
    "CrO-",
    "O-",
]


file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder = "comparable_scans"

output_file_name = "all_Fe_oxide.pdf"


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


DP_file_name = [
    Files.untreated_RT,
    Files.treated_RT_1,
    Files.untreated_200,
    Files.treated_200_1,
    Files.untreated_300,
    Files.treated_300_2,
]


sample_names = [
    FilesNames.untreated_RT,
    FilesNames.treated_RT_1,
    FilesNames.untreated_200,
    FilesNames.treated_200_1,
    FilesNames.untreated_300,
    FilesNames.treated_300_2,
]


class SputterRates:
    RT = 0.11
    T100 = 0.11
    T200 = 0.08
    T300 = 0.08
    T400 = 0.07


rate = SputterRates.T300

rate_test = [
    SputterRates.RT,
    SputterRates.RT,
    SputterRates.T200,
    SputterRates.T200,
    SputterRates.T300,
    SputterRates.T300,
    SputterRates.T300,
]


x_range = [0, 200]
y_range = [0, 1.1]

number_of_plots = len(DP_file_name)

fig_width = 3.5
mass_spectra_ar = 6 / 2
padding = [[0.7, 0.3], [0.2, 0.7]]  # padding = [[left, right], [top, bottom]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.03

spectra_width = fig_width - sum(padding[0])
spectra_height = spectra_width / mass_spectra_ar
fig_height = (
    number_of_plots * spectra_height
    + (number_of_plots - 1) * vertical_gap
    + (sum(padding[1]))
)

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
                    vertical_padding_rel[0]
                    + spectra_height_rel
                    + i * (spectra_height_rel + vertical_gap_rel)
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


for l, data in enumerate(panda_list):
    ax = axs[l]
    for species in mass_list:
        x_values = data[data.columns[0]][3:].values
        y_values = data[species][3:].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        label = re.sub(r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species))
        ax.scatter(x_values, norm_y_values, s=2, label=label)
        ax.set_xlim(x_range[0], x_range[1])
        ax.set_ylim(y_range[0], y_range[1])
        if l < (len(panda_list) - 1):
            ax.set_xticklabels([])
        ax.legend(
            loc="upper right",
            bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.9),
            fancybox=False,
            prop={"size": 5},
            frameon=True,
            ncol=2,
            markerscale=2,
        )
    ax.text(
        0.8,
        0.9,
        sample_names[l],
        size=5,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )


axs[number_of_plots - 1].set_xlabel("Sputter time / s")
axs[math.ceil(number_of_plots / 2)].set_ylabel("Normalised counts")

fig.savefig(os.path.join(file_dir, output_file_name))
