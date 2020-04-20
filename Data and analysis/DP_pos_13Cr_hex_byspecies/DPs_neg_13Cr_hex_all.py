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
from mass_list_pos import mass_list


file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder = "comparable_scans_pos"

output_file_name = "test_all_species_not_norm"


class Files:
    untreated_RT = "a23913_1.txt"
    treated_RT_1 = "a23917_1.txt"
    treated_RT_2 = "a23923_1.txt"
    untreated_200 = "a23909_1.txt"
    treated_200_1 = "a23912_1.txt"
    treated_200_2 = "a23926_1.txt"
    untreated_300 = "a23908_1.txt"
    treated_300_1 = "a23905_2.txt"
    treated_300_2 = "a23927_1.txt"


class FilesNames:
    untreated_RT = "RT untreated"
    treated_RT_1 = "RT treated 1(23917)"
    treated_RT_2 = "RT treated 2 (23923)"
    untreated_200 = "200$^o$C untreated"
    treated_200_1 = "200$^o$C treated 1 (23912)"
    treated_200_2 = "200$^o$C treated 2 (23926)"
    untreated_300 = "300$^o$C untreated"
    treated_300_1 = "300$^o$C treated 1 (23905)"
    treated_300_2 = "300$^o$C treated 2 (23927)"


DP_file_name = [
    Files.untreated_RT,
    Files.treated_RT_1,
    Files.treated_RT_2,
    Files.untreated_200,
    Files.treated_200_1,
    Files.treated_200_2,
    Files.untreated_300,
    Files.treated_300_1,
    Files.treated_300_2,
]


sample_names = [
    FilesNames.untreated_RT,
    FilesNames.treated_RT_1,
    FilesNames.treated_RT_2,
    FilesNames.untreated_200,
    FilesNames.treated_200_1,
    FilesNames.treated_200_2,
    FilesNames.untreated_300,
    FilesNames.treated_300_1,
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

x_range = [0, 600]
y_range = [0, 1000]


fig_width = 7
mass_spectra_ar = 7 / 3
padding = [[0.7, 0.3], [0.5, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.03

panda_list = []

for file_name in DP_file_name:
    dataframe = pd.read_csv(
        os.path.join(data_dir, DP_data_folder, file_name),
        header=1,
        skiprows=[2, 3],
        sep="\t",
    )
    panda_list.append(dataframe)

num_pdfs = math.ceil(len(mass_list) / 40)

for slice_num in range(num_pdfs):
    slice_mass_list = mass_list[40 * slice_num : 40 * (slice_num + 1)]
    slice_length = len(slice_mass_list)
    spectra_width = fig_width - sum(padding[0])
    spectra_height = spectra_width / mass_spectra_ar
    fig_height = slice_length * (spectra_height + sum(padding[1]))
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
            for i in range(slice_length)
        ]
    )

    for j in range(slice_length):
        ax = axs[j]
        species = slice_mass_list[j]
        for l, data in enumerate(panda_list):
            x_values = data[data.columns[0]][3:].values
            y_values = data[species][3:].values
            max_y_value = y_values.max()
            norm_y_values = y_values / max_y_value
            ax.scatter(x_values, y_values, s=2, label=sample_names[l])
            # x.scatter(x_values * rate, norm_y_values, s=2, label=sample_names[l])
            # ax.scatter(
            #     x_values * rate_test[l], norm_y_values, s=2, label=sample_names[l]
            # )
            ax.text(
                0.95,
                0.95,
                re.sub(r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species)),
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
            )
            ax.set_xlim(x_range[0], x_range[1])
            # ax.set_ylim(y_range[0], y_range[1])
            ax.legend(
                loc="upper right",
                bbox_to_anchor=(1 - (0.08 / mass_spectra_ar), 0.92),
                fancybox=False,
                prop={"size": 8},
                frameon=True,
                ncol=1,
                markerscale=2,
            )

    ax.set_xlabel("Sputter time / s")
    ax.set_ylabel("Normalised counts")

    fig.savefig(
        os.path.join(file_dir, output_file_name + "_" + str(slice_num) + ".pdf"), dpi=50
    )
