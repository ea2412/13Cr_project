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

file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder = "depth_profiles_all"
DP_file_name = [
    "a14704_1.txt",
    "a14707_1.txt",
    "a14708_1.txt",
    "a23924_1.TXT",
    "a23925_1.TXT",
    "a23928_1.TXT",
    "a23929_1.TXT",
]
sample_names = [
    "RT",
    "200$^o$C",
    "300$^o$C",
    "RT_diesel",
    "200$^o$C_diesel",
    "300$^o$C_diesel",
    "300$^o$C_diesel",
]


sputter_rates = [0.11, 0.11, 0.08, 0.08, 0.07]
x_range = [0, 25]
y_range = [0, 1.1]

output_file_name = "all_comparisons_measured_different_conditions.pdf"

len_species_list = len(mass_list)

fig_width = 7
mass_spectra_ar = 7 / 4
padding = [[0.7, 0.3], [0.5, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.03

spectra_width = fig_width - sum(padding[0])
spectra_height = spectra_width / mass_spectra_ar
fig_height = len_species_list * (spectra_height + sum(padding[1]))

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
                    + i * (spectra_height_rel + vertical_padding_rel[0])
                ),
                spectra_width_rel,
                spectra_height_rel,
            ]
        )
        for i in range(len_species_list)
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


rate = sputter_rates[3]


for j in range(len_species_list):
    ax = axs[j]
    species = mass_list[j]
    for l, data in enumerate(panda_list):
        x_values = data[data.columns[0]][3:].values
        y_values = data[species][3:].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        ax.scatter(x_values * rate, norm_y_values, s=2, label=sample_names[l])
        ax.text(
            0.95,
            0.95,
            re.sub(r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species)),
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        ax.set_xlim(0, 75)
        ax.legend(
            loc="upper right",
            bbox_to_anchor=(1 - (0.08 / mass_spectra_ar), 0.92),
            fancybox=False,
            prop={"size": 8},
            frameon=True,
            ncol=1,
            markerscale=2,
        )

    # max_y_value = max(y_values)
    # index_max = np.where(y_values == max_y_value)[0][0]
    # x_position = x_values[index_max]
    # x_position_distance = x_position * rate
    # ax.axvline(12.16, color="black", linewidth=0.75, dashes=(5, 5))
# if species < len(species_to_plot) - 1:
#     labels = [item.get_text() for item in ax.get_xticklabels()]
#     empty_string_labels = [""] * len(labels)
#     ax.set_xticklabels(empty_string_labels)

# labels = [item.get_text() for item in ax.get_yticklabels()]
# empty_string_labels = [""] * len(labels)
# ax.set_xlim(x_range)
# ax.set_ylim(y_range)
# if len(species_to_plot[species]) > 3:
#     num_col = 2
# else:
#     num_col = 1


ax.set_xlabel("Distance / nm")
ax.set_ylabel("Intensity / 10$^3$ counts")

fig.savefig(os.path.join(file_dir, output_file_name))
