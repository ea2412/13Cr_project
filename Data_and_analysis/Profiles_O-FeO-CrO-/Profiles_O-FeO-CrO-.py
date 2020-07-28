import matplotlib
import os
import re
import string
import math
import sys

matplotlib.use("TkAgg")  # Fixes some weird error
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(__file__, os.path.pardir, os.path.pardir)))

from mass_list_py_13Cr import mass_list
from Colors import ColorSet, MarkerSet
from FilesNamesetc import Files, FilesNames, Temp, SputterRates2015

mass_list_O = ["O-"]
mass_list_FeOCrO = ["FeO-", "CrO-"]
mass_list_CrO = ["CrO-"]
mass_list_FeO = ["FeO-"]


file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder_2015 = "depth_profiles_temp"
DP_data_folder_2018 = "comparable_scans_neg"

output_file_name = "OFeOCrO_profiles_labelled.pdf"

number_of_plots = 1
number_of_big_plots = 1


DP_file_name_2015 = [
    Files.untreated_2015_RT,
    Files.untreated_2015_100,
    Files.untreated_2015_200,
    Files.untreated_2015_300,
    Files.untreated_2015_400,
]

sample_names_2015 = [
    "RT",
    "100$^\degree$C",
    "200$^\degree$C",
    "300$^\degree$C",
    "400$^\degree$C",
]


def colours(species):
    if species == "FeO-":
        return ColorSet.FeO
    if species == "CrO-":
        return ColorSet.CrO
    else:
        return "black"


cmap = plt.get_cmap("inferno")
colors_sequential = cmap(np.linspace(0, 1.0, len(sample_names_2015) + 1))
colors_sequential = colors_sequential[:-1]

FeO_label = mpatches.Patch(color=ColorSet.FeO, label="FeO$^-$")
CrO_label = mpatches.Patch(color=ColorSet.CrO, label="CrO$^-$")

temperature_2015 = [
    Temp.t_RT,
    Temp.t_100,
    Temp.t_200,
    Temp.t_300,
    Temp.t_400,
]

rate_average_1sf = 0.1

sputter_rates = [
    SputterRates2015.RT,
    SputterRates2015.T100,
    SputterRates2015.T200,
    SputterRates2015.T300,
    SputterRates2015.T400,
]

fig_width = 7
profile_ar = 7 / 3
profile_ar_1 = 1
padding = [[0.5, 0.3], [0.5, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.5
letter_padding = 0.03

spectra_width = fig_width - sum(padding[0])
spectra_height = spectra_width / profile_ar
spectra_height_1 = spectra_width / profile_ar_1
fig_height = (
    sum(padding[1])
    + number_of_plots * (spectra_height)
    + (number_of_plots + number_of_big_plots - 1) * vertical_gap
    + number_of_big_plots * (spectra_height_1)
)

spectra_width_rel = spectra_width / fig_width
spectra_height_rel = spectra_height / fig_height
spectra_height_1_rel = spectra_height_1 / fig_height
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
                    + i * (spectra_height_rel + vertical_gap_rel)
                ),
                spectra_width_rel,
                spectra_height_rel,
            ]
        )
        for i in range(number_of_plots)
    ]
)


axs_big = np.array(
    [
        fig.add_axes(
            [
                horizontal_padding_rel[0],
                1
                - (
                    vertical_padding_rel[0] / 2
                    + spectra_height_1_rel
                    + number_of_plots * (spectra_height_rel + vertical_gap_rel)
                    + i * (spectra_height_1_rel + vertical_gap_rel)
                ),
                spectra_width_rel,
                spectra_height_1_rel,
            ]
        )
        for i in range(number_of_big_plots)
    ]
)


# Plot data
panda_list_2015 = []
panda_list_2018 = []

for file_name in DP_file_name_2015:
    dataframe = pd.read_csv(
        os.path.join(data_dir, DP_data_folder_2015, file_name),
        header=1,
        skiprows=[2, 3],
        sep="\t",
    )
    panda_list_2015.append(dataframe)

# Profiles
ax = axs[0]
for k, data in enumerate(panda_list_2015):
    for m, species in enumerate(mass_list_O):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        # ax.plot(x_values, norm_y_values, label=sample_names[k])
        # ax.scatter(
        #     x_values,
        #     norm_y_values,
        #     s=6,
        #     label=sample_names_2015[k],
        #     color=colors_sequential[k],
        # )
        ax.scatter(
            x_values * rate_average_1sf,
            norm_y_values,
            s=6,
            label=sample_names_2015[k],
            color=colors_sequential[k],
        )
ax.set_xlim([0, 60])
ax.set_ylim([0, 1.2])
# ax.set_xlim([0, 1000])
# ax.set_xlabel("Sputter time (s)")
ax.set_xlabel("Depth from surface (nm)")
ax.get_yaxis().set_ticks([])
ax.set_ylabel("Normalised Counts")
ax.legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / profile_ar), 1),
    fancybox=False,
    prop={"size": 10},
    frameon=False,
    ncol=1,
    markerscale=2,
)

ax = axs_big[number_of_plots - 1]
for k, data in enumerate(panda_list_2015):
    for m, species in enumerate(mass_list_FeOCrO):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        # ax.scatter(
        #     x_values,
        #     norm_y_values + (k * 1.1),
        #     s=6,
        #     label=sample_names_2015[k],
        #     color=colours(species),
        # )
        ax.scatter(
            x_values * rate_average_1sf,
            norm_y_values + (k * 1.1),
            s=6,
            label=sample_names_2015[k],
            color=colours(species),
        )
        # ax.set_xlim([0, 600])
        ax.set_xlim([0, 60])
        ax.set_ylim([0, 5.7])
        ax.legend(
            handles=[FeO_label, CrO_label],
            loc="upper right",
            bbox_to_anchor=(0.5, 0.15),
            fancybox=False,
            prop={"size": 10},
            frameon=False,
            ncol=1,
            markerscale=2,
        )
ax.get_yaxis().set_ticks([])
# ax.set_xlabel("Sputter time (s)")
ax.set_xlabel("Depth from surface (nm)")
ax.set_ylabel("Normalised Counts")
# ax.text(500, 0.1, "RT", fontsize=10)
# ax.text(500, 1.2, r"100$^\degree$C", fontsize=10)
# ax.text(500, 2.3, r"200$^\degree$C", fontsize=10)
# ax.text(500, 3.4, r"300$^\degree$C", fontsize=10)
# ax.text(500, 4.5, r"400$^\degree$C", fontsize=10)
ax.text(53, 0.1, "RT", fontsize=10)
ax.text(53, 1.2, r"100$^\degree$C", fontsize=10)
ax.text(53, 2.3, r"200$^\degree$C", fontsize=10)
ax.text(53, 3.4, r"300$^\degree$C", fontsize=10)
ax.text(53, 4.45, r"400$^\degree$C", fontsize=10)


for i, ax in enumerate(axs):
    ax.text(
        letter_padding,
        1 - letter_padding * profile_ar,
        string.ascii_lowercase[i],
        transform=ax.transAxes,
        ha="center",
        va="center",
    )

for j, ax in enumerate(axs_big):
    ax.text(
        letter_padding,
        1 - letter_padding * profile_ar_1,
        string.ascii_lowercase[len(axs) + j],
        transform=ax.transAxes,
        ha="center",
        va="center",
    )

fig.savefig(os.path.join(file_dir, output_file_name))
