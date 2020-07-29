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

from matplotlib.ticker import FormatStrFormatter

# from palettable.cartocolors.sequential import YlOrRd_6
# from palettable.cartocolors.sequential import Blues_6

import sys

sys.path.append(os.path.abspath(os.path.join(__file__, os.path.pardir, os.path.pardir)))
from test_list import mass_list
from Colors import Settings
from FilesNamesetc import Files, FilesNames, Temp, SputterRates2015

mass_list_Fe = [
    "Fe_3O-",
    "Fe_2O-",
    "Fe_3O_2-",
    "FeO-",
    "Fe_2O_2-",
    # "Fe_3O_3-",
    "FeO_2-",
    "FeO_3-",
]


mass_list_Cr = [
    # "Cr-",
    "CrO_3-",
    "Cr_2O_5-",
    "Cr_2O_4-",
    "Cr_4O_6-",
    "CrO-",
    "Cr_2O-",
    # "CCr-",
    #     "COCr-",
    #     # "CsCr-",
    #     "Na_3Cr_2O_4-",
    #     "Cr_2FeH-",
    #     "Cr_3O_3-",
    #     "KNa_2Cr_2O_4-",
    #     "CrH_3-",
    #     "C_2NO_2Cr-",
]


file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder = "comparable_scans_neg"

output_file_name = "profiles_normalised_plot_test.pdf"

DP_file_name = [
    Files.untreated_2018_RT,
    Files.treated_2018_RT,
    Files.untreated_2018_200,
    Files.treated_2018_200_2,
    Files.untreated_2018_300,
    Files.treated_2018_300_2,
]

sample_names = [
    FilesNames.untreated_RT_T,
    FilesNames.treated_RT,
    FilesNames.untreated_200,
    FilesNames.treated_200,
    FilesNames.untreated_300,
    FilesNames.treated_300,
]

x_range = [0, 150]
y_range = [0, 1.05]

number_of_plots = len(DP_file_name)
number_of_columns = 2
column_labels = ["Fe-containing species", "Cr-containing species"]

fig_width = 3.5
mass_spectra_ar = 7 / 3
padding = [[0.4, 0.3], [0.2, 0.4]]  # padding = [[left, right], [top, bottom]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.03

spectra_width = (fig_width - sum(padding[0]) - horizontal_gap) / 2
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
        [
            fig.add_axes(
                [
                    horizontal_padding_rel[0]
                    + j * (spectra_width_rel + horizontal_gap_rel),
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
        for j in range(number_of_columns)
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

mass_lists = [mass_list_Fe, mass_list_Cr]
cmap = [
    plt.get_cmap("winter_r"),
    plt.get_cmap("autumn_r"),
]

# cmap = [
#     plt.get_cmap("Blues_r"),
#     plt.get_cmap("YlOrRd_r"),
# ]

# cmap = [
#     plt.get_cmap("winter_r"),
#     plt.get_cmap("autumn_r"),
# ]

# other cmap options: jet,Blues_r, YlOrRd_r

for m in range(number_of_columns):
    for l, data in enumerate(panda_list):
        ax = axs[m, l]
        mass_list = mass_lists[m]
        colors = cmap[m](np.linspace(0, 1.0, len(mass_list)))
        for n, species in enumerate(mass_list):
            x_values = data[data.columns[0]][3:].values
            y_values = data[species][3:].values
            y_values_O = data["O-"][3:].values
            max_y_value = y_values.max()
            max_y_value_O = y_values_O.max()
            norm_y_values = y_values / max_y_value
            norm_y_values_O = y_values_O / max_y_value_O
            label = re.sub(
                r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species)
            )
            label_O = re.sub(
                r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", "O-")
            )
            # ax.scatter(x_values, norm_y_values, s=1, label=label, color=colors[n])
            ax.plot(
                x_values, norm_y_values, label=label, color=colors[n], linewidth=1,
            )
            # ax.plot(x_values, y_values, label=label, color=colors[n], linewidth=1)
            ax.set_xlim(x_range[0], x_range[1])
            ax.set_ylim(y_range[0], y_range[1])
            ax.get_yaxis().set_ticks([])
            if l < (len(panda_list) - 1):
                ax.set_xticklabels([])
                ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
            else:
                ax.set_xticklabels(
                    np.arange(x_range[0], x_range[1], 50),
                    fontsize=Settings.axis_fontsize,
                )
        ax.plot(
            x_values, norm_y_values_O, label=label_O, color="black", linewidth=1,
        )


axs[0, 0].legend(
    loc="upper right",
    bbox_to_anchor=(1 - 0.4 / mass_spectra_ar, 1),
    fancybox=False,
    prop={"size": Settings.legend_fontsize},
    frameon=False,
    ncol=2,
    markerscale=2,
)

axs[1, 0].legend(
    loc="upper right",
    bbox_to_anchor=(1 - 0.4 / mass_spectra_ar, 1),
    fancybox=False,
    prop={"size": Settings.legend_fontsize},
    frameon=False,
    ncol=2,
    markerscale=2,
)

for m in range(2):
    axs[m, (number_of_plots - 1)].set_xlabel(
        "Sputter time / s", fontsize=Settings.axis_fontsize
    )

fig.text(
    0.5 * horizontal_padding_rel[0],
    vertical_padding_rel[1] + 3 * spectra_height_rel + 2 * vertical_gap_rel,
    "Normalised counts",
    ha="left",
    va="center",
    rotation=90,
    fontproperties=fm.FontProperties(size=Settings.axis_fontsize),
)

for i, text in enumerate(sample_names):
    fig.text(
        1 - (0.5 * horizontal_padding_rel[1]),
        1
        - (
            vertical_padding_rel[0]
            + 0.5 * spectra_height_rel
            + i * (vertical_gap_rel + spectra_height_rel)
        ),
        text,
        ha="center",
        va="center",
        rotation=270,
        fontproperties=fm.FontProperties(size=Settings.axis_fontsize),
    )

for i, text in enumerate(column_labels):
    fig.text(
        horizontal_padding_rel[0]
        + (i + 0.5) * spectra_width_rel
        + i * horizontal_gap_rel,
        1 - (0.5 * vertical_padding_rel[0]),
        text,
        ha="center",
        va="center",
        rotation=0,
        fontproperties=fm.FontProperties(size=Settings.axis_fontsize),
    )

for i in range(len(panda_list)):
    for j in range(number_of_columns):
        ax = axs[j, i]
        ax.text(
            1 - letter_padding / spectra_width_rel,
            1 - letter_padding / spectra_height_rel,
            string.ascii_lowercase[(i * 2) + j],
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="black",
            fontsize=Settings.axis_fontsize,
        )


fig.savefig(os.path.join(file_dir, output_file_name))
