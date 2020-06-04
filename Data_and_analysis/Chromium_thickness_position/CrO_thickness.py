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


new_mass_list = ["CrO-"]


file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder_2015 = "depth_profiles_temp"
DP_data_folder_2018 = "comparable_scans_neg"

output_file_name = "CrO_2015.pdf"

# used RT_1, 200_2, 300_2


class Files:
    untreated_2015_RT = "a14704_4.txt"
    untreated_2015_100 = "a14705_1.txt"
    untreated_2015_200 = "a14707_1.txt"
    untreated_2015_300 = "a14708_3.txt"
    untreated_2015_400 = "a14710_2.txt"
    untreated_2018_RT = "a23914_1.txt"
    treated_2018_RT = "a23916_1.txt"
    untreated_2018_200 = "a23910_1.txt"
    treated_2018_200_2 = "a23925_1.txt"
    untreated_2018_300 = "a23907_1.txt"
    treated_2018_300_2 = "a23928_1.txt"


class FilesNames:
    untreated_RT = "RT untreated"
    treated_RT = "RT treated"
    untreated_100 = "100$^o$C untreated"
    untreated_200 = "200$^o$C untreated"
    treated_200 = "200$^o$C treated"
    untreated_300 = "300$^o$C untreated"
    treated_300 = "300$^o$C treated"
    untreated_400 = "400$^o$C untreated"


class Temp:
    t_RT = 25
    t_100 = 100
    t_200 = 200
    t_300 = 300
    t_400 = 400


DP_file_name_2015 = [
    Files.untreated_2015_RT,
    Files.untreated_2015_100,
    Files.untreated_2015_200,
    Files.untreated_2015_300,
    Files.untreated_2015_400,
]

DP_file_name_2018 = [
    Files.untreated_2018_RT,
    Files.treated_2018_RT,
    Files.untreated_2018_200,
    Files.treated_2018_200_2,
    Files.untreated_2018_300,
    Files.treated_2018_300_2,
]


sample_names_2015 = [
    FilesNames.untreated_RT,
    FilesNames.untreated_100,
    FilesNames.untreated_200,
    FilesNames.untreated_300,
    FilesNames.untreated_400,
]

sample_names_2018 = [
    FilesNames.untreated_RT,
    FilesNames.treated_RT,
    FilesNames.untreated_200,
    FilesNames.treated_200,
    FilesNames.untreated_300,
    FilesNames.treated_300,
]

colours = [
    "yellow",
    "yellow",
    "yellow",
    "yellow",
    "yellow",
    "blue",
    "red",
    "blue",
    "red",
    "blue",
    "red",
]


temperature_2015 = [
    Temp.t_RT,
    Temp.t_100,
    Temp.t_200,
    Temp.t_300,
    Temp.t_400,
]

temperature_2018 = [
    Temp.t_RT,
    Temp.t_RT,
    Temp.t_200,
    Temp.t_200,
    Temp.t_300,
    Temp.t_300,
]


class SputterRates2015:
    RT = 0.11
    T100 = 0.11
    T200 = 0.08
    T300 = 0.08
    T400 = 0.07


sputter_rates = [
    SputterRates2015.RT,
    SputterRates2015.T100,
    SputterRates2015.T200,
    SputterRates2015.T300,
    SputterRates2015.T400,
]

number_of_plots = 3
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

for file_name in DP_file_name_2018:
    dataframe = pd.read_csv(
        os.path.join(data_dir, DP_data_folder_2018, file_name),
        header=1,
        skiprows=[2, 3],
        sep="\t",
    )
    panda_list_2018.append(dataframe)

panda_both = panda_list_2015 + panda_list_2018
sample_names = sample_names_2015 + sample_names_2018
temperatures = temperature_2015 + temperature_2018

# ## Line profiles of CrO-

for k, data in enumerate(panda_list_2015):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        # axs[0].plot(x_values, norm_y_values, label=sample_names[k])
        axs[0].scatter(x_values, norm_y_values, s=6, label=sample_names_2015[k])
# axs[0].hlines(0.5, 0, 1000)
axs[0].set_xlim([0, 1000])
axs[0].set_xlabel("Sputter time (s)")
axs[0].set_ylabel("Normalised Counts")
# axs[0].title.set_text("Profiles of CrO-")
axs[0].legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)

# ## Line profiles of CrO-

change_crater_size = 350 ** 2 / 300 ** 2

for k, data in enumerate(panda_list_2018):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        # axs[0].plot(x_values, norm_y_values, label=sample_names[k])
        axs[1].scatter(
            x_values * change_crater_size,
            norm_y_values,
            s=6,
            label=sample_names_2018[k],
        )
# axs[0].hlines(0.5, 0, 1000)
axs[1].set_xlim([0, 300])
axs[1].set_xlabel("Sputter time (s)")
axs[1].set_ylabel("Normalised Counts")
# axs[1].title.set_text("Profiles of CrO-")
axs[1].legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)


### CrO- width
CrO_position = []
CrO_lower_value = []
CrO_higher_value = []

for k, data in enumerate(panda_both):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        index_ymax = np.where(y_values == max_y_value)[0][0]
        xvalue_forymax = x_values[index_ymax]
        CrO_position.append(xvalue_forymax)
        data["norm"] = data[species] / max_y_value
        lower_bound = data[:index_ymax]
        higher_bound = data[index_ymax:]
        lower_range = lower_bound[(lower_bound.norm > 0.3) & (lower_bound.norm < 0.7)]
        lower_range["centred"] = lower_range["norm"] - 0.5
        gradient_lower, yintercept_lower = np.polyfit(
            lower_range.iloc[:, 0], lower_range.centred, 1
        )
        xintercept_lower = abs(-yintercept_lower / gradient_lower)
        CrO_lower_value.append(xintercept_lower)
        higher_range = higher_bound[
            (higher_bound.norm > 0.3) & (higher_bound.norm < 0.7)
        ]
        higher_range["centred"] = higher_range["norm"] - 0.5
        gradient_higher, yintercept_higher = np.polyfit(
            higher_range.iloc[:, 0], higher_range.centred, 1
        )
        xintercept_higher = abs(-yintercept_higher / gradient_higher)
        CrO_higher_value.append(xintercept_higher)
CrO_width = [st2 - st1 for st2, st1 in zip(CrO_higher_value, CrO_lower_value)]
lower_difference = [st2 - st1 for st2, st1 in zip(CrO_position, CrO_lower_value)]
upper_difference = [st2 - st1 for st2, st1 in zip(CrO_higher_value, CrO_position)]
print(CrO_position, CrO_width)
for k in range(len(panda_both)):
    axs[2].errorbar(
        temperatures[k],
        CrO_position[k],
        yerr=[[lower_difference[k]], [upper_difference[k]]],
        fmt=" ",
    )
    axs[2].scatter(
        temperatures[k], CrO_position[k], label=sample_names[k], s=7,
    )
axs[2].set_xlim([0, 325])
axs[2].set_ylim([0, 200])
axs[2].set_xlabel("Temperature ($^o$C)")
axs[2].set_ylabel("Sputter time / s")
# axs[1].title.set_text("Time to sputter through CrO layer")
axs[2].legend(
    loc="upper right",
    bbox_to_anchor=(0.3 / mass_spectra_ar, 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)

fig.savefig(os.path.join(file_dir, output_file_name))
