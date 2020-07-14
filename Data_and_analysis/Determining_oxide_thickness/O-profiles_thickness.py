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
DP_data_folder_2015 = "depth_profiles_temp"
DP_data_folder_2018 = "comparable_scans_neg"

output_file_name = "O-_time.pdf"

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

DP_file_name_2018_t = [
    Files.treated_2018_RT,
    Files.treated_2018_200_2,
    Files.treated_2018_300_2,
]

DP_file_name_2018_all = [
    Files.untreated_2018_RT,
    Files.treated_2018_RT,
    Files.untreated_2018_200,
    Files.treated_2018_200_2,
    Files.untreated_2018_300,
    Files.treated_2018_300_2,
]

DP_file_name_2018_u = [
    Files.untreated_2018_RT,
    Files.untreated_2018_200,
    Files.untreated_2018_300,
]


sample_names_2015 = [
    FilesNames.untreated_RT,
    FilesNames.untreated_100,
    FilesNames.untreated_200,
    FilesNames.untreated_300,
    FilesNames.untreated_400,
]

sample_names_2018_all = [
    FilesNames.untreated_RT,
    FilesNames.treated_RT,
    FilesNames.untreated_200,
    FilesNames.treated_200,
    FilesNames.untreated_300,
    FilesNames.treated_300,
]


sample_names_2018_u = [
    FilesNames.untreated_RT,
    FilesNames.untreated_200,
    FilesNames.untreated_300,
]

sample_names_2018_t = [
    FilesNames.treated_RT,
    FilesNames.treated_200,
    FilesNames.treated_300,
]

colours = [
    "yellow",
    "yellow",
    "yellow",
    "yellow",
    "yellow",
    "blue",
    # "red",
    "blue",
    # "red",
    "blue",
    # "red",
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

temperature_2018_u = [
    Temp.t_RT,
    Temp.t_200,
    Temp.t_300,
]


class SputterRates2015:
    RT = 0.11
    T100 = 0.11
    T200 = 0.08
    T300 = 0.08
    T400 = 0.07


rate_2015 = 0.1

sputter_rates = [
    SputterRates2015.RT,
    SputterRates2015.T100,
    SputterRates2015.T200,
    SputterRates2015.T300,
    SputterRates2015.T400,
]

number_of_plots = 4
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

rate = 0.1

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

for file_name in DP_file_name_2018_u:
    dataframe = pd.read_csv(
        os.path.join(data_dir, DP_data_folder_2018, file_name),
        header=1,
        skiprows=[2, 3],
        sep="\t",
    )
    panda_list_2018.append(dataframe)

panda_both = panda_list_2015 + panda_list_2018
sample_names = sample_names_2015 + sample_names_2018_u
temperatures = temperature_2015 + temperature_2018_u

# Profiles
ax = axs[0]
for k, data in enumerate(panda_list_2015):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        # ax.plot(x_values, norm_y_values, label=sample_names[k])
        ax.scatter(
            x_values, norm_y_values, s=6, label=sample_names_2015[k],
        )
# ax.hlines(0.5, 0, 1000)
ax.set_xlim([0, 1000])
ax.set_xlabel("Sputter time (s)")
ax.set_ylabel("Normalised Counts")
# ax.title.set_text("Profiles of CrO-")
ax.legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)

ax = axs[2]
for k, data in enumerate(panda_list_2018):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        # ax.plot(x_values, norm_y_values, label=sample_names[k])
        ax.scatter(
            x_values * (350 ** 2 / 300 ** 2),
            norm_y_values,
            s=6,
            label=sample_names_2018_u[k],
        )
# ax.hlines(0.5, 0, 1000)
ax.set_xlim([0, 1000])
ax.set_xlabel("Sputter time (s)")
ax.set_ylabel("Normalised Counts")
# ax.title.set_text("Profiles of CrO-")
ax.legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)

# plot thicknesses

ax = axs[1]
oxide_sputtertime = []
oxide_sputtertime_2015 = []
oxide_sputtertime_2018 = []

for l, data in enumerate(panda_list_2015):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]][5:].values
        y_values = data[species][5:].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        if l <= 3:
            finding_05 = norm_y_values - 0.5
            y2 = max(finding_05[finding_05 < 0])
            y1 = min(finding_05[finding_05 > 0])
            index_y2 = np.where(finding_05 == y2)[0][0]
            index_y1 = np.where(finding_05 == y1)[0][0]
            x1 = x_values[index_y1]
            x2 = x_values[index_y2]
            gradient, yintercept = np.polyfit([x1, x2], [y1, y2], 1)
            xintercept = abs(-yintercept / gradient)
            # oxide_sputtertime_2015.append(xintercept)
        if l > 3:
            data["norm"] = data[species] / max_y_value
            middle_section = data[(data.norm > 0.4) & (data.norm < 0.6)]
            data["centred"] = data["norm"] - 0.5
            gradient, yintercept = np.polyfit(data.iloc[:, 0], data.centred, 1)
            xintercept = abs(-yintercept / gradient)
        oxide_sputtertime_2015.append(xintercept)
thickness_2015 = [time * rate_2015 for time in oxide_sputtertime_2015]
for k in range(len(panda_list_2015)):
    ax.scatter(
        temperature_2015[k], oxide_sputtertime_2015[k], label=sample_names_2015[k], s=7,
    )
    # ax.scatter(
    #     temperature_2015[k], thickness_2015[k], label=sample_names_2015[k], s=7,
    # )
# ax.set_xlim([0, 410])
# ax.set_ylim([0, 60])
ax.set_ylim([0, 600])
ax.set_xlabel("Temperature ($^o$C)")
# ax.set_ylabel("Oxide thickness / nm")
ax.set_ylabel("Sputter time / s")
ax.legend(
    loc="upper right",
    bbox_to_anchor=(0.3 / mass_spectra_ar, 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)


ax = axs[3]
oxide_sputtertime = []
oxide_sputtertime_2015 = []
oxide_sputtertime_2018 = []

for l, data in enumerate(panda_list_2018):
    for m, species in enumerate(new_mass_list):
        x_values = data[data.columns[0]][5:].values
        y_values = data[species][5:].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        if l <= len(panda_list_2018):
            finding_05 = norm_y_values - 0.5
            y2 = max(finding_05[finding_05 < 0])
            y1 = min(finding_05[finding_05 > 0])
            index_y2 = np.where(finding_05 == y2)[0][0]
            index_y1 = np.where(finding_05 == y1)[0][0]
            x1 = x_values[index_y1]
            x2 = x_values[index_y2]
            gradient, yintercept = np.polyfit([x1, x2], [y1, y2], 1)
            xintercept = abs(-yintercept / gradient)
        # if l > 2:
        #     data["norm"] = data[species] / max_y_value
        #     middle_section = data[(data.norm > 0.3) & (data.norm < 0.7)]
        #     data["centred"] = data["norm"] - 0.5
        #     gradient, yintercept = np.polyfit(data.iloc[:, 0], data.centred, 1)
        #     xintercept_lower = abs(-yintercept / gradient)
        oxide_sputtertime_2018.append(xintercept)
adjusted_oxide_sputtertime_2018 = [
    st * (350 ** 2 / 300 ** 2) for st in oxide_sputtertime_2018
]
for k in range(len(panda_list_2018)):
    ax.scatter(
        temperature_2018_u[k],
        adjusted_oxide_sputtertime_2018[k],
        label=sample_names_2018_u[k],
        s=7,
    )
ax.set_xlim([0, 410])
ax.set_ylim([0, 200])
ax.set_xlabel("Temperature ($^o$C)")
ax.set_ylabel("Sputter time / s")
# ax.title.set_text("Time to sputter through CrO layer")
ax.legend(
    loc="upper right",
    bbox_to_anchor=(0.3 / mass_spectra_ar, 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)


# thickness_2015 = [
#     time * rate for time, rate in zip(oxide_sputtertime_2015, sputter_rates)
# ]

# thickness_2015 = [time * rate_2015 for time in oxide_sputtertime_2015]


# for l, data in enumerate(panda_list_2018):
#     for m, species in enumerate(new_mass_list):
#         x_values = data[data.columns[0]][3:].values
#         y_values = data[species][3:].values
#         max_y_value = y_values.max()
#         norm_y_values = y_values / max_y_value
#         finding_05 = norm_y_values - 0.5
#         y2 = max(finding_05[finding_05 < 0])
#         y1 = min(finding_05[finding_05 > 0])
#
#         sputtertime = ((x2 - x1) / (y2 - y1)) * (-y1) + x1
#         print(sputtertime)
#         oxide_sputtertime_2018.append(sputtertime)

# assumed_thicknesses = [thickness_2015[0]] + [thickness_2015[2]] + [thickness_2015[3]]
# print(assumed_thicknesses)

# rate_2018 = 0.1 * (350 ** 2 / 300 ** 2)

# thickness_2018 = [time * rate_2018 for time in oxide_sputtertime_2018]


# thicknesses = thickness_2015 + thickness_2018

# temp_K = [temp + 273 for temp in temperatures]
# log_temp = np.log(temp_K)
# log_thickness = np.log(thicknesses)

# cmap = plt.get_cmap("jet")
# colors = cmap(np.linspace(0, 1.0, len(sample_names)))


# for k in range(len(sample_names)):
#     ax.scatter(
#         temperatures[k], thicknesses[k], color=colours[k], s=6, label=sample_names[k],
#     )

# # for k in range(len(sample_names)):
# #     ax.scatter(
# #         temperatures[k], thicknesses[k], color=colours[k], s=6, label=sample_names[k],
# #     )
# ax.set_ylabel("Film thickness (nm)")

# ax.legend(
#     loc="upper right",
#     bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.5),
#     fancybox=False,
#     prop={"size": 4},
#     frameon=True,
#     ncol=1,
#     markerscale=2,
# )

# ax.set_xlabel("Temperature (K)")


fig.savefig(os.path.join(file_dir, output_file_name))
