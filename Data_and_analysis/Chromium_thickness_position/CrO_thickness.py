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
        x_values = data[data.columns[0]][3:].values
        y_values = data[species][3:].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        # axs[0].plot(x_values, norm_y_values, label=sample_names[k])
        axs[0].scatter(x_values, norm_y_values, s=6, label=sample_names[k])
axs[0].hlines(0.5, 0, 1000)
axs[0].set_xlim([150, 600])
axs[0].set_xlabel("Sputter time (s)")
axs[0].set_ylabel("Normalised Counts")
axs[0].title.set_text("Profiles of CrO-")
axs[0].legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.9),
    fancybox=False,
    prop={"size": 4},
    frameon=True,
    ncol=1,
    markerscale=2,
)


# ### CrO- position
# CrO_position = []

# for k, data in enumerate(panda_both):
#     for m, species in enumerate(new_mass_list):
#         x_values = data[data.columns[0]][3:].values
#         y_values = data[species][3:].values
#         max_y_value = y_values.max()
#         norm_y_values = y_values / max_y_value
#         index_ymax = np.where(y_values == max_y_value)[0][0]
#         layer_position = x_values[index_ymax]
#         # axs[0].plot(x_values, norm_y_values, label=sample_names[k])
#         # axs[0].scatter(x_values, norm_y_values, s=6, label=sample_names[k])
#         CrO_position.append(layer_position)
# for k in range(len(sample_names)):
#     axs[0].scatter(
#         temperatures[k], CrO_position[k], color=colours[k], s=6, label=sample_names[k],
#     )
# axs[0].set_xlim([0, 400])
# axs[0].set_xlabel("Temperature ($^o$C)")
# axs[0].set_ylabel("Sputter time / s")
# axs[0].title.set_text("Position of CrO- max")
# axs[0].legend(
#     loc="upper right",
#     bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.9),
#     fancybox=False,
#     prop={"size": 4},
#     frameon=True,
#     ncol=1,
#     markerscale=2,
# )


# ### CrO- width
# CrO_position = []
# CrO_lower_value = []
# CrO_higher_value = []

# for k, data in enumerate(panda_list_2015):
#     for m, species in enumerate(new_mass_list):
#         x_values = data[data.columns[0]][3:].values
#         y_values = data[species][3:].values
#         max_y_value = y_values.max()
#         norm_y_values = y_values / max_y_value
#         index_ymax = np.where(y_values == max_y_value)[0][0]
#         norm_higher_y_values = y_values[:index_ymax] / max_y_value
#         norm_lower_y_values = y_values[index_ymax:] / max_y_value
#         finding_05 = norm_lower_y_values - 0.5
#         ly2 = max(finding_05[finding_05 < 0])
#         ly1 = min(finding_05[finding_05 > 0])
#         index_ly2 = np.where(finding_05 == ly2)[0][0]
#         index_ly1 = np.where(finding_05 == ly1)[0][0]
#         lx1 = x_values[index_ly1]
#         lx2 = x_values[index_ly2]
#         sputtertime1 = ((lx2 - lx1) / (ly2 - ly1)) * (-ly1) + lx1
#         print(sputtertime1)
#         CrO_lower_value.append(sputtertime1)
#         finding_05 = norm_higher_y_values - 0.5
#         hy2 = max(finding_05[finding_05 < 0])
#         hy1 = min(finding_05[finding_05 > 0])
#         index_hy2 = np.where(finding_05 == hy2)[0][0]
#         index_hy1 = np.where(finding_05 == hy1)[0][0]
#         hx1 = x_values[index_hy1]
#         hx2 = x_values[index_hy2]
#         sputtertime2 = ((hx2 - hx1) / (hy2 - hy1)) * (-hy1) + hx1
#         print(sputtertime2)
#         CrO_higher_value.append(sputtertime2)
# CrO_width = [st2 - st1 for st2, st1 in zip(CrO_higher_value, CrO_lower_value)]
# print(CrO_width)
# for k in range(len(sample_names)):
#     axs[0].scatter(
#         temperatures[k], CrO_width[k], color=colours[k], s=6, label=sample_names[k],
#     )
# axs[0].set_xlim([0, 410])
# axs[0].set_xlabel("Temperature ($^o$C)")
# axs[0].set_ylabel("Sputter time / s")
# axs[0].title.set_text("Time to sputter through CrO layer")
# axs[0].legend(
#     loc="upper right",
#     bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.5),
#     fancybox=False,
#     prop={"size": 4},
#     frameon=True,
#     ncol=1,
#     markerscale=2,
# )


# ax = axs[0]
# oxide_sputtertime_2015 = []
# oxide_sputtertime_2018 = []

# for l, data in enumerate(panda_list_2015):
#     for m, species in enumerate(new_mass_list):
#         x_values = data[data.columns[0]][3:].values
#         y_values = data[species][3:].values
#         max_y_value = y_values.max()
#         norm_y_values = y_values / max_y_value
#         finding_05 = norm_y_values - 0.5
#         y2 = max(finding_05[finding_05 < 0])
#         y1 = min(finding_05[finding_05 > 0])
#         index_y2 = np.where(finding_05 == y2)[0][0]
#         index_y1 = np.where(finding_05 == y1)[0][0]
#         x1 = x_values[index_y1]
#         x2 = x_values[index_y2]
#         sputtertime = ((x2 - x1) / (y2 - y1)) * (-y1) + x1
#         print(sputtertime)
#         oxide_sputtertime_2015.append(sputtertime)

# rate_2015 = 0.1
# # thickness_2015 = [
# #     time * rate for time, rate in zip(oxide_sputtertime_2015, sputter_rates)
# # ]

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
#         index_y2 = np.where(finding_05 == y2)[0][0]
#         index_y1 = np.where(finding_05 == y1)[0][0]
#         x1 = x_values[index_y1]
#         x2 = x_values[index_y2]
#         sputtertime = ((x2 - x1) / (y2 - y1)) * (-y1) + x1
#         print(sputtertime)
#         oxide_sputtertime_2018.append(sputtertime)

# assumed_thicknesses = [thickness_2015[0]] + [thickness_2015[2]] + [thickness_2015[3]]
# print(assumed_thicknesses)

# rate_2018 = 0.1 * (350 ** 2 / 300 ** 2)

# thickness_2018 = [time * rate_2018 for time in oxide_sputtertime_2018]

# temperatures = temperature_2015 + temperature_2018
# thicknesses = thickness_2015 + thickness_2018

# temp_K = [temp + 273 for temp in temperatures]

# sample_names = sample_names_2015 + sample_names_2018
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
