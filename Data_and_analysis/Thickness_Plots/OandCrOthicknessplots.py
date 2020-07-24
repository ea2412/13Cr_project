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

output_file_name = "test1.pdf"

number_of_plots = 2


DP_file_name_2015 = [
    Files.untreated_2015_RT,
    Files.untreated_2015_100,
    Files.untreated_2015_200,
    Files.untreated_2015_300,
    Files.untreated_2015_400,
]

DP_file_name_2018_u = [
    Files.untreated_2018_RT,
    Files.untreated_2018_200,
    Files.untreated_2018_300,
]

sample_names_2015 = [
    FilesNames.untreated_RT_T,
    FilesNames.untreated_100_T,
    FilesNames.untreated_200_T,
    FilesNames.untreated_300_T,
    FilesNames.untreated_400_T,
]

sample_names_2018_u = [
    FilesNames.untreated_RT_T,
    FilesNames.untreated_200_T,
    FilesNames.untreated_300_T,
]

cmap = plt.get_cmap("inferno")
colors_sequential = cmap(np.linspace(0, 1.0, len(sample_names_2015) + 1))
colors_sequential = colors_sequential[:-1]

temperature_2015 = [
    Temp.t_RT,
    Temp.t_100,
    Temp.t_200,
    Temp.t_300,
    Temp.t_400,
]

temperature_2018_u = [
    Temp.t_RT,
    Temp.t_200,
    Temp.t_300,
]

rate_2015 = 0.1

sputter_rates = [
    SputterRates2015.RT,
    SputterRates2015.T100,
    SputterRates2015.T200,
    SputterRates2015.T300,
    SputterRates2015.T400,
]

fig_width = 7
mass_spectra_ar = 4 / 3
padding = [[0.7, 0.3], [0.5, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.5
letter_padding = 0.03

spectra_width = fig_width - sum(padding[0])
spectra_height = spectra_width / mass_spectra_ar
fig_height = (
    sum(padding[1])
    + number_of_plots * (spectra_height)
    + (number_of_plots - 1) * vertical_gap
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

# Plot data
panda_list_2015 = []
panda_list_2018 = []
panda_list_2018_u = []

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
    panda_list_2018_u.append(dataframe)

panda_both_u = panda_list_2015 + panda_list_2018_u
sample_names = sample_names_2015 + sample_names_2018_u
temperatures_u = temperature_2015 + temperature_2018_u

# Plots
ax = axs[0]

oxide_sputtertime_2015 = []
oxide_sputtertime_2018 = []

for l, data in enumerate(panda_list_2015):
    for m, species in enumerate(mass_list_O):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        norm_y_values_first_points_removed = norm_y_values[3:]
        if l <= 3:
            finding_05 = norm_y_values_first_points_removed - 0.5
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

for l, data in enumerate(panda_list_2018_u):
    for m, species in enumerate(mass_list_O):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        norm_y_values_first_points_removed = norm_y_values[3:]
        if l <= len(panda_list_2018_u):
            finding_05 = norm_y_values_first_points_removed - 0.5
            y2 = max(finding_05[finding_05 < 0])
            y1 = min(finding_05[finding_05 > 0])
            index_y2 = np.where(finding_05 == y2)[0][0]
            index_y1 = np.where(finding_05 == y1)[0][0]
            x1 = x_values[index_y1]
            x2 = x_values[index_y2]
            gradient, yintercept = np.polyfit([x1, x2], [y1, y2], 1)
            xintercept = abs(-yintercept / gradient)
        oxide_sputtertime_2018.append(xintercept)

ax.scatter(
    temperature_2015,
    oxide_sputtertime_2015,
    label="Pressure 1",
    s=15,
    marker="x",
    color=colors_sequential,
    # color="black",
)
ax.scatter(
    temperature_2018_u,
    oxide_sputtertime_2018,
    label="Pressure 2",
    s=15,
    marker="s",
    color=[colors_sequential[0], colors_sequential[2], colors_sequential[3]],
    # color="black",
)
ax.set_xlim([0, 425])
# ax.set_ylim([0, 200])
ax.set_xlabel("Temperature ($^o$C)")
ax.set_ylabel("Sputter time / s")
# ax.title.set_text("Time to sputter through CrO layer")
ax.legend(
    loc="upper right",
    bbox_to_anchor=(0.25, 0.95),
    fancybox=False,
    prop={"size": 10},
    frameon=False,
    ncol=1,
    markerscale=2,
)


ax = axs[1]
### CrO- width
CrO_position_2015 = []
CrO_position_2018 = []
CrO_lower_value_2015 = []
CrO_higher_value_2015 = []
CrO_lower_value_2018 = []
CrO_higher_value_2018 = []

for k, data in enumerate(panda_list_2015):
    for m, species in enumerate(mass_list_CrO):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        index_ymax = np.where(y_values == max_y_value)[0][0]
        xvalue_forymax = x_values[index_ymax]
        CrO_position_2015.append(xvalue_forymax)
        data["norm"] = data[species] / max_y_value
        lower_bound = data[:index_ymax]
        higher_bound = data[index_ymax:]
        lower_range = lower_bound[(lower_bound.norm > 0.3) & (lower_bound.norm < 0.7)]
        lower_range["centred"] = lower_range["norm"] - 0.5
        gradient_lower, yintercept_lower = np.polyfit(
            lower_range.iloc[:, 0], lower_range.centred, 1
        )
        xintercept_lower = abs(-yintercept_lower / gradient_lower)
        CrO_lower_value_2015.append(xintercept_lower)
        higher_range = higher_bound[
            (higher_bound.norm > 0.3) & (higher_bound.norm < 0.7)
        ]
        higher_range["centred"] = higher_range["norm"] - 0.5
        gradient_higher, yintercept_higher = np.polyfit(
            higher_range.iloc[:, 0], higher_range.centred, 1
        )
        xintercept_higher = abs(-yintercept_higher / gradient_higher)
        CrO_higher_value_2015.append(xintercept_higher)

for k, data in enumerate(panda_list_2018_u):
    for m, species in enumerate(mass_list_CrO):
        x_values = data[data.columns[0]].values
        y_values = data[species].values
        max_y_value = y_values.max()
        norm_y_values = y_values / max_y_value
        index_ymax = np.where(y_values == max_y_value)[0][0]
        xvalue_forymax = x_values[index_ymax]
        CrO_position_2018.append(xvalue_forymax)
        data["norm"] = data[species] / max_y_value
        lower_bound = data[:index_ymax]
        higher_bound = data[index_ymax:]
        lower_range = lower_bound[(lower_bound.norm > 0.3) & (lower_bound.norm < 0.7)]
        lower_range["centred"] = lower_range["norm"] - 0.5
        gradient_lower, yintercept_lower = np.polyfit(
            lower_range.iloc[:, 0], lower_range.centred, 1
        )
        xintercept_lower = abs(-yintercept_lower / gradient_lower)
        CrO_lower_value_2018.append(xintercept_lower)
        higher_range = higher_bound[
            (higher_bound.norm > 0.3) & (higher_bound.norm < 0.7)
        ]
        higher_range["centred"] = higher_range["norm"] - 0.5
        gradient_higher, yintercept_higher = np.polyfit(
            higher_range.iloc[:, 0], higher_range.centred, 1
        )
        xintercept_higher = abs(-yintercept_higher / gradient_higher)
        CrO_higher_value_2018.append(xintercept_higher)
CrO_width_2015 = [
    st2 - st1 for st2, st1 in zip(CrO_higher_value_2015, CrO_lower_value_2015)
]
lower_difference_2015 = [
    st2 - st1 for st2, st1 in zip(CrO_position_2015, CrO_lower_value_2015)
]
upper_difference_2015 = [
    st2 - st1 for st2, st1 in zip(CrO_higher_value_2015, CrO_position_2015)
]
CrO_width_2018 = [
    st2 - st1 for st2, st1 in zip(CrO_higher_value_2018, CrO_lower_value_2018)
]
lower_difference_2018 = [
    st2 - st1 for st2, st1 in zip(CrO_position_2018, CrO_lower_value_2018)
]
upper_difference_2018 = [
    st2 - st1 for st2, st1 in zip(CrO_higher_value_2018, CrO_position_2018)
]

for k in range(len(temperature_2015)):
    ax.errorbar(
        temperature_2015[k],
        CrO_position_2015[k],
        yerr=[[lower_difference_2015[k]], [upper_difference_2015[k]]],
        fmt=" ",
        color=colors_sequential[k],
        # color="black",
    )
reduced_colours = [colors_sequential[0], colors_sequential[2], colors_sequential[3]]
for m in range(len(temperature_2018_u)):
    ax.errorbar(
        temperature_2018_u[m],
        CrO_position_2018[m],
        yerr=[[lower_difference_2018[m]], [upper_difference_2018[m]]],
        fmt=" ",
        color=reduced_colours[m],
        # color="black",
    )
ax.scatter(
    temperature_2015,
    CrO_position_2015,
    label="Pressure 1",
    s=15,
    marker="x",
    color=colors_sequential,
    # color="black",
)
ax.scatter(
    temperature_2018_u,
    CrO_position_2018,
    label="Pressure 2",
    s=15,
    marker="s",
    color=[colors_sequential[0], colors_sequential[2], colors_sequential[3]],
    # color="black",
)
ax.set_xlim([0, 425])
ax.set_ylim([0, 600])
ax.set_xlabel("Temperature ($^o$C)")
ax.set_ylabel("Sputter time / s")
# ax.title.set_text("Time to sputter through CrO layer")
ax.legend(
    loc="upper right",
    bbox_to_anchor=(0.25, 0.95),
    fancybox=False,
    prop={"size": 10},
    frameon=False,
    ncol=1,
    markerscale=2,
)


fig.savefig(os.path.join(file_dir, output_file_name))
