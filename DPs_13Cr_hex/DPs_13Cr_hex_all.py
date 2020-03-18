import matplotlib
import os
import re
import string

matplotlib.use("TkAgg")  # Fixes some weird error
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd

file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
DP_data_folder = "depth_profiles_hex"
DP_file_name = "a23928_1.txt"

species_to_plot = ["Fe-", "Si_2-", "Cr-", "C_3O-", "CrO-"]
sputter_rates = [0.11, 0.11, 0.08, 0.08, 0.07]
x_range = [0, 25]
y_range = [0, 1.1]

output_file_name = "test_2.pdf"

fig_width = 7
mass_spectra_ar = 7 / 4
padding = [[0.7, 0.3], [0.5, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.03

spectra_width = fig_width - sum(padding[0])
spectra_height = spectra_width / mass_spectra_ar
fig_height = spectra_height + sum(padding[1])

spectra_width_rel = spectra_width / fig_width
spectra_height_rel = spectra_height / fig_height
vertical_gap_rel = vertical_gap / fig_height
horizontal_gap_rel = horizontal_gap / fig_width
horizontal_padding_rel = [x / fig_width for x in padding[0]]
vertical_padding_rel = [x / fig_height for x in padding[1]]


# Make figure
fig = plt.figure(figsize=(fig_width, spectra_height))
ax = fig.add_axes(
    [
        horizontal_padding_rel[0],
        vertical_padding_rel[0],
        spectra_width_rel,
        spectra_height_rel,
    ]
)

# Plot data
data = pd.read_csv(
    os.path.join(data_dir, DP_data_folder, DP_file_name),
    header=1,
    skiprows=[2, 3],
    sep="\t",
)

rate = sputter_rates[3]

x_values = data[data.columns[0]][2:].values
for i in range(1, 478):
    y_values = data[data.columns[i]][2:].values
    max_y_value = y_values.max()
    norm_y_values = y_values / max_y_value
    # label = re.sub(r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species))
    ax.scatter(x_values * rate, norm_y_values, s=2)

# index_max = np.where(y_values == max_y_value)[0][0]
# x_position = x_values[index_max]
# x_position_distance = x_position * rate
# ax.axvline(12.16, color="black", linewidth=0.75, dashes=(5, 5))


# for species in species_to_plot:
#     x_values = data[data.columns[0]][2:].values
#     y_values = data[species][2:].values
#     max_y_value = y_values.max()
#     norm_y_values = y_values / max_y_value
#     label = re.sub(r"\-", r"$^{-}$", re.sub(r"(\^|_)(\d+)", r"$\1{\2}$", species))
#     ax.scatter(x_values * rate, norm_y_values, s=2, label=label)
#     max_y_value = max(y_values)
#     index_max = np.where(y_values == max_y_value)[0][0]
#     x_position = x_values[index_max]
#     x_position_distance = x_position * rate
#     ax.axvline(12.16, color="black", linewidth=0.75, dashes=(5, 5))
# if species < len(species_to_plot) - 1:
#     labels = [item.get_text() for item in ax.get_xticklabels()]
#     empty_string_labels = [""] * len(labels)
#     ax.set_xticklabels(empty_string_labels)
labels = [item.get_text() for item in ax.get_yticklabels()]
empty_string_labels = [""] * len(labels)
ax.set_xlim(x_range)
ax.set_ylim(y_range)
# if len(species_to_plot[species]) > 3:
#     num_col = 2
# else:
#     num_col = 1
ax.legend(
    loc="upper right",
    bbox_to_anchor=(1 - (0.05 / mass_spectra_ar), 0.95),
    fancybox=False,
    prop={"size": 8},
    frameon=True,
    ncol=1,
    markerscale=2,
)

ax.set_xlabel("Distance / nm")
ax.set_ylabel("Intensity / 10$^3$ counts")

fig.savefig(os.path.join(file_dir, output_file_name))
