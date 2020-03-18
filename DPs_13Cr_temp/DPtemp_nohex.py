import matplotlib.font_manager as fm
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import os
import numpy as np
import pandas as pd
import string

file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")
output_file_name = "test.pdf"
row_labels = ["RT", "100$^o$C", "200$^o$C", "300$^o$C", "400$^o$C"]
image_labels = ["Position 1", "Position 2"]
sputter_rates = [0.11, 0.11, 0.08, 0.08, 0.07]

DP_data_folder = "depth_profiles_temp"
DP_file_names = [
    "a14704_4.txt",
    "a14705_1.txt",
    "a14707_1.txt",
    "a14708_3.txt",
    "a14710_2.txt",
]
DP_species_to_plot = ["FeO-", "CrO-"]
x_range = [0, 100]
y_range = [0, 1.2]


def colours(species):
    if species == "FeO-":
        color = "blue"
    else:
        color = "orange"
    return color


images_data_folder = "images_at_depth"
FeO_images_file_name = [
    "a14704_5_FeOCrO_BO.bmp",
    "a14705_21_FeOCrO_BO.bmp",
    "a14707_2_FeOCrO_BO.bmp",
    "a14708_4_FeOCrO_BO.bmp",
    "a14710_3_FeOCrO_BO.bmp",
]
CrO_images_file_name = [
    "a14704_6_FeOCrO_BO.bmp",
    "a14705_22_FeOCrO_BO.bmp",
    "a14707_3_FeOCrO_BO.bmp",
    "a14708_5_FeOCrO_BO.bmp",
    "a14710_4_FeOCrO_BO.bmp",
]
images_pixel_per_units = 192 / 100
images_scalebar_sizes = 20
images_scalebar_units = "$\mu$m"


fig_width = 7
depth_profile_width_proportion = 1 / 2
image_ar = 1
padding = [[0.5, 0.3], [0.5, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.02

depth_profile_width = fig_width * depth_profile_width_proportion
image_width = (
    fig_width - depth_profile_width - 2 * horizontal_gap - sum(padding[0])
) / 2
image_height = image_width / image_ar

fig_height = (
    len(CrO_images_file_name) * image_height
    + (len(CrO_images_file_name) - 1) * vertical_gap
    + sum(padding[1])
)
depth_profile_height = image_height

# Calculate relative dimensions
depth_profile_width_rel = depth_profile_width / fig_width
depth_profile_height_rel = depth_profile_height / fig_height
image_width_rel = image_width / fig_width
image_height_rel = image_height / fig_height
vertical_gap_rel = vertical_gap / fig_height
horizontal_gap_rel = horizontal_gap / fig_width
horizontal_padding_rel = [x / fig_width for x in padding[0]]
vertical_padding_rel = [x / fig_height for x in padding[1]]

# Prepare axis
fig = plt.figure(figsize=(fig_width, fig_height))
all_ax = np.array(
    [
        [
            fig.add_axes(
                [
                    horizontal_padding_rel[0],
                    1
                    - vertical_padding_rel[1]
                    - image_height_rel
                    - i * (image_height_rel + vertical_gap_rel),
                    depth_profile_width_rel,
                    depth_profile_height_rel,
                ]
            ),
            fig.add_axes(
                [
                    horizontal_padding_rel[0]
                    + depth_profile_width_rel
                    + horizontal_gap_rel,
                    1
                    - vertical_padding_rel[1]
                    - image_height_rel
                    - i * (image_height_rel + vertical_gap_rel),
                    image_width_rel,
                    image_height_rel,
                ]
            ),
            fig.add_axes(
                [
                    horizontal_padding_rel[0]
                    + depth_profile_width_rel
                    + 2 * horizontal_gap_rel
                    + image_width_rel,
                    1
                    - vertical_padding_rel[1]
                    - image_height_rel
                    - i * (image_height_rel + vertical_gap_rel),
                    image_width_rel,
                    image_height_rel,
                ]
            ),
        ]
        for i in range(len(CrO_images_file_name))
    ]
)

# Add depth profile


for i, files in enumerate(DP_file_names):
    data = pd.read_csv(
        os.path.join(data_dir, DP_data_folder, files),
        header=1,
        skiprows=[2, 3],
        sep="\t",
    )

    rate = sputter_rates[i]

    for species in DP_species_to_plot:
        x_values = data[data.columns[0]][2:].values
        y_values = data[species][2:].values
        x_values = data[data.columns[0]][2:].values
        y_values = data[species][2:].values
        max_y_value = y_values.max()
        norm_y_value = y_values / max_y_value
        all_ax[i, 0].scatter(x_values * rate, norm_y_value, s=2, color=colours(species))
        all_ax[i, 0].set_xlim(x_range)
        all_ax[i, 0].set_ylim(y_range)
        max_y_value = max(y_values)
        index_max = np.where(y_values == max_y_value)[0][0]
        x_position = x_values[index_max]
        x_position_distance = x_position * rate
        all_ax[i, 0].plot(
            [x_position_distance, x_position_distance],
            [-1, 4],
            color="black",
            linewidth=0.75,
        )
        if species == "FeO-":
            line_label = 1
        else:
            line_label = 2
        all_ax[i, 0].text(
            x_position_distance + 1,
            1.1,
            line_label,
            ha="center",
            va="center",
            fontproperties=fm.FontProperties(size=6),
        )
    if i < len(DP_file_names) - 1:
        labels = [item.get_text() for item in all_ax[i, 0].get_xticklabels()]
        empty_string_labels = [""] * len(labels)
        all_ax[i, 0].set_xticklabels(empty_string_labels)
    labels = [item.get_text() for item in all_ax[i, 0].get_yticklabels()]
    empty_string_labels = [""] * len(labels)
    all_ax[i, 0].set_yticklabels(empty_string_labels)


fig.text(
    horizontal_padding_rel[0] + 0.5 * depth_profile_width_rel,
    0.3 * vertical_padding_rel[0],
    "Distance / nm",
    ha="center",
    va="center",
)

fig.text(
    horizontal_padding_rel[0] * 0.5,
    vertical_padding_rel[0] + 2.5 * depth_profile_height_rel + 2 * vertical_gap_rel,
    "Normalised intensity",
    ha="center",
    va="center",
    rotation=90,
)


# Add images
images = [
    plt.imread(os.path.join(data_dir, images_data_folder, filename))
    for filename in [*FeO_images_file_name, *CrO_images_file_name]
]

for j in range(2):
    for i in range(len(FeO_images_file_name)):
        ax = all_ax[i, 1 + j]
        ax.imshow(images[len(FeO_images_file_name) * j + i], vmin=0, vmax=1)
        ax.axis("off")
        scalebar = AnchoredSizeBar(
            ax.transData,
            images_scalebar_sizes * images_pixel_per_units,
            str(images_scalebar_sizes) + " " + images_scalebar_units,
            "lower right",
            pad=0.1,
            color="white",
            frameon=False,
            size_vertical=1,
            fontproperties=fm.FontProperties(size=8),
        )
        ax.add_artist(scalebar)


# Add subplot labels
for j in range(3):
    for i in range(5):
        if j == 0:
            color = "black"
            ax_width = depth_profile_width_rel
            ax_height = depth_profile_height_rel
            x_position = 1 - letter_padding / ax_width
        else:
            color = "white"
            ax_width = image_width_rel
            ax_height = image_height_rel
            x_position = letter_padding / ax_width
        all_ax[i, j].text(
            x_position,
            1 - letter_padding / ax_height,
            string.ascii_lowercase[i * 3 + j],
            transform=all_ax[i, j].transAxes,
            ha="center",
            va="center",
            color=color,
        )

for i, text in enumerate(row_labels):
    fig.text(
        1 - 0.5 * horizontal_padding_rel[1],
        1
        - (
            vertical_padding_rel[1]
            + (i + 0.5) * image_height_rel
            + i * vertical_gap_rel
        ),
        text,
        ha="center",
        va="center",
        rotation=270,
        fontproperties=fm.FontProperties(size=10),
    )

fig.text(
    horizontal_padding_rel[0] + 0.5 * depth_profile_width_rel,
    1 - (0.5 * vertical_padding_rel[1]),
    "Depth Profiles",
    ha="center",
    va="center",
    rotation=0,
    fontproperties=fm.FontProperties(size=10),
)

for i, text in enumerate(image_labels):
    fig.text(
        horizontal_padding_rel[0]
        + depth_profile_width_rel
        + horizontal_gap_rel
        + (i + 0.5) * image_width_rel
        + i * horizontal_gap_rel,
        1 - (0.5 * vertical_padding_rel[1]),
        text,
        ha="center",
        va="center",
        rotation=0,
        fontproperties=fm.FontProperties(size=10),
    )

fig.savefig(os.path.join(file_dir, output_file_name))
