import os
import string
import sys
import numpy as np
import pandas as pd
import matplotlib.font_manager as fm
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

sys.path.append(os.path.abspath(os.path.join(__file__, os.path.pardir, os.path.pardir)))
from Colors import ColorSet

file_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(file_dir), "data")

output_file_name = "slices_test.pdf"

slices_data_folder = "slices"
images_data_folder = "images_at_depth"
axis_data_folder = "axis_schematic"

number_of_rows = 4

Fe_slices_file_name = [
    "a14704_4_FeOslice1.bmp",
    "a14705_1_FeOslice1.bmp",
    "a14707_1_FeOslice1.bmp",
    "a14708_3_FeOslice1.bmp",
    "a14710_2FeOslice1.bmp",
]
Cr_slices_file_name = [
    "a14704_4_CrOslice1.bmp",
    "a14705_1_CrOslice1.bmp",
    "a14707_1_CrOslice1.bmp",
    "a14708_3_CrOslice1.bmp",
    "a14710_2_CrOslice1.bmp",
]

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


scale_bar_file_name = "slice example scale bar.bmp"
cube_pos1_file_name = "cubepos1.png"
cube_pos2_file_name = "cubepos2.png"
xz_schematic_file_name = "axisxz.png"
xy_schematic_file_name = "axisxy.png"
xyz_schematic_file_name = "axisxyz.png"

slices_list = [*Fe_slices_file_name, *Cr_slices_file_name]
images_list = [*FeO_images_file_name, *CrO_images_file_name]

images_pixel_per_unit = 192 / 100
images_scalebar_size = 20
images_scalebar_unit = "$\mu$m"

fig_width = 7
image_ar = 1
scale_bar_ar = 19 / 195
scale_bar_proportion = 1 / 20
scale_bar_width = scale_bar_proportion * fig_width
padding = [[0.3, 0.2], [0.2, 0.3]]  # padding = [[left, right], [bottom, top]]
horizontal_gap = 0.1
vertical_gap = 0.1
letter_padding = 0.18

image_width = (
    fig_width
    - sum(padding[0])
    - len(Fe_slices_file_name) * horizontal_gap
    - 2 * scale_bar_width
) / (len(Fe_slices_file_name) + 1)
image_height = image_width / image_ar
scale_bar_height = 2 * image_height + vertical_gap

fig_height = (
    sum(padding[1])
    + number_of_rows * image_height
    + vertical_gap * (number_of_rows - 1)
)

# calculate relative dimensions
image_width_rel = image_width / fig_width
image_height_rel = image_height / fig_height
scale_bar_width_rel = scale_bar_width / fig_width
scale_bar_height_rel = scale_bar_height / fig_height
axis_image_width_rel = image_width_rel / 2
axis_image_height_rel = image_height_rel / 2
vertical_gap_rel = vertical_gap / fig_height
horizontal_gap_rel = horizontal_gap / fig_width
horizontal_padding_rel = [x / fig_width for x in padding[0]]
vertical_padding_rel = [x / fig_height for x in padding[1]]

# Prepare axis
fig = plt.figure(figsize=(fig_width, fig_height))

axs = np.array(
    [
        [
            fig.add_axes(
                [
                    horizontal_padding_rel[0]
                    + i * (image_width_rel + horizontal_gap_rel),
                    1
                    - (
                        vertical_padding_rel[1]
                        + image_height_rel
                        # + vertical_gap_rel
                        + j * (image_height_rel + vertical_gap_rel)
                    ),
                    image_width_rel,
                    image_height_rel,
                ]
            )
            for i in range(len(Fe_slices_file_name))
        ]
        for j in range(number_of_rows)
    ]
)

scale_ax = fig.add_axes(
    [
        horizontal_padding_rel[0]
        + len(Fe_slices_file_name) * (image_width_rel + horizontal_gap_rel),
        1 - (vertical_padding_rel[1] + image_height_rel * 2 + vertical_gap_rel),
        scale_bar_width_rel,
        scale_bar_height_rel,
    ]
)


axis_images_axs = np.array(
    [
        fig.add_axes(
            [
                horizontal_padding_rel[0]
                + len(Fe_slices_file_name) * (image_width_rel + horizontal_gap_rel)
                + horizontal_gap_rel
                + axis_image_width_rel,
                1
                - (
                    vertical_padding_rel[1]
                    + image_height_rel * 0.5
                    + p * (image_height_rel + vertical_gap_rel)
                ),
                axis_image_width_rel,
                axis_image_height_rel,
            ]
        )
        for p in range(number_of_rows)
    ]
)

cube_axs = np.array(
    [
        fig.add_axes(
            [
                horizontal_padding_rel[0]
                + len(Fe_slices_file_name) * (image_width_rel + horizontal_gap_rel),
                1
                - (
                    vertical_padding_rel[1]
                    + (3 * image_height_rel)
                    + (2 * vertical_gap_rel)
                    + q * (image_height_rel + vertical_gap_rel)
                ),
                axis_image_width_rel,
                axis_image_height_rel,
            ]
        )
        for q in range(2)
    ]
)

slices_FeO = [
    np.divide(plt.imread(os.path.join(data_dir, slices_data_folder, filename)), 255)
    for filename in Fe_slices_file_name
]

slices_CrO = [
    np.divide(plt.imread(os.path.join(data_dir, slices_data_folder, filename)), 255)
    for filename in Cr_slices_file_name
]

images_FeO = [
    np.divide(plt.imread(os.path.join(data_dir, images_data_folder, filename)), 255)
    for filename in FeO_images_file_name
]

images_CrO = [
    np.divide(plt.imread(os.path.join(data_dir, images_data_folder, filename)), 255)
    for filename in CrO_images_file_name
]

scale_bar = plt.imread(os.path.join(data_dir, slices_data_folder, scale_bar_file_name))
scale_ax.imshow(scale_bar, vmin=0, vmax=1)
scale_ax.axis("off")

axis_schematic_list = [
    xz_schematic_file_name,
    xz_schematic_file_name,
    xy_schematic_file_name,
    xy_schematic_file_name,
]

cubepos1 = mpimg.imread(os.path.join(data_dir, axis_data_folder, cube_pos1_file_name))
cube_axs[0].imshow(cubepos1, vmin=0, vmax=1)
cube_axs[0].axis("off")

# for l in range(len(axis_schematic_list)):
#     ax = axis_images_axs[l]
#     axis_schematic = plt.imread(
#         os.path.join(data_dir, axis_data_folder, axis_schematic_list[l])
#     )
#     ax.imshow(axis_schematic, vmin=0, vmax=1)
#     # ax.axis("off")

combined_list = [slices_FeO, slices_CrO, images_FeO, images_CrO]

for i in range(number_of_rows):
    for j in range(len(FeO_images_file_name)):
        ax = axs[i, j]
        ax.imshow(combined_list[i][j], vmin=0, vmax=1)
        ax.axis("off")
        scalebar = AnchoredSizeBar(
            ax.transData,
            images_scalebar_size * images_pixel_per_unit,
            str(images_scalebar_size) + " " + images_scalebar_unit,
            "lower right",
            pad=0.1,
            color="white",
            frameon=False,
            size_vertical=1,
            fontproperties=fm.FontProperties(size=8),
        )
    ax.add_artist(scalebar)


# Add subplot labels
for i in range(number_of_rows):
    for j in range(len(FeO_images_file_name)):
        ax = axs[i, j]
        ax.text(
            letter_padding / image_width,
            letter_padding / image_height,
            string.ascii_lowercase[(i * 5) + j],
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="white",
        )

row_labels = ["FeO$^-$", "CrO$^-$", "Position 1", "Position 2"]
column_labels = ["RT", "100$^o$C", "200$^o$C", "300$^o$C", "400$^o$C"]

for i, text in enumerate(row_labels):
    fig.text(
        0.5 * horizontal_padding_rel[0],
        1
        - (
            vertical_padding_rel[1]
            + (i + 0.5) * image_height_rel
            + i * vertical_gap_rel
        ),
        text,
        ha="center",
        va="center",
        rotation=90,
        fontproperties=fm.FontProperties(size=10),
    )

for i, text in enumerate(column_labels):
    fig.text(
        horizontal_padding_rel[0]
        + (i + 0.5) * image_width_rel
        + i * horizontal_gap_rel,
        1 - (0.5 * vertical_padding_rel[1]),
        text,
        ha="center",
        va="center",
        rotation=0,
        fontproperties=fm.FontProperties(size=10),
    )

FeO_label = mpatches.Patch(color="blue", label="FeO$^-$")
CrO_label = mpatches.Patch(color="orange", label="CrO$^-$")


scale_ax.legend(
    handles=[FeO_label, CrO_label],
    loc="upper right",
    bbox_to_anchor=(3, -0.15),
    fancybox=False,
    prop={"size": 8},
    frameon=False,
    ncol=1,
    markerscale=2,
)

fig.savefig(os.path.join(file_dir, output_file_name))
