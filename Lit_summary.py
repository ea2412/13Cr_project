from matplotlib import pyplot as plt

temperature = [
    50,
    50,
    157,
    164,
    275,
    300,
    445,
    481,
    484,
    726,
    758,
    800,
    800,
    800,
    25,
    100,
    200,
    300,
    400,
]
composition = [20, 40, 9, 9, 40, 20, 40, 9, 9, 9, 9, 17, 19, 22, 13, 13, 13, 13, 13]
# blue is iron, orange is chromium, black is sublayer, purple is mix

surface_elemental_composition = [
    "blue",
    "purple",
    "blue",
    "black",
    "black",
    "orange",
    "orange",
    "black",
    "orange",
    "orange",
    "orange",
    "orange",
    "orange",
    "orange",
    "black",
    "black",
    "black",
    "black",
    "black",
]


fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_ylabel("Temperature ($^o$C)")
ax1.set_title("Composition (at. %)")


ax1.scatter(
    composition, temperature, color=surface_elemental_composition,
)

plt.savefig("Litreview_surface.pdf")
