# def create_matrix(delta_x):
#
#     length = int(1 / delta_x + 1)
#     mesh = np.zeros(shape=(length, length))
#     no_nodes = length * length
#
#     for row in range(length):
#         for column in range(length):
#             mesh[row, column] = length * (length - row - 1) + column
#     # (length * (length - row)) - (length - column)
#
#     return mesh

# arr = np.array([7, 4])
# arr2 = np.array([3, 6])
#
# result = np.array([arr, arr2])
#
# print(result)
# final = np.flip(result)
# final = np.fliplr(final)
# print(final)

import numpy as np
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# def x_y(delta_x):
#     columns = int(1 / delta_x) + 1
#     x, y = np.meshgrid(np.linspace(0, 1, columns), np.linspace(0, 1, columns))
#
#     return x, y
#
#
# coords = x_y(1 / 3)
#
# plt.scatter(coords[0], coords[1])
#
# segs1 = np.stack((coords[0], coords[1]), axis=2)
# print(segs1)
# segs2 = segs1.transpose(1, 0, 2)
# plt.gca().add_collection(LineCollection(segs1, color='#000000'))
#
# plt.show()

# A = np.array([15, 25, 60])
# B = np.sum(A)
# print(B)


rand_dict = {1: 5.0, 2: 9.0, 128: 213.0, 23: 1230.0}
print(rand_dict)
print(rand_dict.keys())
print(rand_dict.values())
