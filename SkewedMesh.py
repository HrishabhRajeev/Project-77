import csv

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from matplotlib.collections import LineCollection
import seaborn as sns


class Elements:

    def __init__(self, delta_x, alpha=0.0):
        self.MeshSpacing = delta_x
        self.Domain = 1
        self.Row_Col_Length = int(1 / delta_x + 1)
        self.cell_size = self.Row_Col_Length - 1
        self.NumberOfElementCenters = int(np.power((1 / delta_x), 2))
        self.NumberOfVertices = int(np.power((1 / delta_x + 1), 2))
        self.alpha = alpha
        self.z_vector = np.array([0, 0, 1])

    @property
    def squareVertices(self):
        mesh_spacing = self.MeshSpacing
        row_col_length = self.Row_Col_Length
        vertices = np.zeros(shape=(row_col_length, row_col_length, 2))

        for row in range(row_col_length):
            for column in range(row_col_length):
                vertices[row][column] = np.array([column * mesh_spacing, row * mesh_spacing])

        return vertices

    @property
    def skewedVertices(self):
        mesh_spacing = self.MeshSpacing
        vertices = self.squareVertices
        row_col_length = self.Row_Col_Length
        alpha = self.alpha
        for row_counter, row in enumerate(vertices[1:row_col_length]):
            currentrow = row_counter + 1
            columnccounter = 0

            for columnccounter, column in enumerate(row[1:row_col_length - 1]):
                currentcolumn = columnccounter + 1

                if (currentrow + 1) % 2 == 0:
                    vertices[currentrow][currentcolumn][0] = column[0] - (alpha * mesh_spacing)

        return vertices

    @property
    def internal_vertex_coordinates(self):

        length = self.Row_Col_Length
        vertices = self.internal_vertices[0]
        coords = np.flipud(self.skewedVertices)
        array_element_coordinates = coords[1:(length - 1), 1:(length - 1)]
        element_coordinates = {}

        for row, vertex_row in enumerate(vertices):
            for column, vertex in enumerate(vertex_row):
                element_coordinates[vertex] = array_element_coordinates[row][column]

        return element_coordinates

    @property
    def nodes(self):

        vertices = self.skewedVertices
        row_col_length = self.Row_Col_Length
        nodes = {}

        for row_number, row_values in enumerate(vertices):
            for column, coordinate in enumerate(row_values):
                index = column + (row_number * row_col_length)
                nodes[index] = coordinate

        return nodes

    @property
    def node_connect(self):

        nodes = self.nodes
        row_col_length_1 = self.Row_Col_Length - 1
        row_col_length = self.Row_Col_Length
        node_connect = {}

        for row in range(row_col_length_1):
            for col in range(row_col_length_1):
                current_element_center = col + (row * row_col_length_1)

                sw = current_element_center + row
                se = current_element_center + row + 1
                ne = current_element_center + row_col_length + 1 + row
                nw = current_element_center + row_col_length + row

                node_connect[current_element_center] = np.array([nodes[sw], nodes[se], nodes[ne], nodes[nw]])

        return node_connect

    @property
    def element_centers(self):

        connected_nodes = self.node_connect
        centroid = {}
        for key in connected_nodes:
            polygon = Polygon(connected_nodes[key])
            centroid[key] = np.array(polygon.centroid.coords).flatten()

        return centroid

    def drawGraph(self):

        vertices = self.skewedVertices  # this is what you need to change to plot the new skewed graph
        unadjusted_vertices = self.squareVertices
        row_col_length = self.Row_Col_Length
        number_of_vertices = self.NumberOfVertices
        edge_centers = self.edge_connect
        element_centers = self.element_centers
        element_center_dict = self.element_center_connect
        hybrid_volume_centroids = self.hybrid_volume[1]
        standard_volume_centroids = self.standard_volume[1]
        x_coords_temp = []
        y_coords_temp = []
        centroid_x_temp = []
        centroid_y_temp = []
        hybrid_volume_centroids_x_temp = []
        hybrid_volume_centroids_y_temp = []
        standard_volume_centroids_x_temp = []
        standard_volume_centroids_y_temp = []
        squarey = []
        squarex = []
        edge_centers_x = []
        edge_centers_y = []
        for entry in vertices:
            x_coords_temp.append(entry[:, 0])
            y_coords_temp.append(entry[:, 1])

        for key in element_centers:
            centroid_x_temp.append(element_centers[key][0])
            centroid_y_temp.append(element_centers[key][1])

        for enter in unadjusted_vertices:
            squarey.append(enter[:, 1])
            squarex.append(enter[:, 0])

        for key in edge_centers:
            edge_centers_x.append(edge_centers[key][:, 0])
            edge_centers_y.append(edge_centers[key][:, 1])
            hybrid_volume_centroids_x_temp.append((hybrid_volume_centroids[key][0]))
            hybrid_volume_centroids_y_temp.append((hybrid_volume_centroids[key][1]))
            standard_volume_centroids_x_temp.append((standard_volume_centroids[key][0]))
            standard_volume_centroids_y_temp.append((standard_volume_centroids[key][1]))

        edge_centers_y = np.array(edge_centers_y)
        edge_centers_x = np.array(edge_centers_x)
        x_coords = np.array(x_coords_temp)
        y_coords = np.array(y_coords_temp)
        centroid_y = np.array(centroid_y_temp)
        centroid_x = np.array(centroid_x_temp)
        hybrid_volume_centroids_x = np.array(hybrid_volume_centroids_x_temp)
        hybrid_volume_centroids_y = np.array(hybrid_volume_centroids_y_temp)
        standard_volume_centroids_x = np.array(standard_volume_centroids_x_temp)
        standard_volume_centroids_y = np.array(standard_volume_centroids_y_temp)
        squarey = np.array(squarey)
        squarex = np.array(squarex)

        fig = plt.figure(figsize=plt.figaspect(0.5))

        axs1 = fig.add_subplot(1, 2, 1)
        axs2 = fig.add_subplot(1, 2, 2)

        axs1.scatter(x_coords, y_coords, color='black')
        axs1.scatter(centroid_x, centroid_y, marker='x', color='r')
        axs1.scatter(standard_volume_centroids_x, standard_volume_centroids_y, marker='D', color='#FFC000')
        axs1.scatter(edge_centers_x, edge_centers_y, facecolor='none', edgecolors='blue')

        axs2.scatter(x_coords, y_coords, color='black')
        axs2.scatter(centroid_x, centroid_y, marker='x', color='r')
        axs2.scatter(hybrid_volume_centroids_x, hybrid_volume_centroids_y, marker='D', color='#FFC000')
        axs2.scatter(edge_centers_x, edge_centers_y, facecolor='none', edgecolors='blue')

        for index in element_center_dict:
            x_coords_element_centers = element_center_dict[index][:, 0]
            y_coords_element_centers = element_center_dict[index][:, 1]
            x_coords_edge_centers = edge_centers[index][:, 0]
            y_coords_edge_centers = edge_centers[index][:, 1]

            x_path = np.stack((x_coords_element_centers, x_coords_edge_centers), axis=1).reshape(8)
            y_path = np.stack((y_coords_element_centers, y_coords_edge_centers), axis=1).reshape(8)
            x_path = np.append(x_path, x_coords_element_centers[0])
            y_path = np.append(y_path, y_coords_element_centers[0])
            axs1.plot(x_path, y_path, color='black', alpha=0.5, linestyle=':')
            axs1.fill(x_path, y_path, facecolor='#061A40', alpha=0.1)

            axs2.plot(x_coords_element_centers, y_coords_element_centers, linestyle=':', color='black')
            axs2.fill(x_coords_element_centers, y_coords_element_centers, facecolor='#061A40', alpha=0.1)
        axs2.plot([element_center_dict[6][0][0],
                   element_center_dict[16][3][0]],
                  [element_center_dict[6][0][1],
                   element_center_dict[16][3][1]], color='black', linestyle=':')

        axs1.plot(x_coords, y_coords, color='black')
        axs1.plot(squarey, squarex, color='black')

        axs2.plot(x_coords, y_coords, color='black')
        axs2.plot(squarey, squarex, color='black')

        return

    @property
    def internal_vertices(self):

        number_of_vertices = self.NumberOfVertices
        row_col_length = self.Row_Col_Length
        vertices = np.zeros(number_of_vertices)
        for vertex in range(number_of_vertices):
            vertices[vertex] = vertex

        vertices = np.flipud(np.reshape(vertices, (row_col_length, row_col_length)))
        internal_vertices = vertices[1:(row_col_length - 1), 1:(row_col_length - 1)]
        internal_vertices = np.array(internal_vertices)

        return internal_vertices, vertices

    @property
    def edge_connect(self):

        edge_connect = {}
        internal_vertices = self.internal_vertices[0]
        vertex_coordinates = self.skewedVertices
        vertex_coordinates = np.flipud(vertex_coordinates)

        for row_count, row in enumerate(internal_vertices):
            for column_count, column in enumerate(row):
                dict_key = internal_vertices[row_count][column_count]

                current_node_coordinate = vertex_coordinates[row_count + 1][column_count + 1]

                n = (current_node_coordinate + vertex_coordinates[row_count][column_count + 1]) / 2
                e = (current_node_coordinate + vertex_coordinates[row_count + 1][column_count + 2]) / 2
                s = (current_node_coordinate + vertex_coordinates[row_count + 2][column_count + 1]) / 2
                w = (current_node_coordinate + vertex_coordinates[row_count + 1][column_count]) / 2

                edge_connect[dict_key] = np.array([s, e, n, w])

        return edge_connect

    @property
    def element_center_connect(self):

        internal_vertices = self.internal_vertices[0]
        element_center_coordinates = self.element_centers
        cell_size = self.Row_Col_Length - 1
        element_center_connect = {}

        for count, row in enumerate(internal_vertices):
            for col in row:
                sw = int(col) - (cell_size * 2) + count
                se = int(col) - (cell_size * 2) + count + 1
                ne = int(col) - cell_size + 1 + count
                nw = int(col) - cell_size + count

                element_center_connect[col] = np.array([element_center_coordinates[sw],
                                                        element_center_coordinates[se],
                                                        element_center_coordinates[ne],
                                                        element_center_coordinates[nw]])

        return element_center_connect

    @property
    def hybrid_volume(self):

        element_centers_per_node = self.element_center_connect
        volume = {}
        volume_centroid = {}
        for key in element_centers_per_node:
            polygon = Polygon(element_centers_per_node[key])
            volume[key] = polygon.area
            volume_centroid[key] = np.array(polygon.centroid.coords).flatten()

        return volume, volume_centroid

    @property
    def standard_volume(self):

        element_centers_per_node = self.element_center_connect
        edge_centers_per_node = self.edge_connect
        volume = {}
        volume_centroid = {}

        for key in element_centers_per_node:
            polygon_coordinates = np.stack((element_centers_per_node[key],
                                            edge_centers_per_node[key]),
                                           axis=1).reshape(8, 2)
            polygon = Polygon(polygon_coordinates)
            volume[key] = polygon.area
            volume_centroid[key] = np.array(polygon.centroid.coords).flatten()

        return volume, volume_centroid

    # this bit of code generates the graphs and plots them for three separate meshes and alpha values ------------------

    @property
    def dphi_dxdy_hybrid(self):

        node_connect = self.element_center_connect
        edge_connect = self.edge_connect

        internal_vertex_size = int(np.power((self.cell_size - 1), 2))
        z_vector = self.z_vector
        volume = self.hybrid_volume[0]
        dphidx = np.zeros(shape=internal_vertex_size)
        dphidy = np.zeros(shape=internal_vertex_size)
        Cfx = np.zeros(shape=internal_vertex_size)
        Cfy = np.zeros(shape=internal_vertex_size)

        for node_number, node in enumerate(node_connect):
            coords = node_connect[node]
            edge_coords = edge_connect[node]

            for edge in range(4):

                if edge == 3:
                    vector_l = coords[0] - coords[edge]

                else:
                    vector_l = coords[edge + 1] - coords[edge]

                normal = np.cross(vector_l, z_vector)
                unit_normal_x = (normal / np.linalg.norm(normal))[0]
                unit_normal_y = (normal / np.linalg.norm(normal))[1]
                area = np.linalg.norm(vector_l)

                Cfx[node_number] += (area * unit_normal_x) * np.sin(edge_coords[edge][0])
                Cfy[node_number] += (area * unit_normal_y) * np.cos(edge_coords[edge][1]) * 1.5

            dphidx[node_number] = Cfx[node_number] / volume[node]
            dphidy[node_number] = Cfy[node_number] / volume[node]
        return dphidx, dphidy

    @property
    def dphi_dxdy_hybrid_plus(self):

        node_connect = self.element_center_connect
        edge_connect = self.edge_connect

        internal_vertex_size = int(np.power((self.cell_size - 1), 2))
        z_vector = self.z_vector
        volume = self.hybrid_volume[0]
        volume_centroid = self.hybrid_volume[1]
        dphidx = np.zeros(shape=internal_vertex_size)
        dphidy = np.zeros(shape=internal_vertex_size)
        Cfx = np.zeros(shape=internal_vertex_size)
        Cfy = np.zeros(shape=internal_vertex_size)

        for node_number, node in enumerate(node_connect):
            coords = node_connect[node]
            volume_centroid_coord = volume_centroid[node]
            edgeee = edge_connect[node]
            for edge in range(4):

                if edge == 3:
                    vector_l = coords[0] - coords[edge]


                else:
                    vector_l = coords[edge + 1] - coords[edge]
                    edge_coords = (coords[edge + 1] + coords[edge]) / 2

                normal = np.cross(vector_l, z_vector)
                unit_normal_x = (normal / np.linalg.norm(normal))[0]
                unit_normal_y = (normal / np.linalg.norm(normal))[1]
                area = np.linalg.norm(vector_l)

                Cfx[node_number] += (area * unit_normal_x) * np.sin(volume_centroid_coord[0])
                Cfy[node_number] += (area * unit_normal_y) * np.cos(volume_centroid_coord[1]) * 1.5
                # print('edgecenter:', edgeee[edge])
                # print('edgecoord:', edge_coords)
                # print('volumecoord:', volume_centroid_coord)
                # print('\n\n')
            dphidx[node_number] = (Cfx[node_number] / volume[node])
            dphidy[node_number] = (Cfy[node_number] / volume[node])
        return dphidx, dphidy

    @property
    def dphi_dxdy_standard(self):

        node_connect = self.element_center_connect
        edge_connect = self.edge_connect

        internal_vertex_size = int(np.power((self.cell_size - 1), 2))
        z_vector = self.z_vector
        volume = self.standard_volume[0]
        dphidx = np.zeros(shape=internal_vertex_size)
        dphidy = np.zeros(shape=internal_vertex_size)
        Cfx = np.zeros(shape=internal_vertex_size)
        Cfy = np.zeros(shape=internal_vertex_size)

        index_edge_center = np.array([0, 0, 1, 1, 2, 2, 3, 3])
        index_element_center = np.array([0, 1, 1, 2, 2, 3, 3, 0])

        for vertex_number, vertex in enumerate(node_connect):

            edge_center_coordinates = edge_connect[vertex]
            element_center_coordinates = node_connect[vertex]

            for i in range(8):
                l_vector = (edge_center_coordinates[index_edge_center[i]]
                            - element_center_coordinates[index_element_center[i]]) * np.power(-1, i)
                normal = np.cross(l_vector, z_vector)
                unit_normal_x = (normal / np.linalg.norm(normal))[0]
                unit_normal_y = (normal / np.linalg.norm(normal))[1]
                area = np.linalg.norm(l_vector)

                Cfx[vertex_number] += (area * unit_normal_x) * np.sin(edge_center_coordinates
                                                                      [index_edge_center[i]][0])
                Cfy[vertex_number] += (area * unit_normal_y) * 1.5 * np.cos(edge_center_coordinates
                                                                            [index_edge_center[i]][1])

            dphidx[vertex_number] = Cfx[vertex_number] / volume[vertex]
            dphidy[vertex_number] = Cfy[vertex_number] / volume[vertex]

        return dphidx, dphidy

    @property
    def hybrid_analytical(self):

        internal_vertex_size = int(np.power((self.cell_size - 1), 2))
        dphidx_analytical = np.zeros(shape=internal_vertex_size)
        dphidy_analytical = np.zeros(shape=internal_vertex_size)
        internal_vertex_coordinates = self.hybrid_volume[1]

        for counter, coordinates in enumerate(internal_vertex_coordinates):
            x_coordiate = internal_vertex_coordinates[coordinates][0]
            y_coordinate = internal_vertex_coordinates[coordinates][1]
            dphidx_analytical[counter] = np.cos(x_coordiate)
            dphidy_analytical[counter] = -1.5 * np.sin(y_coordinate)

        return dphidx_analytical, dphidy_analytical

    @property
    def standard_analytical(self):

        internal_vertex_size = int(np.power((self.cell_size - 1), 2))
        dphidx_analytical = np.zeros(shape=internal_vertex_size)
        dphidy_analytical = np.zeros(shape=internal_vertex_size)
        internal_vertex_coordinates = self.standard_volume[1]

        for counter, coordinates in enumerate(internal_vertex_coordinates):
            x_coordiate = internal_vertex_coordinates[coordinates][0]
            y_coordinate = internal_vertex_coordinates[coordinates][1]
            dphidx_analytical[counter] = np.cos(x_coordiate)
            dphidy_analytical[counter] = -1.5 * np.sin(y_coordinate)

        return dphidx_analytical, dphidy_analytical

    @property
    def error_dxdy_hybrid(self):

        dphidxdy = self.dphi_dxdy_hybrid
        analytical_dphi = self.hybrid_analytical
        volume = self.hybrid_volume[0]
        volume_array = np.array(list(volume.values()))

        phi_difference_x = dphidxdy[0] - analytical_dphi[0]
        phi_difference_y = dphidxdy[1] - analytical_dphi[1]
        phi_difference_x_absolute = np.abs(phi_difference_x)
        phi_difference_y_absolute = np.abs(phi_difference_y)

        numerator_x = np.dot(phi_difference_x_absolute, volume_array)
        numerator_y = np.dot(phi_difference_y_absolute, volume_array)
        avg_grad = (np.sin(1) + (1.5 * np.cos(0))) / (1 - 0)
        denominator = np.sum(volume_array)

        error_metric_dx = numerator_x / (denominator * avg_grad)
        error_metric_dy = numerator_y / (denominator * avg_grad)

        return error_metric_dx, error_metric_dy

    @property
    def error_dxdy_standard(self):
        dhpidxdy = self.dphi_dxdy_standard
        analytical_dphi = self.standard_analytical
        volume = self.standard_volume[0]
        volume_array = np.array(list(volume.values()))

        phi_difference_x = dhpidxdy[0] - analytical_dphi[0]
        phi_difference_y = dhpidxdy[1] - analytical_dphi[1]
        phi_difference_x_absolute = np.abs(phi_difference_x)
        phi_difference_y_absolute = np.abs(phi_difference_y)

        numerator_x = np.dot(phi_difference_x_absolute, volume_array)
        numerator_y = np.dot(phi_difference_y_absolute, volume_array)
        avg_grad = (np.sin(1) + (1.5 * np.cos(0))) / (1 - 0)
        denominator = np.sum(volume_array)

        error_metric_dx = numerator_x / (denominator * avg_grad)
        error_metric_dy = numerator_y / (denominator * avg_grad)

        return error_metric_dx, error_metric_dy

    def heatmap(self):

        dphi_dxdy_hybrid = self.dphi_dxdy_hybrid_plus
        dphi_dxdy_standard = self.dphi_dxdy_standard
        dphi_dx_hybrid = dphi_dxdy_hybrid[0]
        dphi_dy_hybrid = dphi_dxdy_hybrid[1]
        dphi_dx_standard = dphi_dxdy_standard[0]
        dphi_dy_standard = dphi_dxdy_standard[1]
        internal_vertex_length = int(self.cell_size - 1)
        phi_val_dictionary = {}

        avg = self.internal_vertices
        internal_vertex_coordinates = self.internal_vertex_coordinates
        for key in internal_vertex_coordinates:
            phi_val_dictionary[key] = np.sin(internal_vertex_coordinates[key][0]) + \
                                      (1.5 * np.cos(internal_vertex_coordinates[key][1]))

        phi_val = np.array(list(phi_val_dictionary.values()))
        phi_val = np.reshape(phi_val, (internal_vertex_length, internal_vertex_length))
        # print(phi_val_dictionary,'\n\n',phi_val)

        x = np.linspace(0, 1, 10)
        y = np.linspace(0, 1, 10)
        X, Y = np.meshgrid(x, y)
        phi1 = np.sin(X) + (1.5 * np.cos(Y))

        x2 = np.linspace(0, 1, 20)
        y2 = np.linspace(0, 1, 20)
        X2, Y2 = np.meshgrid(x2, y2)
        phi2 = np.sin(X2) + (1.5 * np.cos(Y2))

        x3 = np.linspace(0, 1, 40)
        y3 = np.linspace(0, 1, 40)
        X3, Y3 = np.meshgrid(x3, y3)
        phi3 = np.sin(X3) + (1.5 * np.cos(Y3))

        x4 = np.linspace(0, 1, 80)
        y4 = np.linspace(0, 1, 80)
        X4, Y4 = np.meshgrid(x4, y4)
        phi4 = np.sin(X4) + (1.5 * np.cos(Y4))

        fig = plt.figure(figsize=plt.figaspect(0.5))

        axs = fig.add_subplot(2, 2, 1, projection='3d')
        surface1 = axs.plot_surface(X, Y, phi1, rstride=1, cstride=1,
                                    cmap='viridis', edgecolor='none')
        fig.colorbar(surface1, shrink=0.5, aspect=20)

        axs = fig.add_subplot(2, 2, 2, projection='3d')
        surface2 = axs.plot_surface(X2, Y2, phi2, rstride=1, cstride=1,
                                    cmap='viridis', edgecolor='none')
        fig.colorbar(surface2, shrink=0.5, aspect=20)

        axs = fig.add_subplot(2, 2, 3, projection='3d')
        surface3 = axs.plot_surface(X3, Y3, phi3, rstride=1, cstride=1,
                                    cmap='viridis', edgecolor='none')
        fig.colorbar(surface3, shrink=0.5, aspect=20)

        axs = fig.add_subplot(2, 2, 4, projection='3d')
        surface4 = axs.plot_surface(X4, Y4, phi4, rstride=1, cstride=1,
                                    cmap='viridis', edgecolor='none')
        fig.colorbar(surface4, shrink=0.5, aspect=20)

        # axs = fig.add_subplot(1, 3, 2)
        # line = axs.plot(dx, self.analytical_standard[0])

        # axs = fig.add_subplot(1, 3, 3)
        # line = axs.plot(dx, self.analytical_standard[1])
        # axs = fig.add_subplot(2, 2, 2, projection='3d')
        # surface = axs.plot_surface(X2, Y2, phi2, rstride=1, cstride=1,
        #                            cmap='viridis', edgecolor='none')
        #
        # axs = fig.add_subplot(2, 2, 3, projection='3d')
        # surface = axs.plot_surface(X3, Y3, phi3, rstride=1, cstride=1,
        #                            cmap='viridis', edgecolor='none')
        #
        # axs = fig.add_subplot(2, 2, 4, projection='3d')
        # surface = axs.plot_surface(X4, Y4, phi4, rstride=10, cstride=10)

        plt.show()
        return

# mesh1 = Elements(1 / 4, alpha=0.26)
# mesh3 = Elements(1 / 4, alpha=0.46)
# print(mesh1.node_connect)
# print('\n\n', mesh1.element_centers)
# print('\n\n', mesh2.node_connect)
# print('\n\n', mesh2.element_centers)
# mesh1.drawGraph()
# print(mesh1.dphi_dxdy_hybrid[0])
# print(mesh1.edge_connect)
# print('Hybrid vol:\n', mesh1.hybrid_volume, '\n\nStandard vol:\n', mesh1.standard_volume)
# print(mesh1.internal_vertex_coordinates)
# with open('Volumes.csv', 'w', encoding='UTF8') as w:
#
#     writer = csv.writer(w)
#     hybrid_vol = mesh1.hybrid_volume
#     standard_vol = mesh1.standard_volume
#     writer.writerow(hybrid_vol.keys())
#     writer.writerow(hybrid_vol.values())
#     writer.writerow(standard_vol.keys())
#     writer.writerow(standard_vol.values())
# mesh1.drawGraph()
# print(mesh1.dphi_dxdy_hybrid_plus)
# mesh3.drawGraph()
# plt.show()
# print(mesh1.node_connect)

# ----------------------------------------------------------------------------------------------------------------------
