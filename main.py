# Author: Hrishabh Rajeev ----------------------------------------------------------------------------------------------
# Institution: University of Cape Town ---------------------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as ntx
import pandas as pd
import scipy as spy
from sympy import pretty_print as pp, latex
import SkewedMesh as skw
import shapely as shy
import csv


def main():
    class Elements:

        def __init__(self, delta_x):
            self.MeshSpacing = delta_x
            self.mesh_size = int(1 / delta_x + 1)
            self.cell_size = int(1 / delta_x)
            self.z_vector = np.array([0, 0, 1])
            self.alpha_0 = 0.0
            self.alpha_1 = 0.1
            self.alpha_2 = 0.2
            self.alpha_3 = 0.46
            self.alpha_4 = 0.05

            self.alpha_0_element = skw.SkewedElements(self.MeshSpacing, alpha=self.alpha_0)
            self.alpha_1_element = skw.SkewedElements(self.MeshSpacing, alpha=self.alpha_1)
            self.alpha_2_element = skw.SkewedElements(self.MeshSpacing, alpha=self.alpha_2)
            self.alpha_3_element = skw.SkewedElements(self.MeshSpacing, alpha=self.alpha_3)
            self.alpha_4_element = skw.SkewedElements(self.MeshSpacing, alpha=self.alpha_4)

            self.alpha_0_error_hybrid = self.alpha_0_element.error_dxdy_hybrid
            self.alpha_0_error_standard = self.alpha_0_element.error_dxdy_standard
            self.alpha_1_error_hybrid = self.alpha_1_element.error_dxdy_hybrid
            self.alpha_1_error_standard = self.alpha_1_element.error_dxdy_standard
            self.alpha_2_error_hybrid = self.alpha_2_element.error_dxdy_hybrid
            self.alpha_2_error_standard = self.alpha_2_element.error_dxdy_standard
            self.alpha_3_error_hybrid = self.alpha_3_element.error_dxdy_hybrid
            self.alpha_3_error_standard = self.alpha_3_element.error_dxdy_standard
            self.alpha_4_error_hybrid = self.alpha_4_element.error_dxdy_hybrid
            self.alpha_4_error_standard = self.alpha_4_element.error_dxdy_standard

        @property
        def node_matrix(self):
            """This function creates a matrix of the vertices and the corresponding coordinates depending on the
            provided delta x """

            # Defining parameters and mesh size ------------------------------------------------------------------------

            length = self.mesh_size
            mesh_spacing = self.MeshSpacing
            mesh = np.zeros(shape=(length, length))
            coordinates = np.zeros(shape=(length, length, 2))

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            # Looping through mesh and storing node and coordinate values ----------------------------------------------
            for row in range(length):
                for column in range(length):
                    mesh[row, column] = length * (length - row - 1) + column
                    coordinates[row, column] = np.array([column * mesh_spacing, row * mesh_spacing])

            # ----------------------------------------------------------------------------------------------------------

            return mesh, np.flipud(coordinates)

        @property
        def column_connect(self):
            """This function creates a dictionary of the column vectors according
            to each column labelled from left to right"""

            # Initialing  parameters -----------------------------------------------------------------------------------

            length = self.mesh_size
            node_vals = self.node_matrix[0]
            column_node_dict = {}  # note that the column vector is a dictionary it can be changed to a normal matrix
            # class instances with vector calculations?
            # also replace with pandas dataframe - this is difficult to do for some reason idk why?

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            # Storing vertical nodes in vectors in a dictionary --------------------------------------------------------
            for column in range(length):
                column_node_dict[column] = node_vals[:, column]

            # ----------------------------------------------------------------------------------------------------------

            return column_node_dict

        @property
        def node_connect(self):
            """This function calculates and stores the surrounding vertex numbers in a dictionary"""

            # Initialing  parameters -----------------------------------------------------------------------------------

            internal_vertex_matrix = self.internal_vertex_matrix
            element_center_coordinates = self.element_center_coordinates[0]
            length = self.mesh_size
            num_cells_row = self.cell_size
            node_connect = {}

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            # a mathematical relationship was found to relate the element center coordinates ---------------------------
            # to the node center coordinates ---------------------------------------------------------------------------
            # the coordinates in the cardinal directions SW, SE, NE, NW corresponds to the current row and the column --
            # the following code enforces the relationship -------------------------------------------------------------

            for count, row in enumerate(internal_vertex_matrix):
                for col in row:
                    sw = int(col) - (num_cells_row * 2) + count
                    se = int(col) - (num_cells_row * 2) + count + 1
                    ne = int(col) - num_cells_row + 1 + count
                    nw = int(col) - num_cells_row + count
                    node_connect[col] = np.array([element_center_coordinates[sw],
                                                  element_center_coordinates[se],
                                                  element_center_coordinates[ne],
                                                  element_center_coordinates[nw]])

            # ----------------------------------------------------------------------------------------------------------

            return node_connect

        @property
        def edge_connect(self):
            """This functions returns a dictionary of corresponding adjacent edge centers"""

            # Initialing  parameters -----------------------------------------------------------------------------------

            # We first need the internal vertex matrix that acts as a reference to the adjacent vertices ---------------
            internal_vertex_matrix = self.internal_vertex_matrix
            vertex_coordinates = self.node_matrix[1]
            vertex_connect = {}

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            for row_count, row in enumerate(internal_vertex_matrix):
                for column_count, column in enumerate(row):
                    dict_key = internal_vertex_matrix[row_count][column_count]

                    # the dict key is the vertex number counting from left to right of the node matrix -----------------

                    current_node_coordinate = vertex_coordinates[row_count + 1][column_count + 1]

                    n = (current_node_coordinate + vertex_coordinates[row_count][column_count + 1]) / 2
                    e = (current_node_coordinate + vertex_coordinates[row_count + 1][column_count + 2]) / 2
                    s = (current_node_coordinate + vertex_coordinates[row_count + 2][column_count + 1]) / 2
                    w = (current_node_coordinate + vertex_coordinates[row_count + 1][column_count]) / 2

                    # [s e n w] are the cardinal directions; they refer to the coordinates of the edge centers in the
                    # cardinal direction to the vertex -----------------------------------------------------------------

                    vertex_connect[dict_key] = np.array([s, e, n, w])

                    # DEBUG --------------------------------------------------------------------------------------------

                    # print(dict_key, '\t', current_node_coordinate)

                    # --------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------

            return vertex_connect

        @property
        def element_center_coordinates(self):
            """This function returns the coordinates of the element centers"""
            # This function is required to relate the vertices to element centers

            # Initialing  parameters -----------------------------------------------------------------------------------
            node_matrix = self.node_matrix[0]
            cell_length = self.cell_size
            delta_x = self.MeshSpacing
            x_coords_relative = np.zeros(cell_length * cell_length)
            y_coords_relative = np.zeros(cell_length * cell_length)
            x_coords_absolute = np.zeros(cell_length * cell_length)
            y_coords_absolute = np.zeros(cell_length * cell_length)
            counter = 0

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            for row in range(cell_length):
                for col in range(cell_length):
                    x_coords_relative[counter] = 0.5 + (1 * col)
                    y_coords_relative[counter] = 0.5 + (1 * row)
                    # the relative coordinates are used for plotting but I can probably leave it out
                    x_coords_absolute[counter] = (delta_x / 2) + (delta_x * col)
                    y_coords_absolute[counter] = (delta_x / 2) + (delta_x * row)
                    counter += 1

            coords = np.vstack((x_coords_absolute, y_coords_absolute)).T

            # ----------------------------------------------------------------------------------------------------------

            return [coords, x_coords_relative, y_coords_relative]

        def draw_graph(self):
            """This function uses networkx to draw a 2d grid of the mesh. Purely a
            visualisation package, doing computations with this may be tough"""

            # Initialing  parameters -----------------------------------------------------------------------------------

            length = self.mesh_size
            cell_length = self.cell_size
            x_coords = self.element_center_coordinates[1]
            y_coords = self.element_center_coordinates[2]

            # ----------------------------------------------------------------------------------------------------------

            # Defining grid --------------------------------------------------------------------------------------------
            grid_nodes = ntx.grid_2d_graph(length, length)
            # grid_cells = ntx.grid_2d_graph(cell_length, cell_length)

            position_nodes = dict((n, n) for n in grid_nodes.nodes())
            # position_cells = dict((m, m) for m in grid_cells.nodes())

            labels_nodes = dict(((column, row), (row * length + column)) for row, column in grid_nodes.nodes())
            # labels_cells = dict(((col, row_1), (row_1 * cell_length + col)) for row_1, col in grid_cells.nodes())
            # ----------------------------------------------------------------------------------------------------------

            # Drawing grid ---------------------------------------------------------------------------------------------
            ntx.draw(grid_nodes, pos=position_nodes, node_color='#D9C5B2', labels=labels_nodes, font_color='#1C1C1C')
            # ntx.draw(grid_cells, pos=position_cells, node_color='#000000', labels=labels_cells, font_color='#FFFFFF')
            # ntx.draw(grid, pos=position, node_color='#FFE8C2', labels=labels, font_color='#6A8D73')

            plt.scatter(x_coords, y_coords)

            return

        @property
        def internal_vertex_matrix(self):
            """This function returns the internal node numbers"""

            # Initialing  parameters -----------------------------------------------------------------------------------

            length = self.mesh_size
            elements = self.node_matrix[0]

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            element_matrix = elements[1:(length - 1), 1:(length - 1)]

            # ----------------------------------------------------------------------------------------------------------

            return element_matrix

        @property
        def internal_vertex_coordinates(self):
            """This function returns the internal node coordinates"""

            # Initialing  parameters -----------------------------------------------------------------------------------

            length = self.mesh_size
            vertices = self.internal_vertex_matrix
            coords = self.node_matrix[1]
            array_element_coordinates = coords[1:(length - 1), 1:(length - 1)]
            element_coordinates = {}

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            for row, vertex_row in enumerate(vertices):

                for column, vertex in enumerate(vertex_row):
                    element_coordinates[vertex] = array_element_coordinates[row][column]

            # ----------------------------------------------------------------------------------------------------------

            return element_coordinates

        @property
        def dphi_dxdy_hybrid(self):
            """This function calculates the change in phi with respect to the x and y coordinate for the element
            center to element center scheme"""

            # Initialising parameters ----------------------------------------------------------------------------------

            node_connect = self.node_connect

            # node connect is a dictionary that contains the vertex as a key and the
            # coordinates of the surrounding element centers as the values

            edge_connect = self.edge_connect

            internal_vertex_size = int(np.power((self.cell_size - 1), 2))
            z_vector = self.z_vector
            volume = self.volume_calc
            dphidx = np.zeros(shape=internal_vertex_size)
            dphidy = np.zeros(shape=internal_vertex_size)
            Cfx = np.zeros(shape=internal_vertex_size)
            Cfy = np.zeros(shape=internal_vertex_size)

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

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

                    # the assumption is that phi acts through the center of the edge

                dphidx[node_number] = Cfx[node_number] / volume[node]
                dphidy[node_number] = Cfy[node_number] / volume[node]

            # ----------------------------------------------------------------------------------------------------------

            return dphidx, dphidy

        @property
        def dphi_dxdy_standard(self):
            """This function calculates the change in phi with respect to the x and y coordinate for the element
                center to edge center scheme"""

            # Initialising parameters ----------------------------------------------------------------------------------

            node_connect = self.node_connect
            edge_connect = self.edge_connect

            # node connect and edge connect are dictionaries containing the vertex-surrounding element centers and
            # edge centers

            internal_vertex_size = int(np.power((self.cell_size - 1), 2))
            z_vector = self.z_vector
            volume = self.volume_calc
            dphidx = np.zeros(shape=internal_vertex_size)
            dphidy = np.zeros(shape=internal_vertex_size)
            Cfx = np.zeros(shape=internal_vertex_size)
            Cfy = np.zeros(shape=internal_vertex_size)

            index_edge_center = np.array([0, 0, 1, 1, 2, 2, 3, 3])
            index_element_center = np.array([0, 1, 1, 2, 2, 3, 3, 0])

            # The indices of the coordinates are required to reference the array and calculate the l vector
            # the indices have been hard-coded, this will be much quicker than doing a conditional statement
            # to check the indices

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

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
                    Cfy[vertex_number] += (area * unit_normal_y) * np.cos(edge_center_coordinates
                                                                          [index_edge_center[i]][1]) * 1.5

                dphidx[vertex_number] = Cfx[vertex_number] / volume[vertex]
                dphidy[vertex_number] = Cfy[vertex_number] / volume[vertex]

            # ----------------------------------------------------------------------------------------------------------

            return dphidx, dphidy

        @property
        def volume_calc(self):
            """This code currently calculates the volume for a square element"""

            # Initialising parameters ----------------------------------------------------------------------------------

            # NOTE # This needs to be changed to a general volume calc when skewing the mesh
            node_connect = self.node_connect
            volume = {}

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            for entry in node_connect:
                xy_1 = node_connect[entry][1] - node_connect[entry][0]
                xy_2 = node_connect[entry][3] - node_connect[entry][0]
                volume[entry] = np.cross(xy_1, xy_2)

            # ----------------------------------------------------------------------------------------------------------

            return volume

        @property
        def analytical_standard(self):
            """This function will calculate the analytical values of dphidxdy at the vertex coordinate and return a
            list with the corresponding 'flow' values"""

            # Initialising parameters ----------------------------------------------------------------------------------

            internal_vertex_size = int(np.power((self.cell_size - 1), 2))
            dphidx_analytical = np.zeros(shape=internal_vertex_size)
            dphidy_analytical = np.zeros(shape=internal_vertex_size)
            internal_vertex_coordinates = self.internal_vertex_coordinates

            # Body -----------------------------------------------------------------------------------------------------

            for counter, coordinates in enumerate(internal_vertex_coordinates):
                x_coordinate = internal_vertex_coordinates[coordinates][0]
                y_coordinate = internal_vertex_coordinates[coordinates][1]
                dphidx_analytical[counter] = np.cos(x_coordinate)
                dphidy_analytical[counter] = np.sin(y_coordinate) * -1.5

            # ----------------------------------------------------------------------------------------------------------

            return dphidx_analytical, dphidy_analytical

        @property
        def analytical_hybrid(self):
            """This function will calculate the analytical values of dphidxdy at the vertex coordinate and return a
            list with the corresponding 'flow' values"""

            # Initialising parameters ----------------------------------------------------------------------------------

            internal_vertex_size = int(np.power((self.cell_size - 1), 2))
            dphidx_analytical = np.zeros(shape=internal_vertex_size)
            dphidy_analytical = np.zeros(shape=internal_vertex_size)
            internal_vertex_coordinates = self.internal_vertex_coordinates

            # Body -----------------------------------------------------------------------------------------------------

            for counter, coordinates in enumerate(internal_vertex_coordinates):
                x_coordinate = internal_vertex_coordinates[coordinates][0]
                y_coordinate = internal_vertex_coordinates[coordinates][1]
                dphidx_analytical[counter] = np.cos(x_coordinate)
                dphidy_analytical[counter] = np.sin(y_coordinate) * -1.5

            # ----------------------------------------------------------------------------------------------------------

            return dphidx_analytical, dphidy_analytical

        @property
        def error_dxdy_hybrid(self):
            """This function returns the average error metric for a mesh with the element center to element center
            scheme """

            # Initialising parameters ----------------------------------------------------------------------------------

            dphidxdy = self.dphi_dxdy_hybrid
            analytical_dphi = self.analytical_hybrid
            volume = self.volume_calc
            volume_array = np.array(list(volume.values()))

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            # very straight forward, simply applying the error metric equation :
            # sum(abs(phi_i - phi_i) * V_i) / sum(V_i)

            phi_difference_x = dphidxdy[0] - analytical_dphi[0]
            phi_difference_y = dphidxdy[1] - analytical_dphi[1]
            phi_difference_x_absolute = np.abs(phi_difference_x)
            phi_difference_y_absolute = np.abs(phi_difference_y)

            numerator_x = np.dot(phi_difference_x_absolute, volume_array)
            numerator_y = np.dot(phi_difference_y_absolute, volume_array)
            denominator = np.sum(volume_array)

            error_metric_dx = numerator_x / denominator
            error_metric_dy = numerator_y / denominator

            # ----------------------------------------------------------------------------------------------------------

            return error_metric_dx, error_metric_dy

        @property
        def error_dxdy_standard(self):

            # Initialising parameters ----------------------------------------------------------------------------------

            dhpidxdy = self.dphi_dxdy_standard
            analytical_dphi = self.analytical_standard
            volume = self.volume_calc
            volume_array = np.array(list(volume.values()))

            # ----------------------------------------------------------------------------------------------------------

            # Body -----------------------------------------------------------------------------------------------------

            phi_difference_x = dhpidxdy[0] - analytical_dphi[0]
            phi_difference_y = dhpidxdy[1] - analytical_dphi[1]
            phi_difference_x_absolute = np.abs(phi_difference_x)
            phi_difference_y_absolute = np.abs(phi_difference_y)

            numerator_x = np.dot(phi_difference_x_absolute, volume_array)
            numerator_y = np.dot(phi_difference_y_absolute, volume_array)
            denominator = np.sum(volume_array)

            error_metric_dx = numerator_x / denominator
            error_metric_dy = numerator_y / denominator

            # ----------------------------------------------------------------------------------------------------------

            return error_metric_dx, error_metric_dy

    mesh_1 = Elements(1 / 10)
    print('Done computing:', mesh_1.MeshSpacing )
    mesh_2 = Elements(1 / 20)
    print('Done computing:', mesh_2.MeshSpacing)
    mesh_3 = Elements(1 / 30)
    print('Done computing:', mesh_3.MeshSpacing)
    mesh_4 = Elements(1 / 40)
    print('Done computing:', mesh_4.MeshSpacing)
    print('Finished')

    # This shows the current node matrix based on the mesh size that you specify ---------------------------------------

    # print('The node matrix is: \n', mesh_1.node_matrix[0], '\n')
    # print('The node matrix is: \n', mesh_2.node_matrix[0], '\n')
    # print('The node matrix is: \n', mesh_3.node_matrix[0], '\n')

    # ------------------------------------------------------------------------------------------------------------------

    # This is general code for checking the entire procedure, haven't updated it in a bit so could be wrong ------------

    # print("\nThis is for a delta x of ", mesh_1.delta_x, '\n')
    # print("The nodes are:\n", mesh_1.node_matrix[0], '\n\nThe node coordinates are:\n', mesh_1.node_matrix[1],
    #       '\n')
    # print("The internal element matrix: \n", mesh_1.internal_vertex_matrix, '\n')
    # print("The element center coordinates: \n", mesh_1.element_center_coordinates[2:4], '\n')

    # print("These are the volumes at each vertex: \n")
    # for val in mesh_1.volume_calc:
    #     print(val, ':', mesh_1.volume_calc[val])

    # ------------------------------------------------------------------------------------------------------------------

    # These are the dphidx and dphidy values for the element center to element center scheme ---------------------------

    # print('For element center to element center scheme:')
    # print('\ndphidx values are: ', mesh_1.dphidxdy_element_to_element[0], '\n\n', 'dphidy values are: ',
    #       mesh_1.dphidxdy_element_to_element[1], '\n\n', sep='')

    # ------------------------------------------------------------------------------------------------------------------

    # This checks the functionality of the edge connect function -------------------------------------------------------

    # for entry in mesh_1.edge_connect:
    #     print('\n', entry, ':\n', mesh_1.edge_connect[entry])

    # ------------------------------------------------------------------------------------------------------------------

    # This checks the functionality of the dhpidxdy function for element center to edge center -------------------------

    # print('For element center to edge center scheme:')
    # print('\ndphidx values are: ', mesh_1.dphidxdy_element_to_edge[0], '\n\n', 'dphidy values are: ',
    #       mesh_1.dphidxdy_element_to_edge[1], '\n\n', sep='')

    # ------------------------------------------------------------------------------------------------------------------

    # checking the analytical calculation ------------------------------------------------------------------------------

    # print('Analytical results element center to element center:')
    # print('\ndphidx values are: ', mesh_1.analytical_element_to_element[0], '\n\n', 'dphidy values are: ',
    #       mesh_1.analytical_element_to_element[1], '\n\n', sep='')
    #
    # print('Analytical results element center to edge center:')
    # print('\ndphidx values are: ', mesh_1.analytical_element_to_edge[0], '\n\n', 'dphidy values are: ',
    #       mesh_1.analytical_element_to_edge[1], '\n\n', sep='')

    # ------------------------------------------------------------------------------------------------------------------

    # error metric printouts -------------------------------------------------------------------------------------------

    # print('\n\nThe avg ε for element - element & mesh size', mesh_1.MeshSpacing, ' :',
    #       mesh_1.alpha_0_error_hybrid)
    # print('\n')
    # print('The avg ε for element - edge & mesh size', mesh_1.MeshSpacing, ' :', mesh_1.alpha_0_error_standard,
    #       '\n\n')
    # print('The avg ε for element - element & mesh size', mesh_2.MeshSpacing, ' :', mesh_2.alpha_0_error_hybrid)
    # print('\n')
    # print('The avg ε for element - edge & mesh size', mesh_2.MeshSpacing, ' :', mesh_2.alpha_0_error_standard,
    #       '\n\n')
    # print('The avg ε for element - element & mesh size', mesh_3.MeshSpacing, ' :', mesh_3.alpha_0_error_hybrid)
    # print('\n')
    # print('The avg ε for element - edge & mesh size', mesh_3.MeshSpacing, ' :', mesh_3.alpha_0_error_standard,
    #       '\n\n')
    # print('The avg ε for element - element & mesh size', mesh_4.MeshSpacing, ' :', mesh_4.alpha_0_error_hybrid)
    # print('\n')
    # print('The avg ε for element - edge & mesh size', mesh_4.MeshSpacing, ' :', mesh_4.alpha_0_error_standard,
    #       '\n\n')

    # ------------------------------------------------------------------------------------------------------------------

    # subplot setup ----------------------------------------------------------------------------------------------------

    # fig, (plt, plt) = plt.subplots(1, 2)
    # fig.suptitle('Second order convergence test Φ = sin(x) + 1.5∙cos(y)')

    # ------------------------------------------------------------------------------------------------------------------

    # log error graphs -------------------------------------------------------------------------------------------------

    log_errors_hybrid_alpha_0 = np.array([[np.log(mesh_1.alpha_0_error_hybrid[0]),
                                           np.log(mesh_2.alpha_0_error_hybrid[0]),
                                           np.log(mesh_3.alpha_0_error_hybrid[0]),
                                           np.log(mesh_4.alpha_0_error_hybrid[0])],
                                          [np.log(mesh_1.alpha_0_error_hybrid[1]),
                                           np.log(mesh_2.alpha_0_error_hybrid[1]),
                                           np.log(mesh_3.alpha_0_error_hybrid[1]),
                                           np.log(mesh_4.alpha_0_error_hybrid[1])]
                                          ])

    log_errors_standard_alpha_0 = np.array([[np.log(mesh_1.alpha_0_error_standard[0]),
                                             np.log(mesh_2.alpha_0_error_standard[0]),
                                             np.log(mesh_3.alpha_0_error_standard[0]),
                                             np.log(mesh_4.alpha_0_error_standard[0])],
                                            [np.log(mesh_1.alpha_0_error_standard[1]),
                                             np.log(mesh_2.alpha_0_error_standard[1]),
                                             np.log(mesh_3.alpha_0_error_standard[1]),
                                             np.log(mesh_4.alpha_0_error_standard[1])]
                                            ])

    log_errors_hybrid_alpha_1 = np.array([[np.log(mesh_1.alpha_1_error_hybrid[0]),
                                           np.log(mesh_2.alpha_1_error_hybrid[0]),
                                           np.log(mesh_3.alpha_1_error_hybrid[0]),
                                           np.log(mesh_4.alpha_1_error_hybrid[0])],
                                          [np.log(mesh_1.alpha_1_error_hybrid[1]),
                                           np.log(mesh_2.alpha_1_error_hybrid[1]),
                                           np.log(mesh_3.alpha_1_error_hybrid[1]),
                                           np.log(mesh_4.alpha_1_error_hybrid[1])]
                                          ])

    log_errors_standard_alpha_1 = np.array([[np.log(mesh_1.alpha_1_error_standard[0]),
                                             np.log(mesh_2.alpha_1_error_standard[0]),
                                             np.log(mesh_3.alpha_1_error_standard[0]),
                                             np.log(mesh_4.alpha_1_error_standard[0])],
                                            [np.log(mesh_1.alpha_1_error_standard[1]),
                                             np.log(mesh_2.alpha_1_error_standard[1]),
                                             np.log(mesh_3.alpha_1_error_standard[1]),
                                             np.log(mesh_4.alpha_1_error_standard[1])]
                                            ])

    log_errors_hybrid_alpha_2 = np.array([[np.log(mesh_1.alpha_2_error_hybrid[0]),
                                           np.log(mesh_2.alpha_2_error_hybrid[0]),
                                           np.log(mesh_3.alpha_2_error_hybrid[0]),
                                           np.log(mesh_4.alpha_2_error_hybrid)[0]],
                                          [np.log(mesh_1.alpha_2_error_hybrid[1]),
                                           np.log(mesh_2.alpha_2_error_hybrid[1]),
                                           np.log(mesh_3.alpha_2_error_hybrid[1]),
                                           np.log(mesh_4.alpha_2_error_hybrid)[1]]
                                          ])

    log_errors_standard_alpha_2 = np.array([[np.log(mesh_1.alpha_2_error_standard[0]),
                                             np.log(mesh_2.alpha_2_error_standard[0]),
                                             np.log(mesh_3.alpha_2_error_standard[0]),
                                             np.log(mesh_4.alpha_2_error_standard[0])],
                                            [np.log(mesh_1.alpha_2_error_standard[1]),
                                             np.log(mesh_2.alpha_2_error_standard[1]),
                                             np.log(mesh_3.alpha_2_error_standard[1]),
                                             np.log(mesh_4.alpha_2_error_standard[1])]
                                            ])

    log_errors_hybrid_alpha_3 = np.array([[np.log(mesh_1.alpha_3_error_hybrid[0]),
                                           np.log(mesh_2.alpha_3_error_hybrid[0]),
                                           np.log(mesh_3.alpha_3_error_hybrid[0]),
                                           np.log(mesh_4.alpha_3_error_hybrid)[0]],
                                          [np.log(mesh_1.alpha_3_error_hybrid)[1],
                                           np.log(mesh_2.alpha_3_error_hybrid)[1],
                                           np.log(mesh_3.alpha_3_error_hybrid)[1],
                                           np.log(mesh_4.alpha_3_error_hybrid)[1]]
                                          ])

    log_errors_standard_alpha_3 = np.array([[np.log(mesh_1.alpha_3_error_standard[0]),
                                             np.log(mesh_2.alpha_3_error_standard[0]),
                                             np.log(mesh_3.alpha_3_error_standard[0]),
                                             np.log(mesh_4.alpha_3_error_standard[0])],
                                            [np.log(mesh_1.alpha_3_error_standard[1]),
                                             np.log(mesh_2.alpha_3_error_standard[1]),
                                             np.log(mesh_3.alpha_3_error_standard[1]),
                                             np.log(mesh_4.alpha_3_error_standard[1])]
                                            ])

    log_errors_hybrid_alpha_4 = np.array([[np.log(mesh_1.alpha_4_error_hybrid[0]),
                                           np.log(mesh_2.alpha_4_error_hybrid[0]),
                                           np.log(mesh_3.alpha_4_error_hybrid[0]),
                                           np.log(mesh_4.alpha_4_error_hybrid[0])],
                                          [np.log(mesh_1.alpha_4_error_hybrid[1]),
                                           np.log(mesh_2.alpha_4_error_hybrid[1]),
                                           np.log(mesh_3.alpha_4_error_hybrid[1]),
                                           np.log(mesh_4.alpha_4_error_hybrid[1])]
                                          ])

    log_errors_standard_alpha_4 = np.array([[np.log(mesh_1.alpha_4_error_standard[0]),
                                             np.log(mesh_2.alpha_4_error_standard[0]),
                                             np.log(mesh_3.alpha_4_error_standard[0]),
                                             np.log(mesh_4.alpha_4_error_standard[0])],
                                            [np.log(mesh_1.alpha_4_error_standard[1]),
                                             np.log(mesh_2.alpha_4_error_standard[1]),
                                             np.log(mesh_3.alpha_4_error_standard[1]),
                                             np.log(mesh_4.alpha_4_error_standard[1])]
                                            ])

    delta_x_log = np.array([np.log(mesh_1.MeshSpacing),
                            np.log(mesh_2.MeshSpacing),
                            np.log(mesh_3.MeshSpacing),
                            np.log(mesh_4.MeshSpacing)])

    second_order_line = np.array([[-4.8, -2.0], [-14.2, -8.6]])
    first_order_line = np.array([[-4.8, -2.0], [-10, -7.2]])

    fig = plt.figure(figsize=plt.figaspect(0.5))

    axs = fig.add_subplot(2, 2, 1)

    axs.plot(delta_x_log, log_errors_hybrid_alpha_0[0], color='#2E86AB',
             label='Hybrid scheme: α = 0.0', marker='D')  # Blue NCS
    axs.plot(delta_x_log, log_errors_standard_alpha_0[0], color='blue',
             label='Standard scheme: α = 0.0', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_4[0], color='#74591B',
             label='Hybrid scheme: α = 0.05', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_4[0], color='#D5AC4E',
             label='Standard scheme: α = 0.05', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_1[0], color='#D64933',
             label='Hybrid scheme: α = 0.1', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_1[0], color='#96031A',
             label='Standard scheme: α = 0.1', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_2[0], color='#26262C',
             label='Hybrid scheme: α = 0.2', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_2[0], color='#050505',
             label='Standard scheme: α = 0.2', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_3[0], color='#A1E8AF',
             label='Hybrid scheme: α = 0.46', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_3[0], color='#3BC14A',
             label='Standard scheme: α = 0.46', linestyle=':')

    axs.plot(second_order_line[0], second_order_line[1], color='red', linestyle=(0, (5, 1)), label='Second order '
                                                                                                   'convergence')
    axs.plot(first_order_line[0], first_order_line[1], color='red', linestyle=(0, (3, 5, 1, 5, 1, 5)),
             label='First order convergence')
    axs.set_xlim(-4.8, -2)
    axs.set_xlabel('log(∆x)')
    axs.set_ylabel('log(ε)')
    axs.set_title('dhpidx Log-log graph')
    # axs.legend()
    axs.grid(color='black', linestyle='--', linewidth=0.5, alpha=0.3)

    axs = fig.add_subplot(2, 2, 2)

    axs.plot(delta_x_log, log_errors_hybrid_alpha_0[1], color='#2E86AB',
             label='Hybrid scheme: α = 0.0', marker='D')  # Blue NCS
    axs.plot(delta_x_log, log_errors_standard_alpha_0[1], color='blue',
             label='Standard scheme: α = 0.0', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_4[1], color='#74591B',
             label='Hybrid scheme: α = 0.05', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_4[1], color='#D5AC4E',
             label='Standard scheme: α = 0.05', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_1[1], color='#D64933',
             label='Hybrid scheme: α = 0.1', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_1[1], color='#96031A',
             label='Standard scheme: α = 0.1', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_2[1], color='#26262C',
             label='Hybrid scheme: α = 0.2', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_2[1], color='#050505',
             label='Standard scheme: α = 0.2', linestyle=':')
    axs.plot(delta_x_log, log_errors_hybrid_alpha_3[1], color='#A1E8AF',
             label='Hybrid scheme: α = 0.46', marker='D')
    axs.plot(delta_x_log, log_errors_standard_alpha_3[1], color='#3BC14A',
             label='Standard scheme: α = 0.46', linestyle=':')

    # plt.plot(second_order_line[0], second_order_line[1], color='black', linestyle=(0, (5, 1)))
    # plt.plot(first_order_line[0], first_order_line[1], color='black', linestyle=(0, (5, 1)))

    axs.plot(second_order_line[0], second_order_line[1], color='red', linestyle=(0, (5, 1)), label='Second order '
                                                                                                   'convergence')
    axs.plot(first_order_line[0], first_order_line[1], color='red', linestyle=(0, (3, 5, 1, 5, 1, 5)),
             label='First order convergence')
    axs.set_xlim(-4.8, -2)
    axs.set_xlabel('log(∆x)')
    axs.set_ylabel('log(ε)')
    axs.set_title('dhpidy Log-log graph')
    # axs.legend()
    handles, labels = axs.get_legend_handles_labels()
    axs.grid(color='black', linestyle='--', linewidth=0.5, alpha=0.3)

    # ------------------------------------------------------------------------------------------------------------------

    # standard error graphs --------------------------------------------------------------------------------------------

    errors_hybrid_alpha_0 = np.array([[mesh_1.alpha_0_error_hybrid[0],
                                       mesh_2.alpha_0_error_hybrid[0],
                                       mesh_3.alpha_0_error_hybrid[0],
                                       mesh_4.alpha_0_error_hybrid[0]],
                                      [mesh_1.alpha_0_error_hybrid[1],
                                       mesh_2.alpha_0_error_hybrid[1],
                                       mesh_3.alpha_0_error_hybrid[1],
                                       mesh_4.alpha_0_error_hybrid[1]]
                                      ])

    errors_standard_alpha_0 = np.array([[mesh_1.alpha_0_error_standard[0],
                                         mesh_2.alpha_0_error_standard[0],
                                         mesh_3.alpha_0_error_standard[0],
                                         mesh_4.alpha_0_error_standard[0]],
                                        [mesh_1.alpha_0_error_standard[1],
                                         mesh_2.alpha_0_error_standard[1],
                                         mesh_3.alpha_0_error_standard[1],
                                         mesh_4.alpha_0_error_standard[1]]
                                        ])

    errors_hybrid_alpha_1 = np.array([[mesh_1.alpha_1_error_hybrid[0],
                                       mesh_2.alpha_1_error_hybrid[0],
                                       mesh_3.alpha_1_error_hybrid[0],
                                       mesh_4.alpha_1_error_hybrid[0]],
                                      [mesh_1.alpha_1_error_hybrid[1],
                                       mesh_2.alpha_1_error_hybrid[1],
                                       mesh_3.alpha_1_error_hybrid[1],
                                       mesh_4.alpha_1_error_hybrid[1]]
                                      ])

    errors_standard_alpha_1 = np.array([[mesh_1.alpha_1_error_standard[0],
                                         mesh_2.alpha_1_error_standard[0],
                                         mesh_3.alpha_1_error_standard[0],
                                         mesh_4.alpha_1_error_standard[0]],
                                        [mesh_1.alpha_1_error_standard[1],
                                         mesh_2.alpha_1_error_standard[1],
                                         mesh_3.alpha_1_error_standard[1],
                                         mesh_4.alpha_1_error_standard[1]]
                                        ])

    errors_hybrid_alpha_2 = np.array([[mesh_1.alpha_2_error_hybrid[0],
                                       mesh_2.alpha_2_error_hybrid[0],
                                       mesh_3.alpha_2_error_hybrid[0],
                                       mesh_4.alpha_2_error_hybrid[0]],
                                      [mesh_1.alpha_2_error_hybrid[1],
                                       mesh_2.alpha_2_error_hybrid[1],
                                       mesh_3.alpha_2_error_hybrid[1],
                                       mesh_4.alpha_2_error_hybrid[1]]
                                      ])

    errors_standard_alpha_2 = np.array([[mesh_1.alpha_2_error_standard[0],
                                         mesh_2.alpha_2_error_standard[0],
                                         mesh_3.alpha_2_error_standard[0],
                                         mesh_4.alpha_2_error_standard[0]],
                                        [mesh_1.alpha_2_error_standard[1],
                                         mesh_2.alpha_2_error_standard[1],
                                         mesh_3.alpha_2_error_standard[1],
                                         mesh_4.alpha_2_error_standard[1]]
                                        ])

    errors_hybrid_alpha_3 = np.array([[mesh_1.alpha_3_error_hybrid[0],
                                       mesh_2.alpha_3_error_hybrid[0],
                                       mesh_3.alpha_3_error_hybrid[0],
                                       mesh_4.alpha_3_error_hybrid[0]],
                                      [mesh_1.alpha_3_error_hybrid[1],
                                       mesh_2.alpha_3_error_hybrid[1],
                                       mesh_3.alpha_3_error_hybrid[1],
                                       mesh_4.alpha_3_error_hybrid[1]]
                                      ])

    errors_standard_alpha_3 = np.array([[mesh_1.alpha_3_error_standard[0],
                                         mesh_2.alpha_3_error_standard[0],
                                         mesh_3.alpha_3_error_standard[0],
                                         mesh_4.alpha_3_error_standard[0]],
                                        [mesh_1.alpha_3_error_standard[1],
                                         mesh_2.alpha_3_error_standard[1],
                                         mesh_3.alpha_3_error_standard[1],
                                         mesh_4.alpha_3_error_standard[1]]
                                        ])

    errors_hybrid_alpha_4 = np.array([[mesh_1.alpha_4_error_hybrid[0],
                                       mesh_2.alpha_4_error_hybrid[0],
                                       mesh_3.alpha_4_error_hybrid[0],
                                       mesh_4.alpha_4_error_hybrid[0]],
                                      [mesh_1.alpha_4_error_hybrid[1],
                                       mesh_2.alpha_4_error_hybrid[1],
                                       mesh_3.alpha_4_error_hybrid[1],
                                       mesh_4.alpha_4_error_hybrid[1]]
                                      ])

    errors_standard_alpha_4 = np.array([[mesh_1.alpha_4_error_standard[0],
                                         mesh_2.alpha_4_error_standard[0],
                                         mesh_3.alpha_4_error_standard[0],
                                         mesh_4.alpha_4_error_standard[0]],
                                        [mesh_1.alpha_4_error_standard[1],
                                         mesh_2.alpha_4_error_standard[1],
                                         mesh_3.alpha_4_error_standard[1],
                                         mesh_4.alpha_4_error_standard[1]]
                                        ])

    delta_x_squared = np.array([np.power(mesh_1.MeshSpacing, 2),
                                np.power(mesh_2.MeshSpacing, 2),
                                np.power(mesh_3.MeshSpacing, 2),
                                np.power(mesh_4.MeshSpacing, 2)])

    axs = fig.add_subplot(2, 2, 3)

    axs.plot(delta_x_squared, errors_hybrid_alpha_0[0], color='#2E86AB',
             label='Hybrid scheme: α = 0.0', marker='D')  # Cinnabar - D64933
    axs.plot(delta_x_squared, errors_standard_alpha_0[0], color='blue',
             label='Standard scheme: α = 0.0', linestyle=':')  # Carmine - 96031A
    axs.plot(delta_x_squared, errors_hybrid_alpha_4[0], color='#74591B',
             label='Hybrid scheme: α = 0.02', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_4[0], color='#D5AC4E',
             label='Standard scheme: α = 0.02', linestyle=':')
    axs.plot(delta_x_squared, errors_hybrid_alpha_1[0], color='#D64933',
             label='Hybrid scheme: α = 0.1', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_1[0], color='#96031A',
             label='Standard scheme: α = 0.1', linestyle=':')
    axs.plot(delta_x_squared, errors_hybrid_alpha_2[0], color='#26262C',
             label='Hybrid scheme: α = 0.2', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_2[0], color='#050505',
             label='Standard scheme: α = 0.2', linestyle=':')
    axs.plot(delta_x_squared, errors_hybrid_alpha_3[0], color='#A1E8AF',
             label='Hybrid scheme: α = 0.46', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_3[0], color='#3BC14A',
             label='Standard scheme: α = 0.46', linestyle=':')

    axs.set_xlabel('∆x²')
    axs.set_ylabel('ε')
    axs.set_title('dphidx Standard error graph')
    # axs.legend()
    axs.grid(color='black', linestyle='--', linewidth=0.5, alpha=0.3)

    axs = fig.add_subplot(2, 2, 4)

    axs.plot(delta_x_squared, errors_hybrid_alpha_0[1], color='#2E86AB',
             label='Hybrid scheme: α = 0.0', marker='D')  # Cinnabar - D64933
    axs.plot(delta_x_squared, errors_standard_alpha_0[1], color='blue',
             label='Standard scheme: α = 0.0', linestyle=':')  # Carmine - 96031A
    axs.plot(delta_x_squared, errors_hybrid_alpha_4[1], color='#3E2F5B',
             label='Hybrid scheme: α = 0.02', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_4[1], color='#814B77',
             label='Standard scheme: α = 0.02', linestyle=':')
    axs.plot(delta_x_squared, errors_hybrid_alpha_1[1], color='#D64933',
             label='Hybrid scheme: α = 0.1', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_1[1], color='#96031A',
             label='Standard scheme: α = 0.1', linestyle=':')
    axs.plot(delta_x_squared, errors_hybrid_alpha_2[1], color='#26262C',
             label='Hybrid scheme: α = 0.2', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_2[1], color='#050505',
             label='Standard scheme: α = 0.2', linestyle=':')
    axs.plot(delta_x_squared, errors_hybrid_alpha_3[1], color='#A1E8AF',
             label='Hybrid scheme: α = 0.46', marker='D')
    axs.plot(delta_x_squared, errors_standard_alpha_3[1], color='#3BC14A',
             label='Standard scheme: α = 0.46', linestyle=':')

    axs.set_xlabel('∆x²')
    axs.set_ylabel('$\epsilon$')
    axs.set_title('dphidy Standard error graph')
    # axs.legend
    axs.grid(color='black', linestyle='--', linewidth=0.5, alpha=0.3)
    fig.legend(handles, labels, loc="lower center", ncol=(int(len(labels)/4)))
    plt.show()

    with open('errors.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerow(['h', 'h', 'h', 'h'])
        writer.writerow(errors_hybrid_alpha_0[0])
        writer.writerow(errors_hybrid_alpha_0[1])
        writer.writerow(errors_hybrid_alpha_4[0])
        writer.writerow(errors_hybrid_alpha_4[1])
        writer.writerow(errors_hybrid_alpha_1[0])
        writer.writerow(errors_hybrid_alpha_1[1])
        writer.writerow(errors_hybrid_alpha_2[0])
        writer.writerow(errors_hybrid_alpha_2[1])
        writer.writerow(errors_hybrid_alpha_3[0])
        writer.writerow(errors_hybrid_alpha_3[1])

        writer.writerow(['s', 's', 's', 's'])
        writer.writerow(errors_standard_alpha_0[0])
        writer.writerow(errors_standard_alpha_0[1])
        writer.writerow(errors_standard_alpha_4[0])
        writer.writerow(errors_standard_alpha_4[1])
        writer.writerow(errors_standard_alpha_1[0])
        writer.writerow(errors_standard_alpha_1[1])
        writer.writerow(errors_standard_alpha_2[0])
        writer.writerow(errors_standard_alpha_2[1])
        writer.writerow(errors_standard_alpha_3[0])
        writer.writerow(errors_standard_alpha_3[1])

    # ------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
