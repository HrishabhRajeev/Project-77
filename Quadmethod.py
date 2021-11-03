import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()


def main():
    class Quadrilateral:

        # This object uses a cellular layout to compute nodes and grids based on adjacent edges and vertices -----------
        def __init__(self, delta_x):
            self.delta_x = delta_x
            self.mesh_size = int(1 / delta_x)
            self.x_scatter = np.zeros(shape=self.mesh_size)
            self.y_scatter = np.zeros(shape=self.mesh_size)
            self.x1y1 = np.array([0, 0])
            self.x2y2 = np.array([0, 0])
            self.x3y3 = np.array([0, 0])
            self.x4y4 = np.array([0, 0])

        def draw_graph(self):

            # This is fucking useless theres already a matplotlib thing that makes this already ------------------------
            # However I might be able to more easily adjust the delta_x value with this code ---------------------------
            # Another benefit is that you are essentially doing node_connect in each iteration of the loop -------------
            # You can also easily store the vertical node vals ---------------------------------------------------------
            # This code is also doing things sequentially one row at a time, which means its slow ----------------------
            delta_x = self.delta_x
            init_x1y1 = self.x1y1
            init_x2y2 = self.x2y2
            init_x3y3 = self.x3y3
            init_x4y4 = self.x4y4

            for row in range(self.mesh_size):

                for column in range(self.mesh_size):
                    x1y1 = np.add(init_x1y1, np.array([delta_x * column, delta_x * row]))
                    x2y2 = np.add(init_x2y2, np.array([delta_x + (delta_x * column), delta_x * row]))
                    x3y3 = np.add(init_x3y3, np.array([delta_x + (delta_x * column), delta_x + (delta_x * row)]))
                    x4y4 = np.add(init_x4y4, np.array([delta_x * column, delta_x + (delta_x * row)]))

                    matrix = np.array([x1y1, x2y2, x3y3, x4y4, x1y1])

                    xvals = matrix[:, 0]
                    yvals = matrix[:, 1]

                    # np.append(x_scatter, xvals)
                    # np.append(y_scatter, yvals)

                    plt.xlim(-0.1, 1.1)
                    plt.ylim(-0.1, 1.1)
                    plt.scatter(xvals[0:4], yvals[0:4], color="#000000")
                    plt.plot(xvals, yvals, color='#000000')
                    # plt.show()

            return

        # def init(self):

    mesh = Quadrilateral(1 / 10)
    mesh.draw_graph()
    plt.show()


if __name__ == '__main__':
    main()
