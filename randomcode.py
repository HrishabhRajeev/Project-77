import numpy as np
from numpy import pi as pi
from scipy.spatial import Delaunay
import matplotlib.pylab as plt
from scipy.optimize import fmin


def ktrimesh(p, bars, pflag=0):
    # create the (x,y) data for the plot
    xx1 = p[bars[:, 0], 0]; yy1 = p[bars[:, 0], 1]
    xx2 = p[bars[:, 1], 0]; yy2 = p[bars[:, 1], 1]
    xmin = np.min(p[:, 0])
    xmax = np.max(p[:, 0])
    ymin = np.min(p[:, 1])
    ymax = np.max(p[:, 1])
    xmin = xmin - 0.05 * (xmax - xmin)
    xmax = xmax + 0.05 * (xmax - xmin)
    ymin = ymin - 0.05 * (ymax - ymin)
    ymax = ymax + 0.05 * (ymax - ymin)

    plt.figure()
    for i in range(len(xx1)):
        xp = np.array([xx1[i], xx2[i]])
        yp = np.array([yy1[i], yy2[i]])
        plt.plot(xmin, ymin, '.', xmax, ymax, '.', markersize=0.1)
        plt.plot(xp, yp, 'k')
    plt.axis('equal')
    if pflag == 0:
        stitle = 'Triangular Mesh'
    if pflag == 1:
        stitle = 'Visual Boundary Integrity Check'
    # plt.title('Triangular Mesh')
    plt.title(stitle)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    return 1


def ccw_tri(p, t):
    """
    orients all the triangles counterclockwise
    """
    # vector A from vertex 0 to vertex 1
    # vector B from vertex 0 to vertex 2
    A01x = p[t[:, 1], 0] - p[t[:, 0], 0]
    A01y = p[t[:, 1], 1] - p[t[:, 0], 1]
    B02x = p[t[:, 2], 0] - p[t[:, 0], 0]
    B02y = p[t[:, 2], 1] - p[t[:, 0], 1]
    # if vertex 2 lies to the left of vector A the component z of
    # their vectorial product A^B is positive
    Cz = A01x * B02y - A01y * B02x
    a = t[np.where(Cz < 0)]
    b = t[np.where(Cz >= 0)]
    a[:, [1, 2]] = a[:, [2, 1]]
    t = np.concatenate((a, b))
    return t


def triqual_flag(p, t):
    # a(1,0), b(2,0), c(2,1)
    a = np.sqrt((p[t[:, 1], 0] - p[t[:, 0], 0]) ** 2 + (p[t[:, 1], 1] - p[t[:, 0], 1]) ** 2)
    b = np.sqrt((p[t[:, 2], 0] - p[t[:, 0], 0]) ** 2 + (p[t[:, 2], 1] - p[t[:, 0], 1]) ** 2)
    c = np.sqrt((p[t[:, 2], 0] - p[t[:, 1], 0]) ** 2 + (p[t[:, 2], 1] - p[t[:, 1], 1]) ** 2)
    A = 0.25 * np.sqrt((a + b + c) * (b + c - a) * (a + c - b) * (a + b - c))
    R = 0.25 * (a * b * c) / A
    r = 0.5 * np.sqrt((a + b - c) * (b + c - a) * (a + c - b) / (a + b + c))
    q = 2.0 * (r / R)
    min_edge = np.minimum(np.minimum(a, b), c)
    min_angle_deg = (180.0 / np.pi) * np.arcsin(0.5 * min_edge / R)

    min_q = np.min(q)
    min_ang = np.min(min_angle_deg)
    return min_q, min_ang


def triqual(p, t, fh, qlim=0.2):
    # a(1,0), b(2,0), c(2,1)
    a = np.sqrt((p[t[:, 1], 0] - p[t[:, 0], 0]) ** 2 + (p[t[:, 1], 1] - p[t[:, 0], 1]) ** 2)
    b = np.sqrt((p[t[:, 2], 0] - p[t[:, 0], 0]) ** 2 + (p[t[:, 2], 1] - p[t[:, 0], 1]) ** 2)
    c = np.sqrt((p[t[:, 2], 0] - p[t[:, 1], 0]) ** 2 + (p[t[:, 2], 1] - p[t[:, 1], 1]) ** 2)
    A = 0.25 * np.sqrt((a + b + c) * (b + c - a) * (a + c - b) * (a + b - c))
    R = 0.25 * (a * b * c) / A
    r = 0.5 * np.sqrt((a + b - c) * (b + c - a) * (a + c - b) / (a + b + c))
    q = 2.0 * (r / R)
    pmid = (p[t[:, 0]] + p[t[:, 1]] + p[t[:, 2]]) / 3.0
    hmid = fh(pmid)
    Ah = A / hmid
    Anorm = Ah / np.mean(Ah)
    min_edge = np.minimum(np.minimum(a, b), c)
    min_angle_deg = (180.0 / np.pi) * np.arcsin(0.5 * min_edge / R)

    plt.figure()
    plt.subplot(3, 1, 1)
    plt.hist(q)
    plt.title('Histogram;Triangle Statistics:q-factor,Minimum Angle and Area')
    plt.subplot(3, 1, 2)
    plt.hist(min_angle_deg)
    plt.ylabel('Number of Triangles')
    plt.subplot(3, 1, 3)
    plt.hist(Anorm)
    plt.xlabel('Note: for equilateral triangles q = 1 and angle = 60 deg')
    plt.show()

    indq = np.where(q < qlim)  # indq is a tuple: len(indq) = 1
    if list(indq[0]) != []:
        print('List of triangles with q < %5.3f and the (x,y) location of their nodes' % qlim)
        print('')
        print('q     t[i]      t[nodes]         [x,y][0]       [x,y][1]       [x,y][2]')
        for i in indq[0]:
            print('%.2f  %4d  [%4d,%4d,%4d]     [%+.2f,%+.2f]  [%+.2f,%+.2f]  [%+.2f,%+.2f]' % \
                  (q[i], i, t[i, 0], t[i, 1], t[i, 2], p[t[i, 0], 0], p[t[i, 0], 1], p[t[i, 1], 0], p[t[i, 1], 1],
                   p[t[i, 2], 0], p[t[i, 2], 1]))
        print('')
        # end of detailed data on worst offenders
    return q, min_angle_deg, Anorm


class Circle:
    def __init__(self, xc, yc, r):
        self.xc, self.yc, self.r = xc, yc, r

    def __call__(self, p):
        xc, yc, r = self.xc, self.yc, self.r
        d = np.sqrt((p[:, 0] - xc) ** 2 + (p[:, 1] - yc) ** 2) - r
        return d


class Rectangle:
    def __init__(self, x1, x2, y1, y2):
        self.x1, self.x2, self.y1, self.y2 = x1, x2, y1, y2

    def __call__(self, p):
        x1, x2, y1, y2 = self.x1, self.x2, self.y1, self.y2
        d1 = p[:, 1] - y1  # if p inside d1 > 0
        d2 = y2 - p[:, 1]  # if p inside d2 > 0
        d3 = p[:, 0] - x1  # if p inside d3 > 0
        d4 = x2 - p[:, 0]  # if p inside d4 > 0
        d = -np.minimum(np.minimum(np.minimum(d1, d2), d3), d4)
        return d


class Polygon:
    def __init__(self, verts):
        self.verts = verts

    def __call__(self, p):
        verts = self.verts
        # close the polygon
        cverts = np.zeros((len(verts) + 1, 2))
        cverts[0:-1] = verts
        cverts[-1] = verts[0]
        # initialize
        inside = np.zeros(len(p))
        dist = np.zeros(len(p))
        Cz = np.zeros(len(verts))  # z-components of the vectorial products
        dist_to_edge = np.zeros(len(verts))
        in_ref = np.ones(len(verts))
        # if np.sign(Cz) == in_ref then point is inside
        for j in range(len(p)):
            Cz = (cverts[1:, 0] - cverts[0:-1, 0]) * (p[j, 1] - cverts[0:-1, 1]) - \
                 (cverts[1:, 1] - cverts[0:-1, 1]) * (p[j, 0] - cverts[0:-1, 0])
            dist_to_edge = Cz / np.sqrt( \
                (cverts[1:, 0] - cverts[0:-1, 0]) ** 2 + \
                (cverts[1:, 1] - cverts[0:-1, 1]) ** 2)

            inside[j] = int(np.array_equal(np.sign(Cz), in_ref))
            dist[j] = (1 - 2 * inside[j]) * np.min(np.abs(dist_to_edge))
        return dist


class Union:
    def __init__(self, fd1, fd2):
        self.fd1, self.fd2 = fd1, fd2

    def __call__(self, p):
        fd1, fd2 = self.fd1, self.fd2
        d = np.minimum(fd1(p), fd2(p))
        return d


class Diff:
    def __init__(self, fd1, fd2):
        self.fd1, self.fd2 = fd1, fd2

    def __call__(self, p):
        fd1, fd2 = self.fd1, self.fd2
        d = np.maximum(fd1(p), -fd2(p))
        return d


class Intersect:
    def __init__(self, fd1, fd2):
        self.fd1, self.fd2 = fd1, fd2

    def __call__(self, p):
        fd1, fd2 = self.fd1, self.fd2
        d = np.maximum(fd1(p), fd2(p))
        return d


class Protate:
    def __init__(self, phi):
        self.phi = phi

    def __call__(self, p):
        phi = self.phi
        c = np.cos(phi)
        s = np.sin(phi)
        temp = np.copy(p[:, 0])
        rp = np.copy(p)
        rp[:, 0] = c * p[:, 0] - s * p[:, 1]
        rp[:, 1] = s * temp + c * p[:, 1]
        return rp


class Pshift:
    def __init__(self, x0, y0):
        self.x0, self.y0 = x0, y0

    def __call__(self, p):
        x0, y0 = self.x0, self.y0
        p[:, 0] = p[:, 0] + x0
        p[:, 1] = p[:, 1] + y0
        return p


def Ellipse_dist_to_minimize(t, p, xc, yc, a, b):
    x = xc + a * np.cos(t)  # coord x of the point on the ellipse
    y = yc + b * np.sin(t)  # coord y of the point on the ellipse
    dist = (p[0] - x) ** 2 + (p[1] - y) ** 2
    return dist


class Ellipse:
    def __init__(self, xc, yc, a, b):
        self.xc, self.yc, self.a, self.b = xc, yc, a, b
        self.t, self.verts = self.pick_points_on_shape()

    def pick_points_on_shape(self):
        xc, yc, a, b = self.xc, self.yc, self.a, self.b
        c = np.array([xc, yc])
        t = np.linspace(0, (7.0 / 4.0) * pi, 8)
        verts = np.zeros((8, 2))
        verts[:, 0] = c[0] + a * np.cos(t)
        verts[:, 1] = c[1] + b * np.sin(t)
        return t, verts

    def inside_ellipse(self, p):
        xc, yc, a, b = self.xc, self.yc, self.a, self.b
        c = np.array([xc, yc])
        r, phase = self.rect_to_polar(p - c)
        r_ellipse = self.rellipse(phase)
        in_ref = np.ones(len(p))
        inside = 0.5 + 0.5 * np.sign(r_ellipse - r)
        return inside

    def rect_to_polar(self, p):
        r = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2)
        phase = np.arctan2(p[:, 1], p[:, 0])
        # note: np.arctan2(y,x) order; phase in +/- pi (+/- 180deg)
        return r, phase

    def rellipse(self, phi):
        a, b = self.a, self.b
        r = a * b / np.sqrt((b * np.cos(phi)) ** 2 + (a * np.sin(phi)) ** 2)
        return r

    def find_closest_vertex(self, point):
        t, verts = self.t, self.verts

        dist = np.zeros(len(t))
        for i in range(len(t)):
            dist[i] = (point[0] - verts[i, 0]) ** 2 + (point[1] - verts[i, 1]) ** 2
        ind = np.argmin(dist)
        t0 = t[ind]
        return t0

    def __call__(self, p):
        xc, yc, a, b = self.xc, self.yc, self.a, self.b
        t, verts = self.t, self.verts
        dist = np.zeros(len(p))
        inside = self.inside_ellipse(p)
        for j in range(len(p)):
            t0 = self.find_closest_vertex(p[j])  # initial guess to minimizer
            opt = fmin(Ellipse_dist_to_minimize, t0, \
                       args=(p[j], xc, yc, a, b), full_output=1, disp=0)
            # add full_output=1 so we can retrieve the min dist(squared)
            # (2nd argument of opt array, 1st argument is the optimum t)
            min_dist = np.sqrt(opt[1])
            dist[j] = min_dist * (1 - 2 * inside[j])
        return dist


def distmesh(fd, fh, h0, xmin, ymin, xmax, ymax, pfix, ttol=0.1, dptol=0.001, Iflag=1, qmin=1.0):
    geps = 0.001 * h0; deltat = 0.2; Fscale = 1.2
    deps = h0 * np.sqrt(np.spacing(1))

    random_seed = 17

    h0x = h0; h0y = h0 * np.sqrt(3) / 2  # to obtain equilateral triangles
    Nx = int(np.floor((xmax - xmin) / h0x))
    Ny = int(np.floor((ymax - ymin) / h0y))
    x = np.linspace(xmin, xmax, Nx)
    y = np.linspace(ymin, ymax, Ny)
    # create the grid in the (x,y) plane
    xx, yy = np.meshgrid(x, y)
    xx[1::2] = xx[1::2] + h0x / 2.0  # shifts even rows by h0x/2
    p = np.zeros((np.size(xx), 2))
    p[:, 0] = np.reshape(xx, np.size(xx))
    p[:, 1] = np.reshape(yy, np.size(yy))

    p = np.delete(p, np.where(fd(p) > geps), axis=0)

    np.random.seed(random_seed)
    r0 = 1.0 / fh(p) ** 2
    p = np.concatenate((pfix, p[np.random.rand(len(p)) < r0 / max(r0), :]))

    pold = np.inf
    Num_of_Delaunay_triangulations = 0
    Num_of_Node_movements = 0  # dp = F*dt

    while (1):
        Num_of_Node_movements += 1
        if Iflag == 1 or Iflag == 3:  # Newton flag
            print('Num_of_Node_movements = %3d' % (Num_of_Node_movements))
        if np.max(np.sqrt(np.sum((p - pold) ** 2, axis=1))) > ttol:
            Num_of_Delaunay_triangulations += 1
            if Iflag == 1 or Iflag == 3:  # Delaunay flag
                print('Num_of_Delaunay_triangulations = %3d' % \
                      (Num_of_Delaunay_triangulations))
            pold = p
            tri = Delaunay(p)  # instantiate a class
            t = tri.vertices
            pmid = (p[t[:, 0]] + p[t[:, 1]] + p[t[:, 2]]) / 3.0
            t = t[np.where(fd(pmid) < -geps)]
            bars = np.concatenate((t[:, [0, 1]], t[:, [0, 2]], t[:, [1, 2]]))
            bars = np.unique(np.sort(bars), axis=0)
            if Iflag == 4:
                min_q, min_angle_deg = triqual_flag(p, t)
                print('Del iter: %3d, min q = %5.2f, min angle = %3.0f deg' \
                      % (Num_of_Delaunay_triangulations, min_q, min_angle_deg))
                if min_q > qmin:
                    break
            if Iflag == 2 or Iflag == 3:
                ktrimesh(p, bars)

        # move mesh points based on bar lengths L and forces F
        barvec = p[bars[:, 0], :] - p[bars[:, 1], :]
        L = np.sqrt(np.sum(barvec ** 2, axis=1))
        hbars = 0.5 * (fh(p[bars[:, 0], :]) + fh(p[bars[:, 1], :]))
        L0 = hbars * Fscale * np.sqrt(np.sum(L ** 2) / np.sum(hbars ** 2))
        F = np.maximum(L0 - L, 0)
        Fvec = np.column_stack((F, F)) * (barvec / np.column_stack((L, L)))
        Ftot = np.zeros((len(p), 2))
        n = len(bars)
        for j in range(n):
            Ftot[bars[j, 0], :] += Fvec[j, :]  # the : for the (x,y) components

            Ftot[bars[j, 1], :] -= Fvec[j, :]

        # force = 0 at fixed points, so they do not move:
        Ftot[0: len(pfix), :] = 0

        # update the node positions
        p = p + deltat * Ftot

        # bring outside points back to the boundary
        d = fd(p); ix = d > 0  # find points outside (d > 0)
        dpx = np.column_stack((p[ix, 0] + deps, p[ix, 1]))
        dgradx = (fd(dpx) - d[ix]) / deps
        dpy = np.column_stack((p[ix, 0], p[ix, 1] + deps))
        dgrady = (fd(dpy) - d[ix]) / deps
        p[ix, :] = p[ix, :] - np.column_stack((dgradx * d[ix], dgrady * d[ix]))

        # termination criterium: all interior nodes move less than dptol:

        if max(np.sqrt(np.sum(deltat * Ftot[d < -geps, :] ** 2, axis=1)) / h0) < dptol:
            break

    final_tri = Delaunay(p)  # another instantiation of the class
    t = final_tri.vertices
    pmid = (p[t[:, 0]] + p[t[:, 1]] + p[t[:, 2]]) / 3.0
    # keep the triangles whose geometrical center is inside the shape
    t = t[np.where(fd(pmid) < -geps)]
    bars = np.concatenate((t[:, [0, 1]], t[:, [0, 2]], t[:, [1, 2]]))
    # delete repeated bars
    # bars = unique_rows(np.sort(bars))
    bars = np.unique(np.sort(bars), axis=0)
    # orient all the triangles counterclockwise (ccw)
    t = ccw_tri(p, t)
    # graphical output of the current mesh
    ktrimesh(p, bars)
    triqual(p, t, fh)
    return p, t, bars


def boundary_bars(t):
    # create the bars (edges) of every triangle
    bars = np.concatenate((t[:, [0, 1]], t[:, [0, 2]], t[:, [1, 2]]))
    # sort all the bars
    data = np.sort(bars)
    # find the bars that are not repeated
    Delaunay_bars = dict()
    for row in data:
        row = tuple(row)
        if row in Delaunay_bars:
            Delaunay_bars[row] += 1
        else:
            Delaunay_bars[row] = 1
    # return the keys of Delaunay_bars whose value is 1 (non-repeated bars)
    bbars = []
    for key in Delaunay_bars:
        if Delaunay_bars[key] == 1:
            bbars.append(key)
    bbars = np.asarray(bbars)
    return bbars


def plot_shapes(xc, yc, r):
    # circle for plotting
    t_cir = np.linspace(0, 2 * pi)
    x_cir = xc + r * np.cos(t_cir)
    y_cir = yc + r * np.sin(t_cir)

    plt.figure()
    plt.plot(x_cir, y_cir)
    plt.grid()
    plt.title('Shapes')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    # plt.show()
    return


plt.close('all')

xc = 0; yc = 0; r = 1.0

x1, y1 = -1.0, -2.0
x2, y2 = 2.0, 3.0

plot_shapes(xc, yc, r)

xmin = -1.5; ymin = -1.5
xmax = 1.5; ymax = 1.5
h0 = 0.4

pfix = np.zeros((0, 2))

fd = Circle(xc,yc,r)

fh = lambda p: np.ones(len(p))

p,t,bars = distmesh(fd,fh,h0,xmin,ymin,xmax,ymax,pfix,Iflag=4)