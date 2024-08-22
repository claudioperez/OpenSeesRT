# ------------------------------------------------------------------------
# The following Python code is implemented by Professor Terje Haukaas at
# the University of British Columbia in Vancouver, Canada. It is made
# freely available online at terje.civil.ubc.ca together with notes,
# examples, and additional Python code. Please be cautious when using
# this code; it may contain bugs and comes without warranty of any form.
# ------------------------------------------------------------------------
# https://gist.githubusercontent.com/terjehaukaas/f633c4afc001badb4473d422ccc146e7/raw/2e7d09dbc850dc800c60e1751fb21f2f76615509/SolidCrossSectionAnalysis.py
# https://civil-terje.sites.olt.ubc.ca/files/2020/02/Screenshot-Solid-Cross-section.pdf
#
import numpy as np
from functools import partial

class TorsionAnalysis:
    pass

def assemble_matrix(Ka, ke, conn, nde, ndf):
    nne = len(conn)
    for j in range(nne):
        for k in range(j + 1):
            for m in range(nde):
                for n in range(nde):
                    Ka[(conn[j] - 1)*ndf + m, (conn[k] - 1)*ndf + n] += ke[j*nde + m, k*nde + n]

                    if j != k:
                        Ka[(conn[k] - 1)*ndf + n, (conn[j] - 1)*ndf + m] += ke[j*nde + m, k*nde + n]

    return Ka


def assemble_vector(Fa, fe, nodes, nde, ndf):
    nen = len(nodes)
    for j in range(nen):
        for m in range(nde):
            Fa[(nodes[j]-1) * ndf + m] += fe[j * nde + m]
    return Fa


class TorsionTriangle:
    _class_mass = 1/12*(np.eye(3) + np.ones((3,3)))

    def __init__(self, xyz):
        ((y1, y2, y3), (z1, z2, z3)) = xyz
        self.xyz = xyz
        self.area = area = 0.5 * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))

        z12 = z1 - z2
        z23 = z2 - z3
        z31 = z3 - z1
        y32 = y3 - y2
        y13 = y1 - y3
        y21 = y2 - y1

        area = 0.5 * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))

        k11 = ( y32**2 +  z23**2)
        k12 = (y13*y32 + z23*z31)
        k13 = (y21*y32 + z12*z23)
        k22 = ( y13**2 +  z31**2)
        k23 = (y13*y21 + z12*z31)
        k33 = ( y21**2 +  z12**2)
        ke = 1/(4.0*area)*np.array([[k11, k12, k13],
                                    [k12, k22, k23],
                                    [k13, k23, k33]])

        fe = -1/6.*np.array([
             ((y1*y32 - z1*z23) + (y2*y32 - z2*z23) + (y3*y32 - z3*z23)),
             ((y1*y13 - z1*z31) + (y2*y13 - z2*z31) + (y3*y13 - z3*z31)),
             ((y1*y21 - z1*z12) + (y2*y21 - z2*z12) + (y3*y21 - z3*z12))])

        self.load = fe
        self.stiffness = ke

    def getMass(self):
        rho = 1.0
        return rho*self.area*self._class_mass

    def getStiffness(self):
        return self.stiffness

    def getLoad(self):
        return self.load

def torsion_element(xyz):
    ((y1, y2, y3), (z1, z2, z3)) = xyz

    z12 = z1 - z2
    z23 = z2 - z3
    z31 = z3 - z1
    y32 = y3 - y2
    y13 = y1 - y3
    y21 = y2 - y1

    area = 0.5 * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))

    k11 = ( y32**2 +  z23**2)
    k12 = (y13*y32 + z23*z31)
    k13 = (y21*y32 + z12*z23)
    k22 = ( y13**2 +  z31**2)
    k23 = (y13*y21 + z12*z31)
    k33 = ( y21**2 +  z12**2)

    me = area/12*(np.eye(3) + np.ones((3,3)))

    ke = 1/(4.0*area)*np.array([[k11, k12, k13],
                                [k12, k22, k23],
                                [k13, k23, k33]])

    fe = -1/6.*np.array([
         ((y1*y32 - z1*z23) + (y2*y32 - z2*z23) + (y3*y32 - z3*z23)),
         ((y1*y13 - z1*z31) + (y2*y13 - z2*z31) + (y3*y13 - z3*z31)),
         ((y1*y21 - z1*z12) + (y2*y21 - z2*z12) + (y3*y21 - z3*z12))])

    return me, ke, fe

def _wrap_elem(centroid, nodes, conn):
    return conn, torsion_element((nodes[conn] - centroid).T)

import multiprocessing
def solve_torsion(section, mesh):
    ndf = 1
    nde = 1
    nodes = mesh.points
    nf = ndf*len(nodes)
    Ka = np.zeros((nf, nf))
    Ma = np.zeros((nf, nf))
    Fa = np.zeros(nf)

    connectivity = mesh.cells[1].data

    threads = 6
    chunk = 200
    with multiprocessing.Pool(threads) as pool:
        for conn, (me, ke, fe) in pool.imap_unordered(
                    partial(_wrap_elem, section.centroid, nodes),
                    mesh.cells[1].data,
                    chunk
            ):
            Ka  = assemble_matrix(Ka, ke, conn+1, nde, ndf)
            Ma  = assemble_matrix(Ma, me, conn+1, nde, ndf)
            Fa  = assemble_vector(Fa, fe, conn+1, nde, ndf)

#   queue = ((conn, torsion_element((nodes[conn] - section.centroid).T))
#            for conn in connectivity)

#   for conn, (me, ke, fe) in queue:
#       Ka  = assemble_matrix(Ka, ke, conn+1, nde, ndf)
#       Ma  = assemble_matrix(Ma, me, conn+1, nde, ndf)
#       Fa  = assemble_vector(Fa, fe, conn+1, nde, ndf)

    # Lock the solution at one node and solve for the others
    Pf = Fa[:nf-1]
    for i in range(nf-1):
        Pf[i] -= Ka[i, nf-1]

    Kf = Ka[:nf-1, :nf-1]
    Uf = np.linalg.solve(Kf, Pf)
    ua = np.append(Uf, 1.0)
    #return ua


#def shear_center(section, mesh, ua):
    # Determine shear center coordinates
    A  = section.area
    Iy = section.iyc
    Iz = section.ixc
    warpIntegral = 0
    yscIntegral = 0
    zscIntegral = 0
    for conn in connectivity:
        ((y1, y2, y3), (z1, z2, z3)) = (nodes[conn] - section.centroid).T

        # Warping values
        u1, u2, u3 = ua[conn]
        area = 0.5 * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))
        warpIntegral += (u1+u2+u3)/3.0 * area
        yscIntegral  += area/12.0*(u1*(2*z1+z2+z3) + u2*(z1+2*z2+z3) + u3*(z1+z2+2*z3))
        zscIntegral  += area/12.0*(u1*(2*y1+y2+y3) + u2*(y1+2*y2+y3) + u3*(y1+y2+2*y3))

#def normalize_shear():
    thetaPrincipal = 0.0
    # Normalizing constant and shear center coordinates
    # normalizingConstant = -warpIntegral / A
    if np.abs(thetaPrincipal) > 1e-3:
        ysc = -yscIntegral / IyRotated
        zsc =  zscIntegral / IzRotated

    else:
        ysc = -yscIntegral / Iy
        zsc =  zscIntegral / Iz

    # Finalize
    finalWarp =   ua[:nf] - warpIntegral / A \
                + np.array([ysc, -zsc])@(mesh.points[:nf] - section.centroid).T

    return np.array(finalWarp)


def warping_constant(section, mesh, warp, sc):
    # ------------------------------------------------------------------------
    # WARPING CONSTANT
    # ------------------------------------------------------------------------

    nodes        = mesh.points
    connectivity = mesh.cells[1].data

    Cw = 0.0
    for conn in connectivity:
        ((y1, y2, y3), (z1, z2, z3)) = (nodes[conn] - section.centroid).T

        # Warping values
        u1, u2, u3 = warp[conn]

        # Warping constant
        dA  = 0.5 * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))
        Cw += dA / 6.0 * (u1**2 + u2**2 + u3**2 + u1*u2 + u1*u3 + u2*u3)

    return Cw

def torsion_constant(section, mesh, warp, sc):
    # ------------------------------------------------------------------------
    # ST. VENANT CONSTANT
    # ------------------------------------------------------------------------

    # for i in range(nf):
    #     yValue = (nodes[i])[0]
    #     zValue = (nodes[i])[1]
    #     warp.append(ua[i] + normalizingConstant + ysc*(zValue-zCentroid) - zsc*(yValue-yCentroid))

    nodes        = mesh.points
    connectivity = mesh.cells[1].data

    # Shear center coordinates in original axis system
    zCentroid, yCentroid = section.centroid

    ysc, zsc = sc
    ysc = yCentroid + ysc
    zsc = zCentroid + zsc

    J  = 0.0
    for conn in connectivity:
        ((y1, y2, y3), (z1, z2, z3)) = (nodes[conn] - section.centroid).T

        z23 = z2 - z3
        z31 = z3 - z1
        z12 = z1 - z2
        y32 = y3 - y2
        y13 = y1 - y3
        y21 = y2 - y1

        # Warping values
        u1, u2, u3 = warp[conn]

        # Element area
        dA  = 0.5 * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))

        # St. Venant constant
        Czeta1  = ( u2*y1 * y13 + u3 *  y1 * y21 + u1 * y1*y32 - u3 * z1 * z12 - u1*z1 * z23 - u2*z1*z31)/(2*dA)
        Czeta2  = (u2*y13 *  y2 + u3 *  y2 * y21 + u1 * y2*y32 - u3 * z12 * z2 - u1*z2 * z23 - u2*z2*z31)/(2*dA)
        Czeta3  = (u2*y13 *  y3 + u3 * y21 *  y3 + u1 * y3*y32 - u3 * z12 * z3 - u1*z23 * z3 - u2*z3*z31)/(2*dA)
        Czeta12 = 2*y1*y2 + 2*z1*z2
        Czeta13 = 2*y1*y3 + 2*z1*z3
        Czeta23 = 2*y2*y3 + 2*z2*z3
        Czeta1s =   y1**2 +   z1**2
        Czeta2s =   y2**2 +   z2**2
        Czeta3s =   y3**2 +   z3**2
        J += ((Czeta1+Czeta2+Czeta3)/3. \
            + (Czeta12+Czeta13+Czeta23)/12. \
            + (Czeta1s+Czeta2s+Czeta3s)/6.)*dA

    return J

#
#  # ------------------------------------------------------------------------
#  # ROTATE AXES BACK
#  # ------------------------------------------------------------------------
#  
#  if np.abs(thetaPrincipal) > 1e-3:
#      yCentroid = storeYcentroid
#      zCentroid = storeZcentroid
#      yscOff = ysc
#      zscOff = zsc
#      ysc = yscOff * np.cos(thetaPrincipal) - zscOff * np.sin(thetaPrincipal) + yCentroid
#      zsc = yscOff * np.sin(thetaPrincipal) + zscOff * np.cos(thetaPrincipal) + zCentroid
#  
#      for i in range(len(nodes)):
#          yOld = (nodes[i])[0]
#          zOld = (nodes[i])[1]
#          (nodes[i])[0] = yOld * np.cos(thetaPrincipal) - zOld * np.sin(thetaPrincipal) + yCentroid
#          (nodes[i])[1] = yOld * np.sin(thetaPrincipal) + zOld * np.cos(thetaPrincipal) + zCentroid
#
#
#  # ------------------------------------------------------------------------
#  # PRINT RESULTS
#  # ------------------------------------------------------------------------
#  
#  print('\n'"Area:", A)
#  print('\n'"Centroid y-coordinate: %.2f" % yCentroid)
#  print('\n'"Centroid z-coordinate: %.2f" % zCentroid)
#  print('\n'"Moment of inertia Iy: %.2f" % Iy)
#  print('\n'"Moment of inertia Iz: %.2f" % Iz)
#  if np.abs(thetaPrincipal) > 1e-3:
#      print('\n'"Principal axis rotation (degrees): %.2f" % (thetaPrincipal/np.pi*180.0))
#      print('\n'"Product of inertia Iyz: %.2f" % Iyz)
#      print('\n'"New principal moment of inertia Iy: %.2f" % IyRotated)
#      print('\n'"New principal moment of inertia Iz: %.2f" % IzRotated)
#      print('\n'"Shear center y-coordinate in original coordinate system: %.2f" % ysc)
#      print('\n'"Shear center z-coordinate in original coordinate system: %.2f" % zsc)
#  else:
#      print('\n'"Shear center y-coordinate: %.2f" % ysc)
#      print('\n'"Shear center z-coordinate: %.2f" % zsc)
#  print('\n'"St. Venant torsion constant J: %.2f" % J)
#  print('\n'"Warping torsion constant Cw: %.2f" % Cw)
#  
#  # ------------------------------------------------------------------------
#  # PLOT
#  # ------------------------------------------------------------------------
#  
#  # Determine warping bounds
#  maxWarp = np.amax(finalWarp)
#  minWarp = np.amin(finalWarp)
#  
#  # Plot the edges of the triangles
#  alreadyPlotted = []
#  for i in range(len(connectivity)):
#  
#      # Node numbers
#      node1 = (connectivity[i])[0]
#      node2 = (connectivity[i])[1]
#      node3 = (connectivity[i])[2]
#  
#      # y-coordinates of vertices
#      y1 = (nodes[node1])[0]
#      y2 = (nodes[node2])[0]
#      y3 = (nodes[node3])[0]
#  
#      # y-coordinates of vertices
#      z1 = (nodes[node1])[1]
#      z2 = (nodes[node2])[1]
#      z3 = (nodes[node3])[1]
#  
#      # Warping value at nodes
#      warp1 = finalWarp[node1]
#      warp2 = finalWarp[node2]
#      warp3 = finalWarp[node3]
#  
#      # Plot first edge
#      if [node1, node2] not in alreadyPlotted:
#          warpAverage = (warp1+warp2)/2.0
#          colorCode = (warpAverage - minWarp) / (maxWarp - minWarp)
#          plt.plot(np.array([y1, y2]), np.array([z1, z2]), color=(colorCode, 0, 1-colorCode), zorder=1)
#          alreadyPlotted.append([node1, node2])
#  
#      # Plot second edge
#      if [node1, node3] not in alreadyPlotted:
#          warpAverage = (warp1+warp3)/2.0
#          colorCode = (warpAverage - minWarp) / (maxWarp - minWarp)
#          plt.plot(np.array([y1, y3]), np.array([z1, z3]), color=(colorCode, 0, 1-colorCode), zorder=1)
#          alreadyPlotted.append([node1, node3])
#  
#      # Plot third edge
#      if [node2, node3] not in alreadyPlotted:
#          warpAverage = (warp2+warp3)/2.0
#          colorCode = (warpAverage - minWarp) / (maxWarp - minWarp)
#          plt.plot(np.array([y2, y3]), np.array([z2, z3]), color=(colorCode, 0, 1-colorCode), zorder=1)
#          alreadyPlotted.append([node2, node3])
#  
#  # Create the plot
#  max = 0
#  min = 0
#  for i in range(len(y)):
#      theMax = np.max(y[i])
#      theMin = np.min(y[i])
#      if theMax > max:
#          max = theMax
#      if theMin < min:
#          min = theMin
#  for i in range(len(z)):
#      theMax = np.max(z[i])
#      theMin = np.min(z[i])
#      if theMax > max:
#          max = theMax
#      if theMin < min:
#          min = theMin
#  yMin = min - (max-min)*0.4
#  border = (max - min) * 0.1
#  plt.plot(np.array([yCentroid, max+border]), np.array([zCentroid, zCentroid+(np.abs(max)+border-yCentroid)*np.tan(thetaPrincipal)]), 'k-',linewidth=1.0)
#  plt.plot(np.array([yCentroid, yMin-border]), np.array([zCentroid, zCentroid-(np.abs(yMin)+border+yCentroid)*np.tan(thetaPrincipal)]), 'k-',linewidth=1.0)
#  plt.plot(np.array([yCentroid, yCentroid-(np.abs(max)+border-zCentroid)*np.tan(thetaPrincipal)]), np.array([zCentroid, max+border]), 'k-',linewidth=1.0)
#  plt.plot(np.array([yCentroid, yCentroid+(np.abs(min)+border+zCentroid)*np.tan(thetaPrincipal)]), np.array([zCentroid, min-border]), 'k-',linewidth=1.0)
#  plt.scatter([yCentroid], [zCentroid], s=100, c='k', marker='o', zorder=2)
#  plt.scatter([ysc], [zsc], s=100, c='k', marker='o', zorder=2)
#  plt.axis([yMin-border, max+border, min-border, max+border])
#  plt.title("Centroid + Shear Centre + Warping")
#  plt.show()
#  
def plot(mesh, values=None, scale=1.0, show_edges=None, savefig:str=None,**kwds):
    from matplotlib import cm
    import pyvista as pv
    pv.set_jupyter_backend("panel")


    pv.start_xvfb(wait=0.05)
    mesh = pv.utilities.from_meshio(mesh)
    if values is not None:
        point_values = scale*values
        mesh.point_data["u"] = point_values
        mesh = mesh.warp_by_scalar("u", factor=scale)
        mesh.set_active_scalars("u")
    if show_edges is None:
        show_edges = True #if sum(len(c.data) for c in mesh.cells) < 1000 else False

    if not pv.OFF_SCREEN:
        plotter = pv.Plotter(notebook=True)
        plotter.add_mesh(mesh,
           show_edges=show_edges,
           cmap=cm.get_cmap("RdYlBu_r"),
           lighting=False,
           **kwds
        )
        # if len(values) < 1000:
        #     plotter.add_mesh(
        #        pv.PolyData(mesh.points), color='red',
        #        point_size=5, render_points_as_spheres=True)
        if savefig:
            plotter.show(screenshot=savefig)
        else:
            plotter.show()


