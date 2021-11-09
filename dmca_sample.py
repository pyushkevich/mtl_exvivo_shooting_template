import pymeshlab
import os
import sys
import argparse
import numpy as np
import trimesh
from scipy.spatial import KDTree, procrustes
from dmcalib import export_mesh

# Add the arguments
parse = argparse.ArgumentParser(
    description="Simplify mesh with corredpondences")
parse.add_argument('Source', metavar='source', type=str, help='input mesh')
parse.add_argument('Result', metavar='result', type=str, help='output samples')
parse.add_argument('-f', type=int, help='Target number of faces', default=2000, dest='FaceNum')
parse.add_argument('-s', type=int, help='number of saved smoothing steps', default=0, dest='TaubinSteps')
parse.add_argument('-l', type=float, help='Taubin lambda parameter', default=0.6, dest='TaubinLambda')
parse.add_argument('-m', type=float, help='Taubin mu parameter', default=-0.4, dest='TaubinMu')
parse.add_argument('-i', dest='InvertNormals', help='Invert normal orientation', action='store_true')

# Do the parsing
args = parse.parse_args()

# Load the mesh
ms = pymeshlab.MeshSet()
ms.load_new_mesh(args.Source)

mesh_stat = lambda ms : (ms.current_mesh().vertex_number(), ms.current_mesh().face_number())
print("Input mesh: %d vertices, %d faces" % mesh_stat(ms))

# Store the vertices and faces of the dense mesh
vmd = ms.current_mesh().vertex_matrix()
fmd = ms.current_mesh().face_matrix()

# Simplify the mesh
ms2 = pymeshlab.MeshSet()
ms2.add_mesh(ms.current_mesh())

ms2.simplification_quadric_edge_collapse_decimation(targetfacenum=args.FaceNum,
                                                    preserveboundary=True,
                                                    preservenormal=True,
                                                    preservetopology=True,
                                                    planarquadric=True)

# Optionally invert the face orientation
if args.InvertNormals:
    ms2.invert_faces_orientation()

# Store the vertices and faces of the current (simplified) mesh
vms = ms2.current_mesh().vertex_matrix()
fms = ms2.current_mesh().face_matrix()

print("Sparse mesh without smoothing: %d vertices, %d faces" % mesh_stat(ms))

if args.TaubinSteps == 0:

    # Just save the mesh as it is
    export_mesh(args.Result, 
                ms2.current_mesh().vertex_matrix(),
                ms2.current_mesh().face_matrix())

else:

    # Smooth the dense mesh
    ms.taubin_smooth(lambda_=args.TaubinLambda, mu=args.TaubinMu, stepsmoothnum=args.TaubinSteps)
    vmr = ms.current_mesh().vertex_matrix()
    fmr = ms.current_mesh().face_matrix()
    print("Smoothed input mesh: %d vertices, %d faces" % mesh_stat(ms))

    # Create a VTK cell locator
    tm_dense = trimesh.Trimesh(vertices=vmd, faces=fmd)
    (tm_closest, tm_dist, tm_tid) = trimesh.proximity.closest_point(tm_dense, vms)

    # Remap the triangles
    vmo = vms * 0.0
    for i in range(vms.shape[0]):
        X = np.zeros((3, 3))
        X[:, 0] = vmd[int(fmd[tm_tid[i], 0]), :]
        X[:, 1] = vmd[int(fmd[tm_tid[i], 1]), :]
        X[:, 2] = vmd[int(fmd[tm_tid[i], 2]), :]
        b = tm_closest[i, :]
        w = np.linalg.solve(X, b)

        X[:, 0] = vmr[int(fmd[tm_tid[i], 0]), :]
        X[:, 1] = vmr[int(fmd[tm_tid[i], 1]), :]
        X[:, 2] = vmr[int(fmd[tm_tid[i], 2]), :]
        z = np.dot(X, w)
        # print(vms[i, :], ' -> ', tm_tid[i], z, tm_dist[i])

        vmo[i] = z

    # Add a tiny amount of randomness to prevent vertices being merged on save
    vmo = vmo + np.random.randn(vmo.shape[0], 3) * 1.0e-5

    # Save the sparse mesh matched to reference
    print("Smoothed output mesh: %d vertices, %d faces" % (vmo.shape[0], fms.shape[0]))
    export_mesh(args.Result, vmo, fms)

# Perform smoothing on the dense mesh and generate corresponding samples
#
#
#
# # If reference specified, find corresponding vertices relative to it
# if args.Reference is not None:
#
#     # Store the vertices and faces of the current (simplified) mesh
#     vms = ms.current_mesh().vertex_matrix()
#     fms = ms.current_mesh().face_matrix()
#
#     # Create a VTK cell locator
#     tm_dense = trimesh.Trimesh(vertices=vmd, faces=fmd)
#     (tm_closest, tm_dist, tm_tid) = trimesh.proximity.closest_point(tm_dense, vms)
#
#     # Load the reference dense mesh
#     msref = pymeshlab.MeshSet()
#     msref.load_new_mesh(args.Reference[0])
#     vmr = msref.current_mesh().vertex_matrix()
#     print("Reference mesh: %d vertices, %d faces" % mesh_stat(msref))
#
#     # Remap the triangles
#     vmo = vms * 0.0
#     for i in range(vms.shape[0]):
#         X=np.zeros((3,3))
#         X[:,0] = vmd[int(fmd[tm_tid[i],0]),:]
#         X[:,1] = vmd[int(fmd[tm_tid[i],1]),:]
#         X[:,2] = vmd[int(fmd[tm_tid[i],2]),:]
#         b = tm_closest[i,:]
#         w = np.linalg.solve(X, b)
#
#         X[:,0] = vmr[int(fmd[tm_tid[i],0]),:]
#         X[:,1] = vmr[int(fmd[tm_tid[i],1]),:]
#         X[:,2] = vmr[int(fmd[tm_tid[i],2]),:]
#         z = np.dot(X, w)
#         print(vms[i,:], ' -> ', tm_tid[i], z, tm_dist[i])
#
#         vmo[i] = z
#
#     # v_pts = vtk.vtkPoints()
#     # v_pts.SetData(numpy_support.numpy_to_vtk(vmd))
#     # v_poly = vtk.vtkPolyData()
#     # v_poly.SetPoints(v_pts)
#     # v_cells = vtk.vtkCellArray()
#     # for i in range(fmd.shape[0]):
#     #     tri = vtk.vtkTriangle()
#     #     tri.GetPointIds().SetId(0, fmd[i,0])
#     #     tri.GetPointIds().SetId(1, fmd[i,1])
#     #     tri.GetPointIds().SetId(2, fmd[i,2])
#     #     v_cells.InsertNextCell(tri)
#     # v_poly.SetPolys(v_cells)
#     # v_poly.BuildCells()
#     # v_poly.BuildLinks()
#     # v_writer = vtk.vtkSTLWriter()
#     # v_writer.SetFileName('/tmp/test.stl')
#     # v_writer.SetInputData(v_poly)
#     # # v_writer.SetFileTypeToASCII()
#     # v_writer.Update()
#     #
#     # # Create the actual locator
#     # # locator = vtk.vtkCellLocator()
#     # locator = vtk.vtkOBBTree()
#     # locator.SetDataSet(v_poly)
#     # locator.BuildLocator()
#     #
#     # # For each node, localize it in the sparse mesh
#     # vmo = vms * 0.0
#     # tree = KDTree(vmd)
#     # for i in range(vms.shape[0]):
#     #     cellId = vtk.reference(0)
#     #     c = [0.0, 0.0, 0.0]
#     #     subId = vtk.reference(0)
#     #     d = vtk.reference(0.0)
#     #     qid = locator.FindCell([vms[i,0],vms[i,1],vms[i,2]])
#     #     locator.FindClosestPoint([vms[i,0],vms[i,1],vms[i,2]], c, cellId, subId, d)
#     #     X=np.zeros((3,3))
#     #     X[:,0] = vmd[int(fmd[cellId,0]),:]
#     #     X[:,1] = vmd[int(fmd[cellId,1]),:]
#     #     X[:,2] = vmd[int(fmd[cellId,2]),:]
#     #     b = np.array(c)
#     #     w = np.linalg.solve(X, b)
#     #
#     #     X[:,0] = vmr[int(fmd[cellId,0]),:]
#     #     X[:,1] = vmr[int(fmd[cellId,1]),:]
#     #     X[:,2] = vmr[int(fmd[cellId,2]),:]
#     #     z = np.dot(X, w)
#     #     print(vms[i,:], ' -> ', cellId, z)
#     #
#     #     vmo[i] = z
#
#     # Save the sparse mesh matched to reference
#     m_new = pymeshlab.Mesh(vmo, fms)
#     msref.add_mesh(m_new)
#     print("Sparse in reference space: %d vertices, %d faces" % mesh_stat(msref))
#     msref.save_current_mesh(args.Reference[1])
#
#
