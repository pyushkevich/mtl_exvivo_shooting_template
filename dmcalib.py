import pymeshlab
import os
import vtk
import trimesh
import pathlib
from vtk.util import numpy_support

def export_mesh(filename, vertices, faces=None):
    """
    Export a mesh (vertices, faces) to a STL or VTK file

    :param filename: Output filename
    :param vertices: Nx3 numpy array of vertex coordinates
    :param faces: Nx3 numpy integer array of face indices
    """
    if pathlib.Path(filename).suffix == '.ply':
        np.savetxt(filename, vertices)

    elif pathlib.Path(filename).suffix == '.vtk':

        # Set VTK points
        v_pts = vtk.vtkPoints()
        v_pts.SetData(numpy_support.numpy_to_vtk(vertices))
        v_poly = vtk.vtkPolyData()
        v_poly.SetPoints(v_pts)

        # Set VTK faces
        if faces is not None:
            v_cells = vtk.vtkCellArray()
            for i in range(faces.shape[0]):
                tri = vtk.vtkTriangle()
                for j in range(3):
                    tri.GetPointIds().SetId(j, faces[i,j])
                v_cells.InsertNextCell(tri)
            v_poly.SetPolys(v_cells)

        # Export to VTK
        v_writer = vtk.vtkPolyDataWriter()
        v_writer.SetFileVersion(vtk.vtkPolyDataWriter.VTK_LEGACY_READER_VERSION_4_2)
        v_writer.SetFileName(filename)
        v_writer.SetInputData(v_poly)
        v_writer.Update()

    else:
        tm = trimesh.Trimesh(vertices=vertices, faces=faces)
        tm.export(filename)

