import pymeshlab
import os
import sys
import argparse
import numpy as np
import pathlib
from dmcalib import export_mesh
from scipy.spatial import KDTree, procrustes

# Add the arguments
parse = argparse.ArgumentParser(
    description="Simplify mesh with corredpondences")
parse.add_argument('Source', metavar='source', type=str, help='input mesh')
parse.add_argument('Result', metavar='result', type=str, help='output text file')
parse.add_argument('-n', type=int, help='number of samples', default=1000, dest='NumSamples')

# Do the parsing
args = parse.parse_args()

# Load the mesh
ms = pymeshlab.MeshSet()
ms.load_new_mesh(args.Source)

# Sample a relatively small number of samples from the final smoothed mesh
print("Sampling with %d points" % args.NumSamples)
ms.poisson_disk_sampling(samplenum=args.NumSamples, subsample=True, exactnumflag=True)

# Save depending on format
vms = ms.current_mesh().vertex_matrix()
export_mesh(args.Result, vms)
