#!/bin/bash
set -x -e

PATH=/Users/pauly/tk/cmrep/xc64rel:/Users/pauly/tk/pointset/xc64rel:$PATH

function set_specimen_vars()
{
  # Read the specimen ID
  local id=${1?}

  # Read the optional prefix
  local PF=${2}

  # Input dir and work dir
  INPUT_DIR=.
  WORK_DIR=.

  # Input mesh
  eval ${PF}INPUT_ADJUSTED_MESH_IMG=$INPUT_DIR/${id}_adjusted_mesh.nii.gz
  eval ${PF}INPUT_MRI=$INPUT_DIR/${id}_axisalign_img.nii.gz
  eval ${PF}INPUT_PHGSEG_IMG=$INPUT_DIR/${id}_axisalign_phgseg_multilabel.nii.gz

  # Input mesh to STL
  eval ${PF}MESH_ADJUSTED=$WORK_DIR/${id}_adjusted_mesh.stl

  # Sparse adjusted mesh (smoothed and unsmoothed)
  eval ${PF}MESH_SPARSE_UNSMOOTH=$WORK_DIR/${id}_sparse_unsmooth.vtk
  eval ${PF}MESH_SPARSE_SMOOTHED=$WORK_DIR/${id}_sparse_smooth.vtk

  # Padded PHGSEG with missing label removed
  eval ${PF}PHGSEG_PADDED=$WORK_DIR/${id}_phgseg_padded.nii.gz

  # Momenta for the smooth/unsmooth lmshoot
  eval ${PF}MOMENTA_SMOOTHED_TO_UNSMOOTH=$WORK_DIR/${id}_fx_smooth_mv_unsmooth_momenta.vtk
  eval ${PF}WARP_SMOOTHED_TO_UNSMOOTH=$WORK_DIR/${id}_fx_smooth_mv_unsmooth_warp.vtk
}

function set_specimen_pair_vars()
{
  local fid=${1?}
  local mid=${2?}

  # The variables will have prfixes FX_ and MV_, e.g., FX_PHGSEG_PADDED
  set_specimen_vars $fid "FX_"
  set_specimen_vars $mid "MV_"

  # Id for this pair
  PAIR_ID="fx_${fid}_mv_${mid}"

  # Pairwise directory for this registration
  PAIRWISE_DIR=${WORK_DIR}/reg_${PAIR_ID}

  # Affine transform
  PAIRWISE_SMOOTH_AFFINE_MATRIX=${WORK_DIR}/${PAIR_ID}_smoothed_affine.mat

  # Fixed smooth mesh transformed into the space of the moving mesh
  PAIRWISE_SMOOTH_AFFINE_FITTED_MESH=${WORK_DIR}/${PAIR_ID}_smoothed_fx_affine_to_mv.vtk

  # Momenta between the affine transformed fixed mesh and moving mesh
  PAIRWISE_SMOOTH_AFFINE_MOMENTA=${WORK_DIR}/${PAIR_ID}_smoothed_fx_affine_to_mv_momenta.vtk

}

function preprocess_specimen()
{
  local id=${1?}
  set_specimen_vars $id

  # Extract the mesh from input image
  vtklevelset $INPUT_ADJUSTED_MESH_IMG $MESH_ADJUSTED 0.5

  # Sample the mesh down to 2000 faces
  python dmca_sample.py -f 2000 -i $MESH_ADJUSTED $MESH_SPARSE_UNSMOOTH

  # Apply padding and missing label removal to the segmentation
  c3d $INPUT_PHGSEG_IMG -replace 2 0 -pad 5x5x5 5x5x5 0 -o $PHGSEG_PADDED

  # Sample the label into mesh space
  mesh_image_sample -V $MESH_SPARSE_UNSMOOTH $PHGSEG_PADDED $MESH_SPARSE_UNSMOOTH Label

  # Generate a smoothed mesh
  python dmca_sample.py -s 320 -l 0.8 -m -0.2 -f 2000 -i $MESH_ADJUSTED $MESH_SPARSE_SMOOTHED

  # Copy the anatomical labeled into smoothed mesh space
  mesh_merge_arrays -r $MESH_SPARSE_SMOOTHED $MESH_SPARSE_SMOOTHED Label $MESH_SPARSE_UNSMOOTH

  # Perform geodesic shooting from smoothed into raw space
  lmshoot -m $MESH_SPARSE_SMOOTHED $MESH_SPARSE_UNSMOOTH -o $MOMENTA_SMOOTHED_TO_UNSMOOTH \
    -s 2.0 -l 200 -d 3 -n 40 -i 600 0 -a L \
    -O /tmp/path%04d.vtk

  # Generate the corresponding warp
  lmtowarp -r $INPUT_MRI -o $WARP_SMOOTHED_TO_UNSMOOTH -m $MOMENTA_SMOOTHED_TO_UNSMOOTH -s 2.0 -d 3 -n 40
}

function pairwise_match_lmshoot()
{
  local fid=${1?}
  local mid=${2?}
  set_specimen_pair_vars $fid $mid

  # Create the working directory
  mkdir -p $PAIRWISE_DIR

  # Perform affine registration between fixed and moving
  ml_affine -m Label $FX_MESH_SPARSE_SMOOTHED $MV_MESH_SPARSE_SMOOTHED $PAIRWISE_SMOOTH_AFFINE_MATRIX

  # Apply the affine transformation to the fixed mesh
  warpmesh $FX_MESH_SPARSE_SMOOTHED $PAIRWISE_SMOOTH_AFFINE_FITTED_MESH $PAIRWISE_SMOOTH_AFFINE_MATRIX

  # Perform geodesic shooting between the latter and the moving mesh
  lmshoot -m $PAIRWISE_SMOOTH_AFFINE_FITTED_MESH $MV_MESH_SPARSE_SMOOTHED \
    -o $PAIRWISE_SMOOTH_AFFINE_MOMENTA \
    -s 3.0 -l 0.1 -d 3 -n 40 -i 300 0 -a C -S 1.0 \
    -O /tmp/pathx%04d.vtk
}

# Main entrypoint
CMD=$1
shift
$CMD "$@"



