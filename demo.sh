#!/usr/bin/env bash
# Both demonstration and integration testing for eg.
# new exceptions that might be thrown.
#
# This script runs the utils over the included demo data,
# then again over the reduced-size output, confirming
# that multiple .ply dialects (Meshlab, Pix4D) are handled.
#
# This data is very, very small (to fit in the git repo);
# multibillion-point clouds work - but make impractical demos.

echo "Running with Meshlab-style header and .xyz offset file"
echo
python main.py test_data/test_point_cloud.ply --savetrees test_data/trees
echo
# Move/rename output to disable recognition and reading of preprocessed file
mv test_point_cloud_sparse.ply test_data/second_point_cloud.ply
echo
echo "Running with first output (Pix4D-style header plus location comment)"
echo
forestutils test_data/second_point_cloud.ply
