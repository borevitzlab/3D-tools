#!/usr/bin/env python3
"""Convert image collections to pointclouds with Pix4D, https://pix4d.com/

Somewhat specialised for use in the ANU Borevitz lab, but easily modified.
"""

import argparse
import os
import subprocess
import sys

PIX4D_PATH = r'C:\Program Files\Pix4Dmapper\pix4dmapper.exe'

# TODO:  put correct paths here
NAT_ARB_OUTPUT_DIR = '.'
OTHER_OUTPUT_DIR = '.'

# Image types recognised by Pix4D
IMAGE_EXTS = ['.jpg', '.jpeg', '.tif', '.tiff']

def get_args():
    """Get the few arguments used by this script."""
    parser = argparse.ArgumentParser(
        description=('Run Pix4D batch jobs in the Boreviz Lab.'))
    parser.add_argument('input',
        help='A list of input paths.  Must include one or more image'
        'directories, and may contain one location file (eg flight log).')
    parser.add_argument('-a', '--arboretum', help='location is at the national'
        'arboretum.  Default false', default=False, action='store_true')
    parser.add_argument('--locationfile', help='image georeference file, '
        'eg flight log.')
    return parser.parse_args()

def validate_locat_file(fname):
    """Return a tuple of Pix4D georeference type and filename.  Return None if
    filetype is invalid or unrecognised."""
    #3drobotics: for "3D Robotics" flight log.
    #sensefly: for "senseFly" flight log.
    #pix4d-lat-long: for "Pix4D Latitude - Longitude" file format.
    #pix4d-long-lat: for "Pix4D Longitude - Latitude" file format.
    #pix4d-x-y: For "Pix4D X - Y" file format.
    #pix4d-y-x: For "Pix4D Y - X" file format.
    return None

def validate_img_dir(dirname):
    """Return True if the directory contains images."""
    for p in os.listdir(dirname):
        if any(p.ends_with(ext) for ext in IMAGE_EXTS):
            return True

def consolidate_images(dirs):
    """Return a single directory with all images.  If multiple dirs are passed,
    a tempdir is created and images copied across."""
    if not dirs:
        raise RuntimeError('No image dir supplied')
    if len(dirs) == 1:
        return dirs[0]
    # TODO: finish this function
    raise NotImplementedError("Can't handle multiple image dirs yet")

def parse_input_paths(input_paths):
    """Return a list of image dirs, and a location file or None if the latter
    does not exist."""
    dirs = []
    locat_file = None
    for p in input_paths:
        if os.path.isfile(p):
            locat_file = validate_locat_file(p)
        elif os.path.isdir(p):
            if validate_img_dir(p):
                dirs.append(p)
            else:
                print('No images in {}'.format(p))
        else:
            print('Nothing at input path {}'.format(p))
    return consolidate_images(dirs), locat_file

def call_pix4d(proj_file):
    """Work out all the arguments and call Pix4D, using our usual settings
    in background mode."""
    p_args = [
        PIX4D_PATH, '--cmdline', '--full', # background, full processing
        '--initial', '--dense', # setup project and create dense pointcloud
        '--feature-scale', '1', # compute features at original image size
        '--dense-image-scale', '0',   # densify with full images
        '--dense-point-density', '1', # highest pointcloud density
        '--dense-use-gpu', 'yes',
        '--mesh-generate', 'no', # do not generate mesh
        '--max-cpus', '11',
        '--image-dir'
        ]
    image_dir, locat_file = parse_input_paths(args.input)
    p_args.append(image_dir)
    if locat_file is not None:
        p_args.extend(['--geolocation-format', locat_file[0],
                       '--geolocation-file', locat_file[1]])
    # Finally, decide where to save the output
    proj_file = os.path.join(
        NAT_ARB_OUTPUT_DIR if args.arboretum else OTHER_OUTPUT_DIR, proj_file)
    if os.path.isfile(proj_file) or not proj_file.endswith('.p4d'):
        raise RuntimeError('Project file is invalid or already exists.')
    else:
        p_args.append(proj_file)
    #... and run the program!
    subprocess.call(p_args)

if __name__ == '__main__':
    args = get_args()
    print(args)
    if not os.path.isfile(PIX4D_PATH):
        sys.exit("Pix4D not found - check installation and try again.")
    call_pix4d('Test_project.p4d')
