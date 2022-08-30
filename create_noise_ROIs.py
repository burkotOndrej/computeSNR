#!/usr/bin/python3
#
# Create 4 cubic ROIs in the input image corners. The ROIs can be used for computation of background noise SD.
#
# Example:
#       python3 create_noise_ROIs.py -i <input_image> -shiftx 10 -shifty 10 -shiftunits px -size 15
#
# Authors: Ondrej Burkot, Jan Valosek
#

import os
import sys
import argparse
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

from argparse import RawDescriptionHelpFormatter


def get_parser():
    parser = argparse.ArgumentParser(description='Create 4 cubic ROIs in the input image corners.'
                                                 '\nThe ROIs can be used for computation of background noise SD.',
                                     formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-i', help='Input image', type=str, required=True)
    parser.add_argument('-shiftx', help='Shift of the ROIs in x axis from the image borders.', type=int,
                        required=False, default=10)
    parser.add_argument('-shifty', help='Shift of the ROIs in y axis from the image borders.', type=int,
                        required=False, default=10)
    parser.add_argument('-shiftunits', help='Shift units in pixels (pix) or percentage (per)', type=str,
                        required=False, default='pix', choices=['pix', 'per'])
    parser.add_argument('-size', help='Size of the ROIs in pixels. Example: 10', type=int, required=False, default=10)
    return parser


def main():

    # Parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # Fetch input arguments
    path_input_file = os.path.abspath(args.i)

    # Check if input file exists
    if os.path.exists(path_input_file):
        testImg1 = nib.load(path_input_file)
    else:
        print('ERROR: Input file {} does not exist.'.format(path_input_file))
        sys.exit()

    size = args.size  # size of the square
    shift_x = args.shiftx
    shift_y = args.shifty
    shift_units = args.shiftunits

    # testImg1 = testImg1.get_fdata()[:, :, 2]
    positionShift = 5                                   # value in %, should cover noise in background
    slice = 2
    color = 'Blues'

    testImg1Hdr = testImg1.header
    nx, ny, nz = testImg1Hdr.get_data_shape()           # data shape, number of pixels in x, y, z

    if shift_units == 'pix':
        # center of ROI [x, y, z]
        center = [np.ceil(nx - shift_x),
                  np.ceil(ny - shift_y),
                  nz]
    elif shift_units == 'per':
        # center of ROI [x, y, z]
        center = [np.ceil(nx*shift_x/100),
                  np.ceil(ny*shift_y/100),
                  nz]

    # mask for 3D volume
    if nz > 1:

        # initialize 2d grid
        xx, yy, zz = np.mgrid[:nx, :ny, :nz]
        xc = center[0]
        yc = center[1]
        zc = center[2]

        radius_x = np.ceil((int(size) - 1) / 2)         # radius of square - 5 to left, 5 to right
        radius_y = radius_x                             # square --> radius left-right and top-bottom are same
        radius_z = nz

        # * 1 to convert F or T to 0 or 1
        mask3d = ((abs(xx - xc) <= radius_x) &
                  (abs(yy - yc) <= radius_y) &
                  (abs(zz - zc) <= radius_z)) * 1

        plt.figure()
        plt.imshow(testImg1.get_fdata()[:, :, slice].T, cmap='gray')
        plt.imshow(np.transpose(mask3d[:, :, slice]), cmap='Blues', alpha=0.3)
        plt.show()

    # mask for single 2D slice
    else:

        # initialize 2d grid
        xx, yy = np.mgrid[:nx, :ny]
        xc = center[0]
        yc = center[1]

        radius_x = np.ceil((int(size) - 1) / 2)         # radius of rectangular - 5 to left, 5 to right
        radius_y = radius_x                             # rectangular --> radius left-right and top-bottom are same

        # * 1 to convert F or T to 0 or 1
        mask2d = ((abs(xx - xc) <= radius_x) &
                  (abs(yy - yc) <= radius_y)) * 1

        plt.figure()
        plt.imshow(testImg1.get_fdata()[:, :, slice].T, cmap='gray')
        plt.imshow(np.transpose(mask2d[:, :]), cmap=color, alpha=0.3)
        plt.show()

    # testImg2 = nib.load('sub-5147B_T2star.nii.gz')
    # testImg1Affine = testImg1.affine.copy()
    # testImg1Header = testImg1.header.copy()


if __name__ == "__main__":
    main()
