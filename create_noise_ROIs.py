#!/usr/bin/python3
#
# Create 4 cubic ROIs in the input image corners. The ROIs can be used for computation of background noise SD.
# Input image should have same number of rows and columns, otherwise script will fail.
#
# Example:
#       python3 create_noise_ROIs.py -i <input_image> -shiftx 10 -shifty 10 -shiftunits pix -size 15
#
# Authors: Ondrej Burkot, Jan Valosek
#

import os
import sys
import argparse
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import operator
from functools import reduce

from argparse import RawDescriptionHelpFormatter

# prejmenovat ulozeny obrazek na jmeno vstupniho souboru.png

filtration_dict = {
    'max': 0,
    'iqr': 1
}

def get_parser():
    parser = argparse.ArgumentParser(description='Create 4 cubic ROIs in the input image corners.'
                                                 '\nThe ROIs can be used for computation of background noise SD.',
                                     formatter_class=RawDescriptionHelpFormatter,
                                     add_help=False)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-i', help='Input image', type=str, required=True)
    optional.add_argument('-shiftx', help='Shift of the ROIs in x axis from the image borders.', type=int,
                        required=False, default=10)
    optional.add_argument('-shifty', help='Shift of the ROIs in y axis from the image borders.', type=int,
                        required=False, default=10)
    optional.add_argument('-shiftunits', help='Shift units in pixels (pix) or percentage (per)', type=str,
                        required=False, default='pix', choices=['pix', 'per'])
    optional.add_argument('-size', help='Size of the ROIs in pixels. Example: 10', type=int, required=False, default=10)
    optional.add_argument('-filter', help='Filter method', default='max', choices=['max', 'iqr', 'None'])
    optional.add_argument('-visualise', help='Visualisation of created ROI for debug. 0 - do not visualise, 1 - visualise',
                        type=int, required=False, default=0, choices=[0, 1])
    optional.add_argument('-outpath', help='Path for saving nii file with created ROIs. Example: /home/<usr>/<directory>'
                        , type=str, required=False, default='/')
    optional.add_argument('-h', '--help', action='help',
                          help='Show this message and exit.')
    return parser

def main():
    # Parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # Fetch input arguments
    path_input_file = os.path.abspath(args.i)

    # Check if input file exists
    try:
        file_data = nib.load(path_input_file)
    except:
        print('ERROR: Input image not found, exitting script!')
        sys.exit()

    size = args.size  # size of the square
    shift_x = args.shiftx
    shift_y = args.shifty
    shift_units = args.shiftunits
    visualise = args.visualise
    path_output = args.outpath

    if os.path.exists(path_output) is False:
        print('ERROR: Output path does not exist!')
        sys.exit()

    if args.filter != 'None':
        filter_method = filtration_dict[args.filter]

    test_img1_hdr = file_data.header
    test_img1_affine = file_data.affine

    # Get input image size
    im_size = test_img1_hdr.get_data_shape()  # data shape, number of pixels in x, y, z

    if im_size[0] != im_size[1]:
        print('ERROR: Input data does not have same number of rows and columns!')
        sys.exit()

    # 2D - if the input image is just single slice (i.e., 2D), set the third dimension to 1
    if len(im_size) == 2:
        nz = 1
    # 3D
    else:
        nz = im_size[2]
    nx, ny = im_size[0], im_size[1]  # data shape, number of pixels in x, y, z

    # fetch image data (numpy ndarray)
    file_data = file_data.get_fdata()

    if shift_units == 'pix':
        # center of ROI [x, y, z]
        center = [np.ceil(nx - shift_x),
                  np.ceil(ny - shift_y),
                  nz]
    elif shift_units == 'per':
        # center of ROI [x, y, z]
        center = [np.ceil(nx * shift_x / 100),
                  np.ceil(ny * shift_y / 100),
                  nz]

    # initialize 3d grid
    xx, yy, zz = np.mgrid[:nx, :ny, :nz]
    mask3d_final = np.zeros((nx, ny, nz))

    # initialize mask list
    mask3d_list = list()
    pix_intensity_list = list()
    median_dict = dict()
    stat_param = dict()

    xc = center[0]
    yc = center[1]
    zc = center[2]

    radius_x = np.ceil((int(size) - 1) / 2)  # radius of square - 5 to left, 5 to right
    radius_y = radius_x  # square --> radius left-right and top-bottom are same
    radius_z = nz

    # TODO - replace rotation with something else, so it'll be possible
    # TODO - to create ROIs on images of different number of rows and columns
    for it in range(4):
        mask3d = ((abs(xx - xc) <= radius_x) &
                  (abs(yy - yc) <= radius_y) &
                  (abs(zz - zc) <= radius_z))
        mask3d = mask3d.astype('int')
        mask3d = np.rot90(mask3d, k=it)

        # Append current ROI to list
        mask3d_list.append(mask3d)

    pix_intensity_list = list()
    # Loop across ROIs
    for it in range(4):
        pix_intensity_list.append(np.where(mask3d_list[it] == 1, file_data, np.NaN))
        # Loop across slices
        for slice in range(nz):
            stat_param[it, slice] = np.nanmax(pix_intensity_list[it][:, :, slice]), \
                                    (np.nanquantile(pix_intensity_list[it][:, :, slice], 0.75) -
                                     np.nanquantile(pix_intensity_list[it][:, :, slice], 0.25))
                                    # [it, slice] - [ROI, slice]
                                    # maxVal, minVal, std, median, IQR

    if args.filter != 'None':
        # Identify outlier values
        for it in range(4):
            for slice in range(nz):
                if it == 0:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 1, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 2, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 3, slice][filter_method]):
                       pix_intensity_list[it][:,:,slice] = np.NaN
                       mask3d_list[it][:,:,slice] = 0
                if it == 1:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 1, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 1, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 2, slice][filter_method]):
                       pix_intensity_list[it][:,:,slice] = np.NaN
                       mask3d_list[it][:, :, slice] = 0
                if it == 2:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 2, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 1, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 1, slice][filter_method]):
                       pix_intensity_list[it][:,:,slice] = np.NaN
                       mask3d_list[it][:, :, slice] = 0
                if it == 3:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 3, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 2, slice][filter_method]) or \
                       (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 1, slice][filter_method]):
                        pix_intensity_list[it][:, :, slice] = np.NaN
                        mask3d_list[it][:, :, slice] = 0

    mask3d_final = reduce(operator.add, mask3d_list)          # Same as mask3d_list[0] + mask3d_list[1] + mask3d_list[2] + mask3d_list[3]

    # Possible visualisation, takes a lot of time
    if visualise == 1:
        print('Creating, please wait...')
        ax = plt.figure().add_subplot(projection='3d')
        ax.voxels(mask3d_final, facecolors='red')

        # Save figure as png file
        # Path to save isn't specified:
        if len(path_output) == 1:
            path = os.path.dirname(os.path.normpath(path_input_file))
            name = os.path.basename(path_input_file.strip('.nii.gz')) + '_noise_mask.png'
            plt.savefig(os.path.join(path, name))

        # Path to save is specified by user:
        else:
            name = os.path.basename(path_input_file.strip('.nii.gz')) + '_noise_mask.png'
            plt.savefig(os.path.join(path_output, name))

    # save selected ROI to NIfTI
    # We have to check, whether output path is
    if len(path_output) == 1:
        output_filename = os.path.normpath(path_input_file.strip('.nii.gz') + '_noise_mask.nii.gz')
    else:
        name = os.path.basename(path_input_file.strip('.nii.gz')) + '_noise_mask.nii.gz'
        output_filename = os.path.join(path_output, name)

    mask3d_final = mask3d_final.astype('uint8')
    mask3d_final = nib.Nifti1Image(mask3d_final, affine=test_img1_affine, header=test_img1_hdr)
    nib.save(mask3d_final, output_filename)
    print('Created {}'.format(output_filename))


if __name__ == "__main__":
    main()
