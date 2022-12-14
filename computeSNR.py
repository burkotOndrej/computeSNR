#!/usr/bin/python3
#
# Script for SNR calculation.
# Create 4 cubic ROIs in the input image corners.
#
# Example:
#       python3 computeSNR.py -i <input_image> -m <roi_mask> -shiftx 10 -shifty 10 -shiftunits pix -size 15
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

filtration_dict = {
    'max': 0,
    'iqr': 1
}


def get_parser():
    """ Define input parameters for script. """

    parser = argparse.ArgumentParser(description='Script for calculation of SNR from input image'
                                                 '\nCreate 4 cubic ROIs in the input image corners.'
                                                 '\nThe ROIs can be used for computation of background noise SD.',
                                     formatter_class=RawDescriptionHelpFormatter,
                                     add_help=False)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-i', help='Input image', type=str, required=True)
    required.add_argument('-m', help='Binary mask of ROI for SNR calculation', required=True)
    optional.add_argument('-shiftx', help='Shift of the ROIs in x axis from the image borders.', type=int,
                          required=False, default=10)
    optional.add_argument('-shifty', help='Shift of the ROIs in y axis from the image borders.', type=int,
                          required=False, default=10)
    optional.add_argument('-shiftunits', help='Shift units in pixels (pix) or percentage (per)', type=str,
                          required=False, default='pix', choices=['pix', 'per'])
    optional.add_argument('-size', help='Size of the ROIs in pixels. Example: 10', type=int, required=False, default=10)
    optional.add_argument('-filter', help='Filter method', default='max', choices=['max', 'iqr', 'None'])
    optional.add_argument('-visualise',
                          help='Visualisation of created ROI for debug. 0 - do not visualise, 1 - visualise',
                          type=int, required=False, default=0, choices=[0, 1])
    optional.add_argument('-outpath',
                          help='Path for saving nii file with created ROIs. Example: /home/<usr>/<directory>'
                          , type=str, required=False, default='')
    optional.add_argument('-h', '--help', action='help',
                          help='Show this message and exit.')
    return parser


# Define function for loading input data
def load_image(image_in_path):
    """ Load input image.

    Parameters
    ----------
    image_in_path : string
        Path to data you want to load.

    Returns
    -------
    image_in : 2D/3D NIfTI image
        Loaded input image in NIfTI format.
    """

    # Check if input file exists
    # TODO - consider to replace if-else condition by try-except
    if os.path.exists(image_in_path):
        image_in = nib.load(image_in_path)
    else:
        print('ERROR: Input file {} does not exist.'.format(image_in_path))
        sys.exit()

    return image_in


# Define function for creation of 4 cubic noise ROIs in background
def create_noise_mask(image):
    """Create 4 cubic ROIs in image corners.

    Parameters
    ----------
    image : NIfTI image
        NIfTI image where noise ROIs will be created.

    Returns
    -------
    mask3d_final : ndarray
        Array with created noise ROIs.
    """

    parser = get_parser()
    args = parser.parse_args()

    size = args.size  # size of the square
    shift_x = args.shiftx
    shift_y = args.shifty
    shift_units = args.shiftunits
    if args.filter != 'None':
        filter_method = filtration_dict[args.filter]

    # Get input image size
    im_size = image.shape  # data shape, number of pixels in x, y, z
    # 2D - if the input image is just single slice (i.e., 2D), set the third dimension to 1
    if len(im_size) == 2:
        nz = 1

    # 3D
    else:
        nz = im_size[2]
    nx, ny = im_size[0], im_size[1]  # data shape, number of pixels in x, y, z

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
    stat_param = dict()

    xc = center[0]
    yc = center[1]
    zc = center[2]

    radius_x = np.ceil((int(size) - 1) / 2)  # radius of square - 5 to left, 5 to right
    radius_y = radius_x  # square --> radius left-right and top-bottom are same
    radius_z = nz

    for it in range(4):
        mask3d = ((abs(xx - xc) <= radius_x) &
                  (abs(yy - yc) <= radius_y) &
                  (abs(zz - zc) <= radius_z))
        mask3d = mask3d.astype('int')
        mask3d = np.rot90(mask3d, k=it)
        # testImg1 = np.rot90(testImg1, k=it)

        # Append current ROI to list
        mask3d_list.append(mask3d)
        mask3d_final = mask3d_final + mask3d

    pix_intensity_list = list()

    if args.filter != 'None':
        # Identify outlier values
        for it in range(4):
            for slice in range(nz):
                if it == 0:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 1, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 2, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 3, slice][filter_method]):
                        pix_intensity_list[it][:, :, slice] = np.NaN
                        mask3d_list[it][:, :, slice] = 0
                if it == 1:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 1, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 1, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 2, slice][filter_method]):
                        pix_intensity_list[it][:, :, slice] = np.NaN
                        mask3d_list[it][:, :, slice] = 0
                if it == 2:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 2, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 1, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it + 1, slice][filter_method]):
                        pix_intensity_list[it][:, :, slice] = np.NaN
                        mask3d_list[it][:, :, slice] = 0
                if it == 3:
                    if (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 3, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 2, slice][filter_method]) or \
                            (stat_param[it, slice][filter_method] > 1.5 * stat_param[it - 1, slice][filter_method]):
                        pix_intensity_list[it][:, :, slice] = np.NaN
                        mask3d_list[it][:, :, slice] = 0

    mask3d_final = reduce(operator.add,
                          mask3d_list)  # Same as mask3d_list[0] + mask3d_list[1] + mask3d_list[2] + mask3d_list[3]

    return mask3d_final


# Define function for saving created noise ROIs
def save_image(image_in, mask_to_save):
    """ Save created noise ROIs.

    Parameters
    ----------
    image_in : NIfTI image
        Input image in NIfTI format.

    mask_to_save : ndarray
        Binary mask which has to be saved (noise ROIs in most cases).

    Returns
    -------
    This function does not return anything, only saves image to output path.
    """

    # Fetch necessary input arguments
    parser = get_parser()
    args = parser.parse_args()

    image_in_path = os.path.abspath(args.i)
    path_to_save = args.outpath

    image_in_affine = image_in.affine
    image_in_hdr = image_in.header

    # save selected ROI to NIfTI
    if len(path_to_save) == 0:
        output_filename = os.path.normpath(image_in_path.strip('.nii.gz') + '_noise_mask.nii.gz')
    else:
        output_filename = os.path.join(path_to_save, (image_in_path.strip('.nii.gz') + '_noise_mask.nii.gz'))

    mask_to_save = mask_to_save.astype('uint8')
    mask_to_save_nifti = nib.Nifti1Image(mask_to_save, affine=image_in_affine, header=image_in_hdr)
    nib.save(mask_to_save_nifti, output_filename)
    print('Created {}'.format(output_filename))


# Define function for SNR calculation
def calculate_snr(image_mask_in, noise_mask_in, input_image):
    """ Calculate SNR of input image using ROI binary mask and noise binary mask.

    Parameters
    ----------
    image_mask_in : ndarray
        Binary mask of ROI (f. e. brain).

    noise_mask_in : ndarray
        Binary mask of background noise.

    input_image : ndarray
        Pixel intensity values of whole image for SNR estimation

    Returns
    -------
    snr_final : float
        Calculated SNR.
    """

    dim = input_image.shape
    pix_intensity_roi = np.zeros(dim)
    pix_intensity_noise = np.zeros(dim)

    for val1 in range(dim[0]):
        for val2 in range(dim[1]):
            for val3 in range(dim[2]):

                # Concentrate indexing values
                vals = val1, val2, val3

                # Copy pixel intensity values for ROI
                if image_mask_in[vals] == 1:
                    pix_intensity_roi[vals] = input_image[vals]
                else:
                    pix_intensity_roi[vals] = 0

                # Copy pixel intensity values for noise
                if noise_mask_in[vals] == 1:
                    pix_intensity_noise[vals] = input_image[vals]
                else:
                    pix_intensity_noise[vals] = 0

    # Calculate standard deviation of noise (sigma_noise)
    sigma_noise = np.std(pix_intensity_noise)

    # Calculate mean value of ROI and decrease it to 2/3 of total value
    mean_roi = 0.66*(np.mean(pix_intensity_roi))

    # Calculate SNR - for formula see Totally accessible MRI, p.131 (DOI - 10.1007/978-0-387-48896-7)
    snr_final = mean_roi/sigma_noise

    return snr_final


def main():
    # Parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # Fetch input arguments
    image_in_path = os.path.abspath(args.i)
    roi_mask_path = os.path.abspath(args.m)
    visualise = args.visualise
    path_to_save = args.outpath  # TODO - check if this path exists

    # Load image
    image_in = load_image(image_in_path)
    image_in_mask = load_image(roi_mask_path)
    image_in_data = image_in.get_fdata()
    image_in_mask = image_in_mask.get_fdata()

    mask3d_final = create_noise_mask(image_in)

    # Possible visualisation, takes a lot of time
    if visualise == 1:
        print('Creating, please wait...')
        ax = plt.figure().add_subplot(projection='3d')
        ax.voxels(mask3d_final, facecolors='red')

        # Save figure as png file
        # Path to save isn't specified:
        if len(path_to_save) == 0:
            # TODO - move the following line at the beginning of the script - to check if the output path exists
            path = os.path.dirname(os.path.normpath(image_in_path))
            plt.savefig(os.path.join(path, 'created_3d_roi.png'))

        # Path to save is specified by user:
        else:
            plt.savefig(os.path.join(path_to_save, (image_in_path.strip('.nii.gz') + '_noise_mask.png')))

    # Save created binary mask
    save_image(image_in, mask3d_final)

    snr = calculate_snr(image_in_mask, mask3d_final, image_in_data)

    out_str = 'SNR of input image is: {snr:.2f}'
    print(out_str.format(snr=snr))


if __name__ == "__main__":
    main()
