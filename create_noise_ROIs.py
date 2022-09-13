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
from matplotlib.ticker import PercentFormatter

from argparse import RawDescriptionHelpFormatter


def get_parser():
    # TODO - make -i flag mandatory argument, now it appears among optional args when you print help
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
    parser.add_argument('-visualise', help='Visualisation of created ROI for debug. 0 - do not visualise, 1 - visualise',
                        type=int, required=False, default=0, choices=[0, 1])
    parser.add_argument('-outpath', help='Path for saving nii file with created ROIs. Example: /home/<usr>/<directory>'
                        , type=str, required=False, default='')
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
    visualise = args.visualise
    path_to_save = args.outpath

    test_img1_hdr = testImg1.header
    test_img1_affine = testImg1.affine
    slice = 3

    # Converting to 2D for testing 2D ROI selection
    # testImg1 = testImg1.get_fdata()[:, :, slice]

    color = 'Blues'                              # setting color for visualisation

    nx, ny, nz = test_img1_hdr.get_data_shape()  # data shape, number of pixels in x, y, z
    # convert NiftiImage class to ndarray
    testImg1 = testImg1.get_fdata()

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

    # mask for 3D volume
    if testImg1.ndim > 2:

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

        for it in range(4):
            mask3d = ((abs(xx - xc) <= radius_x) &
                      (abs(yy - yc) <= radius_y) &
                      (abs(zz - zc) <= radius_z))
            mask3d = mask3d.astype('int')
            mask3d = np.rot90(mask3d, k=it)
            testImg1 = np.rot90(testImg1, k=it)

            # Append current ROI to list
            mask3d_list.append(mask3d)
            mask3d_final = mask3d_final + mask3d

        pix_intensity_list = list()
        # Loop across ROIs
        for it in range(4):
            pix_intensity_list.append(np.where(mask3d_list[it] == 1, testImg1, np.NaN))
            # Loop across slices
            for slice in range(testImg1.shape[2]):
                stat_param[it, slice] = np.nanmax(pix_intensity_list[it][:, :, slice]), \
                                        np.nanmin(pix_intensity_list[it][:, :, slice]), \
                                        np.nanstd(pix_intensity_list[it][:, :, slice]), \
                                        np.nanmedian(pix_intensity_list[it][:, :, slice])
                                        # [it, slice] - [ROI, slice]
                                        # maxVal, minVal, std, median
                                        # slicewise inspection of each slice is now possible
                                        # so we can extract slice containing aliasing artifact

        # Identify outlier values
        for it in range(4):
            for slice in range(testImg1.shape[2]):
                if it == 0:
                    if (stat_param[it, slice][0] > 1.5 * stat_param[it + 1, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it + 2, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it + 3, slice][0]):
                       pix_intensity_list[it][:,:,slice] = np.NaN
                       mask3d_list[it][:,:,slice] = 0
                if it == 1:
                    if (stat_param[it, slice][0] > 1.5 * stat_param[it - 1, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it + 1, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it + 2, slice][0]):
                       pix_intensity_list[it][:,:,slice] = np.NaN
                       mask3d_list[it][:, :, slice] = 0
                if it == 2:
                    if (stat_param[it, slice][0] > 1.5 * stat_param[it - 2, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it - 1, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it + 1, slice][0]):
                       pix_intensity_list[it][:,:,slice] = np.NaN
                       mask3d_list[it][:, :, slice] = 0
                if it == 3:
                    if (stat_param[it, slice][0] > 1.5 * stat_param[it - 3, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it - 2, slice][0]) or \
                       (stat_param[it, slice][0] > 1.5 * stat_param[it - 1, slice][0]):
                        pix_intensity_list[it][:, :, slice] = np.NaN
                        mask3d_list[it][:, :, slice] = 0

        # TODO - Correctly reorient created mask so it'll correspond with image
        # created mask is somehow rotated compared to input image
        # firstly created ROI is in bottom right corner, it should be in upper left!

        #      NOW                     EXPECTED
        #  _____ _____               _____ _____
        # |     |     |             |     |     |
        # |  3  |  4  |             |  1  |  2  |
        # |_____|_____|             |_____|_____|
        #  _____ _____       -->     _____ _____
        # |     |     |             |     |     |
        # |  2  |  1  |             |  3  |  4  |
        # |_____|_____|             |_____|_____|

        # pls suggestions for correction
        # mask3d_final = mask3d_final + np.where(mask3d_list[it] == 1, testImg1, False)

        # Compare histograms of noise
        # plt.hist(pix_intensity_vec, bins=256)
        # plt.ylim([0,1000])

        # Show histograms of all ROIs
        # plt.legend(['data0', 'data1', 'data2', 'data3'])
        # plt.show()

        # Possible visualisation but for some reason takes a lot of time
        if visualise == 1:
            ax = plt.figure().add_subplot(projection='3d')
            ax.voxels(mask3d_final, facecolors='red')

            # Save figure as png file
            # Path to save isn't specified:
            if len(path_to_save) == 0:
                path = os.path.dirname(os.path.normpath(path_input_file))
                plt.savefig(os.path.join(path, 'created_3d_roi.png'))

            # Path to save is specified by user:
            else:
                plt.savefig(os.path.join(path_to_save, 'created_3d_roi.png'))

            # Visualisation in 2D plane
            # plt.figure()
            # plt.imshow(testImg1[:, :, slice].T, cmap='gray')
            # plt.imshow(mask3d_final[:, :, slice].T, cmap=color, alpha=0.3)
            # plt.savefig(os.path.join(os.getcwd(), 'created_2d_roi.png'))

        # save selected ROI to NIfTI
        if len(path_to_save) == 0:
            path = os.path.dirname(os.path.normpath(path_input_file))
            output_filename = os.path.join(path, 'noise_mask.nii')
        else:
            output_filename = os.path.join(path_to_save, 'noise_mask.nii')

        mask3d_final = mask3d_final.astype('uint8')
        mask3d_final = nib.Nifti1Image(mask3d_final, affine=test_img1_affine, header=test_img1_hdr)
        nib.save(mask3d_final, output_filename)
        print('Created {}'.format(output_filename))

    # mask for single 2D slice
    else:

        # initialize 2d grid
        xx, yy = np.mgrid[:nx, :ny]
        mask2d_final = np.zeros((nx, ny))
        pix_intensity = np.zeros((nx, ny))
        xc = center[0]
        yc = center[1]

        radius_x = np.ceil((int(size) - 1) / 2)         # radius of square - 5 to left, 5 to right
        radius_y = radius_x                             # square --> radius left-right and top-bottom are same

        for it in range(4):
            mask2d = ((abs(xx - xc) <= radius_x) &
                      (abs(yy - yc) <= radius_y)) * 1   # * 1 to convert F or T to 0 or 1
            mask2d = np.rot90(mask2d, k=it)
            mask2d_final = mask2d_final + mask2d
            testImg1 = np.rot90(testImg1, k=it)

        # Visualisation in 2D plane
        if visualise == 1:
            plt.figure()
            plt.imshow(np.rot90(testImg1, k=1), cmap='gray')
            plt.imshow(mask2d_final, cmap=color, alpha=0.3)

            if len(path_to_save) == 0:
                path = os.path.dirname(os.path.normpath(path_input_file))
                plt.savefig(os.path.join(path, 'created_2d_roi.png'))
            else:
                plt.savefig(os.path.join(path_to_save, 'created_2d_roi.png'))

        # getting pixel intensities from mask
        pix_intensity = np.where(mask2d_final == 1, testImg1, False)
        std = np.std(pix_intensity)

        # save mask to NIfTI
        # I think it makes no sense to save 2D mask to a nifti file, but may be useful...
        # mask2d_final = nib.Nifti1Image(mask2d_final, affine=test_img1_affine)
        # nib.save(mask2d_final, 'noise_mask_2d.nii')

        # appending values of mask2d to mask2d_final if
        # count of rows and columns of input image are not equal!!!

        txt = 'Standard deviation value is: {0:.5f}'
        print(txt.format(std))


if __name__ == "__main__":
    main()
