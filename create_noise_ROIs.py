#!/usr/bin/python3

import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

testImg1 = nib.load('sub-5147B_acq-02_inv-01_part-mag_IRT1.nii.gz')
# testImg1 = testImg1.get_fdata()[:, :, 2]
positionShift = 5                                   # value in %, should cover noise in background
slice = 2
color = 'Blues'

testImg1Hdr = testImg1.header
nx, ny, nz = testImg1Hdr.get_data_shape()           # data shape, number of pixels in x, y, z

size = 10                                           # size of the square
center = [np.ceil(nx*positionShift/100),
          np.ceil(ny*positionShift/100), nz]        # center of ROI [x, y, z]

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











