#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script used for performing the morphometry analysis
    during the BME Labs

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: March 2023
    """

import os
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import SimpleITK as sitk
from skimage import measure
from Utils import RotationMatrix, Resample, Read, Morphometry


def segment_bone(image, sigma=0.02, threshold=None, nThresholds=2, mask=True, closeSize=None):

    """
    Perform segmentation of bone form gray value image
    Step 1: gaussian filter to smooth values
    Step 2: multiple otsu's algorithm for segmentation
    Step 3: Crop image at bone limits to reduce memory print
    Step 4: pad to avoid border contacts
    Optional
    Step 5: close contour for further hole filling
    Step 6: Mask generation by filling holes in consecutive z slices

    :param image: Gray value image
                    - Type: sitkImage
    :param sigma: Filter width
                    - Type: float
    :param nThresholds: Number of otsu's threshold
                    - Type: int
    :param mask: Generate mask or not
                    - Type: bool
    :param closeSize: Radius used to close the contour
                    - Type: int

    :return grayCrop: Cropped gray value image
            binCrop: Cropped binary segmented image
            Mask: Generated mask
    """

    # Filter image to reduce noise
    gauss = sitk.SmoothingRecursiveGaussianImageFilter()
    gauss.SetSigma(sigma)
    smooth  = gauss.Execute(image)

    if threshold:
        # Segment image using single threshold
        binarize = sitk.BinaryThresholdImageFilter()
        binarize.SetUpperThreshold(threshold)
        binarize.SetOutsideValue(255)
        binarize.SetInsideValue(0)
        bin = binarize.Execute(smooth)

    else:
        # Segment image by thresholding using otsu's method
        otsu = sitk.OtsuMultipleThresholdsImageFilter()
        otsu.SetNumberOfThresholds(nThresholds)
        seg = otsu.Execute(smooth)

        # binarize image to keep bone only
        binarize = sitk.BinaryThresholdImageFilter()
        binarize.SetUpperThreshold(nThresholds - 1)
        binarize.SetOutsideValue(255)
        binarize.SetInsideValue(0)
        bin = binarize.Execute(seg)

    # Remove unconnected components
    labels = sitk.ConnectedComponent(bin)
    sorted = sitk.RelabelComponent(labels, sortByObjectSize=True)
    bin = sorted == 1

    # Crop image to bone
    array = sitk.GetArrayFromImage(bin)
    z, y, x = np.where(array > 0)
    x1, x2 = int(x.min()), int(x.max())
    y1, y2 = int(y.min()), int(y.max())
    z1, z2 = int(z.min()), int(z.max())
    binCrop = sitk.Slice(bin, (x1, y1, z1), (x2, y2, z2))
    grayCrop = sitk.Slice(image, (x1, y1, z1), (x2, y2, z2))

    # pad images to avoid contact with border
    binCrop = sitk.ConstantPad(binCrop, (1, 1, 1), (1, 1, 1))
    grayCrop = sitk.ConstantPad(grayCrop, (1, 1, 1), (1, 1, 1))

    if mask:
        # close contour
        close = sitk.BinaryMorphologicalClosingImageFilter()
        close.SetForegroundValue(255)
        close.SetKernelRadius(closeSize)
        closed = close.Execute(binCrop)

        # Generate mask slice by slice
        size = binCrop.GetSize()
        mask = binCrop
        for start in range(size[2]):
            slice = sitk.Slice(closed, (0, 0, start), (size[0], size[1], start + 1))

            # pad slice in z direction to "close" holes
            pad = sitk.ConstantPad(slice, (0, 0, 1), (0, 0, 1), 255)

            # fill holes
            fill = sitk.BinaryFillholeImageFilter()
            fill.SetForegroundValue(255)
            filled = fill.Execute(pad)

            # Get center slice
            slice = sitk.Slice(filled, (0, 0, 1), (size[0], size[1], 2))

            # Paste slice into original image
            mask = sitk.Paste(mask, slice, slice.GetSize(), destinationIndex=(0, 0, start + 1))

        return grayCrop, binCrop, mask
    
    else:
        return grayCrop, binCrop


def align_cylinder(bin):

    # Compute cylinder main axis
    array = sitk.GetArrayFromImage(bin)
    coords = np.where(array)
    table = pd.DataFrame()
    table['x'] = coords[2]
    table['y'] = coords[1]
    table['z'] = coords[0]

    xm = table.groupby('z')['x'].mean()
    ym = table.groupby('z')['y'].mean()

    x = np.matrix(xm).T
    y = np.matrix(ym).T
    z = table['z'].unique()
    z = np.matrix([np.ones(len(z)), z]).T

    xi, xc = np.linalg.inv(z.T * z) * z.T * x
    yi, yc = np.linalg.inv(z.T * z) * z.T * y

    # Rotate cylinder
    print('\nRotate cylinder')
    phi = np.arctan(float(yc))
    print('phi : %.3f °' % (phi / np.pi * 180))
    theta = np.arctan(float(xc))
    print('theta : %.3f °' % (theta / np.pi * 180))
    r = RotationMatrix(Phi=phi, Theta=theta)

    padSize = 10
    pad = sitk.ConstantPad(bin, (padSize, padSize, padSize), (padSize, padSize, padSize))
    pad.SetOrigin((0, 0, 0))
    center = np.array(pad.GetSize()) * np.array(pad.GetSpacing()) / 2

    t = sitk.VersorRigid3DTransform()
    t.SetCenter(center)
    t.SetMatrix(r.ravel())

    rotated = sitk.Resample(pad, t, sitk.sitkNearestNeighbor)

    return rotated


# Main code

def main(arguments):

    # Read Scan
    isq = Read.ISQ(arguments.sample)[0]

    # Segment scan by using Otsu's threshold
    bin = segment_bone(isq, nThresholds=1, mask=False)[1]
    
    # Align cylinder main axis with image Z axis (for uFE simulation)
    rotated = align_cylinder(bin)

    # Bone morphometry on full sample
    array = sitk.GetArrayFromImage(rotated)
    
    # Get middle slice to compute cirle
    slice = array[rotated.GetSize()[2] // 2, :, :]
    props = measure.regionprops(slice)[0]
    
    # Create corresponding disk mask
    xc = (props.bbox[0] + props.bbox[2]) / 2
    yc = (props.bbox[1] + props.bbox[3]) / 2
    xr = (props.bbox[2] - props.bbox[0]) / 2
    yr = (props.bbox[3] - props.bbox[1]) / 2
    area = np.pi * ((xr + yr) / 2 * rotated.GetSpacing()[0]) ** 2
    x, y = np.ogrid[-int(xc):rotated.GetSize()[0] - int(xc),
           -int(yc):rotated.GetSize()[1] - int(yc)]
    mask = x ** 2 + y ** 2 <= ((xr + yr) / 2) ** 2
    
    # Create corresponding cylinder
    cylinder = np.repeat(mask, rotated.GetSize()[2]) * 255
    cylinder = np.reshape(cylinder, rotated.GetSize()[::-1], order='F').astype('uint')
    
    # Crop external cylinder boundaries to sample limits
    cylinder[array.sum(axis=(1, 2)) == 0] = 0
    
    # Compute height and BV/TV of the sample
    coords = np.argwhere(cylinder)
    height = (coords[-1,0] - coords[0,0]) * rotated.GetSpacing()[2]
    bvtv = array.sum() / cylinder.sum()

    return


if __name__ == '__main__':

    # Initiate the parser with a description
    fc = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=Description, formatter_class=fc)

    # Add long and short argument
    sv = parser.prog + ' version ' + Version
    parser.add_argument('-v', '--version', help='Show script version', action='version', version=sv)
    parser.add_argument('sample', help='Sample number', type=str)

    # Read arguments from the command line
    arguments = parser.parse_args()

    main(arguments)
