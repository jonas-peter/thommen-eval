# %% #!/usr/bin/env python3
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

# %% Imports
# Modules import

import os
import struct
import argparse
import numpy as np
import sympy as sp
import pandas as pd
from pathlib import Path
import SimpleITK as sitk
from skimage import measure


# %% Functions
# Define functions
def ReadISQ(File, BMD=False, Info=False):
    """
    This function read an ISQ file from Scanco and return an ITK image and additional data.

    Adapted from https://github.com/mdoube/BoneJ/blob/master/src/org/bonej/io/ISQReader.java

    Little endian byte order (the least significant bit occupies the lowest memory position.
    00   char    check[16];              // CTDATA-HEADER_V1
    16   int     data_type;
    20   int     nr_of_bytes;
    24   int     nr_of_blocks;
    28   int     patient_index;          //p.skip(28);
    32   int     scanner_id;				//p.skip(32);
    36   int     creation_date[2];		//P.skip(36);
    44   int     dimx_p;					//p.skip(44);
    48   int     dimy_p;
    52   int     dimz_p;
    56   int     dimx_um;				//p.skip(56);
    60   int     dimy_um;
    64   int     dimz_um;
    68   int     slice_thickness_um;		//p.skip(68);
    72   int     slice_increment_um;		//p.skip(72);
    76   int     slice_1_pos_um;
    80   int     min_data_value;
    84   int     max_data_value;
    88   int     mu_scaling;             //p.skip(88);  /* p(x,y,z)/mu_scaling = value [1/cm]
    92	 int     nr_of_samples;
    96	 int     nr_of_projections;
    100  int     scandist_um;
    104  int     scanner_type;
    108  int     sampletime_us;
    112  int     index_measurement;
    116  int     site;                   //coded value
    120  int     reference_line_um;
    124  int     recon_alg;              //coded value
    128  char    name[40]; 		 		//p.skip(128);
    168  int     energy;        /* V     //p.skip(168);
    172  int     intensity;     /* uA    //p.skip(172);
    ...
    508 int     data_offset;     /* in 512-byte-blocks  //p.skip(508);
    * So the first 16 bytes are a string 'CTDATA-HEADER_V1', used to identify
    * the type of data. The 'int' are all 4-byte integers.
    *
    * dimx_p is the dimension in pixels, dimx_um the dimensions in micrometer
    *
    * So dimx_p is at byte-offset 40, then dimy_p at 44, dimz_p (=number of
    * slices) at 48.
    *
    * The microCT calculates so called 'x-ray linear attenuation' values. These
    * (float) values are scaled with 'mu_scaling' (see header, e.g. 4096) to
    * get to the signed 2-byte integers values that we save in the .isq file.
    *
    * e.g. Pixel value 8192 corresponds to lin. att. coeff. of 2.0 [1/cm]
    * (8192/4096)
    *
    * Following to the headers is the data part. It is in 2-byte short integers
    * (signed) and starts from the top-left pixel of slice 1 to the left, then
    * the next line follows, until the last pixel of the last sclice in the
    * lower right.
    """

    try:
        f = open(File, 'rb')
    except IOError:
        print("\n **ERROR**: ISQReader: intput file ' % s' not found!\n\n" % File)
        print('\n E N D E D  with ERRORS \n\n')

    for Index in np.arange(0, 200, 4):
        f.seek(Index)
        # print('Index %s :          %s' % (Index, struct.unpack('i', f.read(4))[0]))
        f.seek(Index)

    f.seek(32)
    CT_ID = struct.unpack('i', f.read(4))[0]

    if CT_ID != 6020:
        print('!!! unknown muCT -> no Slope and Intercept known !!!')

    f.seek(28)
    #    sample_nb = struct.unpack('i', f.read(4))[0]

    f.seek(108)
    Scanning_time = struct.unpack('i', f.read(4))[0] / 1000

    f.seek(168)
    Energy = struct.unpack('i', f.read(4))[0] / 1000.

    f.seek(172)
    Current = struct.unpack('i', f.read(4))[0]

    f.seek(44)
    X_pixel = struct.unpack('i', f.read(4))[0]

    f.seek(48)
    Y_pixel = struct.unpack('i', f.read(4))[0]

    f.seek(52)
    Z_pixel = struct.unpack('i', f.read(4))[0]

    f.seek(56)
    Res_General_X = struct.unpack('i', f.read(4))[0]
    # print('Resolution general X in mu: ', Res_General_X)

    f.seek(60)
    Res_General_Y = struct.unpack('i', f.read(4))[0]
    # print('Resolution general Y in mu: ', Res_General_Y)

    f.seek(64)
    Res_General_Z = struct.unpack('i', f.read(4))[0]
    # print('Resolution general Z in mu: ', Res_General_Z)

    Res_X = Res_General_X / float(X_pixel)
    Res_Y = Res_General_Y / float(Y_pixel)
    Res_Z = Res_General_Z / float(Z_pixel)

    Header_Txt = ['scanner ID:                 %s' % CT_ID,
                  'scaning time in ms:         %s' % Scanning_time,
                  'scaning time in ms:         %s' % Scanning_time,
                  'Energy in keV:              %s' % Energy,
                  'Current in muA:             %s' % Current,
                  'nb X pixel:                 %s' % X_pixel,
                  'nb Y pixel:                 %s' % Y_pixel,
                  'nb Z pixel:                 %s' % Z_pixel,
                  'resolution general X in mu: %s' % Res_General_X,
                  'resolution general Y in mu: %s' % Res_General_Y,
                  'resolution general Z in mu: %s' % Res_General_Z,
                  'pixel resolution X in mu:   %.2f' % Res_X,
                  'pixel resolution Y in mu:   %.2f' % Res_Y,
                  'pixel resolution Z in mu:   %.2f' % Res_Z]
    #    np.savetxt(inFileName.split('.')[0]+'.txt', Header_Txt)

    if Info:
        Write_File = open(File.split('.')[0] + '_info.txt', 'w')
        for Item in Header_Txt:
            Write_File.write("%s\n" % Item)
        Write_File.close()

    f.seek(44)
    Header = np.zeros(6)
    for i in range(0, 6):
        Header[i] = struct.unpack('i', f.read(4))[0]
    # print(Header)

    ElementSpacing = [Header[3] / Header[0] / 1000, Header[4] / Header[1] / 1000, Header[5] / Header[2] / 1000]
    f.seek(508)

    HeaderSize = 512 * (1 + struct.unpack('i', f.read(4))[0])
    f.seek(HeaderSize)

    VoxelModel = np.fromfile(f, dtype='i2')
    # VoxelModel = np.fromfile(f, dtype=np.float)

    NDim = [int(Header[0]), int(Header[1]), int(Header[2])]
    LDim = [float(ElementSpacing[0]), float(ElementSpacing[1]), float(ElementSpacing[2])]

    AdditionalData = {'-LDim': LDim,
                      '-NDim': NDim,
                      'ElementSpacing': LDim,
                      'DimSize': NDim,
                      'HeaderSize': HeaderSize,
                      'TransformMatrix': [1, 0, 0, 0, 1, 0, 0, 0, 1],
                      'CenterOfRotation': [0.0, 0.0, 0.0],
                      'Offset': [0.0, 0.0, 0.0],
                      'AnatomicalOrientation': 'LPS',
                      'ElementType': 'int16',
                      'ElementDataFile': File}

    # print('\nReshape data')
    # Tic = time.time()

    try:
        VoxelModel = VoxelModel.reshape((NDim[2], NDim[1], NDim[0]))
        f.close()
        del f

    except:
        # if the length does not fit the dimensions (len(VoxelModel) != NDim[2] * NDim[1] * NDim[0]),
        # add an offset with seek to reshape the image -> actualise length, delta *2 = seek

        Offset = (len(VoxelModel) - (NDim[2] * NDim[1] * NDim[0]))
        f.seek(0)
        VoxelModel = np.fromfile(f, dtype='i2')

        f.seek((len(VoxelModel) - (NDim[2] * NDim[1] * NDim[0])) * 2)
        VoxelModel = np.fromfile(f, dtype='i2')
        f.close()
        del f

        VoxelModel = VoxelModel.reshape((NDim[2], NDim[1], NDim[0]))
        # the image is flipped by the Offset --> change the order to obtain the continuous image:
        VoxelModel = np.c_[VoxelModel[:, :, -Offset:], VoxelModel[:, :, :(VoxelModel.shape[2] - Offset)]]

    if CT_ID == 6020 and BMD is True:
        # BE CAREFULL, THIS IS FOR BMD CONVERSION:
        Slope = 369.154  # ! ATTENTION, dependent on voltage, Current and time!!!
        Intercept = -191.56
        try:
            VoxelModel = VoxelModel.astype('i4')
            VoxelModel *= Slope
            VoxelModel += Intercept
        except:
            print('\n********* memory not sufficient for BMD values ************\n')

    # Convert numpy array to image
    Image = sitk.GetImageFromArray(VoxelModel)
    Image.SetSpacing(LDim[::-1])
    Image.SetOrigin([0.0, 0.0, 0.0])

    return Image, AdditionalData


def SegmentBone(Image, Sigma=0.02, Threshold=None, nThresholds=2, Mask=True, CloseSize=None):
    """
    Perform segmentation of bone form gray value image
    Step 1: gaussian filter to smooth values
    Step 2: multiple Otsu's algorithm for segmentation
    Step 3: Crop image at bone limits to reduce memory print
    Step 4: Pad to avoid border contacts
    Optional
    Step 5: Close contour for further hole filling
    Step 6: Mask generation by filling holes in consecutive z slices

    :param Image: Gray value image
                    - Type: sitkImage
    :param Sigma: Filter width
                    - Type: float
    :param nThresholds: Number of Otsu's threshold
                    - Type: int
    :param Mask: Generate mask or not
                    - Type: bool
    :param CloseSize: Radius used to close the contour
                    - Type: int

    :return GrayCrop: Cropped gray value image
            BinCrop: Cropped binary segmented image
            Mask: Generated mask
    """

    # Filter image to reduce noise
    Gauss = sitk.SmoothingRecursiveGaussianImageFilter()
    Gauss.SetSigma(Sigma)
    Smooth = Gauss.Execute(Image)

    if Threshold:
        # Segment image using single threshold
        Binarize = sitk.BinaryThresholdImageFilter()
        Binarize.SetUpperThreshold(Threshold)
        Binarize.SetOutsideValue(255)
        Binarize.SetInsideValue(0)
        Bin = Binarize.Execute(Smooth)

    else:
        # Segment image by thresholding using Otsu's method
        Otsu = sitk.OtsuMultipleThresholdsImageFilter()
        Otsu.SetNumberOfThresholds(nThresholds)
        Seg = Otsu.Execute(Smooth)

        # Binarize image to keep bone only
        Binarize = sitk.BinaryThresholdImageFilter()
        Binarize.SetUpperThreshold(nThresholds - 1)
        Binarize.SetOutsideValue(255)
        Binarize.SetInsideValue(0)
        Bin = Binarize.Execute(Seg)

    # Remove unconnected components
    Labels = sitk.ConnectedComponent(Bin)
    Sorted = sitk.RelabelComponent(Labels, sortByObjectSize=True)
    Bin = Sorted == 1

    # Crop image to bone
    Array = sitk.GetArrayFromImage(Bin)
    Z, Y, X = np.where(Array > 0)
    X1, X2 = int(X.min()), int(X.max())
    Y1, Y2 = int(Y.min()), int(Y.max())
    Z1, Z2 = int(Z.min()), int(Z.max())
    BinCrop = sitk.Slice(Bin, (X1, Y1, Z1), (X2, Y2, Z2))
    GrayCrop = sitk.Slice(Image, (X1, Y1, Z1), (X2, Y2, Z2))

    # Pad images to avoid contact with border
    BinCrop = sitk.ConstantPad(BinCrop, (1, 1, 1), (1, 1, 1))
    GrayCrop = sitk.ConstantPad(GrayCrop, (1, 1, 1), (1, 1, 1))

    if Mask:
        # Close contour
        Close = sitk.BinaryMorphologicalClosingImageFilter()
        Close.SetForegroundValue(255)
        Close.SetKernelRadius(CloseSize)
        Closed = Close.Execute(BinCrop)

        # Generate mask slice by slice
        Size = BinCrop.GetSize()
        Mask = BinCrop
        for Start in range(Size[2]):
            Slice = sitk.Slice(Closed, (0, 0, Start), (Size[0], Size[1], Start + 1))

            # Pad slice in z direction to "close" holes
            Pad = sitk.ConstantPad(Slice, (0, 0, 1), (0, 0, 1), 255)

            # Fill holes
            Fill = sitk.BinaryFillholeImageFilter()
            Fill.SetForegroundValue(255)
            Filled = Fill.Execute(Pad)

            # Get center slice
            Slice = sitk.Slice(Filled, (0, 0, 1), (Size[0], Size[1], 2))

            # Paste slice into original image
            Mask = sitk.Paste(Mask, Slice, Slice.GetSize(), destinationIndex=(0, 0, Start + 1))

        return GrayCrop, BinCrop, Mask

    else:
        return GrayCrop, BinCrop


def RotationMatrix(Phi=0.0, Theta=0.0, Psi=0.0, V=np.zeros(3), A=0):
    if (V != 0).any():
        a = np.cos(A) * np.eye(3)
        b = np.sin(A) * np.array([[0, -V[2], V[1]], [V[2], 0, -V[0]], [-V[1], V[0], 0]])
        c = (1 - np.cos(A)) * np.outer(V, V)
        R = np.round(a + b + c, 15)

    else:

        # if list of angles, use numpy for speed
        try:
            len(Phi)
            Phi, Theta, Psi = np.array(Phi), np.array(Theta), np.array(Psi)
            Rx = np.array([[np.ones(len(Phi)), np.zeros(len(Phi)), np.zeros(len(Phi))],
                           [np.zeros(len(Phi)), np.cos(Phi), -np.sin(Phi)],
                           [np.zeros(len(Phi)), np.sin(Phi), np.cos(Phi)]])

            Ry = np.array([[np.cos(Theta), np.zeros(len(Theta)), np.sin(Theta)],
                           [np.zeros(len(Theta)), np.ones(len(Theta)), np.zeros(len(Theta))],
                           [-np.sin(Theta), np.zeros(len(Theta)), np.cos(Theta)]])

            Rz = np.array([[np.cos(Psi), -np.sin(Psi), np.zeros(len(Psi))],
                           [np.sin(Psi), np.cos(Psi), np.zeros(len(Psi))],
                           [np.zeros(len(Psi)), np.zeros(len(Psi)), np.ones(len(Psi))]])

            R = np.einsum('ijl,jkl->lik', Rz, np.einsum('ijl,jkl->ikl', Ry, Rx))

        # if only float angles, use sympy for more accuracy
        except:
            Rx = sp.Matrix([[1, 0, 0],
                            [0, sp.cos(Phi), -sp.sin(Phi)],
                            [0, sp.sin(Phi), sp.cos(Phi)]])

            Ry = sp.Matrix([[sp.cos(Theta), 0, sp.sin(Theta)],
                            [0, 1, 0],
                            [-sp.sin(Theta), 0, sp.cos(Theta)]])

            Rz = sp.Matrix([[sp.cos(Psi), -sp.sin(Psi), 0],
                            [sp.sin(Psi), sp.cos(Psi), 0],
                            [0, 0, 1]])

            R = Rz * Ry * Rx

    return np.array(R, dtype='float')


def AlignCylinder(Bin):
    # Compute cylinder main axis
    Array = sitk.GetArrayFromImage(Bin)
    Coords = np.where(Array)
    Table = pd.DataFrame()
    Table['X'] = Coords[2]
    Table['Y'] = Coords[1]
    Table['Z'] = Coords[0]

    Xm = Table.groupby('Z')['X'].mean()
    Ym = Table.groupby('Z')['Y'].mean()

    X = np.matrix(Xm).T
    Y = np.matrix(Ym).T
    Z = Table['Z'].unique()
    Z = np.matrix([np.ones(len(Z)), Z]).T

    Xi, Xc = np.linalg.inv(Z.T * Z) * Z.T * X
    Yi, Yc = np.linalg.inv(Z.T * Z) * Z.T * Y

    # Rotate cylinder
    print('\nRotate cylinder')
    Phi = np.arctan(float(Yc))
    print('Phi : %.3f °' % (Phi / np.pi * 180))
    Theta = np.arctan(float(Xc))
    print('Theta : %.3f °' % (Theta / np.pi * 180))
    R = RotationMatrix(Phi=Phi, Theta=Theta)

    PadSize = 10
    Pad = sitk.ConstantPad(Bin, (PadSize, PadSize, PadSize), (PadSize, PadSize, PadSize))
    Pad.SetOrigin((0, 0, 0))
    Center = np.array(Pad.GetSize()) * np.array(Pad.GetSpacing()) / 2

    T = sitk.VersorRigid3DTransform()
    T.SetCenter(Center)
    T.SetMatrix(R.ravel())

    Rotated = sitk.Resample(Pad, T, sitk.sitkNearestNeighbor)

    return Rotated


# %% Main
# Main code

def Main(Arguments):
    # Read Scan
    ISQ = ReadISQ(Arguments.Sample)[0]

    # Segment scan by using Otsu's threshold
    Bin = SegmentBone(ISQ, Mask=False, nThresholds=1)[1]

    # Align cylinder main axis with image Z axis (for uFE simulation)
    Rotated = AlignCylinder(Bin)

    # Bone morphometry on full sample
    Array = sitk.GetArrayFromImage(Rotated)

    # Get middle slice to compute cirle
    Slice = Array[Rotated.GetSize()[2] // 2, :, :]
    Props = measure.regionprops(Slice)[0]

    # Create corresponding disk mask
    Xc = (Props.bbox[0] + Props.bbox[2]) / 2
    Yc = (Props.bbox[1] + Props.bbox[3]) / 2
    Xr = (Props.bbox[2] - Props.bbox[0]) / 2
    Yr = (Props.bbox[3] - Props.bbox[1]) / 2
    Area = np.pi * ((Xr + Yr) / 2 * Rotated.GetSpacing()[0]) ** 2
    X, Y = np.ogrid[-int(Xc):Rotated.GetSize()[0] - int(Xc),
           -int(Yc):Rotated.GetSize()[1] - int(Yc)]
    Mask = X ** 2 + Y ** 2 <= ((Xr + Yr) / 2) ** 2

    # Create corresponding cylinder
    Cylinder = np.repeat(Mask, Rotated.GetSize()[2]) * 255
    Cylinder = np.reshape(Cylinder, Rotated.GetSize()[::-1], order='F').astype('uint')

    # Crop external cylinder boundaries to sample limits
    Cylinder[Array.sum(axis=(1, 2)) == 0] = 0

    # Compute height and BV/TV of the sample
    Coords = np.argwhere(Cylinder)
    Height = (Coords[-1, 0] - Coords[0, 0]) * Rotated.GetSpacing()[2]
    BVTV = Array.sum() / Cylinder.sum()

    return


# %% Execution part
# Execution as main
if __name__ == '__main__':
    # Initiate the parser with a description
    FC = argparse.RawDescriptionHelpFormatter
    Parser = argparse.ArgumentParser(description=Description, formatter_class=FC)
    # Add long and short argument
    SV = Parser.prog + ' version ' + Version
    Parser.add_argument('-V', '--Version', help='Show script version', action='version', version=SV)
    Parser.add_argument('Sample', help='Sample number', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)
