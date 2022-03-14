from __future__ import print_function

import SimpleITK as sitk

import sys, time, os
import numpy as np
from pydicom.uid import generate_uid

def convert(patienID: str,patienName: str, input_nifti_path: str,  output_dicom_path: str):
    img = sitk.ReadImage(input_nifti_path)
    # vol = sitk.GetArrayFromImage(img)
    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")
    # print('modification_date',modification_date)
    # sop_instance_uid = generate_uid()
    # print('sop_instance_uid',sop_instance_uid)
    FrameOfReferenceUID='1.2.826.0.1.3680043.2.1125.'+ modification_date + ".2"
    direction = img.GetDirection()
    series_tag_values = [("0008|0005", "ISO_IR 192"),
                         ("0008|0031", modification_time),  # Series Time
                         ("0008|0021", modification_date),  # Series Date
                         ("0008|0008", "DERIVED\\SECONDARY"),  # Image Type
                         ("0020|000d", "1.2.826.0.1.3680043.2.1125." + modification_date + ".1"),
                         ("0020|000e", "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time),
                         ("0020|0052", "1.2.826.0.1.3680043.2.1125." + modification_date + ".2"),
                         ("0020|0011", "1"),
                         # Series Instance UID
                         ("0020|0037",
                          '\\'.join(map(str, (direction[0], direction[3], direction[6],  # Image Orientation (Patient)
                                              direction[1], direction[4], direction[7])))),
                         ("0008|103e", "datu")
                         ]
    list(map(lambda i: write_slices(patienID, patienName,output_dicom_path, series_tag_values, img, i), range(img.GetDepth())))
    # list(map(lambda i: write_slices(output_dicom_path, series_tag_values, img, i), range(0,1)))
def write_slices(patienID, patienName, data_directory, series_tag_values, new_img, i):
    spacing = new_img.GetSpacing()
    origin = new_img.GetOrigin()
    PixelSpacing = [spacing[0], spacing[1]]
    SliceThickness = spacing[2]
    SliceLocation = origin[2]+i*SliceThickness
    ImagePositionPatient=[origin[0], origin[1],SliceLocation]

    image_slice = new_img[:, :, i]

    # Tags shared by the series.
    list(map(lambda tag_value: image_slice.SetMetaData(tag_value[0], tag_value[1]), series_tag_values))

    # Slice specific tags.
    image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))  # Instance Creation Date
    image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))  # Instance Creation Time

    # Setting the type to CT preserves the slice location.
    image_slice.SetMetaData("0008|0060", "CT")  # set the type to CT so the thickness is carried over
    image_slice.SetMetaData("0018|0050", str(SliceThickness))  # Instance Number
    image_slice.SetMetaData("0018|5100",'HFS')
    image_slice.SetMetaData("0010|0010",patienName)
    image_slice.SetMetaData("0010|0020",patienID)

    # (0020, 0032) image position patient determines the 3D spacing between slices.
    image_slice.SetMetaData("0020|0032", '\\'.join( map(str, new_img.TransformIndexToPhysicalPoint((0, 0, i)))))  # Image Position (Patient)
    image_slice.SetMetaData("0020|0013", str(i))  # Instance Number
    image_slice.SetMetaData("0020|1041", str(new_img.TransformIndexToPhysicalPoint((0, 0, i))[2]))  # Instance Number

    # Write to the output directory and add the extension dcm, to force writing in DICOM format.
    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    writer.SetFileName(os.path.join(data_directory, str(i) + '.dcm'))
    writer.Execute(image_slice)

if __name__ == "__main__":
    input_nifti_path = r'../data/nifti/img.nii.gz'
    output_dicom_path = r'../data/out/q'
    patienName = 'GengYuYun'
    patienID = '04539'

    convert(patienID, patienName, input_nifti_path, output_dicom_path)