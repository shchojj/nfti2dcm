import random
import nibabel as nib
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import SimpleITK as sitk
from skimage import measure
from _collections import defaultdict
from loguru import logger
import sys

import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
import time
import datetime
import warnings

def filter_rtss(input_dicom_path, dcm_file_names):
    fn_list=[]
    for fn in dcm_file_names:
        ds = pydicom.dcmread(os.path.join(input_dicom_path, fn), force=True)
        if not hasattr(ds, 'Modality'):
            logger.error(f'object has no attribute Modality: {fn}')
            continue
        if ds.Modality != "RTSTRUCT" and ds.Rows == ds.Columns:
            fn_list.append(fn)
    return fn_list

def ImagePositionPatientOrdering(slices):
    # image position patient ordering
    first = True
    min = 0
    max = 0
    normal = np.empty(3)
    distmultimap = defaultdict(list) # Use a multimap to sort the distances from 0,0,0
    for slice in slices:
        if first:
            cosines = slice.ImageOrientationPatient

            # You only have to do this once for all slices in the volume.Next,
            # for each slice, calculate the distance along the slice normal
            # using the IPP("Image Position Patient") tag.
            # ("dist" is initialized to zero before reading the first slice):
            normal[0] = cosines[1] * cosines[5] - cosines[2] * cosines[4]
            normal[1] = cosines[2] * cosines[3] - cosines[0] * cosines[5]
            normal[2] = cosines[0] * cosines[4] - cosines[1] * cosines[3]

            ipp = slice.ImagePositionPatient

            dist = 0
            for i in range(3):
                dist += normal[i] * ipp[i]

            distmultimap[dist].append(slice)

            max = min = dist
            first = False
        else:
            ipp = slice.ImagePositionPatient

            dist = 0
            for i in range(3):
                dist += normal[i] * ipp[i]

            distmultimap[dist].append(slice)

            if min > dist: min = dist
            if max < dist: max = dist

    # Find out if min/max are coherent
    if min == max:
        logger.error("\nLooks like all images have the exact same image position. No PositionPatientOrdering sort performed\n")
        sys.exit(1)

    # Check to see if image shares a common position
    skip_file_count = 0
    repeated_dist = []
    for dist in distmultimap:
        if len(distmultimap[dist]) != 1:
            skip_file_count += len(distmultimap[dist])-1
            repeated_dist.append(dist)
    if skip_file_count > 0 :
        logger.warning(f"\nDistance: {repeated_dist} position is not unique, and will skip {skip_file_count} files totally\n")

    keys = list(distmultimap.keys())
    keys.sort()
    return [distmultimap[key][0] for key in keys]


def concatenate_coordinates(coordinates_x, coordinates_y, coordinates_z):
    vector = np.zeros((len(coordinates_x) * 3, 1))

    for i in range(len(coordinates_x)):
        vector[i * 3 + 0] = coordinates_x[i]
        vector[i * 3 + 1] = coordinates_y[i]
        vector[i * 3 + 2] = coordinates_z[i]

    return vector


def copy_attributes(rt_ds, im_ds):
    """ copy the image attributes to the struct dataset
        many of these should be the same.
    """
    for attribute in [
        'AccessionNumber',
        'FrameOfReferenceUID',
        'InstitutionAddress',
        'InstitutionName',
        'InstitutionalDepartmentName',
        'Manufacturer',
        'ManufacturerModelName',
        'PatientBirthDate',
        'PatientID',
        'PatientName',
        'PatientSex',
        'PositionReferenceIndicator',
        'ReferringPhysicianName',
        'SeriesDate',
        'SeriesNumber',
        'SeriesTime',
        'SeriesInstanceUID',
        # consider modifying later as structure may be created by different software
        'SoftwareVersions',
        'SpecificCharacterSet',
        'StudyDate',
        'StudyDescription',
        'StudyID',
        'StudyInstanceUID',
        'StudyTime']:
        if hasattr(im_ds, attribute):
            setattr(rt_ds, attribute, getattr(im_ds, attribute))

    return rt_ds


def add_referenced_frame_of_reference(rt_ds, im_ds):
    """ tell the rtstruct to reference the image series """
    # Referenced Frame of Reference Sequence (link structure to sequence)
    frame_of_ref_sequence = Sequence()
    rt_ds.ReferencedFrameOfReferenceSequence = frame_of_ref_sequence
    frame_of_ref = Dataset()
    frame_of_ref.FrameOfReferenceUID = im_ds[0].FrameOfReferenceUID
    referenced_study_sequence = Sequence()
    frame_of_ref.RTReferencedStudySequence = referenced_study_sequence
    referenced_study = Dataset()
    referenced_study.ReferencedSOPClassUID = '1.2.840.10008.3.1.2.3.2'

    # I'm not sure what this should be.
    # I think rt_referenced series instances can be used to store
    # the UIDs of all the corresponding serires files. I'm hoping it
    # isn't essential as it seems useless and overkill to add this here.
    # The series is already referenced by it's SeriesInstanceUID
    referenced_study.ReferencedSOPInstanceUID = im_ds[0].SOPInstanceUID

    referenced_series_sequence = Sequence()
    referenced_study.RTReferencedSeriesSequence = referenced_series_sequence
    referenced_series = Dataset()
    # referenced_series.SeriesInstanceUID = im_ds[0].SeriesInstanceUID

    contour_image_sequence = Sequence()
    referenced_series.ContourImageSequence = contour_image_sequence
    for cidx in range(len(im_ds)):
        contour_im = Dataset()
        contour_im.ReferencedSOPClassUID = im_ds[cidx].SOPClassUID
        contour_im.ReferencedSOPInstanceUID = im_ds[cidx].SOPInstanceUID
        contour_image_sequence.append(contour_im)



    referenced_series_sequence.append(referenced_series)
    referenced_study_sequence.append(referenced_study)
    frame_of_ref_sequence.append(frame_of_ref)


def create_rt_struct_ds(im_ds, filename, label='label', name='name'):
    """ create rt struct dicom dataset """
    SOPClassUIDForRTSTRUCT = '1.2.840.10008.5.1.4.1.1.481.3'
    # SOPInstanceUID is a globally unique identifier for a DICOM file
    sop_instance_uid = generate_uid()
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = SOPClassUIDForRTSTRUCT
    file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    file_meta.MediaStorageSOPInstanceUID = sop_instance_uid

    rt_ds = FileDataset(filename, {}, file_meta=file_meta, preamble=b"\0" * 128)
    # copy patient related attributes etc from an existing image dicom file.
    rt_ds = copy_attributes(rt_ds, im_ds[0])
    rt_ds.SOPInstanceUID = sop_instance_uid
    rt_ds.InstanceCreationTime = str(time.time())

    # Time at which the structures were last modified.
    rt_ds.StructureSetTime = rt_ds.InstanceCreationTime
    # YYYYMMDD
    rt_ds.InstanceCreationDate = datetime.datetime.now().strftime('%Y%m%d')
    rt_ds.StructureSetDate = rt_ds.InstanceCreationDate

    rt_ds.StructureSetLabel = label
    rt_ds.StructureSetName = name
    rt_ds.Modality = 'RTSTRUCT'
    rt_ds.SOPClassUID = SOPClassUIDForRTSTRUCT
    rt_ds.is_implicit_VR = False
    rt_ds.is_little_endian = True

    referenced_study_sequence = Sequence()
    rt_ds.ReferencedStudySequence = referenced_study_sequence
    referenced_study = Dataset()
    referenced_study.ReferencedSOPClassUID = '1.2.840.10008.3.1.2.3.2'
    referenced_study.ReferencedSOPInstanceUID = im_ds[0].SOPInstanceUID
    referenced_study_sequence.append(referenced_study)

    add_referenced_frame_of_reference(rt_ds, im_ds)

    rt_ds.ROIContourSequence = Sequence()
    rt_ds.RTROIObservationsSequence = Sequence()
    rt_ds.StructureSetROISequence = Sequence()
    return rt_ds


def add_contour(rt_ds, contours, name, slice_ds, color=[255, 0, 0]):
    """
    rt_ds is an existing rt struct dataset
    contour is a list of lists
    each list in contours is a list of contour_points
    with the format [x,y,z,x,y,z....]
    """
    # get the max contour number found so far
    # largest ROI number may be greater than the length of the sequence (and some numbers might be missing)
    contour_num = 1
    if len(rt_ds.StructureSetROISequence):
        max_contour_num = max([s.ROINumber for s in rt_ds.StructureSetROISequence])
        contour_num = max_contour_num + 1

    roi_contour = Dataset()
    roi_contour.ROIDisplayColor = color
    contour_sequence = Sequence()
    roi_contour.ContourSequence = contour_sequence

    # each list of contour points_mm represent a closed polygon
    # in a single slice
    # contours do not exist accross multiple axial slices
    # maybe they should?
    # print('contours',len(contours))
    # for points_mm in contours:
    for cidx in range(len(contours)):
        points_mm = contours[cidx]
        for contour in points_mm:
            # if len(contour)>0:
            contour_ds = Dataset()
            contour_ds.ContourGeometricType = 'CLOSED_PLANAR'
            contour_ds.NumberOfContourPoints = str(contour.shape[0] // 3)
            # print('contour',contour.shape)
            contour_ds.ContourData = contour.tolist()

            # Sequence of images containing the contour.
            contour_ds.ContourImageSequence = Sequence()
            contour_im = Dataset()
            contour_im.ReferencedSOPClassUID = slice_ds[cidx].SOPClassUID
            contour_im.ReferencedSOPInstanceUID = slice_ds[cidx].SOPInstanceUID
            # contour_im.ReferencedSOPInstanceUID = rt_ds[cidx].SOPInstanceUID
            contour_ds.ContourImageSequence.append(contour_im)
            contour_sequence.append(contour_ds)

    roi_contour.ReferencedROINumber = contour_num
    rt_ds.ROIContourSequence.append(roi_contour)

    roi_observations = Dataset()
    roi_observations.ObservationNumber = contour_num
    roi_observations.ReferencedROINumber = contour_num
    roi_observations.ROIObservationLabel = name
    roi_observations.RTROIInterpretedType = ''
    roi_observations.ROIInterpreter = ''

    rt_ds.RTROIObservationsSequence.append(roi_observations)
    structure_set_roi = Dataset()
    structure_set_roi.ROIName = name
    structure_set_roi.ROINumber = contour_num
    structure_set_roi.ROIGenerationAlgorithm = ''
    structure_set_roi.ReferencedFrameOfReferenceUID = rt_ds.FrameOfReferenceUID
    rt_ds.StructureSetROISequence.append(structure_set_roi)


def mask_to_contours(mask, spacing, origin):
    label_arr = np.unique(mask)
    ridx_list = label_arr[np.nonzero(label_arr)].tolist()
    # print('ridx_list', ridx_list)
    AllCoordinates = []
    for ridx in ridx_list:
        # print(ridx)
        AllCoordinatesThisRoi = []
        rvolume = mask.copy()
        rvolume[mask != ridx] = 0
        # Loop over slices in volume, get contours for each slice
        for slice in range(rvolume.shape[2]):
            AllCoordinatesThisSlice = []
            image = rvolume[:, :, slice]
            # Get contours in this slice using scikit-image
            contours = measure.find_contours(image, 0.5)
            # Save contours for later use
            for n, contour in enumerate(contours):
                # print("n is ",n,"for slice ",slice)
                nCoordinates = len(contour[:, 0])
                # print("number of coordinates is ",len(contour[:,0])*3," for contour ",n," for slice ",slice)
                zcoordinates = slice * np.ones((nCoordinates, 1))
                # Add patient position offset
                reg_contour = np.append(contour, zcoordinates, -1)
                # Assume no other orientations for simplicity
                reg_contour[:, 0] = reg_contour[:, 0] * spacing[0] + origin[0]
                reg_contour[:, 1] = reg_contour[:, 1] * spacing[1] + origin[1]
                reg_contour[:, 2] = reg_contour[:, 2] * spacing[2] + origin[2]
                # Storing coordinates as mm instead of as voxels
                # coordinates = concatenate_coordinates(contour[:,0] * xPixelSize, contour[:,1] * yPixelSize, zcoordinates * zPixelSize)
                coordinates = concatenate_coordinates(*reg_contour.T)
                coordinates = np.squeeze(coordinates)
                AllCoordinatesThisSlice.append(coordinates)
            AllCoordinatesThisRoi.append(AllCoordinatesThisSlice)
        # print('AllCoordinatesThisRoi',len(AllCoordinatesThisRoi),rvolume.shape[2])
        AllCoordinates.append(AllCoordinatesThisRoi)
    return AllCoordinates


def convert(input_nifti_path: str, input_dicom_path: str, output_dicom_path: str,struct_name: str='SS_1',roi_list: list =None):
    #####read mask
    nii = nib.load(input_nifti_path)
    mask = nii.get_fdata()
    #####read dicom
    dicomFiles = sorted(os.listdir(input_dicom_path))
    dicomFiles = filter_rtss(input_dicom_path, dicomFiles)
    # print(len(dicomFiles))
    image_series_files = []
    for filename in dicomFiles:
        data = pydicom.read_file(os.path.join(input_dicom_path, "%s" % filename) , force=True)
        data.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian  # or whatever is the correct transfer syntax for the file
        image_series_files.append(data)
    image_series_files = ImagePositionPatientOrdering(image_series_files)
    #### set dicom
    RTDCM_name = os.path.join(output_dicom_path, "RTSTRUCT"+image_series_files[0].file_meta.MediaStorageSOPInstanceUID+".dcm")
    new_struct_ds = create_rt_struct_ds(image_series_files, RTDCM_name, label=struct_name, name=struct_name)
    ##### get geometry
    xPixelSize = image_series_files[0].PixelSpacing[0]
    yPixelSize = image_series_files[0].PixelSpacing[1]
    zPixelSize = image_series_files[0].SliceThickness
    patientPosition = image_series_files[0].ImagePositionPatient
    patientStartingZ = image_series_files[0].ImagePositionPatient[2]
    #####
    spacing = [xPixelSize, yPixelSize, zPixelSize]
    origin = [patientPosition[0], patientPosition[1], patientStartingZ]
    AllCoordinates = mask_to_contours(mask, spacing, origin)
    for ridx in range(len(AllCoordinates)):
        if roi_list is not None and len(AllCoordinates) == len(roi_list):
            roi_name = roi_list[ridx]
        else:
            print('roi_list',roi_list,len(AllCoordinates))
            warnings.warn('The list of roi names is incorrect, so we use the default names such as ROI_0, ROI_1...')
            roi_name = 'ROI_'+str(ridx)
        roi = AllCoordinates[ridx]
        color = list(np.random.choice(range(256), size=3))
        add_contour(new_struct_ds, roi, roi_name, image_series_files,color)
    new_struct_ds.SeriesInstanceUID = generate_uid()
    new_struct_ds.save_as(RTDCM_name, write_like_original=False)


if __name__ == "__main__":
    # dicom_series_path = r'../data/dcm/04539_GengYuYun'
    # output_dicom_path = r'../data/result/q'

    dicom_series_path = r'../data/out/q'
    input_nifti_path = r'../data/nifti/label.nii.gz'
    output_dicom_path = r'../data/out/q'

    convert(input_nifti_path,dicom_series_path, output_dicom_path)
