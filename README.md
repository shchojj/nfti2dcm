# nfti2dcm
nifti-formatted image/labels are converted to dicom/rtstruct
In fact, there are many scripts of dicom to Nifti on the Internet, but when I needed to convert NIFti to DICom form, I found that I could not find the right one for a long time. Especially when I need to import some public data set into a supported program, especially if the program input may only support DICOM form, so have this program。
Image is converted into dicom code reference https://simpleitk.readthedocs.io/en/v1.1.0/Examples/DicomSeriesReadModifyWrite/Documentation.html
The label into dicom code reference https://github.com/Abe404/mask_to_dicom_rtstruct
Thanks to these authors.
