import time,os
from src.img2dcm import convert as img_convert
from src.msk2dcm import convert as msk_convert
def nifti2dcm(input_nifti_dir, output_dicom_dir, struct_name, roi_list):
    head, tail = os.path.split(input_nifti_dir)
    patienName = tail
    patienID = time.strftime("%Y%m%d%H%M%S")+'_'+patienName
    img_path = os.path.join(input_nifti_dir,'img.nii.gz')
    msk_path = os.path.join(input_nifti_dir,'label.nii.gz')
    img_convert(patienID, patienName, img_path, output_dicom_dir)
    print('Image conversion completed. ')
    msk_convert(msk_path, output_dicom_dir, output_dicom_dir, struct_name, roi_list)
    print('Label conversion completed. ')

if __name__ == "__main__":
    input_nifti_dir = r'F:\workspace\python\nii2dcm\data\nifti'
    output_dicom_dir = r'C:\Users\shcho\Desktop\1'
    struct_name = 'datu'
    roi_list =['Eye-Right', 'Eye-Left','Mandible', 'Optic-Nerves-Left', 'Optic-Nerves-Right',
               'Parotid-Gland-Left', 'Parotid-Gland-Right', 'Brain']
    # roi_list = ['Eye-Right', 'Eye-Left', 'Mandible', 'Optic-Nerves-Left']
    nifti2dcm(input_nifti_dir, output_dicom_dir, struct_name, roi_list)