import time,os
from src.img2dcm import convert as img_convert
from src.msk2dcm import convert as msk_convert
from src.msk2dcm import convert_list as msk_list_convert

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

def nifti_multiple_label_convert(input_nifti_dir, output_dir, struct_name):
    head, tail = os.path.split(input_nifti_dir)
    patienName = tail
    patienID = time.strftime("%Y%m%d%H%M%S") + '_' + patienName
    ids_ = os.listdir(input_nifti_dir)
    img_path = None
    msk_fn_list = []
    for i, pid in enumerate(ids_):
        # print('pid',pid)
        if 'img' in pid:
            img_path = os.path.join(input_nifti_dir,pid)
        elif 'label' not in pid:
            fn = pid[:-7]
            msk_fn_list.append(fn)

    print(img_path)
    output_dicom_dir = os.path.join(output_dir, patienName)
    os.makedirs(output_dicom_dir, exist_ok=True)
    img_convert(patienID, patienName, img_path, output_dicom_dir)
    # msk_fn_list.remove('bone')
    # msk_fn_list.remove('skin')
    # msk_fn_list.remove('artery')
    # msk_fn_list.remove('leftkidney')
    # msk_fn_list.remove('leftsurrenalgland')
    # msk_fn_list.remove('leftsurretumor')
    # msk_fn_list.remove( 'liver')
    # msk_fn_list.remove('lungs')
    # msk_fn_list.remove('pancreas')
    # msk_fn_list.remove('portalvein')
    # msk_fn_list.remove('rightkidney')
    # msk_fn_list.remove('rightsurrenalgland')
    # msk_fn_list.remove('rightsurretumor')
    # msk_fn_list.remove('spleen')
    # msk_fn_list.remove('stomach')
    # msk_fn_list.remove('venoussystem')

    msk_list_convert(input_nifti_dir, output_dicom_dir, output_dicom_dir, struct_name, msk_fn_list)
if __name__ == "__main__":
    # input_nifti_dir = r'F:\workspace\python\nii2dcm\data\nifti'
    # output_dicom_dir = r'C:\Users\shcho\Desktop\1'
    # struct_name = 'datu'
    # roi_list =['Eye-Right', 'Eye-Left','Mandible', 'Optic-Nerves-Left', 'Optic-Nerves-Right',
    #            'Parotid-Gland-Left', 'Parotid-Gland-Right', 'Brain']
    # nifti2dcm(input_nifti_dir, output_dicom_dir, struct_name, roi_list)

    # input_nifti_dir = r'C:\Users\shcho\Desktop\out\data'
    # output_dicom_dir = r'C:\Users\shcho\Desktop\out\value'
    # struct_name = 'datu'
    # roi_list=None
    # nifti2dcm(input_nifti_dir, output_dicom_dir, struct_name, roi_list)

    input_dir = r'C:\Users\shcho\Desktop\3Dircadb1.5'
    output_dir=r'C:\Users\shcho\Desktop\out\ircadb'
    struct_name = 'datu'
    nifti_multiple_label_convert(input_dir, output_dir, struct_name)