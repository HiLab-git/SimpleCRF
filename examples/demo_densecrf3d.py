import nibabel
import denseCRF3D
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def demo_densecrf3d():
    data_root  = '../dependency/densecrf3d/applicationAndExamples/example/'
    I1Nii = nibabel.load(data_root + 'Flair_normalized.nii.gz')
    I2Nii = nibabel.load(data_root + 'DWI_normalized.nii.gz')
    PNii  = nibabel.load(data_root + 'lesionProbMap.nii.gz')
    I1 = I1Nii.get_data()
    I2 = I2Nii.get_data()
    P  = PNii.get_data()

    # convert input normalized image to intenstiy range of [0, 255]
    I = np.asarray([I1, I2], np.float32)
    I = np.transpose(I, [1, 2, 3, 0])
    I = (I + 3)/6.0 * 255
    I[I < 0] = 0
    I[I > 255] = 255
    I = np.asarray(I, np.uint8)

    # probability map for each class
    P = np.asarray([1.0 - P, P], np.float32)
    P = np.transpose(P, [1, 2, 3, 0])

    dense_crf_param = {}
    dense_crf_param['MaxIterations'] = 2.0
    dense_crf_param['PosW'] = 2.0
    dense_crf_param['PosRStd'] = 5
    dense_crf_param['PosCStd'] = 5
    dense_crf_param['PosZStd'] = 5
    dense_crf_param['BilateralW'] = 3.0
    dense_crf_param['BilateralRStd'] = 5.0
    dense_crf_param['BilateralCStd'] = 5.0
    dense_crf_param['BilateralZStd'] = 5.0
    dense_crf_param['ModalityNum'] = 2
    dense_crf_param['BilateralModsStds'] = (5.0,5.0)

    lab = denseCRF3D.densecrf3d(I, P, dense_crf_param)
    labNii = nibabel.Nifti1Image(lab, np.eye(4))
    nibabel.save(labNii, data_root + 'results/lesionSegMap.nii.gz')

if __name__ == "__main__":
    demo_densecrf3d()



