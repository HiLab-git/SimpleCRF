import dense_crf3d
import numpy as np
import nibabel

import dense_crf
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def test_dense_crf3d():
    I1Nii = nibabel.load('densecrf_3d/applicationAndExamples/example/Flair_normalized.nii.gz')
    I2Nii = nibabel.load('densecrf_3d/applicationAndExamples/example/DWI_normalized.nii.gz')
    PNii = nibabel.load('densecrf_3d/applicationAndExamples/example/lesionProbMap.nii.gz')
    I1 = I1Nii.get_data()
    I2 = I2Nii.get_data()
    P  = PNii.get_data()

    I = np.asarray([I1, I2], np.float32)
    I = np.transpose(I, [1, 2, 3, 0])
    I = (I + 3)/6.0 * 255
    I[I < 0] = 0
    I[I > 255] = 255
    I = np.asarray(I, np.uint8)

    P = np.asarray([1.0 - P, P], np.float32)
    P = np.transpose(P, [1, 2, 3, 0])

    dense_crf_param = {}
    dense_crf_param['MaxIterations'] = 2.0
    dense_crf_param['PosRStd'] = 3.0
    dense_crf_param['PosCStd'] = 3.0
    dense_crf_param['PosZStd'] = 3.0
    dense_crf_param['PosW'] = 1.0
    dense_crf_param['BilateralRStd'] = 5.0
    dense_crf_param['BilateralCStd'] = 5.0
    dense_crf_param['BilateralZStd'] = 5.0
    dense_crf_param['ModalityNum'] = 2
    dense_crf_param['BilateralW'] = 3.0
    dense_crf_param['BilateralModsStds'] = (3.0,3.0)

    lab = dense_crf3d.dense_crf3d(I, P, dense_crf_param)
    labNii = nibabel.Nifti1Image(lab, np.eye(4))
    nibabel.save(labNii, 'densecrf_3d/applicationAndExamples/example/results/lesionSegMap.nii.gz')

if __name__ == "__main__":
    test_dense_crf3d()



