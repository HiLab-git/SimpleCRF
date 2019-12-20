
import numpy as np
import SimpleITK as sitk
import maxflow
from PIL import Image
import matplotlib.pyplot as plt

def demo_maxflow():
    I = Image.open('../data/brain.png')
    Iq = np.asarray(I.convert('L'), np.float32)
    P = np.asarray(Image.open('../data/brain_mask.png').convert('L'), np.float32) / 255

    fP = 0.5 + (P - 0.5) * 0.8
    bP = 1.0 - fP
    lamda = 20.0  
    sigma = 10.0
    param = (lamda, sigma)
    lab = maxflow.maxflow2d(Iq, fP, bP, param)

    plt.subplot(1,3,1); plt.axis('off'); plt.imshow(I);  plt.title('input image')
    plt.subplot(1,3,2); plt.axis('off'); plt.imshow(fP);   plt.title('initial \n segmentation')
    plt.subplot(1,3,3); plt.axis('off'); plt.imshow(lab); plt.title('CRF result')
    plt.show()

def demo_interactive_maxflow():
    I = Image.open('../data/brain.png')
    Iq = np.asarray(I.convert('L'), np.float32)
    P = np.asarray(Image.open('../data/brain_mask.png').convert('L'), np.float32) / 255

    fP = 0.5 + (P - 0.5) * 0.8
    bP = 1.0 - fP

    S = np.asarray(Image.open('../data/brain_scrb.png').convert('L'))
    scrb = np.zeros_like(S)
    scrb[S==170] = 127
    scrb[S==255] = 255

    lamda = 30.0  
    sigma = 8.0
    param = (lamda, sigma)
    lab = maxflow.interactive_maxflow2d(Iq, fP, bP, scrb, param)

    plt.subplot(1,3,1); plt.axis('off'); plt.imshow(I);  plt.title('input image')
    plt.subplot(1,3,2); plt.axis('off'); plt.imshow(fP);   plt.title('initial \n segmentation')
    plt.subplot(1,3,3); plt.axis('off'); plt.imshow(lab); plt.title('CRF result')
    plt.show()

def demo_maxflow3d():
    img_name   = "../data/2013_12_1_img.nii.gz"
    prob_name  = "../data/2013_12_1_init.nii.gz"
    save_name  = "../data/seg_auto.nii.gz"
    img_obj  = sitk.ReadImage(img_name)
    img_data = sitk.GetArrayFromImage(img_obj)
    img_data = np.asarray(img_data, np.float32)
    prob_obj = sitk.ReadImage(prob_name)
    prob_data = sitk.GetArrayFromImage(prob_obj)
    prob_data = np.asarray(prob_data, np.float32)

    fP = 0.5 + (prob_data - 0.5) * 0.8
    bP = 1.0-fP

    lamda = 10.0
    sigma = 15.0
    param = (lamda, sigma)
    lab = maxflow.maxflow3d(img_data, fP, bP, param)
    lab_obj = sitk.GetImageFromArray(lab)
    lab_obj.CopyInformation(img_obj)
    sitk.WriteImage(lab_obj, save_name)
    print('the segmentation has been saved to {0:}'.format(save_name))

def test_interactive_max_flow3d():
    img_name   = "../data/2013_12_1_img.nii.gz"
    prob_name  = "../data/2013_12_1_init.nii.gz"
    seed_name  = "../data/2013_12_1_scrb.nii.gz"
    save_name  = "../data/seg_interact.nii.gz"
    img_obj  = sitk.ReadImage(img_name)
    img_data = sitk.GetArrayFromImage(img_obj)
    img_data = np.asarray(img_data, np.float32)
    prob_obj = sitk.ReadImage(prob_name)
    prob_data = sitk.GetArrayFromImage(prob_obj)
    prob_data = np.asarray(prob_data, np.float32)

    fP = 0.5 + (prob_data - 0.5) * 0.8
    bP = 1.0-fP

    seed_obj  = sitk.ReadImage(seed_name)
    seed_data = sitk.GetArrayFromImage(seed_obj)
    seed_data[seed_data == 1] = 127
    seed_data[seed_data == 2] = 255
    seed_data = np.asarray(seed_data, np.uint8)

    lamda = 10.0
    sigma = 15.0
    param = (lamda, sigma)
    lab = maxflow.interactive_maxflow3d(img_data, fP, bP, seed_data, param)
    lab_obj = sitk.GetImageFromArray(lab)
    lab_obj.CopyInformation(img_obj)
    sitk.WriteImage(lab_obj, save_name)
    print('the segmentation has been saved to {0:}'.format(save_name))

if __name__ == '__main__':
    print("example list")
    print(" 0 -- 2D max flow without interactions")
    print(" 1 -- 2D max flow with interactions")
    print(" 2 -- 3D max flow without interactions")
    print(" 3 -- 3D max flow with interactions")
    print("please enter the index of an example:")
    method = input()
    method = "{0:}".format(method)
    if(method == '0'):
        demo_maxflow()
    elif(method == '1'):
        demo_interactive_maxflow()
    elif(method == '2'):
        demo_maxflow3d()
    elif(method == '3'):
        test_interactive_max_flow3d()
    else:
        print("invalid number : {0:}".format(method))
