
import numpy as np
import SimpleITK as sitk
import maxflow
from PIL import Image
import matplotlib.pyplot as plt

def demo_maxflow():
    I = Image.open('../data/img2d.png')
    Iq = np.asarray(I.convert('L'), np.float32)
    P = np.asarray(Image.open('../data/prob2d.png'), np.float32)

    fP = P/255.0
    bP = 1.0-fP
    lamda = 10.0  
    sigma = 5.0
    param = (lamda, sigma)
    lab = maxflow.maxflow2d(Iq, fP, bP, param)

    plt.subplot(1,4,1); plt.axis('off'); plt.imshow(I);  plt.title('input image')
    plt.subplot(1,4,2); plt.axis('off'); plt.imshow(fP);   plt.title('probability map')
    plt.subplot(1,4,3); plt.axis('off'); plt.imshow(fP > 0.5); plt.title('initial \n segmentation')
    plt.subplot(1,4,4); plt.axis('off'); plt.imshow(lab); plt.title('CRF result')
    plt.show()

def demo_interactive_maxflow():
    I = Image.open('../data/img2d.png')
    Iq = np.asarray(I.convert('L'), np.float32)
    P = np.asarray(Image.open('../data/prob2d.png'), np.float32)
    
    fP = P/255.0
    bP = 1.0-fP
    S = np.zeros_like(fP, np.uint8)
    seed1_x = 27
    seed1_y = 32
    seed2_x = 40
    seed2_y = 52
    S[seed1_y - 2 : seed1_y + 2, seed1_x - 2 : seed1_x + 2] = 255
    S[seed2_y - 2 : seed2_y + 2, seed2_x - 2 : seed2_x + 2] = 127

    lamda = 10.0
    sigma = 5.0
    param = (lamda, sigma)
    lab   = maxflow.interactive_maxflow2d(Iq, fP, bP, S,param)
    plt.subplot(1,3,1); plt.axis('off'); plt.imshow(Iq);  plt.title('input image')
    plt.subplot(1,3,2); plt.axis('off'); plt.imshow(P); 
    plt.plot([seed1_x], [seed1_y], 'bo', markersize = 2)
    plt.plot([seed2_x], [seed2_y], 'ro', markersize = 2)
    plt.title('probability map and \n user interaction')
    plt.subplot(1,3,3); plt.axis('off'); plt.imshow(lab); plt.title('CRF result')
    plt.show()

def demo_maxflow3d():
    img_name   = "../data/atp.nii.gz"
    prob_name  = "../data/prob.nii.gz"
    save_name  = "../data/seg_auto.nii.gz"
    img_obj  = sitk.ReadImage(img_name)
    img_data = sitk.GetArrayFromImage(img_obj)
    img_data = np.asarray(img_data, np.float32)
    prob_obj = sitk.ReadImage(prob_name)
    prob_data = sitk.GetArrayFromImage(prob_obj)
    prob_data = np.asarray(prob_data, np.float32)

    fP = prob_data
    bP = 1.0-fP

    lamda = 10.0
    sigma = 5.0
    param = (lamda, sigma)
    lab = maxflow.maxflow3d(img_data, fP, bP, param)
    lab_obj = sitk.GetImageFromArray(lab)
    lab_obj.CopyInformation(img_obj)
    sitk.WriteImage(lab_obj, save_name)
    print('the segmentation has been saved to {0:}'.format(save_name))

def test_interactive_max_flow3d():
    img_name   = "../data/atp.nii.gz"
    prob_name  = "../data/prob.nii.gz"
    seed_name  = "../data/scrb.nii.gz"
    save_name  = "../data/seg_interact.nii.gz"
    img_obj  = sitk.ReadImage(img_name)
    img_data = sitk.GetArrayFromImage(img_obj)
    img_data = np.asarray(img_data, np.float32)
    prob_obj = sitk.ReadImage(prob_name)
    prob_data = sitk.GetArrayFromImage(prob_obj)
    prob_data = np.asarray(prob_data, np.float32)
    seed_obj  = sitk.ReadImage(seed_name)
    seed_data = sitk.GetArrayFromImage(seed_obj)
    seed_data[seed_data == 1] = 127
    seed_data[seed_data == 2] = 255
    seed_data = np.asarray(seed_data, np.uint8)

    fP = prob_data
    bP = 1.0-fP

    lamda = 20.0
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