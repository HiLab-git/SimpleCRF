import max_flow
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def test_max_flow():
    Iq =np.asarray(Image.open('img.png').convert('L'), np.float32);

    P = np.asarray(Image.open('prob.png'), np.float32);

    fP = P/255.0;
    bP = 1.0-fP;
    param = (5.0,3.5)
    lab = max_flow.max_flow(Iq, fP, bP, param)
    plt.subplot(1,3,1); plt.imshow(Iq)
    plt.subplot(1,3,2); plt.imshow(P)
    plt.subplot(1,3,3); plt.imshow(lab)
    plt.show()
#print lab

def test_interactive_max_flow():
    Iq =np.asarray(Image.open('10.png').convert('L'), np.float32);
    
    P = np.asarray(Image.open('10prob.png'), np.float32);
    
    fP = P/255.0;
    bP = 1.0-fP;
    S = np.asarray(Image.open('10_seed.png'),np.uint8);
    param = (5.0,3.5)
    lab = max_flow.interactive_max_flow(Iq, fP, bP, S,param)
    plt.subplot(1,3,1); plt.imshow(Iq)
    plt.subplot(1,3,2); plt.imshow(P)
    plt.subplot(1,3,3); plt.imshow(lab)
    plt.show()

if __name__ == '__main__':
#    test_max_flow()
    test_interactive_max_flow()