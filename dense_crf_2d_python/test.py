import dense_crf
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def test_dense_crf():
    I = Image.open('../data/img2d.png')
    Iq = np.asarray(I.convert('L'), np.uint8)
    P  = np.asarray(Image.open('../data/prob2d.png'), np.float32)
    fP = P/255.0
    print('img shape:',Iq.shape)
    print('prb shape:',P.shape)

    w1    = 5.0  # weight of bilateral term
    alpha = 10   # spatial std
    beta  = 10   # rgb  std
    w2    = 3.0  # weight of spatial term
    gamma = 10   # spatial std
    it    = 5.0  # iteration
    param = (w1, alpha, beta, w2, gamma, it)
    lab = dense_crf.dense_crf(Iq, fP, param)
    print(lab.min(), lab.max())
    plt.subplot(1,3,1); plt.imshow(I)
    plt.subplot(1,3,2); plt.imshow(P)
    plt.subplot(1,3,3); plt.imshow(lab)
    plt.show()

if __name__ == "__main__":
    test_dense_crf()