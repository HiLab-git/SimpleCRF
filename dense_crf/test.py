import dense_crf
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

Iq =np.asarray(Image.open('10.png'), np.uint8)
P = np.asarray(Image.open('10prob.png'), np.float32)
fP = P/255.0;
print 'img shape:',Iq.shape
print 'prb shape:',P.shape

param = (3.0, 10, 5, 3.0, 5, 5)
lab = dense_crf.dense_crf(Iq, fP, param)
print lab.min(), lab.max()
plt.subplot(1,3,1); plt.imshow(Iq)
plt.subplot(1,3,2); plt.imshow(P)
plt.subplot(1,3,3); plt.imshow(lab)
plt.show()

#print lab
