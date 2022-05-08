import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
img = mpimg.imread('DK.png')
imgplot = plt.imshow(img)
plt.show()

lum_img = img[:, :, 0]
imgplot = plt.imshow(lum_img)
plt.colorbar()
plt.show()

#https://matplotlib.org/users/image_tutorial.html

