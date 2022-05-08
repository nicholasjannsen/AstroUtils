from PIL import Image
from PIL import ImageFilter
import numpy as np
im = Image.open('DK.png')

im1 = im.filter(ImageFilter.BLUR)
im2 = im.filter(ImageFilter.MinFilter(3))
im3 = im.filter(ImageFilter.MinFilter)

im1.show()

#http://pillow.readthedocs.io/en/3.4.x/reference/ImageFilter.html
