from PIL import ImageEnhance
from PIL import Image
from PIL import ImageFilter
import numpy as np


im = Image.open('DSC_0492.JPG')

enhancer = ImageEnhance.Sharpness(im)

for i in range(5):
    factor = i / 4.0
    enhancer.enhance(factor).show("Sharpness %f" % factor)

#http://pillow.readthedocs.io/en/3.4.x/reference/ImageEnhance.html

