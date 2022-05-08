from PIL import Image
import numpy as np
im = Image.open('DK_resized.png')
pixelMap = im.load()

print pixelMap[100,100] 

img = Image.new(im.mode, im.size)
pixelsNew = img.load()
for i in range(img.size[0]):
    for j in range(img.size[1]):

        if (pixelMap[i,j] == (100, 200, 100, 255)):
            pixelsNew[i,j] = (200, 203, 30, 255)
        else:
            pixelsNew[i,j] = pixelMap[i,j]

img.save("out.jpg")

# HTML color codes.
#http://htmlcolorcodes.com/
