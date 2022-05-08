from PIL import Image

img = Image.open("DK.png")
img2 = img.rotate(45, expand=True)
img2.save("DK_45.png")

# x1, y1, x2, y2
img3 = img2.crop((430, 900, 900, 1300))
img3.save("DK_45_crop.jpg")

