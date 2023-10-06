#program to read in a JPG or PNG file and output to a PPM file
from PIL import Image
import sys

def readImage(filename):
    try:
        image = Image.open(filename)
        #rotate the image 90 degrees right
        # image = image.rotate(90)
    except IOError:
        print("Error: file not found")
        sys.exit(1)
    return image

def writePPM(image, filename):
    #open filename for writing
    with open(filename, 'w') as f:
        #write header
        f.write('P3 ')
        f.write(str(image.size[0]) + ' ')
        f.write(str(image.size[1]) + ' ')
        f.write('255 \n')
         
        for x in range(image.size[1]):
            for y in range(image.size[0]):
                #get pixel color
                r, g, b = image.getpixel((y, x))
                #write pixel color
                f.write(str(r) + ' ' + str(g) + ' ' + str(b) + ' ')
            f.write('\n')
    

if __name__ == "__main__":
    image = readImage("image.jpg")
    writePPM(image, "image.ppm")