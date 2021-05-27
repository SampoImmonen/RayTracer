import sys

import imageio
from PIL import Image

if __name__ == '__main__':
    #default imagename
    filename = "image.ppm"
    if len(sys.argv) >= 2:
        filename = sys.argv[1]

    im = imageio.imread(filename)
    image = Image.fromarray(im)
    image.show()