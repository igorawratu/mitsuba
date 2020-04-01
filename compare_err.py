import cv2
import sys
import math
import numpy as np

def findrange(v):
	return 1 - v, v


def float2RGB(v):
	r = 0
	g = 0
	b = 0

	if(v < 0.5):
		b, g = findrange(v * 2)
	else:
		g, r = findrange((v - 0.5) * 2)

	return r, g, b


if len(sys.argv) < 4:
	print("Not enough arguments, format expected is compare_err.py <file1> <file2> <outputfile>")
	sys.exit()

filename1 = sys.argv[1]
filename2 = sys.argv[2]
output = sys.argv[3]


inv = 1 / 2.2

img1 = cv2.imread(filename1, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
img2 = cv2.imread(filename2, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)

i1height, i1width, i1channels = img1.shape
i2height, i2width, i2channels = img2.shape

if i1height != i2height or i1width != i2width or i1channels != i2channels:
	print("Input images do not match")
	sys.exit()

accumulated_square_err = 0

diff_img = np.zeros(img1.shape, img1.dtype)
output_img = np.copy(img1)
cv2.pow(img1, inv, img1)
cv2.pow(img2, inv, img2)

for rowidx, row in enumerate(img1):
	for colidx, pix in enumerate(row):
		pixel_err = 0
		for idx, channel in enumerate(pix):
			v1 = img1[rowidx, colidx, idx]
			v2 = img2[rowidx, colidx, idx]
			dif = v1 - v2
			accumulated_square_err += dif * dif
			pixel_err += abs(dif)
		pixel_err /= i1channels

		col = float2RGB(pixel_err)
		diff_img[rowidx, colidx, 0] = col[2]
		diff_img[rowidx, colidx, 1] = col[1]
		diff_img[rowidx, colidx, 2] = col[0]

rmse = math.sqrt(accumulated_square_err / (i1height * i1width * i1channels))

f = open(output + "_errors.txt", "a+")
f.write(str(rmse) + "\n")
f.close()

cv2.imwrite(output + "_dif_" + str(rmse) + ".exr", diff_img)
cv2.imwrite(output + "_" + str(rmse) + ".exr", output_img)