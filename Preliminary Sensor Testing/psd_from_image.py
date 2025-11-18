import cv2
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt


#read image and convert to greyscale
img = cv2.imread('microscope_panorama.jpg')
grey = cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)

# Thresholding
_, thresh = cv2.threshold(grey, 100, 255, cv2.THRESH_BINARY_INV)

# Morphological closing to get whole particles; opening to get rid of noise
img_mop = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (1, 1)))
img_mop = cv2.morphologyEx(img_mop, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2, 2)))

# skip morphological operation
#img_mop = thresh.copy()

# contours
cnts, _ = cv2.findContours(img_mop, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

# scale info
#15mm is image width
scale = (img.shape[1], img.shape[0], 15e3)  # pixels, pixels, microns

# extract particles
particles = [cv2.boundingRect(cnt) for cnt in cnts]

img_cp = img.copy()

cv2.drawContours(img_cp, cnts, -1, (0, 255, 0), 2)

cv2.imshow('img', cv2.resize(img, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.imshow('gray', cv2.resize(grey, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.imshow('thresh', cv2.resize(thresh, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.imshow('img_mop',  cv2.resize(img_mop, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.imshow('img_cp',  cv2.resize(img_cp, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.waitKey(0)
cv2.destroyAllWindows()

# save smallest size distribution
sizes = [min(p[2], p[3]) / scale[2] * 15e3 for p in particles]
pd.Series(sizes).to_csv('psd_from_image.csv', index=False,header=['size_um'])