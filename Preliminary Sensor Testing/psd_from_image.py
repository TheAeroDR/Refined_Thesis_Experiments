import cv2
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt


#read image and convert to greyscale
img = cv2.imread('microscope_panorama.jpg')
grey = cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)

# Thresholding
_, thresh = cv2.threshold(grey, 120, 255, cv2.THRESH_BINARY_INV)

_,th2 = cv2.threshold(grey,0,255,cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
th3 = cv2.adaptiveThreshold(grey,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY_INV,501,2)
# Morphological closing to get whole particles; opening to get rid of noise
#img_mop = cv2.morphologyEx(th2, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (1, 1)))
#img_mop = cv2.morphologyEx(th2, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2, 2)))

# skip morphological operation
img_mop = th2.copy()

# contours
cnts, _ = cv2.findContours(img_mop, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

# scale info
#15mm is image width
scale = (15e3 / img.shape[1])  # micron per pixel

# extract particles
particles = [cv2.boundingRect(cnt) for cnt in cnts]

img_cp = img.copy()

cv2.drawContours(img_cp, cnts, -1, (0, 255, 0), 2)

cv2.imshow('img', cv2.resize(img, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.imshow('gray', cv2.resize(grey, dsize=(0, 0), fx=0.5, fy=0.5))
#cv2.imshow('thresh', cv2.resize(thresh, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.imshow('thresh2', cv2.resize(th2, dsize=(0, 0), fx=0.5, fy=0.5))
#cv2.imshow('thresh3', cv2.resize(th3, dsize=(0, 0), fx=0.5, fy=0.5))

cv2.imshow('img_mop',  cv2.resize(img_mop, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.imshow('img_cp',  cv2.resize(img_cp, dsize=(0, 0), fx=0.5, fy=0.5))
cv2.waitKey(0)
cv2.destroyAllWindows()

cv2.imwrite('thresholds.jpg', img_mop)
cv2.imwrite('contours.jpg', img_cp)

# save smallest size distribution
sizes = [min(p[2], p[3]) * scale for p in particles]
pd.Series(sizes).to_csv('psd_from_image.csv', index=False,header=['size_um'])