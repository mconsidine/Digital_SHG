"""
@author: Andrew Smith
Version 6 August 2021

"""
from DSHG_solex_util import *

import skimage
import skimage.feature
import sys

import math
import numpy as np
import matplotlib.pyplot as plt

from skimage import data
from skimage import transform
from skimage import filters
from skimage.transform import downscale_local_mean
import cv2

import scipy
#from ellipse import LsqEllipse
from matplotlib.patches import Ellipse
#from DSHG_fit_ellipse_ALT import * #MattC
from DSHG_ellipse_JFP import * #MattC
from scipy import ndimage #MattC

NUM_REG = 1 #6 # include biggest NUM_REG regions in fit


def rot(x):
    return np.array([[np.cos(x), np.sin(x)], [-np.sin(x), np.cos(x)]])

def get_correction_matrix(phi, r):
    """
    IN: phi, ellipse axes ratio (height / width)
    OUT: correction matrix
    """
    stretch_matrix = rot(phi) @ np.array([[r, 0], [0, 1]]) @  rot(-phi)
    theta = np.arctan(-stretch_matrix[0, 1] / stretch_matrix[1, 1]) 
    unrotation_matrix = rot(theta)
    correction_matrix = unrotation_matrix @ stretch_matrix
    return np.linalg.inv(correction_matrix), theta


def get_new_correction_matrix(shear_angle, ratio_corr):      # Modification Jean-Francois, pur shear and scaling without rotation
    """
    IN: phi, ellipse axes ratio (height / width)
    OUT: correction matrix
    """
    #phi = phi - np.pi/2.0 #Jean-Francois
    correction_matrix = np.array([[ratio_corr, -ratio_corr*np.tan(shear_angle)], [0, 1]])
    theta = 0
    return np.linalg.inv(correction_matrix), theta

def dofit(points, options): #MattC
    """IN : numpy points coordinates
    OUT : center, width, height, phi, fit informations
    """
    
    #print("Value of hflip : ", options['hflip'])  #MattC
    
    #if options['hflip']==False:
    if 1:
        reg = LsqEllipse().fit(points)
        center, width, height, phi = reg.as_parameters()
        return center, width, height, phi, reg.return_fit(n_points=100)
    '''
    else:
        reg = ls_ellipse(points[0],points[1])
        verbose=False #MattC
        params = polyToParams(reg,verbose)
        center = (params[0],params[1])
        width = params[2]
        height = params[3]
        phi = params[4]
        print(center)

        n_points=100
        t = np.linspace(0, 2 * np.pi, n_points)

        x = (center[0] + width * np.cos(t) * np.cos(phi) - height * np.sin(t) * np.sin(phi))
        y = (center[1] + width * np.cos(t) * np.sin(phi) + height * np.sin(t) * np.cos(phi))
        return center, width, height, phi, np.c_[x, y]
   '''

def two_step(points, options): #MattC
    """Launch twice an ellipse fit. One with all edge points, one with only tresholded values.
    IN : numpy array of edge points.
    OUT : np.array(center), height, phi, ratio, points_tresholded, ellipse_points
    """
    center, width, height, phi, _ = dofit(points, options) #MattC
    mat, _ = get_correction_matrix(phi, height / width)
    Xr = mat @ (points - np.array(center)).T * height
    values = np.linalg.norm(Xr, axis = 0) - 1
    #print(np.mean(values), np.std(values), max(values), min(values))
    anomaly_threshold = max(values)
    points_tresholded = points[values > -max(values)]
    center, width, height, phi, ellipse_points = dofit(points_tresholded, options) #MattC

    logme('Width : ' + str(width))
    logme('Height : ' + str(height))
    w = 1
    h = height / width
    m = np.tan(0.5*np.pi - phi)
    A = 1/(w*w) + (m*m)/(h*h)
    B = 2*m*np.sqrt(w*w*m*m + h*h)/(h*h)
    x1 = -B/(2*A)
    #x1 = -h*h*np.abs(m)/np.sqrt(w*w*m*m + h*h)
    y1 = np.sqrt(h*h*(1 - (x1*x1)/(w*w)))
    logme('m : ' + str(m))
    logme('x1 : ' + str(x1))
    logme('x1 : ' + str(x1*width))
    logme('y1 : ' + str(y1))
    x2 = x1*np.cos(phi) - y1*np.sin(phi)
    y2 = x1*np.sin(phi) + y1*np.cos(phi)
    logme('phi : ' + str(math.degrees(phi)))
    logme('x2 : ' + str(x2))
    logme('y2 : ' + str(y2))
    #shear_angle = np.arctan2(y2,x2)
    shear_angle = np.arctan(y2/x2)
    logme('shear tilt : ' + str(math.degrees(shear_angle)))
    reg = LsqEllipse().fit(points)
    a, b, c, d, f, g  = reg.ellipse_coeff()        
    logme('a : ' + str(a))
    logme('b : ' + str(b))
    logme('c : ' + str(c))
    logme('d : ' + str(d))
    logme('f : ' + str(f))
    logme('g : ' + str(g))
    x0 = (c*d - b*f)/(b*b - a*c)
    y0 = (a*f - b*d)/(b*b - a*c)
    logme('x0 : ' + str(x0))
    logme('y0 : ' + str(y0))
    x_corr = np.sqrt(height*height*np.cos(phi)*np.cos(phi) + width*width*np.sin(phi)*np.sin(phi)) 
    y_corr = np.sqrt(height*height*np.sin(phi)*np.sin(phi) + width*width*np.cos(phi)*np.cos(phi))
    #y_corr = np.abs(y1*width)
    logme('x_corr : ' + str(x_corr))
    logme('y_corr : ' + str(y_corr))
    logme('ratio : ' + str(y_corr/x_corr))

    mat, _ = get_correction_matrix(phi, height / width)
    Xr = mat @ (points_tresholded - np.array(center)).T * height
    values = np.linalg.norm(Xr, axis = 0) - 1
    #print(np.mean(values), np.std(values), max(values), min(values))
    #ratio = width / height                                                                         # original
    ratio_corr = y_corr/x_corr                                                                      # Modification Jean-Francois
    #if options['hflip'] == True: #MattC kludge?
    #    scale_multiplier = (height+width)/(width)
    #    ratio = ratio*scale_multiplier
    #return np.array(center), height, phi, ratio, points_tresholded, ellipse_points                 # original
    return np.array(center), height, shear_angle, ratio_corr, points_tresholded, ellipse_points     # Modification Jean-Francois

#def correct_image(image, phi, ratio, center, print_log = False):   # original
def correct_image(image, shear_angle, ratio_corr, center, print_log = False):
    """correct image geometry. TODO : a rotation is made instead of a tilt
    IN : numpy array, float, float, numpy array (2 elements)
    OUT : numpy array, numpy array (2 elements)
    """
    
    #mat, theta = get_correction_matrix(phi, ratio)                     # original
    mat, theta = get_new_correction_matrix(shear_angle, ratio_corr)     # Modification Jean-Francois

    if print_log:
        print('unrotation angle theta = ' + "{:.3f}".format(math.degrees(theta)) + " degrees")
        np.set_printoptions(suppress=True)
        #logme('Y/X ratio : ' + "{:.3f}".format(ratio))                                     # original
        logme('Y/X ratio : ' + "{:.3f}".format(ratio_corr))
        #logme('Tilt angle : ' + "{:.3f}".format(math.degrees(phi)) + " degrees")           # original
        logme('Tilt angle : ' + "{:.3f}".format(math.degrees(shear_angle)) + " degrees")
        logme('Linear transform correction matrix: \n' + str(mat))
        np.set_printoptions(suppress=False)
    mat3 = np.zeros((3, 3))
    mat3[:2, :2] = mat
    mat3[2, 2] = 1
    corners = np.array([[0,0], [0, image.shape[0]], [image.shape[1], 0], [image.shape[1], image.shape[0]]])
    new_corners = (np.linalg.inv(mat) @ corners.T).T # use inverse because we represent mat3 as inverse of transform
    new_h = np.max(new_corners[:, 1]) - np.min(new_corners[:, 1])
    new_w = np.max(new_corners[:, 0]) - np.min(new_corners[:, 0])
    mat3 = mat3 @ np.array([[1, 0, np.min(new_corners[:, 0])], [0, 1, np.min(new_corners[:, 1])], [0,   0,    1]]) # apply translation to prevent clipping
    my_transform = transform.ProjectiveTransform(matrix=mat3)
    corrected_img = transform.warp(image, my_transform, output_shape = (np.ceil(new_h),np.ceil(new_w)), cval  = image[0, 0])
    corrected_img = (2**16*corrected_img).astype(np.uint16) # note : 16-bit output
    new_center = (np.linalg.inv(mat) @ center.T).T - np.array([np.min(new_corners[:, 0]), np.min(new_corners[:, 1])])
    return corrected_img, new_center

def get_flood_image(image):
    """
    Return an image, where all the pixels brighter than a threshold
    are made saturated, and all those below average are zeroed.
    the threshhold is chosen as the local minimum of the pixel-brightness
    histogram of the image. As a backup, the average brightness is used if
    a local minimum cannot be found.
    IN: original image
    OUT: modified image
    """
    
    
    thresh = .9 * np.sum(image) / (image.shape[0] * image.shape[1])
    print('thresh=', thresh)
    img_blurred = cv2.blur(image, ksize=(5, 5))
    n, bins = np.histogram(img_blurred.flatten(), bins=20)
    bottom = -1
    for i in range(19, 1, -1):
        if n[i-1] < n [i]:
            tip = i
            break
    for i in range(tip, 1, -1):
        if n[i-1] > n [i]:
            bottom = i
            break
    thresh2 = thresh if bottom == -1 else bins[bottom]
    print('thresh2=', thresh2)
            
    img_blurred[img_blurred < thresh2] = 0
    img_blurred[img_blurred >= thresh2] = 65000
    return img_blurred

def get_edge_list(image, sigma = 2):
    """from a picture, return a numpy array containing edge points
    IN : frame as numpy array, integer
    OUT : numpy array
    TODO: simplify this function?
    """
    if sigma <= 0:
        logme('ERROR: could not find any edges')
        return image, (-1, -1, -1)
    
    low_threshold = np.median(cv2.blur(image, ksize=(5, 5))) / 10 
    high_threshold = low_threshold*1.5
    print('using thresholds:', low_threshold, high_threshold)
    #image = get_flood_image(image) #MattC
    #cv2.namedWindow('test images', cv2.WINDOW_NORMAL) #MattC
    #cv2.moveWindow('test images', 0, 0)
    #cv2.resizeWindow('test images',int(image.shape[1] * 1), int(image.shape[0] * 1))
    #cv2.imshow('test images',image)
    #cv2.waitKey(4000)  # affiche et continue
    #cv2.destroyAllWindows()    
    edges = skimage.feature.canny(
        image=image,
        sigma=sigma,
        low_threshold=low_threshold,
        high_threshold=high_threshold,
    )
    raw_X = np.argwhere(edges)
    labelled, nf = scipy.ndimage.measurements.label(edges, structure  = [[1,1,1], [1,1,1],[1,1,1]])
    if nf == 0:
        return get_edge_list(image, sigma = sigma - 0.5) # try again with less blur, hope it will work
    region_sizes = [-1] + [np.sum(labelled == i) for i in range(1, nf+1)]
    filt = np.zeros(edges.shape)
    for label in sorted(region_sizes, reverse = True)[:min(nf, NUM_REG)]:
        filt[labelled == region_sizes.index(label)] = 1
        
    X = np.argwhere(filt) # find the non-zero pixels
    
    x_min, y_min, x_max, y_max = np.min(X[:, 0]), np.min(X[:, 1]), np.max(X[:, 0]), np.max(X[:, 1])
    dx = x_max - x_min
    dy = y_max - y_min
    crop = 0.015

    mask = np.zeros(filt.shape)
    mask[int(x_min+dx*crop):int(x_max-dx*crop), :] = 1
    filt *= mask
    X = np.argwhere(filt) # find the non-zero pixels again

    x_min, y_min, x_max, y_max = np.min(X[:, 0]), np.min(X[:, 1]), np.max(X[:, 0]), np.max(X[:, 1])
            
    X = np.array(X, dtype='float')
    return np.array([X, raw_X], dtype=object) 

def ellipse_to_circle(image, options):
    """from an entire sun frame, compute ellipse fit and return a circularise picture and center coordinates
    IN : numpy array, dictionnayr of options
    OUt :numpy array, numpy array (2 elements)
    """
    image = image / 65536 # assume 16 bit
    factor = 4
    processed = get_edge_list(downscale_local_mean(image, (factor,factor))) * factor# down-scaled, then upscaled back
    X, raw_X = processed[0], processed[1]
    #center, height, phi, ratio, X_f, ellipse_points = two_step(X)                  # original
    center, height, shear_angle, ratio_corr, X_f, ellipse_points = two_step(X,options)      # Modification Jean-Francois #MattC
    center = np.array([center[1], center[0]])
    
    #fix_img, center = correct_image(image, phi, ratio, center, print_log = True)               # original
    fix_img, center = correct_image(image, shear_angle, ratio_corr, center, print_log = True)   # Modification Jean-Francois
    #fix_img = np.fliplr(np.copy(fix_img)) #MattC  kludge
    #print('Shear angle : ', shear_angle,np.degrees(shear_angle))
    #fix_img = ndimage.rotate(fix_img, np.degrees(shear_angle), reshape=False) #MattC kludge
    
    if options['flag_display']:
        fig, ax = plt.subplots(ncols=2, nrows = 2)
        ax[0][0].imshow(image, cmap=plt.cm.gray) #MattC fix_img?
        ax[0][0].set_title('uncorrected image', fontsize = 11)
        ax[0][1].imshow(image, cmap=plt.cm.gray) #MattC fix_img?
        ax[0][1].plot(raw_X[:, 1], raw_X[:, 0], 'ro', label = 'edge detection')
        #ax[0][1].set_xlim([0, image.shape[1]])     # Modification Jean-Francois, for having the image same orientation as ax[0][0]
        #ax[0][1].set_ylim([0, image.shape[0]])     # Modification Jean-Francois, for having the image same orientation as ax[0][0]
        ax[0][1].legend()
        #ax[1][1].plot(X_f[:, 1], image.shape[0] - X_f[:, 0], 'ro', label = 'filtered edges')       # Modification Jean-Francois
        ax[1][1].plot(X_f[:, 1], X_f[:, 0], 'ro', label = 'filtered edges')
        #ax[1][1].plot(ellipse_points[:, 1], image.shape[0] - ellipse_points[:, 0], color='b', label = 'ellipse fit')   # Modification Jean-Francois
        ax[1][1].plot(ellipse_points[:, 1], ellipse_points[:, 0], color='b', label = 'ellipse fit')
        ax[1][1].set_xlim([0, image.shape[1]]) #MattC fix_img?
        #ax[1][1].set_ylim([0, image.shape[0]])     # Modification Jean-Francois
        ax[1][1].set_ylim([image.shape[0], 0])      # Modification Jean-Francois, for having the scale pointing "down" #MattC fix_img?
        ax[1][1].legend()
        ax[1][0].imshow(fix_img, cmap=plt.cm.gray)
        ax[1][0].set_title('geometrically corrected image', fontsize=11)
        ax[0][1].set_title('remember to close this window \n by pressing the "X"', color = 'red')

        #creating a timer object to auto-close plot after some time
        def close_event():
            plt.close()
        timer = fig.canvas.new_timer(interval = options['tempo']) 
        timer.add_callback(close_event)
        timer.start()
        plt.gca().set_aspect('equal', adjustable='box')     # Modification Jean-Francois, graphic with XY same scale
        plt.show()
  
    if (ratio_corr < 1.0):
        circle = (center[0], center[1], height)             # Modification Jean-Francois, radius == height
    else:
        circle = (center[0], center[1], height*ratio_corr)  # Modification Jean-Francois, radius == width
    #return fix_img, circle, ratio, phi                     # original
    return fix_img, circle, ratio_corr, shear_angle         # Modification Jean-Francois

    
