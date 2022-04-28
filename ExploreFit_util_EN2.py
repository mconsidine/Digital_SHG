# from skimage.measure import profile_line
from skimage import draw
import sys
# import os
import cv2
import numpy as np
from solex_util import logme
import PySimpleGUI as sg
#  from pynput import keyboard
import scipy.interpolate  # MattC for congrid
import scipy.ndimage  # MattC for congrid

#https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

#https://scipython.com/blog/binning-a-2d-array-in-numpy/
def rebin2(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

# https://scipy-cookbook.readthedocs.io/items/Rebinning.html
def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array(a.shape)
    ndims = len(a.shape)
    if len(newdims) != ndims:
        print("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = np.asarray(newdims, dtype=float)
    dimlist = []

    if method == 'neighbour':
        for i in range(ndims):
            base = np.indices(newdims)[i]
            dimlist.append((old[i] - m1) / (newdims[i] - m1) \
                           * (base + ofs) - ofs)
        cd = np.array(dimlist).round().astype(int)
        newa = a[list(cd)]
        return newa

    elif method in ['nearest', 'linear']:
        # calculate new dims
        for i in range(ndims):
            base = np.arange(newdims[i])
            dimlist.append((old[i] - m1) / (newdims[i] - m1) \
                           * (base + ofs) - ofs)
        # specify old dims
        olddims = [np.arange(i, dtype=np.float) for i in list(a.shape)]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d(olddims[-1], a, kind=method)
        newa = mint(dimlist[-1])

        trorder = [ndims - 1]
        trorder.extend(range(ndims - 1))

        print(trorder)

        for i in range(ndims - 2, -1, -1):
            newa = newa.transpose(trorder)

            mint = scipy.interpolate.interp1d(olddims[i], newa, kind=method)
            newa = mint(dimlist[i])

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose(trorder)

        return newa
    elif method in ['spline']:
        oslices = [slice(0, j) for j in old]
        oldcoords = np.ogrid[oslices]
        nslices = [slice(0, j) for j in list(newdims)]
        newcoords = np.mgrid[nslices]

        newcoords_dims = list(range(np.ndim(newcoords)))
        # make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None


def thebigroutine(image_data, thefit, LineRecal, MattC_basefich):
    # https://stackoverflow.com/questions/70246383/image-profile-intensity
    # https://stackoverflow.com/questions/28340950/opencv-how-to-draw-continously-with-a-mouse
    # https://stackoverflow.com/questions/58405119/how-to-resize-the-window-obtained-from-cv2-imshow
    global image, image_copy, img, img_copy, column_selected, fitted_column, calcd_adjustment, \
        SETTINGBWL, SETTINGBWR, BWleft, BWright, bwcalcd_ladjustment, bwcalcd_radjustment, y1temp, y2temp

    profilewindowsize = 400

    SHGINFOADDED = False
    TargetWavelength = 0
    AngPerPixel = 0

    SETTINGY1 = False
    SETTINGY2 = False

    y1temp = -999999
    y2temp = -999999

    SETTINGBWL = False
    SETTINGBWR = False

    BWleft = -999999
    BWright = -999999

    bwcalcd_ladjustment = -999999
    bwcalcd_radjustment = -999999

    def UI_SHGinfo():
        sg.theme('Dark2')
        sg.theme_button_color(('white', '#500000'))

        layout = [

            [sg.Text('Target wavelength in Ang (blank to ignore)', size=(35, 1)),
            # sg.Input(default_text='6562.8', size=(8, 1), key='-TgtLambda-', enable_events=True)],
             sg.Combo(['6562.8', '3933.692'], key='-TgtLambda-')],
            [sg.Text('Ang per pixel (blank to ignore)', size=(35, 1)),
            # sg.Input(default_text='0.034', size=(8, 1), key='-AperP-', enable_events=True)],
             sg.Combo(['0.034', '0.078'], key='-AperP-')],
            [sg.Button('OK'), sg.Cancel()]
        ]

        window = sg.Window('Processing', layout, finalize=True)
        window.BringToFront()

        while True:
            event, values = window.read()
            if event == sg.WIN_CLOSED or event == 'Cancel':
                sys.exit()

            if event == 'OK':
                break
        window.close()

        return values['-TgtLambda-'], values['-AperP-']

    TargetWavelength, AngPerPixel = UI_SHGinfo()
    try:
        targetlambda = float(TargetWavelength) if not TargetWavelength is None else 0
    except ValueError:
        targetlambda = 0
    try:
        angstromperpixel = float(AngPerPixel) if not TargetWavelength is None else 0
    except ValueError:
        angstromperpixel = 0

    if (angstromperpixel > 0):
        SHGINFOADDED = True
        logme('Initial targeted wavelength : ' + str(targetlambda))
        logme('Estimate angstroms per pixel : ' + str(angstromperpixel))

    print("Click in 'image' window to select a slice.  Profile will appear")
    print("  in  smaller window, which can be resized.  Click in smaller")
    print("  window to select a point in the profile.  Status info is in ")
    print("  terminal window.  Pressing 'c' continues with selection. Pressing")
    print("  'x' cancels selection and continues with original algorithm.")

    cv2.namedWindow("profile", cv2.WINDOW_KEEPRATIO)
    cv2.setWindowProperty('profile', 1, cv2.WINDOW_NORMAL)
    cv2.moveWindow("profile", 20, 20)
    cv2.namedWindow("image", cv2.WINDOW_NORMAL) #,  cv2.WINDOW_KEEPRATIO)
    cv2.moveWindow("image", 400, 10)
    cv2.resizeWindow('image', int(image_data.shape[0]), int(image_data.shape[1]))
    image = image_data[..., np.newaxis]
    imgdatatype = np.uint8
    imgmaxlimit = 255
    if np.max(image) > 255:
        imgdatatype = np.uint16
        imgmaxlimit = 65535
    image = (((image - 0 * np.min(image)) / (imgmaxlimit + 0 * np.max(image) - 0 * np.min(image))) * 255).astype(
        np.uint8)

    image = cv2.merge([image, image, image])
    image_copy = image.copy()
    img = np.zeros((profilewindowsize, image.shape[1], 3), np.uint8)  # for profile
    img_copy = img
    cv2.resizeWindow("profile", int(img.shape[1]), int(img.shape[0]) + 15)
    cv2.resizeWindow("image", int(image.shape[1]), int(image.shape[0]) + 15)
    slice_selected = image.shape[0] // 2  # initialize to mid point of slit image
    origl_fit_column = (int(np.asarray(thefit)[slice_selected, 0]) + LineRecal)
    column_selected = origl_fit_column
    fitted_column = origl_fit_column
    calcd_adjustment = 0

    def select_upper_y(event, x, y, flags, param):
        return

    def select_lower_y(event, x, y, flags, param):
        return

    # mouse callback function to select column on profile chart
    def select_column(event, x, y, flags, param):
        global column_selected, img, calcd_adjustment, fitted_column, img_copy, \
            SETTINGBWL, SETTINGBWR, BWleft, BWright, bwcalcd_ladjustment, bwcalcd_radjustment

        bwcalcd_adjustment = 0

        # If left mouse button is clicked, start of line
        if (event == cv2.EVENT_LBUTTONDOWN):
            if ((SETTINGBWL == False) & (SETTINGBWR == False) ) == True:
                # print("not setting either")
                # print(SETTINGBWL, SETTINGBWR)
                column_selected = int(round(x, 0))
                calcd_adjustment = column_selected - origl_fit_column
                if (BWleft != -999999):
                    BWleft = column_selected + bwcalcd_ladjustment
                    BWright = column_selected + bwcalcd_radjustment

                print("Current column selected : ", column_selected)
                print(" fitted column : ", origl_fit_column)

                if SHGINFOADDED:
                    print("Adjustment based on selection : ", calcd_adjustment, " columns; ",
                          calcd_adjustment * angstromperpixel, " estd Angstroms")
                    if (targetlambda > 0):
                        print("Estimated selected wavelength : ", targetlambda + calcd_adjustment * angstromperpixel)
                else:
                    print("Adjustment based on selection : ", calcd_adjustment, " columns ")

            elif ((SETTINGBWL == True) & (SETTINGBWR == False) ) == True:
                print("Setting bandwidth")
                # print(SETTINGBWL, SETTINGBWR)
                BWleft = int(round(x, 0))
                bwcalcd_ladjustment = BWleft - column_selected
                bwcalcd_radjustment = -bwcalcd_ladjustment
                BWright = column_selected + bwcalcd_radjustment
                if SHGINFOADDED:
                    print("Bandwidth based on selection : ", np.abs(2 * (bwcalcd_ladjustment)) + 1, " columns; ",
                          (np.abs(2 * (bwcalcd_ladjustment)) + 1) * angstromperpixel, " estd Angstroms")
                    if (targetlambda > 0):
                        print("left Estimated selected wavelength : ",
                              targetlambda + bwcalcd_ladjustment * angstromperpixel)
                        print("right Estimated selected wavelength : ",
                              targetlambda - bwcalcd_ladjustment * angstromperpixel)
                else:
                    print("Bandwidth based on selection : ", np.abs(2 * bwcalcd_ladjustment) + 1, " columns ")

            '''
            elif ((SETTINGBWL==True) & (SETTINGBWR==False))==True:
                print("setting left")
                print(SETTINGBWL, SETTINGBWR)
                BWleft = int(round(x,0))
                bwcalcd_ladjustment = BWleft-column_selected
                BWright = column_selected-bwcalcd_ladjustment #because its a neg number
                if SHGINFOADDED:
                    print("left Adjustment based on selection : ", bwcalcd_ladjustment," columns; ",
                    bwcalcd_ladjustment*angstromperpixel, " estd Angstroms")
                    if (targetlambda>0):
                        print("left Estimated selected wavelength : ",targetlambda+bwcalcd_ladjustment*angstromperpixel)
                        print("right Estimated selected wavelength : ",targetlambda-bwcalcd_ladjustment*angstromperpixel)                        
                else:
                    print("left Adjustment based on selection : ", bwcalcd_ladjustment, " columns ")
               
            elif ((SETTINGBWL==False) & (SETTINGBWR==True))==True:
                print('setting right')
                print(SETTINGBWL, SETTINGBWR)
                BWright = int(round(x,0))
                bwcalcd_radjustment = BWright-column_selected
                BWleft = column_selected-bwcalcd_radjustment #because its a pos number
                if SHGINFOADDED:
                    print("right Adjustment based on selection : ", bwcalcd_radjustment," columns; ",
                    bwcalcd_radjustment*angstromperpixel, " estd Angstroms")
                    if (targetlambda>0):
                        print("left Estimated selected wavelength : ",targetlambda-bwcalcd_radjustment*angstromperpixel)                         
                        print("right Estimated selected wavelength : ",targetlambda+bwcalcd_radjustment*angstromperpixel)                   
                else:
                    print("right Adjustment based on selection : ", bwcalcd_radjustment, " columns ")
            '''
            '''                     
            elif ((SETTINGBWL==True) & (SETTINGBWR==True))==True:
                print("shouldn't be here")
            else:
                print("really shouldn;t be here")
            '''
        elif (event == cv2.EVENT_MOUSEMOVE):  # LBUTTONUP
            x_start = x
            y_start = 0
            x_end = x_start
            y_end = profilewindowsize
            img = cv2.line(img_copy.copy(), (x_start, y_start), (x_end, y_end), (0, 0, 255), 1)

            if (column_selected >= 0):
                img = cv2.line(img, (column_selected, 0), (column_selected, profilewindowsize), (255, 0, 0), thickness=1)
                if (BWleft >= 0):
                    img = cv2.line(img, (BWleft, 0), (BWleft, profilewindowsize), (255, 127, 63), thickness=1)
                if (BWright >= 0):
                    img = cv2.line(img, (BWright, 0), (BWright, profilewindowsize), (63, 127, 255), thickness=1)

            cv2.imshow("profile", img)

        return

    # mouse callback function to show slice selected from profile
    def print_coords(event, x, y, flags, param):

        # Global variables needed
        global image, image_copy, x_start, y_start, slice_selected, img, img_copy, fitted_column, y1temp ,y2temp

        # If left mouse button is clicked, start of line
        if (event == cv2.EVENT_LBUTTONDOWN):
            if (SETTINGY1 or SETTINGY2) == True:
                if SETTINGY1 == True:
                    y1temp = y
                    print("Setting y1 to : ", y1temp)
                elif SETTINGY2 == True:
                    y2temp = y
                    print("Setting y2 to : ", y2temp)
            '''
            x_start = 0
            y_start = y
            slice_selected = y_start
            x_end = image.shape[1] #x
            y_end = y_start #y
            fitted_column = int(np.asarray(thefit)[slice_selected,0])+LineRecal
            image = cv2.line(image_copy.copy(), (x_start, y_start), (x_end, y_end), (255, 255, 255), 1)
            image = cv2.line(image, (fitted_column+calcd_adjustment, y_start-10), 
                          (fitted_column+calcd_adjustment,y_start+15), (255,0,0), 1)
            line = np.transpose(np.array(draw.line(x_start, y_start, x_end-1, y_end)))
            data = image_copy.copy()[line[:, 1], line[:, 0]] 

            print("Slice : ",y_start)
            print("Calcd column for fit : ",fitted_column)
            #show profile of slice selected.  It's profilewindowsize pixels high and as wide as the selection
            #should be a better way to do this to get larger chart?
            img = np.zeros((profilewindowsize,image.shape[1],3), imgdatatype)
            draw_x = np.linspace(0, len(data), len(data))
            draw_y = profilewindowsize-((data-np.min(data))/(np.max(data)-np.min(data)))*profilewindowsize #0,0 is in upper left
            draw_points = (np.asarray([draw_x, draw_y[:,0]]).T).astype(np.int32)   # needs to be int32 and transposed
            cv2.polylines(img, [draw_points], False, (255,255,255), thickness=1)  # args: image, points, closed, color
            cv2.line(img, (fitted_column,0), (fitted_column,profilewindowsize), (0,255,0), thickness=1)
            img_copy = img
            cv2.imshow("profile", img)
            '''
        elif (event == cv2.EVENT_MOUSEMOVE):
            x_start = 0
            y_start = y
            x_end = image.shape[1]
            y_end = y_start

            image = cv2.line(image_copy.copy(), (x_start, y_start), (x_end, y_end), (255, 255, 255),
                             thickness=3)

            if (y1temp != -999999):
                image = cv2.line(image, (x_start, y1temp), (x_end, y1temp), (0, 127, 127), 3)
            if (y2temp != -999999):
                image = cv2.line(image, (x_start, y2temp), (x_end, y2temp), (0, 127, 255), 3)

            if (column_selected > 0):
                image = cv2.line(image,
                    #(int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment, y_start - 25),
                    #(int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment, y_start + 25),
                    (int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment, y_start - 25),
                    (int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment, y_start + 25),
                    (255, 0, 0), 3)
                if (BWleft > 0):
                    image = cv2.line(image,
                    #int(np.asarray(thefit)[y_start, 0]) + LineRecal + bwcalcd_ladjustment, y_start - 25),
                    #(int(np.asarray(thefit)[y_start, 0]) + LineRecal + bwcalcd_ladjustment,
                    (int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment + bwcalcd_ladjustment,
                    y_start - 25),
                    (int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment + bwcalcd_ladjustment,
                    y_start + 25), (255, 127, 63), 3)
                if (BWright > 0):
                    image = cv2.line(image,
                    #int(np.asarray(thefit)[y_start, 0]) + LineRecal + bwcalcd_radjustment, y_start - 25),
                    #(int(np.asarray(thefit)[y_start, 0]) + LineRecal + bwcalcd_radjustment,
                    (int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment + bwcalcd_radjustment,
                     y_start - 25),
                    (int(np.asarray(thefit)[y_start, 0]) + LineRecal + calcd_adjustment + bwcalcd_radjustment,
                    y_start + 25), (63, 127, 255), 3)
            cv2.imshow("image", image)

    # Set up mouse callback function
    cv2.setMouseCallback("image", print_coords)
    cv2.setMouseCallback("profile", select_column)

    print(image.shape[0], image.shape[1])
    fitted_column = origl_fit_column
    line = np.transpose(np.array(draw.line(0, slice_selected, image.shape[1] - 1, slice_selected)))
    data = image[line[:, 1], line[:, 0]]
    # weights = 1-((data-np.min(data))/(np.max(data)-np.min(data)))
    draw_x = np.linspace(0, len(data), len(data))
    draw_y = profilewindowsize - ((data - np.min(data)) / (np.max(data) - np.min(data))) * profilewindowsize  # 0,0 is in upper left
    draw_points = (np.asarray([draw_x, draw_y[:, 0]]).T).astype(np.int32)  # needs to be int32 and transposed
    cv2.polylines(img, [draw_points], False, (127, 127, 127), thickness=1)  # args: image, points, closed, color
    cv2.line(img, (fitted_column, 0), (fitted_column, profilewindowsize), (0, 255, 0), thickness=1)
    img_copy = img
    cv2.imshow("profile", img)
    '''
    #https://nitratine.net/blog/post/how-to-make-hotkeys-in-python/

    def function_0():
        print('c was pressed')
        return False

    def function_1():
        print('Function 1 activated')
        return True

    def function_2():
        print('Function 2 activated')
        return True

    def on_press(key):
        if key == keyboard.Key.esc:
            return False  # stop listener
        try:
            k = key.char  # single-char keys
        except:
            k = key.name  # other keys
        if k in ['shift+s', '1', '2', 'left', 'right']:  # keys of interest
            # self.keys.append(k)  # store it in global-like variable
            print('Key pressed: ' + k)
            #return False  # stop listener; remove this if want more keys
        else:
            if k == 's' :
                print('lower key2 : ', k)
            elif k == 'S':
                print('upper key2 : ', k)

    listener = keyboard.Listener(on_press=on_press)
    listener.start()  # start to listen on a separate thread
    listener.join()  # remove if main thread is polling self.keys
    '''

    # Loop until the 'c' key is pressed
    while True:
    

        # Display image; wait for keypress
        key = cv2.waitKey(1) & 0xFF

        # If 'c' key is pressed, break from loop
        if key == ord('c'):
            break
        elif key == ord('0'):
            SETTINGY2 = False
            SETTINGY1 = False
            print("Deselecting borders = ", SETTINGY2)
        elif key == ord('2'):
            SETTINGY2 = True
            SETTINGY1 = False
            print("Selecting second border = ", SETTINGY2)
        elif key == ord('1'):
            SETTINGY2 = False
            SETTINGY1 = True
            print("Selecting first border = ", SETTINGY1)
        elif key == ord('b'):
            SETTINGBWL = not SETTINGBWL
            SETTINGBWR = False
            print("Selecting band=", SETTINGBWL)
            print(BWleft, BWright)
        elif key == ord('m'):
            SETTINGBWL = False
            SETTINGBWR = False
            print("Selecting main band")
        elif key == ord('n'):
            SETTINGBWL = False
            SETTINGBWR = False
            print("Resetting band selection to none")
            BWleft = -999999
            BWright = -999999
            # break
        elif key == ord('x'):
            column_selected = -999999
            fitted_column = origl_fit_column
            calcd_adjustment = 0
            SETTINGBWL = False
            SETTINGBWR = False
            SETTINGY2 = False
            SETTINGY1 = False
            BWleft = -999999
            BWright = -999999
            print("Selection discarded.  Original fit will be used.")
            break

    cv2.destroyAllWindows()

    if (column_selected > 0):
        if (BWleft > 0) & (BWright > 0):
            # adjust range so that bandwidth selected is close to FWHM
            FWHMleft = BWleft
            FWHMright = BWright
            print("FWHM boundaries : ", FWHMleft, FWHMright)
            # BWleft = int(column_selected-2*(FWHMright-FWHMleft+1)/(2.354820045)) #2*sqrt(2*ln(2))
            # BWright = int(column_selected+2*(FWHMright-FWHMleft+1)/(2.354820045)) #2*sqrt(2*ln(2))

            bwcalcd_radjustment = int(2*round((FWHMright-FWHMleft+1)/2.354820045)) # for gaussian
            #bwcalcd_radjustment = int(round(2 * (FWHMright - FWHMleft + 1) / (2.0)))
            bwcalcd_ladjustment = -bwcalcd_radjustment
            BWleft = column_selected + bwcalcd_ladjustment
            BWright = column_selected + bwcalcd_radjustment

            print('Left band selected : ', BWleft)
            print('Left adjustment : ', bwcalcd_ladjustment)
            print('Right band selected : ', BWright)
            print('Right adjustment : ', bwcalcd_radjustment)
            print('FWHM Bandwidth : ', (FWHMright - FWHMleft + 1) * angstromperpixel)
            # print("BW sum : ", np.sum(data[BWleft:BWright,0]))
            # print("BW mean : ",np.mean(data[BWleft:BWright,0]))
            # print("BW sd : ",np.std(data[BWleft:BWright,0]))
            # print("BW max : ", np.max(data[BWleft:BWright,0]))
            # print("BW : ",data[BWleft:BWright,0])
            #gsigma = (FWHMright - FWHMleft + 1) / 2 #for Lorentzian ?
            ##print( (2*np.sqrt(2*np.log(2))) ) #for gaussian
            gsigma = (FWHMright-FWHMleft+1)/(2.354820045) # for gaussian
            # print(gsigma)
            gmean = column_selected
            # print(gmean)
            grange = range(BWleft, BWright + 1)  # add one because of python range function
            grange = np.asarray(grange)
            # print(grange)
            # for Lorentzian:
            ##gweights = 1.0/(1.0+(grange-gmean)**2/gsigma**2)
            # for Gaussian,pdf = :
            gweights = 1.0 / (gsigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((grange - gmean) / gsigma) ** 2)
            # print("first gweights : ",gweights)
            ''' # maybe not needed to scale?? #MattC
            gsumweights = np.sum(gweights)
            print("first sum gweight : ", gsumweights)
            gweights = gweights / gsumweights
            print("second gweights : ",gweights)
            print("second sum gweights : ", np.sum(gweights))
            '''
            # print("data : ",data)
            # print("drawy : ", draw_y)
            # print("Weights in band : ", weights[BWleft:BWright,0])
            # print("Weights in band : ", weights[BWleft:BWright,0]/np.sum(weights[BWleft:BWright,0]))

        # print("Slice chosen : ", slice_selected)
        logme('Slice chosen : ' + str(slice_selected))
        # print(type(column_selected))
        # print("Column selected : ", column_selected)
        logme('Column selected : ' + str(column_selected))
        # fitted_column = int(np.asarray(thefit)[slice_selected,0])+LineRecal #add back LineRecal to compare to image data
        # print("  vs fitted value : ",fitted_column)
        logme('  vs fitted value : ' + str(fitted_column))
        calcd_adjustment = column_selected - fitted_column
        # print("Calcd shift vs fitted value : ",calcd_adjustment)
        logme('Calcd shift vs fitted value : ' + str(calcd_adjustment))

        if SHGINFOADDED == True:
            # print("Adjustment based on selection : ", calcd_adjustment," columns ",
            # calcd_adjustment*angstromperpixel, "estd Angstroms")
            logme('Adjustment based on selection : ' + str(calcd_adjustment) + ' columns; ' +
                  str(calcd_adjustment * angstromperpixel) + ' estd Angstroms')
            # print("Estimated selected wavelength : ",targetlambda+calcd_adjustment*angstromperpixel)
            if (targetlambda > 0):
                logme('Estimated selected wavelength : ' + str(targetlambda + calcd_adjustment * angstromperpixel))

        if (calcd_adjustment != 0):
            newfitvals = np.asarray(thefit)[:, 0] + calcd_adjustment  # already accounts for LineRecal
            newfit = list(zip(newfitvals.astype('int'), np.asarray(thefit)[:, 1], np.asarray(thefit)[:, 2]))

            # create 3-channel grayscale of original for overlay purposes only
            gray = image_data[..., np.newaxis]
            gray = (((gray - 0 * np.min(gray)) / (imgmaxlimit + 0 * np.max(gray) - 0 * np.min(gray))) * 255).astype(
                np.uint8)

            # create overlay of original fit
            im2_data = cv2.merge([gray, gray, gray])
            pts = np.asarray(thefit)[:, (0, 2)].astype(np.int32)
            # offset=5 #positive offsets are to right in the image; concave direction
            # pluspts = np.asarray(thefit)[:,(0,2)].astype(np.int32)+offset
            # minuspts = np.asarray(thefit)[:,(0,2)].astype(np.int32)-offset

            # cv2.polylines(im2_data,[minuspts],isClosed=False, color = (255, 0, 0), thickness=1) #BGR
            cv2.polylines(im2_data, [pts], isClosed=False, color=(0, 255, 0), thickness=1)  # BGR
            # cv2.polylines(im2_data,[pluspts],isClosed=False, color = (0, 0, 255), thickness=1) #BGR

            cv2.namedWindow('Original fit', cv2.WINDOW_NORMAL)
            cv2.moveWindow('Original fit', 0, 0)
            scale = 2
            cv2.resizeWindow('Original fit', int(im2_data.shape[1] * scale), int(im2_data.shape[0] * scale))

            while (1):
                cv2.imshow('Original fit', im2_data)
                key = cv2.waitKey(1) & 0xFF

                # If 'c' key is pressed, break from loop
                if key == ord("c"):
                    break

            cv2.destroyAllWindows()

            cv2.imwrite(MattC_basefich + 'origfitmap.png', im2_data)

            ##gray=image_data[...,np.newaxis]
            ##gray=(((gray-np.min(gray))/(np.max(gray)-np.min(gray)))*255).astype(imgdatatype)

            # create overlay of adjusted fit
            im2_data = cv2.merge([gray, gray, gray])
            pts = np.asarray(newfit)[:, (0, 2)].astype(np.int32)
            # offset=5 #positive offsets are to right in the image; concave direction
            # pluspts = np.asarray(newfit)[:,(0,2)].astype(np.int32)+offset
            # minuspts = np.asarray(newfit)[:,(0,2)].astype(np.int32)-offset

            # cv2.polylines(im2_data,[minuspts],isClosed=False, color = (255, 0, 0), thickness=1) #BGR
            cv2.polylines(im2_data, [pts], isClosed=False, color=(255, 0, 0), thickness=1)  # BGR
            # cv2.polylines(im2_data,[pluspts],isClosed=False, color = (0, 0, 255), thickness=1) #BGR

            cv2.namedWindow('Adjusted fit', cv2.WINDOW_NORMAL)
            cv2.moveWindow('Adjusted fit', 0, 0)
            scale = 2
            cv2.resizeWindow('Adjusted fit', int(im2_data.shape[1]), int(im2_data.shape[0]))

            while (1):
                cv2.imshow('Adjusted fit', im2_data)
                key = cv2.waitKey(1) & 0xFF

                # If 'c' key is pressed, break from loop
                if key == ord("c"):
                    break

            cv2.destroyAllWindows()

            cv2.imwrite(MattC_basefich + 'newfitmap.png', im2_data)

            print("Guide images saved ...")
            if (BWleft > 0):
                print("New fit used, bandwidth added")
                print(BWleft, bwcalcd_ladjustment, BWright, bwcalcd_radjustment, calcd_adjustment, column_selected)
                return newfit, True, BWleft, bwcalcd_ladjustment, BWright, bwcalcd_radjustment, column_selected, gweights, y1temp, y2temp
            else:
                print("New fit, no bandwidth added")
                return newfit, True, 0, 0, 0, 0, column_selected, 0, y1temp, y2temp
        else:
            if (BWleft > 0) & (BWright > 0):
                print("Origl fit used, bandwidth added")
                print(BWleft, bwcalcd_ladjustment, BWright, bwcalcd_radjustment, calcd_adjustment, column_selected)
                return thefit, True, BWleft, bwcalcd_ladjustment, BWright, bwcalcd_radjustment, column_selected, gweights, y1temp, y2temp
            else:
                print("Selection = origl fit, no bandwidth added")
                return thefit, False, 0, 0, 0, 0, origl_fit_column, 0, y1temp, y2temp  # newfit True column_selected
    else:
        logme('No new column selected.  Original fit will be used.')
        return thefit, False, 0, 0, 0, 0, origl_fit_column, 0, 0, 0
