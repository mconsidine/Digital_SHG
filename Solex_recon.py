# -*- coding: utf-8 -*-
"""
@author: Valerie Desnoux
with improvements by Andrew Smith
contributors: Jean-Francois Pittet, Jean-Baptiste Butet, Pascal Berteau, Matt Considine
Version 8 September 2021

------------------------------------------------------------------------
reconstruction on an image from the deviations between the minimum of the line and a reference line

calcul sur une image des ecarts simples entre min de la raie et une ligne de reference
-------------------------------------------------------------------------

"""

from solex_util import *
from ser_read_video import *
from ellipse_to_circle import ellipse_to_circle, correct_image
from ExploreFit_util_EN2 import rebin, rebin2, congrid, thebigroutine #MattC
from scipy.optimize import curve_fit #MattC

# read video and return constructed image of sun using fit and LineRecal
def read_video_improved(serfile, fit, LineRecal, options, theflattener):
    rdr = ser_reader(serfile)
    ih, iw = rdr.ih, rdr.iw
    FrameMax = rdr.FrameCount
    disk_list = [np.zeros((ih, FrameMax), dtype=rdr.infiledatatype)  # MattC avi
                 for _ in options['shift']]
    #disk_list_d1 = np.copy(disk_list)  # MattC for doppler
    #disk_list_d2 = np.copy(disk_list)  # MattC for doppler
    #disk_list_c = np.copy(disk_list)  # MattC for doppler
    disk_list_bw = np.zeros((ih, FrameMax), dtype=rdr.infiledatatype)  # MattC for bandwidth

    if options['flag_display']:
        cv2.namedWindow('disk', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('disk', FrameMax // 3, ih // 3)
        cv2.moveWindow('disk', 200, 0)
        cv2.namedWindow('image', cv2.WINDOW_NORMAL)
        cv2.moveWindow('image', 0, 0)
        cv2.resizeWindow('image', int(iw), int(ih))

    col_indeces = []

    for shift in options['shift']:
        ind_l = (np.asarray(fit)[:, 0] + np.ones(ih)
                 * (LineRecal + shift)).astype(int)

        # CLEAN if fitting goes too far
        ind_l[ind_l < 0] = 0
        ind_l[ind_l > iw - 2] = iw - 2
        ind_r = (ind_l + np.ones(ih)).astype(int)
        col_indeces.append((ind_l, ind_r))

    left_weights = np.ones(ih) - np.asarray(fit)[:, 1]
    right_weights = np.ones(ih) - left_weights

    # start of mods for doppler, bandwidth  # MattC

    if m1 > 0:  # MattC
        the_weights = (np.ones((ih, (m3 - m1 + 1))).astype(float)) * m6  # +1
        test_coords = []

        for i, shift in enumerate(options['shift']):
            #ind_bw = (np.asarray(fit)[:, 0] + np.ones(ih) * (LineRecal + 0 * shift)).astype(int)  # MattC
            #ind_bw[ind_bw < 0] = 0
            #ind_bw[ind_bw > iw - 2] = iw - 2
            ind_bw, _ = col_indeces[i]
            startcol = np.add(ind_bw, m2)
            endcol = np.add(ind_bw, m4)
            # img_suff = (['.fits', '_d1_' + str(shift_dop) + '.fits', '_d2_' + str(shift_dop) + '.fits',
            #             '_cont_' + str(shift_cont) + '.fits', '_bw_' + str(np.abs(m2) + np.abs(m4) + 1) + '.fits'])
            testcoords = np.linspace(startcol, endcol, np.abs(m2) + np.abs(m4) + 1).astype(int).T  # +1
            testcoords[testcoords < 0] = 0
            testcoords[testcoords > iw - 2] = iw - 2
            test_coords.append((testcoords))

    # end of mods for doppler, bandwidth  # MattC

    # lance la reconstruction du disk a partir des trames
    print('reader num frames:', rdr.FrameCount)
    while rdr.has_frames():
        img = rdr.next_frame()
        if options['flag_display'] and rdr.FrameIndex % 10 == 0:
            cv2.imshow('image', img)
            if cv2.waitKey(1) == 27:
                cv2.destroyAllWindows()
                sys.exit()

        for i in range(len(options['shift'])):
            ind_l, ind_r = col_indeces[i]
            left_col = img[np.arange(ih), ind_l]
            right_col = img[np.arange(ih), ind_r]
            IntensiteRaie = left_col * left_weights + right_col * right_weights
            #disk_list[i][:, rdr.FrameIndex] = np.divide(IntensiteRaie, theflattener).astype(rdr.infiledatatype) #MattC
            disk_list[i][:, rdr.FrameIndex] = IntensiteRaie.astype(rdr.infiledatatype) #MattC avi
            # start of bandwidth mods MattC

            if m1 > 0:  # MattC
                imgsubset = np.take_along_axis(img, test_coords[i], axis=1)
                newvalues = imgsubset * the_weights
                newvalues = np.sum(newvalues, axis=1)
                IntensiteRaie_bw = np.round(newvalues, 0)
                #disk_list_bw[:, rdr.FrameIndex] = np.divide(IntensiteRaie_bw, theflattener).astype(rdr.infiledatatype) #MattC
                disk_list_bw[:, rdr.FrameIndex] = IntensiteRaie_bw.astype(rdr.infiledatatype)  # MattC
            # end of bandwidth mods MattC
        if options['flag_display'] and rdr.FrameIndex % 10 == 0:
            # disk_list[1] is always shift = 0
            cv2.imshow('disk', disk_list[1])
            if cv2.waitKey(
                    1) == 27:                     # exit if Escape is hit
                cv2.destroyAllWindows()
                sys.exit()
    if m1 > 0:  # MattC
        options['shift'] = options['shift'] + [str(int(m4))+'FWHM']  #MattC adding bw 2*m4+1
        disk_list.append(disk_list_bw)

    print("rdr infiletype : ",rdr.infiledatatype)

    for i in range(len(options['shift'])):  # MattC this should just change 8bit to 16bit
        imgmax = np.max(disk_list[i])
        imgmin = np.min(disk_list[i])
        print('img min max dtype : ',i, imgmin, imgmax, disk_list[i].dtype)
        ##disk_list[i] = ((disk_list[i]-imgmin)*64000/(imgmax-imgmin)).astype('uint16')
        #disk_list[i] = (disk_list[i]*(64000/255)).astype('uint16')

    if options['binsize'] > 1 :  # MattC
        binsize=options['binsize']
        print('old ih and iw : ',ih, iw)
        print('length disk_list : ',len(disk_list))
        print('length options shift : ',len(options['shift']))
        print('options binsize : ',options['binsize'])
        print('binsize : ',binsize)
        disk_list2 = [np.zeros((ih//binsize,FrameMax), dtype='uint16') for _ in options['shift']]
        print("here is range len options shift ",range(len(options['shift'])))
        print("here is length disklist ", len(disk_list))
        for i in range(len(options['shift'])):
            img = rebin2(disk_list[i], (binsize, binsize))
            #origmax = np.max(disk_list[i])
            #img = congrid(disk_list[i],(ih // binsize, FrameMax // binsize), method='spline', centre=True)
            #img = ((img-0)*origmax/(np.max(img)-0)).astype('uint16') #spline can give value > 65535 so reset
            disk_list2[i] = img
        ih,iw = disk_list2[1].shape
        print('new ih and iw : ',ih, iw)
        return disk_list2, ih, iw, rdr.FrameCount
    else:
        return disk_list, ih, iw, rdr.FrameCount


def make_header(rdr):
    # initialisation d'une entete fits (etait utilisé pour sauver les trames
    # individuelles)
    hdr = fits.Header()
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = 32
    hdr['NAXIS'] = 2
    hdr['NAXIS1'] = rdr.iw
    hdr['NAXIS2'] = rdr.ih
    hdr['BZERO'] = 0
    hdr['BSCALE'] = 1
    hdr['BIN1'] = 1
    hdr['BIN2'] = 1
    hdr['EXPTIME'] = 0
    return hdr

def compute_flattener(animg, imgh, afit, ay1, ay2):  # MattC
    corrintprofile = np.zeros((imgh,1))
    MinRayfit = animg[np.arange(0, imgh), (np.asarray(afit)[0:imgh, 0]).astype(np.int)]
    popt, __ = curve_fit(func2,np.asarray(afit)[0:imgh, 2],MinRayfit,p0=np.ones(3,))  # np.ones sets poly
    v2 = func2(np.asarray(afit)[0:imgh, 2],*popt) #.astype('int')

    mawindow = 15 #5
    corrintprofile = np.divide(MinRayfit, v2)
    corrintprofile = moving_average(corrintprofile, mawindow)

    profbuffer = 100 #18
    #corrintprofile[np.abs(corrintprofile) > 1.5 ] = 0
    #now fix the tails ...
    a=[corrintprofile[min(ay2+1+profbuffer,len(corrintprofile)-1)]]*(mawindow-1)
    corrintprofile=np.concatenate((corrintprofile,a))
    corrintprofile[0:max(ay1-1-profbuffer,0)] = corrintprofile[max(ay1-1-profbuffer,0)] #18 was 25
    print(min(ay2+1+profbuffer,ih))
    print(ih)
    print(imgh)
    print("length of array", len(corrintprofile))
    corrintprofile[min(ay2+1+profbuffer,ih-1):ih-1] = corrintprofile[min(ay2+1+profbuffer,ih-1)]

    # plt.plot(corrintprofile)
    # plt.show()
    return corrintprofile

# compute mean image of video
def compute_mean(serfile):
    """IN : serfile path"
    OUT :numpy array
    """
    rdr = ser_reader(serfile)
    logme('Width, Height : ' + str(rdr.Width) + ' ' + str(rdr.Height))
    logme('Number of frames : ' + str(rdr.FrameCount))
    my_data = np.zeros((rdr.ih, rdr.iw), dtype='uint64')
    my_datacube = [np.zeros((rdr.ih, rdr.iw), dtype=rdr.infiledatatype)  # MattC
                 for _ in range(rdr.FrameCount)]
    i = 0  # MattC
    while rdr.has_frames():
        img = rdr.next_frame()
        my_data += img
        my_datacube[i] = img  # MattC
        i += 1  # MattC
    return (my_data / rdr.FrameCount).astype(rdr.infiledatatype), my_datacube  # MattC avi

def func2(x, *p):  # MattC
    """Polynomial fitting function of arbitrary degree."""
    poly = 0.
    for i, n in enumerate(p):
       poly += n * x**i
    return poly

def moving_average(x, w):  # MattC
        return np.convolve(x, np.ones(w), 'valid') / w

def compute_mean_return_fit(serfile, options, LineRecal=1): 
    global hdr, ih, iw
    global m1,m2,m3,m4,m5,m6,m7,m8
    """
    ----------------------------------------------------------------------------
    Reconstuit l'image du disque a partir de l'image moyenne des trames et
    des trames extraite du fichier ser avec un fit polynomial
    Corrige de mauvaises lignes et transversallium

    basefich: nom du fichier de base de la video sans extension, sans repertoire
    shift: ecart en pixel par rapport au centre de la raie pour explorer longueur d'onde decalée
    ----------------------------------------------------------------------------
    """
    flag_display = options['flag_display']
    # first compute mean image
    # rdr is the ser_reader object
    mean_img, thedatacube = compute_mean(serfile) #MattC

    """
    ----------------------------------------------------------------------------
    Calcul polynome ecart sur une image au centre de la sequence
    ----------------------------------------------------------------------------
    """

    
    #savefich=basefich+'_mean'

    if options['save_fit']:
        DiskHDU = fits.PrimaryHDU(mean_img, header=hdr)
        DiskHDU.writeto(basefich0 + '_mean.fits', overwrite='True')

    # affiche image moyenne
    if flag_display:
        cv2.namedWindow('Video mean', cv2.WINDOW_NORMAL)  # MattC avi
        cv2.resizeWindow('Video mean', iw, ih) # MattC avi
        cv2.moveWindow('Video mean', 100, 0) # MattC avi
        cv2.imshow('Video mean', mean_img) # MattC avi
        if cv2.waitKey(2000) == 27:                     # exit if Escape is hit
            cv2.destroyAllWindows()
            sys.exit()

        cv2.destroyAllWindows()

    y1, y2 = detect_bord(mean_img, axis=1, offset=0)  # MattC  5
    logme('Vertical limits y1, y2 : ' + str(y1) + ' ' + str(y2))

    if (y1>0.2*ih): #MattC
        y1 = 5
    if (y2<0.8*ih): #MattC
        y2 = ih-5

    PosRaieHaut = y1

    PosRaieBas = y2

    """
    -----------------------------------------------------------
    Trouve les min intensité de la raie
    -----------------------------------------------------------
    """
    # construit le tableau des min de la raie a partir du haut jusqu'en bas
    MinOfRaie = []

    for i in range(PosRaieHaut, PosRaieBas):
        line_h = mean_img[i, :]
        MinX = line_h.argmin()
        MinOfRaie.append([MinX, i])
        #print('MinOfRaie x,y', MinX,i)

    # best fit d'un polynome degre 2, les lignes y sont les x et les colonnes
    # x sont les y
    np_m = np.asarray(MinOfRaie)
    xm, ym = np_m.T
    # LineRecal=xm.min()

    p = np.polyfit(ym, xm, 2)
    
    # calcul des x colonnes pour les y lignes du polynome
    a = p[0]
    b = p[1]
    c = p[2]
    fit = []
    # ecart=[]
    for y in range(0, ih):
        x = a * y**2 + b*y + c
        deci = x - int(x)
        fit.append([int(x) - LineRecal, deci, y])

    adjustflag=False #MattC
    thenewfit, adjustflag, m1, m2, m3, m4, m5, m6, m7, m8 = thebigroutine(mean_img, fit, LineRecal, basefich0) #MattC
    print(adjustflag,"m1 : ",m1,m2,m3,m4,m5,m6,m7,m8,"y1 : ",y1,y2,"Haut : ",PosRaieHaut,PosRaieBas,"ih :",ih)
    if m7>0:
        print("m7 ",m7)
        y1 = m7
        PosRaieHaut = y1
    if m8>0:
        print("m8 ",m8)
        y2 = m8
        PosRaieBas = y2

    if adjustflag : #MattC
        print("Origl fit : ",np.asarray(fit)[1:5,]) #MattC
        print("Adjusted fit : ", np.asarray(thenewfit)[1:5,]) #MattC
        fit=thenewfit #MattC

    corrintprofile = compute_flattener(mean_img, ih, fit, y1, y2)
    # plt.plot(corrintprofile)
    # plt.show()
    framenums = len(thedatacube)
    disk_flat = np.zeros((ih, framenums), dtype=np.float) #, dtype=np.uint16)
    for i in range(0,framenums):
        # disk_flat[:, i] = compute_flattener(thedatacube[i], ih, fit, y1, y2)
        disk_flat[:, i] = corrintprofile
        # if i == (framenums // 2):
        #     print("disk f ranges ",np.min(disk_flat[:,i]), np.max(disk_flat[:,i]),i)
        #     plt.plot(disk_flat[:, i])
        #     plt.show()

    # print("disk f ranges ",np.min(disk_flat), np.max(disk_flat))
    # disk_flat[disk_flat>3]=1
    # disk_flat[disk_flat < 0] = 1

    # diskf = (disk_flat-np.min(disk_flat))*64000/(np.max(disk_flat)-np.min(disk_flat))
    cv2.imwrite(basefich0 + '_flatdisk_MattC.png', disk_flat.astype(np.uint16))

    if 2==1:
        MinRayfit = mean_img[np.arange(0, ih), (np.asarray(fit)[0:ih, 0]).astype(np.int)] #
        popt, __ = curve_fit(func2,np.asarray(fit)[0:ih, 2],MinRayfit,p0=np.ones(3,))  # np.ones sets poly
        v2 = func2(np.asarray(fit)[0:ih, 2],*popt) #.astype('int')

        mawindow = 15 #5
        corrintprofile = np.divide(MinRayfit, v2)
        corrintprofile = moving_average(corrintprofile, mawindow)

        profbuffer = 100 #18
        #corrintprofile[np.abs(corrintprofile) > 1.5 ] = 0
        corrintprofile[0:(y1-1-profbuffer)] = corrintprofile[(y1-1-profbuffer)] #18 was 25
        corrintprofile[(y2+1+profbuffer):ih] = corrintprofile[(y2+1+profbuffer)]
        a=[corrintprofile[(y2+1+profbuffer)]]*(mawindow-1)
        corrintprofile=np.concatenate((corrintprofile,a))
        plt.plot(corrintprofile)
        plt.show()

    return fit, a, b, c, corrintprofile, disk_flat

def correct_bad_lines_and_geom(Disk, options, not_fake):
    global hdr, basefich

    iw = Disk.shape[1]
    ih = Disk.shape[0]
    img = Disk

    if (m7>0) and (m8>0):
        y1 = m7
        y2 = m8
    else:
        y1, y2 = detect_bord(img, axis=1,offset=0)  # MattC   # bords verticaux

    # detection de mauvaises lignes

    # somme de lignes projetées sur axe Y
    ysum = np.mean(img, 1)

    # ne considere que les lignes du disque avec marge de 15 lignes
    ymargin = 0  # MattC was 15
    ysum = ysum[y1 + ymargin:y2 - ymargin]

    # filtrage sur fenetre de 31 pixels, polynome ordre 3 (etait 101 avant)
    yc = savgol_filter(ysum, 31, 3)

    # divise le profil somme par le profil filtré pour avoir les hautes
    # frequences
    hcol = np.divide(ysum, yc)

    # met à zero les pixels dont l'intensité est inferieur à 1.03 (3%)
    hcol[abs(hcol - 1) <= 0.03] = 0

    # tableau de zero en debut et en fin pour completer le tableau du disque
    a = [0] * (y1 + ymargin)
    b = [0] * (ih - y2 + ymargin)
    hcol = np.concatenate((a, hcol, b))

    # creation du tableau d'indice des lignes a corriger
    l_col = np.where(hcol != 0)
    listcol = l_col[0]

    # correction de lignes par filtrage median 13 lignes, empririque
    img_copy = np.copy(img)
    for c in listcol:
        m = img[c - 7:c + 6, ]
        s = np.median(m, 0)
        img_copy[c - 1:c, ] = s
    '''
    if options['save_fit']:
        DiskHDU=fits.PrimaryHDU(img_copy,header=hdr)
        DiskHDU.writeto(basefich+'_corr.fits', overwrite='True')
    '''
    return img_copy

def correct_transversalium(img, flag_nobords, options, not_fake, theflattener, theflatdisk):  # MattC
    global hdr, ih, basefich
    frame = img

    newiw=img.shape[1]
    ih=img.shape[0]
    flag_nobords = False
    # on cherche la projection de la taille max du soleil en Y
    if (m7>0) and (m8>0): #MattC
        y1 = m7
        y2 = m8
    else:
        y1, y2 = detect_bord(frame, axis=1,offset=0)  # MattC

    #y1,y2=detect_bord(frame, axis=1,offset=0) #MattC comment out if using above lines
    #print ('flat ',y1,y2)
    # si mauvaise detection des bords en x alors on doit prendre toute l'image
    if flag_nobords:
        ydisk=np.median(img,1)
    else:
        #plt.hist(frame.ravel(),bins=1000,)
        #plt.show()
        #plt.hist(frame.ravel(),bins=1000,cumulative=True)
        #plt.show()
        ##seuil_bas=np.percentile(frame,25)
        seuil_haut=np.percentile(frame,97) 
        #print ('Seuils de flat: ',seuil_bas, seuil_haut)
        #print ('Seuils bas x: ',seuil_bas*4)
        #print ('Seuils haut x: ',seuil_haut*0.25)
        #myseuil=seuil_haut*0.2
        myseuil=seuil_haut*0.5
        # filtre le profil moyen en Y en ne prenant que le disque
        ydisk=np.empty(ih+1)
        for j in range(0,ih):
            temp=np.copy(frame[j,:])
            temp=temp[temp>myseuil]
            if len(temp)!=0:
                ydisk[j]=np.median(temp)
            else:
                ydisk[j]=1
    ##y1=y1
    ##y2=y2
    ToSpline= ydisk[y1:y2]

    Smoothed2=savgol_filter(ToSpline,301, 3) # window size, polynomial order

    ''' #origly commented out
    #best fit d'un polynome degre 4
    np_m=np.asarray(ToSpline)
    ym=np_m.T
    xm=np.arange(y2-y1)
    p=np.polyfit(xm,ym,4)
    
    #calcul des x colonnes pour les y lignes du polynome
    a=p[0]
    b=p[1]
    c=p[2]
    d=p[3]
    e=p[4]
    x = np.arange(y2-y1)
    Smoothed=a*x**4+b*x**3+c*x**2+d*x+e
    '''

    # divise le profil reel par son filtre ce qui nous donne le flat
    hf=np.divide(ToSpline,Smoothed2)
       
    # elimine possible artefact de bord
    hf=hf[5:-5]
    
    #reconstruit le tableau du pofil complet 
    a=[1]*(y1+5)
    b=[1]*(ih-y2+5)
    hf=np.concatenate((a,hf,b))
    
    
    ##Smoothed=np.concatenate((a,Smoothed,b))
    ToSpline=np.concatenate((a,ToSpline,b))
    Smoothed2=np.concatenate((a,Smoothed2,b))
    '''
    # MattC test
    corrintprofile = np.zeros((y2-y1,1))
    popt, __ = curve_fit(func2,range(y1,y2),ydisk,p0=np.ones(3,))  # np.ones sets poly
    v2 = func2(range(y1,y2),*popt) #.astype('int')

    mawindow = 15 #5
    corrintprofile = np.divide(ydisk, v2)
    corrintprofile = moving_average(corrintprofile, mawindow)

    profbuffer = 0 # 100 18
    #corrintprofile[np.abs(corrintprofile) > 1.5 ] = 0
    corrintprofile[0:(y1-1-profbuffer)] = corrintprofile[(y1-1-profbuffer)] #18 was 25
    corrintprofile[(y2+1+profbuffer):ih] = corrintprofile[(y2+1+profbuffer)]
    a=[corrintprofile[(y2+1+profbuffer)]]*(mawindow-1)
    Smoothed2=np.concatenate((corrintprofile,a))
    # plt.plot(Smoothed2)
    # plt.show()
    # MattC end test
    '''
    # genere tableau image de flat 
    flat=[]

    hf = np.array(hf) / max(0.9, min(hf)) # don't make things bigger
    hf[hf==0] = 1

    for i in range(0,newiw):
        flat.append(hf)
    '''
    BelleImage=np.divide(frame,theflatdisk)  # MattC
    framef=np.array(BelleImage, dtype='uint16')
    DiskHDU = fits.PrimaryHDU(framef, header=hdr)
    DiskHDU.writeto(basefich + '_flat_MattC.fits', overwrite='True')
    '''
    np_flat=np.asarray(flat)
    flat = np_flat.T
    #print(hf, sum(hf)/len(hf), max(hf), min(hf))    
    # divise image par le flat
    BelleImage=np.divide(frame,flat)
    frame=np.array(BelleImage, dtype='uint16')

    # sauvegarde de l'image deflattée
    if options['save_fit'] and not_fake:
        DiskHDU = fits.PrimaryHDU(frame, header=hdr)
        DiskHDU.writeto(basefich + '_flat.fits', overwrite='True')
    return frame


def solex_proc(serfile, options):
    global hdr, ih, iw, basefich0, basefich
    clearlog()
    # plt.gray()              #palette de gris si utilise matplotlib pour visu
    # debug
    logme('Using pixel shift : ' + str(options['shift']))
    options['shift'] = [10, 0] + options['shift'] # 10, 0 are "fake"
    WorkDir = os.path.dirname(serfile) + "/"
    os.chdir(WorkDir)
    base = os.path.basename(serfile)
    basefich0 = os.path.splitext(base)[0]
    LineRecal = 1
    rdr = ser_reader(serfile)
    hdr = make_header(rdr)
    ih = rdr.ih
    iw = rdr.iw

    fit, a, b, c, flattener, flatdisk = compute_mean_return_fit(serfile, options, LineRecal)

    # Modification Jean-Francois: correct the variable names: A0, A1, A2
    logme('Coeff A0, A1, A2 :  ' + str(a) + '  ' + str(b) + '  ' + str(c))

    disk_list, ih, iw, FrameCount = read_video_improved(serfile, fit, LineRecal, options, flattener)

    hdr['NAXIS1'] = iw  # note: slightly dodgy, new width
    #if m1 > 0:
    #    options['shift'] = options['shift'] + [str(int(m4))+'FWHM']  #MattC adding bw 2*m4+1
    # sauve fichier disque reconstruit

    if options['flag_display']:
        cv2.destroyAllWindows()

    cercle = (-1, -1, -1)
    frames_circularized = []
    for i in range(len(disk_list)):
        basefich = basefich0 + '_shift=' + str(options['shift'][i])
        if options['save_fit'] and i >= 2:
            DiskHDU = fits.PrimaryHDU(disk_list[i], header=hdr)
            DiskHDU.writeto(basefich + '_img.fits', overwrite='True')

        """
        --------------------------------------------------------------------
        --------------------------------------------------------------------
        Badlines and geometry
        --------------------------------------------------------------------
        --------------------------------------------------------------------
        """
        try:
            img = correct_bad_lines_and_geom(disk_list[i], options, i >= 2)

            """
            --------------------------------------------------------------
            transversallium correction
            --------------------------------------------------------------
            """
            flag_nobords = False
            frame_flatted =  correct_transversalium(img, flag_nobords, options, i >= 2, flattener, flatdisk)
            # frame_flatted = img  # MattC
        except Exception:
            logme('WARNING: correct_bad_lines / correct_transversalium FAILED')
            frame_flatted = disk_list[i]

        """
        We now apply ellipse_fit to apply the geometric correction

        """
        # disk_list[0] is always shift = 10, for more contrast for ellipse fit
        if options['ratio_fixe'] is None and options['slant_fix'] is None:
            frame_circularized, cercle, options['ratio_fixe'], phi = ellipse_to_circle(
                frame_flatted, options, basefich)
            # in options angles are stored as degrees for some reason
            options['slant_fix'] = math.degrees(phi)
            frames_circularized.append(frame_circularized)
        else:
            ratio = options['ratio_fixe'] if not options['ratio_fixe'] is None else 1.0
            phi = math.radians(
                options['slant_fix']) if not options['slant_fix'] is None else 0.0
            frames_circularized.append(correct_image(frame_flatted / 65536, phi, ratio, np.array(
                [-1.0, -1.0]), -1.0, print_log=i == 0)[0])  # Note that we assume 16-bit

        # sauvegarde en fits de l'image finale

        if options['save_fit'] and i >= 2:  # first two shifts are not user specified
            DiskHDU = fits.PrimaryHDU(frames_circularized[-1], header=hdr)
            DiskHDU.writeto(basefich + '_recon.fits', overwrite='True')

    with open(basefich0 + '_log.txt', "w") as logfile:
        logfile.writelines(mylog)

    return frames_circularized[2:], hdr, cercle
