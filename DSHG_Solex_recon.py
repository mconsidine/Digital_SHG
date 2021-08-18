# -*- coding: utf-8 -*-
"""
@author: Valerie Desnoux
with improvements by Andrew Smith
contributors: Jean-Francois Pittet, Jean-Baptiste Butet, Pascal Berteau, Matt Considine
Version 13 August 2021

------------------------------------------------------------------------
reconstruction on an image from the deviations between the minimum of the line and a reference line

calcul sur une image des ecarts simples entre min de la raie et une ligne de reference
-------------------------------------------------------------------------

"""

from DSHG_solex_util import *
from DSHG_ser_read_video import *
from DSHG_ellipse_to_circle_JFP import ellipse_to_circle, correct_image
#import numpy_indexed as npi #MattC
from scipy import signal, fftpack #MattC
import astropy.units as u #MattC
from astropy.time import Time #MattC
from astropy.coordinates import SkyCoord, EarthLocation, get_sun, AltAz #MattC
import xlsxwriter #MattC

# read video and return constructed image of sun using fit and LineRecal
'''
def read_video_improved(serfile, fit, LineRecal, options):
    rdr = ser_reader(serfile)
    ih, iw = rdr.ih, rdr.iw
    
    if options['flag_display']:
        cv2.namedWindow('disk', cv2.WINDOW_NORMAL)
        FrameMax=rdr.FrameCount
        cv2.resizeWindow('disk', FrameMax//3, ih//3)
        cv2.moveWindow('disk', 200, 0)
        #initialize le tableau qui va recevoir la raie spectrale de chaque trame
        Disk=np.zeros((ih,FrameMax), dtype=rdr.infiledatatype) #MattC
        
        cv2.namedWindow('image', cv2.WINDOW_NORMAL)
        cv2.moveWindow('image', 0, 0)
        cv2.resizeWindow('image', int(iw), int(ih))
    else:
        #Disk=np.zeros((ih,1), dtype='uint16')
        FrameMax=rdr.FrameCount
        Disk=np.zeros((ih,FrameMax), dtype=rdr.infiledatatype) #MattC
        
    shift = options['shift']
    ind_l = (np.asarray(fit)[:, 0] + np.ones(ih) * (LineRecal + shift )).astype(int) #MattC

    #CLEAN if fitting goes too far
    ind_l[ind_l < 0] = 0
    ind_l[ind_l > iw - 2] = iw - 2
    #print('ind l ',ind_l[0:49]) #MattC
    #print('ind l shape ',ind_l.shape)    #MattC
    ind_r = (ind_l + np.ones(ih)).astype(int)
    #print('ind r ',ind_r[0:49]) #MattC
    #print('ind r shape ',ind_r.shape) #MattC
    left_weights = np.ones(ih) - np.asarray(fit)[:, 1]
    right_weights = np.ones(ih) - left_weights
    # lance la reconstruction du disk a partir des trames
    #print('L Weights : ', left_weights.shape) #MattC
    #print('R Weights : ', right_weights.shape) #MattC
    #print('reader num frames:', rdr.FrameCount) #MattC
    ##an_array = np.arange(0,iw)
    ##repetitions = ih
    ##testmask = np.tile(an_array, (repetitions))#.reshape(ih,iw)
    ##print('testmask shape',testmask.shape)
    col_mask = np.tile(range(0,options['pixel_bandwidth']),ih).reshape(ih,options['pixel_bandwidth'])  #MattC
    #print('col mask shape',col_mask.shape) #MattC
    #print('col mask data ',col_mask[0:49,:])
 
    while rdr.has_frames():
        img = rdr.next_frame()               
        #img2 = img.reshape(ih*iw)
        if options['flag_display'] and rdr.FrameIndex % 10 == 0 :
            cv2.imshow('image', img)
            if cv2.waitKey(1)==27:
                cv2.destroyAllWindows()
                sys.exit()

        if options['pixel_bandwidth'] >= 0: #MattC TODO: need to test this against prior default
            left_banda = np.repeat(ind_l, options['pixel_bandwidth'],axis=0).reshape(ih,options['pixel_bandwidth'])
            left_band = left_banda - col_mask
            mx1 = np.take_along_axis(img,left_band.reshape(ih,options['pixel_bandwidth']),1)
            left_col = np.mean(mx1,axis=1) 

            right_banda = np.repeat(ind_r, options['pixel_bandwidth'],axis=0).reshape(ih,options['pixel_bandwidth'])
            right_band = right_banda + col_mask
            mx2 = np.take_along_axis(img,right_band.reshape(ih,options['pixel_bandwidth']),1)
            right_col = np.mean(mx2,axis=1)              
            
            #
            #if rdr.FrameIndex == 2000: #MattC
            #    print('FRAMEINDEX ',rdr.FrameIndex)
            #    #print('mx1 shape ',mx1.shape)            
            #    #print('mx2 shape ',mx2.shape)    
            #    #print('mx1 type ',type(mx1))
            #    #print('mx2 type ',type(mx2)) 
            #    print('img data ',img[0:49,:])
            #    print('left banda data ',left_banda[0:49,:])
            #    print('left band data ',left_band[0:49,:])
            #    print('mx1 data ',mx1[0:49,:])          
            #    print('left col data ',left_col[0:49])
            #    print('right banda data ',right_banda[0:49,:])
            #    print('right band data ',right_band[0:49,:])                 
            #    print('mx2 data ',mx2[0:49,:])
            #    print('right col data ',right_col[0:49])                        
            #
            mx3 = np.hstack((mx1,mx2))
            #
            #if rdr.FrameIndex == 2000: #MattC
            #    print('mx3 shape ',mx3.shape)
            #    print('mx3 data ',mx3[0:49,:])                 
            #
            IntensiteRaie = np.mean(mx3, axis=1) #MattC could be any function, techically
            #
            #if rdr.FrameIndex == 2000: #MattC
            #    print('Intens shape ',IntensiteRaie.shape)
            #    print('Intens data ',IntensiteRaie[0:49].astype(rdr.infiledatatype))
            #
        #else:
        #    left_col = img[np.arange(ih), ind_l]
        #    right_col = img[np.arange(ih), ind_r]

        #    IntensiteRaie = left_col*left_weights + right_col*right_weights
        
        #ajoute au tableau disk 
        Disk[:,rdr.FrameIndex]=IntensiteRaie.astype(rdr.infiledatatype) #MattC
        
        if options['flag_display'] and rdr.FrameIndex % 10 ==0:
            cv2.imshow ('disk', Disk)
            if cv2.waitKey(1) == 27:                     # exit if Escape is hit
                cv2.destroyAllWindows()    
                sys.exit()
    if rdr.infiledatatype == 'uint8': #MattC deal with AVI scale
        Disk = (Disk*(65535/255))
    return Disk, ih, iw, rdr.FrameCount
'''
def read_video_improved(serfile, fit, LineRecal, options):
    rdr = ser_reader(serfile)
    ih, iw = rdr.ih, rdr.iw
    FrameMax=rdr.FrameCount
    disk_list = [np.zeros((ih,FrameMax), dtype='uint16') for _ in options['shift']]

    
    if options['flag_display']:
        cv2.namedWindow('disk', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('disk', FrameMax//3, ih//3)
        cv2.moveWindow('disk', 200, 0)
        #initialize le tableau qui va recevoir la raie spectrale de chaque trame
        #Disk=np.zeros((ih,FrameMax), dtype=rdr.infiledatatype) #MattC
        
        cv2.namedWindow('image', cv2.WINDOW_NORMAL)
        cv2.moveWindow('image', 0, 0)
        cv2.resizeWindow('image', int(iw), int(ih))
    #else: #MattC
    #    #Disk=np.zeros((ih,1), dtype='uint16')
    #    FrameMax=rdr.FrameCount
    #    Disk=np.zeros((ih,FrameMax), dtype=rdr.infiledatatype) #MattC

    col_indeces = []
          
    #shift = options['shift']
    #ind_l = (np.asarray(fit)[:, 0] + np.ones(ih) * (LineRecal + shift )).astype(int) #MattC

    for shift in options['shift']:
        ind_l = (np.asarray(fit)[:, 0] + np.ones(ih) * (LineRecal + shift)).astype(int)
        
        #CLEAN if fitting goes too far
        ind_l[ind_l < 0] = 0
        ind_l[ind_l > iw - 2] = iw - 2
        ind_r = (ind_l + np.ones(ih)).astype(int)
        col_indeces.append((ind_l, ind_r))
        
    left_weights = np.ones(ih) - np.asarray(fit)[:, 1]
    right_weights = np.ones(ih) - left_weights

    # lance la reconstruction du disk a partir des trames
    print('reader num frames:', rdr.FrameCount)
    
    col_mask = np.tile(range(0,options['pixel_bandwidth']),ih).reshape(ih,options['pixel_bandwidth'])  #MattC
    if rdr.scalemax == 255: #MattC
        ascalefactor = 255
        options['scalemax'] = 255
    else:
        ascalefactor = 1
        options['scalemax'] = 64000
        
    while rdr.has_frames():
        img = rdr.next_frame()               
        if options['flag_display'] and rdr.FrameIndex % 10 == 0 :
            cv2.imshow('image', img)
            if cv2.waitKey(1)==27:
                cv2.destroyAllWindows()
                sys.exit()

        for i in range(len(options['shift'])):
            ind_l, ind_r = col_indeces[i]
            if options['pixel_bandwidth'] > 1: #MattC TODO: need to test this against prior default
                left_banda = np.repeat(ind_l, options['pixel_bandwidth'],axis=0).reshape(ih,options['pixel_bandwidth'])
                left_band = left_banda - col_mask
                mx1 = np.take_along_axis(img,left_band.reshape(ih,options['pixel_bandwidth']),1)
                left_col = np.mean(mx1,axis=1) 

                right_banda = np.repeat(ind_r, options['pixel_bandwidth'],axis=0).reshape(ih,options['pixel_bandwidth'])
                right_band = right_banda + col_mask
                mx2 = np.take_along_axis(img,right_band.reshape(ih,options['pixel_bandwidth']),1)
                right_col = np.mean(mx2,axis=1)              

                mx3 = np.hstack((mx1,mx2))

                IntensiteRaie = np.mean(mx3, axis=1) #MattC could be any function, techically

            else:
                left_col = img[np.arange(ih), ind_l]
                right_col = img[np.arange(ih), ind_r]
                IntensiteRaie = left_col*left_weights + right_col*right_weights
        
            #ajoute au tableau disk 
            #Disk[:,rdr.FrameIndex]=IntensiteRaie.astype(rdr.infiledatatype) #MattC
            #if options['hflip']:
            #    disk_list[i][:,rdr.FrameCount-rdr.FrameIndex-1]=IntensiteRaie.astype(rdr.infiledatatype) #MattC
            #else:
            if 1:
                disk_list[i][:,rdr.FrameIndex]=IntensiteRaie.astype(rdr.infiledatatype) #MattC
                    
        if options['flag_display'] and rdr.FrameIndex % 10 ==0:
            cv2.imshow ('disk', disk_list[1]*ascalefactor) # disk_list[1] is always shift = 0 #MattC
            if cv2.waitKey(1) == 27:                     # exit if Escape is hit
                cv2.destroyAllWindows()    
                sys.exit()
    if rdr.infiledatatype == 'uint8': #MattC deal with AVI scale
        for i in range(len(options['shift'])):
            disk_list[i] = (disk_list[i]*(64000/255)).astype('uint16')
    return disk_list, ih, iw, rdr.FrameCount

def make_header(rdr):        
    # initialisation d'une entete fits (etait utilisé pour sauver les trames individuelles)
    hdr= fits.Header()
    hdr['SIMPLE']='T'
    hdr['BITPIX']=32
    hdr['NAXIS']=2
    hdr['NAXIS1'] = rdr.iw
    hdr['NAXIS2'] = rdr.ih
    hdr['BZERO']=0
    hdr['BSCALE']=1
    hdr['BIN1']=1
    hdr['BIN2']=1
    hdr['EXPTIME']=0
    return hdr

# compute mean image of video
def compute_mean(serfile, options): #MattC
    """IN : serfile path"
    OUT :numpy array
    """
    rdr = ser_reader(serfile)
    logme('Width, Height : '+str(rdr.Width)+' '+str(rdr.Height)) 
    logme('Number of frames : '+str(rdr.FrameCount))
    my_data = np.zeros((rdr.ih, rdr.iw),dtype='uint64')
        
    if options['phasecorr'] == True: #MattC  I'm pretty sure this isn't even close to working.  Shifts seem reasonable
        #but I haven't been successful in reflecting them in the image.  Getting tripped up by ... something.
        my_shift_data = np.zeros((rdr.FrameCount),dtype='int16') #NOTE: need a signed int here
        testcount = 0    
        while rdr.has_frames():
            img = rdr.next_frame()
            #phase correlation for shifts: https://stackoverflow.com/questions/4688715/find-time-shift-between-two-similar-waveforms
            if testcount == 0: #MattC https://stackoverflow.com/questions/35567906/how-to-apply-phase-correlation-in-1d-signal
                #establish the slice from which movement is judged
                
                baseimg = img[:,rdr.iw//2]                
                fft_sig1 = fftpack.fft(img[:,rdr.iw//2])
                Ar = -fft_sig1.conjugate()
                #fft_sig1 = np.fft.fft(img[:,rdr.iw//2])                
                my_shift_data[testcount] = 0
                shift_img = img
            else:
                '''
                fft_sig2 = np.fft.fft(img[:,rdr.iw//2])
                fft_sig2_conj = np.conj(fft_sig2)
                R = (fft_sig1 * fft_sig2_conj) / abs(fft_sig1 * fft_sig2_conj)
                r = np.fft.ifft(R)
                pix_shift = np.argmax(r)
                if testcount % 500 == 0:
                    print('pixel shift = ', testcount, pix_shift)
                my_shift_data[testcount] = pix_shift
                '''
                
                fft_sig2 = fftpack.fft(img[:,rdr.iw//2])
                Br = -fft_sig2.conjugate()
                pix_shift = np.argmax(np.abs(fftpack.ifft(Ar*fft_sig2)))
                print(pix_shift)
                #pix_shift2 = np.argmax(np.abs(fftpack.ifft(fft_sig1*Br)))
                my_shift_data[testcount] = pix_shift    
                           
                #NOTE: I have no idea yet if this is really working; shift data still needs to go to Disk

                shift_img = np.zeros_like(img) #We're going to shift the entire image
                shift_img = np.roll(img,-pix_shift,axis=0) #shift up/down #TODO: preserve rolled rows
                
                fft_sig1 = fft_sig2 #fftpack.fft(shift_img[:,rdr.iw//2]) #????
                Ar = -fft_sig1.conjugate()
                
            testcount += 1 #MattC               
            my_data += shift_img #MattC shift data should be passed back so that smile calc/fit can accomodate    
    else:
        while rdr.has_frames():
            img = rdr.next_frame()
            my_data += img
    #workbook = xlsxwriter.Workbook('arrays.xlsx')
    #worksheet = workbook.add_worksheet()
    #row = 0
    #for col, data in enumerate(my_shift_data.T):
    #    worksheet.write_column(row, col, data)
    #workbook.close()   

    return (my_data / rdr.FrameCount).astype(rdr.infiledatatype) #MattC

def compute_mean_return_fit(serfile, options, LineRecal = 1): 
    global hdr, ih, iw

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
    mean_img= compute_mean(serfile, options) #MattC

    """
    ----------------------------------------------------------------------------
    Calcul polynome ecart sur une image au centre de la sequence
    ----------------------------------------------------------------------------
    """

    
    #savefich=basefich+'_mean'

    if options['save_fit']:
        DiskHDU=fits.PrimaryHDU(mean_img,header=hdr)
        DiskHDU.writeto(basefich0+'_mean.fits', overwrite='True')


    #affiche image moyenne
    if flag_display:
        cv2.namedWindow('Video mean', cv2.WINDOW_NORMAL) #MattC
        cv2.resizeWindow('Video mean', iw, ih) #MattC
        cv2.moveWindow('Video mean', 100, 0)  #MattC
        cv2.imshow ('Video mean', mean_img)  #MattC
        if cv2.waitKey(2000) == 27:                     # exit if Escape is hit
            cv2.destroyAllWindows()
            sys.exit()
        
        cv2.destroyAllWindows()
    
    y1,y2=detect_bord(mean_img, axis=1, offset=5)
    logme('Vertical limits y1, y2 : '+str(y1)+' '+str(y2))
    
    PosRaieHaut=y1
    PosRaieBas=y2
    
    """
    -----------------------------------------------------------
    Trouve les min intensité de la raie
    -----------------------------------------------------------
    """
    # construit le tableau des min de la raie a partir du haut jusqu'en bas
    MinOfRaie=[]
    
    for i in range(PosRaieHaut,PosRaieBas):
        line_h=mean_img[i,:]
        MinX=line_h.argmin()
        MinOfRaie.append([MinX,i])
        #print('MinOfRaie x,y', MinX,i)
    
    #best fit d'un polynome degre 2, les lignes y sont les x et les colonnes x sont les y
    np_m=np.asarray(MinOfRaie)
    xm,ym=np_m.T
    #LineRecal=xm.min()
    
    p=np.polyfit(ym,xm,2)
    
    #calcul des x colonnes pour les y lignes du polynome
    a=p[0]
    b=p[1]
    c=p[2]
    fit=[]
    #ecart=[]
    for y in range(0,ih):
        x=a*y**2+b*y+c
        deci=x-int(x)
        fit.append([int(x)-LineRecal,deci,y])
    return fit, a, b, c

def correct_bad_lines_and_geom(Disk, options, not_fake):
    global hdr, basefich
    
    iw=Disk.shape[1]
    ih=Disk.shape[0]
    img=Disk
    
    y1,y2=detect_bord (img, axis=1,offset=5)    # bords verticaux
    
    #detection de mauvaises lignes
    
    # somme de lignes projetées sur axe Y
    ysum=np.mean(img,1)
    #plt.plot(ysum)
    #plt.show()
    # ne considere que les lignes du disque avec marge de 15 lignes 
    ysum=ysum[y1+15:y2-15]
    
    # filtrage sur fenetre de 31 pixels, polynome ordre 3 (etait 101 avant)
    yc=savgol_filter(ysum,31, 3)

    # divise le profil somme par le profil filtré pour avoir les hautes frequences
    hcol=np.divide(ysum,yc)

    # met à zero les pixels dont l'intensité est inferieur à 1.03 (3%)
    hcol[abs(hcol-1)<=0.03]=0

    
    # tableau de zero en debut et en fin pour completer le tableau du disque
    a=[0]*(y1+15)
    b=[0]*(ih-y2+15)
    hcol=np.concatenate((a,hcol,b))
    #plt.plot(hcol)
    #plt.show()
    
    # creation du tableau d'indice des lignes a corriger
    l_col=np.where(hcol!=0)
    listcol=l_col[0]
    

    # correction de lignes par filtrage median 13 lignes, empririque
    img_copy = np.copy(img)
    for c in listcol:
        m=img[c-7:c+6,]
        s=np.median(m,0)
        img_copy[c-1:c,]=s
    
    #sauvegarde le fits

    if options['save_fit'] and not_fake:
        DiskHDU=fits.PrimaryHDU(img_copy,header=hdr)
        DiskHDU.writeto(basefich+'_corr.fits', overwrite='True')
        
        
    return img_copy

def correct_transversalium(img, flag_nobords, options, not_fake):
    global hdr, ih, basefich
    frame = img
    newiw=img.shape[1]
    ih=img.shape[0]
    flag_nobords = False
    # on cherche la projection de la taille max du soleil en Y
    y1,y2=detect_bord(frame, axis=1,offset=0)
    #print ('flat ',y1,y2)
    # si mauvaise detection des bords en x alors on doit prendre toute l'image
    if flag_nobords:
        ydisk=np.median(img,1)
    else:
        #plt.hist(frame.ravel(),bins=1000,)
        #plt.show()
        #plt.hist(frame.ravel(),bins=1000,cumulative=True)
       # plt.show()
        seuil_bas=np.percentile(frame,25)
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
    y1=y1
    y2=y2
    ToSpline= ydisk[y1:y2]
 
    
    Smoothed2=savgol_filter(ToSpline,301, 3) # window size, polynomial order
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
    '''
    Smoothed=[]
    for x in range(0,y2-y1):
        y=a*x**4+b*x**3+c*x**2+d*x+e
        Smoothed.append(y)
    '''
    x = np.arange(y2-y1)
    Smoothed=a*x**4+b*x**3+c*x**2+d*x+e
    
    # divise le profil reel par son filtre ce qui nous donne le flat
    hf=np.divide(ToSpline,Smoothed2)
       
    # elimine possible artefact de bord
    hf=hf[5:-5]
    
    #reconstruit le tableau du pofil complet 
    a=[1]*(y1+5)
    b=[1]*(ih-y2+5)
    hf=np.concatenate((a,hf,b))
    
    
    Smoothed=np.concatenate((a,Smoothed,b))
    ToSpline=np.concatenate((a,ToSpline,b))
    Smoothed2=np.concatenate((a,Smoothed2,b))

    
    # genere tableau image de flat 
    flat=[]
    hf = np.array(hf) / max(0.9, min(hf)) # don't make things bigger
    hf[hf==0] = 1
    for i in range(0,newiw):
        flat.append(hf)
        
    np_flat=np.asarray(flat)
    flat = np_flat.T
    #print(hf, sum(hf)/len(hf), max(hf), min(hf))    
    # divise image par le flat
    BelleImage=np.divide(frame,flat)
    frame=np.array(BelleImage, dtype='uint16')
    # sauvegarde de l'image deflattée
    if options['save_fit'] and not_fake:
        DiskHDU=fits.PrimaryHDU(frame,header=hdr)
        DiskHDU.writeto(basefich+'_flat.fits', overwrite='True')
    return frame

def solex_proc(serfile, options):
    global hdr, ih, iw, basefich0, basefich
    clearlog()
    #plt.gray()              #palette de gris si utilise matplotlib pour visu debug
    logme('Using pixel shift : ' + str(options['shift']))
    options['shift'] = [10, 0] + options['shift'] # 10, 0 are "fake"
    WorkDir=os.path.dirname(serfile)+"/"
    os.chdir(WorkDir)
    base=os.path.basename(serfile)
    basefich0=os.path.splitext(base)[0]
    LineRecal=1
    rdr = ser_reader(serfile)
    hdr = make_header(rdr)
    ih = rdr.ih
    iw = rdr.iw
    
    fit, a, b, c = compute_mean_return_fit(serfile, options, LineRecal)
    
    # Modification Jean-Francois: correct the variable names: A0, A1, A2
    logme('Coeff A0, A1, A2 :  '+str(a)+'  '+str(b)+'  '+str(c))
    
    disk_list, ih, iw, FrameCount = read_video_improved(serfile, fit, LineRecal, options)

    hdr['NAXIS1']=iw # note: slightly dodgy, new width
   
    #sauve fichier disque reconstruit

    '''
    Start of angle/shear test ... MattC
    '''
    '''
    print(rdr.DTime)
    print(rdr.DTimeUTC)
    
    the_location = EarthLocation.from_geodetic(lat=45.2*u.deg, lon=-70*u.deg, height=330*u.m)
    utcoffset = -4*u.hour
    
    startobstime = Time('2021-07-22 10:58:23') - utcoffset
    timeframe = AltAz(obstime=startobstime, location=the_location)
    startsunaltaz = get_sun(startobstime).transform_to(timeframe)
    print("Sun's altitude = {0.alt:.8}".format(startsunaltaz))
    print("Sun's alzimuth = {0.az:.8}".format(startsunaltaz))

    endobstime = startobstime + 2*u.min
    timeframe = AltAz(obstime=endobstime, location=the_location)
    endsunaltaz = get_sun(endobstime).transform_to(timeframe)
    print("Sun's altitude = {0.alt:.8}".format(endsunaltaz))
    print("Sun's azimuth = {0.az:.8}".format(endsunaltaz))

    print(endsunaltaz.alt-startsunaltaz.alt)
    print(endsunaltaz.az-startsunaltaz.az)
    the_slope = ((endsunaltaz.alt-startsunaltaz.alt)/(endsunaltaz.az-startsunaltaz.az))/2
    
    slope_rad=math.atan(the_slope)
    slope_deg=math.degrees(slope_rad)
    print(the_slope)
    print(slope_deg)

    print("strip width in pixels : ",disk_list[0].shape[1])
    print("strip height in pixels : ",disk_list[0].shape[0])
    
    plt.imshow(disk_list[0])
    plt.show()        

    rows, cols = disk_list[0].shape    
    angle_shs = 2*slope_deg #rotation about center; undo shear then add angle
    shearx_shs = 0*the_slope #positive leans to left; anchors on bottom
    sheary_shs = the_slope #positive skews down from left; anchors on left
    transl_shs = 0
    type_border = cv2.BORDER_CONSTANT
    color_border = (255,255,255)
    
    ###Calc space needed for rotation
    ##M=cv2.getRotationMatrix2D((cols/2,rows/2), 45, 1) #45 would be max
    ##cos_part = np.abs(M[0,0])
    ##sin_part = np.abs(M[0,1])
    ##new_cols = int((rows * sin_part) + (cols * cos_part))
    ##new_rows = int((rows * cos_part) + (cols * sin_part))
    ###Calc space needed for shear
    ##new_cols = new_cols + (shearx_shs*new_cols)
    ##new_rows = new_rows + (sheary_shs*new_rows)
    
    new_rows = max(rows,cols)
    new_cols = max(rows,cols)+2*rows
    #Calc space needed for border
    up_down = int((new_rows-rows)/2)
    left_right = int((new_cols-cols)/2)

    sheared_Disk = cv2.copyMakeBorder(disk_list[0], up_down, up_down, left_right, left_right, type_border, value=color_border)
    
    plt.axis('off')
    plt.imshow(sheared_Disk)
    plt.show()
    
    rows, cols = sheared_Disk.shape
    #Apply transform
    M_rot = cv2.getRotationMatrix2D((cols/2, rows/2),angle_shs,1)
   
    translat_center_x = 0
    translat_center_y = 0
    translat_center_x = -(shearx_shs*cols)/2
    translat_center_y = -(sheary_shs*rows)/2
    
    M = M_rot + np.float64([[0,shearx_shs,transl_shs + translat_center_x], [sheary_shs,0,transl_shs + translat_center_y]])

    sheared_Disk = cv2.warpAffine(sheared_Disk, M, (cols,rows), borderMode=type_border, borderValue=color_border)
    
    sheared_Disk = cv2.flip(sheared_Disk, 1)
    
    plt.axis('off')
    plt.imshow(sheared_Disk)
    plt.show()
    #options['slant_fix']=slope_deg
    # ... end of angle test MattC   
    '''
    
    if options['flag_display']:
        cv2.destroyAllWindows()

    cercle = (-1, -1, -1)
    frames_circularized = []
    for i in range(len(disk_list)):
        basefich = basefich0 + '_hbw='+str(options['pixel_bandwidth']) + '_shift='+str(options['shift'][i]) #MattC

        if options['save_fit'] and i >= 2:
            DiskHDU=fits.PrimaryHDU(disk_list[i],header=hdr)
            DiskHDU.writeto(basefich+'_img.fits', overwrite='True')
    
        """
        --------------------------------------------------------------------
        --------------------------------------------------------------------
        Badlines and geometry
        --------------------------------------------------------------------
        --------------------------------------------------------------------
        """
        img = correct_bad_lines_and_geom(disk_list[i], options, i >= 2)
            
        """
        --------------------------------------------------------------
        transversallium correction
        --------------------------------------------------------------
        """
        flag_nobords = False
        frame_flatted = correct_transversalium(img,flag_nobords, options, i >= 2)

        """
        We now apply ellipse_fit to apply the geometric correction

        """
        # disk_list[0] is always shift = 10, for more contrast for ellipse fit
        if options['ratio_fixe'] is None and options['slant_fix'] is None:
            frame_circularized, cercle, options['ratio_fixe'], phi = ellipse_to_circle(frame_flatted, options)
            options['slant_fix'] = math.degrees(phi) # in options angles are stored as degrees for some reason
            frames_circularized.append(frame_circularized)
        else:
            ratio = options['ratio_fixe']             if not options['ratio_fixe'] is None else 1.0
            phi = math.radians(options['slant_fix'])  if not options['slant_fix'] is None else 0.0
            frames_circularized.append(correct_image(frame_flatted / 65536, phi, ratio, np.array([-1.0, -1.0]), print_log = i == 0)[0]) # Note that we assume 16-bit

        # sauvegarde en fits de l'image finale
        
        if options['save_fit'] and i >= 2: # first two shifts are not user specified
            DiskHDU=fits.PrimaryHDU(frames_circularized[-1],header=hdr)
            DiskHDU.writeto(basefich + '_recon.fits', overwrite='True')
            
    with  open(basefich0+'_log.txt', "w") as logfile:
        logfile.writelines(mylog)
    
    return frames_circularized[2:], hdr, cercle
