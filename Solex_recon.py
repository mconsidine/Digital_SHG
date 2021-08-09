# -*- coding: utf-8 -*-
"""
@author: Valerie Desnoux
with improvements by Andrew Smith
contributors: Jean-Francois Pittet, Jean-Baptiste Butet, Pascal Berteau, Matt Considine
Version 1 August 2021

------------------------------------------------------------------------
reconstruction on an image from the deviations between the minimum of the line and a reference line

calcul sur une image des ecarts simples entre min de la raie et une ligne de reference
-------------------------------------------------------------------------

"""

from solex_util import *
from ser_read_video import *
from ellipse_to_circle import ellipse_to_circle, correct_image
import numpy_indexed as npi

# read video and return constructed image of sun using fit and LineRecal
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
            '''
            if rdr.FrameIndex == 2000: #MattC
                print('FRAMEINDEX ',rdr.FrameIndex)
                #print('mx1 shape ',mx1.shape)            
                #print('mx2 shape ',mx2.shape)    
                #print('mx1 type ',type(mx1))
                #print('mx2 type ',type(mx2)) 
                print('img data ',img[0:49,:])
                print('left banda data ',left_banda[0:49,:])
                print('left band data ',left_band[0:49,:])
                print('mx1 data ',mx1[0:49,:])          
                print('left col data ',left_col[0:49])
                print('right banda data ',right_banda[0:49,:])
                print('right band data ',right_band[0:49,:])                 
                print('mx2 data ',mx2[0:49,:])
                print('right col data ',right_col[0:49])                        
            '''
            mx3 = np.hstack((mx1,mx2))
            '''
            if rdr.FrameIndex == 2000: #MattC
                print('mx3 shape ',mx3.shape)
                print('mx3 data ',mx3[0:49,:])                 
             '''
            IntensiteRaie = np.mean(mx3, axis=1) #MattC could be any function, techically
            '''
            if rdr.FrameIndex == 2000: #MattC
                print('Intens shape ',IntensiteRaie.shape)
                print('Intens data ',IntensiteRaie[0:49].astype(rdr.infiledatatype))
            '''
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


def make_header(rdr):        
    # initialisation d'une entete fits (etait utilisé pour sauver les trames individuelles
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
def compute_mean(serfile):
    """IN : serfile path"
    OUT :numpy array
    """
    rdr = ser_reader(serfile)
    logme('Width, Height : '+str(rdr.Width)+' '+str(rdr.Height)) 
    logme('Number of frames : '+str(rdr.FrameCount))
    my_data = np.zeros((rdr.ih, rdr.iw),dtype='uint64')
    while rdr.has_frames():
        img = rdr.next_frame()
        my_data += img
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
    mean_img= compute_mean(serfile)
    
    """
    ----------------------------------------------------------------------------
    Calcul polynome ecart sur une image au centre de la sequence
    ----------------------------------------------------------------------------
    """

    
    #savefich=basefich+'_mean'

    if options['save_fit']:
        DiskHDU=fits.PrimaryHDU(mean_img,header=hdr)
        DiskHDU.writeto(basefich+'_mean.fits', overwrite='True')


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

def correct_bad_lines_and_geom(Disk, options):
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
    for c in listcol:
        m=img[c-7:c+6,]
        s=np.median(m,0)
        img[c-1:c,]=s
    
    #sauvegarde le fits

    if options['save_fit']:
        DiskHDU=fits.PrimaryHDU(img,header=hdr)
        DiskHDU.writeto(basefich+'_corr.fits', overwrite='True')
        
        
    return img

def correct_transversalium(img, flag_nobords, options):
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
    Smoothed=[]
    for x in range(0,y2-y1):
        y=a*x**4+b*x**3+c*x**2+d*x+e
        Smoothed.append(y)


    
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
    for i in range(0,newiw):
        flat.append(hf)
        
    np_flat=np.asarray(flat)
    flat = np_flat.T
    #evite les divisions par zeros...
    flat[flat==0]=1
    
    # divise image par le flat
    BelleImage=np.divide(frame,flat)
    frame=np.array(BelleImage, dtype='uint16')
    # sauvegarde de l'image deflattée
    if options['save_fit']:
        DiskHDU=fits.PrimaryHDU(frame,header=hdr)
        DiskHDU.writeto(basefich+'_flat.fits', overwrite='True')
    return frame

def solex_proc(serfile, options):
    global hdr, ih, iw, basefich
    clearlog()
    #plt.gray()              #palette de gris si utilise matplotlib pour visu debug
    logme('Using pixel shift : ' + str(options['shift']))
    WorkDir=os.path.dirname(serfile)+"/"
    os.chdir(WorkDir)
    base=os.path.basename(serfile)
    basefich=os.path.splitext(base)[0]
    LineRecal=1
    rdr = ser_reader(serfile)
    hdr = make_header(rdr)
    ih = rdr.ih
    iw = rdr.iw
    
    fit, a, b, c = compute_mean_return_fit(serfile, options, LineRecal)
    
    # Modification Jean-Francois: correct the variable names: A0, A1, A2
    logme('Coeff A0, A1, A2 :  '+str(a)+'  '+str(b)+'  '+str(c))
    
    Disk, ih, iw, FrameCount = read_video_improved(serfile, fit, LineRecal, options)

    hdr['NAXIS1']=iw # note: slightly dodgy, new width
   
    #sauve fichier disque reconstruit


    if options['save_fit']:
        DiskHDU=fits.PrimaryHDU(Disk,header=hdr)
        DiskHDU.writeto(basefich+'_img.fits', overwrite='True')
    
    if options['flag_display']:
        cv2.destroyAllWindows()
    
    """
    --------------------------------------------------------------------
    --------------------------------------------------------------------
    Badlines and geometry
    --------------------------------------------------------------------
    --------------------------------------------------------------------
    """
    img = correct_bad_lines_and_geom(Disk, options)
        
    """
    --------------------------------------------------------------
    transversallium correction
    --------------------------------------------------------------
    """
    flag_nobords = False
    frame_flatted = correct_transversalium(img,flag_nobords, options)

    """
    We now apply ellipse_fit to apply the geometric correction

    """
    if options['ratio_fixe'] is None and options['slant_fix'] is None:
        frame_circularized, cercle = ellipse_to_circle(frame_flatted, options)
    else:
        ratio = options['ratio_fixe']             if not options['ratio_fixe'] is None else 1.0
        phi = math.radians(options['slant_fix'])  if not options['slant_fix'] is None else 0.0
        frame_circularized, cercle = correct_image(frame_flatted / 65536, phi, ratio, np.array([-1.0, -1.0]))[0], (-1, -1, -1) # Note that we assume 16-bit

    # sauvegarde en fits de l'image finale
    
    if options['save_fit']:
        DiskHDU=fits.PrimaryHDU(frame_circularized,header=hdr)
        DiskHDU.writeto(basefich+'_recon.fits', overwrite='True')
            
    with  open(basefich+'_log.txt', "w") as logfile:
        logfile.writelines(mylog)
    
    return frame_circularized, hdr, cercle
    


