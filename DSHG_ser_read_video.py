"""
@author: Valerie Desnoux
with improvements by Andrew Smith
Version 1 August 2021
additional changes by Matt Considine

"""
import numpy as np
import cv2 #MattC

class ser_reader:

    def __init__(self, serfile):
    
        self.serfile = serfile
        
        if (self.serfile.upper().endswith('.SER')==True): #MattC 20210726
            self.SER_flag=True
            self.AVI_flag=False
            self.infiledatatype='uint16'
        elif (self.serfile.upper().endswith('.AVI')==True):
            self.SER_flag=False
            self.AVI_flag=True
            self.infiledatatype='uint8'
        else:
            self.SER_flag=False
            self.AVI_flag=False
        
        #ouverture et lecture de l'entete du fichier ser

        if (self.SER_flag==True): #MattC
            self.FileID=np.fromfile(serfile, dtype='int8',count=14)
            offset=14

            self.LuID=np.fromfile(serfile, dtype=np.uint32, count=1, offset=offset)
            offset=offset+4
        
            self.ColorID=np.fromfile(serfile, dtype='uint32', count=1, offset=offset)
            offset=offset+4
        
            self.littleEndian=np.fromfile(serfile, dtype='uint32', count=1,offset=offset)
            offset=offset+4
        
            self.Width=np.fromfile(serfile, dtype='uint32', count=1,offset=offset)[0]
            offset=offset+4
        
            self.Height=np.fromfile(serfile, dtype='uint32', count=1,offset=offset)[0]
            offset=offset+4
        
            PixelDepthPerPlane=np.fromfile(serfile, dtype='uint32', count=1,offset=offset)
            self.PixelDepthPerPlane=PixelDepthPerPlane[0]
            offset=offset+4
        
            FrameCount=np.fromfile(serfile, dtype='uint32', count=1,offset=offset)
            self.FrameCount=FrameCount[0]
        
            if self.PixelDepthPerPlane==8:
                self.infiledatatype='uint8'
                self.count=self.Width*self.Height       # Nombre d'octet d'une trame
                self.infilebytes=1
                self.scalemax=255
            else:
                self.infiledatatype='uint16'
                self.count=self.Width*self.Height      # Nombre d'octet d'une trame
                self.infilebytes=2
                self.scalemax=65535                
            self.FrameIndex=-1             # Index de trame, on evite les deux premieres
            self.offset=178               # Offset de l'entete fichier ser
            self.fileoffset=178 #MattC to avoid stomping on offset accumulator
            
        elif (self.AVI_flag==True): #MattC 
    	    #deal with avi file
            self.serfile = cv2.VideoCapture(serfile)

            self.Width = int(self.serfile.get(cv2.CAP_PROP_FRAME_WIDTH))
            self.Height = int(self.serfile.get(cv2.CAP_PROP_FRAME_HEIGHT))
            self.PixelDepthPerPlane=1*8
            self.FrameCount = int(self.serfile.get(cv2.CAP_PROP_FRAME_COUNT))            
            self.count=self.Width*self.Height
            self.infilebytes=1
            self.scalemax=255
            self.FrameIndex=0 #MattC cv2 starts indexing at 1?
            self.offset = 0
            self.fileoffset = 0 #MattC to avoid stomping on offset accumulator
        else : #MattC
    	    ok_flag=False
     
        if self.Width>self.Height:
            self.flag_rotate=True
            self.ih=self.Width
            self.iw=self.Height
        else:
            self.flag_rotate=False
            self.iw=self.Width
            self.ih=self.Height
        

    def next_frame(self):
        self.FrameIndex += 1
        self.offset=self.fileoffset+self.FrameIndex*self.count*self.infilebytes #MattC track offset
      
        if (self.SER_flag==True): #MattC
            img=np.fromfile(self.serfile, dtype=self.infiledatatype,count=self.count, offset=self.offset)
        elif (self.AVI_flag==True):
            ret, img = self.serfile.read()
            img=cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        else:
            img=np.fromfile(self.serfile, dtype=self.infiledatatype,count=self.count, offset=self.offset)

        img=np.reshape(img,(self.Height,self.Width))
        
        if self.flag_rotate:
            img=np.rot90(img)
        return img

    def has_frames(self):
        return self.FrameIndex+1 < self.FrameCount




