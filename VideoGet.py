from threading import Thread
import cv2

class VideoGet:
    """
    Class that continuously gets frames from a VideoCapture object
    with a dedicated thread.
    """

    def __init__(self, src=0):
        self.stream = cv2.VideoCapture(src)
        (self.grabbed, self.frame) = self.stream.read()
        self.stopped = False
        self.Width = int(self.src.get(cv2.CAP_PROP_FRAME_WIDTH))
        self.Height = int(self.src.get(cv2.CAP_PROP_FRAME_HEIGHT))
        self.PixelDepthPerPlane=1*8
        self.FrameCount = int(self.src.get(cv2.CAP_PROP_FRAME_COUNT))            
        self.count=self.Width*self.Height
        self.FrameIndex=-1
        self.offset = 0
        self.fileoffset = 0 #MattC to avoid stomping on offset accumulator

        if self.Width>self.Height:
            self.flag_rotate=True
            self.ih=self.Width
            self.iw=self.Height
        else:
            self.flag_rotate=False
            self.iw=self.Width
            self.ih=self.Height  
                      
    def start(self):    
        Thread(target=self.get, args=()).start()
        return self

    def get(self):
        while not self.stopped:
            if not self.grabbed:
                self.stop()
            else:
                self.FrameIndex += 1
                self.offset=self.fileoffset+self.FrameIndex*self.count*2 #MattC track offset
               
                (self.grabbed, self.frame) = self.stream.read()
                self.frame=cv2.cvtColor(self.frame, cv2.COLOR_BGR2GRAY)
                self.frame=np.reshape(self.frame,(self.Height,self.Width))
        
                if self.flag_rotate:
                    self.frame=np.rot90(self.frame)
            
    def stop(self):
        self.stopped = True
        
    def has_frames(self):
        return self.FrameIndex+1 < self.FrameCount
