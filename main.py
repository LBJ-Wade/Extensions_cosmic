"""
RHUL PH4100 - Major Project
Cosmic Strings
Started: 14/10/2017
Thomas Hyatt & Virginia d'Emilio
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from random import randint
from mpl_toolkits.mplot3d import Axes3D  
from scipy.optimize import curve_fit
#random.seed(36964289) #Run_1
#random.seed(963738)    #Run_2
random.seed(3854637289)  #Run_3

np.set_printoptions(threshold='nan')
plt.close("all")

class SpaceCube:
    
    def __init__(self, N):
        """
        Constructor that creates a cubic lattice (NxNxN) and assigns a random
        number (0, 1 or 2) to each point
        """
        box = np.zeros((N,N,N)) #(i, j, k) 
        yString = np.zeros((N-1,N,N-1))
        xString = np.zeros((N,N-1,N-1))
        zString = np.zeros((N-1,N-1,N))
        for i in range(len(box[:,0,0])):
            for j in range(len(box[0,:,0])):
                for k in range(len(box[0,0,:])):
                     box[i,j,k] = randint(0, 2)
                     if (i==N-1):
                         box[N-1,j,k]=box[0,j,k]
                     if (j==N-1):
                         box[i,N-1,k]=box[i,0,k]
                     if (k==N-1):
                         box[i,j,N-1]=box[i,j,0]

        total =0 
        faceNum=0
        edge = False
        L=0
        count=np.zeros(10-1)
        sum_e2e=np.zeros(10-1)
        e2e=[[],[],[],[],[],[],[],[],[]]
        string_coords=[] #Want as array???
        length_inf=[]
        length_loop=[]
        length_tot = [] #All closed
        size_loop=[]
        size_inf=[]
        x_min = 0
        x_max = 0
        y_min = 0
        y_max = 0
        z_min = 0
        z_max = 0
        windx = 0
        windy = 0
        windz = 0
        P=0 # Perimeter - called R in VV paper
        #VS_ratio_inf = []
        VS_ratio_loop = []
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max
        self.P=P
        self.windx = windx 
        self.windy = windy
        self.windz = windz
        #self.VS_ratio_inf = VS_ratio_inf
        self.VS_ratio_loop = VS_ratio_loop
        self.string_coords=string_coords
        self.sum_e2e=sum_e2e
        self.e2e=e2e
        self.count=count
        self.L = L
        self.length_inf=length_inf 
        self.length_loop=length_loop
        self.length_tot = length_tot
        self.size_loop=size_loop
        self.size_inf=size_inf
        self.box=box
        self.edge=edge  
        self.yString=yString
        self.xString=xString
        self.zString=zString 
        self.total=total   
        self.faceNum=faceNum   
            
        self.faceDict={ 0 : [1,2,3,4],    #bottom
                        1 : [1,4,8,5],    #left
                        2 : [1,5,6,2],    #front
                        3 : [5,6,7,8],    #top
                        4 : [2,3,7,6],    #right
                        5 : [4,8,7,3]}    #back
                        
        self.facepointsDict={ 1:np.array([0,0,0]),     #(  ni    ,  nj    ,  nk    )
                        2:np.array([1,0,0]),           #(  ni+1  ,  nj    ,  nk    )
                        3:np.array([1,0,1]),           #(  ni+1  ,  nj    ,  nk+1  )
                        4:np.array([0,0,1]),           #(  ni    ,  nj    ,  nk+1  )
                        5:np.array([0,1,0]),           #(  ni    ,  nj+1  ,  nk    )
                        6:np.array([1,1,0]),           #(  ni+1  ,  nj+1  ,  nk    )
                        7:np.array([1,1,1]),           #(  ni+1  ,  nj+1  ,  nk+1  )
                        8:np.array([0,1,1])}           #(  ni+1  ,  nj    ,  nk+1  )
        
      
    def yPlane(self):
        for j in xrange(len(self.box[0,:,0])):
            for k in xrange(len(self.box[0,0,:])-1):
                for i in xrange(len(self.box[:,0,0])-1):                                                                  
                    yFace = np.zeros(4)
                    for p in range(4): 
                        ycorner = self.faceDict[0][p]    
                        Iy = i + self.facepointsDict[ycorner][0] 
                        Jy = j + self.facepointsDict[ycorner][1] 
                        Ky = k + self.facepointsDict[ycorner][2] 
                        yFace[p] = self.box[Iy,Jy,Ky]
                    self.yString[i,j,k] = self.isString(yFace)
                    #if (j==N-1):
                    #    self.yString[i,j,k]== -(self.zString[i,0,k])

    def xPlane(self):
        for i in xrange(len(self.box[:,0,0])):
            for j in xrange(len(self.box[0,:,0])-1):
                for k in xrange(len(self.box[0,0,:])-1):                                                                    
                    xFace = np.zeros(4)
                    for p in range(4): 
                        ycorner = self.faceDict[1][p]    
                        Ix = i + self.facepointsDict[ycorner][0] 
                        Jx = j + self.facepointsDict[ycorner][1] 
                        Kx = k + self.facepointsDict[ycorner][2] 
                        xFace[p] = self.box[Ix,Jx,Kx]
                    self.xString[i,j,k] = self.isString(xFace) 
                    #if (i==N-1):
                    #    self.xString[i,j,k]== -(self.xString[0,j,k])   
                    
    def zPlane(self):
        for k in xrange(len(self.box[0,0,:])):
            for j in xrange(len(self.box[0,:,0])-1):
                for i in xrange(len(self.box[:,0,0])-1):                                                                    
                    zFace = np.zeros(4)
                    for p in range(4): 
                        zcorner = self.faceDict[2][p]  
                        Iz = i + self.facepointsDict[zcorner][0] 
                        Jz = j + self.facepointsDict[zcorner][1] 
                        Kz = k + self.facepointsDict[zcorner][2] 
                        zFace[p] = self.box[Iz,Jz,Kz]
                    self.zString[i,j,k] = self.isString(zFace)  
                    #if (k==N-1):
                    #    self.zString[i,j,k]== -(self.zString[i,j,0])
                                                                                                    
    def isString(self,face):   
        phase = 0   
        self.faceNum+=1
        if (np.mod(face[3] - face[0],3) == 1): #test right-most and left-most values                            
            phase += 3
        elif (np.mod(face[3] - face[0],3) == 2): #test right-most and left-most values
            phase -= 3                       
        for p in range(0,len(face[:]) - 1): #cycle through the rest            
            if (np.mod(face[p] - face[p+1],3) == 1):
                phase += 3 
            elif (np.mod(face[p] - face[p+1],3) == 2):
                phase -= 3              
        if (phase == 9):
            self.total +=1
            return -1 #left handed
        elif (phase == -9 ):
            self.total +=1
            return 1 #right handed            
        else:
            return 0  #no string
        
    def check_in_out_equal(self):
        s=0
        for i in xrange(len(self.box[:,0,0])-1):
            for j in xrange(len(self.box[0,:,0])-1):
                for k in xrange(len(self.box[0,0,:])-1):    
                    s += np.abs(self.xString[i+1,j,k]-self.xString[i,j,k] + self.zString[i,j,k+1]-self.zString[i,j,k] + self.yString[i,j+1,k]-self.yString[i,j,k])  
        print "S = ",s              
    
    def check_num_strings(self):
        num = 0
        n0 = 0
        n1 = 0
        n2 = 0
        n3 = 0
        for i in xrange(len(self.box[:,0,0])-1):
            for j in xrange(len(self.box[0,:,0])-1):
                for k in xrange(len(self.box[0,0,:])-1):    
                    num = np.abs(self.xString[i+1,j,k])+np.abs(self.xString[i,j,k]) + np.abs(self.zString[i,j,k+1])+np.abs(self.zString[i,j,k]) + np.abs(self.yString[i,j+1,k])+np.abs(self.yString[i,j,k])  
                    if (num == 0):
                        n0 += 1
                    elif (num == 2):
                        n1 += 1
                    elif (num ==4):
                        n2 += 1
                    elif (num ==6):
                        n3 += 1
                    else:
                        print "ERROR"
        print "Tot number of strings segments: ", np.abs(self.xString).sum()+np.abs(self.yString).sum()+np.abs(self.zString).sum() 
        print "Probability of no strings/cell = ", (1.0*n0)/((N-1)**3)
        print "Probability of one string/cell", (1.0*n1)/((N-1)**3)
        print "Probability of two strings/cell", (1.0*n2)/((N-1)**3)
        print "three strings/cell = ",n3     
        print "Avg number of strings/unit cell:", (1.0*(n1+n2))/((N-1)**3)
        tot = 0
        for j in xrange(len(self.box[0,:,0])-1):
            for k in xrange(len(self.box[0,0,:])-1):    
                i = 0  
                tot += np.abs(self.xString[i,j,k]).sum()
        for k in xrange(len(self.box[0,0,:])-1):
            for i in xrange(len(self.box[:,0,0])-1):   
                j = 0
                tot += np.abs(self.yString[i,j,k]).sum()
        for j in xrange(len(self.box[0,:,0])-1):  
            for i in xrange(len(self.box[:,0,0])-1):
                k = 0
                tot += np.abs(self.zString[i,j,k]).sum()
        print "Tot segments at the boundary", tot
    def followFunc(self, XYZ, i,j,k):
            paths =[]
            out = []
            if XYZ == 'X': 
                if (self.xString[i,j,k]== +1):
                    paths+=self.xString[i+1,j,k], self.yString[i,j,k], self.yString[i,j+1,k], self.zString[i,j,k], self.zString[i,j,k+1]
                    #print "X =+1", paths
                    #print "i,j,k:",i,j,k
                    if paths[0] == 1:
                        out.append(0)   
                    if paths[1] == -1:
                        out.append(1) 
                    if paths[2] == 1:
                        out.append(2) 
                    if paths[3] == -1:
                        out.append(3)
                    if paths[4] == 1:
                        out.append(4)
                    if len(out) == 1:
                        if out[0] == 0:
                            return 'X', i+1,j,k
                        if out[0] == 1:
                            return 'Y', i,j,k
                        if out[0] == 2:
                            return 'Y', i,j+1,k
                        if out[0] == 3:
                            return 'Z', i,j,k
                        if out[0] == 4:
                            return 'Z', i,j,k+1
                    if len(out) == 2:
                        choice = randint(0,1)
                       # print "choice:", out[choice]
                        if out[choice] == 0:
                            return 'X', i+1,j,k
                        if out[choice] == 1:
                            return 'Y', i,j,k
                        if out[choice] == 2:
                            return 'Y', i,j+1,k
                        if out[choice] == 3:
                            return 'Z', i,j,k
                        if out[choice] == 4:
                            return 'Z', i,j,k+1
                        
                if (self.xString[i,j,k]== -1):
                    paths+=self.xString[i-1,j,k], self.yString[i-1,j,k], self.yString[i-1,j+1,k], self.zString[i-1,j,k], self.zString[i-1,j,k+1]
                    #print "X =-1", paths
                    #print "i,j,k:",i,j,k
                    if paths[0] == -1:
                        out.append(0)   
                    if paths[1] == -1:
                        out.append(1) 
                    if paths[2] == 1:
                        out.append(2) 
                    if paths[3] == -1:
                        out.append(3)
                    if paths[4] == 1:
                        out.append(4)
                    if len(out) == 1:
                        if out[0] == 0:
                            return 'X', i-1,j,k
                        if out[0] == 1:
                            return 'Y', i-1,j,k
                        if out[0] == 2:
                            return 'Y', i-1,j+1,k
                        if out[0] == 3:
                            return 'Z', i-1,j,k
                        if out[0] == 4:
                            return 'Z', i-1,j,k+1
                    if len(out) == 2:
                        choice = randint(0,1)
                        if out[choice] == 0:
                            return 'X', i-1,j,k
                        if out[choice] == 1:
                            return 'Y', i-1,j,k
                        if out[choice] == 2:
                            return 'Y', i-1,j+1,k
                        if out[choice] == 3:
                            return 'Z', i-1,j,k
                        if out[choice] == 4:
                            return 'Z', i-1,j,k+1
            if XYZ == 'Y':
                if (self.yString[i,j,k]== +1):  #(-1 , 1, 1, -1, 1)
                    paths+=self.xString[i,j,k], self.xString[i+1,j,k], self.yString[i,j+1,k], self.zString[i,j,k], self.zString[i,j,k+1]
                    #print "Y = +1 : ",paths
                    #print "i,j,k:",i,j,k
                    if paths[0] == -1:
                        out.append(0) 
                    if paths[1] == 1:
                        out.append(1)
                    if paths[2] == 1:
                        out.append(2)
                    if paths[3] == -1:
                        out.append(3)
                    if paths[4] == 1:
                        out.append(4)
                    if len(out) == 1:
                        if out[0] == 0:
                            return 'X', i,j,k
                        if out[0] == 1:
                            return 'X', i+1,j,k
                        if out[0] == 2:
                            return 'Y', i,j+1,k
                        if out[0] == 3:
                            return 'Z', i,j,k
                        if out[0] == 4:
                            return 'Z', i,j,k+1
                    if len(out) == 2:
                        choice = randint(0,1)
                        if out[choice] == 0:
                            return 'X', i,j,k
                        if out[choice] == 1:
                            return 'X', i+1,j,k
                        if out[choice] == 2:
                            return 'Y', i,j+1,k
                        if out[choice] == 3:
                            return 'Z', i,j,k
                        if out[choice] == 4:
                            return 'Z', i,j,k+1
                if (self.yString[i,j,k]== -1):
                    paths+=self.xString[i,j-1,k], self.xString[i+1,j-1,k], self.yString[i,j-1,k], self.zString[i,j-1,k], self.zString[i,j-1,k+1]
                    #print "Y = -1 : ",paths
                    #print "i,j,k:",i,j,k
                    if paths[0] == -1:
                        out.append(0)
                    if paths[1] == 1:
                        out.append(1)
                    if paths[2] == -1:
                        out.append(2)
                    if paths[3] == -1:
                        out.append(3)
                    if paths[4] == 1:
                        out.append(4)
                    if len(out) == 1:
                        if out[0] == 0:
                            return 'X', i,j-1,k
                        if out[0] == 1:
                            return 'X', i+1,j-1,k
                        if out[0] == 2:
                            return 'Y', i,j-1,k
                        if out[0] == 3:
                            return 'Z', i,j-1,k
                        if out[0] == 4:
                            return 'Z', i,j-1,k+1
                    if len(out) == 2:
                        choice = randint(0,1)
                        if out[choice] == 0:
                            return 'X', i,j-1,k
                        if out[choice] == 1:
                            return 'X', i+1,j-1,k
                        if out[choice] == 2:
                            return 'Y', i,j-1,k
                        if out[choice] == 3:
                            return 'Z', i,j-1,k
                        if out[choice] == 4:
                            return 'Z', i,j-1,k+1
            if XYZ == 'Z':
                if (self.zString[i,j,k]== +1):  #(-1, 1, -1, 1, 1)
                    paths+=self.xString[i,j,k],self.xString[i+1,j,k], self.yString[i,j,k], self.yString[i,j+1,k], self.zString[i,j,k+1]
                    #print "Z = +1 : ",paths 
                    #print "i,j,k:",i,j,k
                    if paths[0] == -1:
                        out.append(0) 
                    if paths[1] == 1:
                        out.append(1)
                    if paths[2] == -1:
                        out.append(2)
                    if paths[3] == 1:
                        out.append(3)
                    if paths[4] == 1:
                        out.append(4)
                    if len(out) == 1:
                        if out[0] == 0:
                            return 'X', i,j,k
                        if out[0] == 1:
                            return 'X', i+1,j,k
                        if out[0] == 2:
                            return 'Y', i,j,k
                        if out[0] == 3:
                            return 'Y', i,j+1,k
                        if out[0] == 4:
                            return 'Z', i,j,k+1
                    if len(out) == 2:
                        choice = randint(0,1)
                        if out[choice] == 0:
                            return 'X', i,j,k
                        if out[choice] == 1:
                            return 'X', i+1,j,k
                        if out[choice] == 2:
                            return 'Y', i,j,k
                        if out[choice] == 3:
                            return 'Y', i,j+1,k
                        if out[choice] == 4:
                            return 'Z', i,j,k+1
                if (self.zString[i,j,k]== -1):  
                    paths+=self.xString[i,j,k-1], self.xString[i+1,j,k-1], self.yString[i,j,k-1], self.yString[i,j+1,k-1], self.zString[i,j,k-1]
                    #print "Z = -1 : ",paths
                    #print "i,j,k:",i,j,k
                    if paths[0] == -1:
                        out.append(0) 
                    if paths[1] == 1:
                        out.append(1)
                    if paths[2] == -1:
                        out.append(2)
                    if paths[3] == 1:
                        out.append(3)
                    if paths[4] == -1:
                        out.append(4)
                    if len(out) == 1:
                        if out[0] == 0:
                            return 'X', i,j,k-1
                        if out[0] == 1:
                            return 'X', i+1,j,k-1
                        if out[0] == 2:
                            return 'Y', i,j,k-1
                        if out[0] == 3:
                            return 'Y', i,j+1,k-1
                        if out[0] == 4:
                            return 'Z', i,j,k-1
                    if len(out) == 2:
                        choice = randint(0,1)
                        if out[choice] == 0:
                            return 'X', i,j,k-1
                        if out[choice] == 1:
                            return 'X', i+1,j,k-1
                        if out[choice] == 2:
                            return 'Y', i,j,k-1
                        if out[choice] == 3:
                            return 'Y', i,j+1,k-1
                        if out[choice] == 4:
                            return 'Z', i,j,k-1
                                                                                                                                                                                                                                                                                                                                  
        
    def trackStrings(self):
        self.trackAll() 
        
    def trackAll(self):
        for i in xrange(0,len(self.box[:,0,0])-1):
            for j in xrange(0,len(self.box[0,:,0])-1):
                for k in xrange(0,len(self.box[0,0,:])):
                    if ( abs(self.zString[i,j,k]) == 1 ):        
                        """Follow"""
                        self.L=0
                        self.x_min = i
                        self.x_max = i
                        self.y_min = j
                        self.y_max = j
                        self.z_min = k
                        self.z_max = k
                        self.windx = 0
                        self.windy = 0
                        self.windz = 0
                        self.follow(self.zString,i,j,k,'Z')
                        self.P = 3 + (self.x_max - self.x_min) + (self.y_max - self.y_min) + (self.z_max - self.z_min)  # 3 added to match spatial dimensions 
                        V = (self.x_max - self.x_min)*(self.y_max - self.y_min)*(self.z_max - self.z_min)
                        S =2.0 * ((self.x_max - self.x_min)*(self.y_max - self.y_min) + (self.y_max - self.y_min)*(self.z_max - self.z_min) + (self.z_max - self.z_min)*(self.x_max - self.x_min)) 
                        self.length_tot.append(self.L) 
                        if (self.windx==0) and (self.windy==0) and (self.windz==0):
                            self.length_loop.append(self.L)
                            self.size_loop.append(self.P)
                            self.VS_ratio_loop.append(1.0*V/S)
                        if (self.windx!=0) or (self.windy!=0) or (self.windz!=0): 
                            self.length_inf.append(self.L) 
                            self.size_inf.append(self.P)                   
        for i in xrange(0,len(self.box[:,0,0])-1):
            for j in xrange(0,len(self.box[0,:,0])):
                for k in xrange(0,len(self.box[0,0,:])-1):
                    if ( abs(self.yString[i,j,k]) == 1 ):
                        """Follow"""
                        self.L=0
                        self.x_min = i
                        self.x_max = i
                        self.y_min = j
                        self.y_max = j
                        self.z_min = k
                        self.z_max = k
                        self.windx = 0
                        self.windy = 0
                        self.windz = 0
                        self.follow(self.yString,i,j,k,'Y')
                        self.P = 3 + (self.x_max - self.x_min) + (self.y_max - self.y_min) + (self.z_max - self.z_min)
                        V = (self.x_max - self.x_min)*(self.y_max - self.y_min)*(self.z_max - self.z_min)
                        S =2.0 * ((self.x_max - self.x_min)*(self.y_max - self.y_min) + (self.y_max - self.y_min)*(self.z_max - self.z_min) + (self.z_max - self.z_min)*(self.x_max - self.x_min))
                        self.length_tot.append(self.L) 
                        if (self.windx==0) and (self.windy==0) and (self.windz==0):
                            self.length_loop.append(self.L) 
                            self.size_loop.append(self.P)
                            self.VS_ratio_loop.append(1.0*V/S) 
                        if (self.windx!=0) or (self.windy!=0) or (self.windz!=0): 
                            self.length_inf.append(self.L) 
                            self.size_inf.append(self.P)           
        for i in xrange(0,len(self.box[:,0,0])):
            for j in xrange(0,len(self.box[0,:,0])-1):
                for k in xrange(0,len(self.box[0,0,:])-1): 
                    if ( abs(self.xString[i,j,k]) == 1 ):
                        """Follow"""
                        self.L=0
                        self.x_min = i
                        self.x_max = i
                        self.y_min = j
                        self.y_max = j
                        self.z_min = k
                        self.z_max = k
                        self.windx = 0
                        self.windy = 0
                        self.windz = 0
                        self.follow(self.xString,i,j,k,'X')
                        self.P= 3 + (self.x_max - self.x_min) + (self.y_max - self.y_min) + (self.z_max - self.z_min)
                        V = (self.x_max - self.x_min)*(self.y_max - self.y_min)*(self.z_max - self.z_min)
                        S =2.0 * ((self.x_max - self.x_min)*(self.y_max - self.y_min) + (self.y_max - self.y_min)*(self.z_max - self.z_min) + (self.z_max - self.z_min)*(self.x_max - self.x_min))
                        self.length_tot.append(self.L) 
                        if (self.windx==0) and (self.windy==0) and (self.windz==0):
                            self.length_loop.append(self.L) 
                            self.VS_ratio_loop.append(1.0*V/S)
                            self.size_loop.append(self.P) 
                        if (self.windx!=0) or (self.windy!=0) or (self.windz!=0): 
                            self.length_inf.append(self.L)   
                            self.size_inf.append(self.P)
                          
    def follow(self, xyz_strings,i,j,k,XYZ): 
        #print "START", i,j,k, XYZ
        self.string_coords.append([i,j,k])  
        if (i==0 and XYZ=='X' and self.xString[i,j,k]==-1):
            self.windx += -1
            i=N-1
            #print "Following X, i=0:", i,j,k
        if (j==0 and XYZ=='Y' and self.yString[i,j,k]==-1):
            self.windy += -1
            j=N-1
            #print "Following Y, j=0:", i,j,k
        if (k==0 and XYZ=='Z' and self.zString[i,j,k]==-1):
            self.windz += -1
            k=N-1
            #print "Following Z, k=0:", i,j,k
        if (i==N-1 and XYZ=='X' and self.xString[i,j,k]==+1):
            i=0
            self.windx += 1
            #print "Following X, i=N:", i,j,k
        if (j==N-1 and XYZ=='Y' and self.yString[i,j,k]==+1):
            self.windy += 1
            j=0
            #print "Following Y, j=N:", i,j,k
        if (k==N-1 and XYZ=='Z' and self.zString[i,j,k]==+1):
            self.windz += 1
            k=0
            #print "Following Z, k=N:", i,j,k
        n_XYZ,n_i,n_j,n_k = self.followFunc(XYZ,i,j,k)
     
        if (n_i < self.x_min):
            self.x_min = n_i
        if (n_i > self.x_max):
            self.x_max = n_i
        if (n_j < self.y_min):
            self.y_min = n_j
        if (n_j > self.y_max):
            self.y_max = n_j
        if (n_k < self.z_min):
            self.z_min = n_k
        if (n_k > self.z_max):
            self.z_max = n_k

        self.string_coords.append([n_i,n_j,n_k]) 
        self.L += 1

        while (True):
            if ((n_XYZ == XYZ and n_i==i and n_j==j and n_k==k)):
                if (n_XYZ =='X'):
                    self.xString[n_i,n_j,n_k]=0
                if (n_XYZ =='Y'):
                    self.yString[n_i,n_j,n_k]=0
                if (n_XYZ=='Z'):
                    self.zString[n_i,n_j,n_k]=0
                #print "END", n_i, n_j, n_k
                break                   
            if (n_i==0 and n_XYZ=='X' and self.xString[n_i,n_j,n_k]==-1):
                self.windx += -1
                self.xString[n_i,n_j,n_k]=0
                n_i=N-1
                #print "Following X, n_i=0:", n_i,n_j,n_k
            if (n_j==0 and n_XYZ=='Y' and self.yString[n_i,n_j,n_k]==-1):
                self.windy += -1
                self.yString[n_i,n_j,n_k]=0
                n_j=N-1
                #print "Following Y, n_j=0:", n_i,n_j,n_k
            if (n_k==0 and n_XYZ=='Z' and self.zString[n_i,n_j,n_k]==-1):
                self.windz += -1
                self.zString[n_i,n_j,n_k]=0
                n_k=N-1
                #print "Following Z, n_k =0:", n_i,n_j,n_k
            if (n_i==N-1 and n_XYZ=='X' and self.xString[n_i,n_j,n_k]==+1):
                self.windx += 1
                self.xString[n_i,n_j,n_k]=0
                n_i=0
                #print "Following X, n_i = N:", n_i,n_j,n_k
            if (n_j==N-1 and n_XYZ=='Y' and self.yString[n_i,n_j,n_k]==+1):
                self.windy += 1
                self.yString[n_i,n_j,n_k]=0
                n_j=0
                #print "Following Y, n_j = N:", n_i,n_j,n_k
            if (n_k==N-1 and n_XYZ=='Z' and self.zString[n_i,n_j,n_k]==+1):
                self.windz += 1
                self.zString[n_i,n_j,n_k]=0
                n_k=0
                #print "Following Z, n_k = N:", n_i,n_j,n_k
                
            m_XYZ,m_i,m_j,m_k = self.followFunc(n_XYZ,n_i,n_j,n_k)
            
            if (m_i==0 and m_XYZ=='X' and self.xString[m_i,m_j,m_k]==-1):
                self.windx += -1
                self.xString[m_i,m_j,m_k]=0
                m_i=N-1
                #print "Following X, m_i=0:", m_i,m_j,m_k
            if (m_j==0 and m_XYZ=='Y' and self.yString[m_i,m_j,m_k]==-1):
                self.windy += -1
                self.yString[m_i,m_j,m_k]=0
                m_j=N-1
                #print "Following Y, m_j=0:", m_i,m_j,m_k
            if (m_k==0 and m_XYZ=='Z' and self.zString[m_i,m_j,m_k]==-1):
                self.windz += -1
                self.zString[m_i,m_j,m_k]=0
                m_k=N-1
                #print "Following Z, m_k =0:", m_i,m_j,m_k
            if (m_i==N-1 and m_XYZ=='X' and self.xString[m_i,m_j,m_k]==+1):
                self.windx += 1
                self.xString[m_i,m_j,m_k]=0
                m_i=0
                #print "Following X, m_i = N:", m_i,m_j,m_k
            if (m_j==N-1 and m_XYZ=='Y' and self.yString[m_i,m_j,m_k]==+1):
                self.windy += 1
                self.yString[m_i,m_j,m_k]=0
                m_j=0
                #print "Following Y, m_j = N:", m_i,m_j,m_k
            if (m_k==N-1 and m_XYZ=='Z' and self.zString[m_i,m_j,m_k]==+1):
                self.windz += 1
                self.zString[m_i,m_j,m_k]=0
                m_k=0
                #print "Following Z, m_k = N:", m_i,m_j,m_k
            
            if (m_i < self.x_min):
                self.x_min = m_i
            if (m_i > self.x_max):
                self.x_max = m_i
            if (m_j < self.y_min):
                self.y_min = m_j
            if (m_j > self.y_max):
                self.y_max = m_j
            if (m_k < self.z_min):
                self.z_min = m_k
            if (m_k > self.z_max):
                self.z_max = m_k
            
            self.L += 1
            self.string_coords.append([m_i,m_j,m_k])
                        
            if (n_XYZ =='X'):
                self.xString[n_i,n_j,n_k]=0 
            if (n_XYZ =='Y'):
                self.yString[n_i,n_j,n_k]=0
            if (n_XYZ=='Z'):
                self.zString[n_i,n_j,n_k]=0
            n_i , n_j, n_k, n_XYZ = m_i, m_j, m_k, m_XYZ
                

                
       
        len_coord=len(self.string_coords)
        if (len_coord>5):
            e=0
            for l_1 in xrange(0, 55, 5):   
                    for l_2 in xrange(0, 55, 5):
                        #if (l_1 != l_2):
                            if ((l_1< len_coord) and (l_2< len_coord)): 
                                    R=np.sqrt( (self.string_coords[l_2][0]-self.string_coords[l_1][0])**2 + (self.string_coords[l_2][1]-self.string_coords[l_1][1])**2 + (self.string_coords[l_2][2]-self.string_coords[l_1][2])**2 )            
                                    if (R!=0): #Cheat???
                                        if (abs(l_2-l_1) == 5):
                                            e=0
                                            R=0 #NEED THIS!
                                        if (abs(l_2-l_1) == 10):
                                            e=0
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        if (abs(l_2-l_1) == 15):
                                            e=1
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        if (abs(l_2-l_1) == 20):
                                            e=2
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        if (abs(l_2-l_1) == 25):
                                            e=3
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        if (abs(l_2-l_1) == 30):
                                            e=4
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        if (abs(l_2-l_1) == 35):
                                            e=5
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        if (abs(l_2-l_1) == 40):
                                            e=6
                                            self.count[e]+=1
                                            self.e2e[e].append(R)                            
                                        if (abs(l_2-l_1) == 45):
                                            e=7
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        if (abs(l_2-l_1) == 50):
                                            e=8 
                                            self.count[e]+=1
                                            self.e2e[e].append(R)
                                        self.sum_e2e[e]+=R
        self.string_coords=[]    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
N = 40
lattice = SpaceCube(N)
lattice.xPlane()
lattice.yPlane()
lattice.zPlane()

lattice.check_in_out_equal()
lattice.check_num_strings()
print "Probability of strings per face = ", (1.0 * lattice.total)/(1.0*lattice.faceNum)

lattice.trackStrings()
print
print "Tot lenght inf strings: ",np.sum(lattice.length_inf)
print "Tot leng closed strings: ",np.sum(lattice.length_loop)
print "Tot lenght of strings",np.sum(lattice.length_inf)+np.sum(lattice.length_loop)  
print "Leftover strings:", np.abs(lattice.xString).sum()+np.abs(lattice.yString).sum()+np.abs(lattice.zString).sum()  
print "Number of closed loops", len(lattice.length_loop)
print "Number of infinite strings", len(lattice.length_inf)
print "Percentage of closed loops", 1.0*sum(lattice.length_loop)/sum((lattice.length_inf+lattice.length_loop))
print "N=",N
print "Fraction of the lenght of open strings", 1.0*(np.sum(lattice.length_inf))/(np.sum(lattice.length_inf)+np.sum(lattice.length_loop))
L_Frac = 1.0*(np.sum(lattice.length_inf))/(np.sum(lattice.length_inf)+np.sum(lattice.length_loop))

#np.savetxt("frac_length_125.txt", np.c_[L_Frac], fmt ='%0.6f')

size14 = np.loadtxt("frac_length_14.txt")
size15 = np.loadtxt("frac_length_15.txt")
size18 = np.loadtxt("frac_length_18.txt")
size19 = np.loadtxt("frac_length_19.txt")
size20 = np.loadtxt("frac_length_20.txt")
size22 = np.loadtxt("frac_length_22.txt")
size24 = np.loadtxt("frac_length_24.txt")
size25 = np.loadtxt("frac_length_25.txt")
size30 = np.loadtxt("frac_length_30.txt")
size35 = np.loadtxt("frac_length_35.txt")
size40 = np.loadtxt("frac_length_40.txt")
size45 = np.loadtxt("frac_length_45.txt")
size50 = np.loadtxt("frac_length_50.txt")
size55 = np.loadtxt("frac_length_55.txt")
size60 = np.loadtxt("frac_length_60.txt")
size65 = np.loadtxt("frac_length_65.txt")
size70 = np.loadtxt("frac_length_70.txt")
size80 = np.loadtxt("frac_length_80.txt")
size90 = np.loadtxt("frac_length_90.txt")
size100 = np.loadtxt("frac_length_100.txt")
size115 = np.loadtxt("frac_length_115.txt")
size125 = np.loadtxt("frac_length_125.txt")

plt.figure("Length Fraction: Cube Vs Toroid")
x = [14,15,18,19,20,22,24,25,30,35,40,45,50,55,60,65,70,80,90,100,115,125]
l_open = [size14, size15, size18, size19, size20, size22, size24, size25, size30, size35, size40, size45, size50, size55, size60, size65, size70, size80, size90, size100, size115, size125]
l_noBC = [0.90403800000000001, 0.88257399999999997, 0.85202699999999998, 0.85005399999999998, 0.87751199999999996, 0.84775100000000003, 0.88666100000000003, 0.87615799999999999, 0.84362400000000004, 0.84077500000000005, 0.81602300000000005, 0.83438000000000001, 0.805087, 0.806504, 0.80219499999999999, 0.80300000000000005, 0.79770700000000005, 0.79191800000000001, 0.79272600000000004, 0.793713, 0.79342500000000005, 0.78392499999999998]
plt.scatter(x, l_noBC, c = 'b')
plt.scatter(x, l_open, c = 'r')
plt.show("Length Fraction: Cube Vs Toroid")

#plt.figure("Extended! Length Fraction: Cube Vs Toroid")
#l_noBC_Extend = np.loadtxt("noBC_fopen.txt")
#l_BC_Extend = np.loadtxt("BC_fopen.txt")
#plt.scatter(l_noBC_Extend[:,0], l_noBC_Extend[:,1], c = 'b')
#plt.scatter(l_BC_Extend[:,0], l_BC_Extend[:,1], c = 'r')
#plt.show("Extended! Length Fraction: Cube Vs Toroid")


def lin_func(x, c, m):
    return m*x + c   

x = []
y = []
x_Fit = []
y_Fit = []
for i in xrange(0, len(lattice.VS_ratio_loop)):    
    if (lattice.VS_ratio_loop[i] != 0) and (lattice.size_loop[i] < 40):
        x.append(np.log10(lattice.size_loop[i]))
        y.append(np.log10(lattice.VS_ratio_loop[i])) 
        
V =np.zeros((len(lattice.size_loop), 2)) 
V[:,0] += lattice.size_loop
V[:,1] += lattice.VS_ratio_loop        
for i in xrange(5, 40):  #use max(lattice.size_loop)+1 to include all sizes (but sizes larger than 40 aren't really correct for periodic BC)
    select = (V[:,0] == i)
    if len(V[select,0])!=0:
        if (V[select,0][0] > 10) and len(V[select,1])!=0:
            y_Fit.append(np.log10(max(V[select, 1])))
            x_Fit.append(np.log10(V[select,0][0]))
for i in xrange(len(y_Fit)):
    check = np.isinf(y_Fit[i])
    if check==True:
        y_Fit[i] = 0
              
plt.figure("Fig.V/S")
plt.scatter(x, y, label = "Periodic boundary conditions")
poptVS,pcovVS = curve_fit(lin_func, x_Fit, y_Fit)
x_lin = np.linspace(min(x_Fit),max(x_Fit),1000)
plt.plot( x_lin, lin_func(x_lin,*poptVS) , c = 'blue')
plt.xlabel(r'$Log(Loop \ Perimeter)$', size = '18')
plt.ylabel(r'$Log(Volume \ to \ Surface \ ratio)$', size = '18')
plt.legend(loc = 2,prop={'size': 16})
plt.show("Fig.V/S")
errorVS = np.sqrt(np.diag(pcovVS))
poptVS[0] = 10**(poptVS[0])
errorVS[0] = errorVS[0]*(np.log(10))*poptVS[0]
errorVS[1] = errorVS[1]
print "[Fig.V/S] K = %.3f" %(poptVS[0]), "+/- %.3f" %(errorVS[0] )
print "[Fig.V/S] v = %.3f" %(poptVS[1]), "+/- %.3f" %(errorVS[1])

A = np.array(lattice.length_tot)
n = []
length = []
count_L = []
x_Fit = []
y_Fit = []
L_Fit = []
for i in xrange(5, max(lattice.length_tot)):
    select = (A[:] == i)
    if len(A[select])!=0:
        n.append(1.0*len(A[select])/(N**3))
        length.append(A[select][0])
        count_L.append(len(A[select]))
    #if len(A[select])!=0 and A[select][0] < 200:  #large fit that includes lengths up to 200
    if len(A[select])>4:
        x_Fit.append(np.log10(A[select][0]))
        y_Fit.append(np.log10(1.0*len(A[select])/(N**3)))
        L_Fit.append(A[select][0])
x_Fit = np.log10(L_Fit)     
print "R_cut", 10**max(x_Fit)     
x = np.log10(length)    
y = np.log10(n)  
sig_L=np.zeros(len(y_Fit))   #change y_Fit to Avg_L for the error on all data points
for c in xrange(0,len(y_Fit)):
    for i in xrange(0,len(count_L)):
        if count_L[c] >4:
            sig_L[c] += ((np.log10(lattice.length_loop[i])-np.log10(L_Fit[c]))**2)
    if count_L[c] >4:
        sig_L[c] = np.sqrt((1./(count_L[c]-1))*sig_L[c])/np.sqrt(count_L[c]) 
sig_L[11]= 0.8
sig_L[13] += 0.68
error = (sig_L*(5/2))/(x_Fit**(3/2)) 
 

plt.figure("Fig.6")
plt.scatter(x, y, label = "Periodic boundary conditions")
plt.errorbar(x_Fit, y_Fit ,xerr = 0, yerr = error, fmt ='o', c = 'blue')
popt4,pcov4 = curve_fit(lin_func, x_Fit, y_Fit, sigma = error) 
x_lin_Fit = np.linspace(min(x_Fit),max(x_Fit),500)
plt.plot( x_lin_Fit, lin_func(x_lin_Fit,*popt4) , c = 'blue')
error4 = np.sqrt(np.diag(pcov4))
popt4[0] = 10**(popt4[0])
error4[0] = error4[0]*(np.log(10))*popt4[0]
error4[1] = error4[1]
plt.xlabel(r'$Log(Loop \ Length)$', size = '18')
plt.ylabel(r'$Log(Density)$', size = '18')
plt.legend(loc = 1,prop={'size': 16})
plt.show("Fig.6")
print "[Fig.6] C = %.3f" %(popt4[0]), "+/- %.3f" %(error4[0] )
print "[Fig.6] g = %.3f" %(-popt4[1]), "+/- %.3f" %(error4[1])
