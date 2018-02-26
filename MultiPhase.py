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

def PrintPnF(i,j,k):
    '''
    Prints the points of the cube around the bottom left identifiter ijk
    and also prints the string state of each face of that cube
    '''
    print
    print "P1:(0,0,0) = ",lattice.box[i,j,k]
    print "P2:(1,0,0) = ",lattice.box[i+1,j,k]
    print "P3:(1,0,1) = ",lattice.box[i+1,j,k+1]
    print "P4:(0,0,1) = ",lattice.box[i,j,k+1]
    print "P5:(0,1,0) = ",lattice.box[i,j+1,k]
    print "P6:(1,1,0) = ",lattice.box[i+1,j+1,k]
    print "P7:(1,1,1) = ",lattice.box[i+1,j+1,k+1]
    print "P8:(0,1,1) = ",lattice.box[i,j+1,k+1]
    print
    print "X(0,0,0) = ",lattice.yString[i,j,k]
    print "X(0,1,0) = ",lattice.yString[i,j+1,k]
    print "Y(0,0,0) = ",lattice.xString[i,j,k]
    print "Y(1,0,0) = ",lattice.xString[i+1,j,k]
    print "Z(0,0,0) = ",lattice.zString[i,j,k]
    print "Z(0,0,1) = ",lattice.zString[i,j,k+1]
    
def PlotLengthHist():
    figHist=plt.figure("Histogram", figsize=(16,9))   
    bins = range(min(lattice.length_inf), max(lattice.length_loop))
    plt.hist(lattice.length_inf, bins, histtype= 'bar', color ='r', label = r'$Infinite \ strings$', alpha=0.5)
    plt.hist(lattice.length_loop, bins, histtype= 'bar', color = 'b', label = r'$Closed \ strings$', alpha=0.5)
    plt.xlabel(r'$Length \ of \ Strings$', fontsize=18)
    plt.ylabel(r'$Number \ of \ Strings$', fontsize=18)
    plt.title(r'$Histogram \ of \ String \ Lengths$', fontsize=18, y=1.025)
    plt.legend(loc='upper right', fontsize=18)
    plt.annotate(r'$Size: \ N\xi \ = \ {0}$'.format(N), xy=(750, 325), xycoords='figure points', fontsize=18,
    bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.25'))
    plt.show("Histogram")

def Plot3DStrings():
    fig=plt.figure("3DFig", figsize=(8.0*1.2,6.0*1.2))
    ax = Axes3D(fig)  
    for n in xrange(0,len(lattice.tot_loop_coord_i)-1):
        ax.plot3D(lattice.tot_loop_coord_i[n],lattice.tot_loop_coord_j[n],lattice.tot_loop_coord_k[n], color='blue', label=r'$Closed \ Strings$'  if n == 0 else "", linewidth = 3)
    for n in xrange(0,len(lattice.tot_inf_coord_i)-1):
        ax.plot3D(lattice.tot_inf_coord_i[n],lattice.tot_inf_coord_j[n],lattice.tot_inf_coord_k[n], color='r', label=r'$Infinite \ Strings$'  if n == 0 else "", linewidth = 3)        
    ax.set_title(r"$3D \ Plot \ of \ Cosmic \ Strings $", fontsize=25, y=1.025, x=0.32)
    ax.set_xlabel(r"$X/ \xi$", fontsize=18, labelpad=10)
    ax.set_ylabel(r"$Y/ \xi$", fontsize=18, labelpad=10)
    ax.set_zlabel(r"$Z/ \xi$", fontsize=18, labelpad=10)
    ax.legend(loc = 'upper right', fontsize=22)
    ax.annotate(r'$Size: \ N\xi \ = \ {0}$'.format(N), xy=(525, 20), xycoords='figure points', fontsize=22,
    bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.25'))
    plt.show("3DFig")  
    fig1=plt.figure("3DFigInf", figsize=(8.0*1.2,6.0*1.2))
    ax1=Axes3D(fig1) 
    for n in xrange(0,len(lattice.tot_inf_coord_i)-1):
        ax1.plot3D(lattice.tot_inf_coord_i[n],lattice.tot_inf_coord_j[n],lattice.tot_inf_coord_k[n], color='r', label=r'$Infinite \ Strings$' if n == 0 else "", linewidth = 1.5)      
    ax1.set_title(r"$3D \ Plot \ of \ Infinite \ Strings $", fontsize=25, y=1.025, x=0.32)
    ax1.set_xlabel(r"$X/ \xi$", fontsize=18, labelpad=10)
    ax1.set_ylabel(r"$Y/ \xi$", fontsize=18, labelpad=10)
    ax1.set_zlabel(r"$Z/ \xi$", fontsize=18, labelpad=10)
    ax1.legend(loc = 'upper right', fontsize=22)
    ax1.annotate(r'$Size: \ N\xi \ = \ {0}$'.format(N), xy=(525, 20), xycoords='figure points', fontsize=22,
    bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.25'))
    plt.show("3DFig1")   
    fig2=plt.figure("3DFigLoop", figsize=(8.0*1.2,6.0*1.2))
    ax2=Axes3D(fig2)    
    for n in xrange(0,len(lattice.tot_loop_coord_i)-1):
        ax2.plot3D(lattice.tot_loop_coord_i[n],lattice.tot_loop_coord_j[n],lattice.tot_loop_coord_k[n], color='b', label=r'$Closed \ Strings$' if n == 0 else "", linewidth = 1.5)       
    ax2.set_title(r"$3D \ Plot \ of \ Closed \ Strings $", fontsize=25, y=1.025, x=0.32)
    ax2.set_xlabel(r"$X/ \xi$", fontsize=18, labelpad=10)
    ax2.set_ylabel(r"$Y/ \xi$", fontsize=18, labelpad=10)
    ax2.set_zlabel(r"$Z/ \xi$", fontsize=18, labelpad=10)
    ax2.legend(loc = 'upper right', fontsize=22)
    ax2.annotate(r'$Size: \ N\xi \ = \ {0}$'.format(N), xy=(525, 20), xycoords='figure points', fontsize=22,
    bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.25'))
    plt.show("3DFig1")

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
        size_loop=[]
        loop_coord_i=[]
        loop_coord_j=[]
        loop_coord_k=[]
        tot_loop_coord_i=[]
        tot_loop_coord_j=[]
        tot_loop_coord_k=[]
        inf_coord_i=[]
        inf_coord_j=[]
        inf_coord_k=[]
        tot_inf_coord_i=[]
        tot_inf_coord_j=[]
        tot_inf_coord_k=[]
        x_min = 0
        x_max = 0
        y_min = 0
        y_max = 0
        z_min = 0
        z_max = 0
        P=0 # Perimeter - called R in VV paper
        VS_ratio = []
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max
        self.P=P
        self.VS_ratio = VS_ratio 
        self.string_coords=string_coords
        self.sum_e2e=sum_e2e
        self.e2e=e2e
        self.count=count
        self.L = L
        self.length_inf=length_inf 
        self.length_loop=length_loop
        self.size_loop=size_loop
        self.loop_coord_i=loop_coord_i
        self.loop_coord_j=loop_coord_j
        self.loop_coord_k=loop_coord_k
        self.tot_loop_coord_i=tot_loop_coord_i
        self.tot_loop_coord_j=tot_loop_coord_j
        self.tot_loop_coord_k=tot_loop_coord_k
        self.inf_coord_i=inf_coord_i
        self.inf_coord_j=inf_coord_j
        self.inf_coord_k=inf_coord_k
        self.tot_inf_coord_i=tot_inf_coord_i
        self.tot_inf_coord_j=tot_inf_coord_j
        self.tot_inf_coord_k=tot_inf_coord_k
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
        print "Tot number of strings: ", np.abs(self.xString).sum()+np.abs(self.yString).sum()+np.abs(self.zString).sum()  
        print "Probability of no strings/cell = ", (1.0*n0)/((N-1)**3)
        print "Probability of one string/cell", (1.0*n1)/((N-1)**3)
        print "Probability of two strings/cell", (1.0*n2)/((N-1)**3)
        print "three strings/cell = ",n3     
        print "Avg number of strings/unit cell:", (1.0*(n1+n2))/((N-1)**3)
        
    def followFunc(self, XYZ, i,j,k):
           # print "In ijk: ", i, j, k
            paths =[]
            out = []
            if XYZ == 'X': 
              #  print "In X"
                if (self.xString[i,j,k]== +1):
                    paths+=self.xString[i+1,j,k], self.yString[i,j,k], self.yString[i,j+1,k], self.zString[i,j,k], self.zString[i,j,k+1]
                   # print paths
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
                  #  print paths
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
               # print "In Y"
                if (self.yString[i,j,k]== +1):  #(-1 , 1, 1, -1, 1)
                    paths+=self.xString[i,j,k], self.xString[i+1,j,k], self.yString[i,j+1,k], self.zString[i,j,k], self.zString[i,j,k+1]
                   # print paths
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
                   # print 
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
                #print "In Z"
                if (self.zString[i,j,k]== +1):  #(-1, 1, -1, 1, 1)
                    paths+=self.xString[i,j,k],self.xString[i+1,j,k], self.yString[i,j,k], self.yString[i,j+1,k], self.zString[i,j,k+1]
                   # print 
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
                    #print "found: ",self.zString[i,j,k]
                    paths+=self.xString[i,j,k-1], self.xString[i+1,j,k-1], self.yString[i,j,k-1], self.yString[i,j+1,k-1], self.zString[i,j,k-1]
                    #print "Prob", paths
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
        self.edge = True
        self.trackEdge() #Infinite Strings
        self.edge = False
        self.trackCentre() #Closed Strings
        
    def trackEdge(self):
        """Z-Edges"""
        for j in xrange(len(self.box[0,:,0])-1):  
            for i in xrange(len(self.box[:,0,0])-1):
                k = 0
                if ( self.zString[i,j,k] == 1 ):
                    """Follow"""
                    self.L=1
                    self.inf_coord_i=[]
                    self.inf_coord_j=[]
                    self.inf_coord_k=[]
                    self.follow(self.zString,i,j,k,'Z')  
                    self.tot_inf_coord_i.append(self.inf_coord_i)
                    self.tot_inf_coord_j.append(self.inf_coord_j)
                    self.tot_inf_coord_k.append(self.inf_coord_k)  
                    self.length_inf.append(self.L)                                  
                k = N-1
                if ( self.zString[i,j,k] == -1 ):
                    """Follow"""
                    self.L=1
                    self.inf_coord_i=[]
                    self.inf_coord_j=[]
                    self.inf_coord_k=[]
                    self.follow(self.zString,i,j,k,'Z')
                    self.tot_inf_coord_i.append(self.inf_coord_i)
                    self.tot_inf_coord_j.append(self.inf_coord_j)
                    self.tot_inf_coord_k.append(self.inf_coord_k) 
                    self.length_inf.append(self.L) 
        """Y-Edges"""    
        for k in xrange(len(self.box[0,0,:])-1):
            for i in xrange(len(self.box[:,0,0])-1):   
                j = 0    
                if ( self.yString[i,j,k] == 1 ):
                    """Follow"""
                    self.L=1
                    self.inf_coord_i=[]
                    self.inf_coord_j=[]
                    self.inf_coord_k=[]
                    self.follow(self.yString,i,j,k,'Y')
                    self.tot_inf_coord_i.append(self.inf_coord_i)
                    self.tot_inf_coord_j.append(self.inf_coord_j)
                    self.tot_inf_coord_k.append(self.inf_coord_k)  
                    self.length_inf.append(self.L)                                   
                j = N-1
                if ( self.yString[i,j,k] == -1 ):
                    """Follow"""
                    self.L=1
                    self.inf_coord_i=[]
                    self.inf_coord_j=[]
                    self.inf_coord_k=[]
                    self.follow(self.yString,i,j,k,'Y')
                    self.tot_inf_coord_i.append(self.inf_coord_i)
                    self.tot_inf_coord_j.append(self.inf_coord_j)
                    self.tot_inf_coord_k.append(self.inf_coord_k) 
                    self.length_inf.append(self.L)  
        """X-Edges"""    
        for j in xrange(len(self.box[0,:,0])-1):
            for k in xrange(len(self.box[0,0,:])-1):    
                i = 0    
                if ( self.xString[i,j,k] == 1 ):
                    """Follow"""
                    self.L=1
                    self.inf_coord_i=[]
                    self.inf_coord_j=[]
                    self.inf_coord_k=[]
                    self.follow(self.xString,i,j,k,'X')
                    self.tot_inf_coord_i.append(self.inf_coord_i)
                    self.tot_inf_coord_j.append(self.inf_coord_j)
                    self.tot_inf_coord_k.append(self.inf_coord_k)   
                    self.length_inf.append(self.L)                                   
                i = N-1
                if ( self.xString[i,j,k] == -1 ):
                    """Follow"""
                    self.L=1
                    self.inf_coord_i=[]
                    self.inf_coord_j=[]
                    self.inf_coord_k=[]
                    self.follow(self.xString,i,j,k,'X')
                    self.tot_inf_coord_i.append(self.inf_coord_i)
                    self.tot_inf_coord_j.append(self.inf_coord_j)
                    self.tot_inf_coord_k.append(self.inf_coord_k) 
                    self.length_inf.append(self.L)  
                                                                                                              
    def trackCentre(self):
        for i in xrange(0,len(self.box[:,0,0])-1):
            for j in xrange(0,len(self.box[0,:,0])-1):
                for k in xrange(1,len(self.box[0,0,:])-2):
                    if ( abs(self.zString[i,j,k]) == 1 ):
                        """Follow"""
                        self.L=0
                        self.x_min = i
                        self.x_max = i
                        self.y_min = j
                        self.y_max = j
                        self.z_min = k
                        self.z_max = k
                        self.loop_coord_i=[]
                        self.loop_coord_j=[]
                        self.loop_coord_k=[]
                        self.follow(self.zString,i,j,k,'Z')
                        self.tot_loop_coord_i.append(self.loop_coord_i)
                        self.tot_loop_coord_j.append(self.loop_coord_j)
                        self.tot_loop_coord_k.append(self.loop_coord_k)
                        self.P = 3 + (self.x_max - self.x_min) + (self.y_max - self.y_min) + (self.z_max - self.z_min)  # 3 added to match spatial dimensions
                        V = (self.x_max - self.x_min)*(self.y_max - self.y_min)*(self.z_max - self.z_min)
                        S =2.0 * ((self.x_max - self.x_min)*(self.y_max - self.y_min) + (self.y_max - self.y_min)*(self.z_max - self.z_min) + (self.z_max - self.z_min)*(self.x_max - self.x_min))
                        self.VS_ratio.append(1.0*V/S)
                        self.length_loop.append(self.L)
                        self.size_loop.append(self.P)                        
        for i in xrange(0,len(self.box[:,0,0])-1):
            for j in xrange(1,len(self.box[0,:,0])-2):
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
                        self.loop_coord_i=[]
                        self.loop_coord_j=[]
                        self.loop_coord_k=[]
                        self.follow(self.yString,i,j,k,'Y')
                        self.tot_loop_coord_i.append(self.loop_coord_i)
                        self.tot_loop_coord_j.append(self.loop_coord_j)
                        self.tot_loop_coord_k.append(self.loop_coord_k)
                        self.P = 3 + (self.x_max - self.x_min) + (self.y_max - self.y_min) + (self.z_max - self.z_min)
                        V = (self.x_max - self.x_min)*(self.y_max - self.y_min)*(self.z_max - self.z_min)
                        S =2.0 * ((self.x_max - self.x_min)*(self.y_max - self.y_min) + (self.y_max - self.y_min)*(self.z_max - self.z_min) + (self.z_max - self.z_min)*(self.x_max - self.x_min))
                        self.VS_ratio.append(1.0*V/S)
                        self.length_loop.append(self.L)
                        self.size_loop.append(self.P)           
        for i in xrange(1,len(self.box[:,0,0])-2):
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
                        self.loop_coord_i=[]
                        self.loop_coord_j=[]
                        self.loop_coord_k=[]
                        self.follow(self.xString,i,j,k,'X')
                        self.tot_loop_coord_i.append(self.loop_coord_i)
                        self.tot_loop_coord_j.append(self.loop_coord_j)
                        self.tot_loop_coord_k.append(self.loop_coord_k)
                        self.P= 3 + (self.x_max - self.x_min) + (self.y_max - self.y_min) + (self.z_max - self.z_min)
                        V = (self.x_max - self.x_min)*(self.y_max - self.y_min)*(self.z_max - self.z_min)
                        S =2.0 * ((self.x_max - self.x_min)*(self.y_max - self.y_min) + (self.y_max - self.y_min)*(self.z_max - self.z_min) + (self.z_max - self.z_min)*(self.x_max - self.x_min))
                        self.VS_ratio.append(1.0*V/S)
                        self.length_loop.append(self.L)
                        self.size_loop.append(self.P)  
                
                          
    def follow(self, xyz_strings,i,j,k,XYZ): 
        #Edge == True means looking for infinite strings
        #Edge == False means looking for closed strings
        self.string_coords.append([i,j,k])
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
        self.loop_coord_i.append(i)
        self.loop_coord_j.append(j)
        self.loop_coord_k.append(k)
        self.inf_coord_i.append(i)
        self.inf_coord_j.append(j)
        self.inf_coord_k.append(k)
        self.loop_coord_i.append(n_i)
        self.loop_coord_j.append(n_j)
        self.loop_coord_k.append(n_k)
        self.inf_coord_i.append(n_i)
        self.inf_coord_j.append(n_j)
        self.inf_coord_k.append(n_k)
        if (self.edge == False):
            while (True):
                if (XYZ == n_XYZ and n_i==i and n_j==j and n_k==k):
                    if (n_XYZ =='X'):
                        self.xString[n_i,n_j,n_k]=0
                    if (n_XYZ =='Y'):
                        self.yString[n_i,n_j,n_k]=0
                    if (n_XYZ=='Z'):
                        self.zString[n_i,n_j,n_k]=0
                    break                              
                m_XYZ,m_i,m_j,m_k = self.followFunc(n_XYZ,n_i,n_j,n_k)
                
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
                               
                self.loop_coord_i.append(m_i)
                self.loop_coord_j.append(m_j)
                self.loop_coord_k.append(m_k)
                if (n_XYZ =='X'):
                    self.xString[n_i,n_j,n_k]=0 
                if (n_XYZ =='Y'):
                    self.yString[n_i,n_j,n_k]=0
                if (n_XYZ=='Z'):
                    self.zString[n_i,n_j,n_k]=0
                n_i , n_j, n_k, n_XYZ = m_i, m_j, m_k, m_XYZ
                
        if (self.edge == True):
            if (XYZ =='X'):
                self.xString[i,j,k]=0
            if (XYZ =='Y'):
                self.yString[i,j,k]=0
            if (XYZ=='Z'):
                self.zString[i,j,k]=0
            while (True):               
                if (n_XYZ == 'X'):
                    if (n_i==N-1 or n_i==0):
                     self.xString[n_i,n_j,n_k]=0
                     #print "Break X" 
                     break
                if (n_XYZ == 'Y'):
                    if (n_j==N-1 or n_j==0):
                     self.yString[n_i,n_j,n_k]=0
                     #print "Break Y"
                     break
                if (n_XYZ == 'Z'):
                    if (n_k==N-1 or n_k==0):
                     self.zString[n_i,n_j,n_k]=0
                     #print "Break Z"
                     break
                self.L += 1 
                m_XYZ,m_i,m_j,m_k = self.followFunc(n_XYZ,n_i,n_j,n_k)
                self.string_coords.append([m_i,m_j,m_k])                   
                self.inf_coord_i.append(m_i)
                self.inf_coord_j.append(m_j)
                self.inf_coord_k.append(m_k)
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
#Plot3DStrings()
#PlotLengthHist()


Avg_e2e = lattice.sum_e2e/lattice.count
segment_length=[10,15,20,25,30,35,40,45,50]
sig_e2e = np.zeros(9)
summation = np.zeros(9)

for i in xrange(0,len(lattice.count)):
    for j in xrange(0,int(lattice.count[i])):
        summation[i] += (np.log10(Avg_e2e[i]) - np.log10(lattice.e2e[i][j]))**2 
        #if ((np.log10(Avg_e2e[i]) - np.log10(lattice.e2e[i][j]))**2 == float("inf")):
            #print "check1: ", np.log10(Avg_e2e[i]), " check2: ",np.log10(lattice.e2e[i][j])
    sig_e2e[i] = np.sqrt(1.0/(lattice.count[i]-1) * summation[i])/np.sqrt(lattice.count[i])


#np.savetxt("test_multirun_e2e_2.txt", np.c_[Avg_e2e,sig_e2e], fmt ='%0.6f')  #change seed and change file name, then run 
Run_1 = np.loadtxt("test_multirun_e2e.txt")
Run_2 = np.loadtxt("test_multirun_e2e_1.txt")            #RETAKE DATA for all runs
Run_3 = np.loadtxt("test_multirun_e2e_2.txt")
Avg_e2e_run = (Run_1[:,0]+Run_2[:,0]+Run_3[:,0])/3
segment_length = np.array(10-1)
segment_length=[10,15,20,25,30,35,40,45,50]
Error_e2e_run = np.sqrt((1./3)**2 * (Run_1[:,1]**2 + Run_2[:,1]**2 + Run_3[:,1]**2 ))

def lin_func(x, c, m):
    return m*x + c                
#y = np.log10(Avg_e2e_run)
y = np.log10(Avg_e2e)  
x = np.log10(segment_length)     
yerr = Error_e2e_run*10

plt.figure("Fig.2")
#plt.scatter(x, y)
plt.errorbar(x, y, xerr = 0, yerr = yerr, fmt ='o', c = 'blue', label = r"$Error \times 10$")
popt1,pcov1 = curve_fit(lin_func, x, y, sigma = Error_e2e_run)
x_lin = np.linspace(min(x)-0.1,max(x)+0.1,500)
plt.plot( x_lin, lin_func(x_lin,*popt1) , c = 'blue')
plt.xlabel(r'$Log(segment \ lenght)$', size = '18')
plt.ylabel(r'$Log(end \ to \ end \ distance)$', size = '18')
plt.title(r'$Estimation \ of \ the \ fractal \ dimension$', size = '18')
plt.legend(loc = 2,prop={'size': 18})
plt.show("Fig.2")
error = np.sqrt(np.diag(pcov1))
popt1[0] = 10**(popt1[0])
error[0] = error[0]*(np.log(10))*popt1[0]
error[1] = error[1]
popt1[1]=1.0/popt1[1]
print "[Fig.2] A = %.3f" %(popt1[0]), "+/- %.3f" %(error[0])
print "[Fig.2] d = %.3f" %(popt1[1]), "+/- %.3f" %(error[1])


Error_Range = []
Avg_L = []
P_Range = []
count_L=[]
P_Fit=[]
L_Fit=[]
A = np.zeros((len(lattice.size_loop), 2))
A[:,0] += lattice.length_loop
A[:,1] += lattice.size_loop
#i = 0
for i in xrange(5, max(lattice.size_loop)):
    select = (A[:,1] == i)
    if (len(A[select,0])!=0) and (A[select,1][0]>5):
    #if (len(A[select,0])!=0):
        Avg_L.append(sum(A[select,0])/len(A[select,0]))
        P_Range.append(A[select,1][0])
        count_L.append(len(A[select,0]))
        if(len(A[select,0])>4):
            P_Fit.append(A[select,1][0])
            L_Fit.append(sum(A[select,0])/len(A[select,0]))

x = np.log10(P_Range)            
y = np.log10(Avg_L)
x_Fit = np.log10(P_Fit)
y_Fit = np.log10(L_Fit)

sig_L=np.zeros(len(y_Fit))   #change y_Fit to Avg_L for the error on all data points
for c in xrange(0,len(y_Fit)):
    for i in xrange(0,len(count_L)):
        if count_L[c] >3:
            sig_L[c] += ((np.log10(lattice.length_loop[i])-np.log10(L_Fit[c]))**2)
    if count_L[c] >3:
        sig_L[c] = np.sqrt((1./(count_L[c]-1))*sig_L[c])/np.sqrt(count_L[c])
#sig_L[11]+= 0.8   # because at box size N = 40, sig_L has size (y_Fit -1) and this would be zero
#sig_L[14]+=0.9      #for N=55
#sig_L[26]+= 1.2        #for N=90
                                                                                                                                    
plt.figure("Fig.4")
plt.scatter(x, y)
plt.errorbar(x_Fit, y_Fit,xerr = 0, yerr = sig_L, fmt ='o', c = 'blue')
popt2,pcov2 = curve_fit(lin_func, x_Fit, y_Fit, sigma = sig_L)
x_lin_Fit = np.linspace(min(x_Fit),max(x_Fit),500)
plt.plot( x_lin_Fit, lin_func(x_lin_Fit,*popt2) , c = 'blue')
error2 = np.sqrt(np.diag(pcov2))
popt2[0] = 10**(popt2[0])
error2[0] = error2[0]*(np.log(10))*popt2[0]
error2[1] = error2[1]
plt.xlabel(r'$Log(Loop \ Perimeter)$', size = '18')
plt.ylabel(r'$Log(< Length >)$', size = '18')
plt.title(r'$Estimation \ of \ the \ fractal \ dimension$', size = '18')
plt.show("Fig.4")
print "[Fig.4] A = %.3f" %(popt2[0]), "+/- %.3f" %(error2[0] )
print "[Fig.4] d = %.3f" %(popt2[1]), "+/- %.3f" %(error2[1])

n = []
L_Range = []
P_Range = []  #redefine  P_Range if cutting at R_c for in the plot above
count_P = []
P_Fit = []
L_Fit = []
n_Fit = []
for i in xrange(5, max(lattice.size_loop)):
    select = (A[:,1] == i)
    if (len(A[select,0])!=0) and (A[select,1][0]>5):
        n.append(1.0*len(A[select,1])/(N**3))
        L_Range.append(A[select,0][0])
        P_Range.append(A[select,1][0])
        count_P.append(len(A[select,1]))
        if(len(A[select,0])>4):
            P_Fit.append(A[select,1][0])
            L_Fit.append(A[select,0][0])
            n_Fit.append(1.0*len(A[select,1])/(N**3))

                        
#l = P_Range            #Alternative method for power law parameter estimation
#beta = 0
#for i in xrange(0, len(l)):
#    beta += (np.log(l[i]/min(l)))
#beta = 1 + (len(l)/beta)

x = np.log10(P_Range)
y = np.log10(n)
x_Fit = np.log10(P_Fit)
y_Fit = np.log10(n_Fit)

sig_P=np.zeros(len(y_Fit))  
for c in xrange(0,len(y_Fit)):
    for i in xrange(0,len(count_P)):
        if count_P[c] >3:
            sig_P[c] += ((np.log10(lattice.size_loop[i])-np.log10(P_Fit[c]))**2) 
    if count_P[c] >3:
        sig_P[c] = np.sqrt((1./(count_P[c]-1))*sig_P[c])/np.sqrt(count_P[c])
sig_n = (4.0*sig_P)/(x_Fit**5)
#sig_n[11]+= 0.7  # because at box size N = 40, this would be zero
#sig_n[14]+=0.9      #for N=55
#sig_n[20]+=0.4         #for N=90

plt.figure("Fig.5")
plt.scatter(x, y)
plt.errorbar(x_Fit, y_Fit ,xerr = 0, yerr = sig_n, fmt ='o', c = 'blue')
popt3,pcov3 = curve_fit(lin_func, x_Fit, y_Fit, sigma = sig_n)
x_lin_Fit = np.linspace(min(x_Fit),max(x_Fit),500)
plt.plot( x_lin_Fit, lin_func(x_lin_Fit,*popt3) , c = 'blue')
error3 = np.sqrt(np.diag(pcov3))
popt3[0] = 10**(popt3[0])
error3[0] = error3[0]*(np.log(10))*popt3[0]
error3[1] = error3[1]
plt.xlabel(r'$Log(Loop \ Perimeter)$', size = '18')
plt.ylabel(r'$Log(Density)$', size = '18')
plt.show("Fig.5")
print "[Fig.5] B = %.3f" %(popt3[0]), "+/- %.3f" %(error3[0] )
print "[Fig.5] b = %.3f" %(-popt3[1]), "+/- %.3f" %(error3[1])

x = np.log10(L_Range)
y = np.log10(n)
x_Fit = np.log10(L_Fit)
y_Fit = np.log10(n_Fit)

sig_n2 = (sig_L*(5/2))/(x_Fit**(3/2))   

plt.figure("Fig.6")
plt.scatter(x, y)
plt.errorbar(x_Fit, y_Fit ,xerr = 0, yerr = sig_n2, fmt ='o', c = 'blue')
popt4,pcov4 = curve_fit(lin_func, x_Fit, y_Fit, sigma = sig_n2)
x_lin_Fit = np.linspace(min(x_Fit),max(x_Fit),500)
plt.plot( x_lin_Fit, lin_func(x_lin_Fit,*popt4) , c = 'blue')
error4 = np.sqrt(np.diag(pcov4))
popt4[0] = 10**(popt4[0])
error4[0] = error4[0]*(np.log(10))*popt4[0]
error4[1] = error4[1]
plt.xlabel(r'$Log(Loop \ Length)$', size = '18')
plt.ylabel(r'$Log(Density)$', size = '18')
plt.show("Fig.6")
print "[Fig.6] C = %.3f" %(popt4[0]), "+/- %.3f" %(error4[0] )
print "[Fig.6] g = %.3f" %(-popt4[1]), "+/- %.3f" %(error4[1])

n_open = []
Inf_Range = []
M =np.zeros((len(lattice.length_inf),2))
M[:,0] += lattice.length_inf
nopen_Fit =[]
Inf_Fit = []
count_inf =[]
for i in xrange(4, max(lattice.length_inf)):
    select = (M[:,0] == i)
    if (len(M[select,0])!=0):
        n_open.append(1.0*len(M[select,0])/(N**3))
        Inf_Range.append(M[select,0][0])
        count_inf.append(len(M[select,0]))
        if(len(M[select,0])>4):
                Inf_Fit.append(M[select,0][0])
                nopen_Fit.append(1.0*len(M[select,0])/(N**3))
            
            
x = np.log10(Inf_Range)
y = np.log10(n_open)
y_Fit = np.log10(nopen_Fit)
x_Fit = np.log10(Inf_Fit)

sig_inf=np.zeros(len(y_Fit))  
for c in xrange(0,len(y_Fit)):
    for i in xrange(0,len(count_inf)):
        if count_inf[c] >3:
            sig_inf[c] += ((np.log10(lattice.length_inf[i])-np.log10(Inf_Fit[c]))**2) 
    if count_inf[c] >3:
        sig_inf[c] = np.sqrt((1./(count_inf[c]-1))*sig_inf[c])/np.sqrt(count_inf[c])
#sig_inf = ??

plt.figure("Fig.7")
plt.scatter(x, y)
#plt.errorbar(x_Fit, y_Fit ,xerr = 0, yerr = sig_inf, fmt ='o', c = 'blue')
popt5,pcov5 = curve_fit(lin_func, x_Fit, y_Fit) #, sigma = sig_inf)
x_lin_Fit = np.linspace(min(x_Fit),max(x_Fit),500)
plt.plot( x_lin_Fit, lin_func(x_lin_Fit,*popt5) , c = 'blue')
error5 = np.sqrt(np.diag(pcov5))
popt5[0] = 10**(popt5[0])
error5[0] = error5[0]*(np.log(10))*popt5[0]
error5[1] = error5[1]
plt.xlabel(r'$Log(Segment \ Length)$', size = '18')
plt.ylabel(r'$Log(Density)$', size = '18')
plt.show("Fig.7")
print "[Fig.7] D = %.3f" %(popt5[0]), "+/- %.3f" %(error5[0] )
print "[Fig.7] d = %.3f" %(-popt5[1]), "+/- %.3f" %(error5[1])

x = []
y = []
x_Fit = []
y_Fit = []
for i in xrange(0, len(lattice.VS_ratio)):    
    if lattice.VS_ratio[i] != 0:
        x.append(np.log10(lattice.size_loop[i]))
        y.append(np.log10(lattice.VS_ratio[i]))
    #if lattice.VS_ratio[i] != 0 and (lattice.size_loop[i] > 20):   #Comment out if want to include only max value of VS_ratio
    #    x_Fit.append(np.log10(lattice.size_loop[i]))
    #    y_Fit.append(np.log10(lattice.VS_ratio[i]))  
        
V =np.zeros((len(lattice.size_loop), 2)) 
V[:,0] += lattice.size_loop
V[:,1] += lattice.VS_ratio        
for i in xrange(5, max(lattice.size_loop)+1):  #Including only max values of VS_ratio, comment out for previous x_Fit, y_Fit
    select = (V[:,0] == i)
    if len(V[select,0])!=0:
        if (V[select,0][0] > 10) and len(V[select,1])!=0:
            y_Fit.append(np.log10(max(V[select, 1])))
            x_Fit.append(np.log10(V[select,0][0]))
  
#VS_data_150 = np.loadtxt("VS_data_150.txt")   #VS plot for box size=150
#x = []
#y = []
#x_Fit = []
#y_Fit = []
#size_150 = VS_data_150[:,0]
#VS_150 = VS_data_150[:,1]
#for i in xrange(0, len(VS_150)):           
#    if VS_150[i] != 0:
#        x.append(np.log10(size_150[i]))
#        y.append(np.log10(VS_150[i]))
##    if VS_150[i] != 0 and (size_150[i] > 30):      #Comment out if want to include only max value of VS_ratio
##        x_Fit.append(np.log10(size_150[i]))
##        y_Fit.append(np.log10(VS_150[i]))    
#for i in xrange(5, int(max(size_150)+1)):  #Including only max values of VS_ratio, comment out for previous x_Fit, y_Fit
#    select = (VS_data_150[:,0] == i)
#    if len(VS_data_150[select,0])!=0:
#        if (VS_data_150[select,0][0] > 30) and len(VS_data_150[select,1])!=0:
#            y_Fit.append(np.log10(max(VS_data_150[select, 1])))
#            x_Fit.append(np.log10(VS_data_150[select,0][0])) 

plt.figure("Fig.V/S")
plt.scatter(x, y)
poptVS,pcovVS = curve_fit(lin_func, x_Fit, y_Fit)
x_lin = np.linspace(min(x_Fit),max(x_Fit),1000)
plt.plot( x_lin, lin_func(x_lin,*poptVS) , c = 'blue')
plt.xlabel(r'$Log(Loop \ Perimeter)$', size = '18')
plt.ylabel(r'$Log(Volume \ to \ Surface \ ratio)$', size = '18')
plt.show("Fig.V/S")
errorVS = np.sqrt(np.diag(pcovVS))
poptVS[0] = 10**(poptVS[0])
errorVS[0] = errorVS[0]*(np.log(10))*poptVS[0]
errorVS[1] = errorVS[1]
print "[Fig.V/S] K = %.3f" %(poptVS[0]), "+/- %.3f" %(errorVS[0] )
print "[Fig.V/S] v = %.3f" %(poptVS[1]), "+/- %.3f" %(errorVS[1])
#np.savetxt("VS_data_150.txt", np.c_[lattice.size_loop,lattice.VS_ratio], fmt ='%0.6f') #For when N=150 -> save this raw data

print "Fraction of the lenght of open strings", 1.0*(np.sum(lattice.length_inf))/(np.sum(lattice.length_inf)+np.sum(lattice.length_loop))
L_Frac = 1.0*(np.sum(lattice.length_inf))/(np.sum(lattice.length_inf)+np.sum(lattice.length_loop)) 

#np.savetxt("gradients_14.txt", np.c_[popt1[1],popt2[1],popt3[1],popt4[1],poptVS[1],popt5[1],L_Frac], fmt ='%0.6f')

Size_14 = np.loadtxt("gradients_14.txt") #new
Size_15 = np.loadtxt("gradients_15.txt")
Size_18 = np.loadtxt("gradients_18.txt") #new
Size_19 = np.loadtxt("gradients_19.txt") #new
Size_20 = np.loadtxt("gradients_20.txt")
Size_22 = np.loadtxt("gradients_22.txt") #new
Size_24 = np.loadtxt("gradients_24.txt") #new
Size_25 = np.loadtxt("gradients_25.txt")
Size_30 = np.loadtxt("gradients_30.txt")
Size_35 = np.loadtxt("gradients_35.txt")
Size_40 = np.loadtxt("gradients_40.txt")
Size_45 = np.loadtxt("gradients_45.txt")
Size_50 = np.loadtxt("gradients_50.txt")
Size_55 = np.loadtxt("gradients_55.txt")
Size_60 = np.loadtxt("gradients_60.txt")
Size_65 = np.loadtxt("gradients_65.txt")
Size_70 = np.loadtxt("gradients_70.txt")
Size_80 = np.loadtxt("gradients_80.txt")
Size_90 = np.loadtxt("gradients_90.txt")
Size_100 = np.loadtxt("gradients_100.txt")
Size_115 = np.loadtxt("gradients_115.txt") #new
Size_125 = np.loadtxt("gradients_125.txt") #new

grad_d1 = [Size_14[0],Size_15[0],Size_18[0],Size_19[0],Size_20[0],Size_22[0],Size_24[0],Size_25[0],Size_30[0],Size_35[0],Size_40[0],Size_45[0],Size_50[0],Size_55[0],Size_60[0],Size_65[0],Size_70[0],Size_80[0],Size_90[0],Size_100[0],Size_115[0],Size_125[0]]
grad_d2 = [Size_14[1],Size_15[1],Size_18[1],Size_19[1],Size_20[1],Size_22[1],Size_24[1],Size_25[1],Size_30[1],Size_35[1],Size_40[1],Size_45[1],Size_50[1],Size_55[1],Size_60[1],Size_65[1],Size_70[1],Size_80[1],Size_90[1],Size_100[1],Size_115[1],Size_125[1]]
grad_b = [Size_14[2],Size_15[2],Size_18[2],Size_19[2],Size_20[2],Size_22[2],Size_24[2],Size_25[2],Size_30[2],Size_35[2],Size_40[2],Size_45[2],Size_50[2],Size_55[2],Size_60[2],Size_65[2],Size_70[2],Size_80[2],Size_90[2],Size_100[2],Size_115[2],Size_125[2]]
grad_g = [Size_14[3],Size_15[3],Size_18[3],Size_19[3],Size_20[3],Size_22[3],Size_24[3],Size_25[3],Size_30[3],Size_35[3],Size_40[3],Size_45[3],Size_50[3],Size_55[3],Size_60[3],Size_65[3],Size_70[3],Size_80[3],Size_90[3],Size_100[3],Size_115[3],Size_125[3]]
grad_vs = [Size_14[4],Size_15[4],Size_18[4],Size_19[4],Size_20[4],Size_22[4],Size_24[4],Size_25[4],Size_30[4],Size_35[4],Size_40[4],Size_45[4],Size_50[4],Size_55[4],Size_60[4],Size_65[4],Size_70[4],Size_80[4],Size_90[4],Size_100[4],Size_115[4],Size_125[4]]
grad_delta = [Size_14[5],Size_15[5],Size_18[5],Size_19[5],Size_20[5],Size_22[5],Size_24[5],Size_25[5],Size_30[5],Size_35[5],Size_40[5],Size_45[5],Size_50[5],Size_55[5],Size_60[5],Size_65[5],Size_70[5],Size_80[5],Size_90[5],Size_100[5],Size_115[5],Size_125[5]]
grad_L_Frac = [Size_14[6],Size_15[6],Size_18[6],Size_19[6],Size_20[6],Size_22[6],Size_24[6],Size_25[6],Size_30[6],Size_35[6],Size_40[6],Size_45[6],Size_50[6],Size_55[6],Size_60[6],Size_65[6],Size_70[6],Size_80[6],Size_90[6],Size_100[6],Size_115[6],Size_125[6]]
print grad_L_Frac
x = [14,15,18,19,20,22,24,25,30,35,40,45,50,55,60,65,70,80,90,100,115,125]
plt.figure("Fig.Param")
plt.scatter(x, grad_vs) #, label = 'average length vs \n loop perimeter')
plt.title("Parameter $\it{f}_{open}$ as a function of box size")
#plt.legend(loc = 4,prop={'size': 16})
ax = plt.axes()
ax.set_xticks([15,20,25,30,35,40,45,50,55,60,65,70,80,90,100, 125])
plt.show("Fig.Param")
