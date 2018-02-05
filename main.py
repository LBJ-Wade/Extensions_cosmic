"""
RHUL PH4100 - Major Project - EXTENSIONS!!!
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
    plt.xlabel(r'$Length \ of \ Strings$', fontsize=22)
    plt.ylabel(r'$Number \ of \ Strings$', fontsize=22)
    plt.title(r'$Histogram \ of \ String \ Lengths$', fontsize=25, y=1.025)
    plt.legend(loc='upper right', fontsize=25)
    plt.annotate(r'$Size: \ N\xi \ = \ {0}$'.format(N), xy=(800, 325), xycoords='figure points', fontsize=22,
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
                     box[i,j,k] = randint(0, 3) ###Change
        total =0 
        faceNum=0
        edge = False
        L=0
        count=np.zeros(10-1)
        sum_e2e=np.zeros(10-1)
        e2e=[[],[],[],[],[],[],[],[],[]]
        string_coords=[] 
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
        if (np.mod(face[3] - face[0],3) == 1): #test right-most and left-most values      ###Change                      
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
                
                          
    def follow(self,xyz_string,i,j,k,XYZ): 
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

