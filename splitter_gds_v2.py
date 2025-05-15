# -*- coding: utf-8 -*-
"""
Created on Tue May 13 12:20:53 2025

@author: Dell
"""

import numpy as np
import matplotlib.pyplot as plt
import gdspy
import cv2
import time
from shapely.geometry import Polygon
import os

timestamp=time.strftime("%Y%m%d_%H%M%S", time.localtime())

class splitter_3dB():
    def __init__(self, Lwg,Wg, Ltaper,Wtaper,L_MMI,W_MMI, WG_gap,Edge_gap,Lcoupler,Bend,
                 wfiber=125,Ltotal=518):
        self.Lwg=Lwg
        self.Wg=Wg
        self.Ltaper=Ltaper
        self.Wtaper=Wtaper
        self.L_MMI=L_MMI
        self.W_MMI=W_MMI
        self.WG_gap=WG_gap
        self.Edge_gap=Edge_gap  
        self.Lcoupler=Lcoupler  
        self.Bend=Bend
        self.wfiber=wfiber
        self.Ltotal=Ltotal
    def MMI_Footprint(self,shift_x,shift_y):
        ## Parameters:
        wcoupler=0.18
        Lcoupler=(self.Ltotal-self.L_MMI)/2 - (2*self.Bend+self.Lwg+self.Ltaper)
        Wc_MMI=np.max([self.W_MMI,2*self.Edge_gap+self.WG_gap+2*self.Wtaper])        
        y0_up=self.WG_gap/2+self.Wtaper/2+shift_y
        y0_down=-(y0_up-shift_y)+shift_y
        y0_middle=(y0_up+y0_down)/2
        x0=shift_x  - (self.Lwg+self.Ltaper+self.L_MMI/2)
        y0_coupler=y0_middle+self.wfiber/2
        dy_bend=(y0_coupler-self.Bend)-(y0_up+self.Bend)
                
        WG_Left_Up=gdspy.Path(self.Wtaper, (x0+self.Lwg+self.Ltaper, y0_up))
        WG_Left_Up.segment(self.Ltaper,"-x",layer=gds_layer, datatype=0,final_width=self.Wg)
        WG_Left_Up.segment(self.Lwg,"-x",layer=gds_layer, datatype=0)
        WG_Left_Up.turn(self.Bend, angle=-np.pi/2,layer=gds_layer, datatype=0)
        WG_Left_Up.segment(dy_bend,"+y",layer=gds_layer, datatype=0)
        WG_Left_Up.turn(self.Bend, angle=np.pi/2,layer=gds_layer, datatype=0)
        WG_Left_Up.segment(Lcoupler,"-x",layer=gds_layer, datatype=0,final_width=wcoupler)

        WG_Left_Down=gdspy.Path(self.Wtaper, (x0+self.Lwg+self.Ltaper, y0_down))
        WG_Left_Down.segment(self.Ltaper,"-x",layer=gds_layer, datatype=0,final_width=self.Wg)
        WG_Left_Down.segment(self.Lwg,"-x",layer=gds_layer, datatype=0)
        WG_Left_Down.turn(self.Bend, angle=np.pi/2,layer=gds_layer, datatype=0)
        WG_Left_Down.segment(dy_bend,"-y",layer=gds_layer, datatype=0)
        WG_Left_Down.turn(self.Bend, angle=-np.pi/2,layer=gds_layer, datatype=0)
        WG_Left_Down.segment(Lcoupler,"-x",layer=gds_layer, datatype=0,final_width=wcoupler)

        x0_right= x0+self.Lwg+self.Ltaper+self.L_MMI
        
        WG_right_Up=gdspy.Path(self.Wtaper, (x0_right, y0_up))
        WG_right_Up.segment(self.Ltaper,"+x",layer=gds_layer, datatype=0,final_width=self.Wg)
        WG_right_Up.segment(self.Lwg,"+x",layer=gds_layer, datatype=0)
        WG_right_Up.turn(self.Bend, angle=np.pi/2,layer=gds_layer, datatype=0)
        WG_right_Up.segment(dy_bend,"+y",layer=gds_layer, datatype=0)
        WG_right_Up.turn(self.Bend, angle=-np.pi/2,layer=gds_layer, datatype=0)
        WG_right_Up.segment(Lcoupler,"+x",layer=gds_layer, datatype=0,final_width=wcoupler)
        
        WG_right_Down=gdspy.Path(self.Wtaper, (x0_right, y0_down))
        WG_right_Down.segment(self.Ltaper,"+x",layer=gds_layer, datatype=0,final_width=self.Wg)
        WG_right_Down.segment(self.Lwg,"+x",layer=gds_layer, datatype=0)
        WG_right_Down.turn(self.Bend, angle=-np.pi/2,layer=gds_layer, datatype=0)
        WG_right_Down.segment(dy_bend,"-y",layer=gds_layer, datatype=0)
        WG_right_Down.turn(self.Bend, angle=np.pi/2,layer=gds_layer, datatype=0)
        WG_right_Down.segment(Lcoupler,"+x",layer=gds_layer, datatype=0,final_width=wcoupler)
        
        MMI=gdspy.Rectangle((shift_x-self.L_MMI/2,y0_middle-self.W_MMI/2),(shift_x+self.L_MMI/2,y0_middle+self.W_MMI/2),layer=gds_layer, datatype=0)
        Wtotal=self.wfiber+self.Wg
        
        print('total length %G and width %G of the device' %(Ltotal,Wtotal))
        
        return WG_Left_Up,WG_Left_Down,WG_right_Up,WG_right_Down,MMI
    
    def Straight_WG(self,shift_x,shift_y):
        wcoupler=0.18
        Lcoupler=(self.Ltotal-self.L_MMI)/2 - (2*self.Bend+self.Lwg+self.Ltaper)
        Wc_MMI=np.max([self.W_MMI,2*self.Edge_gap+self.WG_gap+2*self.Wtaper])        
        y0_up=self.WG_gap/2+self.Wtaper/2+shift_y
        y0_down=-(y0_up-shift_y)+shift_y
        y0_middle=(y0_up+y0_down)/2
        x0=shift_x  - (self.Lwg+self.Ltaper+self.L_MMI/2)
        y0_coupler=y0_middle+self.wfiber/2
        dy_bend=(y0_coupler-self.Bend)-(y0_up+self.Bend)
        
        
        xlef=shift_x-self.L_MMI/2-self.Ltaper-self.Lwg-2*self.Bend-Lcoupler
        xright=shift_x-(-self.L_MMI/2-self.Ltaper-self.Lwg-2*self.Bend-Lcoupler)
        
        straight_wg_left=gdspy.Path(wcoupler, (xlef, y0_middle))
        straight_wg_left.segment(Lcoupler,"+x",layer=gds_layer, datatype=0,final_width=self.Wg)
        
        straight_wg_right=gdspy.Path(wcoupler, (xright, y0_middle))
        straight_wg_right.segment(Lcoupler,"-x",layer=gds_layer, datatype=0,final_width=self.Wg)
        
        
        dL=(xright-Lcoupler)-(xlef+Lcoupler)
        straight_wg_mid=gdspy.Path(Wg, (xlef+Lcoupler, y0_middle))
        straight_wg_mid.segment(dL,"+x",layer=gds_layer, datatype=0)
 
        return straight_wg_left,straight_wg_right,straight_wg_mid
    
    
    def Bend_WG(self,shift_x,shift_y):
        wcoupler=0.18
        Lcoupler=(self.Ltotal-self.L_MMI)/2 - (2*self.Bend+self.Lwg+self.Ltaper)
        Wc_MMI=np.max([self.W_MMI,2*self.Edge_gap+self.WG_gap+2*self.Wtaper])        
        y0_up=self.WG_gap/2+self.Wtaper/2+shift_y
        y0_down=-(y0_up-shift_y)+shift_y
        y0_middle=(y0_up+y0_down)/2
        x0=shift_x  - (self.Lwg+self.Ltaper+self.L_MMI/2)
        y0_coupler=y0_middle+self.wfiber/2
        dy_bend=(y0_coupler-self.Bend)-(y0_middle+self.Bend)
        
        WG_Left_Up=gdspy.Path(self.Wg, (shift_x, y0_middle))
        WG_Left_Up.segment(self.Lwg+self.Ltaper+self.L_MMI/2,"-x",layer=gds_layer, datatype=0)
        WG_Left_Up.turn(self.Bend, angle=-np.pi/2,layer=gds_layer, datatype=0)
        WG_Left_Up.segment(dy_bend,"+y",layer=gds_layer, datatype=0)
        WG_Left_Up.turn(self.Bend, angle=np.pi/2,layer=gds_layer, datatype=0)
        WG_Left_Up.segment(Lcoupler, "-x",layer=gds_layer, datatype=0,final_width=wcoupler)
        
            
        WG_right_Up=gdspy.Path(self.Wg, (shift_x, y0_middle))
        WG_right_Up.segment(self.Lwg+self.Ltaper+self.L_MMI/2+2*self.Bend+Lcoupler,"+x",layer=gds_layer, datatype=0,final_width=wcoupler)
        

        

        return WG_Left_Up,WG_right_Up




                                        
# GDS file setting:
cell_name = 'Splitter_3dB'
lib_name=cell_name
gds_file_name ='Splitter.gds'
folder_name="Splitter_MaskFile"  
gds_layer = 10
box_layer = 82
  
#### Preparing the Folder:
if not os.path.exists(folder_name):
    os.makedirs(folder_name)   
os.chdir(folder_name)     
gdspy.current_library = gdspy.GdsLibrary()  
cell = gdspy.Cell(cell_name)  
###############################################################################
################################ Parameters ###################################
Lwg=10
Wg=0.5
Ltaper=20
Wtaper=1
L_MMI=18
W_MMI=4
WG_gap=Wg
Edge_gap=0.8
Lcoupler=200
Bend=10
circuit_gap=20
Ltotal=1000
wfiber=125
FP_xgap=2
FP_ygap=50
# length variations:   
n=0
for i in np.arange(L_MMI-2,L_MMI+2,1):     
    Splitter=splitter_3dB(Lwg,Wg, Ltaper,Wtaper,i,W_MMI, WG_gap,Edge_gap
                       ,Lcoupler,Bend,wfiber,Ltotal).MMI_Footprint(shift_x=0,shift_y=(wfiber+circuit_gap)*n) 
    cell.add(Splitter)  
    n+=1
# gap variations:   
m=0
for i in np.arange(WG_gap-0.3,WG_gap+0.3,0.1):     
    Splitter=splitter_3dB(Lwg,Wg, Ltaper,Wtaper,L_MMI,W_MMI, i,Edge_gap
                       ,Lcoupler,Bend,wfiber,Ltotal).MMI_Footprint(shift_x=0,shift_y=(wfiber+circuit_gap)*(n+m)) 
    cell.add(Splitter)  
    m+=1
   
Test_Devices1 = splitter_3dB(Lwg,Wg, Ltaper,Wtaper,L_MMI,W_MMI, WG_gap,Edge_gap
                   ,Lcoupler,Bend,wfiber,Ltotal).Straight_WG(shift_x=0,shift_y=(wfiber+circuit_gap)*(n+m-1)+wfiber/2+circuit_gap)

cell.add(Test_Devices1)

Test_Devices2 = splitter_3dB(Lwg,Wg, Ltaper,Wtaper,L_MMI,W_MMI, WG_gap,Edge_gap
                   ,Lcoupler,Bend,wfiber,Ltotal).Bend_WG(shift_x=0,shift_y=(wfiber+circuit_gap)*(n+m-1)+wfiber/2+circuit_gap*2)
 
cell.add(Test_Devices2)                     
###############################################################################
#################### add box (design domain) to GDS file ######################

y_up=(wfiber+circuit_gap)*(n+m-1)+wfiber+circuit_gap*2+FP_ygap+Wg/2
y_dowm=-wfiber/2-FP_ygap-Wg/2

FP = gdspy.Rectangle((-FP_xgap-Ltotal/2,y_dowm),(Ltotal/2+FP_xgap,y_up),layer=box_layer, datatype=0)
cell.add(FP)
###############################################################################
############################### write to gds file #############################
lib = gdspy.GdsLibrary(name=lib_name,unit=1e-6,precision=1e-9)
lib.add(cell)
lib.write_gds(gds_file_name)











