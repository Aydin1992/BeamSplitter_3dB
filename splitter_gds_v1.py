# -- coding utf-8 --

Created on Tue May 13 122053 2025

@author Dell


import numpy as np
import matplotlib.pyplot as plt
import gdspy
import cv2
import time
from shapely.geometry import Polygon
import os

timestamp=time.strftime(%Y%m%d_%H%M%S, time.localtime())

class splitter_3dB()
    def __init__(self, Lwg,Wg, Ltaper,Wtaper,L_MMI,W_MMI, WG_gap,Edge_gap,Lcoupler)
        self.Lwg=Lwg
        self.Wg=Wg
        self.Ltaper=Ltaper
        self.Wtaper=Wtaper
        self.L_MMI=L_MMI
        self.W_MMI=W_MMI
        self.WG_gap=WG_gap
        self.Edge_gap=Edge_gap  
        self.Lcoupler=Lcoupler  
    def MMI_Footprint_Polygon(self,shift_x,shift_y)
        Wg_left_up=[]
        Wg_left_down=[]        
        taper_left_up=[]
        taper_left_down=[]        
        MMI=[]
        Wg_right_up=[]
        Wg_right_down=[]
        taper_right_up=[]
        taper_right_down=[]
        
        Wcoupler=0.18
        Wc_MMI=np.max([self.W_MMI,2self.Edge_gap+self.WG_gap+2self.Wtaper])
        
        y0_up=self.WG_gap2+self.Wtaper2+shift_y
        y0_down=-(y0_up-shift_y)+shift_y
        y0_middle=(y0_up+y0_down)2
        x0=0+shift_x 
        
        y0_coupler_up=y0_middle + 1252
        y0_coupler_down=y0_middle - 1252
        
        path1=gdspy.Path(1,(0,0))
        
        taper_left_up.extend([ (x0,y0_coupler_up-Wcoupler2) , (x0,y0_coupler_up+Wcoupler2),
                              (x0+Lcoupler,y0_coupler_up+Wg2),(x0+Lcoupler,y0_coupler_up-Wg2) ,
                              (x0,y0_coupler_up-Wcoupler2)])
        
        taper_left_down.extend([ (x0,y0_coupler_down-Wcoupler2) , (x0,y0_coupler_down+Wcoupler2),
                              (x0+Lcoupler,y0_coupler_down+Wg2),(x0+Lcoupler,y0_coupler_down-Wg2) ,
                              (x0,y0_coupler_down-Wcoupler2)])
        
        
        Wg_left_up.extend([(x0+self.Lcoupler,y0_up-self.Wg2),(x0+self.Lcoupler,y0_up+self.Wg2),
                           (x0+self.Lcoupler+self.Lwg,y0_up+self.Wg2),
                           (x0+self.Lcoupler+self.Lwg+self.Ltaper,y0_up+self.Wtaper2),
                           (x0+self.Lcoupler+self.Lwg+self.Ltaper,y0_up-self.Wtaper2),
                           (x0+self.Lcoupler+self.Lwg,y0_up-self.Wg2),(x0+self.Lcoupler,y0_up-self.Wg2)])
        
        Wg_left_down.extend([(x0+self.Lcoupler,y0_down-self.Wg2),(x0+self.Lcoupler,y0_down+self.Wg2),
                             (x0+self.Lcoupler+self.Lwg,y0_down+self.Wg2),
                   (x0+self.Lcoupler+self.Lwg+self.Ltaper,y0_down+self.Wtaper2),
                   (x0+self.Lcoupler+self.Lwg+self.Ltaper,y0_down-self.Wtaper2),
                   (x0+self.Lcoupler+self.Lwg,y0_down-self.Wg2),(x0+self.Lcoupler,y0_down-self.Wg2)])
        
        MMI.extend([(x0+self.Lcoupler+self.Lwg+self.Ltaper,y0_middle+Wc_MMI2),
                    (x0+self.Lcoupler+self.Lwg+self.Ltaper+self.L_MMI,y0_middle+Wc_MMI2),
                    (x0+self.Lcoupler+self.Lwg+self.Ltaper+self.L_MMI,y0_middle-Wc_MMI2),
                    (x0+self.Lcoupler+self.Lwg+self.Ltaper,y0_middle-Wc_MMI2)])
                
        x0_right= x0+self.Lwg+self.Ltaper+self.L_MMI+self.Lcoupler
        
        Wg_right_up.extend([(x0_right,y0_up-self.Wtaper2),(x0_right,y0_up+self.Wtaper2),
                            (x0_right+self.Ltaper,y0_up+self.Wg2),
                   (x0_right+self.Lwg+self.Ltaper,y0_up+self.Wg2),
                   (x0_right+self.Lwg+self.Ltaper,y0_up-self.Wg2),
                   (x0_right+self.Ltaper,y0_up-self.Wg2),(x0_right,y0_up-self.Wtaper2)])
        
        Wg_right_down.extend([(x0_right,y0_down-self.Wtaper2),(x0_right,y0_down+self.Wtaper2),
                              (x0_right+self.Ltaper,y0_down+self.Wg2),
                   (x0_right+self.Lwg+self.Ltaper,y0_down+self.Wg2),
                   (x0_right+self.Lwg+self.Ltaper,y0_down-self.Wg2),
                   (x0_right+self.Ltaper,y0_down-self.Wg2),(x0_right,y0_down-self.Wtaper2)])
        return taper_left_up,taper_left_down,Wg_left_up, Wg_left_down,MMI,Wg_right_up,Wg_right_down
    
                                    
# GDS file setting
cell_name = '3dB_Splitter'
gds_file_name ='Splitter.gds'
folder_name=Splitter_MaskFile  
gds_layer = 1
box_layer = 82
  
#### Preparing the Folder
if not os.path.exists(folder_name)
    os.makedirs(folder_name)   
os.chdir(folder_name)     
gdspy.current_library = gdspy.GdsLibrary()  
cell = gdspy.Cell(cell_name)  
###############################################################################
################################ Parameters ###################################
Lwg=10
Wg=0.5
Ltaper=30
Wtaper=5
L_MMI=60
W_MMI=1
WG_gap=3
Edge_gap=0
Lcoupler=200

MMI_3dB=splitter_3dB(Lwg,Wg, Ltaper,Wtaper,L_MMI,W_MMI, WG_gap,Edge_gap,Lcoupler)
components = MMI_3dB.MMI_Footprint_Polygon(shift_x=0,shift_y=0);
for pol in components
    polygon = gdspy.Polygon(pol, layer=gds_layer, datatype=0)
    cell.add(polygon)                         
###############################################################################
#################### add box (design domain) to GDS file ######################
# box = [(np.min(x), np.min(y)), (np.min(x), np.max(y)), (np.max(x), np.max(y)), (np.max(x), np.min(y)),
#                (np.min(x), np.min(y))]
# polygon = gdspy.Polygon(box, layer=box_layer, datatype=0)
# cell.add(polygon)
###############################################################################
############################### write to gds file #############################
lib = gdspy.GdsLibrary()
lib.add(cell)
lib.write_gds(gds_file_name)