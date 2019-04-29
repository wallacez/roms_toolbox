#
## 3d_gui.py
## 3D visualization of 2 variables in the Patagonia region 
#

import sys, os, random
import matplotlib
import interp    # f2py generated fortran wrappers from interpolate.f90
import numpy as np
import math
import re

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from netCDF4 import Dataset, Variable

from collections import OrderedDict
from datetime    import *
from time        import *


# Global Variables
global exp_dir, file_list, fhead
global tday_ref, date_ref, date_file_1
global fnum_min, fnum_max
global varnamelist, var_1_combined, var_2_combined
global F1, F2, zr, nt, var_1_name, var_2_name
global it, nt, ii1, jj1, int_i, int_j
global lon_rho, lat_rho, mask_rho, lon_u, lat_u, mask_u, lon_v, lat_v, mask_v
global lon_psi, lat_psi, mask_psi, h
global xi_rho, eta_rho, xi_psi, eta_psi 
global xi_u, eta_u, xi_v, eta_v
global lon_1, lon_2, lat_1, lat_2
global mask_rho_nan
global var_1_loaded_from, var_2_loaded_from

## uncomment the below to get information on fortran routines used 
#
##############################
#                            #
# print im.hindices.__doc__  #
# print im.try_range.__doc__ #
# print im.inside.__doc__    #
#                            #
##############################


##### matplotlibrc settings #####

# tick label size        
matplotlib.rc('xtick', labelsize=7)
matplotlib.rc('ytick', labelsize=7)

# title size
matplotlib.rc('axes', titlesize=10)

# font size
matplotlib.rc('font', size=10.0)

##### End matplotlibrc settings #####


class AppForm(QMainWindow):
    
    def __init__(self, parent=None):
        '''Function to instantiate the GUI window object and widgets inside
        '''
        
        QMainWindow.__init__(self, parent)
        
        self.setWindowTitle('3D Patagonia Model')
        
        ### read grid data in from grid file
        self.read_grid_data()
        
        ### create window widgets
        self.create_plot_window()
        
 
    def create_plot_window(self):
        '''Function to create UI
        '''
        
        global varnamelist, fhead
        
        self.plot_window = QWidget()

        ### create two figures
        self.dpi      = 150
        self.fig_1    = Figure((10.0, 10.0), dpi=self.dpi)
        self.fig_2    = Figure((10.0, 10.0), dpi=self.dpi)
        self.canvas_1 = FigureCanvas(self.fig_1)
        self.canvas_2 = FigureCanvas(self.fig_2)
        self.canvas_1.setParent(self.plot_window)
        self.canvas_2.setParent(self.plot_window)
        
        ### create 2 subplots
        self.axes_1 = self.fig_1.add_subplot(111)
        self.axes_2 = self.fig_2.add_subplot(111)
        
        ### GUI controls ###
        
        #### input fields (choose model, date, lon/lat, variables)
        self.model_select_label = QLabel('Select Model')
        self.model_select       = QComboBox(self)
        self.model_select.addItem('Select Model')
        self.model_select.addItem('Patagonia Surface and Vertical Model')
        
        self.date_start_label = QLabel('Start Date')
        self.date_start       = QLineEdit('26 Jun 2001')
        self.date_start.setMaximumWidth(100)
        
        self.date_end_label = QLabel('End date')
        self.date_end       = QLineEdit('20 Dec 2001')
        self.date_end.setMaximumWidth(100)
        
        self.lon_east_label = QLabel('Lon Limit East')
        self.lon_east       = QLineEdit('-55')
        self.lon_east.setMaximumWidth(100)
        
        self.lon_west_label = QLabel('Lon Limit West')
        self.lon_west       = QLineEdit('-75')
        self.lon_west.setMaximumWidth(100)

        self.lat_north_label = QLabel('Lat Limit North')
        self.lat_north       = QLineEdit('-50')
        self.lat_north.setMaximumWidth(100)
        
        self.lat_south_label = QLabel('Lat Limit South')
        self.lat_south       = QLineEdit('-60')
        self.lat_south.setMaximumWidth(100)
        
        self.variable_1_label = QLabel('Variable 1')
        self.variable_1       = QComboBox(self)
        
        self.variable_2_label = QLabel('Variable 2')
        self.variable_2       = QComboBox(self)
       
        ### Populate variable comboBoxes with netcdf variables
        varnamelist = self.read_var_data()
        
        self.variable_1.addItem('Choose Variable')
        for name in varnamelist:
            self.variable_1.addItem(name)
        self.variable_1.activated.connect(lambda: self.load_var_1_surf(self.date_start.text(), 
                                          self.date_end.text(), self.lon_east.text(), 
                                          self.lon_west.text(), self.lat_north.text(), 
                                          self.lat_south.text(),self.variable_1.currentText()))
        
        self.variable_2.addItem('Choose Variable')
        for name in varnamelist:
            self.variable_2.addItem(name)
        self.variable_2.activated.connect(lambda: self.load_var_2_surf(self.date_start.text(), 
                                          self.date_end.text(), self.lon_east.text(), 
                                          self.lon_west.text(), self.lat_north.text(), 
                                          self.lat_south.text(),self.variable_2.currentText()))
        
        #### Plot buttons 
        self.surf_surf  = QPushButton('Surface-Surface')
        self.surf_surf.setMaximumWidth(120)
        self.connect(self.surf_surf, SIGNAL('clicked()'), self.draw_surf_surf)
        
        self.depth_depth = QPushButton('Depth-Depth')
        self.depth_depth.setMaximumWidth(120)
        self.connect(self.depth_depth, SIGNAL('clicked()'), self.draw_depth_depth)
        
        self.two_plots = QPushButton('Surface-Depth')
        self.two_plots.setMaximumWidth(120)
        self.connect(self.two_plots, SIGNAL('clicked()'), self.on_draw_2_plots) 
        
        #### User-defined section select button
        self.select_section = QPushButton('Select Section')
        self.select_section.setMaximumWidth(120)
        self.connect(self.select_section, SIGNAL('clicked()'), self.on_select_section)
        
        #### Reset button
        self.reset_button = QPushButton('Reset')
        self.reset_button.setMaximumWidth(100)
        self.connect(self.reset_button, SIGNAL('clicked()'), self.on_reset_button)
        
        #### slider
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMaximumWidth(1287)
        self.slider.sliderMoved.connect(self.slider_moved)
        
        #### play button
        self.play_button = QPushButton('Play')
        self.play_button.setCheckable(True)
        self.play_button.setMaximumWidth(75)
        #self.connect(self.play_button, SIGNAL('clicked()'), self.on_play)
        self.play_button.clicked[bool].connect(self.on_play)

        input_dict = OrderedDict()
        input_dict[self.model_select_label] = self.model_select
        input_dict[self.date_start_label  ] = self.date_start
        input_dict[self.date_end_label    ] = self.date_end
        input_dict[self.lon_east_label    ] = self.lon_east
        input_dict[self.lon_west_label    ] = self.lon_west
        input_dict[self.lat_north_label   ] = self.lat_north
        input_dict[self.lat_south_label   ] = self.lat_south
        input_dict[self.variable_1_label  ] = self.variable_1
        input_dict[self.variable_2_label  ] = self.variable_2
        input_dict[self.surf_surf         ] = self.surf_surf
        input_dict[self.depth_depth       ] = self.depth_depth
        input_dict[self.two_plots         ] = self.two_plots
        input_dict[self.select_section    ] = self.select_section
        input_dict[self.reset_button      ] = self.reset_button
        
        ### End GUI controls ###
        
        ### Window Layout
        hbox        = QHBoxLayout()
        vbox        = QVBoxLayout()
        form_layout = QFormLayout()
        
        for i in input_dict:
            form_layout.addRow(i, input_dict[i])
            
        hbox.addWidget(self.canvas_1)
        hbox.addWidget(self.canvas_2)

        hbox.addLayout(form_layout)
        vbox.addLayout(hbox)
        vbox.addWidget(self.slider)
        vbox.addWidget(self.play_button)
        
        self.plot_window.setLayout(vbox)
        self.setCentralWidget(self.plot_window)
        
        
    def load_var_1_surf(self, date_start, date_end, lon_e, lon_w, 
                            lat_n, lat_s, var_1):
        '''Function to load the variable chosen in comboBox 1
        '''
    
        global file_list            
        global tday_ref, date_ref, date_file_1, tday
        global fnum_min, fnum_max
        global varnamelist
        global F1, var_1_combined
        global it, nt, ii1, jj1, int_i, int_j
        global lon_rho, lat_rho, mask_rho, lon_u, lat_u, mask_u, lon_v, lat_v, mask_v
        global lon_psi, lat_psi, mask_psi, h
        global lon_1, lon_2, lat_1, lat_2
        global var_1_loaded_from
    
        print "variable 1:\t" , var_1
        print "start date:\t" , date_start
        print "end date:\t"   , date_end
        print "lon east:\t"   , lon_e
        print "lon west:\t"   , lon_w
        print "lat north:\t"  , lat_n
        print "lat south:\t"  , lat_s
    
        ### ROMS time
        self.t_date_start = self.yeardayAK(date_start, date_ref, tday_ref)
        self.t_date_end   = self.yeardayAK(date_end,   date_ref, tday_ref)
    
        Tan = []
        Tan.append(self.t_date_start * 3600 * 24)
        Tan.append(self.t_date_end   * 3600 * 24)
    
        ### find indices of files to read from
        d0 = self.yeardayAK(date_file_1, date_ref, tday_ref)
        d1 = self.yeardayAK(date_start , date_ref, tday_ref) - 1
        d2 = self.yeardayAK(date_end   , date_ref, tday_ref) + 1
    
        fnum1 = max([math.floor((d1 - d0) / 30), fnum_min])
        fnum2 = min([math.floor((d2 - d0) / 30), fnum_max])
    
        ### Coordinates of view frame
        lon_1 = float(lon_w)
        lon_2 = float(lon_e)
        lat_1 = float(lat_n)
        lat_2 = float(lat_s)
    
        ### Get variable to plot
        valF       = str(var_1)            # get name of variable to use
        var_1_name = str(var_1)
    
        varlist = np.array(varnamelist[valF])
        
        ### get netcdf global attributes
        nc = Dataset(file_list[0], mode='r')

        ### get N, the number of levels
        for dim in nc.dimensions:
            if dim == 's_rho':
                N = str(nc.dimensions[dim])
                [N for N in N.split() if N.isdigit()]  # get '40' from string
                N = int(N)
            
        theta_s = getattr(nc, 'theta_s')
        theta_b = getattr(nc, 'theta_b')
        Tcline  = getattr(nc, 'Tcline' )
        hc      = getattr(nc, 'hc'     )
    
        ### set attributes for variables
        if valF == 'comb_phyto' or valF == 'comb_zoo':
            for n in varnamelist:
                if n == 'LPHY':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_1_units    = varcheck.units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
        elif valF == 'comb_N':
            for n in varnamelist:
                if n == 'NO3':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_1_units    = varcheck.units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
        elif valF == 'comb_Fe':
            for n in varnamelist:
                if n == 'FeD':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_1_units    = varcheck.units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
        else:            # non-combined variable
            dimid_2  = varnamelist[str(valF)].dimensions[2]
            dimid_3  = varnamelist[str(valF)].dimensions[3]
    
            ### combine variable name and units for label in plots    
            var_1_units    = varnamelist[str(valF)].units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
    
        nc.close()
    
        ### Use eta_rho, xi_u, eta_v appropriately for coordinates
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_u':
            lon_r = lon_u
            lat_r = lat_u
        if dimid_2 == 'eta_v'   and dimid_3 == 'xi_rho':
            lon_r = lon_v
            lat_r = lat_v
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_rho':
            lon_r = lon_rho
            lat_r = lat_rho

        ii1 = []
        jj1 = []

        for idx, val in enumerate(lon_r[0, :]):
            if val >= lon_1 and val <= lon_2:
                ii1.append(idx)   
    
        for idx, val in enumerate(lat_r[:, 0]):
            if val <= lat_1 and val >= lat_2:
                jj1.append(idx)
    
        i1 = ii1[0]
        j1 = jj1[0]
        nx = len(ii1)
        ny = len(jj1)    

        start_3D = (i1, j1, N)
        count_3D = (nx, ny, 1)
    
        print 'start_3D', start_3D
        print 'count_3D', count_3D
    
        F1 = 0          # variable with plotting data contained

        print var_1, varnamelist[valF]
    
        ### squeeze variables 
        #### start with combined variables

        if var_1 == 'comb_phyto':
            vartemp = ('SPHY', 'LPHY')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        elif var_1 == 'comb_zoo':
            vartemp = ('SZOO', 'LZOO', 'PZOO')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        elif var_1 == 'comb_N':
            vartemp = ('NO3', 'NH4', 'PON', 'DON')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        elif var_1 == 'comb_Fe':
            vartemp = ('FeD', 'FeSp', 'FeLp')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        else:           # non-combined variable
            vartemp       = (valF,)
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            F1   = np.squeeze(a[valF])
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        int_i = int(str(re.search('\d+',str(np.shape(ii1))).group()))    # get size of ii1 as int
        int_j = int(str(re.search('\d+',str(np.shape(jj1))).group()))    # get size of jj1 as int
        
        ii1 = np.reshape(ii1, (int_i,1))
        jj1 = np.transpose(np.reshape(jj1, (int_j,1)))
    
        F1 = np.transpose(F1)
        
        self.slider.setRange(0, nt - 1)
        
        var_1_loaded_from = 'load_var_1_surf'
    
        # show the first frame
        it = 0
    

    def load_var_2_surf(self, date_start, date_end, lon_e, lon_w, 
                            lat_n, lat_s, var_2):
        '''Function to load the variable chosen in comboBox 2
        '''
                           
        global file_list            
        global tday_ref, date_ref, date_file_1, tday
        global fnum_min, fnum_max
        global varnamelist, var_2_combined
        global F2, zr, nt, var_2_name
        global it, nt, ii1, jj1, int_i, int_j
        global lon_rho, lat_rho, mask_rho, lon_u, lat_u, mask_u, lon_v, lat_v, mask_v
        global lon_psi, lat_psi, mask_psi, h
        global lon_1, lon_2, lat_1, lat_2
        global var_2_loaded_from
    
        print "\nvariable 2:\t" , var_2
        print "start date:\t"   , date_start
        print "end date:\t"     , date_end
        print "lon east:\t"     , lon_e
        print "lon west:\t"     , lon_w
        print "lat north:\t"    , lat_n
        print "lat south:\t"    , lat_s
    
        ### ROMS time
        self.t_date_start = self.yeardayAK(date_start, date_ref, tday_ref)
        self.t_date_end   = self.yeardayAK(date_end,   date_ref, tday_ref)
    
        Tan = []
        Tan.append(self.t_date_start * 3600 * 24)
        Tan.append(self.t_date_end   * 3600 * 24)
    
        ### find indices of files to read from
        d0 = self.yeardayAK(date_file_1, date_ref, tday_ref)
        d1 = self.yeardayAK(date_start , date_ref, tday_ref) - 1
        d2 = self.yeardayAK(date_end   , date_ref, tday_ref) + 1
    
        fnum1 = max([math.floor((d1 - d0) / 30), fnum_min])
        fnum2 = min([math.floor((d2 - d0) / 30), fnum_max])
    
        ### Coordinates of view frame
        lon_1 = float(lon_w)
        lon_2 = float(lon_e)
        lat_1 = float(lat_n)
        lat_2 = float(lat_s)
    
        ### Get variable to plot
        valF       = str(var_2)            # get name of variable to use
        var_2_name = str(var_2)
    
        varlist = np.array(varnamelist[valF])
        
        ### get netcdf global attributes
        nc = Dataset(file_list[0], mode='r')

        ### get N, the number of levels
        for dim in nc.dimensions:
            if dim == 's_rho':
                N = str(nc.dimensions[dim])
                [N for N in N.split() if N.isdigit()]  # get '40' from string
                N = float(N)
            
        theta_s = getattr(nc, 'theta_s')
        theta_b = getattr(nc, 'theta_b')
        Tcline  = getattr(nc, 'Tcline' )
        hc      = getattr(nc, 'hc'     )
    
        ### set attributes for variables
        if valF == 'comb_phyto' or valF == 'comb_zoo':
            for n in varnamelist:
                if n == 'LPHY':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varcheck.units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
        elif valF == 'comb_N':
            for n in varnamelist:
                if n == 'NO3':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varcheck.units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
        elif valF == 'comb_Fe':
            for n in varnamelist:
                if n == 'FeD':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varcheck.units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
        else:            # non-combined variable
            dimid_2  = varnamelist[str(valF)].dimensions[2]
            dimid_3  = varnamelist[str(valF)].dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varnamelist[str(valF)].units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
    
        nc.close()
    
        ### Use eta_rho, xi_u, eta_v appropriately for coordinates
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_u':
            lon_r = lon_u
            lat_r = lat_u
        if dimid_2 == 'eta_v'   and dimid_3 == 'xi_rho':
            lon_r = lon_v
            lat_r = lat_v
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_rho':
            lon_r = lon_rho
            lat_r = lat_rho
    
        ii1 = []
        jj1 = []

        for idx, val in enumerate(lon_r[0, :]):
            if val >= lon_1 and val <= lon_2:
                ii1.append(idx)   
    
        for idx, val in enumerate(lat_r[:, 0]):
            if val <= lat_1 and val >= lat_2:
                jj1.append(idx)
    
        i1 = ii1[0]
        j1 = jj1[0]
        nx = len(ii1)
        ny = len(jj1)     

        start_3D = (i1, j1, int(N))
        count_3D = (nx, ny, 1)
    
        print 'start_3D', start_3D
        print 'count_3D', count_3D   
    
        F2 = 0          # variable to contain plotting data

        print var_2, varnamelist[valF]
    
        ### squeeze variables 
    
        #### start with combined variables
        if var_2 == 'comb_phyto':
            vartemp = ('SPHY', 'LPHY')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        elif var_2 == 'comb_zoo':
            vartemp = ('SZOO', 'LZOO', 'PZOO')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        elif var_2 == 'comb_N':
            vartemp = ('NO3', 'NH4', 'PON', 'DON')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        elif var_2 == 'comb_Fe':
            vartemp = ('FeD', 'FeSp', 'FeLp')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        else:           # non-combined variable
            vartemp       = (valF,)
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            F2   = np.squeeze(a[valF])
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        int_i = int(str(re.search('\d+',str(np.shape(ii1))).group()))    # get size of ii1 as int
        int_j = int(str(re.search('\d+',str(np.shape(jj1))).group()))    # get size of jj1 as int
        
        ii1 = np.reshape(ii1, (int_i,1))
        jj1 = np.transpose(np.reshape(jj1, (int_j,1)))
  
        F2 = np.transpose(F2)
        
        var_2_loaded_from = 'load_var_2_surf'
        
        self.slider.setRange(0, nt - 1)
    
        ### show the first frame
        it = 0
        

    def load_var_1_depth(self, date_start, date_end, lon_e, lon_w, 
                            lat_n, lat_s, var_1):
        '''Function to load the variable chosen in comboBox 1
        '''
        
        global file_list            
        global tday_ref, date_ref, date_file_1, tday
        global fnum_min, fnum_max
        global varnamelist
        global F1, var_1_combined
        global it, nt, ii1, jj1, int_i, int_j
        global lon_rho, lat_rho, mask_rho, lon_u, lat_u, mask_u, lon_v, lat_v, mask_v
        global lon_psi, lat_psi, mask_psi, h
        global lon_1, lon_2, lat_1, lat_2
        global var_1_loaded_from
        
        print "variable 1:\t" , var_1
        print "start date:\t" , date_start
        print "end date:\t"   , date_end
        print "lon east:\t"   , lon_e
        print "lon west:\t"   , lon_w
        print "lat north:\t"  , lat_n
        print "lat south:\t"  , lat_s
        
        ### ROMS time
        self.t_date_start = self.yeardayAK(date_start, date_ref, tday_ref)
        self.t_date_end   = self.yeardayAK(date_end,   date_ref, tday_ref)
        
        Tan = []
        Tan.append(self.t_date_start * 3600 * 24)
        Tan.append(self.t_date_end   * 3600 * 24)
        
        ### find indices of files to read from
        d0 = self.yeardayAK(date_file_1, date_ref, tday_ref)
        d1 = self.yeardayAK(date_start , date_ref, tday_ref) - 1
        d2 = self.yeardayAK(date_end   , date_ref, tday_ref) + 1
        
        fnum1 = max([math.floor((d1 - d0) / 30), fnum_min])
        fnum2 = min([math.floor((d2 - d0) / 30), fnum_max])
        
        ### Coordinates of view frame
        lon_1 = float(lon_w)
        lon_2 = float(lon_e)
        lat_1 = float(lat_n)
        lat_2 = float(lat_s)
        
        ### Get variable to plot
        valF       = str(var_1)            # get name of variable to use
        var_1_name = str(var_1)
        
        varlist = np.array(varnamelist[valF])
            
        ### get netcdf global attributes
        nc = Dataset(file_list[0], mode='r')
 
        ### get N, the number of levels
        for dim in nc.dimensions:
            if dim == 's_rho':
                N = str(nc.dimensions[dim])
                [N for N in N.split() if N.isdigit()]  # get '40' from string
                N = int(N)
                
        theta_s = getattr(nc, 'theta_s')
        theta_b = getattr(nc, 'theta_b')
        Tcline  = getattr(nc, 'Tcline' )
        hc      = getattr(nc, 'hc'     )
        
        ### set attributes for variables
        if valF == 'comb_phyto' or valF == 'comb_zoo':
            for n in varnamelist:
                if n == 'LPHY':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_1_units    = varcheck.units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
        elif valF == 'comb_N':
            for n in varnamelist:
                if n == 'NO3':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_1_units    = varcheck.units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
        elif valF == 'comb_Fe':
            for n in varnamelist:
                if n == 'FeD':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_1_units    = varcheck.units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
        else:            # non-combined variable
            dimid_2  = varnamelist[str(valF)].dimensions[2]
            dimid_3  = varnamelist[str(valF)].dimensions[3]
        
            ### combine variable name and units for label in plots    
            var_1_units    = varnamelist[str(valF)].units
            var_1_combined = var_1_name + ' ' + '(' + var_1_units + ')'
        
        nc.close()
        
        ### Use eta_rho, xi_u, eta_v appropriately for coordinates
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_u':
            lon_r = lon_u
            lat_r = lat_u
        if dimid_2 == 'eta_v'   and dimid_3 == 'xi_rho':
            lon_r = lon_v
            lat_r = lat_v
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_rho':
            lon_r = lon_rho
            lat_r = lat_rho
        
        ii1 = []
        jj1 = []

        for idx, val in enumerate(lon_r[0, :]):
            if val >= lon_1 and val <= lon_2:
                ii1.append(idx)   
        
        for idx, val in enumerate(lat_r[:, 0]):
            if val <= lat_1 and val >= lat_2:
                jj1.append(idx)
        
        i1 = ii1[0]
        j1 = jj1[0]
        nx = len(ii1)
        ny = len(jj1)    

        # start_3D = (i1, j1, N)
        # count_3D = (nx, ny, 1)
        
        start_3D = (i1, j1, 1)
        count_3D = (nx, ny, N)
        
        print 'start_3D', start_3D
        print 'count_3D', count_3D
        
        int_i = int(str(re.search('\d+',str(np.shape(ii1))).group()))    # get size of ii1 as int
        int_j = int(str(re.search('\d+',str(np.shape(jj1))).group()))    # get size of jj1 as int

        ii1 = np.reshape(ii1, (int_i,1))
        jj1 = np.transpose(np.reshape(jj1, (int_j,1)))

        ### set depth variable zr
        vartemp = ('zeta',)
        scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp,
                                        start_3D, count_3D, Tan)
        zeta1 = a['zeta']
        zeta1 = np.transpose(zeta1)

        Vst = 1
        Vtr = 1

        for i in range(int(nt)):
            ssh = np.squeeze(zeta1[:,:,i])
            zr = self.set_depth_(Vtr, Vst, theta_s, theta_b, hc, N, 1, h[jj1,ii1], ssh, 0)
        
        F1 = 0          # variable with plotting data contained

        print var_1, varnamelist[valF]
        
        ### squeeze variables 
        #### start with combined variables

        if var_1 == 'comb_phyto':
            vartemp = ('SPHY', 'LPHY')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        elif var_1 == 'comb_zoo':
            vartemp = ('SZOO', 'LZOO', 'PZOO')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        elif var_1 == 'comb_N':
            vartemp = ('NO3', 'NH4', 'PON', 'DON')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        elif var_1 == 'comb_Fe':
            vartemp = ('FeD', 'FeSp', 'FeLp')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F1 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        else:           # non-combined variable
            vartemp       = (valF,)
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            F1   = np.squeeze(a[valF])
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
        
        F1 = np.transpose(F1)
        
        var_1_loaded_from = 'load_var_1_depth'
        
        self.slider.setRange(0, nt - 1)
        
        # show the first frame
        it = 0
        

    def load_var_2_depth(self, date_start, date_end, lon_e, lon_w, 
                            lat_n, lat_s, var_2):
        '''Function to load the variable chosen in comboBox 2
        '''
                               
        global file_list            
        global tday_ref, date_ref, date_file_1, tday
        global fnum_min, fnum_max
        global varnamelist, var_2_combined
        global F2, zr, nt, var_2_name
        global it, nt, ii1, jj1, int_i, int_j
        global lon_rho, lat_rho, mask_rho, lon_u, lat_u, mask_u, lon_v, lat_v, mask_v
        global lon_psi, lat_psi, mask_psi, h
        global lon_1, lon_2, lat_1, lat_2
        global var_2_loaded_from
        
        print "\nvariable 2:\t" , var_2
        print "start date:\t"   , date_start
        print "end date:\t"     , date_end
        print "lon east:\t"     , lon_e
        print "lon west:\t"     , lon_w
        print "lat north:\t"    , lat_n
        print "lat south:\t"    , lat_s
        
        ### ROMS time
        self.t_date_start = self.yeardayAK(date_start, date_ref, tday_ref)
        self.t_date_end   = self.yeardayAK(date_end,   date_ref, tday_ref)
        
        Tan = []
        Tan.append(self.t_date_start * 3600 * 24)
        Tan.append(self.t_date_end   * 3600 * 24)
        
        ### find indices of files to read from
        d0 = self.yeardayAK(date_file_1, date_ref, tday_ref)
        d1 = self.yeardayAK(date_start , date_ref, tday_ref) - 1
        d2 = self.yeardayAK(date_end   , date_ref, tday_ref) + 1
        
        fnum1 = max([math.floor((d1 - d0) / 30), fnum_min])
        fnum2 = min([math.floor((d2 - d0) / 30), fnum_max])
        
        ### Coordinates of view frame
        lon_1 = float(lon_w)
        lon_2 = float(lon_e)
        lat_1 = float(lat_n)
        lat_2 = float(lat_s)
        
        ### Get variable to plot
        valF       = str(var_2)            # get name of variable to use
        var_2_name = str(var_2)
        
        varlist = np.array(varnamelist[valF])
            
        ### get netcdf global attributes
        nc = Dataset(file_list[0], mode='r')
 
        ### get N, the number of levels
        for dim in nc.dimensions:
            if dim == 's_rho':
                N = str(nc.dimensions[dim])
                [N for N in N.split() if N.isdigit()]  # get '40' from string
                N = float(N)
                
        theta_s = getattr(nc, 'theta_s')
        theta_b = getattr(nc, 'theta_b')
        Tcline  = getattr(nc, 'Tcline' )
        hc      = getattr(nc, 'hc'     )
        
        ### set attributes for variables
        if valF == 'comb_phyto' or valF == 'comb_zoo':
            for n in varnamelist:
                if n == 'LPHY':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varcheck.units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
        elif valF == 'comb_N':
            for n in varnamelist:
                if n == 'NO3':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varcheck.units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
        elif valF == 'comb_Fe':
            for n in varnamelist:
                if n == 'FeD':
                    varcheck = varnamelist[n]
            dimid_2 = varcheck.dimensions[2]
            dimid_3 = varcheck.dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varcheck.units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
        else:            # non-combined variable
            dimid_2  = varnamelist[str(valF)].dimensions[2]
            dimid_3  = varnamelist[str(valF)].dimensions[3]
            ### combine variable name and units for label in plots    
            var_2_units    = varnamelist[str(valF)].units
            var_2_combined = var_2_name + ' ' + '(' + var_2_units + ')'
        
        nc.close()
        
        ### Use eta_rho, xi_u, eta_v appropriately for coordinates
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_u':
            lon_r = lon_u
            lat_r = lat_u
        if dimid_2 == 'eta_v'   and dimid_3 == 'xi_rho':
            lon_r = lon_v
            lat_r = lat_v
        if dimid_2 == 'eta_rho' and dimid_3 == 'xi_rho':
            lon_r = lon_rho
            lat_r = lat_rho
        
        ii1 = []
        jj1 = []

        for idx, val in enumerate(lon_r[0, :]):
            if val >= lon_1 and val <= lon_2:
                ii1.append(idx)   
        
        for idx, val in enumerate(lat_r[:, 0]):
            if val <= lat_1 and val >= lat_2:
                jj1.append(idx)
        
        i1 = ii1[0]
        j1 = jj1[0]
        nx = len(ii1)
        ny = len(jj1)     

        start_3D = (i1, j1, 1)
        count_3D = (nx, ny, int(N))
        
        print 'start_3D', start_3D
        print 'count_3D', count_3D
        
        int_i = int(str(re.search('\d+',str(np.shape(ii1))).group()))    # get size of ii1 as int
        int_j = int(str(re.search('\d+',str(np.shape(jj1))).group()))    # get size of jj1 as int
        
        ii1 = np.reshape(ii1, (int_i,1))
        jj1 = np.transpose(np.reshape(jj1, (int_j,1)))
        
        ### set depth variable zr
        vartemp = ('zeta',)
        scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp,
                                        start_3D, count_3D, Tan)                               
        zeta1 = a['zeta']
        zeta1 = np.transpose(zeta1)
                
        Vst = 1
        Vtr = 1

        for i in range(int(nt)):
            ssh = np.squeeze(zeta1[:,:,i])
            zr = self.set_depth_(Vtr, Vst, theta_s, theta_b, hc, N, 1, h[jj1,ii1], ssh, 0)    
        
        F2 = 0          # variable to contain plotting data

        print var_2, varnamelist[valF]
        
        ### squeeze variables 
        
        #### start with combined variables
        if var_2 == 'comb_phyto':
            vartemp = ('SPHY', 'LPHY')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        elif var_2 == 'comb_zoo':
            vartemp = ('SZOO', 'LZOO', 'PZOO')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        elif var_2 == 'comb_N':
            vartemp = ('NO3', 'NH4', 'PON', 'DON')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        elif var_2 == 'comb_Fe':
            vartemp = ('FeD', 'FeSp', 'FeLp')
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            for var in a:
                Ftemp = np.squeeze(a[var])
                F2 += Ftemp 
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
            
        else:           # non-combined variable
            vartemp       = (valF,)
            scrum_time, a = self.read_his_fastAK(fnum1, fnum2, vartemp, 
                                                 start_3D, count_3D, Tan)
            F2   = np.squeeze(a[valF])
            nt   = len(scrum_time)
            tday = scrum_time / 24 / 3600
 
        F2 = np.transpose(F2)

        var_2_loaded_from = 'load_var_2_depth'
        
        self.slider.setRange(0, nt - 1)
        
        ### show the first frame
        it = 0

         
    def set_depth_(self, Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, zeta, report):
        '''function that implements set_depth.m in python  
        '''
        
        varargin = [Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, zeta, report]
        nargin = len(varargin)
        
        s1, s2 = (np.shape(zeta))
        z = np.zeros((s1, s2, N))       # initialize z
        
        if (nargin < 8):
            disp_(char(' '))
            disp_([setstr_(7),char('*** Error:  SET_DEPTH - too few arguments.'),setstr_(7)])
            disp_([setstr_(7),char('                     number of supplied arguments: '),num2str_(nargin),setstr_(7)])
            disp_([setstr_(7),char('                     number of required arguments: 8'),setstr_(7)])
            disp_(char(' '))
            return z
        if (Vtransform < 1 or Vtransform > 2):
            disp_(char(' '))
            disp_([setstr_(7),char('*** Error:  SET_DEPTH - Illegal parameter Vtransform = '),num2str_(Vtransform),setstr_(7)])
            disp_(char(' '))
            return z
        if (Vstretching < 1 or Vstretching > 4):
            disp_(char(' '))
            disp_([setstr_(7),char('*** Error:  SET_DEPTH - Illegal parameter Vstretching = '),num2str_(Vstretching),setstr_(7)])
            disp_(char(' '))
            return z
        if (hc > h.min() and Vtransform == 1):
            disp_(char(' '))
            disp_([setstr_(7),char('*** Error:  SET_DEPTH - critical depth exceeds minimum'),char(' bathymetry value.'),setstr_(7)])
            disp_([setstr_(7),char('                        Vtranform = '),num2str_(Vtransform),setstr_(7)])
            disp_([setstr_(7),char('                        hc        = '),num2str_(hc),setstr_(7)])
            disp_([setstr_(7),char('                        hmax      = '),num2str_(min_(min_(h))),setstr_(7)])
            disp_(char(' '))
            return z
        if (nargin < 9):
            zeta=zeros_(size_(h))
        if (nargin < 10):
            report=1
            
        Np = N + 1
        Lp, Mp = (np.shape(h))
        L = Lp - 1
        M = Mp - 1
        
        hmin=h.min()
        hmax=h.max()
        
        if (report):
            disp_(char(' '))
            if (Vtransform == 1):
                disp_([char('Vtransform  = '),num2str_(Vtransform),char('   original ROMS')])
            else:
                if (Vtransform == 2):
                    disp_([char('Vtransform  = '),num2str_(Vtransform),char('   ROMS-UCLA')])
            if 1 == (igrid):
                disp_([char('   igrid    = '),num2str_(igrid),char('   at horizontal RHO-points')])
            else:
                if 2 == (igrid):
                    disp_([char('   igrid    = '),num2str_(igrid),char('   at horizontal PSI-points')])
                else:
                    if 3 == (igrid):
                        disp_([char('   igrid    = '),num2str_(igrid),char('   at horizontal U-points')])
                    else:
                        if 4 == (igrid):
                            disp_([char('   igrid    = '),num2str_(igrid),char('   at horizontal V-points')])
                        else:
                            if 5 == (igrid):
                                disp_([char('   igrid    = '),num2str_(igrid),char('   at horizontal RHO-points')])
        if (igrid == 5):
            kgrid = 1
        else:
            kgrid = 0
            
        s, C = self.stretching_(Vstretching, theta_s, theta_b, hc, N, kgrid, report)
        
        if 1 == (igrid):
            hr    = h
            zetar = zeta
        else:
            if 2 == (igrid):
                hp=0.25.dot((h[1:L,1:M] + h[2:Lp,1:M] + h[1:L,2:Mp] + h[2:Lp,2:Mp]))
                zetap=0.25.dot((zeta[1:L,1:M] + zeta[2:Lp,1:M] + zeta[1:L,2:Mp] + zeta[2:Lp,2:Mp]))
            else:
                if 3 == (igrid):
                    hu=0.5.dot((h[1:L,1:Mp] + h[2:Lp,1:Mp]))
                    zetau=0.5.dot((zeta[1:L,1:Mp] + zeta[2:Lp,1:Mp]))
                else:
                    if 4 == (igrid):
                        hv=0.5.dot((h[1:Lp,1:M] + h[1:Lp,2:Mp]))
                        zetav=0.5.dot((zeta[1:Lp,1:M] + zeta[1:Lp,2:Mp]))
                    else:
                        if 5 == (igrid):
                            hr=copy_(h)
                            zetar=copy_(zeta)
        if (Vtransform == 1):
            if 1 == (igrid):
                for k in range(0, int(N)-1):
                    z0 = np.float_((s[k] - C[k]) * hc) + np.float_(C[k] * (hr))
                    z[:,:,k] = z0 + zetar * (1.0 + (z0 / hr))
            else:
                if 2 == (igrid):
                    for k in arange_(1,N).reshape(-1):
                        z0=(s[k] - C[k]) * hc + C[k].dot(hp)
                        z[:,:,k]=z0 + zetap.dot((1.0 + z0 / hp))
                else:
                    if 3 == (igrid):
                        for k in arange_(1,N).reshape(-1):
                            z0=(s[k] - C[k]) * hc + C[k].dot(hu)
                            z[:,:,k]=z0 + zetau.dot((1.0 + z0 / hu))
                    else:
                        if 4 == (igrid):
                            for k in arange_(1,N).reshape(-1):
                                z0=(s[k] - C[k]) * hc + C[k].dot(hv)
                                z[:,:,k]=z0 + zetav.dot((1.0 + z0 / hv))
                        else:
                            if 5 == (igrid):
                                z[:,:,1]=- hr
                                for k in arange_(2,Np).reshape(-1):
                                    z0=(s[k] - C[k]) * hc + C[k].dot(hr)
                                    z[:,:,k]=z0 + zetar.dot((1.0 + z0 / hr))
        else:
            if (Vtransform == 2):
                if 1 == (igrid):
                    for k in arange_(1,N).reshape(-1):
                        z0=(hc.dot(s[k]) + C[k].dot(hr)) / (hc + hr)
                        z[:,:,k]=zetar + (zeta + hr).dot(z0)
                else:
                    if 2 == (igrid):
                        for k in arange_(1,N).reshape(-1):
                            z0=(hc.dot(s[k]) + C[k].dot(hp)) / (hc + hp)
                            z[:,:,k]=zetap + (zetap + hp).dot(z0)
                    else:
                        if 3 == (igrid):
                            for k in arange_(1,N).reshape(-1):
                                z0=(hc.dot(s[k]) + C[k].dot(hu)) / (hc + hu)
                                z[:,:,k]=zetau + (zetau + hu).dot(z0)
                        else:
                            if 4 == (igrid):
                                for k in arange_(1,N).reshape(-1):
                                    z0=(hc.dot(s[k]) + C[k].dot(hv)) / (hc + hv)
                                    z[:,:,k]=zetav + (zetav + hv).dot(z0)
                            else:
                                if 5 == (igrid):
                                    for k in arange_(1,Np).reshape(-1):
                                        z0=(hc.dot(s[k]) + C[k].dot(hr)) / (hc + hr)
                                        z[:,:,k]=zetar + (zetar + hr).dot(z0)
        return z
        

    def stretching_(self, Vstretching, theta_s, theta_b, hc, N, kgrid, report):
        '''function that implements stretching.m in python
        '''
        
        varargin = [Vstretching, theta_s, theta_b, hc, N, kgrid, report]
        nargin   = len(varargin)

        if (nargin < 6):
            disp_(char(' '))
            disp_(char('*** Error:  STRETCHING - too few arguments.'))
            disp_([char('                     number of supplied arguments: '),num2str_(nargin)])
            disp_(char('                     number of required arguments: 6'))
            disp_(char(' '))
            return s,C
        if (Vstretching < 1 or Vstretching > 4):
            disp_(char(' '))
            disp_([char('*** Error:  STRETCHING - Illegal parameter Vstretching = '),num2str_(Vstretching)])
            disp_(char(' '))
            return s,C
        if (nargin < 7):
            report=copy(false)
        Np = N + 1
        if (Vstretching == 1):
            ds=1.0 / N
            if (kgrid == 1):
                Nlev=copy_(Np)
                lev=arange_(0,N)
                s=(lev - N).dot(ds)
            else:
                Nlev = N
                lev  = (np.arange(1, N)) - 0.5
                s    = (lev - N) * (ds)
            if (theta_s > 0):
                Ptheta = np.sinh(theta_s * (s)) / np.sinh(theta_s)
                Rtheta = np.tanh(theta_s * ((s + 0.5))) / (2.0 * np.tanh(0.5 * theta_s)) - 0.5
                C=(1.0 - theta_b) * (Ptheta) + theta_b * (Rtheta)
            else:
                C=copy_(s)
        else:
            if (Vstretching == 2):
                alfa=1.0
                beta=1.0
                ds=1.0 / N
                if (kgrid == 1):
                    Nlev=copy_(Np)
                    lev=arange_(0,N)
                    s=(lev - N).dot(ds)
                else:
                    Nlev=copy_(N)
                    lev=(arange_(1,N)) - 0.5
                    s=(lev - N).dot(ds)
                if (theta_s > 0):
                    Csur=(1.0 - cosh_(theta_s.dot(s))) / (cosh_(theta_s) - 1.0)
                    if (theta_b > 0):
                        Cbot=- 1.0 + sinh_(theta_b * (s + 1.0)) / sinh_(theta_b)
                        weigth=(s + 1.0) ** alfa.dot((1.0 + (alfa / beta).dot((1.0 - (s + 1.0) ** beta))))
                        C=weigth.dot(Csur) + (1.0 - weigth).dot(Cbot)
                    else:
                        C=copy_(Csur)
                else:
                    C=copy_(s)
            else:
                if (Vstretching == 3):
                    ds=1.0 / N
                    if (kgrid == 1):
                        Nlev=copy_(Np)
                        lev=arange_(0,N)
                        s=(lev - N).dot(ds)
                    else:
                        Nlev=copy_(N)
                        lev=(arange_(1,N)) - 0.5
                        s=(lev - N).dot(ds)
                    if (theta_s > 0):
                        exp_s=copy_(theta_s)
                        exp_b=copy_(theta_b)
                        alpha=3
                        Cbot=log_(cosh_(alpha * (s + 1) ** exp_b)) / log_(cosh_(alpha)) - 1
                        Csur=- log_(cosh_(alpha * abs_(s) ** exp_s)) / log_(cosh_(alpha))
                        weight=(1 - tanh_(alpha * (s + 0.5))) / 2
                        C=weight.dot(Cbot) + (1 - weight).dot(Csur)
                    else:
                        C=copy_(s)
                else:
                    if (Vstretching == 4):
                        ds=1.0 / N
                        if (kgrid == 1):
                            Nlev=copy_(Np)
                            lev=arange_(0,N)
                            s=(lev - N).dot(ds)
                        else:
                            Nlev=copy_(N)
                            lev=(arange_(1,N)) - 0.5
                            s=(lev - N).dot(ds)
                        if (theta_s > 0):
                            Csur=(1.0 - cosh_(theta_s.dot(s))) / (cosh_(theta_s) - 1.0)
                        else:
                            Csur=- s ** 2
                        if (theta_b > 0):
                            Cbot=(exp_(theta_b.dot(Csur)) - 1.0) / (1.0 - exp_(- theta_b))
                            C=copy_(Cbot)
                        else:
                            C=copy_(Csur)
        if (report):
            disp_(char(' '))
            if (Vstretching == 1):
                disp_([char('Vstretching = '),num2str_(Vstretching),char('   Song and Haidvogel (1994)')])
            else:
                if (Vstretching == 2):
                    disp_([char('Vstretching = '),num2str_(Vstretching),char('   Shchepetkin (2005)')])
                else:
                    if (Vstretching == 3):
                        disp_([char('Vstretching = '),num2str_(Vstretching),char('   Geyer (2009), BBL')])
                    else:
                        if (Vstretching == 4):
                            disp_([char('Vstretching = '),num2str_(Vstretching),char('   Shchepetkin (2010)')])
            if (kgrid == 1):
                disp_([char('   kgrid    = '),num2str_(kgrid),char('   at vertical W-points')])
            else:
                disp_([char('   kgrid    = '),num2str_(kgrid),char('   at vertical RHO-points')])
            disp_([char('   theta_s  = '),num2str_(theta_s)])
            disp_([char('   theta_b  = '),num2str_(theta_b)])
            disp_([char('   hc       = '),num2str_(hc)])
            disp_(char(' '))
            disp_(char(' S-coordinate curves: k, s(k), C(k)'))
            disp_(char(' '))
            if (kgrid == 1):
                for k in arange_(Nlev,1,- 1).reshape(-1):
                    disp_([char('    '),sprintf_(char('%3g'),k - 1),char('   '),sprintf_(char('%20.12e'),s[k]),char('   '),sprintf_(char('%20.12e'),C[k])])
            else:
                for k in arange_(Nlev,1,- 1).reshape(-1):
                    disp_([char('    '),sprintf_(char('%3g'),k),char('   '),sprintf_(char('%20.12e'),s[k]),char('   '),sprintf_(char('%20.12e'),C[k])])
            disp_(char(' '))
        return (s, C)


    def dateAK(self, day, date_ref, tday_ref):
        '''dateAK function rewritten in Python
        '''
        
        date_ref = datetime.strptime(date_ref, "%d %b %Y")
        date_ref = datetime.date(date_ref)
    
        tdays_jc = day - tday_ref + datetime.toordinal(date_ref)

        ymd = datetime.fromordinal(int(tdays_jc - 365))     # -365 to account for python starting at year 1 instead of 0
    
        date = datetime.date(ymd)
        date = date.strftime("%d %b %Y")
    
        return date
        
 
    def yeardayAK(self, date, date_ref, day_ref):
        '''yeardayAK function rewritten in python
           MATLAB function: day=datenum(date)-datenum(date_ref)+day_ref;
        '''
        
        date = datetime.strptime(str(date), "%d %b %Y")
        date = datetime.date(date) 
        
        date_ref = datetime.strptime(date_ref, "%d %b %Y")
        date_ref = datetime.date(date_ref)
        
        day  = datetime.toordinal(date) - datetime.toordinal(date_ref) + day_ref

        return day
        
   
    def read_his_fastAK(self, fnum1, fnum2, field_list, start_3D, count_3D, Tan):
        '''read_his_fastAK function rewritten in Python
        '''
        
        global exp_dir, file_list, fhead, nt
        
        scrum_time = []
        a          = {}

        fnums = range(int(fnum1), int(fnum2) + 1)   # +1 to include fnum2
        nf    = len(fnums)
        nnt   = np.zeros(nf)
        
        ### Open time files and extract scrum_time
        ### .nc files must be in the format <fhead>.000x.nc
        for kf in range(nf):
            fname   = exp_dir + fhead + '000' + str(fnums[kf]) + '.nc'
            nc      = Dataset(fname, mode='r')
            t       = nc.variables['scrum_time'][:]
            IN      = [idx for idx, val in enumerate(t)]
            nnt[kf] = len(IN)
            nc.close()
        
        nt = sum(nnt)
        
        ### initialize output arrays
        scrum_time = np.zeros((nt, 1))

        if nt > 0:
            fname   = exp_dir + fhead + '000' + str(fnums[kf]) + '.nc'
            nc      = Dataset(fname, mode='r')
            for kv in field_list:
                varid     = nc.variables[kv]
                ndims     = len(varid.dimensions)
                if ndims == 3:      # 2D field
                    a[kv] = np.zeros(((nt, count_3D[1], count_3D[0])))
                elif ndims == 4:    # 3D field
                    a[kv] = np.zeros(((nt, count_3D[2], count_3D[1], count_3D[0])))
        nc.close()      
        
        ### fill in arrays (dictionary a)
        it1 = 0
    
        for kf in range(nf):
            if nnt[kf] > 0:
                fname = exp_dir + fhead + '000' + str(fnums[kf]) + '.nc'
                nc              = Dataset(fname, mode='r')
                t               = nc.variables['scrum_time'][:]
                IN              = [idx for idx, val in enumerate(t)]
                ntk             = len(IN)
                itt             = range(it1, (it1 + ntk))
                t               = np.reshape(t, (ntk,1))
                scrum_time[itt] = t[IN]
                
                for kv in field_list:
                    varid = nc.variables[kv]
                    ndims = len(varid.dimensions)

                    #### grab data from varid according to the lat/lon chosen
                    if ndims == 3:          # 2D Field
                        a[kv][itt,:,:] = varid[IN[0]:ntk,
                                               start_3D[1]:start_3D[1]+count_3D[1],
                                               start_3D[0]:start_3D[0]+count_3D[0]]
                    elif ndims == 4:        # 3D Field
                        a[kv][itt,:,:,:] = varid[IN[0]:ntk,
                                                 start_3D[2]-1:start_3D[2]+count_3D[2],
                                                 start_3D[1]  :start_3D[1]+count_3D[1],
                                                 start_3D[0]  :start_3D[0]+count_3D[0]]
                it1 += ntk
                nc.close()

        return (scrum_time, a)
        
        
    def read_grid_data(self):
        '''Function to read the grid data from the netcdf file
        '''
        
        global fhead
        global lon_rho, lat_rho, mask_rho, lon_u, lat_u, mask_u, lon_v, lat_v, mask_v
        global lon_psi, lat_psi, mask_psi, h
        global tday_ref, date_ref, date_file_1
        global xi_rho, eta_rho, xi_psi, eta_psi
        global xi_u, eta_u, xi_v, eta_v
        global lon_1, lon_2, lat_1, lat_2
        global mask_rho_nan
        
        ### Set up grid and time variables 
        grdname = '/Users/zach/Documents/MATLAB/VData/roms_grd.nc.1'
        fhead   = '2001_NEMURO_avg_p.'

        fnum_min = 1
        fnum_max = 6

        tday_ref    = 0
        date_ref    = '01 Jan 0001' # Python ordinal starts at 1 Jan 1
        date_file_1 = '26 Jun 2001' # <- the date of file with index 1
       
        ### read grid into variables and transpose
        nc       = Dataset(grdname, mode='r')
        
        lon_rho  = nc.variables['lon_rho' ][:]
        lat_rho  = nc.variables['lat_rho' ][:]
        mask_rho = nc.variables['mask_rho'][:]
        lon_u    = nc.variables['lon_u'   ][:]
        lat_u    = nc.variables['lat_u'   ][:]
        mask_u   = nc.variables['mask_u'  ][:]
        lon_v    = nc.variables['lon_v'   ][:]
        lat_v    = nc.variables['lat_v'   ][:]
        mask_v   = nc.variables['mask_v'  ][:]
        lon_psi  = nc.variables['lon_psi' ][:]
        lat_psi  = nc.variables['lat_psi' ][:]
        mask_psi = nc.variables['mask_psi'][:]
        h        = nc.variables['h'       ][:]
        
        nc.close()

        ### assign dimensions
        (xi_rho, eta_rho) = np.shape(lon_rho)
        (xi_u, eta_u)     = np.shape(lon_u)
        (xi_v  , eta_v)   = np.shape(lon_v)
        (xi_psi, eta_psi) = np.shape(lon_psi)
        
        ### lat/lon limits of the grid   
        lon_1 = lon_rho[ 1, 1]
        lon_2 = lon_rho[-1, 1]
        lat_1 = lat_rho[ 1, 1]
        lat_2 = lat_rho[-1, 1]
        
        #### convert 0s to NaN
        mask_rho_nan = mask_rho.astype(float)
        mask_rho_nan[mask_rho_nan == 0] = np.nan
        
        mask_psi_nan = mask_psi.astype(float)
        mask_psi_nan[mask_psi_nan == 0] = np.nan
        
        mask_u_nan   = mask_u.astype(float)
        mask_u_nan[mask_u_nan == 0]     = np.nan
        
        mask_v_nan   = mask_v.astype(float)
        mask_v_nan[mask_v_nan == 0]     = np.nan
        
        ### ensure nans are a different color (masked)
        mask_rho_nan = np.ma.array(mask_rho_nan, mask=np.isnan(mask_rho_nan))
        mask_psi_nan = np.ma.array(mask_psi_nan, mask=np.isnan(mask_psi_nan))
        mask_u_nan   = np.ma.array(mask_u_nan  , mask=np.isnan(mask_u_nan  ))
        mask_v_nan   = np.ma.array(mask_v_nan  , mask=np.isnan(mask_v_nan  ))  


    def read_var_data(self):
        '''Function to read the variables from the netcdf files into the variable_1
           and variable_2 combo boxes
        '''
        global exp_dir, file_list, fhead
        global fnum_min, fnum_max
        global varnamelist
        
        ### access netcdf files
        exp_dir = '/Users/zach/Documents/MATLAB/VData/test_DEMO_GUI/nc_files/'
        file_list = os.listdir(exp_dir)

        fnum_min = int(file_list[0 ][len(fhead):len(fhead)+4])
        fnum_max = int(file_list[-1][len(fhead):len(fhead)+4])
        
        for k in range(len(file_list)):
             file_list[k] = exp_dir + file_list[k]   
    
        ### open first netcdf file and get variables inside
        file_1        = file_list[0]
        nc            = Dataset(file_1, mode='r')
        nc_vars       = nc.variables

        varnamelist = OrderedDict()

        #### get all variables with 4 dimensions
        for i in nc_vars:
            var_dimensions = nc_vars[i].dimensions
            if len(var_dimensions) == 4:
                varnamelist[i] = nc_vars[i]  
        ### add combined variables
        varnamelist['comb_phyto'] = len(varnamelist) + 1
        varnamelist['comb_zoo'  ] = len(varnamelist) + 1
        varnamelist['comb_N'    ] = len(varnamelist) + 1
        varnamelist['comb_Fe'   ] = len(varnamelist) + 1
        
        return varnamelist


    def on_draw(self):
        '''Function to draw one plot
        '''
        
        global lon_rho, lat_rho, mask_rho_nan
        global lon_1, lon_2, lat_1, lat_2
        global F1, var_1_name
        global it, ii1, jj1, int_i, int_j
        global tday, date_ref, tday_ref
        
        self.num_plots = "one_plot"
        
        ### plot variable 
        plot_1 = self.axes_1.pcolormesh(lon_rho[jj1,ii1],
                                        lat_rho[jj1,ii1],
                                        F1[:,:,it] * mask_rho_nan[jj1,ii1],
                                        vmin=F1.min(),   # set lower bound of variable range
                                        vmax=F1.max(),   # set upper bound of variable range
                                        shading='gouraud')
        ### set lat/lon axes
        self.axes_1.axis([lon_rho[jj1,ii1].min(),
                          lon_rho[jj1,ii1].max(),
                          lat_rho[jj1,ii1].min(),
                          lat_rho[jj1,ii1].max()])
                          
        self.axes_1.set_xlabel(var_1_combined)
        
        ### show colorbar
        cb = self.fig_1.colorbar(plot_1, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        
        ### show title
        self.axes_1.set_title('%s' %self.dateAK(tday[it],date_ref,tday_ref))
        
        self.canvas_1.draw()
        
        ### show updated plot as it is drawn
        ### this line is what causes the plot to look animated
        QApplication.processEvents()

        ### clear plot
        self.axes_1.cla()
        cb.remove()


    def on_draw_2_plots(self):
        '''Function to plot the surface variable in plot window 1
           and the depth variable in plot window 2
        '''
        
        global lon_rho, lat_rho, mask_rho_nan
        global lon_1, lon_2, lat_1, lat_2
        global F1, F2, zr, var_1_combined, var_2_combined
        global it, ii1, jj1, int_i, int_j
        global tday, date_ref, tday_ref
        global var_1_loaded_from, var_2_loaded_from
        
        if var_1_loaded_from != 'load_var_1_surf':
            self.load_var_1_surf(self.date_start.text(), 
                                 self.date_end.text(), self.lon_east.text(), 
                                 self.lon_west.text(), self.lat_north.text(), 
                                 self.lat_south.text(),self.variable_1.currentText())
            print var_1_loaded_from
        if var_2_loaded_from != 'load_var_2_depth':
            self.load_var_2_depth(self.date_start.text(), 
                                  self.date_end.text(), self.lon_east.text(), 
                                  self.lon_west.text(), self.lat_north.text(), 
                                  self.lat_south.text(),self.variable_2.currentText())
            print var_2_loaded_from    

        self.plot_type = 'surf-depth'
        
        ### set up vertical grid    
        x = [lon_1, lon_2]
        y = [lat_1, lat_2]
        
        ### set up number of vertical contours
        F2_step = (F2.max() - F2.min()) / 20
        OPT = {
            'contours': np.arange(F2.min(), F2.max(), F2_step),
        }

        ########## plot surface variable ##########
        
        if var_1_combined == var_2_combined:    ### if same variables, ensure colorbar encompasses the same range
            plot_1 = self.axes_1.pcolormesh(lon_rho[jj1, ii1],
                                            lat_rho[jj1, ii1],
                                            F1[:, :, it] * mask_rho_nan[jj1, ii1],
                                            vmin = F2.min(),   # set lower bound of variable range
                                            vmax = F2.max(),   # set upper bound of variable range
                                            shading = 'gouraud')
            ### set lat/lon axes
            self.axes_1.axis([lon_rho[jj1,ii1].min(),
                              lon_rho[jj1,ii1].max(),
                              lat_rho[jj1,ii1].min(),
                              lat_rho[jj1,ii1].max()])
                              
            self.axes_1.set_xlabel(var_1_combined)
            
            ### show colorbar
            cb_1 = self.fig_1.colorbar(plot_1, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        
            ### set line denoting cross-section displayed in second plot
            self.axes_1.plot(x, y, 'k-')
        else:
            plot_1 = self.axes_1.pcolormesh(lon_rho[jj1, ii1],
                                            lat_rho[jj1, ii1],
                                            F1[:, :, it] * mask_rho_nan[jj1, ii1],
                                            vmin = F1.min(),   # set lower bound of variable range
                                            vmax = F1.max(),   # set upper bound of variable range
                                            shading = 'gouraud')
            ### set lat/lon axes
            self.axes_1.axis([lon_rho[jj1,ii1].min(),
                              lon_rho[jj1,ii1].max(),
                              lat_rho[jj1,ii1].min(),
                              lat_rho[jj1,ii1].max()])
                              
            self.axes_1.set_xlabel(var_1_combined)
            
            ### show colorbar
            cb_1 = self.fig_1.colorbar(plot_1, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        
            ### set line denoting cross-section displayed in second plot
            self.axes_1.plot(x, y, 'k-')
            
            ########## end surface variable  ##########
        
        
        ########## plot depth variable ##########
        longr = lon_rho[jj1,ii1]
        latgr = lat_rho[jj1,ii1]

        Xi, Z, SECT, Ipos, Jpos, xcoord, ycoord, I, J = self.rnt_section_(longr, latgr, zr,
                                                                          F2[:, :, :, it],
                                                                          x, y, OPT)                                                               
        plot_2 = self.axes_2.pcolormesh(np.asarray(Xi),
                                        np.asarray(Z),
                                        np.asarray(SECT),
                                        vmin = F2.min(),   # set lower bound of variable range
                                        vmax = F2.max(),   # set upper bound of variable range
                                        shading = 'gouraud')
        self.axes_2.set_xlabel(var_2_combined)
        self.axes_2.set_xlim([np.asarray(Xi).min(), np.asarray(Xi).max()])
        self.axes_2.set_ylim([-200, 0])     # set depth
                
        ### show colorbar
        cb_2 = self.fig_2.colorbar(plot_2, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        
        ########## end depth variable ##########
        
        
        ### set titles
        self.axes_1.set_title('%s' %self.dateAK(tday[it], date_ref, tday_ref))
        self.axes_2.set_title('%s' %self.dateAK(tday[it], date_ref, tday_ref))
        
        self.canvas_1.draw()
        self.canvas_2.draw()
        
        ### show updated plots as they are drawn
        ### this line is what causes the plots to look animated
        QApplication.processEvents()
     
        ### clear plots 
        self.axes_1.cla()
        self.axes_2.cla()
        cb_1.remove()
        cb_2.remove() 
        
        
    def draw_surf_surf(self):
        '''Plot both variables at surface
        '''
        global lon_rho, lat_rho, mask_rho_nan
        global lon_1, lon_2, lat_1, lat_2
        global F1, F2, zr, var_1_combined, var_2_combined
        global it, ii1, jj1, int_i, int_j
        global tday, date_ref, tday_ref
        global var_1_loaded_from, var_2_loaded_from
        
        if var_1_loaded_from != 'load_var_1_surf':
            self.load_var_1_surf(self.date_start.text(), 
                                 self.date_end.text(), self.lon_east.text(), 
                                 self.lon_west.text(), self.lat_north.text(), 
                                 self.lat_south.text(),self.variable_1.currentText())
            print var_1_loaded_from
        if var_2_loaded_from != 'load_var_2_surf':
            self.load_var_2_surf(self.date_start.text(), 
                                 self.date_end.text(), self.lon_east.text(), 
                                 self.lon_west.text(), self.lat_north.text(), 
                                 self.lat_south.text(),self.variable_2.currentText())
            print var_2_loaded_from
        
        self.plot_type = 'surf-surf'
        
        ### set up vertical grid    
        x = [lon_1, lon_2]
        y = [lat_1, lat_2]
        
        ### plot first variable
        if var_1_combined == var_2_combined:    ### if same variables, ensure colorbar encompasses the same range 
            plot_1 = self.axes_1.pcolormesh(lon_rho[jj1, ii1],
                                            lat_rho[jj1, ii1],
                                            F1[:, :, it] * mask_rho_nan[jj1, ii1],
                                            vmin = F2.min(),   # set lower bound of variable range
                                            vmax = F2.max(),   # set upper bound of variable range
                                            shading = 'gouraud')
            ### set lat/lon axes
            self.axes_1.axis([lon_rho[jj1,ii1].min(),
                              lon_rho[jj1,ii1].max(),
                              lat_rho[jj1,ii1].min(),
                              lat_rho[jj1,ii1].max()])
                          
            self.axes_1.set_xlabel(var_1_combined)
        
            ### show colorbar
            cb_1 = self.fig_1.colorbar(plot_1, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        else:
            plot_1 = self.axes_1.pcolormesh(lon_rho[jj1, ii1],
                                            lat_rho[jj1, ii1],
                                            F1[:, :, it] * mask_rho_nan[jj1, ii1],
                                            vmin = F1.min(),   # set lower bound of variable range
                                            vmax = F1.max(),   # set upper bound of variable range
                                            shading = 'gouraud')
            ### set lat/lon axes
            self.axes_1.axis([lon_rho[jj1,ii1].min(),
                              lon_rho[jj1,ii1].max(),
                              lat_rho[jj1,ii1].min(),
                              lat_rho[jj1,ii1].max()])
                          
            self.axes_1.set_xlabel(var_1_combined)
        
            ### show colorbar
            cb_1 = self.fig_1.colorbar(plot_1, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
            
        
        ### plot second variable    
        plot_2 = self.axes_2.pcolormesh(lon_rho[jj1, ii1],
                                        lat_rho[jj1, ii1],
                                        F2[:, :, it] * mask_rho_nan[jj1, ii1],
                                        vmin = F2.min(),   # set lower bound of variable range
                                        vmax = F2.max(),   # set upper bound of variable range
                                        shading = 'gouraud')
        ### set lat/lon axes
        self.axes_2.axis([lon_rho[jj1,ii1].min(),
                          lon_rho[jj1,ii1].max(),
                          lat_rho[jj1,ii1].min(),
                          lat_rho[jj1,ii1].max()])
                          
        self.axes_2.set_xlabel(var_2_combined)
        
        ### show colorbar
        cb_2 = self.fig_2.colorbar(plot_2, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')        
        
        ### set titles
        self.axes_1.set_title('%s' %self.dateAK(tday[it], date_ref, tday_ref))
        self.axes_2.set_title('%s' %self.dateAK(tday[it], date_ref, tday_ref))
        
        self.canvas_1.draw()
        self.canvas_2.draw()
        
        ### show updated plots as they are drawn
        ### this line is what causes the plots to look animated
        QApplication.processEvents()
     
        ### clear plots 
        self.axes_1.cla()
        self.axes_2.cla()
        cb_1.remove()
        cb_2.remove()
        
        
    def draw_depth_depth(self):
        '''Plot both variables at depth 
        '''
        global lon_rho, lat_rho, mask_rho_nan
        global lon_1, lon_2, lat_1, lat_2
        global F1, F2, zr, var_1_combined, var_2_combined
        global it, ii1, jj1, int_i, int_j
        global tday, date_ref, tday_ref
        global var_1_loaded_from, var_2_loaded_from
        
        if var_1_loaded_from != 'load_var_1_depth':
            self.load_var_1_depth(self.date_start.text(), 
                                  self.date_end.text(), self.lon_east.text(), 
                                  self.lon_west.text(), self.lat_north.text(), 
                                  self.lat_south.text(),self.variable_1.currentText())
            print var_2_loaded_from
        if var_2_loaded_from != 'load_var_2_depth':
            self.load_var_2_depth(self.date_start.text(), 
                                  self.date_end.text(), self.lon_east.text(), 
                                  self.lon_west.text(), self.lat_north.text(), 
                                  self.lat_south.text(),self.variable_2.currentText())
            print var_2_loaded_from
        
        self.plot_type = 'depth-depth'
        
        ### set up vertical grid    
        x = [lon_1, lon_2]
        y = [lat_1, lat_2]
        
        ### set up number of vertical contours
        F1_step = (F1.max() - F1.min()) / 20
        OPT = {
            'contours': np.arange(F1.min(), F1.max(), F1_step),
        }
        F2_step = (F2.max() - F2.min()) / 20
        OPT = {
            'contours': np.arange(F2.min(), F2.max(), F2_step),
        }
        
        longr = lon_rho[jj1,ii1]
        latgr = lat_rho[jj1,ii1]
        
        ### plot first variable
        if var_1_combined == var_2_combined:    ### if same variables, ensure colorbar encompasses the same range
            Xi, Z, SECT, Ipos, Jpos, xcoord, ycoord, I, J = self.rnt_section_(longr, latgr, zr,
                                                                              F1[:, :, :, it],
                                                                              x, y, OPT)                                                               
            plot_1 = self.axes_1.pcolormesh(np.asarray(Xi),
                                            np.asarray(Z),
                                            np.asarray(SECT),
                                            vmin = F2.min(),   # set lower bound of variable range
                                            vmax = F2.max(),   # set upper bound of variable range
                                            shading = 'gouraud')

            self.axes_1.set_xlim([np.asarray(Xi).min(), np.asarray(Xi).max()]) 
            self.axes_1.set_ylim([-200, 0])     # set depth
                          
            self.axes_1.set_xlabel(var_1_combined)
        
            ### show colorbar
            cb_1 = self.fig_1.colorbar(plot_1, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        else:
            Xi, Z, SECT, Ipos, Jpos, xcoord, ycoord, I, J = self.rnt_section_(longr, latgr, zr,
                                                                              F1[:, :, :, it],
                                                                              x, y, OPT)                                                               
            plot_1 = self.axes_1.pcolormesh(np.asarray(Xi),
                                            np.asarray(Z),
                                            np.asarray(SECT),
                                            vmin = F1.min(),   # set lower bound of variable range
                                            vmax = F1.max(),   # set upper bound of variable range
                                            shading = 'gouraud')

            self.axes_1.set_xlim([np.asarray(Xi).min(), np.asarray(Xi).max()]) 
            self.axes_1.set_ylim([-200, 0])     # set depth
                          
            self.axes_1.set_xlabel(var_1_combined)
        
            ### show colorbar
            cb_1 = self.fig_1.colorbar(plot_1, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        
        ### Plot second variable
        Xi, Z, SECT, Ipos, Jpos, xcoord, ycoord, I, J = self.rnt_section_(longr, latgr, zr,
                                                                          F2[:, :, :, it],
                                                                          x, y, OPT)                                                               
        plot_2 = self.axes_2.pcolormesh(np.asarray(Xi),
                                        np.asarray(Z),
                                        np.asarray(SECT),
                                        vmin = F2.min(),   # set lower bound of variable range
                                        vmax = F2.max(),   # set upper bound of variable range
                                        shading = 'gouraud')
        
        self.axes_2.set_xlabel(var_2_combined)
        self.axes_2.set_xlim([np.asarray(Xi).min(), np.asarray(Xi).max()]) 
        self.axes_2.set_ylim([-200, 0])     # set depth
                
        ### show colorbar
        cb_2 = self.fig_2.colorbar(plot_2, orientation='horizontal', pad=0.1, fraction=0.05, format='%.1f')
        
        ### set titles
        self.axes_1.set_title('%s' %self.dateAK(tday[it], date_ref, tday_ref))
        self.axes_2.set_title('%s' %self.dateAK(tday[it], date_ref, tday_ref))
        
        self.canvas_1.draw()
        self.canvas_2.draw()
        
        ### show updated plots as they are drawn
        ### this line is what causes the plots to look animated
        QApplication.processEvents()
     
        ### clear plots 
        self.axes_1.cla()
        self.axes_2.cla()
        cb_1.remove()
        cb_2.remove()       
        
        
    def on_select_section(self):
        '''Allow user to select two grid points on the surface plot
           and plot the corresponding section at depth in the second
           plot
        '''
        global lon_1, lon_2, lat_1, lat_2
        global var_1_loaded_from, var_2_loaded_from
        
        ### load appropriate plots 
        if var_1_loaded_from != 'load_var_1_surf':
            self.load_var_1_surf(self.date_start.text(), 
                                  self.date_end.text(), self.lon_east.text(), 
                                  self.lon_west.text(), self.lat_north.text(), 
                                  self.lat_south.text(),self.variable_1.currentText())
        if var_2_loaded_from != 'load_var_2_depth':
            self.load_var_2_depth(self.date_start.text(),
                                  self.date_end.text(), self.lon_east.text(),
                                  self.lon_west.text(), self.lat_north.text(),
                                  self.lat_south.text(),self.variable_2.currentText())
        self.on_draw_2_plots()      # show appropiate plots
        
        ### Get new coordinates from user
        lon_1, lat_1 = self.fig_1.ginput()[0]
        self.on_draw()
        
        lon_2, lat_2 = self.fig_1.ginput()[0]
        self.on_draw_2_plots()
        
        print 'new coordinates', lon_1, lat_1, lon_2, lat_2
        
    
    def on_play(self, pressed):
        '''Play data from date_start to date_end
        '''        
        global it, nt
        
        source = self.sender()
        print 'int pressed', int(pressed)
        if pressed:
            print "Plot type:", self.plot_type
            btn_info = str(source.text())
            self.play_button.setText('Pause')
            pressed_2 = False
        else:
            btn_info = 'not_play'
            self.play_button.setText('Play')
            return
        # #print 'button text:', self.sender().text(), type(self.sender().text())
        # btn_text = str(self.sender().text())
        # if btn_text == 'Play':
        #     self._running = True
        #     btn_text = 'Stop'
        # else:
        #     return
        #
        # while self._running:
        while it < nt and pressed == 1: #btn_info == 'Play':
            self.slider.setValue(it)            # update slider 
            #self.connect(self.play_button, SIGNAL('clicked()'), self.on_pause)
            if pressed:
                pressed_2 = False
            if pressed_2:
                btn_info = 'not_play'
                continue
            
            if   self.plot_type == 'surf-surf':
                self.draw_surf_surf()
            elif self.plot_type == 'depth-depth':
                self.draw_depth_depth()
            elif self.plot_type == 'surf-depth':
                self.on_draw_2_plots()    
    
            it += 1

        if it == nt:
            #self.end_movie()
            #self._running = False
            self.play_button.setText('Play')
            self.connect(self.play_button, SIGNAL('clicked()'), self.on_play)
            it = 0      # allow movie to be replayed
            btn_info == 'stop'
        
    
    def on_pause(self):
        '''Pause movie
        '''
        self._running = False
        self.connect(self.play_button, SIGNAL('clicked()'), self.on_play)
        
    def slider_moved(self, position):
        '''Connect the slider to each frame in the movie
        "position" denotes the frame (it) to display
        '''
        global it
        
        print 'slider position:', position
        
        it = val
        
        if   self.plot_type == 'surf-surf':
            self.draw_surf_surf()
        elif self.plot_type == 'depth-depth':
            self.draw_depth_depth()
        elif self.plot_type == 'surf-depth':
            self.on_draw_2_plots()
            
         
    def rnt_section_(self, lonr, latr, zr, field, x, y, OPT):
        '''function to implement rnt_section.m in python
        '''
        varargin = [lonr, latr, zr, field, x, y, OPT]
        nargin   = len(varargin)

        OPT['interp'] = 'none'
        OPT['res']    = 5
        OPT['iplot']  = 1

        K = 0
        
        if OPT['res'] == 0:
            xg = copy_(x)
            yg = copy_(y)
        else:
            for i in range(len(x) - 1):
                x1 = x[i]
                y1 = y[i]
                x2 = x[i + 1]
                y2 = y[i + 1]
                
                a = x2 - x1
                b = y2 - y1
                c = np.sqrt(a * a + b * b)
                dc = 1.0 / (100 / OPT['res'])
                
                theta = np.sign(b) * np.arccos(a / c)
                
                xg = x1
                yg = y1
                
                for c_tmp in np.arange(dc, c, dc):
                    xg = np.append(np.asarray(xg), x1 + c_tmp * np.cos(theta))
                    yg = np.append(np.asarray(yg), y1 + c_tmp * np.sin(theta))
         
        Ipos, Jpos, triX, triY, Ival, Jval = self.rnt_hindicesTRI_(xg, yg, lonr, latr)
        
        Ipos = Ipos[~np.isnan(Ipos)]
        Jpos = Jpos[~np.isnan(Jpos)]
        I    = np.round(Ipos) 
        J    = np.round(Jpos)
        
        arr        = np.array([I, J], dtype=int)
        SIZ1, SIZ2 = np.shape(lonr)
        IND        = np.ravel_multi_index(arr, dims=(SIZ1, SIZ2), mode='clip', order='C')

        blon = np.zeros((np.shape(IND)))
        for i in range(len(IND)):
            blon[i] = lonr.item(IND[i])
            
        blat = np.zeros((np.shape(IND)))
        for i in range(len(IND)):
            blat[i] = latr.item(IND[i])
        
        distances = self.rnt_earthdist_(blon[0], blat[0], blon, blat)

        if OPT['interp'] == 'none':
            Z    = np.zeros((np.shape(IND)[0], np.shape(field)[2]))
            SECT = np.zeros((np.shape(IND)[0], np.shape(field)[2]))
            X    = np.zeros((np.shape(IND)[0], np.shape(field)[2]))
            
            for k in range(np.shape(field)[2]):
                tmp       = zr[:,:,k]
                for i in range(len(IND)):
                    Z[:,k][i] = tmp.item(IND[i])
                tmp       = field[:,:,k]
                for i in range(len(IND)):
                    SECT[:,k][i] = tmp.item(IND[i])
                X[:,k]    = distances

        else:
            II,JJ=size_(lonr,nargout=2)
            II,JJ=meshgrid_(arange_(1,II),arange_(1,JJ),nargout=2)
            tmp=interp2_(II,JJ,lonr.T,Ipos,Jpos,OPT.interp)
            lonr1=tmp[:]
            tmp=interp2_(II,JJ,latr.T,Ipos,Jpos,OPT.interp)
            latr1=tmp[:]
            distances=rnt_earthdist_(lonr1[1],latr1[1],lonr1,latr1)
            for k in arange_(1,size_(field,3)).reshape(-1):
                tmp=zr[:,:,k]
                tmp=interp2_(II,JJ,tmp.T,Ipos,Jpos,OPT.interp)
                Z[:,k]=tmp[:]
                tmp=field[:,:,k]
                tmp=interp2_(II,JJ,tmp.T,Ipos,Jpos,OPT.interp)
                SECT[:,k]=tmp[:]
                X[:,k]=distances
                
        if OPT['iplot'] == 1:
            X = np.real(X) / 1000   # dividing by 1000 sets x limit of depth plots in km
            
            h_bottom = Z[:, 0].T 
            x_coord  = X[:, 0].T / 1000
            xr       = x_coord
            
            x_coord  = np.append(x_coord , [x_coord[-1], x_coord[0], x_coord[0]])
            h_bottom = np.append(h_bottom, [np.amin(h_bottom) - 10, np.amin(h_bottom) - 10, h_bottom[0]])
        
        i    = np.zeros((np.shape(Ipos)))
        j    = np.zeros((np.shape(Jpos)))

        i[0] = np.round(Ipos[0])
        j[0] = np.round(Jpos[0])
        n    = 0
        
        for k in np.arange(2, len(Ipos)):
            i_tmp = np.round(Ipos[k])
            j_tmp = np.round(Jpos[k])
            if i_tmp != i[n] or j_tmp != j[n]:
                n = n + 1
                i[n] = i_tmp
                j[n] = j_tmp
        I = i
        J = j
        
        return X, Z, SECT, Ipos, Jpos, xg, yg, I, J
           
       
    def rnt_hindicesTRI_(self, Xpos, Ypos, Xgrd, Ygrd):
        '''function to implement rnt_hindicesTRI.m in python 
        '''
        
        I, J   = (np.shape(Xgrd))
        I1, J1 = (np.shape(Ygrd))

        if I1 != I or J != J1:
            print 'RNT_HINDICES - Inconsistent size of grid arrays - STOP'
            return
            
        ### set output matrix according to sizes
        m_in = np.shape(Xpos)[0]
        n_in = 1
        
        Npos = m_in * n_in
        
        ### set dimensions of parameter arrays
        Ipos = np.zeros((Npos))
        Jpos = np.zeros((Npos))
        triX = np.zeros((Npos, 3))
        triY = np.zeros((Npos, 3))
        Ival = np.zeros((Npos, 3))
        Jval = np.zeros((Npos, 3))
        
        ### get size of Xgrd
        Lp, Mp = np.shape(Xgrd)

        Ipos, Jpos, triX, triY, Ival, Jval = interp.hindices(Ipos, Jpos, Xpos, Ypos,
                                                             Xgrd, Ygrd, triX, triY,
                                                             Ival, Jval, Npos, Lp, Mp)
        Ival = Ival + 1
        Jval = Jval + 1
        L    = len(Xpos)
        xi   = Xpos
        yi   = Ypos
        w    = np.zeros((Npos, 3))    
        
        del_ = (triX[:,1] - triX[:,0])*(triY[:,2] - triY[:,0]) - \
               (triX[:,2] - triX[:,0])*(triY[:,1] - triY[:,0])
        
        w[:, 2] = (((triX[:,0] - xi)*(triY[:,1] - yi)) - \
                   ((triX[:,1] - xi)*(triY[:,0] - yi))) / del_
        w[:, 1] = (((triX[:,2] - xi)*(triY[:,0] - yi)) - \
                   ((triX[:,0] - xi)*(triY[:,2] - yi))) / del_
        w[:, 0] = (((triX[:,1] - xi)*(triY[:,2] - yi)) - \
                   ((triX[:,2] - xi)*(triY[:,1] - yi))) / del_
        
        Ipos = np.sum(Ival*w, axis=1)
        Jpos = np.sum(Jval*w, axis=1)
        
        return Ipos, Jpos, triX, triY, Ival, Jval
        
             
    def rnt_earthdist_(self, alon, alat, blon, blat):
        '''rnt_earthdist -- Distance in meters between two lon/lats.
        earthdist(alon, aloat, blon, blat, radius) returns the
        distance in maters between locations (alon, alat)
        and loncations (blon, blat).  The default earth radius is
        assumed to be 6371*1000 meters, the radius for
        a sphere of equal-volume.
    
        from Chuck..
        
        Edited by zwallace 28 October 2015
        '''
    
        varargin = [alon, alat, blon, blat]
        nargin   = len(varargin)

        radius = 6371 * 1000 
        RCF    = 180 / np.pi
        
        alon = alon / RCF
        alat = alat / RCF
        blon = blon / RCF
        blat = blat / RCF
        
        c    = np.cos(alat)
        ax   = c * (np.cos(alon))
        ay   = c * (np.sin(alon))
        az   = np.sin(alat)
        
        c    = np.cos(blat)
        bx   = c * (np.cos(blon))
        by   = c * (np.sin(blon))
        bz   = np.sin(blat)
        
        result = np.arccos(ax*bx + ay*by + az*bz) * radius
        
        return result
        
        
    def on_reset_button(self):
        '''Reset plots
        '''
        self.fig_1.clf()
        self.fig_2.clf()
        self.create_plot_window()
        
        

def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()