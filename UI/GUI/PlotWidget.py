#import gi

import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_gtk3agg import (FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.colors import ListedColormap

from scipy.optimize import curve_fit

class plotwidget(object):
    def __init__(self, data_queue):
        self.plottype = 'normal'
        self.fig = Figure(figsize = (5, 5), dpi = 100)
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_xlabel('X', size = 12)
        self.ax.set_ylabel('Y', size = 12)
        self.fig.subplots_adjust(left = 0.15, top = 0.85)
        self.ax.axis([0, 255, 0, 255])
        self.x_vals = np.empty(0, np.uint16)
        self.y_vals = np.empty(0, np.uint16)
        self.t_vals = np.empty(0, np.uint16)
        self.intensity = np.empty(0, np.uint16)
        self.length = np.empty(0, np.uint16)
        self.i = 0
        self.occupancy_array = np.array([[0, 0, 0]])
        self.occ_length = []
        self.j = 0
        self.old_elements = np.array([[0, 0]])
        self.occupancy = []
        self.integration_length = 500
        self.color_depth = 10
        cmap = self.fading_colormap(5)
        self.data_queue = data_queue

        self.scatter = self.ax.scatter(self.x_vals,
                                       self.y_vals,
                                       c = [],
                                       s = 1,
                                       marker = 's',
                                       cmap = cmap,
                                       vmin = 0,
                                       vmax = 1,
                                       animated=True)

        self.canvas = FigureCanvas(self.fig)
        self.canvas.set_size_request(500, 500)
        #self.ax.plot()
        self.canvas.draw()
        self.fig_bgr = self.canvas.copy_from_bbox(self.fig.bbox)
        # my bayan
        # -----------------------------
        self.ax.draw_artist(self.scatter)
        self.canvas.blit(self.fig.bbox)
        plt.pause(0.1)
        # -----------------------------

    def fading_colormap(self, steps = 5):
        # This creates a fading colormap of 'steps' steps. Each step is more transparent,
        # so if plotted correctly the data appears to be fading out.
        if steps <= 0:
            self.colorsteps = 1
            print('ERROR: Minimum number of colorsteps is 1. Colorsteps have been set to 1.')
        else:
            self.colorsteps = steps

        i = 1
        viridis = cm.get_cmap('viridis', 256)
        newcmap = viridis(np.linspace(0, 1, 256))
        newmap1 = np.tile(newcmap, (self.colorsteps, 1))
        while(i<self.colorsteps):
            newmap1[(i - 1) * 256:(i * 256), -1] = np.linspace(i * 1 / self.colorsteps, i * 1 / self.colorsteps, 256)
            i = i + 1
        cmap = ListedColormap(newmap1)

        return cmap

    def get_new_vals(self):
        #Get values from Chip
        #t need to between 0 and 1 then the calculation 1 - (t / self.colorsteps) needs
        #to be done in order to distribute is correctly over the colormap
        x = np.empty(0, np.uint16)
        y = np.empty(0, np.uint16)
        t = np.empty(0, np.uint16)

        if not self.data_queue.empty():
            #x_new, y_new, t_new = self.data_queue.get()
            x_new, y_new, t_new, _ = self.data_queue.get()
            x = np.append(x, x_new)
            y = np.append(y, y_new)
            t = np.append(t, t_new)
        if len(t) > 0:
            self.t_vals = np.append(self.t_vals, t)
            max_value = np.amax(self.t_vals) # This has to be changed to a more general way
            t = t/max_value

        return x, y, t

    def update_plot(self):
        # some temp checks 
        #print("UI::GUI::PlotWidget::update_plot") 
        print("UPD plot: len-of-queue {}, ".format(self.data_queue.qsize()),flush=True)
        print("UPD plot: arrays {},{},{} ".format(self.x_vals.shape, self.y_vals.shape, self.t_vals.shape),flush=True)
        print("UPD plot: self.length:{}, {} ".format(len(self.length),self.length[0:4]),flush=True)

        #Plot the fading plot with new data.
        new_xvals, new_yvals, new_tvals = self.get_new_vals()
        self.x_vals = np.append(self.x_vals, new_xvals)
        self.y_vals = np.append(self.y_vals, new_yvals)
        self.length = np.append(self.length, new_xvals.size)

        #print("UPD plot: self.new_xvals.size: {}, self.i = {}".format(new_xvals.shape[0], self.i),flush=True)
        #Cut plotting arrays to n_colorsteps Timeblocks
        if self.i < (self.colorsteps):
            #print("UPD plot: self.i({}) < (self.colorsteps)({})".format(self.i,self.colorsteps))
            self.i = self.i + 1

        elif self.i == (self.colorsteps):
            #print("UPD plot: self.i({}) == (self.colorsteps ({}))".format(self.i, self.colorsteps),flush=True)
            number = np.arange(self.length[0])
            self.length = np.delete(self.length, 0)
            self.x_vals = np.delete(self.x_vals, number)
            self.y_vals = np.delete(self.y_vals, number)
            self.t_vals = np.delete(self.t_vals, number)
            self.intensity = np.delete(self.intensity, number)

        elif self.i > (self.colorsteps):
            #print("UPD plot: self.i({}) > (self.colorsteps ({}))".format(self.i, self.colorsteps),flush=True)

            while self.i >= (self.colorsteps):
                number = np.arange(self.length[0])
                self.length = np.delete(self.length, 0)
                self.x_vals = np.delete(self.x_vals, number)
                self.y_vals = np.delete(self.y_vals, number)
                self.t_vals = np.delete(self.t_vals, number)
                self.intensity = np.delete(self.intensity, number)
                self.i = self.i-1

        #print("x->{}".format(self.x_vals[0:16]),flush=True)
        #print("y->{}".format(self.y_vals[0:16]),flush=True)
        #huya = np.c_[self.x_vals, self.y_vals]
        #print("np.c_[x,y]->size={}".format(huya.size),flush=True)
        #print("np.c_[x,y][0:16] - > {}".format(huya[0:16]),flush=True)
        #print("np.c_[x,y]-> size->{}".format(huya.size()),flush=True)

        xy_arr = np.c_[self.x_vals,self.y_vals]

        #if np.c_[self.x_vals, self.y_vals].size != 0:
        #       ^ if x & y are not the same length/size -> exception!
        if xy_arr.size != 0:
            
            #print("UPD plot: np.c_[..].size!=0",flush=True)

            #self.canvas.restore_region(self.fig_bgr)

            #self.scatter.set_offsets(np.c_[self.x_vals, self.y_vals])
            self.scatter.set_offsets(xy_arr)
            self.intensity = np.concatenate((np.array(self.intensity) - (1 / self.colorsteps), new_tvals))
            self.scatter.set_array(self.intensity)

            #self.canvas.draw()
            # -------------------------------
            
            self.canvas.restore_region(self.fig_bgr)
            self.ax.draw_artist(self.scatter)
            self.canvas.blit(self.fig.bbox)
            self.canvas.flush_events()            

            xy_arr = None

        return True

    def update_occupancy_plot(self):
        new_xvals, new_yvals, new_tvals = self.get_new_vals()
        new_elements = np.c_[new_xvals, new_yvals]
        self.occ_length.append(new_elements.shape[0])
        self.old_elements = np.append(self.old_elements, new_elements, axis = 0)

        #count hited pixel
        for new_element in new_elements:
            pos = np.argwhere(np.all(self.occupancy_array[ : , :2] == new_element, axis = 1) == True)
            if pos.size == 0:
                # add new element
                self.occupancy_array = np.append(self.occupancy_array, [np.append(new_element, 1)], axis = 0)

            elif pos.size == 1:
                #increment element at pos
                x = pos[0, 0]
                self.occupancy_array[pos[0, 0], 2] = (self.occupancy_array[pos[0, 0], 2] + 1)

            else:
                print('Error')

        #remove hitted pixel
        if self.j <= (self.integration_length):
            self.j = self.j + 1

        elif self.j == (self.integration_length + 1):
            number = self.occ_length[0]
            self.occ_length.pop(0)
            k = 0
            while k < number:
                k = k + 1
                pos = np.argwhere(np.all(self.occupancy_array[ : , :2] == self.old_elements[1], axis = 1) == True)
                self.old_elements = np.delete(self.old_elements, 1, axis = 0)
                if self.occupancy_array[pos[0, 0], 2] == 1:
                    #Remove item if no count left
                    self.occupancy_array = np.delete(self.occupancy_array, pos[0, 0], axis = 0)

                elif self.occupancy_array[pos[0, 0], 2] > 1:
                    #decrement element at pos
                    self.occupancy_array[pos[0, 0], 2] = (self.occupancy_array[pos[0, 0], 2] - 1)

                else:
                    print('Error')

        elif self.j > (self.integration_length + 1):
            while self.j > (self.integration_length + 1):
                number = self.occ_length[0]
                self.occ_length.pop(0)
                k = 0
                while k < number:
                    k = k + 1
                    pos = np.argwhere(np.all(self.occupancy_array[ : , :2] == self.old_elements[1], axis = 1) == True)
                    self.old_elements = np.delete(self.old_elements, 1, axis = 0)
                    if self.occupancy_array[pos[0, 0], 2] == 1:
                        #Remove item if no count left
                        self.occupancy_array = np.delete(self.occupancy_array,pos[0, 0], axis = 0)

                    elif self.occupancy_array[pos[0, 0], 2] > 1:
                        #Decrement element at pos
                        self.occupancy_array[pos[0, 0], 2] = (self.occupancy_array[pos[0, 0], 2] - 1)

                    else:
                        print('Error')

                self.j = self.j - 1

        if self.occupancy_array.size > 3:
            self.scatter.set_offsets(self.occupancy_array[ : , :2])
            self.occupancy = self.occupancy_array[ : , 2: ]
            self.scatter.set_array(np.squeeze(self.occupancy))
            self.canvas.draw()

        return True

    def reset_occupancy(self):
        self.occupancy_array = np.array([[0, 0, 0]])
        self.occ_length = []
        self.j = 0
        self.old_elements = np.array([[0, 0]])
        self.occupancy = []

        return True

    def change_colormap(self, colormap, vmin = 0, vmax = 1):
        x_vals = []
        y_vals = []
        cmap = colormap
        vmin = vmin
        vmax = vmax
        if self.plottype == 'occupancy':
            self.color_depth = vmax

        self.ax.remove()
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_xlabel('X', size = 12)
        self.ax.set_ylabel('Y', size = 12)
        self.ax.axis([0, 255, 0, 255])
        self.scatter = self.ax.scatter(x_vals, y_vals, c = [], s = 1, marker = 's', cmap = cmap, vmin = vmin, vmax = vmax)
        self.ax.plot()

        return True

    def get_iteration_depth(self,function = 'normal'):
        function = function
        if function == 'normal':
            return self.colorsteps

        elif function == 'occupancy':
            return self.integration_length

        elif function == 'occupancy.color':
            return self.color_depth

        else:
            print('Unknown argument. Use "normal", "occupancy" or "occupancy.color')
            return False

    def get_plottype(self):
        return self.plottype

    def set_plottype(self, plottype):
        self.plottype = plottype
        return True

    def set_occupancy_length(self, length):
        self.integration_length = length
        return True

    def set_color_depth(self, color_depth):
        self.color_depth = color_depth
        return True

    def set_color_steps(self, color_steps):
        self.colorsteps = color_steps
        return True

class TOTplot(object):
    def __init__(self, data_queue, upd_data, upd_rate, readout_interval):
        
        # adding another subplot
        self.figtot = Figure(figsize=(5,5),dpi=100)
        self.ax_tot = self.figtot.add_subplot(111)
        #self.canvas_tot.draw()
        #self.figtot.subplots_adjust(left = 0.2, top = 0.9)
        self.ax_tot.set_xlabel('TOT')
        self.ax_tot.set_ylabel('Entries')
        self.data_queue = data_queue
        # histo parameters
        self.nbins = 51
        self.bin_range = (0,1025)
        self.bin_edges = np.linspace(self.bin_range[0],self.bin_range[1],self.nbins+1)
        self.bin_cnt = np.zeros(self.nbins, dtype=np.int64)
        # other params
        self.logScale = False
        ###################################

        #self.update_callback = upd_data
        self.update_tot_callback = upd_data

        self.update_rate_callback = upd_rate
        self.rate_list = []
        self.t_interval = readout_interval
        #print("len(bin_edges):{} -> first10:{}, bin_cnt[0:9]:{}".format(len(self.bin_edges),self.bin_edges[0:9],self.bin_cnt[0:9]))

        # some test plot here
        self.tot_histo = self.ax_tot.hist(self.bin_edges[:-1],
                                          weights=self.bin_cnt,
                                          range=self.bin_range,
                                          bins=self.nbins,
                                          align='left',
                                          #density = True,# makes y axis to show relative integral values and not absolute ones
                                          histtype='stepfilled',
                                          edgecolor='black',
                                          facecolor='g'
                                          )

        self.canvas_tot = FigureCanvas(self.figtot)
        self.canvas_tot.set_size_request(500,500)
        self.ax_tot.plot()

    def get_tot_val(self):

        tot_val = np.empty(0,np.uint16)

        if not self.data_queue.empty():
            _,_,tot_val, _ = self.data_queue.get()
        return tot_val

    def upd_histo(self):
        
        new_tot = self.get_tot_val()
        print(f"new_tot = {new_tot.shape[0]}")
        #n_tot_hits = len(new_tot)
        avg_hit_rate = self.get_rate(new_tot.shape[0])

        self.update_rate_callback(avg_hit_rate)        

        i_bin_cnt, _ = np.histogram(new_tot, bins=self.bin_edges)
        self.bin_cnt += i_bin_cnt

        self.ax_tot.cla()
        if(len(self.bin_cnt>0) and self.logScale):
            self.ax_tot.set_yscale('log')

        self.ax_tot.set_xlabel('TOT')
        self.ax_tot.set_ylabel('Entries')
        self.tot_histo = self.ax_tot.hist(self.bin_edges[:-1],
                                          weights=self.bin_cnt,
                                          bins=self.nbins,
                                          range=self.bin_range,
                                          align='left',
                                          #density = True,
                                          histtype='stepfilled',
                                          edgecolor='black',
                                          facecolor='g'
                                          )

        #self.canvas_tot.draw()
        dist_mean = 0.0
        if(np.sum(self.bin_cnt)>0):
           dist_mean = np.average(np.linspace(0,1025,self.nbins), weights=self.bin_cnt)
           self.update_tot_callback(dist_mean)

           # works...
           #popt, pcov = curve_fit(self.exponent, np.linspace(0,1025,self.nbins), self.bin_cnt)  
           popt, _ = curve_fit(self.exponent, np.linspace(0,1025,self.nbins), self.bin_cnt)  
           #print(f"Fit for TOT: results:{popt}")
           self.ax_tot.plot(np.linspace(0.1025,self.nbins), self.exponent(np.linspace(0.1025,self.nbins), *popt), 'r-')

        self.canvas_tot.draw()

        return True

    def reset_histo(self):
        
        self.bin_cnt = np.zeros(self.nbins, dtype=np.int64)
        self.update_tot_callback(0.0)
        self.update_rate_callback(0.0)

    def setLogScale(self):
    
        self.logScale = True

    def exponent(self, x, a, b, c):
        return (a * np.exp(-b * x) + c)

    def get_rate(self,nhits):

        maxlen = 10
        rate = 0.0

        if(len(self.rate_list)<maxlen):
            self.rate_list.append(nhits/self.t_interval)
            rate = sum(self.rate_list)/len(self.rate_list)
                
        else:
            self.rate_list.pop(0)
            self.rate_list.append(nhits/self.t_interval)
            rate = sum(self.rate_list)/len(self.rate_list)
        
        #print(f"rate_list={self.rate_list} -> averaged={rate}")

        return rate

class TOAplot(object):
    def __init__(self, data_queue):
        
        # adding another subplot
        self.figtoa = Figure(figsize=(5,5),dpi=100)
        self.ax_toa = self.figtoa.add_subplot(111)
        self.ax_toa.set_xlabel('TOA')
        self.ax_toa.set_ylabel('Entries')
        self.data_queue = data_queue
        # histo parameters
        # for reasons currently unknown i do not get TOA above ~4000 cts...
        #self.bin_range = (0,16400) 
        self.bin_range = (0,4097)
        self.nbins = 50
        self.bin_edges = np.linspace(self.bin_range[0],self.bin_range[1],self.nbins+1)
        self.bin_cnt = np.zeros(self.nbins, dtype=np.int64)
        # other params
        self.logScale = False

        ###################################
        #self.fit_callback = upd_data
        ###################################
        # some test plot here
        self.toa_histo = self.ax_toa.hist(self.bin_edges[:-1],
                                          weights=self.bin_cnt,
                                          range=self.bin_range,
                                          bins=self.nbins,
                                          align='left',
                                          #density = True,# makes y axis to show relative integral values and not absolute ones
                                          #histtype='stepfilled',
                                          histtype='step',
                                          edgecolor='green'
                                          #facecolor='g' # need if using 'stepfilled' type
                                          )

        self.canvas_toa = FigureCanvas(self.figtoa)
        self.canvas_toa.set_size_request(500,500)
        self.ax_toa.plot()

    def get_toa_val(self):

        toa_val = np.empty(0,np.uint16)

        if not self.data_queue.empty():
            _, _, _, toa_val = self.data_queue.get()
            #print(f" got TOA={toa_val} from queue!")
        return toa_val

    def upd_histo(self):
        
        new_toa = self.get_toa_val()

        i_bin_cnt, _ = np.histogram(new_toa, bins=self.bin_edges)
        self.bin_cnt += i_bin_cnt

        self.ax_toa.cla()
        if(len(self.bin_cnt>0) and self.logScale):
            self.ax_toa.set_yscale('log')

        self.ax_toa.set_xlabel('TOA')
        self.ax_toa.set_ylabel('Entries')
        self.toa_histo = self.ax_toa.hist(self.bin_edges[:-1],
                                          weights=self.bin_cnt,
                                          bins=self.nbins,
                                          range=self.bin_range,
                                          align='left',
                                          #density = True,
                                          histtype='stepfilled',
                                          edgecolor='black',
                                          facecolor='b'
                                          )

        self.canvas_toa.draw()

        #if(np.sum(self.bin_cnt)>0):
        
        return True

    def reset_histo(self):
        
        self.bin_cnt = np.zeros(self.nbins, dtype=np.int64)
        

    def setLogScale(self):
    
        self.logScale = True


