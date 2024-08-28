#import gi

import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_gtk3agg import (FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.colors import ListedColormap

class plotwidget(object):
    def __init__(self, data_queue):
        self.plottype = 'normal'
        self.fig = Figure(figsize = (5, 5), dpi = 100)
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_xlabel('X', size = 12)
        self.ax.set_ylabel('Y', size = 12)
        #self.ax.yaxis.set_label_coords()
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

        self.scatter = self.ax.scatter(self.x_vals, self.y_vals, c = [], s = 1, marker = '.', cmap = cmap, vmin = 0, vmax = 1)
        #                                                                                 ^ this one was "s"
        self.canvas = FigureCanvas(self.fig)
        self.canvas.set_size_request(500, 500)
        self.ax.plot()

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
            x_new, y_new, t_new = self.data_queue.get()
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
        print("UI::GUI::PlotWidget::update_plot") 
        print("UPD plot: len-of-queue {}, ".format(self.data_queue.qsize()),flush=True)
        print("UPD plot: arrays {},{},{} ".format(self.x_vals.shape, self.y_vals.shape, self.t_vals.shape),flush=True)
        #print("UPD plot: len-of-queue {}, ".format(self.data_queue.qsize()),flush=True)

        #Plot the fading plot with new data.
        new_xvals, new_yvals, new_tvals = self.get_new_vals()
        self.x_vals = np.append(self.x_vals, new_xvals)
        self.y_vals = np.append(self.y_vals, new_yvals)
        self.length = np.append(self.length, new_xvals.size)

        #Cut plotting arrays to n_colorsteps Timeblocks
        if self.i < (self.colorsteps):
            self.i = self.i + 1

        elif self.i == (self.colorsteps):
            number = np.arange(self.length[0])
            self.length = np.delete(self.length, 0)
            self.x_vals = np.delete(self.x_vals, number)
            self.y_vals = np.delete(self.y_vals, number)
            self.t_vals = np.delete(self.t_vals, number)
            self.intensity = np.delete(self.intensity, number)

        elif self.i > (self.colorsteps):
            while self.i >= (self.colorsteps):
                number = np.arange(self.length[0])
                self.length = np.delete(self.length, 0)
                self.x_vals = np.delete(self.x_vals, number)
                self.y_vals = np.delete(self.y_vals, number)
                self.t_vals = np.delete(self.t_vals, number)
                self.intensity = np.delete(self.intensity, number)
                self.i = self.i-1

        if np.c_[self.x_vals, self.y_vals].size != 0:
            self.scatter.set_offsets(np.c_[self.x_vals, self.y_vals])
            self.intensity = np.concatenate((np.array(self.intensity) - (1 / self.colorsteps), new_tvals))
            self.scatter.set_array(self.intensity)

            self.canvas.draw()

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
    def __init__(self, data_queue):
        
        # adding another subplot
        self.localcnt = 0
        self.figtot = Figure(figsize=(5,5),dpi=100)
        self.ax_tot = self.figtot.add_subplot(111)
        #self.figtot.subplots_adjust(left = 0.2, top = 0.9)
        self.ax_tot.set_xlabel('TOT')
        self.ax_tot.set_ylabel('Entries')
        self.ax_tot.set_yscale('log')
        self.data_queue = data_queue
        self.tot_array = np.array([], dtype=np.uint16)
        # some test plot here
        self.tot_histo = self.ax_tot.hist(self.tot_array,
                                          bins=103,
                                          range=(0,1030), 
                                          #density = True,
                                          histtype='stepfilled', 
                                          facecolor='g'
                                          )
        self.canvas_tot = FigureCanvas(self.figtot)
        self.canvas_tot.set_size_request(500,500)
        self.ax_tot.plot()

    def get_tot_val(self):

        tot_val = np.empty(0,np.uint16)

        if not self.data_queue.empty():
            _,_,tot_val = self.data_queue.get()
        return tot_val

    def upd_histo(self):
        
        new_tot = self.get_tot_val()
        self.tot_array = np.concatenate([self.tot_array, new_tot])
        print("UI::GUI::PlotWidget::TOTplot: TOT array size: {}, last entries-> {}".format(
            self.tot_array.shape[0],
            self.tot_array[len(self.tot_array)-32:len(self.tot_array)])
        ,flush=True)

        self.tot_histo = self.ax_tot.hist(self.tot_array,
                                          bins=103,
                                          range=(0,1030), 
                                          #density = True,
                                          histtype='stepfilled', 
                                          facecolor='g'
                                          )
        self.canvas_tot.draw()
        self.localcnt +=1

        if(self.localcnt == 20):
           np.savetxt("ebala.txt", self.tot_array, delimiter=',')

        return True


