#
# ------------------------------------------------------------
# Copyright (c) All rights reserved
# SiLab, Institute of Physics, University of Bonn
# ------------------------------------------------------------
#

'''
    This script scans over different thresholds for one testpulse height
'''
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from tqdm import tqdm
import numpy as np
import time
import tables as tb
import math

from tpx3.scan_base import ScanBase
import tpx3.analysis as analysis
import tpx3.plotting as plotting
import tpx3.utils as utils
from six.moves import range

local_configuration = {
    # Scan parameters
    'mask_step'        : 16,
    'Vthreshold_start' : 1800,
    'Vthreshold_stop'  : 2800,
    'n_injections'     : 100,
    'thrfile'         : './output_data/20200401_160123_mask.h5'
}


class ThresholdScan(ScanBase):

    scan_id      = "ThresholdScan"
    wafer_number = 0
    y_position   = 0
    x_position   = 'A'

    def scan(self, Vthreshold_start=0, Vthreshold_stop=2911, n_injections=100, tp_period = 1, mask_step=16, progress = None, status = None, **kwargs):
        '''
            Takes data for threshold scan in a range of threshold and a given number of testpulses per iteration
            If progress is None a tqdm progress bar is used else progress should be a Multiprocess Queue which stores the progress as fraction of 1
            If there is a status queue information about the status of the scan are put into it
        '''

        # Check if parameters are valid before starting the scan
        if Vthreshold_start < 0 or Vthreshold_start > 2911:
            raise ValueError("Value {} for Vthreshold_start is not in the allowed range (0-2911)".format(Vthreshold_start))
        if Vthreshold_stop < 0 or Vthreshold_stop > 2911:
            raise ValueError("Value {} for Vthreshold_stop is not in the allowed range (0-2911)".format(Vthreshold_stop))
        if Vthreshold_stop <= Vthreshold_start:
            raise ValueError("Value for Vthreshold_stop must be bigger than value for Vthreshold_start")
        if n_injections < 1 or n_injections > 65535:
            raise ValueError("Value {} for n_injections is not in the allowed range (1-65535)".format(n_injections))
        if mask_step not in {4, 16, 64, 256}:
            raise ValueError("Value {} for mask_step is not in the allowed range (4, 16, 64, 256)".format(mask_step))

        for chip in self.chips[1:]:
            print('Current chip: ' + str(chip))
            # Set general configuration registers of the Timepix3
            self.chips[0].write(chip.write_general_config(write=False))

            # Write to the test pulse registers of the Timepix3
            # Write to period and phase tp registers
            self.chips[0].write(chip.write_tp_period(tp_period, 0, write=False))

            # Write to pulse number tp register
            self.chips[0].write(chip.write_tp_pulsenumber(n_injections, write=False))

        self.logger.info('Preparing injection masks...')
        if status != None:
            status.put("Preparing injection masks")

        # Get the shutter sleep time
        sleep_time = self.get_shutter_sleep_time(tp_period = tp_period, n_injections = n_injections)

        # Create the masks for all steps
        mask_cmds = self.create_scan_masks(mask_step, progress = progress)

        # Start the scan
        self.logger.info('Starting scan...')
        if status != None:
            status.put("Starting scan")
        if status != None:
            status.put("iteration_symbol")
        thresholds = utils.create_threshold_list(utils.get_coarse_jumps(Vthreshold_start, Vthreshold_stop))

        if progress == None:
            # Initialize progress bar
            pbar = tqdm(total=len(mask_cmds) * len(thresholds))
        else:
            # Initialize counter for progress
            step_counter = 0

        scan_param_id = 0
        for threshold in thresholds:
            for chip in self.chips[1:]:
                # Set the threshold
                self.chips[0].write(chip.set_dac("Vthreshold_coarse", int(threshold[0]), write=False))
                self.chips[0].write(chip.set_dac("Vthreshold_fine", int(threshold[1]), write=False))

            with self.readout(scan_param_id=scan_param_id):
                step = 0
                for mask_step_cmd in mask_cmds:
                    for chip in self.chips[1:]:
                        # Only activate testpulses for columns with active pixels
                        self.chips[0].write(chip.write_ctpr(list(range(step//(mask_step//int(math.sqrt(mask_step))), 256, mask_step//int(math.sqrt(mask_step)))), write=False))

                    # Write the pixel matrix for the current step plus the read_pixel_matrix_datadriven command
                    self.chips[0].write(mask_step_cmd)

                    # Open the shutter, take data and update the progress bar
                    with self.shutter():
                        time.sleep(sleep_time)
                        if progress == None:
                            # Update the progress bar
                            pbar.update(1)
                        else:
                            # Update the progress fraction and put it in the queue
                            step_counter += 1
                            fraction = step_counter / (len(mask_cmds) * len(thresholds))
                            progress.put(fraction)
                    for chip in self.chips[1:]:
                        self.chips[0].write(chip.stop_readout(write=False))
                        time.sleep(0.001)
                    step += 1
                for chip in self.chips[1:]:
                    self.chips[0].write(chip.reset_sequential(write=False))
                    time.sleep(0.001)
            scan_param_id += 1

        if progress == None:
            # Close the progress bar
            pbar.close()

        if status != None:
            status.put("iteration_finish_symbol")

        self.logger.info('Scan finished')

    def analyze(self, progress = None, status = None, **kwargs):
        '''
            Analyze the data of the scan
            If progress is None a tqdm progress bar is used else progress should be a Multiprocess Queue which stores the progress as fraction of 1
            If there is a status queue information about the status of the scan are put into it
        '''

        h5_filename = self.output_filename + '.h5'

        self.logger.info('Starting data analysis...')
        if status != None:
            status.put("Performing data analysis")

        # Open the HDF5 which contains all data of the scan
        with tb.open_file(h5_filename, 'r+') as h5_file:
            # Read raw data, meta data and configuration parameters
            raw_data       = h5_file.root.raw_data[:]
            meta_data      = h5_file.root.meta_data[:]
            run_config     = h5_file.root.configuration.run_config[:]
            general_config = h5_file.root.configuration.generalConfig[:]
            op_mode        = [row[1] for row in general_config if row[0]==b'Op_mode'][0]
            vco            = [row[1] for row in general_config if row[0]==b'Fast_Io_en'][0]
            # 'Simulate' more chips
            chip_IDs_new = [b'W18-K7',b'W18-K7',b'W17-D8',b'W17-D8',b'W14-E9', b'W14-E9',b'W15-C5', b'W15-C5']
            for new_Id in range(8):
                h5_file.root.configuration.links.cols.chip_id[new_Id] = chip_IDs_new[new_Id]

            # Get link configuration
            link_config = h5_file.root.configuration.links[:]
            print(link_config)
            chip_IDs    = link_config['chip_id']

            # Create dictionary of Chips and the links they are connected to
            self.chip_links = {}
    
            for link, ID in enumerate(chip_IDs):
                if ID not in self.chip_links:
                    self.chip_links[ID] = [link]
                else:
                    self.chip_links[ID].append(link)
            print('Chip links: ' + str(self.chip_links))

            # Sanity check
            #link_number = 3
            #for link, chipID in enumerate(self.chip_links):
            #    if link_number in self.chip_links[chipID]:
            #        print(link, chipID)

            # Get the number of chips
            self.num_of_chips = len(self.chip_links)

            # Create group to save all data and histograms to the HDF file
            h5_file.create_group(h5_file.root, 'interpreted', 'Interpreted Data')

            self.logger.info('Interpret raw data...')
            # Interpret the raw data (2x 32 bit to 1x 48 bit)
            hit_data = analysis.interpret_raw_data(raw_data, op_mode, vco, self.chip_links, meta_data, progress = progress)
            #raw_data = None

            for chip in range(self.num_of_chips):
                # get chipID of current chip
                chipID = str([ID for number, ID in enumerate(self.chip_links) if chip == number])[3:-2]
                print(chip, chipID)

                # create group for current chip
                h5_file.create_group(h5_file.root.interpreted, name = chipID)

                # get group for current chip
                chip_group = h5_file.root.interpreted._f_get_child(chipID)

                # Select only data which is hit data
                hit_data_chip = hit_data[chip][hit_data[chip]['data_header'] == 1]
                h5_file.create_table(chip_group, 'hit_data', hit_data_chip, filters=tb.Filters(complib='zlib', complevel=5))
                pix_occ       = np.bincount(hit_data_chip['x'] * 256 + hit_data_chip['y'], minlength=256 * 256).astype(np.uint32)
                hist_occ      = np.reshape(pix_occ, (256, 256)).T
                h5_file.create_carray(chip_group, name='HistOcc', obj=hist_occ)
                param_range   = np.unique(meta_data['scan_param_id'])
                #meta_data     = None
                pix_occ       = None
                hist_occ      = None

                # Create histograms for number of detected hits for individual thresholds
                scurve   = analysis.scurve_hist(hit_data_chip, param_range)
                #hit_data = None

                # Read needed configuration parameters
                n_injections     = [int(item[1]) for item in run_config if item[0] == b'n_injections'][0]
                Vthreshold_start = [int(item[1]) for item in run_config if item[0] == b'Vthreshold_start'][0]
                Vthreshold_stop  = [int(item[1]) for item in run_config if item[0] == b'Vthreshold_stop'][0]

                # Fit S-Curves to the histograms for all pixels
                param_range = list(range(Vthreshold_start, Vthreshold_stop + 1))
                thr2D, sig2D, chi2ndf2D = analysis.fit_scurves_multithread(scurve, scan_param_range=param_range, n_injections=n_injections, invert_x=False, progress = progress)

                h5_file.create_carray(chip_group, name='HistSCurve', obj=scurve)
                h5_file.create_carray(chip_group, name='Chi2Map', obj=chi2ndf2D.T)
                h5_file.create_carray(chip_group, name='ThresholdMap', obj=thr2D.T)
                h5_file.create_carray(chip_group, name='NoiseMap', obj=sig2D.T)

    def plot(self, status = None, plot_queue = None, **kwargs):
        '''
            Plot data and histograms of the scan
            If there is a status queue information about the status of the scan are put into it
        '''

        h5_filename = self.output_filename + '.h5'

        self.logger.info('Starting plotting...')
        if status != None:
            status.put("Create Plots")
        with tb.open_file(h5_filename, 'r+') as h5_file:

            with plotting.Plotting(h5_filename) as p:

                # Read needed configuration parameters
                Vthreshold_start = int(p.run_config[b'Vthreshold_start'])
                Vthreshold_stop = int(p.run_config[b'Vthreshold_stop'])
                n_injections = int(p.run_config[b'n_injections'])

                # Plot a page with all parameters
                p.plot_parameter_page()

                # Plot the equalisation bits histograms
                thr_matrix  = h5_file.root.configuration.thr_matrix[:]
                p.plot_distribution(thr_matrix, plot_range=np.arange(-0.5, 16.5, 1), title='Pixel threshold distribution', x_axis_title='Pixel threshold', y_axis_title='# of hits', suffix='pixel_threshold_distribution', plot_queue=plot_queue)
                        
                for chip in range(self.num_of_chips):
                    # get chipID of current chip
                    chipID = str([ID for number, ID in enumerate(self.chip_links) if chip == number])[3:-2]
                    print(chipID)

                    # get group of current chip in H5 file
                    chip_group = h5_file.root.interpreted._f_get_child(chipID)
                    
                    mask = h5_file.root.configuration.mask_matrix[:].T

                    # Plot the occupancy matrix
                    occ_masked  = np.ma.masked_array(chip_group.HistOcc[:], mask)
                    p.plot_occupancy(occ_masked, title='Integrated Occupancy, chip %s' %chipID, z_max='median', suffix='occupancy', plot_queue=plot_queue)
                        
                    # Plot the S-Curve histogram, put title for plot
                    scurve_hist = chip_group.HistSCurve[:].T
                    max_occ     = n_injections * 5
                    p.plot_scurves(scurve_hist, list(range(Vthreshold_start, Vthreshold_stop)), chipID, scan_parameter_name="Vthreshold", max_occ=max_occ, plot_queue=plot_queue)

                    # Do not plot pixels with converged  S-Curve fits 
                    chi2_sel        = chip_group.Chi2Map[:] > 0.  # Mask not converged fits (chi2 = 0)
                    mask[~chi2_sel] = True

                    # Plot the threshold distribution based on the S-Curve fits
                    hist = np.ma.masked_array(chip_group.ThresholdMap[:], mask)
                    p.plot_distribution(hist, plot_range=np.arange(Vthreshold_start-0.5, Vthreshold_stop-0.5, 1), x_axis_title='Vthreshold', title='Threshold distribution %s' %chipID, suffix='threshold_distribution', plot_queue=plot_queue)

                    # Plot the occupancy
                    p.plot_occupancy(hist, z_label='Threshold', title='Threshold, chip %s' %chipID, show_sum=False, suffix='threshold_map', z_min=Vthreshold_start, z_max=Vthreshold_stop, plot_queue=plot_queue)

                    # Plot the noise map
                    hist = np.ma.masked_array(chip_group.NoiseMap[:], mask)
                    p.plot_distribution(hist, plot_range=np.arange(0.1, 20, 0.1), title='Noise distribution, chip %s' %chipID, suffix='noise_distribution', plot_queue=plot_queue)
                    p.plot_occupancy(hist, z_label='Noise', title='Noise, chip %s' %chipID, show_sum=False, suffix='noise_map', z_min=0.1, z_max=20.0, plot_queue=plot_queue)


if __name__ == "__main__":
    scan = ThresholdScan()
    scan.start(**local_configuration)
    scan.analyze()
    scan.plot()
