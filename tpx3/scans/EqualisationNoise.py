#
# ------------------------------------------------------------
# Copyright (c) All rights reserved
# SiLab, Institute of Physics, University of Bonn
# ------------------------------------------------------------
#

'''
    This script performs an equalisation of pixels based on a threshold scan
    with noise.
'''
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from tqdm import tqdm
import numpy as np
import time
import tables as tb
import os
import math

from tpx3.scan_base import ScanBase
import tpx3.analysis as analysis
import tpx3.plotting as plotting
import tpx3.utils as utils

from tables.exceptions import NoSuchNodeError
from six.moves import range

local_configuration = {
    # Scan parameters
    'mask_step'        : 16,
    'Vthreshold_start' : 1000,
    'Vthreshold_stop'  : 1350
}


class EqualisationNoise(ScanBase):

    scan_id = "EqualisationNoise"
    wafer_number = 0
    y_position = 0
    x_position = 'A'

    def scan(self, Vthreshold_start = 1000, Vthreshold_stop = 1350, mask_step = 16, progress = None, status = None, **kwargs):
        '''
            Takes data for equalisation. Therefore a threshold scan is performed for all pixel thresholds at 0 and at 15.
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
        if mask_step not in {4, 16, 64, 256}:
            raise ValueError("Value {} for mask_step is not in the allowed range (4, 16, 64, 256)".format(mask_step))

        # Set general configuration registers of the Timepix3
        self.chip.write_general_config()

        # Write to the test pulse registers of the Timepix3
        # Write to period and phase tp registers
        # This is needed here to open the internal Timepix3 shutter
        data = self.chip.write_tp_period(1, 0)

        self.logger.info('Preparing injection masks...')
        if status != None:
            status.put("Preparing injection masks")

        # Create the masks for all steps for the scan at 0 and at 15
        mask_cmds = self.create_scan_masks(mask_step, pixel_threhsold = 0, progress = progress)
        mask_cmds2 = self.create_scan_masks(mask_step, pixel_threhsold = 15, progress = progress)

        # Scan with pixel threshold 0
        self.logger.info('Starting scan for THR = 0...')
        if status != None:
            status.put("Starting scan for THR = 0")
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
            # Set the threshold
            self.chip.set_dac("Vthreshold_coarse", int(threshold[0]))
            self.chip.set_dac("Vthreshold_fine", int(threshold[1]))

            with self.readout(scan_param_id=scan_param_id):
                for mask_step_cmd in mask_cmds:
                    # Write the pixel matrix for the current step plus the read_pixel_matrix_datadriven command
                    self.chip.write(mask_step_cmd)

                    # Open the shutter, take data and update the progress bar
                    with self.shutter():
                        time.sleep(0.01)
                        if progress == None:
                            # Update the progress bar
                            pbar.update(1)
                        else:
                            # Update the progress fraction and put it in the queue
                            step_counter += 1
                            fraction = step_counter / (len(mask_cmds) * len(thresholds))
                            progress.put(fraction)
                    self.chip.stop_readout()
                    time.sleep(0.001)
                self.chip.reset_sequential()
                time.sleep(0.001)
            scan_param_id += 1

        if progress == None:
            # Close the progress bar
            pbar.close()

        # Scan with pixel threshold 15
        self.logger.info('Starting scan for THR = 15...')
        if status != None:
            status.put("Starting scan for THR = 15")
        if status != None:
            status.put("iteration_symbol")

        if progress == None:
            # Initialize progress bar
            pbar = tqdm(total=len(mask_cmds2) * len(threshold))
        else:
            # Initialize counter for progress
            step_counter = 0

        scan_param_id = 0
        for threshold in thresholds:
            # Set the threshold
            self.chip.set_dac("Vthreshold_coarse", int(threshold[0]))
            self.chip.set_dac("Vthreshold_fine", int(threshold[1]))

            with self.readout(scan_param_id=scan_param_id + len(thresholds)):
                for mask_step_cmd in mask_cmds2:
                    # Only activate testpulses for columns with active pixels
                    self.chip.write(mask_step_cmd)

                    # Open the shutter, take data and update the progress bar
                    with self.shutter():
                        time.sleep(0.01)
                        if progress == None:
                            # Update the progress bar
                            pbar.update(1)
                        else:
                            # Update the progress fraction and put it in the queue
                            step_counter += 1
                            fraction = step_counter / (len(mask_cmds2) * len(thresholds))
                            progress.put(fraction)
                    self.chip.stop_readout()
                    time.sleep(0.001)
                self.chip.reset_sequential()
                time.sleep(0.001)
            scan_param_id += 1

        if progress == None:
            # Close the progress bar
            pbar.close()

        if status != None:
            status.put("iteration_finish_symbol")

        self.logger.info('Scan finished')

    def analyze(self, progress = None, status = None, result_path = None, **kwargs):
        '''
            Analyze the data of the equalisation and calculate the equalisation matrix
            If progress is None a tqdm progress bar is used else progress should be a Multiprocess Queue which stores the progress as fraction of 1
            If there is a status queue information about the status of the scan are put into it
        '''

        h5_filename = self.output_filename + '.h5'

        self.logger.info('Starting data analysis...')
        if status != None:
            status.put("Performing data analysis")

        # Open the HDF5 which contains all data of the equalisation
        with tb.open_file(h5_filename, 'r+') as h5_file:
            # Read raw data, meta data and configuration parameters
            meta_data = h5_file.root.meta_data[:]
            run_config = h5_file.root.configuration.run_config[:]
            general_config = h5_file.root.configuration.generalConfig[:]
            op_mode = [row[1] for row in general_config if row[0]==b'Op_mode'][0]
            vco = [row[1] for row in general_config if row[0]==b'Fast_Io_en'][0]

            self.logger.info('Interpret raw data...')

            # THR = 0
            param_range, index = np.unique(meta_data['scan_param_id'], return_index=True)
            meta_data_th0 = meta_data[meta_data['scan_param_id'] < len(param_range) // 2]
            param_range_th0 = np.unique(meta_data_th0['scan_param_id'])

            # THR = 15
            meta_data_th15 = meta_data[meta_data['scan_param_id'] >= len(param_range) // 2]
            param_range_th15 = np.unique(meta_data_th15['scan_param_id'])

            # shift indices so that they start with zero
            start = meta_data_th15['index_start'][0]
            meta_data_th15['index_start'] = meta_data_th15['index_start']-start
            meta_data_th15['index_stop'] = meta_data_th15['index_stop']-start

            self.logger.info('THR = 0')
            #THR = 0
            raw_data_thr0 = h5_file.root.raw_data[:meta_data_th0['index_stop'][-1]]
            hit_data_thr0 = analysis.interpret_raw_data(raw_data_thr0, op_mode, vco, meta_data_th0, progress = progress)
            raw_data_thr0 = None

            self.logger.info('THR = 15')
            #THR = 15
            raw_data_thr15 = h5_file.root.raw_data[meta_data_th0['index_stop'][-1]:]
            hit_data_thr15 = analysis.interpret_raw_data(raw_data_thr15, op_mode, vco, meta_data_th15, progress = progress)
            raw_data_thr15 = None

        # Read needed configuration parameters
        Vthreshold_start = [int(item[1]) for item in run_config if item[0] == b'Vthreshold_start'][0]
        Vthreshold_stop = [int(item[1]) for item in run_config if item[0] == b'Vthreshold_stop'][0]
        chip_wafer = [int(item[1]) for item in run_config if item[0] == b'chip_wafer'][0]
        chip_x = [item[1].decode() for item in run_config if item[0] == b'chip_x'][0]
        chip_y = [int(item[1]) for item in run_config if item[0] == b'chip_y'][0]

        # Select only data which is hit data
        hit_data_thr0 = hit_data_thr0[hit_data_thr0['data_header'] == 1]
        hit_data_thr15 = hit_data_thr15[hit_data_thr15['data_header'] == 1]

        # Divide the data into two parts - data for pixel threshold 0 and 15
        param_range = np.unique(meta_data['scan_param_id'])
        meta_data = None
        param_range_th0 = np.unique(hit_data_thr0['scan_param_id'])
        param_range_th15 = np.unique(hit_data_thr15['scan_param_id'])

        # Create histograms for number of detected hits for individual thresholds
        self.logger.info('Get the global threshold distributions for all pixels...')
        scurve_th0 = analysis.scurve_hist(hit_data_thr0, param_range_th0)
        hit_data_thr0 = None
        scurve_th15 = analysis.scurve_hist(hit_data_thr15, param_range_th15)
        hit_data_thr15 = None

        # Calculate the mean of the threshold distributions for all pixels
        self.logger.info('Calculate the mean of the global threshold distributions for all pixels...')
        vths_th0 = analysis.vths(scurve_th0, param_range_th0, Vthreshold_start)
        scurve_th0 = None
        vths_th15 = analysis.vths(scurve_th15, param_range_th15, Vthreshold_start)
        scurve_th15 = None

        # Get the treshold distributions for both scan
        self.logger.info('Get the cumulated global threshold distributions...')
        hist_th0 = analysis.vth_hist(vths_th0, Vthreshold_stop)
        hist_th15 = analysis.vth_hist(vths_th15, Vthreshold_stop)
        vths_th15 = None

        # Use the threshold histograms and one threshold distribution to calculate the equalisation
        self.logger.info('Calculate the equalisation matrix...')
        eq_matrix = analysis.eq_matrix(hist_th0, hist_th15, vths_th0, Vthreshold_start, Vthreshold_stop)

        # Don't mask any pixels in the mask file
        mask_matrix = np.zeros((256, 256), dtype=bool)
        mask_matrix[:, :] = 0

        # Write the equalisation matrix to a new HDF5 file
        self.save_thr_mask(eq_matrix, chip_wafer, chip_x ,chip_y)

        if result_path != None:
            result_path.put(self.thrfile)


if __name__ == "__main__":
    scan = EqualisationNoise()
    scan.start(**local_configuration)
    scan.analyze()