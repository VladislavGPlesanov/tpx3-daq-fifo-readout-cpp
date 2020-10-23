import readline
import sys
import os
from multiprocessing import Process
from shutil import copy
from tpx3.scans.ToT_calib import ToTCalib
from tpx3.scans.scan_threshold import ThresholdScan
from tpx3.scans.scan_testpulse import TestpulseScan
from tpx3.scans.PixelDAC_opt import PixelDAC_opt
from tpx3.scans.take_data import DataTake
from tpx3.scans.Threshold_calib import ThresholdCalib
from UI.tpx3_logger import TPX3_datalogger, file_logger  #TODO:check if already opened instance by GUI

#In this part all callable function names should be in the list functions
functions = ['ToT', 'ToT_Calibration', 'tot_Calibration', 'tot', 
                'Threshold_Scan', 'THL_Scan', 'THL', 'threshold_scan', 'thl_scan', 'thl', 
                'Threshold_Calibration', 'THL_Calib', 'threshold_calibration', 'thl_calib',
                'Pixel_DAC_Optimisation', 'Pixel_DAC', 'PDAC', 'pixel_dac_optimisation', 'pixel_dac', 'pdac', 
                'Testpulse_Scan', 'TP_Scan', 'Tp_Scan' 'TP', 'testpulse_scan', 'tp_scan' 'tp', 
                'Run_Datataking', 'Run', 'Datataking', 'R', 'run_datataking', 'run', 'datataking', 'r',
                'Set_DAC', 'set_dac',
                'Load_Equalisation', 'Load_Equal', 'LEQ','load_equalisation', 'load_equal', 'leq',
                'Save_Equalisation', 'Save_Equal', 'SEQ','save_equalisation', 'save_equal', 'seq',
                'Save_Backup', 'Backup','save_backup', 'backup',
                'GUI',
                'Set_Polarity', 'Set_Pol', 'Polarity', 'Pol','set_polarity', 'set_pol', 'polarity','pol',
                'Set_Mask', 'Mask', 'set_mask', 'mask', 
                'Set_operation_mode', 'Set_Op_mode', 'Op_mode', 'set_operation_mode', 'set_Op_mode', 'op_mode',
                'Set_Fast_Io', 'Fast_Io', 'set_fast_io', 'fast_io', 'Fast_Io_en', 'fast_io_en',
                'Expert', 'expert',
                'Help', 'help', 'h', '-h',
                'End', 'end', 'Quit', 'quit', 'q', 'Q', 'Exit', 'exit']
help_functions = ['ToT_Calibration', 'Threshold_Scan', 'Threshold_Calibration', 'Pixel_DAC_Optimisation', 'Testpulse_Scan', 'Run_Datataking', 'Set_DAC','Load_Equalisation', 'Save_Equalisation', 'Set_Polarity', 'Set_operation_mode', 'Set_Fast_Io', 'Save_Backup',  'GUI', 'Help', 'Quit']

def completer(text, state):
    options = [function for function in functions if function.startswith(text)]
    try:
        return options[state]
    except IndexError:
        return None

class TPX3_multiprocess_start(object):
    def process_call(function, **kwargs):
        
        def startup_func(function, **kwargs):
            try:  
                call_func = (function+'()')
                scan = eval(call_func)
                scan.start(**kwargs)
                scan.analyze()
                scan.plot()
            except KeyboardInterrupt:
                sys.exit(1)
            except ValueError as e:
                print(e)
            except NotImplementedError:
                pass

        p = Process(target=startup_func, args=(function,), kwargs=kwargs)
        p.start()
        p.join()

class TPX3_CLI_funktion_call(object):#TODO: change to function_call

    TPX3_multiprocess_start = TPX3_multiprocess_start()

    def ToT_Calibration(object, VTP_fine_start = None, VTP_fine_stop = None, mask_step = None):
        if VTP_fine_start == None:
            print('> Please enter the VTP_fine_start value (0-511):')
            VTP_fine_start = int(input('>> '))
            print('> Please enter the VTP_fine_stop value (0-511):')
            VTP_fine_stop = int(input('>> '))
            print('> Please enter the number of steps(4, 16, 64, 256):')
            mask_step = int(input('>> '))
            
        print ('ToT calibration with VTP_fine_start =', VTP_fine_start, 'VTP_fine_stop =',VTP_fine_stop, 'mask_step =', mask_step)
        TPX3_multiprocess_start.process_call(function = 'ToTCalib', VTP_fine_start = VTP_fine_start, VTP_fine_stop = VTP_fine_stop, mask_step = mask_step, maskfile = TPX3_datalogger.read_value(name = 'Equalisation_path'))

    def Threshold_Scan(object, Vthreshold_start = None, Vthreshold_stop = None, n_injections = None, mask_step = None):
        if Vthreshold_start == None:
            print('> Please enter the Vthreshold_start value (0-2911):')
            Vthreshold_start = int(input('>> '))
            print('> Please enter the Vthreshold_stop value (0-2911):')
            Vthreshold_stop = int(input('>> '))
            print('> Please enter the number of injections (1-65535):')
            n_injections = int(input('>> '))
            print('> Please enter the number of steps(4, 16, 64, 256):')
            mask_step = int(input('>> '))
            
        print ('Threshold scan with Vthreshold_start =', Vthreshold_start, 'Vthreshold_stop =', Vthreshold_stop, 'Number of injections = ', n_injections, 'mask_step = ', mask_step)
        TPX3_multiprocess_start.process_call(function = 'ThresholdScan', Vthreshold_start = Vthreshold_start, Vthreshold_stop = Vthreshold_stop, n_injections = n_injections, mask_step = mask_step, maskfile = TPX3_datalogger.read_value(name = 'Equalisation_path'))

    def Threshold_Calib(object, Vthreshold_start = None, Vthreshold_stop = None, n_injections = None, mask_step = None, n_pulse_heights = None):
        if Vthreshold_start == None:
            print('> Please enter the Vthreshold_start value (0-2911):')
            Vthreshold_start = int(input('>> '))
            print('> Please enter the Vthreshold_stop value (0-2911):')
            Vthreshold_stop = int(input('>> '))
            print('> Please enter the number of injections (1-65535):')
            n_injections = int(input('>> '))
            print('> Please enter the number of steps(4, 16, 64, 256):')
            mask_step = int(input('>> '))
            print('> Please enter the number of pulse height steps(2-100):')
            n_pulse_heights = int(input('>> '))
            
        print ('Threshold scan with Vthreshold_start =', Vthreshold_start, 'Vthreshold_stop =', Vthreshold_stop, 'Number of injections = ', n_injections, 'mask_step = ', mask_step, 'Number of pulse heights = ', n_pulse_heights)
        TPX3_multiprocess_start.process_call(function = 'ThresholdCalib', Vthreshold_start = Vthreshold_start, Vthreshold_stop = Vthreshold_stop, n_injections = n_injections, mask_step = mask_step, n_pulse_heights = n_pulse_heights, maskfile = TPX3_datalogger.read_value(name = 'Equalisation_path'))


    def Testpulse_Scan(object, VTP_fine_start = None, VTP_fine_stop = None, n_injections = None, mask_step = None):
        if VTP_fine_start == None:
            print('> Please enter the VTP_fine_start value (0-511):')
            VTP_fine_start = int(input('>> '))
            print('> Please enter the VTP_fine_stop value (0-511):')
            VTP_fine_stop = int(input('>> '))
            print('> Please enter the number of injections (1-65535):')
            n_injections = int(input('>> '))
            print('> Please enter the number of steps(4, 16, 64, 256):')
            mask_step = int(input('>> '))
            
        print ('Testpulse scan with VTP_fine_start =', VTP_fine_start, 'VTP_fine_stop =',VTP_fine_stop, 'Number of injections = ', n_injections, 'mask_step =', mask_step)
        TPX3_multiprocess_start.process_call(function = 'TestpulseScan', VTP_fine_start = VTP_fine_start, VTP_fine_stop = VTP_fine_stop, n_injections = n_injections, mask_step = mask_step, maskfile = TPX3_datalogger.read_value(name = 'Equalisation_path'))

    def Pixel_DAC_Optimisation(object, Vthreshold_start = None, Vthreshold_stop = None, n_injections = None, mask_step = None):
        if Vthreshold_start == None:
            print('> Please enter the Vthreshold_start value (0-2911):')
            Vthreshold_start = int(input('>> '))
            print('> Please enter the Vthreshold_stop value (0-2911):')
            Vthreshold_stop = int(input('>> '))
            print('> Please enter the number of injections (1-65535):')
            n_injections = int(input('>> '))
            print('> Please enter the number of steps (4, 16, 64, 256):')
            mask_step = int(input('>> '))
        print ('Pixel DAC optimisation with Vthreshold_start =', Vthreshold_start, 'Vthreshold_stop =', Vthreshold_stop, 'Number of injections = ', n_injections, 'mask_step =', mask_step)
        TPX3_multiprocess_start.process_call(function = 'PixelDAC_opt', iteration = 0, Vthreshold_start = Vthreshold_start, Vthreshold_stop = Vthreshold_stop, n_injections = n_injections, mask_step = mask_step)

    def Set_DAC(object, DAC_Name = None, DAC_value = None):
        if DAC_Name == None:
            print('> Please enter the DAC-name from:\n    Ibias_Preamp_ON (0-255)\n    VPreamp_NCAS (0-255)\n    Ibias_Ikrum (0-255)\n    Vfbk (0-255)\n    Vthreshold_fine (0-511)\n    Vthreshold_coarse (0-15)\n    Ibias_DiscS1_ON (0-255)\n    Ibias_DiscS2_ON (0-255)\n    Ibias_PixelDAC (0-255)\n    Ibias_TPbufferIn (0-255)\n    Ibias_TPbufferOut (0-255)\n    VTP_coarse (0-255)\n    VTP_fine (0-511)\n    Ibias_CP_PLL (0-255)\n    PLL_Vcntrl (0-255)')
            DAC_Name = input('>> ')
            if DAC_Name in {'Ibias_Preamp_ON', 'VPreamp_NCAS', 'Ibias_Ikrum', 'Vfbk', 'Ibias_DiscS1_ON', 'Ibias_DiscS2_ON', 'Ibias_PixelDAC', 'Ibias_TPbufferIn', 'Ibias_TPbufferOut', 'VTP_coarse', 'Ibias_CP_PLL', 'PLL_Vcntrl'}:
                print('> Please enter the DAC value (0-255):')
                DAC_value = int(input('>> '))
            elif DAC_Name in {'Vthreshold_coarse'}:
                print('> Please enter the DAC value ( 0-15):')
                DAC_value = int(input('>> '))
            elif DAC_Name in {'Vthreshold_fine', 'VTP_fine'}:
                print('> Please enter the DAC value (0-511):')
                DAC_value = int(input('>> '))

        if DAC_Name in {'Ibias_Preamp_ON', 'VPreamp_NCAS', 'Ibias_Ikrum', 'Vfbk', 'Ibias_DiscS1_ON', 'Ibias_DiscS2_ON', 'Ibias_PixelDAC', 'Ibias_TPbufferIn', 'Ibias_TPbufferOut', 'VTP_coarse', 'Ibias_CP_PLL', 'PLL_Vcntrl'}:
            if DAC_value >= 0 and DAC_value <= 255:
                TPX3_datalogger.write_value(name = DAC_Name, value = DAC_value)
                TPX3_datalogger.write_to_yaml(name = DAC_Name)
                print('> Set ' + DAC_Name + ' to value ' + str(DAC_value) + '.')
            else:
                print('> Value ' + str(DAC_value) + ' is not in range (0-255)')
        elif DAC_Name in {'Vthreshold_coarse'}:
            if DAC_value >= 0 and DAC_value <= 15:
                TPX3_datalogger.write_value(name = DAC_Name, value = DAC_value)
                TPX3_datalogger.write_to_yaml(name = DAC_Name)
                print('> Set ' + DAC_Name + ' to value ' + str(DAC_value) + '.')
            else:
                print('> Value ' + str(DAC_value) + ' is not in range (0-15)')
        elif DAC_Name in {'Vthreshold_fine', 'VTP_fine'}:
            if DAC_value >= 0 and DAC_value <= 511:
                TPX3_datalogger.write_value(name = DAC_Name, value = DAC_value)
                TPX3_datalogger.write_to_yaml(name = DAC_Name)
                print('> Set ' + DAC_Name + ' to value ' + str(DAC_value) + '.')
            else:
                print('> Value ' + str(DAC_value) + ' is not in range (0-511)')
        else:
            print('Unknown DAC-name')

    def Load_Equalisation(object, equal_path = None):
        user_path = '~'
        user_path = os.path.expanduser(user_path)
        user_path = os.path.join(user_path, 'Timepix3')
        user_path = os.path.join(user_path, 'scans')
        user_path = os.path.join(user_path, 'hdf')
        
        if equal_path == None:
            print('> Please enter the name of the equalisation you like to load:')
            equal_path = input('>> ')
        try:
            #look if path exists
            full_path = user_path + os.sep + equal_path
            if os.path.isfile(full_path) == True:
                TPX3_datalogger.write_value(name = 'Equalisation_path', value = full_path)
        except:
            print('Path does not exist')

    def Save_Equalisation(object, file_name = None):
        user_path = '~'
        user_path = os.path.expanduser(user_path)
        user_path = os.path.join(user_path, 'Timepix3')
        user_path = os.path.join(user_path, 'scans')
        user_path = os.path.join(user_path, 'hdf')
        
        if file_name == None:
            print('> Please enter the path of the name you like to save the equalisation under:')
            file_name = input('>> ')
        try:
            #look if path exists
            full_path = user_path + os.sep + file_name + '.h5'
            if os.path.isfile(full_path) == True:
                print('File already exists')
            else:
                current_equal = TPX3_datalogger.read_value(name = 'Equalisation_path')
                copy(current_equal, full_path)
        except:
            print('Could not write file')

    def Save_Backup(object, file_name = None):
        user_path = '~'
        user_path = os.path.expanduser(user_path)
        user_path = os.path.join(user_path, 'Timepix3')
        user_path = os.path.join(user_path, 'backups')
        
        if file_name == None:
            print('> Please enter the path of the equalisation  you like to save the backup under:')
            file_name = input('>> ')
        try:
            #look if path exists
            full_path = user_path + os.sep + file_name + '.TPX3'
            if os.path.isfile(full_path) == True:
                print('File already exists')
            else:
                file.logger.write_backup(file = full_path)
        except:
            print('Could not write file')

    def Set_Polarity(object, polarity = None):
        if polarity == None:
            print('> Please enter the polarity (0 for positive or 1 for negative):')
            polarity = int(input('>> '))
        if polarity == 1 or polarity == 0:
            TPX3_datalogger.write_value(name = 'Polarity', value = polarity)
            TPX3_datalogger.write_to_yaml(name = 'Polarity')
        else:
            print('Unknown polarity')
    
    def Set_operation_mode(object, Op_mode = None):
        if Op_mode == None:
            print('> Please enter the operation mode (0 for ToT and TOA, 1 for only TOA, 2 for Event Count & Integral ToT):')
            Op_mode = int(input('>> '))
        if Op_mode >= 0 and Op_mode <= 2:
            TPX3_datalogger.write_value(name = 'Op_mode', value = Op_mode)
            TPX3_datalogger.write_to_yaml(name = 'Op_mode')
        else:
            print('Unknown operation mode')

    def Set_Fast_Io(object, Fast_Io_en = None):
        if Fast_Io_en == None:
            print('> Please enter the fast IO enable (0 for off or 1 for on):')
            Fast_Io_en = int(input('>> '))
        if Fast_Io_en == 1 or Fast_Io_en == 0:
            TPX3_datalogger.write_value(name = 'Fast_Io_en', value = Fast_Io_en)
            TPX3_datalogger.write_to_yaml(name = 'Fast_Io_en')
        else:
            print('Unknown polarity')

    def Run_Datataking(object, scan_timeout = None):
        if scan_timeout == None:
            print('> Please enter the required run time in seconds (choose 0 for an infinite run):')
            scan_timeout = int(input('>> '))
        
        if scan_timeout == 0:
            print('Infinite data taking run started! You can close the run with "ctrl. c"')
        else:
            print('{} s long data taking run started!'.format(scan_timeout))
            
        TPX3_multiprocess_start.process_call(function = 'DataTake', scan_timeout = scan_timeout, maskfile = TPX3_datalogger.read_value(name = 'Equalisation_path'))


class TPX3_CLI_TOP(object):
    def __init__(self, ext_input_list = None):
        readline.set_completer(completer)
        readline.parse_and_bind("tab: complete")
        funktion_call = TPX3_CLI_funktion_call()
        expertmode = False
        print ('\n Welcome to the Timepix3 control Software\n')

        # Here the main part of the cli starts. Every usercomand needs to be processed here.
        while 1:

            #if
            if ext_input_list == None:
                cmd_input = input('> ')
                #Catch if no input given
                if cmd_input == '':
                    print ('Something enter you must!')
                else:
                    inputlist = cmd_input.split()
            #Input is given
            else:
                inputlist = ext_input_list
                cmd_input = ' '.join(ext_input_list)
                ext_input_list = ['Quit']# To exit the while loop
            if inputlist:
                #Help
                if inputlist[0] in {'Help', 'help', 'h', '-h'}:
                    print('If you need detailed help on a function type [functionname -h].\n Possible options are:')
                    for function in help_functions:
                        print (function)

                #ToT_Calibration
                elif inputlist[0] in {'ToT', 'ToT_Calibration', 'tot_Calibration', 'tot'}:
                    if len(inputlist) == 1:
                        print('ToT_Calibration')
                        try:
                            funktion_call.ToT_Calibration()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the ToT calibration. As arguments you can give the start testpulse value (0-511), the stop testpulse value (0-511) and the number of steps (4, 16, 64, 256).')
                        elif len(inputlist) < 4:
                            print ('Incomplete set of parameters:')
                            try:
                                funktion_call.ToT_Calibration()
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) == 4:
                            try:
                                funktion_call.ToT_Calibration(VTP_fine_start = int(inputlist[1]), VTP_fine_stop = int(inputlist[2]), mask_step = int(inputlist[3]))
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 4:
                            print ('To many parameters! The given function takes only three parameters:\n start testpulse value (0-511),\n stop testpulse value (0-511),\n number of steps (4, 16, 64, 256).')

                #Threshold_Scan
                elif inputlist[0] in {'Threshold_Scan', 'THL_Scan', 'THL', 'threshold_scan', 'thl_scan', 'thl'}:
                    if len(inputlist) == 1:
                        print('Threshold_Scan')
                        try:
                            funktion_call.Threshold_Scan()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the Threshold scan. As arguments you can give the start threshold value (0-2911), the stop threshold value (0-2911), the number of testpulse injections (1-65535) and the number of steps (4, 16, 64, 256).')
                        elif len(inputlist) < 5:
                            print ('Incomplete set of parameters:')
                            try:
                                funktion_call.Threshold_Scan()
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) == 5:
                            try:
                                funktion_call.Threshold_Scan(Vthreshold_start = int(inputlist[1]), Vthreshold_stop = int(inputlist[2]), n_injections = int(inputlist[3]), mask_step = int(inputlist[4]))
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 5:
                            print ('To many parameters! The given function takes only four parameters:\n start testpulse value (0-2911),\n stop testpulse value (0-2911),\n number of injections (1-65535),\n number of steps (4, 16, 64, 256).')
               
                #Threshold_Calib
                elif inputlist[0] in {'Threshold_Calibration', 'THL_Calib', 'threshold_calibration', 'thl_calib',}:
                    if len(inputlist) == 1:
                        print('Threshold_Calibration')
                        try:
                            funktion_call.Threshold_Calib()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the Threshold calibration. As arguments you can give the start threshold value (0-2911), the stop threshold value (0-2911), the number of testpulse injections (1-65535), the number of steps (4, 16, 64, 256) and the number of pulse height steps (2-100).')
                        elif len(inputlist) < 6:
                            print ('Incomplete set of parameters:')
                            try:
                                funktion_call.Threshold_Scan()
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) == 6:
                            try:
                                funktion_call.Threshold_Scan(Vthreshold_start = int(inputlist[1]), Vthreshold_stop = int(inputlist[2]), n_injections = int(inputlist[3]), mask_step = int(inputlist[4]), n_pulse_height = int(inputlist[5]))
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 6:
                            print ('To many parameters! The given function takes only four parameters:\n start testpulse value (0-2911),\n stop testpulse value (0-2911),\n number of injections (1-65535),\n number of steps (4, 16, 64, 256),\n number of pulse height steps (2-100).')

                #Testpulse_Scan
                elif inputlist[0] in {'Testpulse_Scan', 'TP_Scan', 'Tp_Scan' 'TP', 'testpulse_scan', 'tp_scan' 'tp'}:
                    if len(inputlist) == 1:
                        print('Testpulse_Scan')
                        try:
                            funktion_call.Testpulse_Scan()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the Testpulse Scan. As arguments you can give the the start testpulse value (0-511), the stop testpulse value (0-511), the number of testpulse injections (1-65535) and the number of steps (4, 16, 64, 256).')
                        elif len(inputlist) < 5:
                            print ('Incomplete set of parameters:')
                            try:
                                funktion_call.Testpulse_Scan()
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) == 5:
                            try:
                                funktion_call.Testpulse_Scan(VTP_fine_start = int(inputlist[1]), VTP_fine_stop = int(inputlist[2]), n_injections = int(inputlist[3]), mask_step = int(inputlist[4]))
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 5:
                            print ('To many parameters! The given function takes only four parameters:\n start testpulse value (0-511),\n stop testpulse value (0-511),\n number of injections (1-65535),\n number of steps (4, 16, 64, 256).')

                #Pixel_DAC_Optimisation
                elif inputlist[0] in {'Pixel_DAC_Optimisation', 'Pixel_DAC', 'PDAC', 'pixel_dac_optimisation', 'pixel_dac', 'pdac'}:
                    if len(inputlist) == 1:
                        print('Pixel_DAC_Optimisation')
                        try:
                            funktion_call.Pixel_DAC_Optimisation()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the Pixel DAC Optimisation. As arguments you can give the start threshold value (0-2911), the stop threshold value (0-2911), the number of testpulse injections (1-65535) and the number of steps (4, 16, 64, 256).')
                        elif len(inputlist) < 5:
                            print ('Incomplete set of parameters:')
                            try:
                                funktion_call.Pixel_DAC_Optimisation()
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) == 5:
                            try:
                                funktion_call.Pixel_DAC_Optimisation(Vthreshold_start = int(inputlist[1]), Vthreshold_stop = int(inputlist[2]), n_injections = int(inputlist[3]), mask_step = int(inputlist[4]))
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 5:
                            print ('To many parameters! The given function takes only four parameters:\n start testpulse value (0-2911),\n stop testpulse value (0-2911),\n number of injections (1-65535),\n number of steps (4, 16, 64, 256).')

                #Set_DAC
                elif inputlist[0] in {'Set_DAC', 'set_dac'}:
                    if len(inputlist) == 1:
                        print('Set_DAC')
                        try:
                            funktion_call.Set_DAC()
                        except:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the Set DAC function. As arguments you can give the DAC-name/DAC-number  and the new value.\n The following DACs are aviable:\n     1.) Ibias_Preamp_ON (0-255)\n     2.) VPreamp_NCAS (0-255)\n     3.) Ibias_Ikrum (0-255)\n     4.) Vfbk (0-255)\n     5.) Vthreshold_fine (0-511)\n     6.) Vthreshold_coarse (0-15)\n     7.) Ibias_DiscS1_ON (0-255)\n     8.) Ibias_DiscS2_ON (0-255)\n     9.) Ibias_PixelDAC (0-255)\n    10.) Ibias_TPbufferIn (0-255)\n    11.) Ibias_TPbufferOut (0-255)\n    12.) VTP_coarse (0-255)\n    13.) VTP_fine (0-511)\n    14.) Ibias_CP_PLL (0-255)\n    15.) PLL_Vcntrl (0-255)')                
                        elif len(inputlist) < 3:
                            print ('Incomplete set of parameters:')
                            try:
                                funktion_call.Set_DAC()
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) == 3:
                            if inputlist[1] in {'1', 'Ibias_Preamp_ON'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_Preamp_ON', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'2', 'VPreamp_NCAS'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'VPreamp_NCAS', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'3', 'Ibias_Ikrum'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_Ikrum', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'4', 'Vfbk'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Vfbk', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'5', 'Vthreshold_fine'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Vthreshold_fine', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'6', 'Vthreshold_coarse'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Vthreshold_coarse', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'7', 'Ibias_DiscS1_ON'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_DiscS1_ON', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'8', 'Ibias_DiscS2_ON'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_DiscS2_ON', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'9', 'Ibias_PixelDAC'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_PixelDAC', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'10', 'Ibias_TPbufferIn'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_TPbufferIn', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'11', 'Ibias_TPbufferOut'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_TPbufferOut', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'12', 'VTP_coarse'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'VTP_coarse', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'13', 'VTP_fine'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'VTP_fine', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'14', 'Ibias_CP_PLL'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'Ibias_CP_PLL', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'15', 'PLL_Vcntrl'}:
                                try:
                                    funktion_call.Set_DAC(DAC_Name = 'PLL_Vcntrl', DAC_value = int(inputlist[2]))
                                except KeyboardInterrupt:
                                    print('User quit')
                            else:
                                print ('Unknown DAC-name')
                        elif len(inputlist) > 3:
                            print ('To many parameters! The given function takes only two parameters:\n The DAC-name and its value.')

                #Data taking
                elif inputlist[0] in {'Run_Datataking', 'Run', 'Datataking', 'R', 'run_datataking', 'run', 'datataking', 'r'}:
                    if len(inputlist) == 1:
                        print('Run_Datataking')
                        try:
                            funktion_call.Run_Datataking()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the datataking function. As argument you can give the scan timeout (in seconds, if 0 is entered the datataking will run infinitely')
                        elif len(inputlist) == 2:
                            try:
                                funktion_call.Run_Datataking(scan_timeout = int(inputlist[1]))
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 2:
                            print ('To many parameters! The given function takes only one parameters:\n scan timeout (in seconds).')

                #Load equalisation
                elif inputlist[0] in {'Load_Equalisation', 'Load_Equal', 'LEQ','load_equalisation', 'load_equal', 'leq'}:
                    if len(inputlist) == 1:
                        print('Load_Equalisation')
                        try:
                            funktion_call.Load_Equalisation()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the load equalisation function. As argument you can give the path of the equalisation you like to load')
                        elif len(inputlist) == 2:
                            try:
                                funktion_call.Load_Equalisation(equal_path = inputlist[1])
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 2:
                            print ('To many parameters! The given function takes only one parameters:\n equalisation path.')

                #Save equalisation
                elif inputlist[0] in {'Save_Equalisation', 'Save_Equal', 'SEQ','save_equalisation', 'save_equal', 'seq'}:
                    if len(inputlist) == 1:
                        print('Save_Equalisation')
                        try:
                            funktion_call.Save_Equalisation()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the save equalisation function. As argument you can give the name of the equalisation file')
                        elif len(inputlist) == 2:
                            try:
                                funktion_call.Save_Equalisation(file_name = inputlist[1])
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 2:
                            print ('To many parameters! The given function takes only one parameters:\n equalisation file name.')

                #Save backup
                elif inputlist[0] in {'Save_Backup', 'Backup','save_backup', 'backup'}:
                    if len(inputlist) == 1:
                        print('Save_Backup')
                        try:
                            funktion_call.Save_Backup()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the save backup function. As argument you can give the name of the backup file')
                        elif len(inputlist) == 2:
                            try:
                                funktion_call.Save_Backup(file_name = inputlist[1])
                            except KeyboardInterrupt:
                                print('User quit')
                        elif len(inputlist) > 2:
                            print ('To many parameters! The given function takes only one parameters:\n backup file name.')


                #Set polarity
                elif inputlist[0] in {'Set_Polarity', 'Set_Pol', 'Polarity', 'Pol','set_polarity', 'set_pol', 'polarity','pol'}:
                    if len(inputlist) == 1:
                        print('Set_Polarity')
                        try:
                            funktion_call.Set_Polarity()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the set polarity function. As argument you can give the polarity as {negative, neg, -, 1} or {positive, pos, +, 0}')
                        elif len(inputlist) == 2:
                            if inputlist[1] in {'negative', 'neg', '-', '1'}:
                                try:
                                    funktion_call.Set_Polarity(polarity = 1)
                                except KeyboardInterrupt:
                                    print('User quit')
                            elif inputlist[1] in {'positive', 'pos', '+', '0'}:
                                try:
                                    funktion_call.Set_Polarity(polarity = 0)
                                except KeyboardInterrupt:
                                    print('User quit')
                            else:
                                print('Unknown polarity use {negative, neg, -, 1} or {positive, pos, +, 0}')
                        elif len(inputlist) > 2:
                            print ('To many parameters! The given function takes only one parameters:\n polarity.')

                #Set operation mode
                elif inputlist[0] in {'Set_operation_mode', 'Set_Op_mode', 'Op_mode', 'set_operation_mode', 'set_Op_mode', 'op_mode'}:
                    if len(inputlist) == 1:
                        print('Set_operation_mode')
                        try:
                            funktion_call.Set_operation_mode()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the Set operation mode function. As argument you can give the operation mode as 0 for ToT & ToA, 1 for only ToA or 2 for Event Count & Integral ToT')
                        elif len(inputlist) == 2:
                                try:
                                    funktion_call.Set_operation_mode(Op_mode = int(inputlist[1]))
                                except KeyboardInterrupt:
                                    print('User quit')
                        elif len(inputlist) > 2:
                            print ('To many parameters! The given function takes only one parameters:\n polarity.')
                
                #Set Fast Io mode
                elif inputlist[0] in {'Set_Fast_Io', 'Fast_Io', 'Fast_Io_en', 'set_fast_io', 'fast_io', 'fast_io_en'}:
                    if len(inputlist) == 1:
                        print('Set_Fast_Io')
                        try:
                            funktion_call.Set_Fast_Io()
                        except KeyboardInterrupt:
                            print('User quit')
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This is the Fast Io enable function. As argument you can give the enable as 0 (off) or 1 (on)')
                        elif len(inputlist) == 2:
                                try:
                                    funktion_call.Set_Fast_Io(Fast_Io_en = int(inputlist[1]))
                                except KeyboardInterrupt:
                                    print('User quit')
                        elif len(inputlist) > 2:
                            print ('To many parameters! The given function takes only one parameters:\n Fast Io enable.')

                #Start GUI
                elif inputlist[0] in {'GUI'}:
                    if len(inputlist) == 1:
                        #Start GUI
                        print('GUI started')
                        break
                    else:
                        if inputlist[1] in {'Help', 'help', 'h', '-h'}:
                            print('This will start the GUI')
                        elif len(inputlist) > 1:
                            print('GUI takes no parameters')

                #Set expert mode
                elif inputlist[0] in {'Expert', 'expert'}:
                    if expertmode == False:
                        expertmode = True
                    elif expertmode == True:
                        expertmode = False

                #Quit
                elif inputlist[0] in {'End', 'end', 'Quit', 'quit', 'q', 'Q', 'Exit', 'exit'}:
                    print('Goodbye and have a nice day.')
                    break

                #Unknown command
                else:
                    print ('Unknown command: ', cmd_input, 'Use a language I understand.')

if __name__ == "__main__":
    ext_input_list = sys.argv
    ext_input_list.pop(0)
    if ext_input_list == []:
        tpx3_cli = TPX3_CLI_TOP()
    else:
        tpx3_cli = TPX3_CLI_TOP(ext_input_list = ext_input_list)