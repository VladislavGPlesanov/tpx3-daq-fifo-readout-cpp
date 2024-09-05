//
// C++ class for reading data from TPX3 chip located on FPGA FIFO (of basil SiTCP type)
//
// Uni-Bonn, Vlad
// 
// Description of python object structures and interplay are described in how-to-firmware.html
// (ask around for this file)
//
//
#define PY_SSIZE_T_CLEAN

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>
#include <Python.h>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <thread>
#include <string>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <bitset>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL readoutBuffer_ARRAY_API

#include <numpy/arrayobject.h>
#include <tpx3constants.h>

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC PyInit_Mod(void){
    Py_Initialize();
    import_array();
    if(PyErr_Occurred()){
        std::cerr<<"[ERROR] Failed to import numpy Python module(s)"<<std::endl;
    }
    assert(PyArray_API);
    return NULL;
}


// note to any unfortunate successor/developer:
//
// For many finctions Global Interpreter Lock (GIL) states are fixed at the start to the end
// of their run logic to ensure thread safe operation compliant that python operation requires.
// Whenever one deals with a python object one must consider its pointer reference count, which 
// is handled by Py_INCREF(obj), which increases reference count, and Py_DECREF(obj), 
// which decreases reference count.
// This number should not drop to zero.
// At the end of any routine/function/loop object's ref count should be decremented
// In some case incrementation of this counter is welcome to ensure the counter doesn't 
// drop to zero. (Keep track of it using Py_REF_CNT(obj))
// Be carefull...
//

///////// global variables ///////////////////////////

signed long int global_words_recorded = 0; 

//////////////////////////////////////////////////////

extern "C" {

    // used for timing
    chronoTime tick(){
        return std::chrono::high_resolution_clock::now();
    }


    void flush_debug(std::string msg){
        std::cout<<msg<<"\n"<<std::flush;
    };

//    void dumpVector(std::vector<int> vect){
//        int cnt = 0;
//        for(auto &el: vect){
//            std::cout<<"RX["<<cnt<<"]=["<<el<<"]"<<"\t"<<std::flush;
//            cnt++;
//        }
//    };

//    void printVectUint(std::vector<uint32_t> vect){    
//        
//        int vsize = vect.size();
//        int iter_max = 64;
//        if(vsize<iter_max){
//            iter_max = (int)vsize;
//        }
//
//        for(int i = 0; i < iter_max; i++){
//            std::cout<<vect.at(i)<<"|"<<std::flush;
//            if(i !=0 && i % 8 == 0){
//              std::cout<<"\n"<<std::flush;
//            }
//        }
//        int items_left = 0;
//        if(vsize >= iter_max){
//            items_left = vsize - iter_max;
//        }
//        std::cout<<" and "<<items_left<<" more\n"<<std::flush;
//    }

//    void printReverseVectUint(std::vector<uint32_t> vect){    
//
//        int vsize = vect.size();
//        int iter_max = 64;
//
//        if(vsize >= iter_max*2){
//
//            for(int i = vsize-iter_max; i<vsize; i++){
//                std::cout<<vect.at(i)<<"|"<<std::flush;
//                if(i !=vsize-iter_max && i % 8 == 0){
//                  std::cout<<"\n"<<std::flush;
//                }
//            }
//        }
//    }

    // NIU, but can be used for debuging pyobjects directly or as example
    //
    const char* printObjType(PyObject *obj){ 

        PyTypeObject *type = obj->ob_type;
        const char* type_name = type->tp_name;
        std::cout<<"[TYPE CHECK] ("<<type_name<<")\n"<<std::flush;
        return type_name;

    }

    // returns time duration in milliseconds that is left before certain time point
    float time_to_read(chronoTime start_time){
         float millisec =  std::chrono::duration<float, std::milli>(
                           std::chrono::high_resolution_clock::now() - start_time).count();
  
         return millisec;
    };

    // returns time difference in milliseconds for two fixed time points
    float timeDiff(chronoTime t0, chronoTime t1){
         float millisec = std::chrono::duration<float, std::milli>(t1 - t0).count();
         return millisec;
    };

 
    // retuns epoch time in seconds  
    double get_float_time(){

        auto now = std::chrono::system_clock::now(); //current time point
        auto duration = now.time_since_epoch(); // getting epoch time
        // converting epoch time to seconds
        double seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count(); 
        //converting epoch tme to microseconds
        double microseconds = std::chrono::duration_cast<std::chrono::microseconds>(duration).count() % 1000000;

        return seconds + microseconds * 1e-6; 

    };

    // updates timestamps for chunk timing
    // used in readoutToDeque()
    void upd_time(double &curr_time, double &last_time, double &timestamp){
        
        // get current time
        curr_time = get_float_time(); 
        //assign timestamp (strt time of this chunk) to previous chunk stop time 
        last_time = timestamp; 
        //assign stop time of current chunk to current time
        timestamp = curr_time; 

    };

    // returns sum of the error counters of all TPX3_RX
    //used in readoutToDeque
    long int getNerror(PyObject *list){

        // error counter    
        long int n_errors = 0;
        //obtaining length of receiver (RX) list of counters
        Py_ssize_t listlen = PyObject_Length(list);
        // iterating over counter list
        for(Py_ssize_t it = 0; it < listlen; ++it){
            PyObject *item = PyList_GetItem(list, it);
            if(item!=NULL){
                Py_INCREF(item);
                n_errors+=PyLong_AsLong(item);
                Py_DECREF(item);
            }else{
                flush_debug("can not access item in:");
                PyObject_Print(list,stdout,0);
                Py_DECREF(item);
            }
            Py_XDECREF(item);
            item=NULL;
        };
        // return sum of counters
        return n_errors;
    };

    // returns empty PyList
    // used in readFifoStatus()
    PyObject* emptyList(){
        PyObject *emList = PyList_New(8);
        for(int i=0; i<8; i++){
            PyObject *val = PyLong_FromLong(0);
            PyList_SetItem(emList, i, val);
        }
        return emList;
    }

    // puts data of a single readut iteration that is stored in temporary std::vector<uint32_t>
    // and puts it into a python PyArray for transfer to python side via python tuple
    // used in readoutToDeque()    
    PyObject* fillPyArray(std::vector<uint32_t> &data){// works perfectly fine
    
        // get umber of data words
        const std::size_t datasize = data.size();
        // initalize container size to store data in py array
        const npy_intp pyarrsize[] = {(npy_intp)(datasize)};
        
        // initialize empty PyArray with element size pyarrsize[] and of type numpy.uint32
        PyObject* pyarray = PyArray_SimpleNew(1,pyarrsize,NPY_UINT32);
        if(pyarray==NULL){
            flush_debug("ARRAY IS NULL!");
        }
    
        // allocate sufficient memory size for python array that is of size  
        // on the python memory stack
        uint32_t* thisData = reinterpret_cast<uint32_t*>(
                PyArray_DATA(reinterpret_cast<PyArrayObject*>(pyarray)));
    
        // copy c++ data into pyarray 
        //std::copy(data.begin(),data.end(),thisData); 
        // to try:
        memcpy(thisData, &data[0], datasize*sizeof(uint32_t)); // this is actually faster..
                                                               // but not sure why reads less...
        return pyarray;
    
    };
   
    // returns std::vector<uint32_t> of data words recorded from basil deque
    // used in readoutToDeque()
    std::vector<uint32_t> fillVectorFromPyData(PyArrayObject *pyarray){
    
        //get address of the first element of python numpy array in data member
        // of py array structure
        uint32_t* result_data = (uint32_t*)PyArray_DATA(pyarray);
        //adjust get the dimension of numpy array
        npy_intp* arr_dim = PyArray_DIMS(pyarray);
        //get size of each dimension
        npy_intp* arr_size = &arr_dim[0];
        //put data from numpy array into the vector
        std::vector<uint32_t> result_vector(result_data, result_data + *arr_size);
        Py_DECREF(pyarray); // Have to release pyarray reference here
                            // - otherwise 2Gb leak in 5 seconds! 
    
        return result_vector;
    
    };
    //////////////////////////////////////////////////////////////////////////
    ///---------- working numpy C-API functions --------------------------
    //////////////////////////////////////////////////////////////////////////
    
    // returns pyobject of SiTcp class of basil that is used to get data from the socket
    // can be used to manipulate parameters of 
    // SiTcp class (though its should not have been possible...)
    PyObject *get_interface(PyObject *fifo){
       //get pointer to the interface object defined in basil.dut class as "_intf" 
       PyObject *intf = PyObject_GetAttrString(fifo, "_intf");
       if(intf=NULL){
          flush_debug("[ERROR] can not get basil.SiTcp._intf object -> nullptr!");
       }
       Py_DECREF(fifo); //release reference of self.chip['FIFO'] object
    
       return intf;
    
    }
    
    // can be used to check polling interval if select.select function that polls data
    // within basil.SiTcp._tcp_readout loop
    void checkTcpReadoutInt(PyObject *intf){
    
        if(intf != NULL){  
          PyObject *tcp_readout_interval = PyObject_GetAttrString(intf, "_tcp_readout_interval");
          if(tcp_readout_interval!=NULL){
            flush_debug("Got instance of tcp readout interval");
            std::cout<<"basil tcp readout interval is ["
                    <<PyObject_Print(tcp_readout_interval,stdout,0)<<"] s"<<std::flush;
            Py_DECREF(tcp_readout_interval);
          }else{
            flush_debug("Can not access: tcp readout interval");
            Py_DECREF(tcp_readout_interval);
          }
    
        }else{
            flush_debug("Object \"intf\" is NULL!");
        }
    
    }
    
    //can be used to set the desired data polling interval 
    //wihtin basil.SiTcp._tcp_readout object
    void set_tcp_interval(PyObject *interface, float time){ 
    
        PyObject_SetAttrString(interface,
                               "_tcp_readout_interval", 
                               PyFloat_FromDouble((double)time));
    
    }
   
    // calls status/counter check functions in tpx3/fifo_readout.py
    //
    PyObject* readFifoStatus(PyObject *self, const char* option){       
    
    // test case when passing self in fifo_readout 
    // object on the python side should be: <tpx3.fifo_readout.FifoReadout>
    
        PyGILState_STATE gstate = PyGILState_Ensure(); 
        PyObject *status;
    
        /*returns status and counters for FIFO for following options:
            get_rx_sync_status,
            get_rx_en_status,
            get_rx_fifo_discard_count,
            get_rx_decode_error_count
        */
        if(self != NULL){
            Py_INCREF(self);
            status = PyObject_CallMethod(self,option,NULL);// not null
            if(status==NULL){
                std::stringstream ss;
                ss<<"[ERROR] Could not obtain RX status for option"<<option;
                flush_debug(ss.str());
                return emptyList();
            }
        }else{
            flush_debug("readFifoStatus::self object is NULL!");
            return emptyList();
        }
    
        Py_DECREF(self);
        PyGILState_Release(gstate);
    
        return status;
    };
   
    // returns sum of tpx3 receiver fifo sizes (of N receivers available)
    long int get_rx_fifo_size(PyObject *chip){
    
        long int total_fifo_size = 0;
    
        for(int i=0; i<8; i++){
   
            // makeregister name of RX<Nr> [0-8] 
            std::string rxname = "RX"+std::to_string(i);
    
            Py_INCREF(chip);
            // retrieve receiver object from self.chip object
            PyObject *rx = PyObject_GetItem(chip, PyUnicode_FromString(rxname.c_str()));
            Py_INCREF(rx);
            // get pyobject of the fifo size in the register
            PyObject *fifosize = PyObject_GetItem(rx,PyUnicode_FromString("FIFO_SIZE"));
            // convert pyobject to suitable c++ format
            long int ith_fifo_size = PyLong_AsLong(fifosize);
            // add ith RX counts to total
            total_fifo_size += ith_fifo_size;
        
            Py_DECREF(fifosize);
            Py_DECREF(rx);
            rx=NULL;
            fifosize=NULL;
            Py_DECREF(chip);
        }
        // return sum of counts
        return total_fifo_size;
    
    }

    // should it stay ??? 
    void checkIsRunning(PyObject *self){ // not used
    
        PyObject *isRunning = PyObject_GetAttrString(self, "_is_running");
        std::cout<<"\n\t[self._is_runnning]="<<PyObject_IsTrue(isRunning)<<"\n"<<std::flush;
        Py_DECREF(isRunning);
        isRunning=NULL;
    
    }
    
    // checks if self._set_calculate object is set true in tpx3/fifo_readout.py
    bool checkSelfCalculate(PyObject *self){// check if needed []
    
        PyObject *calc = PyObject_GetAttrString(self, "_calculate");
        bool status = false;
    
        if(calc!=NULL){
            PyObject *is_set = PyObject_GetAttrString(calc,"is_set");
            if(is_set!=NULL){
                status = PyObject_IsTrue(is_set);
            }else{
                flush_debug("Could not get attr=\"_calculate.is_set()\"");
            }
            Py_DECREF(is_set); 
        }else{
            flush_debug("Could not get attr=\"_calculate\"");
        }    
    
        Py_DECREF(calc);
        return status;
    
    }
    
    //clears self._set_calculate object in tpx3/fifo_readout.py
    void clearSelfCalculate(PyObject *self){ // not used, check necessity []
        
        PyObject *calc = PyObject_GetAttrString(self, "_calculate");
        bool status = false;
    
        if(calc!=NULL){
            PyObject *clear = PyObject_CallMethod(calc,"clear",NULL); 
            if(clear==NULL){
                flush_debug("Could not call [\"self._calculate.clear()\"]");
            }
            Py_DECREF(clear); 
        }else{
            flush_debug("Could not get attr=\"_calculate\"");
        }    
    
        Py_DECREF(calc);
    
    }
    
    // sets total number of recorded words during the total runtime to corresponding oject 
    // on the python level at self._record_count 
    // (mainly ofr the function tpx3::fifo_readout::print_status)
    void setSelfRecordCount(PyObject *self){
        
        // enable GIL lock 
        PyGILState_STATE gstate = PyGILState_Ensure();
          
        // get pointer to the object of local counter of total words recorded 
        PyObject *total_words = PyLong_FromLong(global_words_recorded);
        // set pyobject referece to the python level object of _record_count
        PyObject_SetAttrString(self,"_record_count", total_words);
        Py_DECREF(total_words);// release local pointer to the counter
    
        PyGILState_Release(gstate);
    
    }
   
    // returns no data timeout from the python level defined in the kwargs of 
    // tpx3::scan_base::start_readout and tpx3::fifo_readout::readout functions
    float getNoDataTimeout(PyObject *self){
      
        // get the pyobject of self.no_data_timeout
        PyObject *nd_timeout = PyObject_GetAttrString(self, "no_data_timeout");
        // return zero if object is NULL
        if(nd_timeout==NULL){
           //flush_debug("could not find no_data_timeout!");
           Py_XDECREF(nd_timeout);
           return 0.0; 
        }  
        // return the actual number if the pointer not nil and object is float 
        if(PyFloat_Check(nd_timeout)){ 
            float timeout = PyFloat_AsDouble(nd_timeout);
            Py_DECREF(nd_timeout);
            return timeout;
        }
        // if pointer type is not of expected type - return zero
        if(!(nd_timeout==nd_timeout)){
            Py_DECREF(nd_timeout);
            return 0.0;
        }
        Py_XDECREF(nd_timeout);
        nd_timeout=NULL;
    
        return 0.0;
    }
    
    // returns boolean flag of the sate of the object self.stop_readout in 
    // tpx3::fifo_readout:__init__
    // necesary to stop the readout in between the iterations of the scans 
    bool local_flagStopReadout(PyObject *self) {

        // get pointer reference to the object self.stop_readout    
        PyObject *stop_readout = PyObject_GetAttrString(self, "stop_readout");
        //check if reference is zero (reference exists)
        if (stop_readout != NULL) { //if pinter exists
            // get object reference ofr the threading.Event status object "is_set"
            PyObject *is_set_method = PyObject_GetAttrString(stop_readout, "is_set");
            if (is_set_method != NULL) {// "is_set exists"
                // calling self.stop_readout.is_set()
                PyObject *is_set_result = PyObject_CallObject(is_set_method, NULL);
                if (is_set_result != NULL) {//if method returns something
                    // check if result obejct is of type true
                    if (PyObject_IsTrue(is_set_result)) {
                        //flush_debug("self.stop_readout.is_set=[TRUE]");
                        Py_DECREF(is_set_result);
                        Py_DECREF(is_set_method);
                        Py_DECREF(stop_readout);
                        //release pointers and return true 
                        return true;
                    } else {
                        //flush_debug("self.stop_readout.is_set=[FALSE]");
                        Py_DECREF(is_set_result);
                        Py_DECREF(is_set_method);
                        Py_DECREF(stop_readout);
                        //release refrences and return false
                        return false;
                    }
                } else {
                    flush_debug("[DEBUG] self.stop_readout.is_set() call failed");
                    Py_DECREF(is_set_method);
                    Py_DECREF(stop_readout);
                    return false;
                }
            } else {
                flush_debug("[DEBUG] self.stop_readout.is_set method is NULL");
                Py_DECREF(stop_readout);
                return false;
            }
        } else {
            flush_debug("[DEBUG] self.stop_readout IS NULL");
            return false;
        }
        // if at any point if this tree the reference is null - return false
    }
    
    // calls tpx3::fifo_readout::get_data method and returns 
    // std:;vector<uint32_t> of data words
    std::vector<uint32_t> localQuerryFifo(PyObject *chipFifo){
        //ensure GIL lock
        PyGILState_STATE gstate = PyGILState_Ensure(); 
        // initialize pointer to python numpy.ndarray recorded by the get_data() method
        PyArrayObject *result; 
    
        if(chipFifo != NULL){// if self.chip['FIFO'] is not null
            // increase and then decrease reference to *chipFifo
            // working version below
            Py_INCREF(chipFifo); 
            result = reinterpret_cast<PyArrayObject*>(PyObject_CallMethod(chipFifo,"get_data",NULL));
            Py_DECREF(chipFifo); 
        }else{ // if can not access self.chip['FIFO'] object 
            flush_debug("localQuerryFifo:[Error] chipFifo is null\nOR Object was not passed correctly");
            result = 0;
            // TODO: maybe put in empty PyArray_DATA..?
        }
    
        // increase reference counter and transfer data from python side to c++ 
        // via fillVectoFromPyData (defined above)
        Py_INCREF(result);
        std::vector<uint32_t> result_vector = fillVectorFromPyData(result);// technically works
        Py_DECREF(result);
        //delete pointer after data is transfered
        result = NULL;
        //release GIL lock
        PyGILState_Release(gstate); 
    
        return result_vector;
    
    };
    
    // function records uint32_t words from basil interface via "localQuerryFifo"
    // within user-defined time of the chunk "interval" in [ms] units
    //
    std::vector<uint32_t> local_getFifoData(PyObject *chipFifo, float interval){
    
        //ensure GIL lock
        PyGILState_STATE gstate = PyGILState_Ensure(); 
    
        // initialize data vectors:
        std::vector<uint32_t> temp; // temp vector to store data from single querry call 
        std::vector<uint32_t> fifo_data; // vector of dat of full chunk
        
        long int cnt = 0;  
     
        ///////////////////////////////////////////////////////////////
        // start measurng time 
        auto start_time = std::chrono::high_resolution_clock::now();
        // while loop breaks when its runtime, checked by time_to_read is leser than t_interval
        while(time_to_read(start_time) < interval){ 
            Py_INCREF(chipFifo);
            //recording data from single querry to basil in a temp vector
            std::vector<uint32_t> temp = localQuerryFifo(chipFifo);
            Py_DECREF(chipFifo);
            if(temp.size()>0){// if querry yielded some data
                // append data of this iteration of querry to the main chunk vector
                fifo_data.insert(fifo_data.end(),
                                 temp.begin(),
                                 temp.end());// should use temp vector here
            }
            cnt++;
        }
        //clear temp vector for good measure
        temp.clear();
   
        // release GIL state
        PyGILState_Release(gstate); 
        // return data chunk
        return fifo_data;
    
    };
   
    // checks size of tpx3 FIFO. IN actuality checks size of the 
    // basil.SiTcp._read_buffer deque object from which 
    // tpx3::fifo_readout::read_data() collects data words as numpy.ndarray
    void checkFifoSize(PyObject *fifo){
    
        PyObject *fsize = PyObject_GetItem(fifo,PyUnicode_FromString("FIFO_SIZE"));
        long int fifoSize = PyLong_AsLong(fsize);
        std::cout<<"[DEBUG] FIFO_SIZE >> "<<fifoSize<<" words\n"<<std::flush;
        Py_DECREF(fsize);
        fsize=NULL;
    
    };
    
    // resets error counter of tpx3 receivers via call to self.rx_error_reset 
    // method in tpx3::fifo_readout
    void resetRXErrorCounters(PyObject *self){

        // call self.rx_error_reset    
        PyObject *reg = PyObject_CallMethod(self,"rx_error_reset",NULL);
        if(reg==NULL){// if return object is null - call failed
            flush_debug("resetRXErrorCounters: failed to reset RX!");
        }
        Py_DECREF(reg);
        reg=NULL;
    
    };
    
//    // dumps status of status registers of all receivers on tpx3
//    // test function
//    PyObject* getStatusAllRx(PyObject *self, const char* REG){ // this one's not used
//    
//        // Can be used to directly obtain:
//        // [ENABLE, DATA_DELAY, INVERT, SAMPLING_EDGE, READY] register values
//        // of RX channels
//    
//        PyGILState_STATE gstate = PyGILState_Ensure(); 
//        PyObject *chip; // top object reference
//        PyObject *STAT; // statusobject reference
//    
//        std::vector<int> rx_status;
//     
//        chip = PyObject_GetAttrString(self, "chip");
//        if(chip!=NULL){
//            for(int i=0; i<tpx3::NUM_RX_CHAN ; i++){ //TODO: put nr of RX channels in a const file.
//                std::string ithRX = "RX"+std::to_string(i);
//                PyObject *irx = PyObject_GetItem(chip, PyUnicode_FromString(ithRX.c_str()));
//                if(irx != NULL){
//                    STAT = PyObject_GetItem(irx, PyUnicode_FromString(REG));           
//                    if(STAT != NULL){
//                        int status = PyLong_AsLong(STAT);
//                        rx_status.push_back(status); 
//                        Py_DECREF(STAT);
//                        Py_DECREF(irx);
//                    }else{
//                       std::stringstream ss;
//                       ss<<"Pointer to (" << REG << ") is NULL!";
//                       flush_debug(ss.str());
//                       return Py_None;
//                    }
//                }else{
//                   flush_debug(ithRX+" is NULL");
//                   return Py_None;    
//                }
//            }
//        }else{
//            flush_debug("chip is NULL");
//            return Py_None;
//        }
//        std::stringstream stream;
//        stream<<"Dumping "<<REG<<"]:\n";
//        flush_debug(stream.str());
//        dumpVector(rx_status);
//        rx_status.clear();
//    
//        Py_DECREF(chip);
//    
//        PyGILState_Release(gstate); 
//        return Py_None;
//    };
    
    std::string checkReadStatus(PyObject *self){ // maybe could be used...
        
        PyGILState_STATE gstate = PyGILState_Ensure();
        const char* status;
        PyObject *stopreadout = PyObject_GetAttrString(self,"stop_readout");
        if(stopreadout!=NULL){
            PyObject_Print(stopreadout,stdout,0);
            /////////////////////
            PyObject *isSet = PyObject_CallMethod(stopreadout,"is_set",NULL);
            //const char* temp = 
            std::cout<<"[DEBUG] CHEKCING self.stop_readout.isSet() = "
                     <<std::atoi(PyUnicode_AsUTF8(PyObject_Str(isSet)))<<"\n"<<std::flush;
            ////////////////////////
            status = PyUnicode_AsUTF8(PyObject_Str(stopreadout));
            //std::cout<<"[DEBUG] CHCEKING readout FLAG <pepelaugh> = "<<status<<"\n"<<std::flush;
        }else{
            flush_debug("[ERROR]: stopreadout object is NULL!");
            return "HUYASE!";
        };
        Py_DECREF(stopreadout);
        stopreadout = NULL;
        PyGILState_Release(gstate);
        return status;
    
    };

    // checks the sate of SHUTTER register of the tpx3
    // returns 1 if shuter is open, 0 if closed    
    int checkSHUTTER(PyObject *self){ // could be used
    
        // ensure GIL lock
        PyGILState_STATE gstate = PyGILState_Ensure();
        // get self.chip object
        PyObject *chip = PyObject_GetAttrString(self, "chip");
        // initialize shutter status number 
        int status;
        if(chip != NULL){// if self.chip exists
            //get pointer to to the item of the self.chip['CONTROL'] in self.chip dictionary
            PyObject *control = PyObject_GetItem(chip, PyUnicode_FromString("CONTROL"));
            if(control != NULL){// if ['CONTROL'] item exists
                // get pointer to value of ['SHUTTER'] keyword register
                PyObject *shutter = PyObject_GetItem(control, PyUnicode_FromString("SHUTTER"));
                // convert and assign number at the pointer to the c++ format
                status = std::atoi(PyUnicode_AsUTF8(PyObject_Str(shutter)));
                Py_DECREF(shutter);
            }else{ // if self.chip['CONTROL'] was not found - return failure (-1)_
                flush_debug("[ERROR] checkSHUTTER(): SHUTTER IS NULL!");
                Py_XDECREF(control);
                Py_XDECREF(chip); //might be redundant
                return -1;
            }
            Py_DECREF(control);
            control = NULL;
        }else{// if self.chip was not found - return failure (-1)_
            flush_debug("[ERROR] checkSHUTTER(): CONTROL IS NULL!");
            Py_XDECREF(chip);
            return -1;
        }
    
        Py_DECREF(chip);
        //release GIL lock
        PyGILState_Release(gstate); 
        //return 1 on open shutter, 0 on closed shutter
        return status;
    
    };
    
    // test/debug function 
    // can be used to set tpx3 registers defined in tpx3.yml for tpx3 object (not tpx3_rx)
    //
    PyObject *setRegister(PyObject *self, const char* SETTING, int value){ //could be used...
    
        PyGILState_STATE gstate = PyGILState_Ensure(); 
        PyObject *chip; // top object reference
        PyObject *REG; // control reg object reference
    
        chip = PyObject_GetAttrString(self, "chip");
    
        if(chip!=NULL){
            for(int i=0; i<tpx3::NUM_RX_CHAN ; i++){
                std::string ithRX = "RX"+std::to_string(i);
                PyObject *irx = PyObject_GetItem(chip, 
                                                 PyUnicode_FromString(ithRX.c_str()));
                if(irx != NULL){
                    Py_INCREF(chip);
                    PyObject_SetItem(chip,REG,PyLong_FromLong(value));
                    std::stringstream ss;
                    ss<<"Set ["<<SETTING<<"]="<<PyLong_FromLong(value);
                    flush_debug(ss.str());
                    ss.clear();
                    Py_DECREF(REG);
                }else{
                    std::stringstream ss;
                    ss<<"Register "<<SETTING<<"["<<REG<<"] is NULL";
                    flush_debug(ss.str());
                    ss.clear();
                    return Py_None;             
                }
            }
        }else{
            flush_debug("chip is NULL");
            return Py_None;
        }
    
        Py_DECREF(chip);
    
        PyGILState_Release(gstate); 
        return Py_None;
    
    };
    
    //---------------------------------------------------------------------------------
    //
    //  Substitute for tpx3/fifo_readout/readout function
    //
    //---------------------------------------------------------------------------------
    // a complete replica of tpx3::fifo_readout::readout function that runs in a separate thread
    // in tpx3::fifo_readout::start.
    // this function reads data words from basil buffer and assembles characteristic data
    // into python tuple and transfers it directly to tpx3::fifo_readout::worker function(thread)
    //  
    //  input:
    //  PyObject *self -> python object of the tpx3::fifo_readout class
    //  PyObject *deque -> python object of the tpx3::fifo_readout::__init__::_data_deque() 
    //  float interval -> readout interval in [ms]
    PyObject* readoutToDeque(PyObject *self, PyObject *deque, float interval){

      //ensure GIL state locked
      PyGILState_STATE gstate = PyGILState_Ensure();
      // importing numpy API functionality
      import_array();
     
      //flush_debug("tpx3::cpp::tpx3redoutcpp: [INFO] starting readout"); 
    
      //if(l_global_start_time == std::chrono::_V2::system_clock::time_point()){
      //  l_global_start_time == tick();
      //}

      // defining last chunk read time, current time, time to wait (not implemented yet)
      double last_time = 0.0, curr_time = 0.0, timestamp = 0.0; 
      float time_wait = 0.0;
    
      // initializing
      int glob_tries = 0; // global number of readout iterations 
      int prev_shutter = 0; // tracker of the SHUTTER state
      float noDataTimeout = 0.0; // no data timeout setting 
      std::vector<float> iter_times;
    
      // getting common global objects
      // self.chip
      PyObject *chip = PyObject_GetAttrString(self,"chip");
      //self.chip['FIFO']
      PyObject *fifo = PyObject_GetItem(chip,PyUnicode_FromString("FIFO"));
    
      // checking user-defined timeout when no data is available
      //
      Py_INCREF(self);
      noDataTimeout = getNoDataTimeout(self);
      Py_DECREF(self);
      // checking if user enabled self._set_calculate option
      Py_INCREF(self);
      bool eba = checkSelfCalculate(self);
      Py_DECREF(self);
    
      //std::cout<<"[debug] is \"calculate\" set? = "
      //         <<ansi_orange<<"["<<eba<<"]"<<ansi_reset<<"\n"<<std::flush;  
    
      // get current time     
      curr_time = get_float_time();
    
      // start of global readout loop where data is assembled in chunks
      // during runtime of t_interval, [milliseconds]
      while (true){/*loop breaks via "stop_readout" signal below*/
   
        // fixing start of the readout  
        auto iter_start = tick();
    
        // checking SHUTTER state 
        Py_INCREF(self);
        int shutter_now = checkSHUTTER(self); // ~6 us
        Py_DECREF(self);
    
        // Allocating memory for the output tuple for a fixed number of items
        PyObject *tuple = PyTuple_New(tpx3::NUM_DEQUE_TUPLE_ITEMS+2);
        if(tuple==NULL){
            flush_debug("Could not instantiate data tuple");
            PyGILState_Release(gstate); 
            return NULL;
        };
    
        // initializin' pointers to the data and rx error counter objects
        //
        PyObject *fifoData, *fifo_decerr, *fifo_discerr;
      
        // recoding  fifo data locally
        //
        Py_INCREF(fifo);
        std::vector<uint32_t> data = local_getFifoData(fifo,interval);// usually conforms to
                                                                      // interval set by user
                                                                      // with acc-cy of O(us)
        Py_DECREF(fifo);
   
        long int n_words = data.size();
    
        // recording number of collected unt32_t data words
        //
        global_words_recorded += n_words;
        //creating python object ref to currently recorded number of data words
        PyObject *py_n_words = PyLong_FromLong(n_words);      
        
        // put recored words from vector.data() to numpy array
        //
        fifoData = fillPyArray(data); // ~ 40 us @ 100ms interval
    
        /////////////////////
        // reading fifo discard error counters
        Py_INCREF(self); 
        fifo_discerr = readFifoStatus(self,tpx3::CNTR_DISCARD_ERR);
        Py_DECREF(self); 
    
        // reading fifo decode error counters
        Py_INCREF(self); 
        fifo_decerr = readFifoStatus(self,tpx3::CNTR_DECODING_ERR);
        Py_DECREF(self); 
    
        ///////////////////////////////
        // counting numbers in error lists for discard errors
        Py_INCREF(fifo_discerr);
        long int n_errdisc = getNerror(fifo_discerr);
        Py_DECREF(fifo_discerr);
       
        // counting numbers in error lists for decode errors
        Py_INCREF(fifo_discerr);
        Py_INCREF(fifo_decerr);
        long int n_errdec = getNerror(fifo_decerr);
        Py_DECREF(fifo_decerr);
        ///////////////////////////////
    
        // updating timestamp
        float t_end = timeDiff(iter_start,tick());//---------------> ~53ms from start
    
        //adjusting timestamps for chunk sart/stop and current time
        upd_time(curr_time, last_time, timestamp);
    
        // assemble python tuple that will be sent to fifo_readout::worker function
        // which later arrives at scan_base::handle_data via self.callback 
        // call on the python side
        PyTuple_SetItem(tuple, 0, fifoData); // uint32_t data
        PyTuple_SetItem(tuple, 1, PyFloat_FromDouble(last_time)); //chunk start time
        PyTuple_SetItem(tuple, 2, PyFloat_FromDouble(curr_time)); // chunk stop time
        PyTuple_SetItem(tuple, 3, PyLong_FromLong(n_errdisc)); // n discard errors
        PyTuple_SetItem(tuple, 4, PyLong_FromLong(n_errdec)); // n decoder err
    
        // Adding python tuple to deque:
        Py_INCREF(deque);
        PyObject *result = PyObject_CallMethod(deque, "append", "(O)", tuple);
        if(result==NULL){// if can not access self._data_deque() return error and nullptr
            flush_debug("[WARNING] CAN NOT access self._data_deque!");
            Py_DECREF(result);
            Py_DECREF(self);
            PyGILState_Release(gstate);
            return NULL;
        }
        // < 1.0 us to set tuple items and append to deque
        
        // releasing pointers related to self._data_deque()
        Py_DECREF(result);
        Py_DECREF(tuple);
        Py_DECREF(deque);
        // setting obj-ts to NULL to save memory (EXCEPT deque)!
        result = NULL;
        tuple = NULL;
        fifoData = NULL;
        fifo_discerr = NULL;
        fifo_decerr = NULL;
    
        //////////////////////////////////////////////////////////////
        // Add number of words recorded in ith iteration of while loop
        // to calculate average reading rate on the python level 
        Py_INCREF(self); //+self 
        // get pointer to self._words_per_read object
        PyObject *wpr = PyObject_GetAttrString(self,"_words_per_read");
        // calling self._words_poer_read.append(py_n_words)
        PyObject *wpread = PyObject_CallMethod(wpr, "append", "O", py_n_words);
    
        if(wpread==NULL){// if method call fails do nothing
            flush_debug("[PyObject ISSUE] Could not append to _words_per_read!\n");
        }
        // release pointers and delete data related to self._words_per_read to save memory
        // and prevent leaks
        Py_DECREF(self); //-self
        Py_DECREF(wpread);// -wpread
        Py_DECREF(wpr);// -wpread
        wpread=NULL;
        wpr=NULL;
        //////////////////////////////////////////////////////////////
        
        // reset rx error counters if detected any discard
        // or decoding errors
        if(n_errdisc != 0 || n_errdec !=0){ //if any decoding or discard error happens
            std::cout<<ansi_red<<"[DEBUG] Detected: (N_DISCARD="<<n_errdisc
                     <<"), (N_DECODE="<<n_errdec<<") errors.\n"<<ansi_reset<<std::flush;
    
            Py_INCREF(self);
            // calling self.rx_error_reset
            resetRXErrorCounters(self);
            Py_DECREF(self);
            // reseting local counters
            n_errdisc = 0;
            n_errdec = 0;
            
        }
    
        // calculating time to wait for befre next iteration of while loop     
        float t_iteration = timeDiff(iter_start,tick());
        time_wait = interval - t_iteration;
   
        // checking state of self.stop_readout flag  
        if(local_flagStopReadout(self)){ //if self.stop_readout is true - break the loop
           // clear std::vector<uibnt32_t> for good measure...
           data.clear();
           break;
    
        }
    
        //assign current shutter state to be previous shutter state since all vital
        //operations are have been performed
        prev_shutter = shutter_now;
    
        //////////////////////////////////////////
    
        // clear std::vector<uint32_t> of recorded data words to avoid duplicates/overwrites
        data.clear();
    
        // increment global tries
        glob_tries++;
           
      };// end the while(true) loop
      
      // set self._record_count parameter for txp3::fifo_readout::print_status
      Py_INCREF(self);
      setSelfRecordCount(self);
      Py_DECREF(self);
    
      // adding Py_None tuple here that will stop tpx3::fifo_readout::worker function
      // and stop readout
      //
      //initialize pytuple with a singe litem of None type
      PyObject *lastTuple = PyTuple_Pack(1,Py_None);
      // appending None data to ypx3::fifo_readout::_data_deque()
      PyObject *append = PyObject_CallMethod(deque, "append", "O", lastTuple);
      if(append == NULL){
        std::stringstream warn_msg;
        warn_msg<<"tpx3::cpp::tpx3readoutcpp::readoutToDeque: [WARNING]"
                <<"Could not append to self._data_deque()";
        flush_debug(warn_msg.str());
      }
      
      // releasing tuple objectobjects
      Py_DECREF(append);
      Py_DECREF(lastTuple);
      append = NULL;
      lastTuple=NULL;
    
      //flush_debug("tpx3::cpp::tpx3readoutcpp:[INFO] - readoutToDeque DONE!"); 
    
      // release GIL state
      PyGILState_Release(gstate);   
      ///// finita la comedia! //////    
      return Py_None;
    
    };

// end of class

}

