
for i in {1..10}
do
    echo " "
done
echo " ================ [COMPILING C++ readout] ================"
echo " "

if which python3; then
    echo "found python3 "
fi

PYVERSION=$(python3 -c 'import sys; ver=sys.version_info[0:3]; print("{0}.{1}".format(*ver))')

# checking if python headers and core file paths are included in env
if env | grep CPATH; then 
    echo "CPATH Path exists"
else
    echo exporting CPATH
    export CPATH=/usr/include/python$PYVERSION/:$CPATH
fi

if env | grep LD_LIBRARY_PATH; then
    echo "LD_LIB exists"
else
    echo "exporting LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH=/usr/lib/:$LD_LIBRARY_PATH
fi

##############################################################################
NUMPY_VER=$(python3 -c 'import numpy; print(numpy.__version__)')
echo "found numpy version $NUMPY_VER"
NUMPY_INCLUDE=$(python3 -c 'import numpy; print(numpy.get_include())') 
NUMPY_LIB=$(python3 -c 'import numpy; print(numpy.__file__.rsplit("/",1)[0] + "/core/lib")')

echo $NUMPY_INCLUDE
echo $NUMPY_LIB

g++ tpx3readout.cpp \
-o tpx3readoutcpp.so \
-shared -std=c++14 \
-fPIC \
-fpermissive \
-I $NUMPY_INCLUDE \
-L $NUMPY_LIB
##############################################################################
echo " "
echo " ================ [END] ================"
for i in {1..10}
do
    echo " "
done

