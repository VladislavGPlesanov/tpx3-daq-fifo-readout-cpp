# tpx3-daq
[![Build Status](https://dev.azure.com/SiLab-Bonn/tpx3-daq/_apis/build/status/SiLab-Bonn.tpx3-daq?branchName=master)](https://dev.azure.com/SiLab-Bonn/tpx3-daq/_build/latest?definitionId=1&branchName=master)

DAQ for the [Timepix3](https://medipix.web.cern.ch/technology-chip/timepix3-chip) chip based on the [Basil](https://github.com/SiLab-Bonn/basil) framework.

### Installation

Description below works on Ubuntu OS of versions > LTS 12.04

- For the first time installation update existing software on the system with `sudo apt update` and `sudo apt upgrade`

Install main libraries for the tpx3 software:

```bash
sudo apt install curl libgirepository1.0-dev gcc libcairo2-dev pkg-config python3-dev gir1.2-gtk-3.0
```

- Install [conda](https://conda.io/miniconda.html) for python

```bash
mkdir ~/miniconda
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh
bash miniconda.sh -u -b -p ~/miniconda
```
then export miniconda executable path to the enviroment:
```bash
export PATH=<PATH_TO_MINICONDA>/miniconda/bin:$PATH
```
Update conda and install dependencies:
```bash
conda update --yes conda
conda install --yes numpy bitarray pytest pytest-cov pyyaml scipy numba pytables pyqt matplotlib tqdm pyzmq blosc psutil setuptools
```
- Instaling the rest of dependencies. 
Basil backend:
```bash
pip install basil_daq==3.2.0
```
and GUI:
```bash
pip install pycairo
pip install PyGObject
```
 - Finally, clone the repository nad install it
```bash
mkdir tpx3-daq
git clone https://github.com/SiLab-Bonn/tpx3-daq.git tpx3-daq/
cd tpx3-daq
pip install -e .
```
In case the installation of `pycairo` fails, then try to additionally install
```bash
sudo apt-get install libcairo2-dev libjpeg-dev libgif-dev
```
and rerun `pip install pycairo`

<!--- Install dependencies and tpx3-daq:
```
conda install numpy bitarray pyyaml scipy numba pytables matplotlib tqdm pyzmq blosc psutil
pip install git+https://github.com/SiLab-Bonn/tpx3-daq.git@master
```
-->

### Usage

- Flash appropriate bit file (TBD)

- Manually configure ethernet connection according to [instructions](https://tpx3-daq.readthedocs.io/en/latest/usage.html#installation)

- Run a scan from CLI (TBD):
```
tpx3 tune_noise
```
- Start GUI via:
```bash
tpx3_gui
```
- For help, run (TBD):
```
tpx3 --help
```
- Use the online monitor:
```
tpx3_monitor
```
In case monitor fails to start up try installing:

```bash
sudo apt install libxcb-xinerama0
```


