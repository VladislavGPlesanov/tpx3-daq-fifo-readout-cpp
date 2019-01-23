#!/usr/bin/env python

from __future__ import absolute_import
from setuptools import setup
from setuptools import find_packages

import tpx3
from symbol import except_clause

version = tpx3.__version__

author = ''
author_email = ''

# Requirements
install_requires = ['basil-daq==2.4.11', 'bitarray>=0.8.1', 'matplotlib',
                    'numpy', 'online_monitor>=0.4.0<0.5',
                    'pixel_clusterizer==3.1.3', 'tables', 'pyyaml', 'pyzmq',
                    'scipy', 'numba', 'tqdm']
setup(
    name='tpx3-daq',
    version=version,
    description='DAQ for TimePix3 ASIC',
    url='https://github.com/SiLab-Bonn/tpx3-daq',
    license='',
    long_description='',
    author=author,
    maintainer=author,
    author_email=author_email,
    maintainer_email=author_email,
    install_requires=install_requires,
    packages=find_packages(),
    setup_requires=['online_monitor>=0.4.0<0.5'],
    include_package_data=True,
    platforms='any',
    entry_points={
        'console_scripts': [
            'tpx3_monitor = tpx3.online_monitor.start_tpx3_monitor:main'
        ]
    },

)

try:
    import os
    from online_monitor.utils import settings
    # Get the absoulte path of this package
    package_path = os.path.dirname(tpx3.__file__)
    # Add online_monitor plugin folder to entity search paths
    settings.add_producer_sim_path(os.path.join(package_path,
                                                'online_monitor'))
    settings.add_converter_path(os.path.join(package_path,
                                             'online_monitor'))
    settings.add_receiver_path(os.path.join(package_path,
                                            'online_monitor'))
except ImportError:
    pass
