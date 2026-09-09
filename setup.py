"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='PPIFold',
    version='1.0.7',
    description=(
        'Automatic pipeline using AlphaPulldown to predict PPI and homo-oligomer'
    ),
    author='Quentin Rouger',
    author_email='quentin.rouger@univ-rennes.fr',
    license='GPL-3.0 license',
    url='https://github.com/Qrouger/PPIFold',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'alphapulldown==2.2.0',
        'seaborn',
        'matplotlib',
        'scipy',
        'graphviz',
        'numpy==1.26.4',
        'gemmi',
        'nvidia-cublas-cu12==12.6.4.1',
        'nvidia-cuda-cupti-cu12==12.6.80',
        'nvidia-cuda-nvcc-cu12==12.8.93',
        'nvidia-cuda-nvrtc-cu12==12.6.77',
        'nvidia-cuda-runtime-cu12==12.6.77',
        'nvidia-cudnn-cu12==9.5.1.17',
        'nvidia-cufft-cu12==11.3.0.4',
        'nvidia-cufile-cu12==1.11.1.6',
        'nvidia-curand-cu12==10.3.7.77',
        'nvidia-cusolver-cu12==11.7.1.2',
        'nvidia-cusparse-cu12==12.5.4.2',
        'nvidia-cusparselt-cu12==0.6.3',
        'nvidia-nccl-cu12==2.26.2',
        'nvidia-nvjitlink-cu12==12.6.85',
        'nvidia-nvtx-cu12==12.6.77'
    ],
    entry_points={'console_scripts': ['PPIFold=PPIFold.PPIFold:main',],}
)
