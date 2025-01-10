"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='PPIFold',
    version='0.1.5',
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
        'alphapulldown==2.0.0',
        'seaborn',
        'urllib3',
        'matplotlib',
        'scipy'
    ],
    entry_points={'console_scripts': ['PPIFold=PPIFold.PPIFold:main',],}
)
