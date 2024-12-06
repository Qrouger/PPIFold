"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='PPIFold',
    version=0.0.1,
    description=(
        'Pipeline that automates AlphaPulldown'
    ),
    author='Quentin Rouger',
    author_email='quentin.rouger@univ-rennes.fr',
    license='GPL-3.0 license',
    url='https://github.com/Qrouger/PPIFold',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'alphapulldown'
    ]
)
