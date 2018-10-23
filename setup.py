from setuptools import setup, find_packages
#import setuptools
import sys
#mport setuptools, sys
if sys.version_info[0] > 2:
    raise "Invalid Python interpreter, must be 2.X"

setup(
    name='quanttb',
    version='0.1dev',
    packages=find_packages(),
    description='Quantification of mixed TB infections',
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.txt').read(),
    install_requires=[
         'numpy', 'subprocess32','sklearn', 'scipy'      ],
    package_dir={'quanttb': 'quanttb/'},
    package_data={'quanttb': ['data/*']},
    #include_package_data = True,
    #scripts=['scripts/classify.py', 'scripts/runnucmer.py', 'scripts/quanttb.py'],
    entry_points={
        'console_scripts': [
            'quanttb = quanttb.scripts.quanttb:main'
        ]},
    

)