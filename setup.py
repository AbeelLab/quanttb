from setuptools import setup, find_packages
import sys
if sys.version_info[0] > 2:
    raise Exception( "Invalid Python interpreter, must be 2.X")

setup(
    name='quanttb',
    version='0.1dev',
    packages=find_packages(),
    description='Quantification of mixed TB infections',
    license='GNU GENERAL PUBLIC LICENSE',
    long_description=open('README.md').read(),
    install_requires=[
         'numpy==1.14.5', 'subprocess32==3.5.2', 'scipy==1.1.0'],
    package_dir={'quanttb': 'quanttb/'},
    package_data={'quanttb': ['data/*']},
    #include_package_data = True,
    #scripts=['scripts/classify.py', 'scripts/runnucmer.py', 'scripts/quanttb.py'],
    entry_points={
        'console_scripts': [
            'quanttb = quanttb.scripts.quanttb:main'
        ]},
    

)