from setuptools import setup, find_packages
import os
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = os.path.join(lib_folder, '/requirements.txt')
install_requires = []
if os.path.isfile(requirement_path):
    with open(requirement_path, 'r', encoding='utf-8') as f:
        install_requires = f.read().splitlines()

setup(
    name='evpfft_gui',
    python_requires='>=3.8',
    author='Sarah Hankins',
    author_email='shankins@lanl.gov',
    version='0.1',
    description='Package for running LANL\'s EVPFFT software with a graphical front end.',
    packages=find_packages(exclude=['ez_setup', 'tests', 'tests.*']),
    include_package_data=True,
    install_requires=install_requires,
    entry_points={
        'console_scripts' : [
            'evpfft-gui = evpfft_gui.gui:main'
        ]
    }
)