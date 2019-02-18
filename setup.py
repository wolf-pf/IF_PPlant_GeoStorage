from setuptools import setup
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='IF_PPlant_GeoStorage',
      version='0.0.1 dev',
      description='Interface for coupled simulation of power plant and geological energy storage',
      url='https://github.com/project-angus/IF_PPlant_GeoStorage',
      author='Wolf Tilmann Pfeiffer, Francesco Witte',
      author_email='francesco.witte@hs-flensburg.de',
      long_description=read('README.rst'),
      license='GPL-3.0',
      packages=['coupled_simulation'],
      python_requires='>=3',
      install_requires=['TESPy >= 0.1.0',
                        'numpy >= 1.13.3',
                        'pandas >= 0.19.2',
                        'scipy >= 0.19.1'])
