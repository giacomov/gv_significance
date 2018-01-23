from setuptools import setup

setup(
    name='gv_significance',
    version='1.0',
    packages=['gv_significance'],
    url='github.com/giacomov/gv_significance',
    license='BSD-3.0',
    author='Giacomo Vianello',
    author_email='giacomov@stanford.edu',
    description='Implement the formulae from Vianello (2018)',
    install_requires=['scipy >= 0.18',
                      'numpy',
                      'ncephes',
                      'numba']
)
