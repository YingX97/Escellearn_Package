from setuptools import setup

setup(
    name='Escellearn',
    version='0.1',
    description='a backend for exploring and visualising single cell data from neonatal lung',
    url='',
    author='YING XU',
    author_email='yingxu0928@gmail.com',
    license='MIT',
    packages=['Escellearn'],
    install_requires=[
        'h5py',
        'numpy',
        'pandas',
        'scipy'
    ],
    include_package_data=True,
    zip_safe=False
)