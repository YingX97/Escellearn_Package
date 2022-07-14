from setuptools import setup

setup(
    name='cell_atlas_portal_2022',
    version='0.1',
    description='a webportal for exploring and visualising single cell data from neonatal lung',
    url='',
    author='YING XU',
    author_email='yingxu0928@gmail.com',
    license='MIT',
    packages=['cell_atlas_portal_2022'],
    install_requires=[
        'h5py',
        'numpy',
        'pandas',
        'scipy'
    ],
    include_package_data=True,
    zip_safe=False
)