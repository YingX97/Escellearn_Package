# Escellearn_Package

The Escellearn package contains serveral functions that perform the same functionality as the Escellearn webportal

To install and used the package locally:

1.1 -- Change directory to Escellearn

```bash
$ cd Escellearn
```



1.2 -- Install the package

```bash
$ pip install .
```



1.3 -- Use the package to explore the condensed neonatal lung cell atlas

```bash
$ python3
```



1.4 -- import the class

```python
from Escellearn import h5Reader
```



1.5 -- Use the class to read in the input file

```python
my_reader = h5Reader('condensed_lung_atlas_ordered.h5')
```



1.6 -- Example of using functions in h5Reader

```python3
my_reader.dataset_unified('Col1a1')
```





### Documentation

documentation of the package is in a static file located in `Escellearn/Escellearn/docs/_build/html/index.html`

