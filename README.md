[![DOI](https://zenodo.org/badge/116432941.svg)](https://zenodo.org/badge/latestdoi/116432941)

# gv_significance
Implements the formulae of Vianello 2018 for computing the significance in 
counting experiments.

# Installation

```bash
> git clone https://github.com/giacomov/gv_significance.git
> cd gv_significance
```

Then you can use `pip`:

```bash
> pip install .
```

or directly the `setup.py` file:

```bash
> python setup.py install
```

# Examples
See Vianello (2018) for details about these methods.

NOTE: all functions accept either single values or arrays as inputs, even 
though the following examples deal with single values only.

## Ideal case
Compute the significance of a measurement of 100 counts with an expected 
background of 90 in the ideal case of no uncertainty on the background:

```python
from gv_significance import ideal_case

ideal_case.significance(n=100, b=90)
```

Compute the the counts that a source must generate in order to be detected 
at 5 sigma with an efficiency of 90% given a background of 100 counts:

```python
from gv_significance import ideal_case

# Note: detection_efficiency must be either 0.5, 0.9 or 0.99
ideal_case.five_sigma_threshold(B=100, detection_efficiency=0.9)
```

## Poisson observation, Poisson background with or without systematics

Compute the significance of a measurement of 100 counts when a side measurement
returned a background of 900 counts with an exposure 10x the one of the observation
on-source (`alpha=0.1`):

```python
from gv_significance import poisson_poisson

poisson_poisson.significance(n=100, b=900, alpha=0.1)
```

Alternative, one can use the $Z_Bi$ estimator described in Cousins 2008 (but be 
careful about the range of usability, see Vianello 2018) by doing:

```python
from gv_significance import poisson_poisson

poisson_poisson.z_bi_significance(n=100, b=900, alpha=0.1)
```

If there is an additional systematic uncertainty on the background measurement of 
at most 10%, then:

```python
from gv_significance import poisson_poisson

poisson_poisson.significance(n=100, b=900, alpha=0.1, k=0.1)
```

If there is an additional systematic uncertainty on the background measurement
that is Gaussian-distributed with a sigma of 10%, then:

```python
from gv_significance import poisson_poisson

poisson_poisson.significance(n=100, b=900, alpha=0.1, sigma=0.1)
```

## Poisson observation, Gaussian background
Compute the significance for an observation of 100 counts when we have a
background estimation procedure that gives us an expected background of
`90 +/- 2.4`:

```python
from gv_significance import poisson_gaussian

poisson_gaussian.significance(n=100, b=90, sigma=2.4)
```
