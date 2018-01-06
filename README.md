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

## Ideal case
Compute the significance of a measurement of 100 counts with an expected 
background of 90 in the ideal case of no uncertainty on the background:

```python
from gv_significance import ideal_case

ideal_case.significance(n=100, b=90)

```

## Poisson observation, Poisson background with or without systematics

Compute the significance of a measurement of 100 counts when a side measurement
returned a background of 900 counts with an exposure 10x the one of the observation
on-source (`alpha=0.1`):

```python
from gv_significance import poisson_poisson

poisson_poisson.significance(n=100, b=900, alpha=0.1)

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

poisson_gaussian.significance(n=100, b=90, sigma=0.4)
```