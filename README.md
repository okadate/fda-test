# FDA

An example implementation of the finite difference-based approximation (FDA) for tangent linear and adjoint models.

This repository contains the code for a simple nutrient-phytoplankton-zooplankton (NPZ) model, including its nonlinear code and FDA-based tangent linear and adjoint code (similar to the pseudo-code in Mattern and Edwards (2018)). It also includes a few examples to test the consistency between models.

## Getting Started

All that is needed is a Python 3 installation with [numpy](https://www.numpy.org/) and (optionally) [matplotlib](http://www.matplotlib.org/).

## Running the tests

To see if it is all working, try running the NPZ nonlinear and tangent linear model (requires matplotlib):
```
python3 NPZmodel.py
```
or perform a series of dot-product tests:
```
python3 dotproduct_test.py
```

## Authors

* Jann Paul Mattern [jpmattern](https://github.com/jpmattern)

## References

Mattern, J. P. and C. A. Edwards (2018), *Journal of Geophysical Research -- Oceans*, A simple difference quotient-based approximation for biogeochemical tangent linear and adjoint models, [doi:10.1029/2018JC014283](https://dx.doi.org/10.1029/2018JC014283).

## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.
