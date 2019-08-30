# The Maryland Equilibrium Solver (MEQ)

This is a fully-featured Free Boundary Grad-Shafranov solver, based on the 
Hybridized Discontinuous Galerkin algorithm of Sanchez-Vizuet & Solano described [here](https://arxiv.org/pdf/1712.04148).

## Getting Started

To run

### Prerequisites

To compile and install MEQ you will need a system with the following

 - A C++17 compliant C++ compiler.
 - An installation of the SUNDIALS library.
 - The Boost C++ Template Library and Boost Unit Testing Framework library.
 - The gmsh library for mesh generation.

Precise dependencies have not been exhaustively tested.

### Installing

To begin, we will clone a copy of the MEQ-patched MFEM library.
```
git clone  <badger badger badger>
```

Then we compile MEQ itself
```
make BODGER
```

As an example, we compute the vacuum field from a pair of Helmholtz coils
```
./meq -c examples/Helmholtz.conf
```

## Running the tests

Various test suites are provided with MEQ.

### Automated Tests

The automated test suite contains both unit tests and functional tests. 
These test independent subsystems of the MEQ code.
To run the test suite

```
make test
```

If the test suite fails, please file a bug.

### Example Configurations

These are fully-fledged example runs of MEQ. As such, they may take considerable computational resoruces to run. 
Each comes with a pre-computed solution against which the output is compared.

To run an example and check its output, for example the 'Helmholtz' example:
```
cd examples/
make Helmholtz
```

Or you can separately run MEQ for the example
```
meq -c examples/Helmholtz.conf
```
and check the output later with
```
cd examples
make Helmholtz-check
```

## Parallelism

Currently none.

## Built With

* [Boost](http://boost.org) - C++ Template library that radically extends the STL
* [MFEM](http://github.com/mfem) - The Finite Element Framework on which MEQ is based
* [Sundials](http://computing.llnl.gov/projects/sundials) - ODE Framework, used for Anderson-accelerated Picard iteration
* [TOML11](http://github.com/toruniina/toml11) - For parsing configuration files written in [TOML](https://github.com/toml-lang/toml)

## Contributing

Contributions to this project are welcome. However, as we are currently in active development, please email the authors if you wish to get involved.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/ianabel/meq/tags). 

## Authors

* **Ian Abel** - *Initial Work* - [Ian Abel at UMD](https://ireap.umd.edu/faculty/abel)

For full copyright attribution, see the [COPYRIGHT](COPYRIGHT) file.
For a summary of contributors, see the [contributors](http://github.com/ianabel/meq/contributors) page.

## License

This project is licensed under the 3-Clause BSD Licence - see the [LICENSE.md](LICENSE.md) file for details

This project links to a patched version of the MFEM library, which licenced under the GNU Lesser General Public Licence version 2.1. 
The edited library is released under the same terms as MFEM itself.

## Acknowledgments

* Many thanks go to the [MFEM](http://girhub.com/mfem) team at [LLNL](http://computing.llnl.gov) for making a great C++ FEM library
* The HDG branch of MFEM used in this project was written by T. Horvath and collaborators at the University of Waterloo
* The HDG Algorithm for solving the Grad-Shafranov equation was developed by T. Sanchez-Vizuet at New York University and M. Solano at Universidad de Concepcion, Chile.

