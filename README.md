# The Maryland Equilibrium Solver (MEQ)

This is a fully-featured Free Boundary Grad-Shafranov solver, based on the 
Hybridized Discontinuous Galerkin algorithm of Sanchez-Vizuet & Solano described [here](https://arxiv.org/pdf/1712.04148).

## Getting Started

To run

### Prerequisites

To compile and install MEQ you will need a system with the following

 - A C++14 compliant C++ compiler
 - An installation of the SUNDIALS library (version 5.0 or newer)
 - 

Precise compiler compatability is currently under test. However, the code is known to compile and
work with g++ and Clang.

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [MFEM](http://github.com/mfem) - The Finite Element Framework on which MEQ is based
* [Sundials](http://computing.llnl.gov/projects/sundials) - ODE Framework, used for Anderson-accelerated Picard iteration
* [TOML11](http://github.com/toruniina/toml11) - For parsing configuration files written in [TOML](https://github.com/toml-lang/toml)

## Contributing

Contributions to this project are welcome. 

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/ianabel/meq/tags). 

## Authors

* **Ian Abel** - *Initial Work* - [Ian Abel at UMD](https://ireap.umd.edu/faculty/abel)

For full copyright attribution, see the [COPYRIGHT](COPYRIGHT) file.
For a summary of contributors, see the [contributors](http://github.com/ianabel/meq/contributors) page.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

This project links to the MFEM library, which licenced under the GNU Lesser General Public Licence version 2.1. 

## Acknowledgments

* Many thanks go to the [MFEM](http://girhub.com/mfem) team at [LLNL](http://computing.llnl.gov) for making a great C++ FEM library
* The HDG branch of MFEM used in this project was written by T. Horvath and collaborators at the University of Waterloo
* The HDG Algorithm for solving the Grad-Shafranov equation was developed by T. Sanchez-Vizuet at New York University and M. Solano at Universidad de Concepcion, Chile.

