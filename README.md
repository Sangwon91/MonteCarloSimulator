# Monte Carlo Simulator
Last update: 2014

Molecular Monte Carlo simulator for the free energy calculations of fluid and solid.

# Features
1. NVT simulations.
2. NPT simulations.
3. Expanded ensemble method simulations.
4. Rigid/flexible molecule simulations.
5. RDF calculations.

# Installation
```bash
$ cd src
$ make
```

# Example

### Inputs
* `system.txt`: Overall simulation options.
* `config.txt`: Initial configurations.
* `molecule.txt`: Molecule information (force field, connectivity).

### Command for simulation
```bash
$ ./mc.x
```
