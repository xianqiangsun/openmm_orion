# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - set up protein:ligand complex and emit OpenMM System
  * `OpenMMSimulation` - run OpenMM simulation.

## Workfloes

* `floes/example_floe.py` - report available OpenMM Platforms
* `floes/openmm_complex_setup.py` - set up protein:ligand complexes as OpenMM System object

## Example

Test setup of OpenMM complex:
```bash
python setup.py develop or pip install -e ./
python floes/openmm_complex_setup.py --protein OpenMMCubes/tests/input/T4-protein.pdb --ligand OpenMMCubes/tests/input/toluene.pdb --molecule_forcefield OpenMMCubes/tests/input/forcefield/smirff99Frosst.ffxml
```
