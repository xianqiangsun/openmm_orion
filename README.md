# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - set up protein:ligand complex and emit OpenMM System
  * `OpenMMSimulation` - run OpenMM simulation.

## Workfloes

* `floes/example_floe.py` - report available OpenMM Platforms
* `floes/openmm_setup_md.py` - setup the protein:ligand complex and run 1000 MD steps.

## Example

Test setup of OpenMM complex:
```bash
python setup.py develop or pip install -e ./

# Setup protein-ligand complex and Run short MD simulation
python floes/openmm_setup_md.py --protein OpenMMCubes/tests/input/T4-protein.pdb --ligand OpenMMCubes/tests/input/smirff_mol.oeb.gz

```
