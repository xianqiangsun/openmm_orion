# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - set up protein:ligand complex and emit OpenMM System

## Workfloes

* `floes/example_floe.py` - report available OpenMM Platforms
* `floes/openmm_complex_setup.py` - set up protein:ligand complexes as OpenMM System object

## Example

Test setup of OpenMM complex:
```python
python floes/openmm_complex_setup.py --complex_setup-receptor OpenMMCubes/tests/input/receptor.pdbfixer.pdb --ifs OpenMMCubes/tests/input/ligand.tripos.mol2
```
