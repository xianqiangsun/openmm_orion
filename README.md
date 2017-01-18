# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - set up protein:ligand complex and emit OpenMM System
  * `OpenMMSimulation` - run OpenMM simulation.

## Workfloes

* `floes/example_floe.py` - report available OpenMM Platforms
* `floes/openmm_complex_setup.py` - set up protein:ligand complexes as OpenMM System object
* `floes/openmm_md.py` - run short md simulation of protein:ligand complex
* `floes/openmm_continue.py` - restart md simulation from saved State

## Example

Test setup of OpenMM complex:
```bash
python setup.py develop or pip install -e ./

# Setup protein-ligand complex
python floes/openmm_complex_setup.py --protein OpenMMCubes/tests/input/T4-protein.pdb --ligand OpenMMCubes/tests/input/toluene.pdb --molecule_forcefield OpenMMCubes/tests/input/forcefield/smirff99Frosst.ffxml

# Run short MD simulation
python floes/openmm_md.py --complex_pdb combined_structure.pdb --system system.xml.xz

# Restart MD simulation from saved State
python floes/openmm_continue.py --complex_pdb combined_structure.pdb --system system.xml.xz --state state.xml
```
