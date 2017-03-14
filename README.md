# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `LigPrepCubes/` - Cubes for preparing molecules
  * `ChargeMCMol` - Assigns partial charges and generates multi-conf molecules with OMEGA
  * `FREDDocking` - Dock MCMols using FRED to a prepared receptor
  * `SMIRNOFFParameterization` - Parametrizes molecule with SMIRNOFF forcefield
  * `GAFFParameterization` - Parametrizes molecule with GAFF forcefield
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - set up protein:ligand complex and emit OpenMM System
  * `OpenMMSimulation` - run OpenMM simulation.

## Workfloes

* `floes/openmm_complex_setup.py` - Setup the protein:ligand complex from **PDBs**.
* `floes/openmm_setup_md.py` - Setup the protein:ligand complex and run MD simulation from **PDBs**
* `floes/openmm_md.py` - Run MD simulation from a prepared *complex.oeb.gz*
* `floes/smiles_ligprep_smirnoff.py` - Parse smiles, dock and parameterize with SMIRNOFF.
* `floes/smiles_ligprep_gaff.py` - Parse smiles, dock, and parameterize with GAFF.
* `floes/smiles_complex_setup.py` - Parse smiles and setup the protein:ligand complex.

## Local Installation
```bash
git clone git@github.com:openeye-private/openmm_orion.git
cd openmm_orion

#Create a new local conda environment and install dependencies
conda create -c omnia -c omnia/label/dev -c mobleylab -n dev python=3.5 openmm==7.1.0rc1 openmoltools==0.7.5 ambermini==16.16.0 smarty==0.1.4 parmed==2.7.1 pdbfixer==1.4 smirff99frosst==1.0.4
source activate dev

#Install the OpenEye-floe package and toolkits
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
pip install OpenEye-floe-0.2.141.tar.gz

#Install the main OpenMM Orion Floes
python setup.py develop

# Run the tests.
python LigPrepCubes/tests/test_cubes.py
python OpenMMCubes/tests/test_cubes.py
```

## Examples
### Starting from PDBs
#### [SetupOpenMMComplex]: Setup the protein:ligand complex.
*Assuming input ligand PDB is in docked position*
Setup the protein-ligand complex.
```
python floes/openmm_complex_setup.py --ligand examples/data/toluene.pdb --protein examples/data/T4-protein.pdb

# Available options
python floes/openmm_complex_setup.py --ligand examples/data/toluene.pdb --protein examples/data/T4-protein.pdb --ffxml smirff99Frosst.ffxml --protein_ff amber99sbildn.xml --solvent_ff tip3p.xml --pH 7.0 --salt_conc 10
```

#### [SetupOpenMMSimulation]: Setup and simulate the protein:ligand complex.
Does the same above and run 10,000 MD steps.
```
python floes/openmm_setup_md.py --ligand examples/data/toluene.pdb --protein examples/data/T4-protein.pdb --steps 10000
```

### Using the prepared complex.oeb.gz files
#### [RunOpenMMSimulation] Run a MD Simulation given a complex.oeb.gz
If you want to just run MD from the prepared complex:
```
python floes/openmm_md.py --complex examples/data/JF6_1-complex.oeb.gz --steps 10000
```


### Starting from SMILES
#### [SmilesLigPrep]: Generate conformers, dock and parameterize molecules
Starting from a SMILE string, generate conformers with OMEGA, dock to a
prepared receptor using FRED, and parameterize the molecules. Writes out the
SMIRFF molecules in their docked poses.
```
python floes/smiles_ligprep.py --ligand examples/data/JF6_1.ism --receptor examples/data/epox_hydrolase_receptor.oeb.gz
```

#### [SmilesComplexPrep]: Setup the protein:ligand complexes
Does the same as above and then prepares the complex from a PDB of the receptor.
Writes out the protein:ligand complex.
```
python floes/smiles_complex_setup.py --ligand examples/data/JF6_1.ism --receptor examples/data/epox_hydrolase_receptor.oeb.gz --protein examples/data/epox_hydrolase_apo-protein.pdb
```
