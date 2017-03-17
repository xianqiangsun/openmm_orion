# OpenMM cubes and workfloes for Orion

## Cube sets

* `PlatformTestCubes/` - simple example cube for testing available OpenMM PlatformTestCubes
* `LigPrepCubes/` - Cubes for preparing molecules
  * `ChargeMCMol` - Assigns partial charges and generates multi-conf molecules with OMEGA
  * `FREDDocking` - Dock MCMols using FRED to a prepared receptor
  * `LigandParameterization` - Parametrize molecules with either GAFF/GAFF2/SMIRNOFF forcefields
* `OpenMMCubes/` - OpenMM utility cubes
  * `OpenMMComplexSetup` - Generate a solvated protein:ligand complex.
  * `OpenMMSimulation` - Runs and OpenMM MD simulation. Minimizes or restarts from saved State.
* `YankCubes/` - YANK cubes
  * `YankHydrationCube` - YANK hydration free energy calculations

## Workfloes

* `floes/openmm_complex_setup.py` - Setup the protein:ligand complex from **PDBs**.
* `floes/openmm_complex_min.py` - Setup the protein:ligand complex, minimize and run short MD.
* `floes/openmm_md.py` - Run MD simulation from a prepared *complex.oeb.gz*
* `floes/smiles_ligprep.py` - Parse SMILES, assign charges, docks, and parameterizes the molecules with GAFF/GAFF2/SMIRNOFF forcefield parameters
* `floes/smiles_complex_setup.py` - Parse SMILES, prepare ligands (above), and generate the protein:ligand complex.

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
## Detailed Walkthrough: From SMILES to Simulation.
### [FLOE] SmilesLigPrep: Preparing a molecule dataset from SMILES strings.
```bash
python floes/smiles_ligprep.py --ligand examples/data/JF6_1.ism --receptor examples/data/epox_hydrolase_receptor.oeb.gz --molecule_forcefield GAFF2 \
--ofs-data_out JF6_1-gaff2.oeb.gz --fail-data_out fail.oeb.gz
```

The first cube uses OpenEye's native `OEMolIStreamCube` to parse a file containing SMILES strings and emit each `oechem.OEMol`.

In the following cube, `LigPrepCubes.cubes.ChargeMCMol` will call functions from the package [openmoltools](https://github.com/choderalab/openmoltools/blob/master/openmoltools/openeye.py) in order to check aromaticity, add explicit hydrogens, and obtain the IUPAC name of the molecule. Then, the cube generates multiple conformers with OMEGA and assign charges using the new `oequacpac.OEAssignCharges` function. By default, the charging engine used is set to `OEAM1BCCCharges`. Before emitting, the IUPAC name and molecule title are stored as SDData using the tags `IUPAC` and `IDTag`, respectively.
###### The IDTag/molecule title is drawn from the string next to the SMILES string from the input file, (i.e. `c1ccc(cc1)c2c(non2)N JF6_1`) or it is randomly generated if not present.

Next, the charged multi-conformer molecules are docked with `LigPrepCubes.cubes.FREDDocking` which uses the FRED or ChemGauss4 docking engine. The pose score and docking method is attached as SDData before emitting.

After docking, the molecules are passed into `LigPrepCubes.cubes.LigandParameterization` where a forcefield parameterized `parmed.Structure` is generated for the respective molecule and attached as generic data. The `LigandParameterization` cube supports parameterizing the molecule with GAFF, GAFF2, or the [SMIRNOFF](https://github.com/open-forcefield-group/smirff99Frosst) forcefields. By default, this cube will use GAFF2.

Molecules (with the `parmed.Structure` attached) are then passed to OpenEye's `OEMolOStreamCube` and written into a dataset file (`oeb`).

### [FLOE] MinimizeComplex: Generate and minimize the protein:ligand complex.
```bash
python floes/openmm_complex_min.py --ligand JF6_1-gaff2.oeb.gz --protein examples/data/epox_hydrolase_apo-protein.pdb \
--ofs-data_out JF6_1-min.oeb.gz --fail-data_out fail.oeb.gz
```


## Other Floe Examples:
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
