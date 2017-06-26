import unittest
import pytest
from floe.test import CubeTestRunner
from openeye import oechem
import OpenMMCubes.utils as utils
from OpenMMCubes.simtools import restraints
from OpenMMCubes.cubes import OpenMMminimizeCube, OpenMMnvtCube, OpenMMnptCube
from simtk import unit, openmm
from simtk.openmm import app
import pickle


class MinimizationCubeTester(unittest.TestCase):
    """
    Test the OpenMM Minimization cube
    Example inputs from `openmm_orion/examples/data`
    """

    def calculate_eng(self, oe_mol):
        # Extract starting MD data from OEMol
        mdData = utils.MDData(oe_mol)
        topology = mdData.topology
        positions = mdData.positions
        structure = mdData.structure
        box = mdData.box

        # OpenMM system
        system = structure.createSystem(nonbondedMethod=eval("app.%s" % self.cube.args.nonbondedMethod),
                                        nonbondedCutoff=self.cube.args.nonbondedCutoff * unit.angstroms,
                                        constraints=eval("app.%s" % self.cube.args.constraints))
        # OpenMM Integrator
        integrator = openmm.LangevinIntegrator(self.cube.args.temperature * unit.kelvin,
                                               1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
        # Set Simulation
        simulation = app.Simulation(topology, system, integrator)

        # Set Positions
        simulation.context.setPositions(positions)

        # Set Box dimensions
        simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

        # Collect the OpenMM state energy info
        state = simulation.context.getState(getEnergy=True)

        # Potential Energy
        peng = state.getPotentialEnergy()

        return peng

    def setUp(self):
        self.cube = OpenMMminimizeCube('minComplex')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def _test_success(self):
        print('Testing cube:', self.cube.name)
        # Complex file name
        complex_fname = utils.get_data_filename('examples', 'data/pbace_lcat13a_complex.oeb.gz')

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(complex_fname) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Calculate starting potential energy
        eng_i = self.calculate_eng(mol)

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()
        # Calculate final potential energy
        eng_f = self.calculate_eng(outmol)

        self.assertLess(eng_f, eng_i)

    @pytest.mark.slow
    def test_success(self):
        self.cube.args.steps = 100000
        self._test_success()

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


class NVTCubeTester(unittest.TestCase):
    """
    Test the OpenMM NVT cube
    Example inputs from `openmm_orion/examples/data`
    """

    def calculate_VT(self, oe_mol):
        # Extract starting MD data from OEMol
        mdData = utils.MDData(oe_mol)
        structure = mdData.structure
        topology = mdData.topology
        positions = mdData.positions
        velocities = mdData.velocities
        box = mdData.box

        volume = box[0][0] * box[1][1] * box[2][2]

        # OpenMM system
        system = structure.createSystem(nonbondedMethod=eval("app.%s" % self.cube.args.nonbondedMethod),
                                        nonbondedCutoff=self.cube.args.nonbondedCutoff * unit.angstroms,
                                        constraints=eval("app.%s" % self.cube.args.constraints))
        # OpenMM Integrator
        integrator = openmm.LangevinIntegrator(self.cube.args.temperature * unit.kelvin,
                                               1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
        # Set Simulation
        simulation = app.Simulation(topology, system, integrator)

        # Set Positions
        simulation.context.setPositions(positions)
        # Set Velocities
        simulation.context.setVelocities(velocities)

        # Set Box dimensions
        simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

        # Collect the OpenMM state energy info
        state = simulation.context.getState(getEnergy=True)

        # Kinetic Energy
        keng = state.getKineticEnergy()

        # Calculate system degrees of freedom:
        dof = 0
        for i in range(system.getNumParticles()):
            if system.getParticleMass(i) > 0 * unit.dalton:
                dof += 3

        dof -= system.getNumConstraints()

        if any(type(system.getForce(i)) == openmm.CMMotionRemover for i in range(system.getNumForces())):
            dof -= 3

        # Calculate the temperature from the equipartition theorem
        temperature = ((2*keng)/(dof*unit.MOLAR_GAS_CONSTANT_R))

        return volume, temperature

    def setUp(self):
        self.cube = OpenMMnvtCube('NVT')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def _test_success(self):
        print('Testing cube:', self.cube.name)
        # Complex file name
        complex_fname = utils.get_data_filename('examples', 'data/pP38_lp38a_2x_complex.oeb.gz')

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(complex_fname) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()

        # Calculate final volume and temperature
        vol_f, temp_f = self.calculate_VT(outmol)

        # Check 3*std volume
        # Average volume and its standard deviation (in nm^3) measured along
        # one 10ns run for the selected system
        avg_volume = 683.492706 * (unit.nanometers ** 3)
        std_volume = 0.000001

        self.assertAlmostEqual(avg_volume/(unit.nanometers ** 3),
                               vol_f.in_units_of(unit.nanometers ** 3) / (unit.nanometers ** 3),
                               delta=3*std_volume)

        # Check temperature
        # Average temperature and its standard deviation (in K) measured along
        # one 10ns run for the selected system
        avg_temperature = 300.0 * unit.kelvin
        std_temperature = 0.9
        self.assertAlmostEqual(avg_temperature/unit.kelvin,
                               temp_f.in_units_of(unit.kelvin)/unit.kelvin,
                               delta=3*std_temperature)

    @pytest.mark.slow
    def test_success(self):
        self.cube.args.time = 50.0  # in picoseconds
        self._test_success()

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


class NPTCubeTester(unittest.TestCase):
    """
    Test the OpenMM NPT cube
    Example inputs from `openmm_orion/examples/data`
    """
    def calculate_VT(self, oe_mol):
        # Extract starting MD data from OEMol
        mdData = utils.MDData(oe_mol)
        structure = mdData.structure
        topology = mdData.topology
        positions = mdData.positions
        velocities = mdData.velocities
        box = mdData.box

        volume = box[0][0] * box[1][1] * box[2][2]

        # OpenMM system
        system = structure.createSystem(nonbondedMethod=eval("app.%s" % self.cube.args.nonbondedMethod),
                                        nonbondedCutoff=self.cube.args.nonbondedCutoff * unit.angstroms,
                                        constraints=eval("app.%s" % self.cube.args.constraints))
        # OpenMM Integrator
        integrator = openmm.LangevinIntegrator(self.cube.args.temperature * unit.kelvin,
                                               1.0 / unit.picoseconds, 0.002 * unit.picoseconds)

        # Add Force Barostat to the system
        system.addForce(openmm.MonteCarloBarostat(self.cube.args.pressure * unit.atmospheres,
                                                  self.cube.args.temperature*unit.kelvin, 25))

        # Set Simulation
        simulation = app.Simulation(topology, system, integrator)

        # Set Positions
        simulation.context.setPositions(positions)

        # Set Velocities
        simulation.context.setVelocities(velocities)

        # Set Box dimensions
        simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

        # Collect the OpenMM state energy info
        state = simulation.context.getState(getEnergy=True)

        # Kinetic Energy
        keng = state.getKineticEnergy()

        # Calculate system degrees of freedom:
        dof = 0
        for i in range(system.getNumParticles()):
            if system.getParticleMass(i) > 0 * unit.dalton:
                dof += 3

        dof -= system.getNumConstraints()

        if any(type(system.getForce(i)) == openmm.CMMotionRemover for i in range(system.getNumForces())):
            dof -= 3

        # Calculate the temperature from the equipartition theorem
        temperature = ((2*keng)/(dof*unit.MOLAR_GAS_CONSTANT_R))

        return volume, temperature

    def setUp(self):
        self.cube = OpenMMnptCube('NPT')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def _test_success(self):
        print('Testing cube:', self.cube.name)

        # Complex file name
        complex_fname = utils.get_data_filename('examples', 'data/pP38_lp38a_2x_complex.oeb.gz')
        # Read OEMol molecule
        mol = oechem.OEMol()
        with oechem.oemolistream(complex_fname) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()

        # Calculate final volume and temperature
        vol_f, temp_f = self.calculate_VT(outmol)

        # Check 3*std volume
        # Average volume and its standard deviation (in nm^3) measured along
        # one 10ns run for the selected system
        avg_volume = 683.0 * (unit.nanometers**3)
        std_volume = 1.0

        self.assertAlmostEqual(avg_volume/(unit.nanometers**3),
                               vol_f.in_units_of(unit.nanometers**3)/(unit.nanometers**3),
                               delta=3*std_volume)

        # Check 3*std temperature
        # Average temperature and its standard deviation (in K) measured along
        # one 10ns run for the selected system
        avg_temperature = 300.0 * unit.kelvin
        std_temperature = 1.0

        self.assertAlmostEqual(avg_temperature/unit.kelvin,
                               temp_f.in_units_of(unit.kelvin)/unit.kelvin,
                               delta=3*std_temperature)

    @pytest.mark.slow
    def test_success(self):
        self.cube.args.time = 50.0  # in picoseconds
        self._test_success()

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


# Check Restraints applications
def test_restraints():

    lig_fname = utils.get_data_filename('examples', 'data/pP38_lp38a_2x_complex.oeb.gz')
    # Read OEMol molecule
    mol = oechem.OEMol()
    with oechem.oemolistream(lig_fname) as ifs:
        oechem.OEReadMolecule(ifs, mol)

    res_dic = {}
        
    mask = 'protein'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'ligand'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'water'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set
    
    mask = 'ions'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'cofactors'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'ca_protein'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'protein or ligand'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'noh (protein or ligand)'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'ca_protein or (noh ligand)'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    mask = 'resid A:4 A:5 A:6 A:7'
    ind_set = restraints(mol, mask=mask)
    res_dic[mask] = ind_set

    dic_fname = utils.get_data_filename('examples', 'data/restraint_test_P38_lp38a_2x.pickle')

    file = open(dic_fname, 'rb')
    res_dic_loaded = pickle.load(file)

    for k in res_dic_loaded:
        if res_dic[k] == res_dic_loaded[k]:
            pass
        else:
            raise ValueError("Restraints checking Errors on mask: {}".format(k))


if __name__ == "__main__":
        unittest.main()
