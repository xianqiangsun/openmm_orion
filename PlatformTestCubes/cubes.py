import time
from floe.api import OEMolComputeCube, parameter
from PlatformTestCubes import testInstallation


# For parallel, import and inherit from ParallelOEMolComputeCube
class PlatformTestCube(OEMolComputeCube):
    """
    Runs a copy of OpenMM's simtk.installation script

    Checks available OpenMM platforms
    """
    title = "Cube Title for view in UI"
    description = """
    *Longform Description*
    Shown in the UI for floe builders. Should explain cube and params
    """
    classification = [
        ["Testing", "OpenMM"],
        ["Testing", "PlatformCheck"],
    ]
    tags = [tag for lists in classification for tag in lists]

    # the string in the first argument has to be the same as the name of the object
    # note that right now, the pdb file has to reside in the same directory as this file
    pdbFileName = parameter.StringParameter("pdbFileName",
            title="pdb File Name",
            description="name of pdb file to use in run_test",
            default='test.pdb',
            )

    def process(self, mol, port):
        output = testInstallation.run_tests( self.args.pdbFileName)
        self.log.info( output)
        ofs = open('openmmPlatformCheck.txt','w')
        ofs.write( output)
        ofs.close()
        self.emit(mol)

