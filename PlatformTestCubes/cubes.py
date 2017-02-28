from floe.api import(
    OEMolComputeCube, parameter, BinaryOutputPort, SourceCube
)

from PlatformTestCubes import testInstallation
from PlatformTestCubes.benchmarking import run_platform_benchmarks


# For parallel, import and inherit from ParallelOEMolComputeCube
class PlatformTestCube(OEMolComputeCube):
    """
    Runs a copy of OpenMM's simtk.installation script

    Checks available OpenMM platforms
    """
    title = "OpenMM Platform Check"
    description = """
    *OpenMM Platform Check*
    Checks to see which OpenMM Platforms are available amongst CPU, OpenCL, and CUDA
    """
    classification = [
        ["OpenMM", "PlatformCheck"],
    ]
    tags = [tag for lists in classification for tag in lists]

    # the string in the first argument has to be the same as the name of the object
    # note that right now, the pdb file has to reside in the same directory as this file
    pdbFileName = parameter.StringParameter(
        "pdbFileName",
        title="pdb File Name",
        description="name of pdb file to use in OpenMM Platform Check",
        default='test.pdb',
    )

    def process(self, mol, port):
        output = testInstallation.run_tests(self.args.pdbFileName)
        self.log.info(output)
        with open('openmmPlatformCheck.txt', 'w') as ofs:
            ofs.write(output)
        self.emit(mol)


class BenchmarkCube(SourceCube):

    title = "OpenMM BenchmarkCube"

    description = """
        Cube that performs a benchmark in a floe on all of the different
        platforms that are available to OpenMM and outputs the byte string
        resulting from the benchmarks to be handed to a FileOutput Cube
    """

    tags = [["OpenMM", "Benchmarking"]]

    success = BinaryOutputPort("success")
    failure = BinaryOutputPort("failure")

    cutoff = parameter.DecimalParameter("cutoff", default=0.9)
    seconds = parameter.IntegerParameter("seconds", default=60)
    polarization = parameter.StringParameter(
        "polarization",
        default="mutual",
        choices=["direct", "extrapolated", "mutual"]
    )
    epsilon = parameter.DecimalParameter("epsilon", default=1e-5)
    heavy = parameter.BooleanParameter("heavy", default=False)
    precision = parameter.StringParameter(
        "precision",
        default="single",
        choices=["single", "mixed", "double"]
    )

    def __iter__(self):
        try:
            run_platform_benchmarks(self.args)
            yield b"Success"
        except Exception as e:
            self.failure.emit(str(e).encode("utf-8"))
