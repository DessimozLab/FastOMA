
import logging
import os
import time

from Bio import SeqIO, AlignIO
from pyparsing import ParseException
import tempfile

from .base_trimmer import MSATrimmer
from ..aligners.base_aligner import AlignmentInput
from ...wrappers import WrapperError
from ..abstract_cli import AbstractCLI
from ..options import OptionSet, StringOption, IntegerOption, FlagOption, FloatOption
from ...seq_utils.utils import iter_seqrecs_from_any
from ...file_utils import TempFile, TempDir

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class TrimAlCLI(AbstractCLI):
    """
       TrimAl low-level command line interface

       :Example:

       ::

           trimal_cli = TrimAlCLI()
           process = trimal_cli(cmd='trimal args...')
           stdout = trimal_cli.get_stdout()
       """
    @property
    def _default_exe(self):
        return 'trimal'

    @property
    def _hyphen_policy(self):
        return 1


def set_representative_options(trimmer, size:int):
    """
    get most representative sequences for msa of size `size`
    """
    trimmer.options = get_default_options()
    trimmer.options['-cluster'].set_value(size)


def get_phylogenetic_ml_options(trimmer):
    """Use an heuristic to decide the optimal method for trimming the alignment.
    (see User Guide for details).
    Optimized for Maximum Likelihood phylogenetic tree reconstruction"""
    trimmer.options = get_default_options()
    trimmer.options['-automated1'].set_value(True)


def get_phylogenetic_neighbour_joining(trimmer):
    """Use automated selection on "strictplus" mode. (see User Guide).
    Optimized for Neighbour Joining phylogenetic tree reconstruction. """
    trimmer.options = get_default_options()
    trimmer.options['-strictplus'].set_value(True)


class TrimAl(MSATrimmer):
    """
    Convenient wrapper for TrimAl MAS Trimmer

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a basic implementation that can be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = TrimAl(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result


    .. note:: There exists an ipython notebook on how to work with wrappers,
         including dealing with non-default parameters.
    """

    def __init__(self, input_, *args, **kwargs):
        super(TrimAl, self).__init__(input_, *args, **kwargs)
        self.options = get_default_options()

    def __call__(self, *args, **kwargs):
        """
        Anything to do with calling Mafft should go here.
        If any extra arguments need to be passed they can
        be specified (listed as *args and **kwargs for now).
        """
        start = time.time()  # time the execution

        with tempfile.TemporaryDirectory() as tmpdir:
            outfn = os.path.join(tmpdir, "output.trimal.fasta")
            if self.input_type == AlignmentInput.OBJECT:  # different operation depending on what it is
                with open(os.path.join(tmpdir, "sequences.fa"), 'wt') as fh:
                    SeqIO.write(iter_seqrecs_from_any(self.input), fh, 'fasta')
                output, error = self._call(fh.name, outfn, *args, **kwargs)
            else:
                output, error = self._call(self.input, outfn, *args, **kwargs)
            self.result = self._read_result(outfn)
        logger.debug('Output of TrimAl: stdout={}; stderr={}'.format(output, error))
        if len(output) == 0 and len(error) > 0:
            raise WrapperError('TrimAl did not compute any alignments. StdErr: {}'.format(error))
        self.stdout = output
        self.stderr = error

        end = time.time()
        self.elapsed_time = end - start
        return self.result
        # End call

    # Any other accessory methods
    def _call(self, filename, outfn, *args, **kwargs):
        """
        Call underlying low level _Mafft wrapper.
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        if self.command() == "-fasta":
            logger.warning("no trimming options specified. will use '-automated1' insted.")
            self.options['-automated1'].set_value(True)

        self.cli('-in {} -out {} {}'.format(filename, outfn, self.command()),
                 wait=True)
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, output):
        """
        Read back the result.
        """
        return AlignIO.read(output, 'fasta')

    def _init_cli(self, binary):
        return TrimAlCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # Use a heuristic selection of the automatic method based on
        # similarity statistics. (see User Guide).
        # (Optimized for Maximum Likelihood phylogenetic tree reconstruction).
        FlagOption('-automated1', False, active=False),

        # Only columns out of internal boundaries (first and last column without gaps) are
        # candidated to be trimmed depending on the applied method
        FlagOption('-terminalonly', False, active=False),

        # Use automated selection on "strictplus" mode. (see User Guide).
        # Optimized for Neighbour Joining phylogenetic tree reconstruction.
        FlagOption('-strictplus', False, active=False),

        # Use automated selection on "gappyout" mode. This method only uses information
        # based on gaps' distribution. (see User Guide).
        FlagOption('-gappyout', False, active=False),

        # Use automated selection on "strict" mode. (see User Guide).
        FlagOption('-strict', False, active=False),

        # Get the relationship between the columns in the old and new alignment.
        FlagOption('-colnumbering', False, active=False),

        # Get the complementary alignment
        FlagOption('-complementary', False, active=False),

        # Minimum column block size to be kept in the trimmed alignment. Available with
        # manual and automatic (gappyout) methods
        IntegerOption('-block', 30, active=False),

        # -gt -gapthreshold <n>    1 - (fraction of sequences with a gap allowed). Range: [0 - 1]
        FloatOption('-gapthreshold', 0.5, active=False),

        # Minimum average similarity allowed. Range: [0 - 1]
        FloatOption('-simthreshold', 0.5, active=False),

        # Minimum percentage of the positions in the original alignment to conserve. Range: [0 - 100]
        IntegerOption('-cons', 50, active=False),

        # Format of output file should be fasta
        FlagOption('-fasta', True, active=True),

        # Use a user-defined scoring matrix.  Default: BLOSUM62
        StringOption('-matrix', '', active=False),
    ])



