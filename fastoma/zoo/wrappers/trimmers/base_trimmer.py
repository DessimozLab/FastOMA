import os, types, itertools
from abc import ABCMeta, abstractmethod
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from ...seq_utils import identify_input
from ...wrappers import WrapperError

import logging
logger = logging.getLogger(__name__)



class MSATrimmer:
    """
    Base class for wrappers of msa trimming software

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a base implementation to be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = ConcreteTrimmer(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result
    """
    __metaclass__ = ABCMeta

    def __init__(self, alignment=None, binary=None):
        """
        Should work the same whether you're working with a Biopython object or a file
            but the implementation differs, e.g. a Biopython object will need
            to be written temporarily to disk for the Trimmer to work on it.

        alignment is one of 4 things:
            a filename
            a Biopython MSA
            a list of Seq objects
            anything else (throw an exception)

        binary is the alignment's executable file, or None
        """

        if alignment is not None:
            self.input_type = identify_input(alignment)  # Figure out what it is - file or object
            self.input = alignment  # store it
        else:
            self.input_type = None
            self.input = None

        self.elapsed_time = None
        self.stdout = None
        self.stderr = None
        try:
            self.cli = self._init_cli(binary)
        except IOError as err:
            raise WrapperError('Error searching for binary: {}'.format(err))
            # End setup

    @abstractmethod
    def __call__(self, *args, **kwargs):
        """
        How to call the underlying aligner
        """
        pass

    @abstractmethod
    def _init_cli(self, binary):
        """
        Set up the command-line interface to the wrapped software
        :param binary: filename of executable binary file
        :return: concrete CLI type inheriting from AbstractCLI
        """
        pass


