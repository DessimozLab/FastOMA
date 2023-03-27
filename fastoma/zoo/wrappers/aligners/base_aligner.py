import itertools
from abc import ABCMeta, abstractmethod
from enum import Enum
from Bio import AlignIO, SeqIO


from ...seq_utils import is_dna, identify_input, AlignmentInput
from .. import WrapperError




DataType = Enum('DataType', 'DNA PROTEIN UNKNOWN')


class Aligner(object):
    """
    Base class for wrappers of Multiple Sequence Aligner software

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a base implementation to be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = ConcreteAligner(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result

    """
    __metaclass__ = ABCMeta

    def __init__(self, input_, datatype=DataType.UNKNOWN, binary=None):
        """
        Should work the same whether you're working with a Biopython object or a file
        but the implementation differs, e.g. a Biopython object will need
        to be written temporarily to disk for the Aligner to work on it.

        :param input_: can be either a filename or a biopython multiple
            sequence alignment (a collection of :class:`Bio.SeqRecord.SeqRecord`)

        :param binary: is the alignment's executable file, or None. If set to
            None, it is assumed to be found in the PATH.

        :param datatype: means is it DNA or protein?
        """
        self.input_type = identify_input(input_)  # Figure out what it is - file or object

        if isinstance(datatype, str):
            try:
                datatype = getattr(DataType, datatype.upper())
            except AttributeError:
                raise ValueError("\"{}\" is an invalid datatype for an Aligner".format(datatype))
        if datatype == DataType.UNKNOWN:
            self.datatype = guess_datatype(input_, from_filename=self.input_type == AlignmentInput.FILENAME)
            if self.input_type == AlignmentInput.OBJECT:
                dup, input_ = itertools.tee(input_)
                self.datatype = guess_datatype(dup, False)
            else:
                self.datatype = guess_datatype(input_, True)
        else:
            self.datatype = datatype

        self.input = input_  # store it
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
        pass

import logging
logger = logging.getLogger()


def guess_datatype(alignment, from_filename=False):
    logger.warning("Guessing is not recommended - specify the sequence type with option datatype={DNA, PROTEIN}, be more confident")
    if from_filename:
        try:
            alignment = SeqIO.parse(alignment, 'fasta')
        except:
            alignment = SeqIO.parse(alignment, 'phylip-relaxed')
    return DataType.DNA if is_dna(alignment) else DataType.PROTEIN


# TODO: Break the identify_input function into two parts - one to work out the datatype, one to work out whether
# this is a file or an object
