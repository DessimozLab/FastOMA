import abc
import re
from ._wrappers import logger


class FastaHeaderTransformer(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def transform(self, header):
        return header


class NoOpFastaHeaderTransformer(FastaHeaderTransformer):
    def transform(self, header):
        return header


class ExtractUniProtAccessionFastaHeaderTransformer(FastaHeaderTransformer):
    def __init__(self):
        self._up_re = re.compile(r"[sptr]{2}\|(?P<acc>[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\|.*")

    def transform(self, header):
        m = self._up_re.match(header)
        if m:
            return m.group('acc')
        logger.warning("cannot extract uniprot accession from header: %s", header)
        return header


def header_transformer(name):
    if name.lower() == "noop":
        return NoOpFastaHeaderTransformer()
    elif name.lower() == 'uniprot':
        return ExtractUniProtAccessionFastaHeaderTransformer()
