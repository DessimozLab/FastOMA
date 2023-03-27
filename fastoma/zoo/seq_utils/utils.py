
import os
import random
import types
from collections import defaultdict
from enum import Enum
from pathlib import Path


from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import ambiguous_dna_letters
from Bio.Seq import Seq # , UnknownSeq
from Bio.SeqRecord import SeqRecord

__all__ = ['bootstrap', 'concatenate', 'is_dna', 'BackTranslator', 'identify_input', 'AlignmentInput',
           'iter_seqrecs_from_any']


def is_dna(obj):
    """check whether a sequence is of type dna"""
    if isinstance(obj, SeqRecord):
        return is_dna(str(obj.seq))
    elif isinstance(obj, Seq):
        return is_dna(str(obj))
    elif isinstance(obj, str):
        chars = [char.upper() for char in obj if char.upper() not in '-?X']
        return set(chars).issubset(set(ambiguous_dna_letters))
    else:
        return all(is_dna(z) for z in obj)


def iter_seqrecs_from_any(obj):
    if isinstance(obj, SeqRecord):
        yield obj
    elif isinstance(obj, str):
        yield SeqRecord(Seq(obj))
    elif isinstance(obj, Seq):
        yield SeqRecord(obj)
    else:
        for sub in obj:
            yield from iter_seqrecs_from_any(sub)

def _sample_wr(population, k):
    _int = int
    _random = random.random
    n = len(population)
    return [population[_int(_random() * n)] for _ in range(k)]


def bootstrap(alignment):
    length = alignment.get_alignment_length()
    columns = [alignment[:, n] for n in range(length)]
    sample = sorted(_sample_wr(columns, length))
    sequences = [''.join(x) for x in zip(*sample)]
    ids = [seqrec.id for seqrec in alignment]
    return _assemble_msa(sequences, ids)


def _assemble_msa(strings, ids):
    seqs = []
    for (s, i) in zip(strings, ids):
        seq = Seq(s)
        seqs.append(SeqRecord(seq, id=i))
    return MultipleSeqAlignment(seqs)


def concatenate(alignments):
    """
    Concatenates a list of multiple sequence alignment objects.

    The alignments are concatenated based on their label, i.e. the
    sequences from the different alignments which have the same id/labels
    will become a single sequence. The order is preserved.

    If any sequences are missing in one or several alignments, these parts
    are padded with unknown data (:py:class:`Bio.Seq.UnknownSeq`).

    :param alignments: the list of alignments objects, i.e. list(:py:class:`Bio.Align.MultipleSeqAlignment`)
    :returns: a single :py:class:`Bio.Align.MultipleSeqAlignment`

    Example::

        >>> sequences = {'aln1': {'seq1': 'acgtca',
        ...                       'seq2': 'acgtt-',
        ...                       'seq3': 'ac-ta-'},
        ...              'aln2': {'seq2': 'ttg-cta',
        ...                       'seq3': 'tcgacta',
        ...                       'seq4': 'ttgacta'}}
        >>> alignments = [MultipleSeqAlignment([SeqRecord(Seq(sequence), id=key, annotations={"molecule_type": "DNA"})
        ...      for (key, sequence) in sequences[aln].items()])
        ...               for aln in ('aln1', 'aln2')]
        >>> con_alignment = concatenate(alignments)
        >>> con_alignment.sort()
        >>> print(con_alignment)
        ExtendedIUPACDNA() alignment with 4 rows and 13 columns
        acgtcaNNNNNNN seq1
        acgtt-ttg-cta seq2
        ac-ta-tcgacta seq3
        NNNNNNttgacta seq4

    :note:

       Limitations: any annotations in the sub-alignments are lost in
       the concatenated alignment.

    """

    # First check to see whether we're inputting filenames of alignments or the Biopython alignments
    # Assume that it's a biopython alignment if it's not a filename
    tmp_aligns = []
    for filename in alignments:
        if identify_input(filename).name == 'FILENAME':
            tmp_aligns.append(AlignIO.read(filename, "fasta"))
        else:
            tmp_aligns.append(filename)

    # Copy back to alignments
    alignments = tmp_aligns

    # Get the full set of labels (i.e. sequence ids) for all the alignments
    all_labels = set(seq.id for aln in alignments for seq in aln)

    # Make a dictionary to store info as we go along
    # (defaultdict is convenient -- asking for a missing key gives back an empty list)
    tmp = defaultdict(list)

    # try to get molecule_type from sequences
    molecule_type = set(seq.annotations.get('molecule_type') for aln in alignments for seq in aln)
    molecule_type.discard(None)
    if len(molecule_type) == 1:
        molecule_type = molecule_type.pop()
        if molecule_type.upper() in ("DNA", "RNA"):
            unknown_char = "N"
        elif molecule_type.lower() == "protein":
            unknown_char = "X"
        else:
            unknown_char = "?"
    else:
        unknown_char = '?'

    for aln in alignments:
        length = aln.get_alignment_length()

        # check if any labels are missing in the current alignment
        these_labels = set(rec.id for rec in aln)
        missing = all_labels - these_labels

        # if any are missing, create unknown data of the right length,
        # stuff the string representation into the tmp dict
        for label in missing:
            new_seq = unknown_char*length  # UnknownSeq(length, character=unknown_char) # deprecate fro biopython
            tmp[label].append(str(new_seq))

        # else stuff the string representation into the tmp dict
        for rec in aln:
            tmp[rec.id].append(str(rec.seq))

    # Stitch all the substrings together using join (most efficient way),
    # and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment
    return MultipleSeqAlignment(SeqRecord(Seq(''.join(v)), id=k)
                                for (k, v) in tmp.items())


AlignmentInput = Enum('AlignmentInput', 'OBJECT FILENAME')


def identify_input(alignment):
    """
    Work out if we're dealing with a Biopython object (return True), a file
    (return False), or invalid input (raise error)
    """
    try:
        if isinstance(alignment, (MultipleSeqAlignment, types.GeneratorType, list)):
            # `alignment` is a Biopython MultipleSequenceAlignment
            return AlignmentInput.OBJECT
        elif hasattr(alignment, '__next__') and hasattr(alignment, '__iter__'):
            # alignmet is an iterator
            return AlignmentInput.OBJECT
        elif isinstance(alignment, (str, Path)) and os.path.exists(alignment):
            # `alignment` is a filepath
            return AlignmentInput.FILENAME
    except:
        # `alignment` is some other thing we can't handle
        raise ValueError('{} is not an alignment object or a valid filename'.format(alignment))


class BackTranslator(object):
    """
    Reverse translate a protein sequence into a DNA sequence using a
    codon usage table
    """

    def __init__(self, table):
        """
        Initialises the amino acid -> codons, codon weights map from
        the supplied codon usage table -- see http://www.kazusa.or.jp/codon
        """
        self.map = self._codon_table_parser(table)

    def _weighted_choice(self, choices):
        total = sum(w for c, w in choices)
        r = random.uniform(0, total)
        upto = 0
        for c, w in choices:
            if upto + w > r:
                return c
            upto += w
        assert False, "Shouldn't get here"

    def _read_entry(self, gen):
        """
        Reads a single entry from the active codon usage table
        """
        codon = next(gen)
        aa = next(gen)
        next(gen)
        next(gen)
        count = next(gen)
        if count == '(':
            count = next(gen)
        count = int(count.strip('()'))
        return codon.replace('U', 'T'), aa, count

    def _codon_table_parser(self, table):
        """
        Parses the table produced by http://www.kazusa.or.jp/codon
        Returns a dict of amino acid -> codons, codon weights
        """
        table = iter(table.replace('\n', '  ').split())
        codons = defaultdict(list)
        counts = defaultdict(list)
        try:
            while table:
                codon, aa, count = self._read_entry(table)
                codons[aa].append(codon)
                counts[aa].append(count)
        except StopIteration:
            pass

        wmap = {}
        for k in counts:
            s = sum(counts[k])
            scaled = [c / float(s) for c in counts[k]]
            wmap[k] = (codons[k], scaled)
        return wmap

    def back_translate(self, seq):
        """
        Back translate a protein sequence into dna. Codons are chosen randomly according
        to the codon usage frequencies in the table the BackTranslator object was 
        initialised with -- see http://www.kazusa.or.jp/codon
        """
        return ''.join(self._weighted_choice(list(zip(*self.map[char.upper()]))) for char in seq)
