import unittest


from Bio import AlignIO
from FastOMA._wrappers import infer_gene_tree
from FastOMA.zoo.wrappers import WrapperError
import pathlib
this_dir = pathlib.Path(__file__).parent


class FastTreeTester(unittest.TestCase):
    def test_failing_tree_building_reports_error_from_fasttree(self):
        msa = AlignIO.read(this_dir / "data" / "failing-msa.fa", "fasta")
        with self.assertLogs("FastOMA", level="ERROR") as cm:
            with self.assertRaises(WrapperError):
                infer_gene_tree(msa, "/tmp")
            self.assertIn("Non-unique name", "\n".join(cm.output))

    def test_treebuilding_with_correct_msa(self):
        msa = AlignIO.read(this_dir / "data" / "correct-msa.fa", "fasta")
        gene_rooting_method= ""
        tree = infer_gene_tree(msa, "/tmp", gene_rooting_method)
        self.assertIn("HUMAN01350", tree)