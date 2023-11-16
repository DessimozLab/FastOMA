import unittest

import FastOMA._config
from FastOMA._wrappers import infer_gene_tree
from FastOMA.zoo.wrappers import WrapperError
import pathlib
this_dir = pathlib.Path(__file__).parent


class FastTreeTester(unittest.TestCase):
    def test_failing_tree_building_reports_error_from_fasttree(self):
        with self.assertLogs("FastOMA", level="ERROR") as cm:
            with self.assertRaises(WrapperError):
                infer_gene_tree(this_dir / "data" / "failing-msa.fa", "/tmp")
            self.assertIn("Non-unique name", "\n".join(cm.output))

    def test_treebuilding_with_correct_msa(self):
        tree = infer_gene_tree(this_dir / "data" / "correct-msa.fa", "/tmp")
        self.assertIn("HUMAN01350", tree)