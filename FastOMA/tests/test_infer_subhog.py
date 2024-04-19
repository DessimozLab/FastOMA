from unittest import TestCase
from ete3 import Tree, TreeNode
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import Namespace
from .._hog_class import HOG, Representative
from .._infer_subhog import LevelHOGProcessor


class TestLevelHogProcessor(TestCase):
    def setUp(self):
        genetree = Tree(
            '(((((G00100_SE001||SE001:153.567,G00100_SE008||SE008:153.567)1:39.499[&&NHX:evoltype=S],(G00100_SE006||SE006:173.507,G00100_SE007||SE007:173.507)1:19.5597[&&NHX:evoltype=S])1:14.0196[&&NHX:evoltype=S],(G00100_SE003||SE003:198.481,((((G00100_SE011||SE011:136.533,G00100_SE012||SE012:136.533)1:7.60673[&&NHX:evoltype=S],(G00100_SE010||SE010:36.1782,G00342_SE010||SE010:36.1782)1:107.961[&&NHX:evoltype=D])1:8.49419[&&NHX:evoltype=S],G00100_SE009||SE009:152.634)1:13.723[&&NHX:evoltype=S],(((G00186_SE004||SE004:143.819,(G00186_SE011||SE011:136.533,(G00186_SE012||SE012:116.411,G00242_SE012||SE012:116.411)1:20.1214[&&NHX:evoltype=D])1:7.28662[&&NHX:evoltype=S])1:0.32011[&&NHX:evoltype=S],(G00186_SE010||SE010:31.4887,G00350_SE010||SE010:31.4887)1:112.651[&&NHX:evoltype=D])1:8.49419[&&NHX:evoltype=S],G00186_SE009||SE009:152.634)1:13.723[&&NHX:evoltype=S])1:32.1245[&&NHX:evoltype=D])1:8.60492[&&NHX:evoltype=S])1:36.2336[&&NHX:evoltype=S],(((G00110_SE001||SE001:153.567,G00110_SE008||SE008:153.567)1:39.499[&&NHX:evoltype=S],(G00110_SE006||SE006:173.507,G00110_SE007||SE007:173.507)1:19.5597[&&NHX:evoltype=S])1:14.0196[&&NHX:evoltype=S],(G00110_SE003||SE003:198.481,(((G00110_SE004||SE004:143.819,(G00110_SE011||SE011:136.533,G00110_SE012||SE012:136.533)1:7.28662[&&NHX:evoltype=S])1:0.32011[&&NHX:evoltype=S],G00110_SE010||SE010:144.139)1:8.49419[&&NHX:evoltype=S],G00110_SE009||SE009:152.634)1:45.8474[&&NHX:evoltype=S])1:8.60492[&&NHX:evoltype=S])1:36.2336[&&NHX:evoltype=S])1:6.68041[&&NHX:evoltype=D],(G00100_SE002||SE002:119.545,(G00100_SE013||SE013:97.4899,(G00100_SE014||SE014:87.2367,G00100_SE015||SE015:87.2367)1:10.2532[&&NHX:evoltype=S])1:22.055[&&NHX:evoltype=S])1:130.455[&&NHX:evoltype=S]);')
        sptree = Tree("dummy;")
        hogs = [HOG(SeqRecord(Seq("AAAAAA"), id=n.name), sptree, "test1") for n in genetree.iter_leaves()]
        conf = Namespace(msa_write=False, gene_tree_write=False, number_of_samples_per_hog=5)
        self.genetree = genetree
        self.lp = LevelHOGProcessor(sptree, hogs, "test1", conf)

    def test_propose_representatives(self):
        rep = self.lp.find_most_divergent_representatives_from_genetree(self.genetree)
        self.assertEqual(len(rep), self.lp.conf.number_of_samples_per_hog)
        self.assertIn(rep[self.lp.conf.number_of_samples_per_hog-1].get_id(), ("G00100_SE013||SE013","G00100_SE015||SE015","G00100_SE014||SE014","G00100_SE002||SE002"))

    def test_reconcilation(self):
        exp = self.genetree.write(features=['evoltype'])
        self.lp.infer_reconciliation(genetree=self.genetree)
        self.assertEqual(exp, self.genetree.write(features=['evoltype']))
        self.assertEqual(self.genetree.sos, 0)
        self.assertEqual(self.genetree.children[0].sos, 1)


