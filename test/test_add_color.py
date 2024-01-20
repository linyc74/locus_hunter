import random
from ngslite import read_genbank, write_genbank
from locus_hunter.add_color import ColorDictGenerator, AddColor
from .setup import TestCase, remove_genbank_date_str


class TestColorDictGenerator(TestCase):

    def test_set_keys(self):
        generator = ColorDictGenerator()
        generator.set_keys(non_unique_keys=[1, 2, 2, 2, 3, 5, 5, 5, 5, 1])
        self.assertListEqual([5, 2, 1, 3], generator.keys)

    def test_main(self):
        random.seed(1)
        color_dict = ColorDictGenerator().main(
            keys=[1, 2, 2, 2, 3, 5, 5, 5, 5, 1])
        expected = {5: '#1F77B4', 2: '#FF7F0E', 1: '#2CA02C', 3: '#D62728'}
        self.assertDictEqual(expected, color_dict)


class TestAddColor(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):

        random.seed(1)

        sorted_loci = read_genbank(file=f'{self.indir}/sorted_loci.gbk')

        colored_loci = AddColor(settings=self.settings).main(
            loci=sorted_loci)

        write_genbank(
            data=colored_loci,
            file=f'{self.outdir}/colored_loci.gbk',
            use_locus_text=False)

        remove_genbank_date_str(f'{self.outdir}/colored_loci.gbk')

        self.assertFileEqual(
            first=f'{self.indir}/colored_loci.gbk',
            second=f'{self.outdir}/colored_loci.gbk'
        )
