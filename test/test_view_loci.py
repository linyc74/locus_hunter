from ngslite import read_genbank, GenericFeature
from locus_hunter.view_loci import GenericToGraphicFeature, ChromosomeToGraphicRecord, ViewLoci, split_list
from .setup import TestCase


class TestViewLoci(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        loci = read_genbank(file=f'{self.indir}/colored_loci.gbk')
        output = f'{self.outdir}/output'

        ViewLoci(self.settings).main(
            loci=loci,
            output=output,
            label_attributes=['gene', 'locus_tag'],
            loci_per_plot=19,
            dpi=300
        )

    def test_very_long_name(self):
        loci = read_genbank(file=f'{self.indir}/colored_loci.gbk')
        for locus in loci:
            locus.seqname = 'x' * 100
        output = f'{self.outdir}/output'

        ViewLoci(self.settings).main(
            loci=loci,
            output=output,
            label_attributes=['gene', 'locus_tag'],
            loci_per_plot=10,
            dpi=300
        )


class TestChromosomeToGraphicRecord(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        chromosome = read_genbank(file=f'{self.indir}/colored_loci.gbk')[0]
        record = ChromosomeToGraphicRecord(self.settings).main(
            chromosome=chromosome,
            label_attributes=['gene', 'locus_tag']
        )

        self.assertEqual(11143, record.sequence_length)

        first_three = str(record.features[0:3])
        expected = '[GF(AK972_RS05835, 2936-4651 (-1)), GF(AK972_RS05840, 5001-6143 (-1)), GF(AK972_RS05845, 6252-7160 (1))]'
        self.assertEqual(expected, first_three)


class TestGenericToGraphicFeature(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):

        generic_feature = GenericFeature(
            seqname='.',
            type_='CDS',
            start=1,
            end=99,
            strand='+')

        graphic_feature = GenericToGraphicFeature(self.settings).main(
            generic_feature=generic_feature,
            label_attributes=['gene', 'locus_tag']
        )

        expected = 'GF(None, 1-99 (1))'

        self.assertTrue(expected, str(graphic_feature))


class TestFunctions(TestCase):

    def test_split_list_dividable(self):
        actual = split_list(ls=[1], size=1)
        self.assertListEqual([[1]], actual)

        actual = split_list(ls=[1, 2], size=1)
        self.assertListEqual([[1], [2]], actual)

        actual = split_list(ls=[1, 2, 3, 4], size=2)
        self.assertListEqual([[1, 2], [3, 4]], actual)

    def test_split_list_non_dividable(self):
        actual = split_list(ls=[1, 2, 3], size=2)
        self.assertListEqual([[1, 2], [3]], actual)

        actual = split_list(ls=[1, 2, 3, 4, 5], size=2)
        self.assertListEqual([[1, 2], [3, 4], [5]], actual)
