from ngslite import read_genbank, GenericFeature
from locus_hunter.view_loci import GenericToGraphicFeature, ChromosomeToGraphicRecord, ViewLoci
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
            label_attributes=['gene', 'locus_tag']
        )

    def test_very_long_name(self):
        loci = read_genbank(file=f'{self.indir}/colored_loci.gbk')
        for locus in loci:
            locus.seqname = 'x' * 100
        output = f'{self.outdir}/output'

        ViewLoci(self.settings).main(
            loci=loci,
            output=output,
            label_attributes=['gene', 'locus_tag']
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

        self.assertEqual(20536, record.sequence_length)

        first_three = str(record.features[0:3])
        expected = '[GF(HP_RS05040, 211-1386 (1)), GF(fur, 1546-1998 (1)), GF(HP_RS05050, 2018-2515 (-1))]'
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
