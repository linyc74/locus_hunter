import shutil
from ngslite import read_genbank, GenericFeature
from locus_hunter.template import Settings
from locus_hunter.view_loci import GenericToGraphicFeature, ChromosomeToGraphicRecord, ViewLoci
from .tools import setup_dir, TestCase


class TestGenericToGraphicFeature(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dir(__file__)
        self.settings = Settings(
            workdir=self.workdir,
            outdir=self.outdir,
            threads=4,
            debug=True)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_main(self):

        generic_feature = GenericFeature(
            seqname='.',
            type_='CDS',
            start=1,
            end=99,
            strand='+')

        graphic_feature = GenericToGraphicFeature(self.settings).main(
            generic_feature=generic_feature)

        expected = 'GF(None, 1-99 (1))'

        self.assertTrue(expected, str(graphic_feature))


class TestChromosomeToGraphicRecord(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dir(__file__)
        self.settings = Settings(
            workdir=self.workdir,
            outdir=self.outdir,
            threads=4,
            debug=True)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_main(self):
        chromosome = read_genbank(file=f'{self.indir}/colored_loci.gbk')[0]
        record = ChromosomeToGraphicRecord(self.settings).main(chromosome)

        self.assertEqual(20536, record.sequence_length)

        first_three = str(record.features[0:3])
        expected = '[GF(None, 211-1386 (1)), GF(None, 1546-1998 (1)), GF(None, 2018-2515 (-1))]'
        self.assertEqual(expected, first_three)


class TestViewLoci(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dir(__file__)
        self.settings = Settings(
            workdir=self.workdir,
            outdir=self.outdir,
            threads=4,
            debug=True)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_main(self):
        loci = read_genbank(file=f'{self.indir}/colored_loci.gbk')
        output = f'{self.outdir}/output'

        ViewLoci(self.settings).main(
            loci=loci,
            output=output)

    def test_very_long_name(self):
        loci = read_genbank(file=f'{self.indir}/colored_loci.gbk')
        for locus in loci:
            locus.seqname = 'x' * 100
        output = f'{self.outdir}/output'

        ViewLoci(self.settings).main(
            loci=loci,
            output=output)
