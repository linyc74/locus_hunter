from ngslite import write_genbank
from locus_hunter.extract_loci import ExtractLoci
from .setup import TestCase, remove_genbank_date_str


class TestExtractLoci(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        loci = ExtractLoci(settings=self.settings).main(
            query_faa=f'{self.indir}/query.faa',
            gbk_dir=f'{self.indir}/gbk_dir',
            extension=5000,
            evalue=1e-3)

        write_genbank(
            data=loci,
            file=f'{self.outdir}/loci.gbk',
            use_locus_text=False)

        remove_genbank_date_str(f'{self.outdir}/loci.gbk')

        self.assertFileEqual(
            first=f'{self.indir}/loci.gbk',
            second=f'{self.outdir}/loci.gbk')

    def test_no_hit(self):
        loci = ExtractLoci(settings=self.settings).main(
            query_faa=f'{self.indir}/wrong_query.faa',
            gbk_dir=f'{self.indir}/gbk_dir',
            extension=5000,
            evalue=1e-20)

        self.assertEqual([], loci)
