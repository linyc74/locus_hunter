import shutil
from ngslite import write_genbank
from locus_hunter.template import Settings
from locus_hunter.extract_loci import ExtractLoci
from .tools import setup_dir, TestCase, remove_genbank_date_str


class TestExtractLoci(TestCase):

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
            file1=f'{self.indir}/loci.gbk',
            file2=f'{self.outdir}/loci.gbk')

    def test_no_hit(self):
        loci = ExtractLoci(settings=self.settings).main(
            query_faa=f'{self.indir}/wrong_query.faa',
            gbk_dir=f'{self.indir}/gbk_dir',
            extension=5000,
            evalue=1e-20)

        self.assertEqual([], loci)
