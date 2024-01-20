import shutil
import random
from locus_hunter.template import Settings
from locus_hunter.locus_hunter import LocusHunter
from .tools import setup_dir, TestCase, remove_genbank_date_str


class TestLocusHunter(TestCase):

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
        random.seed(1)

        LocusHunter(settings=self.settings).main(
            query_faa=f'{self.indir}/query.faa',
            gbk_dir=f'{self.indir}/gbk_dir',
            extension=5000,
            evalue=1e-20,
            ortholog_identity=0.9,
            label_attributes=['gene', 'locus_tag'],
            output=f'{self.outdir}/output'
        )

        remove_genbank_date_str(f'{self.outdir}/output.gbk')

        actual = f'{self.outdir}/output.gbk'
        expected = f'{self.indir}/output.gbk'

        self.assertFileEqual(expected, actual)

    def test_no_hit(self):
        LocusHunter(settings=self.settings).main(
            query_faa=f'{self.indir}/wrong_query.faa',
            gbk_dir=f'{self.indir}/gbk_dir',
            extension=5000,
            evalue=1e-20,
            ortholog_identity=0.9,
            label_attributes=['gene', 'locus_tag'],
            output=f'{self.outdir}/output'
        )
