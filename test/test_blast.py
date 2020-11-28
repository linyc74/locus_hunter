import shutil
import pandas as pd
from ngslite import read_fasta
from locus_hunter.template import Settings
from locus_hunter.blast import Blastp
from .tools import setup_dir, TestCase


class TestBlastp(TestCase):

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

    def test_no_hit(self):

        df = Blastp(settings=self.settings).main(
            query=f'{self.indir}/query.faa',
            library=f'{self.indir}/library.faa',
            evalue=1e-50)

        self.assertTrue(len(df) == 0)

    def test_library_as_faa_data(self):

        library = read_fasta(f'{self.indir}/library.faa')

        df = Blastp(settings=self.settings).main(
            query=f'{self.indir}/query.faa',
            library=library,
            evalue=1e-20)

        expected = pd.read_csv(f'{self.indir}/output.csv')

        self.assertTrue(df.equals(expected))
