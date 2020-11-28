import shutil
from locus_hunter.template import Settings
from locus_hunter.cd_hit import CdHit
from .tools import setup_dir, TestCase


class TestCdHit(TestCase):

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
        actual = CdHit(self.settings).main(
            faa=f'{self.indir}/in.faa',
            sequence_identity=0.9)

        expected = {
            'WP_001262893.1': 1,
            'WP_001175410.1': 2,
            'WP_001175410.2': 2,
            'WP_001175410.3': 2,
            'WP_001862796.1': 3
        }

        self.assertDictEqual(expected, actual)
