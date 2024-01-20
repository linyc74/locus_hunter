import shutil
from locus_hunter.cd_hit import CdHit
from .setup import TestCase


class TestCdHit(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

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
