import random
from locus_hunter.locus_hunter import LocusHunter
from .setup import TestCase, remove_genbank_date_str


class TestLocusHunter(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        random.seed(1)

        LocusHunter(settings=self.settings).main(
            query_faa=f'{self.indir}/query.faa',
            gbk_dir=f'{self.indir}/gbk_dir',
            evalue=1e-20,
            extension=5000,
            min_hits_per_locus=1,
            ortholog_identity=0.9,
            dereplicate_loci=True,
            include_locus_names=[],
            label_attributes=['gene', 'locus_tag'],
            loci_per_plot=3,
            dpi=600,
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
            evalue=1e-20,
            extension=5000,
            min_hits_per_locus=1,
            ortholog_identity=0.9,
            dereplicate_loci=True,
            include_locus_names=[],
            label_attributes=['gene', 'locus_tag'],
            loci_per_plot=3,
            dpi=600,
            output=f'{self.outdir}/output'
        )
