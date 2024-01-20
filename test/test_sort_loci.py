from ngslite import read_genbank, write_genbank
from locus_hunter.sort_loci import SortLoci, SmithWatermanAligner
from .setup import TestCase, remove_genbank_date_str


class TestSortLoci(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):

        loci = read_genbank(file=f'{self.indir}/loci.gbk')

        sorted_loci = SortLoci(settings=self.settings).main(
            loci=loci,
            ortholog_identity=0.9,
            dereplicate_loci=True,
            include_locus_names=['Pseudomonas_aeruginosa_UCBPP-PA14_109']
        )

        write_genbank(
            data=sorted_loci,
            file=f'{self.outdir}/sorted_loci.gbk',
            use_locus_text=False)

        remove_genbank_date_str(f'{self.outdir}/sorted_loci.gbk')

        self.assertFileEqual(
            first=f'{self.indir}/sorted_loci.gbk',
            second=f'{self.outdir}/sorted_loci.gbk')


class TestSmithWatermanAligner(TestCase):

    def test_main(self):
        list1 = [1, 2, 3, 4]
        list2 =    [2, 0, 4, 5]

        score = SmithWatermanAligner().run(list1=list1, list2=list1)
        self.assertEqual(400, score)

        score = SmithWatermanAligner().run(list1=list1, list2=list2)
        self.assertEqual(198, score)
