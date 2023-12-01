import os
from copy import deepcopy
from typing import List, Tuple
from ngslite import Chromosome, read_genbank, get_files, GenericFeature
from .blast import Blastp
from .template import Processor
from .constant import CDS_ID_KEY


JOINER = '___'
INTERVAL_TYPE = List[Tuple[int, int]]


class ExtractLoci(Processor):

    query_faa: str
    gbk_dir: str
    evalue: float
    extension: int

    loci = List[Chromosome]

    def main(
            self,
            query_faa: str,
            gbk_dir: str,
            evalue: float,
            extension: int) -> List[Chromosome]:

        self.query_faa = query_faa
        self.gbk_dir = gbk_dir
        self.evalue = evalue
        self.extension = extension

        self.loci = []

        gbks = self.get_gbks()
        for gbk in gbks:
            chromosomes = self.read_genbank(gbk)
            for chromosome in chromosomes:
                self.extrac_loci_from(chromosome)

        return self.loci

    def get_gbks(self) -> List[str]:
        files = get_files(source=self.gbk_dir, isfullpath=True)
        gbks = [f for f in files
                if not os.path.basename(f).startswith('.')]
        return gbks

    def read_genbank(self, gbk: str) -> List[Chromosome]:
        return ReadGenbank(self.settings).main(gbk=gbk)

    def extrac_loci_from(self, chromosome: Chromosome):
        loci = self.get_loci_from(chromosome)
        self.log(chromosome=chromosome, loci=loci)
        self.loci += loci

    def get_loci_from(self, chromosome: Chromosome) -> List[Chromosome]:
        return GetLociFromChromosome(self.settings).main(
            query_faa=self.query_faa,
            chromosome=chromosome,
            evalue=self.evalue,
            extension=self.extension)

    def log(self,
            chromosome: Chromosome,
            loci: List[Chromosome]):
        msg = f'{chromosome.seqname} -> {len(loci)} loci'
        self.logger.info(msg)


class ReadGenbank(Processor):

    gbk: str

    def main(self, gbk: str) -> List[Chromosome]:
        self.gbk = gbk

        chromosomes = read_genbank(file=self.gbk)

        for chromosome in chromosomes:
            self.modify_one(chromosome)

        return chromosomes

    def modify_one(self, chromosome: Chromosome):
        fname = os.path.basename(self.gbk)
        chromosome.seqname = f'{fname}{JOINER}{chromosome.seqname}'
        for i, feature in enumerate(chromosome.features):
            if feature.type == 'CDS':
                feature.add_attribute(
                    key=CDS_ID_KEY,
                    val=f'{chromosome.seqname}{JOINER}{i+1}')


class GetLociFromChromosome(Processor):

    query_faa: str
    chromosome: Chromosome
    evalue: float
    extension: int

    cds_features: List[GenericFeature]
    cds_ids: List[str]
    intervals: INTERVAL_TYPE
    merged_intervals: INTERVAL_TYPE
    loci: List[Chromosome]

    def main(
            self,
            query_faa: str,
            chromosome: Chromosome,
            evalue: float,
            extension: int) -> List[Chromosome]:

        self.query_faa = query_faa
        self.chromosome = chromosome
        self.evalue = evalue
        self.extension = extension

        self.set_cds_features()

        if len(self.cds_features) == 0:
            return []

        self.set_cds_ids()
        self.set_intervals()
        self.set_merged_intervals()
        self.set_loci()

        return self.loci

    def set_cds_features(self):
        self.cds_features = [
            f for f in self.chromosome.features
            if f.type == 'CDS'
        ]

    def set_cds_ids(self):

        library = self.__get_library_faa_data()

        df = Blastp(self.settings).main(
            query=self.query_faa,
            library=library,
            evalue=self.evalue)

        self.cds_ids = sorted(df['subject'].unique())

    def __get_library_faa_data(self) -> List[Tuple[str, str]]:
        ret = []
        for feature in self.cds_features:
            cds_id = feature.get_attribute(key=CDS_ID_KEY)
            seq = feature.get_attribute(key='translation')
            if seq is not None:
                ret.append((cds_id, seq))
        return ret

    def set_intervals(self):

        cds_dict = {c: True for c in self.cds_ids}

        self.intervals = []
        for feature in self.cds_features:
            cds_id = feature.get_attribute(key=CDS_ID_KEY)
            in_cds_dict = cds_dict.get(cds_id, False)
            if in_cds_dict:
                start = feature.start - self.extension
                end = feature.end + self.extension
                self.intervals.append((start, end))

    def set_merged_intervals(self):

        m = self.intervals[0:1]

        for this in self.intervals[1:]:
            last = m[-1]
            if overlap(this, last):
                m[-1] = (last[0], this[1])
            else:
                m.append(this)

        self.merged_intervals = m

    def set_loci(self):
        self.loci = []
        for start, end in self.merged_intervals:
            locus = self.__get_locus(start=start, end=end)
            self.loci.append(locus)

    def __get_locus(self, start: int, end: int) -> Chromosome:
        locus = deepcopy(self.chromosome)
        locus.genbank_locus_text = ''
        locus.seqname = f'{locus.seqname}{JOINER}{start:,d}-{end:,d}'
        locus.crop(start, end)
        locus.circular = False
        return locus


def overlap(
        a: Tuple[int, int],
        b: Tuple[int, int]) -> bool:

    start_a, end_a = sorted(a)
    start_b, end_b = sorted(b)

    return (end_a >= start_b) and (end_b >= start_a)
