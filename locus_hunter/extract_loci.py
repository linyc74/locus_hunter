import os
from copy import deepcopy
from typing import List, Tuple
from ngslite import Chromosome, read_genbank, get_files, GenericFeature
from locus_hunter.template import Processor, Settings
from .blast import Blastp
from .constant import CDS_ID_KEY


JOINER = '___'


class ReadGenbank(Processor):

    gbk: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

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


INTERVAL_TYPE = List[Tuple[int, int]]


class GetLociFromChromosome(Processor):

    query_faa: str
    chromosome: Chromosome
    evalue: float
    extension: int

    cds_features: List[GenericFeature]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

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

        cds_ids = self.get_cds_ids()
        intervals = self.get_intervals(cds_ids)
        merged_intervals = self.merge(intervals)
        loci = self.get_loci(merged_intervals)

        return loci

    def set_cds_features(self):
        self.cds_features = [
            f for f in self.chromosome.features
            if f.type == 'CDS']

    def get_cds_ids(self) -> List[str]:

        library = self.get_library_faa_data()

        df = Blastp(self.settings).main(
            query=self.query_faa,
            library=library,
            evalue=self.evalue)

        cds_ids = sorted(df['subject'].unique())

        return cds_ids

    def get_library_faa_data(self) -> List[Tuple[str, str]]:
        ret = []
        for feature in self.cds_features:
            cds_id = feature.get_attribute(key=CDS_ID_KEY)
            seq = feature.get_attribute(key='translation')
            if seq is not None:
                ret.append((cds_id, seq))
        return ret

    def get_intervals(self, cds_ids: List[str]) -> INTERVAL_TYPE:

        cds_dict = {c: True for c in cds_ids}

        ret = []
        for feature in self.cds_features:
            cds_id = feature.get_attribute(key=CDS_ID_KEY)
            in_cds_dict = cds_dict.get(cds_id, False)
            if in_cds_dict:
                start = feature.start - self.extension
                end = feature.end + self.extension
                ret.append((start, end))

        return ret

    def merge(self, intervals: INTERVAL_TYPE) -> INTERVAL_TYPE:

        ret = intervals[0:1]

        for this in intervals[1:]:
            last = ret[-1]
            if overlap(this, last):
                ret[-1] = (last[0], this[1])
            else:
                ret.append(this)

        return ret

    def get_loci(self, merged_intervals: INTERVAL_TYPE) -> List[Chromosome]:

        ret = []
        for start, end in merged_intervals:
            locus = self.__get_locus(start=start, end=end)
            ret.append(locus)

        return ret

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


class ExtractLoci(Processor):

    query_faa: str
    gbk_dir: str
    evalue: float
    extension: int

    loci = List[Chromosome]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

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
