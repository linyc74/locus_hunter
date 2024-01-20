from typing import List
from ngslite import write_genbank, Chromosome
from .sort_loci import SortLoci
from .add_color import AddColor
from .view_loci import ViewLoci
from .extract_loci import ExtractLoci
from .template import Settings, Processor
from .constant import CDS_ID_KEY, ORTHOLOG_ID_KEY


class LocusHunter(Processor):

    query_faa: str
    gbk_dir: str
    evalue: float
    extension: int
    ortholog_identity: float
    dereplicate_loci: bool
    label_attributes: List[str]
    output: str

    loci: List[Chromosome]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.loci = []

    def main(
            self,
            query_faa: str,
            gbk_dir: str,
            evalue: float,
            extension: int,
            ortholog_identity: float,
            dereplicate_loci: bool,
            label_attributes: List[str],
            output: str):

        self.query_faa = query_faa
        self.gbk_dir = gbk_dir
        self.evalue = evalue
        self.extension = extension
        self.ortholog_identity = ortholog_identity
        self.dereplicate_loci = dereplicate_loci
        self.label_attributes = label_attributes
        self.output = output

        self.extract_loci()

        if len(self.loci) == 0:
            self.logger.info('No locus found. Abort')
            return

        self.sort_loci()
        self.add_color()
        self.view_loci()
        self.save_genbank()

    def extract_loci(self):
        self.loci = ExtractLoci(self.settings).main(
            query_faa=self.query_faa,
            gbk_dir=self.gbk_dir,
            evalue=self.evalue,
            extension=self.extension)

    def sort_loci(self):
        self.loci = SortLoci(self.settings).main(
            loci=self.loci,
            ortholog_identity=self.ortholog_identity,
            dereplicate_loci=self.dereplicate_loci)

    def add_color(self):
        self.loci = AddColor(self.settings).main(
            loci=self.loci)

    def view_loci(self):
        ViewLoci(self.settings).main(
            loci=self.loci,
            output=self.output,
            label_attributes=self.label_attributes)

    def save_genbank(self):
        self.remove_intermediate_labels()
        write_genbank(
            data=self.loci,
            file=f'{self.output}.gbk',
            use_locus_text=False)

    def remove_intermediate_labels(self):
        intermediate_keys = [CDS_ID_KEY, ORTHOLOG_ID_KEY]
        for locus in self.loci:
            for feature in locus.features:
                for key in intermediate_keys:
                    feature.remove_attribute(key)
