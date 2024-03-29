from typing import List
from ngslite import write_genbank, Chromosome
from .sort_loci import SortLoci
from .add_color import AddColor
from .view_loci import ViewLoci
from .extract_loci import ExtractLoci
from .template import Processor
from .constant import CDS_ID_KEY, ORTHOLOG_ID_KEY


class LocusHunter(Processor):

    query_faa: str
    gbk_dir: str
    evalue: float
    extension: int
    min_hits_per_locus: int
    ortholog_identity: float
    dereplicate_loci: bool
    include_locus_names: List[str]
    label_attributes: List[str]
    loci_per_plot: int
    dpi: int
    output: str

    loci: List[Chromosome]

    def main(
            self,
            query_faa: str,
            gbk_dir: str,
            evalue: float,
            extension: int,
            min_hits_per_locus: int,
            ortholog_identity: float,
            dereplicate_loci: bool,
            include_locus_names: List[str],
            label_attributes: List[str],
            loci_per_plot: int,
            dpi: int,
            output: str):

        self.query_faa = query_faa
        self.gbk_dir = gbk_dir
        self.evalue = evalue
        self.extension = extension
        self.min_hits_per_locus = min_hits_per_locus
        self.ortholog_identity = ortholog_identity
        self.dereplicate_loci = dereplicate_loci
        self.include_locus_names = include_locus_names
        self.label_attributes = label_attributes
        self.loci_per_plot = loci_per_plot
        self.dpi = dpi
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
            extension=self.extension,
            min_hits_per_locus=self.min_hits_per_locus)

    def sort_loci(self):
        self.loci = SortLoci(self.settings).main(
            loci=self.loci,
            ortholog_identity=self.ortholog_identity,
            dereplicate_loci=self.dereplicate_loci,
            include_locus_names=self.include_locus_names)

    def add_color(self):
        self.loci = AddColor(self.settings).main(
            loci=self.loci)

    def view_loci(self):
        ViewLoci(self.settings).main(
            loci=self.loci,
            output=self.output,
            label_attributes=self.label_attributes,
            loci_per_plot=self.loci_per_plot,
            dpi=self.dpi)

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
