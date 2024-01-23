import matplotlib.pyplot as plt
from typing import List, Optional, Dict
from ngslite import GenericFeature, Chromosome
from dna_features_viewer import GraphicFeature, GraphicRecord
from .constant import COLOR_KEY
from .template import Processor


class ViewLoci(Processor):

    loci: List[Chromosome]
    output: str
    label_attributes: List[str]
    loci_per_plot: int
    dpi: int

    graphic_records: List[GraphicRecord]

    def main(
            self,
            loci: List[Chromosome],
            output: str,
            label_attributes: List[str],
            loci_per_plot: int,
            dpi: int):

        self.loci = loci
        self.output = output
        self.label_attributes = label_attributes
        self.loci_per_plot = loci_per_plot
        self.dpi = dpi

        self.set_graphic_records()
        self.plot_graphic_records()

    def set_graphic_records(self):
        self.graphic_records = []
        for locus in self.loci:
            record = ChromosomeToGraphicRecord(self.settings).main(
                chromosome=locus,
                label_attributes=self.label_attributes)
            self.graphic_records.append(record)

    def plot_graphic_records(self):
        PlotGraphicRecords(self.settings).main(
            graphic_records=self.graphic_records,
            loci_per_plot=self.loci_per_plot,
            dpi=self.dpi,
            output=self.output)


#


class ChromosomeToGraphicRecord(Processor):

    FEATURE_TYPES = ['CDS']

    chromosome: Chromosome
    label_attributes: List[str]

    graphic_features: List[GraphicFeature]

    def main(
            self,
            chromosome: Chromosome,
            label_attributes: List[str]) -> GraphicRecord:

        self.chromosome = chromosome
        self.label_attributes = label_attributes

        self.set_graphic_features()

        r = GraphicRecord(
            sequence_length=len(self.chromosome.sequence),
            sequence=None,
            features=self.graphic_features,
            feature_level_height=1,
            first_index=0,
            plots_indexing='biopython',
            labels_spacing=8,
            ticks_resolution='auto')
        r.seqname = self.chromosome.seqname  # attach seqname to the graphic record
        return r

    def set_graphic_features(self):
        self.graphic_features = []
        for feature in self.chromosome.features:
            if feature.type not in self.FEATURE_TYPES:
                continue
            graphic_feature = GenericToGraphicFeature(self.settings).main(
                generic_feature=feature,
                label_attributes=self.label_attributes)
            self.graphic_features.append(graphic_feature)


class GenericToGraphicFeature(Processor):

    STRAND_DICT = {
        '+': 1,
        '-': -1
    }
    DEFAULT_COLOR = '#FFFFFF'
    THICKNESS = 8.0
    LINEWIDTH = 0.5
    LINECOLOR = '#000000'
    OPEN_LEFT = False
    OPEN_RIGHT = False
    FONTDICT = {'fontsize': 6}

    generic_feature: GenericFeature
    label_attributes: List[str]

    strand: int
    label: Optional[str]
    color: str

    graphic_feature: GraphicFeature

    def main(
            self,
            generic_feature: GenericFeature,
            label_attributes: List[str]) -> GraphicFeature:

        self.generic_feature = generic_feature
        self.label_attributes = label_attributes

        self.set_strand()
        self.set_label()
        self.set_color()
        self.set_graphic_feature()

        return self.graphic_feature

    def set_strand(self):
        f = self.generic_feature
        self.strand = self.STRAND_DICT.get(f.strand, 0)

    def set_label(self):
        self.label = None
        for attr in self.label_attributes:
            value = self.generic_feature.get_attribute(key=attr)
            if value is not None:
                self.label = value
                break

    def set_color(self):
        f = self.generic_feature
        self.color = f.get_attribute(key=COLOR_KEY)
        if self.color is None:
            self.color = self.DEFAULT_COLOR

    def set_graphic_feature(self):
        self.graphic_feature = GraphicFeature(
            start=self.generic_feature.start,
            end=self.generic_feature.end,
            strand=self.strand,
            label=self.label,
            color=self.color,
            thickness=self.THICKNESS,
            linewidth=self.LINEWIDTH,
            linecolor=self.LINECOLOR,
            fontdict=self.FONTDICT,
            html=None,
            open_left=self.OPEN_LEFT,
            open_right=self.OPEN_RIGHT,
            box_linewidth=1,
            box_color='auto',
            legend_text=None,
            label_link_color='black')


#


class PlotGraphicRecords(Processor):

    graphic_records: List[GraphicRecord]
    loci_per_plot: int
    dpi: int
    output: str

    output_to_graphic_records: Dict[str, List[GraphicRecord]]

    def main(
            self,
            graphic_records: List[GraphicRecord],
            loci_per_plot: int,
            dpi: int,
            output: str):

        self.graphic_records = graphic_records
        self.loci_per_plot = loci_per_plot
        self.dpi = dpi
        self.output = output

        self.split_graphic_records()
        self.make_plots()

    def split_graphic_records(self):
        self.output_to_graphic_records = SplitGraphicRecords(self.settings).main(
            graphic_records=self.graphic_records,
            loci_per_plot=self.loci_per_plot,
            output_prefix=self.output)

    def make_plots(self):
        for output, graphic_records in self.output_to_graphic_records.items():
            MakeOnePlot(self.settings).main(
                graphic_records=graphic_records,
                dpi=self.dpi,
                output=output)


class SplitGraphicRecords(Processor):

    graphic_records: List[GraphicRecord]
    loci_per_plot: int
    output_prefix: str

    chunk_of_records: List[List[GraphicRecord]]
    output: List[str]

    output_to_graphic_records: Dict[str, List[GraphicRecord]]

    def main(
            self,
            graphic_records: List[GraphicRecord],
            loci_per_plot: int,
            output_prefix: str) -> Dict[str, List[GraphicRecord]]:

        self.graphic_records = graphic_records
        self.loci_per_plot = loci_per_plot
        self.output_prefix = output_prefix

        self.split_graphic_records_into_chunks()
        self.set_output_to_graphic_records()

        return self.output_to_graphic_records

    def split_graphic_records_into_chunks(self):
        self.chunk_of_records = split_list(
            ls=self.graphic_records,
            size=self.loci_per_plot)

    def set_output_to_graphic_records(self):
        if len(self.chunk_of_records) == 1:
            self.output_to_graphic_records = {
                self.output_prefix: self.chunk_of_records[0]
            }
            return

        self.output_to_graphic_records = {}
        for i, records in enumerate(self.chunk_of_records):
            self.output_to_graphic_records[f'{self.output_prefix}_{i + 1}'] = records


class MakeOnePlot(Processor):

    HEIGHT_CM_PER_LOCUS = 2 / 2.54
    WIDTH_CM_PER_KB = 2 / 2.54
    WIDTH_CM_PER_CHAR = 0.2 / 2.54
    SCALE_BAR_KB = 5000

    graphic_records: List[GraphicRecord]
    dpi: int
    output: str

    seqnames: List[str]
    figure: plt.Figure
    axs = List[plt.Axes]

    def main(
            self,
            graphic_records: List[GraphicRecord],
            dpi: int,
            output: str):

        self.graphic_records = graphic_records
        self.dpi = dpi
        self.output = output

        self.set_seqnames()
        self.init_figure()
        self.plot_graphic_records()
        self.config_figure()
        self.save_output()

    def set_seqnames(self):
        self.seqnames = [r.seqname for r in self.graphic_records]

    def init_figure(self):
        self.figure, self.axs = plt.subplots(
            nrows=len(self.graphic_records),
            ncols=1
        )
        if len(self.graphic_records) == 1:
            self.axs = [self.axs]

    def plot_graphic_records(self):
        x_max = self.max_seq_len()
        for i, record in enumerate(self.graphic_records):
            is_last = i == len(self.graphic_records) - 1
            record.plot(
                ax=self.axs[i],
                figure_width=None,
                draw_line=True,
                with_ruler=is_last,
                ruler_color=None,
                plot_sequence=False,
                annotate_inline=True,
                max_label_length=50,
                max_line_length=30,
                level_offset=0,
                strand_in_label_threshold='default',
                elevate_outline_annotations='default',
                x_lim=[0, x_max],
                figure_height=None,
                sequence_params=None)

    def config_figure(self):

        self.annotate_seqname()
        self.set_ylim()
        self.set_scale_bar()
        self.set_figure_size()
        self.figure.tight_layout()

    def annotate_seqname(self):
        for ax, seqname in zip(self.axs, self.seqnames):
            ax.annotate(
                seqname + '  ',
                (0, 0),
                fontsize='small',
                horizontalalignment='right',
                verticalalignment='center')

    def set_ylim(self):
        for ax in self.axs:
            ax.set_ylim(-5, 5)

    def set_scale_bar(self):
        ax = self.axs[-1]
        ax.set_xticks([0, self.SCALE_BAR_KB])
        ax.tick_params(axis='x', which='major', labelsize='small')

    def set_figure_size(self):
        w = self.get_figure_width()
        h = self.HEIGHT_CM_PER_LOCUS * len(self.graphic_records)
        self.figure.set_size_inches(w=w, h=h)

    def get_figure_width(self) -> float:
        locus_width = self.WIDTH_CM_PER_KB * self.max_seq_len() / 1000
        seqname_width = self.WIDTH_CM_PER_CHAR * self.max_seqname_len()
        return locus_width + seqname_width

    def max_seq_len(self) -> int:
        return max([r.sequence_length for r in self.graphic_records])

    def max_seqname_len(self) -> int:
        return max(map(len, self.seqnames))

    def save_output(self):
        for fmt in ['pdf', 'png']:
            self.figure.savefig(f'{self.output}.{fmt}', dpi=self.dpi)


def split_list(ls: list, size: int) -> List[list]:
    quotient = len(ls) // size
    remainder = len(ls) % size

    num_chunks = quotient if remainder == 0 else quotient + 1

    ret = []
    for i in range(num_chunks):
        chunk = ls[i*size:i*size+size]
        ret.append(chunk)
    return ret
