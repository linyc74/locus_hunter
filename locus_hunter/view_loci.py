import matplotlib.pyplot as plt
from typing import Optional, List
from ngslite import GenericFeature, Chromosome
from dna_features_viewer import GraphicFeature, GraphicRecord
from .constant import COLOR_KEY
from .template import Processor, Settings


class ViewLoci(Processor):

    loci: List[Chromosome]
    output: str

    graphic_records: List[GraphicRecord]

    def main(
            self,
            loci: List[Chromosome],
            output: str):

        self.loci = loci
        self.output = output

        self.set_graphic_records()
        self.plot_graphic_records()

    def set_graphic_records(self):
        to_graphic_record = ChromosomeToGraphicRecord(self.settings).main
        self.graphic_records = [to_graphic_record(locus) for locus in self.loci]

    def plot_graphic_records(self):
        seqnames = [locus.seqname for locus in self.loci]
        PlotGraphicRecords(self.settings).main(
            graphic_records=self.graphic_records,
            seqnames=seqnames,
            output=self.output)


class ChromosomeToGraphicRecord(Processor):

    FEATURE_TYPES = ['CDS']

    chromosome: Chromosome

    graphic_features: List[GraphicFeature]

    def main(self, chromosome: Chromosome) -> GraphicRecord:
        self.chromosome = chromosome

        self.set_graphic_features()

        return GraphicRecord(
            sequence_length=len(self.chromosome.sequence),
            sequence=None,
            features=self.graphic_features,
            feature_level_height=1,
            first_index=0,
            plots_indexing='biopython',
            labels_spacing=8,
            ticks_resolution='auto')

    def set_graphic_features(self):
        to_graphic_feature = GenericToGraphicFeature(self.settings).main
        self.graphic_features = [
            to_graphic_feature(f) for f in self.chromosome.features
            if f.type in self.FEATURE_TYPES
        ]


class GenericToGraphicFeature(Processor):

    STRAND_DICT = {
        '+': 1,
        '-': -1
    }
    DEFAULT_COLOR = '#FFFFFF'
    LABEL = None
    THICKNESS = 8.0
    LINEWIDTH = 0.5
    LINECOLOR = '#000000'
    OPEN_LEFT = False
    OPEN_RIGHT = False

    generic_feature: GenericFeature

    graphic_feature: GraphicFeature

    def main(self, generic_feature: GenericFeature) -> GraphicFeature:
        self.generic_feature = generic_feature
        self.set_graphic_feature()
        return self.graphic_feature

    def set_graphic_feature(self):
        self.graphic_feature = GraphicFeature(
            start=self.generic_feature.start,
            end=self.generic_feature.end,
            strand=self.get_strand(),
            label=self.LABEL,
            color=self.get_color(),
            thickness=self.THICKNESS,
            linewidth=self.LINEWIDTH,
            linecolor=self.LINECOLOR,
            fontdict=None,
            html=None,
            open_left=self.OPEN_LEFT,
            open_right=self.OPEN_RIGHT,
            box_linewidth=1,
            box_color='auto',
            legend_text=None,
            label_link_color='black')

    def get_strand(self) -> int:
        f = self.generic_feature
        return self.STRAND_DICT.get(f.strand, 0)

    def get_color(self) -> str:
        f = self.generic_feature
        color = f.get_attribute(key=COLOR_KEY)
        if color is None:
            color = self.DEFAULT_COLOR
        return color


class PlotGraphicRecords(Processor):

    HEIGHT_CM_PER_LOCUS = 1 / 2.54
    WIDTH_CM_PER_KB = 1 / 2.54
    WIDTH_CM_PER_CHAR = 0.2 / 2.54
    SCALE_BAR_KB = 5000
    DPI = 600

    graphic_records: List[GraphicRecord]
    seqnames: List[str]
    output: str

    figure: plt.Figure
    axs = List[plt.Axes]

    def main(
            self,
            graphic_records: List[GraphicRecord],
            seqnames: List[str],
            output: str):

        self.graphic_records = graphic_records
        self.seqnames = seqnames
        self.output = output

        self.init_figure()
        self.plot_graphic_records()
        self.config_figure()
        self.save_output()

    def init_figure(self):
        self.figure, self.axs = plt.subplots(
            nrows=len(self.graphic_records),
            ncols=1
        )

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
            ax.set_ylim(-0.5, 0.5)

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
        for fmt in ['pdf']:
            self.figure.savefig(f'{self.output}.{fmt}', dpi=self.DPI)
