import matplotlib.pyplot as plt
from typing import Optional, List
from ngslite import GenericFeature, Chromosome
from dna_features_viewer import GraphicFeature, GraphicRecord
from .template import Processor, Settings
from .constant import COLOR_KEY


class GenericToGraphicFeature(Processor):

    generic_feature: GenericFeature

    strand_dict = {
        '+': 1,
        '-': -1
    }

    color_key: str = COLOR_KEY
    default_color: str = '#FFFFFF'

    label: Optional[str] = None
    thickness: float = 8.
    linewidth: float = .5
    linecolor: str = '#000000'
    open_left: bool = False
    open_right: bool = False

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, generic_feature: GenericFeature) -> GraphicFeature:
        self.generic_feature = generic_feature
        return self.get_graphic_feature()

    def get_graphic_feature(self) -> GraphicFeature:
        f = self.generic_feature
        return GraphicFeature(
            start=f.start,
            end=f.end,
            strand=self.get_strand(),
            label=self.label,
            color=self.get_color(),
            thickness=self.thickness,
            linewidth=self.linewidth,
            linecolor=self.linecolor,
            fontdict=None,
            html=None,
            open_left=self.open_left,
            open_right=self.open_right,
            box_linewidth=1,
            box_color='auto',
            legend_text=None,
            label_link_color='black')

    def get_strand(self) -> int:
        f = self.generic_feature
        return self.strand_dict.get(f.strand, 0)

    def get_color(self) -> str:
        f = self.generic_feature
        color = f.get_attribute(key=COLOR_KEY)
        if color is None:
            color = self.default_color
        return color


class ChromosomeToGraphicRecord(Processor):

    chromosome: Chromosome

    feature_types = ['CDS']

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.to_graphic_feature = GenericToGraphicFeature(self.settings).main

    def main(self, chromosome: Chromosome) -> GraphicRecord:
        self.chromosome = chromosome

        return GraphicRecord(
            sequence_length=len(self.chromosome.sequence),
            sequence=None,
            features=self.get_graphic_features(),
            feature_level_height=1,
            first_index=0,
            plots_indexing='biopython',
            labels_spacing=8,
            ticks_resolution='auto')

    def get_graphic_features(self) -> List[GraphicFeature]:
        generic_features = self.chromosome.features
        return [self.to_graphic_feature(f) for f in generic_features
                if f.type in self.feature_types]


def to_inch(cm) -> float:
    return cm / 2.54


class PlotGraphicRecords(Processor):

    graphic_records: List[GraphicRecord]
    seqnames: List[str]
    output: str

    height_cm_per_locus: float = to_inch(1)
    width_cm_per_kb: float = to_inch(1)
    width_cm_per_char: float = to_inch(0.2)
    scale_bar_kb: int = 5000
    dpi = 600

    figure: plt.Figure
    axs = List[plt.Axes]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.to_graphic_record = ChromosomeToGraphicRecord(self.settings).main

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
        ax.set_xticks([0, self.scale_bar_kb])
        ax.tick_params(axis='x', which='major', labelsize='small')

    def set_figure_size(self):
        w = self.get_figure_width()
        h = self.height_cm_per_locus * len(self.graphic_records)
        self.figure.set_size_inches(w=w, h=h)

    def get_figure_width(self) -> float:
        locus_width = self.width_cm_per_kb * self.max_seq_len() / 1000
        seqname_width = self.width_cm_per_char * self.max_seqname_len()
        return locus_width + seqname_width

    def max_seq_len(self) -> int:
        return max([r.sequence_length for r in self.graphic_records])

    def max_seqname_len(self) -> int:
        return max(map(len, self.seqnames))

    def save_output(self):
        for fmt in ['pdf']:
            self.figure.savefig(f'{self.output}.{fmt}', dpi=self.dpi)


class ViewLoci(Processor):

    loci: List[Chromosome]
    output: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.to_graphic_record = ChromosomeToGraphicRecord(self.settings).main

    def main(
            self,
            loci: List[Chromosome],
            output: str):

        self.loci = loci
        self.output = output

        graphic_records = self.get_graphic_records()
        self.plot(graphic_records)

    def get_graphic_records(self) -> List[GraphicRecord]:
        return [self.to_graphic_record(locus) for locus in self.loci]

    def plot(self, graphic_records: List[GraphicRecord]):
        seqnames = [locus.seqname for locus in self.loci]
        PlotGraphicRecords(self.settings).main(
            graphic_records=graphic_records,
            seqnames=seqnames,
            output=self.output)
