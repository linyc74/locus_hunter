import os
import shutil
from .locus_hunter import LocusHunter
from .template import Settings
from .tools import get_temp_path


def main(
        query_faa: str,
        gbk_dir: str,
        evalue: float,
        extension: int,
        min_hits_per_locus: int,
        ortholog_identity: float,
        dereplicate_loci: bool,
        include_locus_names: str,
        label_attributes: str,
        loci_per_plot: int,
        dpi: int,
        output: str,
        threads: int,
        debug: bool):

    workdir = get_temp_path(prefix='locus_hunter')

    settings = Settings(
        workdir=workdir,
        outdir='.',
        threads=threads,
        debug=debug)

    os.makedirs(workdir)

    LocusHunter(settings).main(
        query_faa=query_faa,
        gbk_dir=gbk_dir,
        evalue=evalue,
        extension=extension,
        min_hits_per_locus=min_hits_per_locus,
        ortholog_identity=ortholog_identity,
        dereplicate_loci=dereplicate_loci,
        include_locus_names=include_locus_names.split(','),
        label_attributes=label_attributes.split(','),
        loci_per_plot=loci_per_plot,
        dpi=dpi,
        output=output)

    if not settings.debug:
        shutil.rmtree(workdir)
