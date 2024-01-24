import argparse
import locus_hunter


__VERSION__ = '1.1.1'


PROG = 'python locus_hunter'
DESCRIPTION = f'Locus Hunter (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['-q', '--query-faa'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the query fasta file',
        }
    },
    {
        'keys': ['-g', '--gbk-dir'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the genbank folder',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-e', '--evalue'],
        'properties': {
            'type': float,
            'required': False,
            'default': 1e-20,
            'help': 'E-value for blastp (default: %(default)s)',
        }
    },
    {
        'keys': ['-x', '--extension'],
        'properties': {
            'type': int,
            'required': False,
            'default': 5000,
            'help': 'number of base pair extended from each CDS hit by blastp (default: %(default)s)',
        }
    },
    {
        'keys': ['--min-hits-per-locus'],
        'properties': {
            'type': int,
            'required': False,
            'default': 1,
            'help': 'min number of blastp hits per locus (default: %(default)s)',
        }
    },
    {
        'keys': ['--ortholog-identity'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.9,
            'help': 'identity (between 0 and 1) of orthologous genes (default: %(default)s)',
        }
    },
    {
        'keys': ['--dereplicate-loci'],
        'properties': {
            'action': 'store_true',
            'help': 'dereplicate loci with identical sequence of orthologous genes',
        }
    },
    {
        'keys': ['--include-locus-names'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'comma-separated locus names to be included during dereplication (default: %(default)s)',
        }
    },
    {
        'keys': ['--label-attributes'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'gene,locus_tag',
            'help': 'comma-separated attributes to be labeled on each CDS feature (default: %(default)s)',
        }
    },
    {
        'keys': ['--loci-per-plot'],
        'properties': {
            'type': int,
            'required': False,
            'default': 100,
            'help': 'number of loci to be plotted in each output image (default: %(default)s)',
        }
    },
    {
        'keys': ['--dpi'],
        'properties': {
            'type': int,
            'required': False,
            'default': 300,
            'help': 'image resolution (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--output'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'output',
            'help': 'name of output files (default: %(default)s)',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __VERSION__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running Locus Hunter version {__VERSION__}\n', flush=True)
        locus_hunter.main(
            query_faa=args.query_faa,
            gbk_dir=args.gbk_dir,
            evalue=args.evalue,
            extension=args.extension,
            min_hits_per_locus=args.min_hits_per_locus,
            ortholog_identity=args.ortholog_identity,
            dereplicate_loci=args.dereplicate_loci,
            include_locus_names=args.include_locus_names,
            label_attributes=args.label_attributes,
            loci_per_plot=args.loci_per_plot,
            dpi=args.dpi,
            output=args.output,
            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()
