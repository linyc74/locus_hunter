import argparse
import locus_hunter


__version__ = '1.0.5'


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        args = self.parser.parse_args()
        locus_hunter.main(
            query_faa=args.query_faa,
            gbk_dir=args.gbk_dir,
            evalue=args.evalue,
            extension=args.extension,
            ortholog_identity=args.ortholog_identity,
            output=args.output,
            threads=args.threads,
            debug=args.debug)

    def set_parser(self):
        prog = 'python locus_hunter'

        description = f'Locus Hunter (v{__version__}) by Yu-Cheng Lin (yclin.python@gmail.com)'

        dependencies = '\n  '.join([
            'python (>=3.6)',
            'numpy (1.17.2)',
            'pandas (0.25.1)',
            'ngslite (1.1.2)',
            'dna_features_viewer (3.0.1)',
            'blastp (2.5.0+)',
            'CD-HIT (4.8.1)',
        ])

        epilog = f'dependency:\n  {dependencies}'

        self.parser = argparse.ArgumentParser(
            prog=prog,
            description=description,
            epilog=epilog,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')

        group.add_argument(
            '-q', '--query-faa', type=str, required=True,
            help='path to the query fasta file')

        group.add_argument(
            '-g', '--gbk-dir', type=str, required=True,
            help='path to the genbank folder')

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        default = '(default: %(default)s)'

        group.add_argument(
            '-e', '--evalue', type=float, required=False, default=1e-20,
            help=f'E-value for blastp {default}')

        group.add_argument(
            '-x', '--extension', type=int, required=False, default=5000,
            help=f'number of bp extended from each CDS hit {default}')

        group.add_argument(
            '-r', '--ortholog-identity', type=float, required=False, default=0.9,
            help=f'identity (between 0 and 1) of orthologous genes {default}')

        group.add_argument(
            '-o', '--output', type=str, required=False, default='output',
            help=f'name of output files {default}')

        group.add_argument(
            '-t', '--threads', type=int, required=False, default=4,
            help=f'number of CPU threads {default}')

        group.add_argument(
            '-d', '--debug', action='store_true', help=f'debug mode')

        group.add_argument(
            '-h', '--help', action='help',
            # default=SUPPRESS,
            help='show this help message and exit')


if __name__ == '__main__':
    EntryPoint().main()
