import pandas as pd
from ngslite import write_fasta
from typing import List, Tuple, Union
from .tools import get_temp_path
from .template import Processor, Settings


FAA_DATA_TYPE = List[Tuple[str, str]]


class Blastp(Processor):

    query: Union[str, FAA_DATA_TYPE]
    library: Union[str, FAA_DATA_TYPE]
    evalue: float

    query_faa: str
    library_faa: str
    db: str

    output_columns = [
        'query',
        'subject',
        'percent_id',
        'length',
        'mismatch',
        'gapopen',
        'qstart',
        'qend',
        'sstart',
        'send',
        'evalue',
        'bitscore',
    ]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            query: Union[str, FAA_DATA_TYPE],
            library: Union[str, FAA_DATA_TYPE],
            evalue: float) -> pd.DataFrame:

        self.query = query
        self.library = library
        self.evalue = evalue

        self.set_query_faa()
        self.set_library_faa()
        self.set_db()
        df = self.run_blastp()

        return df

    def set_query_faa(self):
        if type(self.query) is str:
            self.query_faa = self.query
        else:
            faa = get_temp_path(prefix=f'{self.workdir}/query', suffix='.faa')
            write_fasta(data=self.query, file=faa)
            self.query_faa = faa

    def set_library_faa(self):
        if type(self.library) is str:
            self.library_faa = self.library
        else:
            faa = get_temp_path(prefix=f'{self.workdir}/library', suffix='.faa')
            write_fasta(data=self.library, file=faa)
            self.library_faa = faa

    def set_db(self):
        logfile = get_temp_path(prefix=f'{self.workdir}/makeblastdb_log')
        self.db = get_temp_path(prefix=f'{self.workdir}/blastp_db')
        lines = [
            'makeblastdb',
            f'-in {self.library_faa}',
            '-dbtype prot',
            f'-logfile {logfile}',
            f'-out {self.db}',
        ]
        cmd = self.CMD_LINEBREAK.join(lines)
        self.call(cmd)

    def run_blastp(self) -> pd.DataFrame:
        blastp_output = get_temp_path(f'{self.workdir}/blastp', '.tsv')

        lines = [
            'blastp',
            f'-query {self.query_faa}',
            f'-db {self.db}',
            f'-evalue {self.evalue}',
            f'-outfmt 6',
            f'-num_threads {self.threads}',
            f'-out {blastp_output}',
        ]
        cmd = self.CMD_LINEBREAK.join(lines)
        self.call(cmd)

        df = pd.read_csv(
            blastp_output, sep='\t', names=self.output_columns)

        return df
