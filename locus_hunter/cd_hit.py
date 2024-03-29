from typing import Dict
from .tools import get_temp_path
from .template import Processor, Settings


class CdHit(Processor):

    faa: str
    sequence_identity: float

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            faa: str,
            sequence_identity: float) -> Dict[str, int]:

        self.faa = faa
        self.sequence_identity = sequence_identity

        clstr = self.run_cd_hit()
        protein_to_cluster_dict = self.read_clstr(file=clstr)

        return protein_to_cluster_dict

    def run_cd_hit(self) -> str:
        output = get_temp_path(prefix=f'{self.workdir}/cd_hit_output')
        lines = [
            'cd-hit',
            f'-i {self.faa}',
            f'-c {self.sequence_identity}',
            '-d 0',
            f'-T {self.threads}',
            f'-o {output}',
            f'1> {self.workdir}/cd-hit.log',
            f'2> {self.workdir}/cd-hit.log',
        ]
        cmd = self.CMD_LINEBREAK.join(lines)
        self.call(cmd)
        return f'{output}.clstr'

    def read_clstr(self, file: str) -> Dict[str, int]:

        ret = {}

        cluster_id = 0
        with open(file) as fh:
            for line in fh:
                if line.startswith('>Cluster'):
                    cluster_id += 1
                else:
                    protein_id = line.split(', >')[1].split('... ')[0]
                    ret[protein_id] = cluster_id

        return ret
