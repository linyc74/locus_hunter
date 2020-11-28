import os
import subprocess
from typing import List, Union
from .template import Processor, Settings


class Caller(Processor):

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def call(self, cmd_or_args: Union[str, List[str]]):

        if type(cmd_or_args) is str:
            cmd = cmd_or_args
        else:
            cmd = ' '.join(map(str, cmd_or_args))

        if self.debug:
            self.logger.info(f'CMD: {cmd}')

        subprocess.check_call(cmd, shell=True)


def get_temp_path(prefix: str = 'temp', suffix: str = '') -> str:
    i = 0
    while True:
        path = f'{prefix}_{i:06d}{suffix}'
        if not os.path.exists(path):
            break
        i += 1
    return path
