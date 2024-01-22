import subprocess
from typing import Any
from datetime import datetime


class Settings:

    workdir: str
    outdir: str
    threads: int
    debug: bool

    def __init__(
            self,
            workdir: str,
            outdir: str,
            threads: int,
            debug: bool):

        self.workdir = workdir
        self.outdir = outdir
        self.threads = threads
        self.debug = debug


class Logger:

    INFO = 'INFO'
    DEBUG = 'DEBUG'

    name: str
    level: str

    def __init__(self, name: str, level: str):
        self.name = name
        assert level in [self.INFO, self.DEBUG]
        self.level = level

    def info(self, msg: Any):
        self.__print(msg_level=self.INFO, msg=msg)

    def debug(self, msg: Any):
        if self.level == self.DEBUG:
            self.__print(msg_level=self.DEBUG, msg=msg)

    def __print(self, msg_level: str, msg: Any):
        if self.level == self.DEBUG:
            print(f'{self.name}\t{msg_level}\t{datetime.now()}\n{msg}\n', flush=True)
        else:
            print(f'{msg}', flush=True)


class Processor:

    CMD_LINEBREAK = ' \\\n  '

    settings: Settings
    workdir: str
    outdir: str
    threads: int
    debug: bool

    logger: Logger

    def __init__(self, settings: Settings):

        self.settings = settings
        self.workdir = settings.workdir
        self.outdir = settings.outdir
        self.threads = settings.threads
        self.debug = settings.debug

        self.logger = Logger(
            name=self.__class__.__name__,
            level=Logger.DEBUG if self.debug else Logger.INFO
        )

    def call(self, cmd: str):
        self.logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)
