class Settings:

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

    def __init__(self, name: str):
        self.name = name

    def info(self, msg: str):
        print(msg, flush=True)


class Processor:

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
        self.logger = Logger(name=self.__class__.__name__)
