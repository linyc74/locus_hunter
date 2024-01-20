import os
import shutil
import unittest
import pandas as pd
from typing import Tuple
from locus_hunter.template import Settings


def get_dirs(py_path: str) -> Tuple[str, str, str]:
    indir = os.path.relpath(path=py_path[:-3], start=os.getcwd())
    basedir = os.path.dirname(indir)
    workdir = os.path.join(basedir, 'workdir')
    outdir = os.path.join(basedir, 'outdir')
    return indir, workdir, outdir


class TestCase(unittest.TestCase):

    def set_up(self, py_path: str):
        self.indir, self.workdir, self.outdir = get_dirs(py_path=py_path)

        for d in [self.workdir, self.outdir]:
            os.makedirs(d, exist_ok=True)

        self.settings = Settings(
            workdir=self.workdir,
            outdir=self.outdir,
            threads=6,
            debug=True)

    def tear_down(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def assertFileEqual(self, first: str, second: str):
        with open(first) as fh1:
            with open(second) as fh2:
                self.assertEqual(fh1.read(), fh2.read())

    def assertDataFrameEqual(self, first: pd.DataFrame, second: pd.DataFrame):
        self.assertListEqual(list(first.columns), list(second.columns))
        self.assertListEqual(list(first.index), list(second.index))
        for c in first.columns:
            for i in first.index:
                a, b = first.loc[i, c], second.loc[i, c]
                if pd.isna(a) and pd.isna(b):
                    continue
                self.assertAlmostEqual(a, b)

    def assertFileExists(self, expected: str, actual: str):
        self.assertEqual(expected, actual)
        assert os.path.exists(actual)


def remove_genbank_date_str(gbk: str):
    with open(gbk) as fh:
        text = fh.read()

    os.remove(gbk)

    with open(gbk, 'w') as fh:
        for line in text.splitlines():
            if line.startswith('LOCUS'):
                pos = line.rfind(' ')
                line = line[:pos]
            fh.write(line + '\n')
