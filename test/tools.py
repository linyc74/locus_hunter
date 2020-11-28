import os
import unittest
import pandas as pd
from typing import Tuple


def setup_dir(__file__: str) -> Tuple[str, str, str]:

    indir = __file__[:-3]

    basedir = os.path.dirname(__file__)
    workdir = os.path.join(basedir, 'workdir')
    outdir = os.path.join(basedir, 'outdir')

    os.makedirs(workdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    return indir, workdir, outdir


class TestCase(unittest.TestCase):

    def assertFileEqual(self, file1: str, file2: str):
        with open(file1) as fh1:
            with open(file2) as fh2:
                self.assertEqual(fh1.read(), fh2.read())

    def assertDataFrameEqual(self, df1: pd.DataFrame, df2: pd.DataFrame):
        self.assertListEqual(list(df1.columns), list(df2.columns))
        self.assertListEqual(list(df1.index), list(df2.index))
        for c in df1.columns:
            for i in df1.index:
                item1 = df1.loc[i, c]
                item2 = df2.loc[i, c]
                if pd.isna(item1) and pd.isna(item2):
                    continue
                self.assertAlmostEqual(item1, item2)


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
