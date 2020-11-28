import random
import pandas as pd
from ngslite import Chromosome
from typing import Dict, Any, List
from .constant import COLOR_KEY, ORTHOLOG_ID_KEY
from .template import Processor, Settings


def sort_by_abundance(items: List[Any]) -> List[Any]:

    item_to_count = {}
    for k in items:
        item_to_count.setdefault(k, 0)
        item_to_count[k] += 1

    def count(item):
        return item_to_count[item]

    return sorted(items, key=count, reverse=True)


def unique(s: List[Any]) -> List[Any]:
    return list(pd.Series(s).unique())


def get_random_color(min_: int, max_: int) -> str:
    assert min_ >= 0
    assert max_ <= 255
    assert min_ <= max_

    ret = '#'
    for _ in range(3):
        c = random.randint(min_, max_)
        ret += f'{c:02X}'
    return ret


class ColorDictGenerator:

    keys: List[Any]

    COLORS = [
        '#1F77B4',
        '#FF7F0E',
        '#2CA02C',
        '#D62728',
        '#9467BD',
        '#8C564B',
        '#E377C2',
        '#7F7F7F',
        '#BCBD22',
        '#17BECF',
    ]

    color_range = (128, 255)

    def main(self, keys: List[Any]) -> Dict[Any, str]:
        self.set_keys(non_unique_keys=keys)
        ret = self.get_color_dict()
        return ret

    def set_keys(self, non_unique_keys: List[Any]):
        self.keys = unique(sort_by_abundance(non_unique_keys))

    def get_color_dict(self):
        ret = {}
        for i, key in enumerate(self.keys):
            if i < len(self.COLORS):
                ret[key] = self.COLORS[i]
            else:
                a, b = self.color_range
                ret[key] = get_random_color(min_=a, max_=b)
        return ret


class AddColor(Processor):

    loci: List[Chromosome]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            loci: List[Chromosome]) -> List[Chromosome]:

        self.loci = loci

        ortholog_ids = self.get_ortholog_ids()
        color_dict = self.get_color_dict(ortholog_ids=ortholog_ids)
        self.add_color(color_dict=color_dict)

        return self.loci

    def get_ortholog_ids(self) -> List[str]:
        ret = []
        for locus in self.loci:
            for feature in locus.features:
                ortholog_id = feature.get_attribute(key=ORTHOLOG_ID_KEY)
                if ortholog_id is not None:
                    ret.append(ortholog_id)
        return ret

    def get_color_dict(self, ortholog_ids: List[str]) -> Dict[str, str]:
        return ColorDictGenerator().main(keys=ortholog_ids)

    def add_color(self, color_dict: Dict[str, str]):
        for locus in self.loci:
            for feature in locus.features:
                ortholog_id = feature.get_attribute(key=ORTHOLOG_ID_KEY)
                color = color_dict.get(ortholog_id, None)
                if color is not None:
                    feature.add_attribute(
                        key=COLOR_KEY,
                        val=color)
