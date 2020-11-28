import numpy as np
from typing import List, Dict, Any
from itertools import combinations
from scipy.cluster import hierarchy
from ngslite import Chromosome, FastaWriter
from .cd_hit import CdHit
from .template import Processor, Settings
from .constant import CDS_ID_KEY, ORTHOLOG_ID_KEY


class SortLoci(Processor):

    loci: List[Chromosome]
    ortholog_identity: float

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            loci: List[Chromosome],
            ortholog_identity: float) -> List[Chromosome]:

        self.loci = loci
        self.ortholog_identity = ortholog_identity

        self.assign_ortholog_id()
        self.sort_loci_by_comparison()

        return self.loci

    def assign_ortholog_id(self):
        self.loci = AssignOrthologId(self.settings).main(
            loci=self.loci,
            ortholog_identity=self.ortholog_identity)

    def sort_loci_by_comparison(self):
        self.loci = SortLociByComparison(self.settings).main(
            loci=self.loci)


class AssignOrthologId(Processor):

    loci: List[Chromosome]
    ortholog_identity: float

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            loci: List[Chromosome],
            ortholog_identity: float) -> List[Chromosome]:

        self.loci = loci
        self.ortholog_identity = ortholog_identity

        faa = self.write_faa()
        cds_ortholog_dict = self.create_orthologs_by_cd_hit(faa=faa)
        self.add_ortholog_id(cds_ortholog_dict=cds_ortholog_dict)

        return self.loci

    def write_faa(self) -> str:

        faa = f'{self.workdir}/loci.faa'

        with FastaWriter(faa) as writer:
            for locus in self.loci:
                for feature in locus.features:
                    cds_id = feature.get_attribute(key=CDS_ID_KEY)
                    seq = feature.get_attribute(key='translation')
                    if seq is not None:
                        writer.write(header=cds_id, sequence=seq)

        return faa

    def create_orthologs_by_cd_hit(self, faa: str) -> Dict[str, int]:
        return CdHit(self.settings).main(
            faa=faa,
            sequence_identity=self.ortholog_identity)

    def add_ortholog_id(self, cds_ortholog_dict: Dict[str, int]):
        for locus in self.loci:
            for feature in locus.features:
                cds_id = feature.get_attribute(CDS_ID_KEY)
                ortholog_id = cds_ortholog_dict.get(cds_id, None)
                if ortholog_id is not None:
                    feature.add_attribute(
                        key=ORTHOLOG_ID_KEY,
                        val=ortholog_id)


class SortLociByComparison(Processor):

    loci: List[Chromosome]

    locus_id_to_ortholog_ids: Dict[str, List[int]]

    linkage_method = 'average'

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.compare = SmithWatermanAligner().run

    def main(
            self,
            loci: List[Chromosome]) -> List[Chromosome]:

        self.loci = loci
        self.set_locus_id_to_ortholog_ids()

        distance_matrix = self.get_condensed_distance_matrix()
        linkage_matrix = self.get_linkage_matrix(distance_matrix)
        idx_order = self.plot_tree_and_get_idx_order(linkage_matrix)
        sorted_loci = self.get_sorted_loci(idx_order)

        return sorted_loci

    def set_locus_id_to_ortholog_ids(self):
        d = {}

        for locus in self.loci:

            orthologs = []
            for feature in locus.features:
                ortholog_id = feature.get_attribute(key=ORTHOLOG_ID_KEY)
                if ortholog_id is not None:
                    orthologs.append(ortholog_id)

            d[locus.seqname] = orthologs

        self.locus_id_to_ortholog_ids = d

    def get_condensed_distance_matrix(self) -> np.ndarray:

        # condensed distance matrix:
        #   1-D array of pairs following the order given by combinations()
        similarity_matrix = []

        for locus1, locus2 in combinations(self.loci, 2):

            orthologs_1 = self.locus_id_to_ortholog_ids[locus1.seqname]
            orthologs_2 = self.locus_id_to_ortholog_ids[locus2.seqname]

            s = self.compare(list1=orthologs_1, list2=orthologs_2)

            similarity_matrix.append(s)

        distance_matrix = np.exp(-np.array(similarity_matrix))

        return distance_matrix

    def get_linkage_matrix(self, distance_matrix: np.ndarray) -> np.ndarray:
        linkage_matrix = hierarchy.linkage(distance_matrix, self.linkage_method)
        return linkage_matrix

    def plot_tree_and_get_idx_order(
            self, linkage_matrix: np.ndarray) -> List[int]:

        idx_order = hierarchy.dendrogram(
            Z=linkage_matrix,
            orientation='top'  # root on top
        )['ivl']  # 'ivl': a list of labels corresponding to the leaf nodes

        idx_order = list(map(int, idx_order))

        return idx_order

    def get_sorted_loci(self, idx_order: List[int]) -> List[Chromosome]:
        return [self.loci[i] for i in idx_order]


class SmithWatermanAligner:

    END_GAP_SCORE = 0.
    GAP_SCORE = -1.
    MATCH_SCORE = 100.
    MISMATCH_SCORE = -1.

    def compare(self, a: Any, b: Any) -> float:
        return self.MATCH_SCORE if a == b else self.MISMATCH_SCORE

    def run(self,
            list1: List[Any],
            list2: List[Any]) -> float:

        L1, L2 = list1, list2

        M = np.zeros((len(L1) + 1, len(L2) + 1), dtype=np.int)

        for r in range(len(L1)):
            M[r + 1, 0] = M[r, 0] + self.END_GAP_SCORE

        for c in range(len(L2)):
            M[0, c + 1] = M[0, c] + self.END_GAP_SCORE

        for r in range(len(L1)):
            for c in range(len(L2)):
                down_right = M[r, c] + self.compare(L1[r], L2[c])
                right = M[r + 1, c] + self.GAP_SCORE
                down = M[r, c + 1] + self.GAP_SCORE

                M[r + 1, c + 1] = max(down_right, right, down)

        return M[-1, -1]
