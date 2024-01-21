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
    dereplicate_loci: bool
    include_locus_names: List[str]

    def main(
            self,
            loci: List[Chromosome],
            ortholog_identity: float,
            dereplicate_loci: bool,
            include_locus_names: List[str]) -> List[Chromosome]:

        self.loci = loci
        self.ortholog_identity = ortholog_identity
        self.dereplicate_loci = dereplicate_loci
        self.include_locus_names = include_locus_names

        self.assign_ortholog_id()
        self.dereplicate_loci_by_identical_ortholog_ids()
        self.sort_loci_by_comparison()

        return self.loci

    def assign_ortholog_id(self):
        self.loci = AssignOrthologId(self.settings).main(
            loci=self.loci,
            ortholog_identity=self.ortholog_identity)

    def dereplicate_loci_by_identical_ortholog_ids(self):
        if self.dereplicate_loci:
            self.loci = DereplicateLociByIdenticalOrthologIDs(self.settings).main(
                loci=self.loci,
                include_locus_names=self.include_locus_names)

    def sort_loci_by_comparison(self):
        self.loci = SortLociByComparison(self.settings).main(
            loci=self.loci)


class AssignOrthologId(Processor):

    loci: List[Chromosome]
    ortholog_identity: float

    faa: str
    cds_id_to_ortholog_id: Dict[str, int]

    def main(
            self,
            loci: List[Chromosome],
            ortholog_identity: float) -> List[Chromosome]:

        self.loci = loci
        self.ortholog_identity = ortholog_identity

        self.write_faa()
        self.create_orthologs_by_cd_hit()
        self.add_ortholog_id_to_features()

        return self.loci

    def write_faa(self):

        self.faa = f'{self.workdir}/loci.faa'

        with FastaWriter(self.faa) as writer:
            for locus in self.loci:
                for feature in locus.features:
                    cds_id = feature.get_attribute(key=CDS_ID_KEY)
                    seq = feature.get_attribute(key='translation')
                    if seq is not None:
                        writer.write(header=cds_id, sequence=seq)

        return self.faa

    def create_orthologs_by_cd_hit(self):
        self.cds_id_to_ortholog_id = CdHit(self.settings).main(
            faa=self.faa,
            sequence_identity=self.ortholog_identity)

    def add_ortholog_id_to_features(self):
        for locus in self.loci:
            for feature in locus.features:
                cds_id = feature.get_attribute(CDS_ID_KEY)
                ortholog_id = self.cds_id_to_ortholog_id.get(cds_id, None)
                if ortholog_id is not None:
                    feature.add_attribute(
                        key=ORTHOLOG_ID_KEY,
                        val=ortholog_id)


class DereplicateLociByIdenticalOrthologIDs(Processor):

    loci: List[Chromosome]
    include_locus_names: List[str]

    ortholog_ids_to_loci: Dict[str, List[Chromosome]]

    dereplicated_loci: List[Chromosome]

    def main(
            self,
            loci: List[Chromosome],
            include_locus_names: List[str]) -> List[Chromosome]:

        self.loci = loci
        self.include_locus_names = include_locus_names

        self.set_ortholog_ids_to_loci()
        self.pick_loci()
        self.log_info()

        return self.dereplicated_loci

    def set_ortholog_ids_to_loci(self):
        self.ortholog_ids_to_loci = {}
        for locus in self.loci:
            ortholog_ids = locus_to_ortholog_ids(locus)
            self.ortholog_ids_to_loci.setdefault(ortholog_ids, []).append(locus)

    def pick_loci(self):
        self.dereplicated_loci = []
        for loci in self.ortholog_ids_to_loci.values():

            self.dereplicated_loci.append(loci[0])  # always include the first locus

            for locus in loci[1:]:  # after the first locus, look for the locus names to be included
                if self.__should_be_included(locus=locus):
                    self.dereplicated_loci.append(locus)

    def __should_be_included(self, locus: Chromosome) -> bool:
        for locus_name in self.include_locus_names:
            if locus_name in locus.seqname:
                return True
        return False

    def log_info(self):
        msg = f'Dereplicated {len(self.loci)} loci to {len(self.dereplicated_loci)} loci'
        self.logger.info(msg)


def locus_to_ortholog_ids(chromosome: Chromosome) -> str:
    ortholog_ids = []
    for feature in chromosome.features:
        orid = str(feature.get_attribute(ORTHOLOG_ID_KEY))  # can be None, so convert to str
        ortholog_ids.append(orid)
    return ','.join(ortholog_ids)


class SortLociByComparison(Processor):

    LINKAGE_METHOD = 'average'

    loci: List[Chromosome]

    locus_id_to_ortholog_ids: Dict[str, List[int]]
    distance_matrix: np.ndarray
    linkage_matrix: np.ndarray
    idx_order: List[int]
    sorted_loci: List[Chromosome]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.compare = SmithWatermanAligner().run

    def main(
            self,
            loci: List[Chromosome]) -> List[Chromosome]:

        self.loci = loci

        self.set_locus_id_to_ortholog_ids()
        self.set_condensed_distance_matrix()
        self.set_linkage_matrix()
        self.plot_tree_and_set_idx_order()
        self.set_sorted_loci()

        return self.sorted_loci

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

    def set_condensed_distance_matrix(self):

        # condensed distance matrix:
        #   1-D array of pairs following the order given by combinations()
        similarity_matrix = []

        for locus1, locus2 in combinations(self.loci, 2):

            orthologs_1 = self.locus_id_to_ortholog_ids[locus1.seqname]
            orthologs_2 = self.locus_id_to_ortholog_ids[locus2.seqname]

            s = self.compare(list1=orthologs_1, list2=orthologs_2)

            similarity_matrix.append(s)

        self.distance_matrix = np.exp(-np.array(similarity_matrix))

    def set_linkage_matrix(self):
        self.linkage_matrix = hierarchy.linkage(self.distance_matrix, self.LINKAGE_METHOD)

    def plot_tree_and_set_idx_order(self):
        idx_order = hierarchy.dendrogram(
            Z=self.linkage_matrix,
            orientation='top'  # root on top
        )['ivl']  # 'ivl': a list of labels corresponding to the leaf nodes

        self.idx_order = list(map(int, idx_order))

    def set_sorted_loci(self):
        self.sorted_loci = [self.loci[i] for i in self.idx_order]


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

        M = np.zeros((len(L1) + 1, len(L2) + 1), dtype=np.int32)

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
