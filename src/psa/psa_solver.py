from abc import ABC, abstractmethod
from typing import List, Tuple, Union, Optional, Callable

from blosum import BLOSUM

from psa.scoring_matrix import ScoringMatrix


class PSASolver(ABC):
    """
    Abstract class for pairwise sequence alignment solvers.
    """

    def __init__(self, config: dict, substitution_matrix: Optional[dict], *args, **kwargs):
        """
        Initialise the PSA solver.
        :param config: Configuration for the PSA solver.
        :param substitution_matrix: Substitution matrix to use for the PSA solver. Defaults to BLOSUM62.
        """
        super().__init__(*args, **kwargs)
        self.config = config
        self.substitution_matrix = substitution_matrix or BLOSUM(62)
        self.scoring_matrix: ScoringMatrix = None
        self.sequence_1: str = None
        self.sequence_2: str = None

    @property
    @abstractmethod
    def scoring_matrix_cls(self) -> Callable:
        """
        Scoring matrix class to use for the PSA solver.
        """
        pass

    def scoring_function(self, item_1: str, item_2: str) -> Union[int, float]:
        """
        Scoring function for the pairwise sequence alignment problem.
        :param item_1: The first item.
        :param item_2: The second item.
        :return: The score for the two items.
        """
        return self.substitution_matrix[item_1][item_2]

    def gap_penalty(self, length: int = 1) -> Union[int, float]:
        """
        Calculate the gap penalty for a given length.
        :param length: Length of the gap.
        :return: Gap penalty.
        """
        return self.config['indel'] * length

    def solve(self, sequence_1: str, sequence_2: str) -> Tuple[int, List[Tuple[str, str]]]:
        """
        Solve the pairwise sequence alignment problem.
        :param sequence_1: The first sequence.
        :param sequence_2: The second sequence.
        :return: The score and all valid alignments.
        """
        self.sequence_1 = sequence_1
        self.sequence_2 = sequence_2
        self.scoring_matrix = self.scoring_matrix_cls(sequence_1, sequence_2)
        self.pre_solve()
        results = self._solve()
        return self.post_solve(results)

    def pre_solve(self):
        """
        Pre-solve hook.
        """
        pass

    def post_solve(self, results: Tuple[int, List[Tuple[str, str]]]):
        """
        Post-solve hook.
        :param results: The results of the PSA solver.
        """
        return results

    @abstractmethod
    def _solve(self) -> Tuple[int, List[Tuple[str, str]]]:
        """
        Solve the pairwise sequence alignment problem.
        :return: The score and all valid alignments.
        """
        pass
