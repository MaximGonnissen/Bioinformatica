from abc import ABC, abstractmethod
from typing import List, Tuple, Optional, Union

from blosum import BLOSUM

from msa.scoring_matrix.scoring_matrix import ScoringMatrix


class MSASolver(ABC):
    """
    Abstract class for multiple sequence alignment solvers.
    """

    def __init__(self, config: dict, substitution_matrix: Optional[dict] = None, *args, **kwargs):
        """
        Initialise the PSA solver.
        :param config: Configuration for the PSA solver.
        :param substitution_matrix: Substitution matrix to use for the PSA solver. Defaults to BLOSUM62.
        """
        super().__init__(*args, **kwargs)
        self.config = config
        self.substitution_matrix = substitution_matrix or BLOSUM(62)
        self.scoring_matrix: ScoringMatrix = None

        if config.get("match") is not None and substitution_matrix is None:
            for key in self.substitution_matrix.keys():
                for other_key in self.substitution_matrix[key].keys():
                    if key == other_key:
                        self.substitution_matrix[key][other_key] = config["match"]

        if config.get("mismatch") is not None and substitution_matrix is None:
            for key in self.substitution_matrix.keys():
                for other_key in self.substitution_matrix[key].keys():
                    if key != other_key:
                        self.substitution_matrix[key][other_key] = config["mismatch"]

    @abstractmethod
    def initialise_scoring_matrix(self, sequences: List[str]) -> ScoringMatrix:
        """
        Initialise the scoring matrix.
        :param sequences: List of sequences to align.
        :return: Initialised scoring matrix.
        """
        pass

    @abstractmethod
    def fill_scoring_matrix(self):
        """
        Fill the scoring matrix.
        """
        pass

    @abstractmethod
    def traceback(self) -> Tuple[Union[int, float], List[str]]:
        """
        Perform the traceback.
        :return: Tuple of the aligned sequences and the alignment score.
        """
        pass

    def pre_solve(self):
        """
        Pre-solve hook.
        """
        pass

    def post_solve(self, results: Tuple[Union[int, float], List[str]]) -> Tuple[Union[int, float], List[str]]:
        """
        Post-solve hook.
        :param results: The results of the PSA solver.
        """
        return results

    def solve(self, sequences: List[str]) -> Tuple[Union[int, float], List[str]]:
        """
        Solve the multiple sequence alignment problem.
        :param sequences: List of sequences to align.
        :return: Tuple of the aligned sequences and the alignment score.
        """
        self.scoring_matrix = self.initialise_scoring_matrix(sequences)
        self.fill_scoring_matrix()
        self.pre_solve()
        return self.post_solve(self.traceback())
