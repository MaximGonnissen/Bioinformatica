from abc import ABC, abstractmethod
from itertools import combinations
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

    def scoring_function(self, *args) -> Union[int, float]:
        """
        Scoring function for the MSA solver.
        :param args: Coordinates in the scoring matrix.
        :return: Score for the given coordinates.
        """
        # Simple version for now, only taking into account matches and mismatches.
        indices = args
        if len(indices) != len(self.scoring_matrix.shape):
            raise IndexError("Incorrect number of indices given.")

        # Get the indices for the sequences.
        sequence_chars = [self.scoring_matrix.sequences[i][indices[i] - 1] for i in range(len(indices))]

        score = 0
        char_combinations = combinations(sequence_chars, 2)
        for char_combination in char_combinations:
            score += self.substitution_matrix[char_combination[0]][char_combination[1]]

        return score

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
    def traceback(self, *args) -> Tuple[Union[int, float], List[str]]:
        """
        Perform the traceback.
        :param args: Coordinates to start the traceback from.
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
