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

    def get_alignment_chars(self, *args, comparison_indices: Tuple[int, ...]) -> Tuple[str, ...]:
        """
        Get the characters to align.
        :param args: Coordinates in the scoring matrix.
        :param comparison_indices: Indices of the sequences to compare.
        :return: Tuple of characters to align.
        """
        indices_aligned = tuple([bool(args[i] - comparison_indices[i]) for i in range(len(args))])

        # Get the indices for the sequences.
        sequence_chars = [self.scoring_matrix.sequences[i][args[i] - 1] if indices_aligned[i] else "-" for i in
                          range(len(args))]

        return tuple(sequence_chars)

    def scoring_function(self, *args, comparison_indices: Tuple[int, ...]) -> Union[int, float]:
        """
        Scoring function for the MSA solver.
        :param args: Coordinates in the scoring matrix.
        :param comparison_indices: Indices of the sequences to compare.
        :return: Pairwise score for the given coordinates.
        """
        indices = args
        if len(indices) != len(self.scoring_matrix.shape):
            raise IndexError("Incorrect number of indices given.")

        sequence_chars = self.get_alignment_chars(*args, comparison_indices=comparison_indices)

        score = 0
        char_combinations = combinations(sequence_chars, 2)
        for char_combination in char_combinations:
            if char_combination == ('-', '-'):
                score += self.config["two gaps"]
            elif '-' in char_combination:
                score += self.config["indel"]
            else:
                score += self.substitution_matrix[char_combination[0]][char_combination[1]]

        return self.scoring_matrix.get_score(*comparison_indices) + score

    def fill_scoring_matrix(self):
        """
        Fill the scoring matrix.
        """
        for index in self.scoring_matrix.iter_non_zero_indices():
            self.update_matrix_position(*index)

    @abstractmethod
    def update_matrix_position(self, *args) -> None:
        """
        Update the position in the scoring matrix, setting the score and traceback.
        :param args: Coordinates in the scoring matrix.
        """
        pass

    @abstractmethod
    def initialise_scoring_matrix(self, sequences: List[str]) -> ScoringMatrix:
        """
        Initialise the scoring matrix.
        :param sequences: List of sequences to align.
        :return: Initialised scoring matrix.
        """
        pass

    @abstractmethod
    def reached_stopping_condition(self, *args) -> bool:
        """
        Check if the stopping condition has been reached.
        :param args: Coordinates to check.
        :return: True if the stopping condition has been reached, False otherwise.
        """
        pass

    def traceback(self, *args) -> List[Tuple[str, ...]]:
        """
        Perform the traceback.
        :param args: Coordinates to start the traceback from.
        :return: List of tuples of aligned sequences.
        """
        if self.reached_stopping_condition(*args):
            return [('',) * len(self.scoring_matrix.sequences)]

        new_alignments = []
        for traceback_direction in self.scoring_matrix.get_traceback(*args):
            traceback_alignments = self.traceback(*traceback_direction)
            alignment_chars = self.get_alignment_chars(*args, comparison_indices=traceback_direction)
            for alignment in traceback_alignments:
                new_alignment = []
                for i in range(len(alignment)):
                    new_alignment.append(alignment[i] + alignment_chars[i])
                new_alignments.append(tuple(new_alignment))

        return new_alignments

    def pre_solve(self):
        """
        Pre-solve hook.
        """
        pass

    def post_solve(self, score: Union[int, float], results: List[Tuple[str, ...]]) -> Tuple[
        Union[int, float], List[Tuple[str, ...]]]:
        """
        Post-solve hook.
        :param score: The alignment score.
        :param results: The results of the PSA solver.
        """
        return score, results

    @abstractmethod
    def get_start_indices(self) -> Tuple[int, ...]:
        """
        Get the starting indices for the traceback.
        :return: The starting indices for the traceback.
        """
        pass

    @abstractmethod
    def get_alignment_score(self) -> Union[int, float]:
        """
        Get the alignment score.
        :return: The alignment score.
        """
        pass

    def solve(self, sequences: List[str]) -> Tuple[Union[int, float], List[Tuple[str, ...]]]:
        """
        Solve the multiple sequence alignment problem.
        :param sequences: List of sequences to align.
        :return: Tuple of the aligned sequences and the alignment score.
        """
        self.scoring_matrix = self.initialise_scoring_matrix(sequences)
        self.fill_scoring_matrix()
        self.pre_solve()
        return self.post_solve(self.get_alignment_score(), self.traceback(*self.get_start_indices()))
