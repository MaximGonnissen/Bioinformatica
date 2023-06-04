from typing import List, Tuple, Union

import numpy as np

from msa.msa_solver import MSASolver
from msa.scoring_matrix.scoring_matrix import ScoringMatrix


class SmithWatermanMSASolver(MSASolver):
    """
    Class for the Smith-Waterman multiple sequence alignment solver.
    """

    def initialise_scoring_matrix(self, sequences: List[str]) -> ScoringMatrix:
        """
        Initialise the scoring matrix.
        :param sequences: List of sequences to align.
        :return: Initialised scoring matrix.
        """
        return ScoringMatrix(sequences)

    def scoring_function(self, *args) -> Union[int, float]:
        """
        Scoring function for the MSA solver.
        :param args: Coordinates for the scoring matrix.
        :return: Score for the given coordinates.
        """
        # Simple version for now, only taking into account matches and mismatches.
        indices = args
        if len(indices) != len(self.scoring_matrix.shape):
            raise IndexError("Incorrect number of indices given.")

        # Get the indices for the sequences.
        sequence_chars = [self.scoring_matrix.sequences[i][indices[i] - 1] for i in range(len(indices))]

        score = 0
        used_chars = []
        for char in sequence_chars:
            for other_char in sequence_chars:
                if char != other_char and char not in used_chars:
                    score += self.substitution_matrix[char][other_char]
            used_chars.append(char)

        return score

    def fill_scoring_matrix(self):
        """
        Fill the scoring matrix.
        """
        for index in self.scoring_matrix.iter_non_zero_indices():
            self.scoring_matrix.set_score(*index, score=self.scoring_function(*index))

    def traceback(self) -> Tuple[Union[int, float], List[str]]:
        pass
