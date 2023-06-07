import itertools
from typing import List, Tuple, Union

from msa.msa_solver import MSASolver
from msa.scoring_matrix.scoring_matrix import ScoringMatrix


class SmithWatermanMSASolver(MSASolver):
    """
    Class for the Smith-Waterman multiple sequence alignment solver.
    """
    add_zero_score = True

    def initialise_scoring_matrix(self, sequences: List[str]) -> ScoringMatrix:
        """
        Initialise the scoring matrix.
        :param sequences: List of sequences to align.
        :return: Initialised scoring matrix.
        """
        return ScoringMatrix(sequences)

    def update_matrix_position(self, *args) -> None:
        """
        Update the position in the scoring matrix, setting the score and traceback.
        :param args: Coordinates in the scoring matrix.
        """
        possible_scores = []
        possible_tracebacks = []

        indices_offset = list(itertools.product((0, -1), repeat=len(args)))
        indices_offset.remove(tuple([0] * len(args)))

        indices_to_check = [tuple([sum(x) for x in zip(args, offset)]) for offset in indices_offset]

        for indices in indices_to_check:
            possible_scores.append(self.scoring_function(*args, comparison_indices=indices))
            possible_tracebacks.append(indices)

        if self.add_zero_score:
            possible_scores.append(0)
            possible_tracebacks.append(None)

        max_score = max(possible_scores)
        self.scoring_matrix.set_score(*args, score=max_score)

        for i in range(len(possible_scores)):
            if possible_scores[i] == max_score:
                self.scoring_matrix.add_traceback(*args, traceback_direction=possible_tracebacks[i])

    def get_start_indices(self) -> Tuple[int, ...]:
        """
        Get the indices to start the traceback from.
        :return: Indices to start the traceback from.
        """
        return self.scoring_matrix.get_max_index()

    def get_alignment_score(self) -> Union[int, float]:
        """
        Get the alignment score.
        :return: Alignment score.
        """
        return self.scoring_matrix.get_max_score()

    def reached_stopping_condition(self, *args) -> bool:
        """
        Check if the stopping condition has been reached.
        :param args: Coordinates to check.
        :return: Whether the stopping condition has been reached.
        """
        if len(self.scoring_matrix.get_traceback(*args)) == 0:
            return True

        return self.scoring_matrix.get_score(*args) == 0


