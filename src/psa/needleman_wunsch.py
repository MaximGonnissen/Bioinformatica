from typing import Union, Optional, Callable

from psa.scoring_matrix.scoring_matrix import ScoringMatrix
from psa.smith_waterman import SmithWatermanPSASolver


class NeedlemanWunschPSASolver(SmithWatermanPSASolver):
    """
    Needleman-Wunsch solver for the PSA problem.
    """

    def __init__(self, config: dict, substitution_matrix: Optional[dict] = None, *args, **kwargs):
        """
        Initialise the solver.
        :param config: Configuration dictionary.
        :param substitution_matrix: Substitution matrix.
        """
        super().__init__(config, substitution_matrix, *args, **kwargs)

    @property
    def scoring_matrix_cls(self) -> Callable:
        """
        Return the scoring matrix class to use.
        :return: Scoring matrix class.
        """
        return ScoringMatrix.needleman_wunsch

    def get_starting_points(self) -> list[tuple[int, int]]:
        """
        Get the starting points for the PSA problem.
        :return: List of starting points.
        """
        return [self.scoring_matrix.corner_index()]

    def calc_matrix_score(self, x: int, y: int) -> list[Union[int, float]]:
        """
        Calculate the scores for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :return: Scores for the matrix entry.
        """
        if x == 0 and y == 0:
            return [0]
        elif x == 0:
            return [self.gap_penalty(y)]
        elif y == 0:
            return [self.gap_penalty(x)]
        return [
            self.scoring_matrix.get_score(x - 1, y - 1) + self.scoring_function(self.scoring_matrix.bottom_char(x),
                                                                                self.scoring_matrix.top_char(y)),
            max([self.scoring_matrix.get_score(x - i, y) + self.gap_penalty(i) for i in range(1, x + 1)]),
            max([self.scoring_matrix.get_score(x, y - j) + self.gap_penalty(j) for j in range(1, y + 1)])
        ]
