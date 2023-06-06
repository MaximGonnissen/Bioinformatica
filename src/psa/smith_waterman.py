from typing import Tuple, Union, Optional, List, Callable

from psa.enums import Direction
from psa.psa_solver import PSASolver
from psa.scoring_matrix.scoring_matrix import ScoringMatrix


class SmithWatermanPSASolver(PSASolver):
    """
    Smith-Waterman solver for the PSA problem.
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
        return ScoringMatrix.smith_waterman

    def get_starting_points(self) -> List[Tuple[int, int]]:
        """
        Get the starting points for the PSA problem.
        :return: List of starting points.
        """
        return self.scoring_matrix.max_score_index_multiple()

    def calculate_scoring_matrix(self) -> None:
        """
        Calculate the scoring matrix.
        """
        for i in range(self.scoring_matrix.height()):
            for j in range(self.scoring_matrix.width()):
                scores = self.calc_matrix_score(i, j)
                max_score = max(scores)
                self.scoring_matrix.set_score(i, j, max_score)
                if len(scores) > 1:
                    if max_score == scores[0]:
                        self.scoring_matrix.add_traceback(i, j, Direction.DIAGONAL)
                    if max_score == scores[1]:
                        self.scoring_matrix.add_traceback(i, j, Direction.UP)
                    if max_score == scores[2]:
                        self.scoring_matrix.add_traceback(i, j, Direction.LEFT)

    def _solve(self) -> Tuple[int, List[Tuple[str, str]]]:
        """
        Solve the PSA problem.
        :return: Tuple containing the score and a list of valid alignments.
        """
        # Calculate the scoring matrix
        self.calculate_scoring_matrix()

        # Find the starting points
        potential_starting_points = self.get_starting_points()

        score = max([self.scoring_matrix.get_score(x, y) for x, y in potential_starting_points])

        # Run traceback for each potential starting point
        results = []
        for starting_point in potential_starting_points:
            results += self.traceback(starting_point[0], starting_point[1])

        return score, results

    def calc_matrix_score(self, x: int, y: int) -> list[Union[int, float]]:
        """
        Calculate the scores for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :return: Scores for the matrix entry.
        """
        if x == 0 or y == 0:
            return [0]
        return [
            self.scoring_matrix.get_score(x - 1, y - 1) + self.scoring_function(self.scoring_matrix.bottom_char(x),
                                                                                self.scoring_matrix.top_char(y)),
            max([self.scoring_matrix.get_score(x - i, y) - self.gap_penalty(i) for i in range(1, x + 1)]),
            max([self.scoring_matrix.get_score(x, y - j) - self.gap_penalty(j) for j in range(1, y + 1)]),
            0
        ]

    def traceback(self, x: int, y: int) -> list[tuple[str, str]]:
        """
        Recursively find all valid paths.
        :param x: Row index.
        :param y: Column index.
        :return: A list of tuples with alignments.
        """
        if x == 0 or y == 0:
            return [('', '')]

        if len(self.scoring_matrix.get_traceback(x, y)) == 0:
            return [('', '')]

        # Follow the traceback path
        new_paths = []
        for direction in self.scoring_matrix.get_traceback(x, y):
            if direction == Direction.DIAGONAL:
                paths = self.traceback(x - 1, y - 1)
                for alignment in paths:
                    new_paths.append((alignment[0] + self.scoring_matrix.top_char(y),
                                      alignment[1] + self.scoring_matrix.bottom_char(x)))

            elif direction == Direction.UP:
                paths = self.traceback(x - 1, y)
                for alignment in paths:
                    new_paths.append((alignment[0] + '-', alignment[1] + self.scoring_matrix.bottom_char(x)))

            elif direction == Direction.LEFT:
                paths = self.traceback(x, y - 1)
                for alignment in paths:
                    new_paths.append((alignment[0] + self.scoring_matrix.top_char(y), alignment[1] + '-'))

        # TODO: This is broken. Undo & just invert width & height in scoring matrix?

        return list(set(new_paths))
