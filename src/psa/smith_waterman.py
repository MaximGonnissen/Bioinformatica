from typing import Tuple, Union

from psa.scoring_matrix import ScoringMatrix
from src.enums import Direction


def smith_waterman(bottom_sequence: str, top_sequence: str, config: dict, substitution_matrix: dict,
                   print_matrix: bool = False) -> Tuple[int, list[Tuple[str, str]]]:
    """
    Smith-Waterman algorithm for local sequence alignment.

    :param bottom_sequence: First sequence to be aligned.
        Can be a string representing the sequence, or a tuple with [sequence_id, sequence_string].
    :param top_sequence: Second sequence to be aligned.
        Can be a string representing the sequence, or a tuple with [sequence_id, sequence_string].
    :param config: Configuration dictionary.
    :param substitution_matrix: Substitution matrix to be used for scoring.
    :param print_matrix: Whether to print the scoring matrix in the console.
    :return: Tuple containing the score and a list of tuples with the aligned sequences matching the score.
    """

    gap_penalty = int(config['indel'])

    def score_function(char1: str, char2: str) -> Union[int, float]:
        """
        Calculate the score for two characters.

        :param char1: First character.
        :param char2: Second character.
        :return: Score for the two characters.
        """
        return substitution_matrix[char1][char2]

    x_length = len(bottom_sequence)
    y_length = len(top_sequence)

    # Initialise scoring matrix
    matrix = ScoringMatrix(top_sequence, bottom_sequence)

    def calc_matrix_scores(_x: int, _y: int) -> list[Union[int, float]]:
        """
        Calculate the scores for a matrix entry.
        :param _x: Row index.
        :param _y: Column index.
        :return: Scores for the matrix entry.
        """
        if _x == 0 or _y == 0:
            return [0]
        return [
            matrix.get_score(_x - 1, _y - 1) + score_function(bottom_sequence[_x - 1], top_sequence[_y - 1]),
            max([matrix.get_score(_x - _i, _y) - gap_penalty * _i for _i in range(1, _x + 1)]),
            max([matrix.get_score(_x, _y - _j) - gap_penalty * _j for _j in range(1, _y + 1)]),
            0
        ]

    # Calculate scores and traceback directions for each matrix entry
    for i in range(x_length + 1):
        for j in range(y_length + 1):
            scores = calc_matrix_scores(i, j)
            matrix.set_score(i, j, max(scores))
            if len(scores) > 1:
                if matrix.get_score(i, j) == scores[0]:
                    matrix.add_traceback(i, j, Direction.DIAGONAL)
                if matrix.get_score(i, j) == scores[1]:
                    matrix.add_traceback(i, j, Direction.UP)
                if matrix.get_score(i, j) == scores[2]:
                    matrix.add_traceback(i, j, Direction.LEFT)

    if print_matrix:
        matrix.print_matrix()

    def traceback(_x: int, _y: int) -> tuple[int, list[tuple[str, str]]]:
        """
        Recursively find the paths with the highest score.
        :param _x: Row index.
        :param _y: Column index.
        :return: Tuple containing the score and a list of tuples with alignments matching the score.
        """
        if _x == 0 or _y == 0:
            return matrix.get_score(_x, _y), [('', '')]

        if matrix.get_score(_x, _y) == 0 or len(matrix.get_traceback(_x, _y)) == 0:
            return matrix.get_score(_x, _y), [('', '')]

        # Follow the traceback path
        potential_paths_scored = []
        for direction in matrix.get_traceback(_x, _y):
            if direction == Direction.DIAGONAL:
                score, paths = traceback(_x - 1, _y - 1)
                new_score = score + matrix.get_score(_x, _y)
                new_alignments = []
                for alignment in paths:
                    new_alignments.append((alignment[0] + bottom_sequence[_x - 1], alignment[1] + top_sequence[_y - 1]))
                potential_paths_scored.append((new_score, new_alignments))

            elif direction == Direction.UP:
                score, paths = traceback(_x - 1, _y)
                new_score = score + matrix.get_score(_x, _y)
                new_alignments = []
                for alignment in paths:
                    new_alignments.append((alignment[0] + bottom_sequence[_x - 1], alignment[1] + '-'))
                potential_paths_scored.append((new_score, new_alignments))

            elif direction == Direction.LEFT:
                score, paths = traceback(_x, _y - 1)
                new_score = score + matrix.get_score(_x, _y)
                new_alignments = []
                for alignment in paths:
                    new_alignments.append((alignment[0] + '-', alignment[1] + top_sequence[_y - 1]))
                potential_paths_scored.append((new_score, new_alignments))

            # Return the paths with the highest score
            max_score = max([path[0] for path in potential_paths_scored])  # Find the max score
            max_paths = [path[1] for path in potential_paths_scored if
                         path[0] == max_score]  # Find the paths with the max score

            # Merge all max score results
            actual_paths = []
            for path in max_paths:
                actual_paths += path

            # Remove duplicates
            actual_paths = list(set(actual_paths))

            return max_score, actual_paths

    # Find the starting points
    potential_starting_points = matrix.max_score_index_multiple()

    # Run traceback for each potential starting point
    max_score = 0
    results = []
    for starting_point in potential_starting_points:
        results.append(traceback(starting_point[0], starting_point[1]))

    # Merge all max score results
    max_score = max([result[0] for result in results])
    max_results = [result[1] for result in results if result[0] == max_score]

    # Merge all max score results
    actual_results = []
    for result in max_results:
        actual_results += result

    # Remove duplicates
    actual_results = list(set(actual_results))

    return max_score, actual_results
