from typing import Tuple, Union


def smith_waterman(sequence1: Union[str, Tuple], sequence2: Union[str, Tuple], config: dict,
                   substitution_matrix: dict) -> Tuple[int, Tuple[str, str], Tuple[str, str]]:
    """
    Smith-Waterman algorithm for local sequence alignment.

    :param sequence1: First sequence to be aligned.
        Can be a string representing the sequence, or a tuple with [sequence_id, sequence_string].
    :param sequence2: Second sequence to be aligned.
        Can be a string representing the sequence, or a tuple with [sequence_id, sequence_string].
    :param config: Configuration dictionary.
    :param substitution_matrix: Substitution matrix to be used for scoring.
    :return: Tuple containing the score and the aligned sequences.
    """
    gap_penalty = config['indel']

    def score_function(char1: str, char2: str) -> int:
        """
        Calculate the score for two characters.

        :param char1: First character.
        :param char2: Second character.
        :return: Score for the two characters.
        """
        return substitution_matrix[char1][char2]

    sequence_id1 = "s1"
    sequence_id2 = "s2"

    # If sequences are given as tuples, extract the sequence strings
    if isinstance(sequence1, tuple):
        sequence_id1 = sequence1[0]
        sequence1 = sequence1[1]

    if isinstance(sequence2, tuple):
        sequence_id2 = sequence2[0]
        sequence2 = sequence2[1]

    x_length = len(sequence1)
    y_length = len(sequence2)

    # Initialize scoring matrix
    matrix = [[0 for _ in range(y_length + 1)] for _ in range(x_length + 1)]

    def calc_matrix_score(x: int, y: int) -> int:
        """
        Calculate the score for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :return: Score for the matrix entry.
        """
        if x == 0 or y == 0:
            return 0
        return max(
            matrix[x - 1][y - 1] + score_function(sequence1[x - 1], sequence2[y - 1]),
            max([matrix[x - i][y] - gap_penalty * i for i in range(1, x_length)]),
            max([matrix[x][y - j] - gap_penalty * j for j in range(1, y_length)]),
            0
        )

    # Fill scoring matrix
    for i in range(x_length + 1):
        for j in range(y_length + 1):
            matrix[i][j] = calc_matrix_score(i, j)

    # Find maximum score to start traceback
    max_score = 0
    max_score_index = (0, 0)
    for i in range(x_length + 1):
        for j in range(y_length + 1):
            if matrix[i][j] >= max_score:
                max_score = matrix[i][j]
                max_score_index = (i, j)

    def traceback(x: int, y: int) -> tuple[int, tuple[str, str]]:
        """
        Recursively find the path with the highest score.
        :param x: Row index.
        :param y: Column index.
        :return: Tuple containing the score and the aligned sequences.
        """
        if matrix[x][y] == 0:
            return 0, (sequence1[x], sequence2[y])

        # Follow the path(s) with the highest score
        left_score = matrix[x - 1][y]
        up_score = matrix[x][y - 1]
        diagonal_score = matrix[x - 1][y - 1]
        max_score = max(left_score, up_score, diagonal_score)

        if max_score == diagonal_score:
            score, strings = traceback(x - 1, y - 1)
            diagonal_score, diagonal_string = score + max_score, (
            strings[0] + sequence1[x - 1], strings[1] + sequence2[y - 1])

        if max_score == left_score:
            score, strings = traceback(x - 1, y)
            left_score, left_string = score + max_score, (strings[0] + sequence1[x - 1], strings[1] + "-")

        if max_score == up_score:
            score, strings = traceback(x, y - 1)
            up_score, up_string = score + max_score, (strings[0] + "-", strings[1] + sequence2[y - 1])

        max_score = max(diagonal_score, left_score, up_score)

        if max_score == diagonal_score:
            return diagonal_score, diagonal_string

        if max_score == left_score:
            return left_score, left_string

        return up_score, up_string

    score, strings = traceback(max_score_index[0], max_score_index[1])

    return score, (sequence_id1, strings[0]), (sequence_id2, strings[1])
