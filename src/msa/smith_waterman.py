from typing import Tuple, Union


def smith_waterman(bottom_sequence: Union[str, Tuple], top_sequence: Union[str, Tuple], config: dict,
                   substitution_matrix: dict, print_matrix: bool = False) -> Tuple[
    int, Tuple[str, str], Tuple[str, str]]:
    """
    Smith-Waterman algorithm for local sequence alignment.

    :param bottom_sequence: First sequence to be aligned.
        Can be a string representing the sequence, or a tuple with [sequence_id, sequence_string].
    :param top_sequence: Second sequence to be aligned.
        Can be a string representing the sequence, or a tuple with [sequence_id, sequence_string].
    :param config: Configuration dictionary.
    :param substitution_matrix: Substitution matrix to be used for scoring.
    :param print_matrix: Whether to print the scoring matrix.
    :return: Tuple containing the score and the aligned sequences.
    """
    gap_penalty = config['indel']

    gap_penalty = int(config['indel'])

    def score_function(char1: str, char2: str) -> Union[int, float]:
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
    if isinstance(bottom_sequence, tuple):
        sequence_id1 = bottom_sequence[0]
        bottom_sequence = bottom_sequence[1]

    if isinstance(top_sequence, tuple):
        sequence_id2 = top_sequence[0]
        top_sequence = top_sequence[1]

    x_length = len(bottom_sequence)
    y_length = len(top_sequence)

    # Initialize scoring matrix
    matrix = [[0 for _ in range(y_length + 1)] for _ in range(x_length + 1)]

    def calc_matrix_score(_x: int, _y: int) -> Union[int, float]:
        """
        Calculate the score for a matrix entry.
        :param _x: Row index.
        :param _y: Column index.
        :return: Score for the matrix entry.
        """
        if _x == 0 or _y == 0:
            return 0
        return max(
            matrix[_x - 1][_y - 1] + score_function(bottom_sequence[_x - 1], top_sequence[_y - 1]),
            max([matrix[_x - _i][_y] - gap_penalty * _i for _i in range(1, _x + 1)]),
            max([matrix[_x][_y - _j] - gap_penalty * _j for _j in range(1, _y + 1)]),
            0
        )

    # Fill scoring matrix
    for i in range(x_length + 1):
        for j in range(y_length + 1):
            matrix[i][j] = calc_matrix_score(i, j)

    if print_matrix:
        top_bar = '\t' * 2 + '\t'.join([f"{char}" for char in top_sequence])
        print(top_bar)
        for i in range(x_length + 1):
            if i == 0:
                print('\t' + '\t'.join([f"{score}" for score in matrix[i]]))
            else:
                print(f"{bottom_sequence[i - 1]}\t" + '\t'.join([f"{int(score)}" for score in matrix[i]]))

    # Find maximum score to start traceback
    max_score = 0
    max_score_index = (0, 0)
    for i in range(x_length + 1):
        for j in range(y_length + 1):
            if matrix[i][j] >= max_score:
                max_score = matrix[i][j]
                max_score_index = (i, j)

    def traceback(_x: int, _y: int) -> tuple[int, tuple[str, str]]:
        """
        Recursively find the path with the highest score.
        :param _x: Row index.
        :param _y: Column index.
        :return: Tuple containing the score and the aligned sequences.
        """
        if matrix[_x][_y] == 0:
            return 0, (bottom_sequence[_x], top_sequence[_y])

        # Follow the path(s) with the highest score
        _up_score = matrix[_x - 1][_y]
        _left_score = matrix[_x][_y - 1]
        _diagonal_score = matrix[_x - 1][_y - 1]
        _max_score = max(_up_score, _left_score, _diagonal_score)

        if _max_score == _diagonal_score:
            _score, _strings = traceback(_x - 1, _y - 1)
            _diagonal_score, diagonal_string = _score + _max_score, (
                _strings[0] + bottom_sequence[_x - 1], _strings[1] + top_sequence[_y - 1])

        if _max_score == _up_score:
            _score, _strings = traceback(_x - 1, _y)
            _up_score, left_string = _score + _max_score, (_strings[0] + bottom_sequence[_x - 1], _strings[1] + "-")

        if _max_score == _left_score:
            _score, _strings = traceback(_x, _y - 1)
            _left_score, up_string = _score + _max_score, (_strings[0] + "-", _strings[1] + top_sequence[_y - 1])

        _max_score = max(_diagonal_score, _up_score, _left_score)

        if _max_score == _diagonal_score:
            return _diagonal_score, diagonal_string

        if _max_score == _up_score:
            return _up_score, left_string

        return _left_score, up_string

    score, strings = traceback(max_score_index[0], max_score_index[1])

    return score, (sequence_id1, strings[0]), (sequence_id2, strings[1])
