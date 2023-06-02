from typing import Union, Tuple, List, Any, Iterable

import numpy as np


class ScoringMatrix:
    """
    N-dimensional scoring matrix for MSA algorithms.
    """

    def __init__(self, sequences: List[str], *args, **kwargs):
        """
        Initialise the scoring matrix.
        :param sequences: List of sequences.
        """
        self.sequences = sequences
        self.shape = tuple(len(sequence) + 1 for sequence in sequences)
        self.matrix = np.zeros(self.shape, dtype=object)

    def __getitem__(self, item: Union[int, Tuple[int, ...]]) -> Any:
        """
        Get an item from the matrix.
        :param item: Index of the item.
        :return: The item.
        """
        return self.matrix[item]

    def __setitem__(self, key: Union[int, Tuple[int, ...]], value) -> None:
        """
        Set an item in the matrix.
        :param key: Index of the item.
        :param value: Value to set the item to.
        """
        self.matrix[key] = value

    def __iter__(self) -> Iterable:
        """
        Iterate over the matrix.
        :return: Iterator for the matrix.
        """
        return iter(self.matrix)

    def __str__(self) -> str:
        """
        Convert the matrix to a string.
        :return: String representation of the matrix.
        """
        return str(self.matrix)

    def __repr__(self) -> str:
        """
        Convert the matrix to a string.
        :return: String representation of the matrix.
        """
        return str(self.matrix)
