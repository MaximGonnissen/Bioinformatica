from abc import ABC, abstractmethod
from typing import List, Tuple, Union


class MSASolver(ABC):
    """
    Abstract class for multiple sequence alignment solvers.
    """

    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def solve(self, sequences: List[str]) -> Tuple[int, List[Tuple[str, str]]]:
        """
        Solves the multiple sequence alignment problem for the given sequences.
        :param sequences: The sequences to align.
        :return: The score of the optimal alignment and the alignments.
        """
        pass

    @abstractmethod
    def scoring_function(self, items: List[str]) -> Union[int, float]:
        """
        Scoring function for the multiple sequence alignment problem.
        :param items: The items to score.
        :return: The score.
        """
        pass
