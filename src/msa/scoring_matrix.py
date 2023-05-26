from typing import Union, Tuple

from msa.enums import Direction


class ScoringMatrix(list):
    """
    Class for the scoring matrix.
    """

    def __init__(self, top_sequence: str, bottom_sequence: str, *args, **kwargs):
        """
        Initialise the scoring matrix.
        :param top_sequence: Top sequence string.
        :param bottom_sequence: Bottom sequence string.
        """
        super().__init__(*args, **kwargs)
        self.top_sequence = top_sequence
        self.bottom_sequence = bottom_sequence
        self.width = len(bottom_sequence) + 1
        self.height = len(top_sequence) + 1

        self.initialise()

    def initialise(self):
        """
        Initialise the scoring matrix.
        """
        self.clear()

        for i in range(self.width):
            self.append([])
            for j in range(self.height):
                self[i].append([0, []])

    def get_score(self, x: int, y: int) -> Union[int, float]:
        """
        Get the score for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :return: Score for the matrix entry.
        """
        return self[x][y][0]

    def get_traceback(self, x: int, y: int) -> list[Direction]:
        """
        Get the traceback for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :return: Traceback for the matrix entry.
        """
        return self[x][y][1]

    def set_score(self, x: int, y: int, score: Union[int, float]):
        """
        Set the score for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :param score: Score for the matrix entry.
        """
        self[x][y][0] = score

    def set_traceback(self, x: int, y: int, traceback: list[Direction]):
        """
        Set the traceback for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :param traceback: Traceback for the matrix entry.
        """
        self[x][y][1] = traceback

    def add_traceback(self, x: int, y: int, traceback: Direction):
        """
        Add a traceback to the traceback list for a matrix entry.
        :param x: Row index.
        :param y: Column index.
        :param traceback: Traceback to be added.
        """
        self[x][y][1].append(traceback)

    def max_score_index(self) -> Tuple[int, int]:
        """
        Find the highest scoring matrix entry.
        :return: Indices of the highest scoring matrix entry.
        """
        max_score = 0
        max_score_index = (0, 0)
        for i in range(self.width):
            for j in range(self.height):
                if self.get_score(i, j) > max_score:
                    max_score = self.get_score(i, j)
                    max_score_index = (i, j)
        return max_score_index

    def max_score(self) -> Union[int, float]:
        """
        Find the highest scoring matrix entry.
        :return: Score of the highest scoring matrix entry.
        """
        return self.get_score(*self.max_score_index())

    def max_score_index_multiple(self) -> list[Tuple[int, int]]:
        """
        Find the indices of the highest scoring matrix entries.
        :return: Indices of the highest scoring matrix entries.
        """
        max_score = self.max_score()
        max_score_indices = []
        for i in range(self.width):
            for j in range(self.height):
                if self.get_score(i, j) == max_score:
                    max_score_indices.append((i, j))
        return max_score_indices

    def print_matrix(self) -> None:
        """
        Print the scoring matrix.
        """
        top_bar = '\t' * 2 + '\t'.join([f"{char}" for char in self.top_sequence])
        print(top_bar)
        for i in range(self.width):
            row = '\t' * 2 + '\t'.join([f"{self.get_score(i, j)}" for j in range(self.height)])
            print(row)
