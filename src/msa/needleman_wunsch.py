from typing import Tuple

import numpy as np

from msa.smith_waterman import SmithWatermanMSASolver


class NeedlemanWunschMSASolver(SmithWatermanMSASolver):
    """
    Class for the Needleman-Wunsch multiple sequence alignment solver.
    """
    add_zero_score = False

    def get_start_indices(self) -> Tuple[int, ...]:
        """
        Get the indices to start the traceback from.
        :return: Indices to start the traceback from.
        """
        return self.scoring_matrix.get_corner_index()

    def get_alignment_score(self) -> int:
        """
        Get the alignment score.
        :return: Alignment score.
        """
        return self.scoring_matrix.get_score(*self.get_start_indices())

    def fill_scoring_matrix(self):
        """
        Fill the scoring matrix.
        """
        for index in self.scoring_matrix.iter_zero_indices():
            distance_from_origin = np.sum(index)
            self.scoring_matrix.set_score(*index, score=self.config["indel"] * distance_from_origin)

        super().fill_scoring_matrix()
