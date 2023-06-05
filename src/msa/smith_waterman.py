from typing import List, Tuple

from msa.msa_solver import MSASolver
from msa.scoring_matrix.scoring_matrix import ScoringMatrix


class SmithWatermanMSASolver(MSASolver):
    """
    Class for the Smith-Waterman multiple sequence alignment solver.
    """

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
        score = self.scoring_matrix.get_score(*args)

    def traceback(self, *args) -> List[Tuple[str, ...]]:
        """
        Perform the traceback.
        :param args: Coordinates to start the traceback from.
        :return: List of tuples of aligned sequences.
        """
        if self.scoring_matrix.is_zero_index(*args):
            return []

        if len(self.scoring_matrix.get_traceback(*args)) == 0:
            return []  # TODO: May need to return sequence chars instead of nothing, in this case?

        new_alignments = []
        for traceback_direction in self.scoring_matrix.get_traceback(*args):
            traceback_alignments = self.traceback(*traceback_direction)
            for alignment in traceback_alignments:
                new_alignment = []
                for i in range(len(alignment)):
                    new_alignment.append(alignment[i] + self.scoring_matrix.sequences[i][args[i] - 1])
                    # TODO: Compare index to previous index, for each index that's the same, add a gap. ~ Something like this.
                new_alignments.append(tuple(new_alignment))

        return new_alignments
