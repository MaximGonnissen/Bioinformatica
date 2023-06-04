from typing import List, Tuple, Union

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

    def scoring_function(self, *args) -> Union[int, float]:
        """
        Scoring function for the MSA solver.
        :param args: Coordinates for the scoring matrix.
        :return: Score for the given coordinates.
        """
        # Simple version for now, only taking into account matches and mismatches.
        indices = args
        if len(indices) != len(self.scoring_matrix.shape):
            raise IndexError("Incorrect number of indices given.")

        # Get the indices for the sequences.
        sequence_chars = [self.scoring_matrix.sequences[i][indices[i] - 1] for i in range(len(indices))]

        score = 0
        used_chars = []
        for char in sequence_chars:
            for other_char in sequence_chars:
                if char != other_char and char not in used_chars:
                    score += self.substitution_matrix[char][other_char]
            used_chars.append(char)

        return score

    def fill_scoring_matrix(self):
        """
        Fill the scoring matrix.
        """
        for index in self.scoring_matrix.iter_non_zero_indices():
            self.scoring_matrix.set_score(*index, score=self.scoring_function(*index))

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
