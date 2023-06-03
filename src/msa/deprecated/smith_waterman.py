from typing import Callable

from msa.deprecated.msa_solver import MSASolver
from psa.psa_solver import PSASolver
from psa.smith_waterman import SmithWatermanPSASolver


class SmithWatermanMSASolver(MSASolver):
    """
    Class for the Smith-Waterman multiple sequence alignment solver.
    """

    @property
    def solver_cls(self) -> Callable[[dict, dict], PSASolver]:
        """
        Returns the solver class.
        :return: The solver class.
        """
        return SmithWatermanPSASolver
