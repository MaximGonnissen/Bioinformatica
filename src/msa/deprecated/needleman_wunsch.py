from typing import Callable

from msa.deprecated.msa_solver import MSASolver
from psa.needleman_wunsch import NeedlemanWunschPSASolver
from psa.psa_solver import PSASolver


class NeedlemanWunschMSASolver(MSASolver):
    """
    Class for the Needleman-Wunsch multiple sequence alignment solver.
    """

    @property
    def solver_cls(self) -> Callable[[dict, dict], PSASolver]:
        """
        Returns the solver class.
        :return: The solver class.
        """
        return NeedlemanWunschPSASolver
