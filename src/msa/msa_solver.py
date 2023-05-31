from abc import ABC, abstractmethod
from typing import List, Tuple, Callable, Optional

from blosum import BLOSUM

from psa.psa_solver import PSASolver
from scoring_matrix.scoring_matrix import ScoringMatrix


class MSASolver(ABC):
    """
    Abstract class for multiple sequence alignment solvers.
    """

    def __init__(self, config: dict, substitution_matrix: Optional[dict] = None, *args, **kwargs):
        """
        Initialise the PSA solver.
        :param config: Configuration for the PSA solver.
        :param substitution_matrix: Substitution matrix to use for the PSA solver. Defaults to BLOSUM62.
        """
        super().__init__(*args, **kwargs)
        self.config = config
        self.substitution_matrix = substitution_matrix or BLOSUM(62)
        self.scoring_matrix: ScoringMatrix = None

        if None not in [config.get("match"), config.get("mismatch")]:
            for key in self.substitution_matrix.keys():
                for other_key in self.substitution_matrix[key].keys():
                    if key == other_key:
                        self.substitution_matrix[key][other_key] = config["match"]
                    else:
                        self.substitution_matrix[key][other_key] = config["mismatch"]

    @property
    @abstractmethod
    def solver_cls(self) -> Callable[[dict, dict], PSASolver]:
        """
        Returns the solver class.
        :return: The solver class.
        """
        pass

    def solve(self, sequences: List[str]) -> List[Tuple[str, str]]:
        """
        Solves the multiple sequence alignment problem for the given sequences.
        :param sequences: The sequences to align.
        :return: The aligned sequences.
        """
        sequence_combinations = []
        for sequence in sequences:
            for other_sequence in sequences:
                if sequence != other_sequence:
                    sequence_combinations.append((sequence, other_sequence))

        solvers = []
        for combination in sequence_combinations:
            combination = sorted(combination,
                                 key=lambda x: len(x))  # Sort by length to help consistency in matrix dimensions
            solver = self.solver_cls(self.config, self.substitution_matrix)
            solver.solve(sequence_1=combination[0], sequence_2=combination[1])
            solvers.append(solver)

        new_solver_matrices = {}
        for target_solver in solvers:
            # First, we create a new matrix of the same type and dimensions as the target solver's matrix
            new_solver_matrix = target_solver.scoring_matrix_cls(target_solver.sequence_1, target_solver.sequence_2)

            # Then, we set all scores to 0 to prepare for the summation
            for i in range(new_solver_matrix.width()):
                for j in range(new_solver_matrix.height()):
                    new_solver_matrix.set_score(i, j, 0)

            for solver in solvers:
                solver_matrix = solver.scoring_matrix
                for i in range(solver_matrix.width()):
                    for j in range(solver_matrix.height()):
                        if i <= new_solver_matrix.width() and j <= new_solver_matrix.height():
                            new_solver_matrix.set_score(i, j,
                                                        new_solver_matrix.get_score(i, j) + solver_matrix.get_score(i,
                                                                                                                    j))

            new_solver_matrices[target_solver] = new_solver_matrix

        for solver in solvers:
            solver.scoring_matrix = new_solver_matrices[solver]

        results = []
        for solver in solvers:
            results += solver.solve(sequence_1=solver.sequence_1, sequence_2=solver.sequence_2,
                                    scoring_matrix=solver.scoring_matrix)[1]

        return results
