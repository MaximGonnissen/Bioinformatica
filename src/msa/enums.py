from enum import Enum


class Direction(Enum):
    """
    Enum for the traceback matrix.
    """
    DIAGONAL = (-1, -1)
    UP = (-1, 0)
    LEFT = (0, -1)
    D = DIAGONAL
    U = UP
    L = LEFT
