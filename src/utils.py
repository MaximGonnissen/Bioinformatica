import os


def check_and_create_dir(path: str) -> bool:
    """
    Check if directory exists, create it if not
    :param path: Path to directory
    :return: True if directory was created, False if directory already exists
    """
    path = os.path.dirname(path)
    if not os.path.exists(path):
        os.makedirs(path)
        return True
    return False


def print_lined_up_alignments(alignments: list[tuple[str, ...]]) -> None:
    """
    Print a list of alignments lined up.

    Legend: * to indicate a match, : to indicate a mismatch, space to indicate a gap

    :param alignments: The alignments to print.
    """
    for index, alignment in enumerate(alignments):
        print(f"Alignment {index + 1}:")
        prev_alignment = alignment[0]
        print(f"\t{prev_alignment}")
        for i in range(1, len(alignment)):
            current_alignment = alignment[i]
            alignment_comparison = ""
            for j in range(len(current_alignment)):
                if current_alignment[j] == prev_alignment[j]:
                    alignment_comparison += "*"
                elif current_alignment[j] == "-" or prev_alignment[j] == "-":
                    alignment_comparison += " "
                else:
                    alignment_comparison += ":"
            print(f"\t{alignment_comparison}\n\t{current_alignment}")
            prev_alignment = current_alignment
