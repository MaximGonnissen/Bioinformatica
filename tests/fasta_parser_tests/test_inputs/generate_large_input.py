from src.fasta_parser.fasta_parser import parse


def generate_large_input(source_input: str, length_multiplier: int, sequence_multiplier: int) -> None:
    """
    Generate a large version of the provided file.

    Sequences are repeated `length_multiplier` times, and are copied `sequence_multiplier` times
    (for a total of `sequence_multiplier` sequences of `length_multiplier` * original characters)

    The resulting file can be large, so be careful.
    Using this on the valid_small.fasta file with length_multiplier=100 and sequence_multiplier=1000 results in a ~5MB file.
    Using length_multiplier=1000 and sequence_multiplier=10000 results in a ~500MB (0.5GB) file. For example.

    :param source_input: The path to the source input file
    :param length_multiplier: The number of times to repeat each sequence
    :param sequence_multiplier: The number of times to copy each sequence
    :return: None
    """
    original_result = parse(source_input)

    with open("./valid_large.fasta", "w") as f:
        identifier = 0
        for i in range(sequence_multiplier):
            for seq_id, seq in original_result.items():
                f.write(f">{seq_id}_{identifier}\n")
                f.write(seq * length_multiplier + "\n")
                identifier += 1


if __name__ == '__main__':
    # Generate a large version of the valid_small.fasta file
    # Sequences are repeated 100 times, and are copied 1000 times (for a total of 1000 sequences of 100 * original characters)
    # The resulting file is ~5MB
    generate_large_input("./valid_small.fasta", 100, 1000)
