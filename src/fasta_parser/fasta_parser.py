from typing import Generator, Dict, Tuple


def parse_generator(fasta_file: str) -> Generator[Tuple[str, str], None, None]:
    """
    Parse a fasta file and return a generator that yields tuples of sequence id and sequence.
    :param fasta_file: path to fasta file
    :return: generator that yields tuples of sequence id and sequence
    """
    seq_id = None
    with open(fasta_file, "r") as file:
        for line in file:  # iterate over lines in file, more memory efficient than readlines()
            if line in ["\n", ""]:  # skip empty lines
                continue
            if line.startswith(">"):  # if line starts with >, it is a sequence id
                seq_id = line[1:].strip()  # remove > and trailing whitespaces
                continue
            if seq_id is None:  # if no sequence id was found yet, raise error
                raise ValueError("No sequence id found before sequence was read")
            yield seq_id, line.strip()  # yield sequence id and sequence
            seq_id = None  # reset sequence id


def parse(fasta_file: str) -> Dict[str, str]:
    """
    Parse a fasta file and return a dictionary with the sequence id as key and the sequence as value.
    :param fasta_file: path to fasta file
    :return: dictionary with sequence id as key and sequence as value
    """
    sequences = {}

    for seq_id, line in parse_generator(fasta_file):  # iterate over tuples of sequence id and sequence
        if seq_id in sequences:  # if sequence id was already found, raise error
            raise ValueError("Sequence id found twice: " + seq_id)
        sequences[seq_id] = line  # add sequence id and sequence to dictionary

    return sequences


def parse_str_generator(fasta_str: str) -> Generator[Tuple[str, str], None, None]:
    """
    Parse a fasta string and return a generator that yields tuples of sequence id and sequence.
    :param fasta_str: fasta string
    :return: generator that yields tuples of sequence id and sequence
    """
    seq_id = None
    for line in fasta_str.split("\n"):  # iterate over lines in string
        if line == "":  # skip empty lines
            continue
        if line.startswith(">"):  # if line starts with >, it is a sequence id
            seq_id = line[1:].strip()  # remove > and trailing whitespaces
            continue
        if seq_id is None:  # if no sequence id was found yet, raise error
            raise ValueError("No sequence id found before sequence was read")
        yield seq_id, line.strip()  # yield sequence id and sequence
        seq_id = None  # reset sequence id


def parse_str(fasta_str: str) -> Dict[str, str]:
    """
    Parse a fasta string and return a dictionary with the sequence id as key and the sequence as value.
    :param fasta_str: fasta string
    :return: dictionary with sequence id as key and sequence as value
    """
    sequences = {}

    for seq_id, line in parse_str_generator(fasta_str):  # iterate over tuples of sequence id and sequence
        if seq_id in sequences:  # if sequence id was already found, raise error
            raise ValueError("Sequence id found twice: " + seq_id)
        sequences[seq_id] = line  # add sequence id and sequence to dictionary

    return sequences
