import argparse
import json
from typing import List

from fasta_parser.fasta_parser import parse
from msa.needleman_wunsch import NeedlemanWunschMSASolver
from msa.smith_waterman import SmithWatermanMSASolver
from psa.needleman_wunsch import NeedlemanWunschPSASolver
from psa.smith_waterman import SmithWatermanPSASolver
from utils import check_and_create_dir


def main(args: List[str] = None) -> int:
    """
    Main function of the program.
    :param args: Commandline arguments, if None, parse from sys.argv
    :param print_output: Print output to stdout
    :return: Exit code
    """
    # Create argument parser
    parser = argparse.ArgumentParser(description='Commandline arguments for the program', prog='Bioinformatics')

    # Add arguments
    parser.add_argument('-c', '--config', help='Path to config file', type=str, default='./data/config/config.json')
    parser.add_argument('-i', '--input', help='Path to input file', type=str,
                        default='./data/input/cs_assignment.fasta')
    parser.add_argument('-o', '--output', help='Path to output file', type=str, default='./data/output/output.txt')
    parser.add_argument('-v', '--verbose', help='Print output to stdout', action='store_true')

    # Add subparsers for pairwise and multiple sequence alignment
    subparsers = parser.add_subparsers(dest='mode', help='Alignment mode', required=True)
    pairwise_parser = subparsers.add_parser('psa', help='Pairwise sequence alignment')
    msa_parser = subparsers.add_parser('msa', help='Multiple sequence alignment')

    # Add subparsers for pairwise alignment
    pairwise_subparsers = pairwise_parser.add_subparsers(dest='pairwise_mode', help='Pairwise alignment mode',
                                                         required=True)
    needleman_wunsch_parser = pairwise_subparsers.add_parser('needleman_wunsch',
                                                             help='Needleman-Wunsch pairwise alignment')
    smith_waterman_parser = pairwise_subparsers.add_parser('smith_waterman', help='Smith-Waterman pairwise alignment')

    # Add subparsers for multiple sequence alignment
    msa_subparsers = msa_parser.add_subparsers(dest='msa_mode', help='Multiple sequence alignment mode', required=True)
    needleman_wunsch_msa_parser = msa_subparsers.add_parser('needleman_wunsch',
                                                            help='Needleman-Wunsch multiple sequence alignment')
    smith_waterman_msa_parser = msa_subparsers.add_parser('smith_waterman',
                                                          help='Smith-Waterman multiple sequence alignment')

    # Parse arguments
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    # Check if output directory exists, create it if not
    check_and_create_dir(args.output)

    # Read config file
    with open(args.config, 'r') as f:
        config = json.load(f)

    # Run program
    solver = None

    if args.mode in ['pairwise', 'psa']:
        if args.pairwise_mode == 'needleman_wunsch':
            solver = NeedlemanWunschPSASolver(config)
        elif args.pairwise_mode == 'smith_waterman':
            solver = SmithWatermanPSASolver(config)
        else:
            raise ValueError('Invalid pairwise alignment mode')

    elif args.mode == 'msa':
        if args.msa_mode == 'needleman_wunsch':
            solver = NeedlemanWunschMSASolver(config)
        elif args.msa_mode == 'smith_waterman':
            solver = SmithWatermanMSASolver(config)
        else:
            raise ValueError('Invalid multiple sequence alignment mode')

    else:
        raise ValueError('Invalid alignment mode')

    sequence_info = parse(args.input)
    sequence_values = [sequence_info[key] for key in sequence_info.keys()]

    score, alignments = None, []

    if args.mode == 'psa':
        score, alignments = solver.solve(sequence_values[0], sequence_values[1])
    else:
        score, alignments = solver.solve(sequence_values)

    if args.verbose:
        print('Alignments score: {}'.format(score))
        print('Alignments written to {}'.format(args.output))

    if len(alignments) == 0:
        if args.verbose:
            print('No alignments found')
        exit(0)

    alignments = sorted(alignments)

    primary_alignment = alignments[0]
    with open(args.output, 'w') as f:
        for i in range(len(primary_alignment)):
            alignment_id = list(sequence_info.keys())[i]
            f.write(f"{alignment_id}: {primary_alignment[i]}\n")

    if len(alignments) > 1:
        with open(args.output.replace('.txt', '_all.txt'), 'w') as f:
            for alignment in alignments:
                for i in range(len(alignment)):
                    alignment_id = list(sequence_info.keys())[i]
                    f.write(f"{alignment_id}: {alignment[i]}\n")
                f.write('\n')

    return 0


if __name__ == '__main__':
    main()
