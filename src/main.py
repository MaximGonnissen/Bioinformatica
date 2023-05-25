import argparse
from utils import check_and_create_dir


if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(description='Commandline arguments for the program', prog='Bioinformatics')

    # Add arguments
    parser.add_argument('-c', '--config', help='Path to config file', type=str, default='./data/config/config.json')
    parser.add_argument('-i', '--input', help='Path to input file', type=str, default='./data/input/input.fasta')
    parser.add_argument('-o', '--output', help='Path to output file', type=str, default='./data/output/output.txt')

    # Parse arguments
    args = parser.parse_args()

    # Check if output directory exists, create it if not
    check_and_create_dir(args.output)

    # Run program
    # TODO: Run program
