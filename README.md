# Bioinformatics project

> Author: Maxim Gonnissen (S0175700)
>
> Year: 2022-2023
>
> Course: Bioinformatics
> 
> University: University of Antwerp
> 
> Faculty: Faculty of Science

## Description

This project determines the optimal [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
local alignment, the optimal [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) global
alignment, and their alignment scores for an arbitrary number of DNA/protein sequences.

## Configuration

Score is calculated based on parameters provided in a config file.

### Example config file

```json
{
  "match": 5,
  "mismatch": -2,
  "indel": -4,
  "two gaps": 0
}
```

## Input files

Input files are expected to be in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.

## Output files

Output files are in a custom format, with the following structure:

```
s<sequence ID from FASTA input file>: <sequence>
```

### Example output file

```
s1: ACTG.GT.CA
s2: .CAGGGT.CA
s3: CCAGGGACCA
```

## Requirements

- [Python 3.11](https://www.python.org/downloads/release/python-3110/)
- Requirements listed in [requirements.txt](requirements.txt)
  - blosum 2.0.2 -- [documentation](https://pypi.org/project/blosum/) `Used for fetching the BLOSUM62 matrix`

## Usage

Run from the root of the project (above the `src` folder):

```shell
python main.py -c <config file> -i <input file> -o <output file>
```

Where `<config file>` is the path to the config file, `<input file>` is the path to the input file, and `<output file>`
is
the path to the output file - all relative to the root of the project.

### Example usage

```shell
python main.py -c ./data/config.json -i ./data/input/input.fasta -o ./data/output/output.txt
```

## Testing

Tests are provided in the `tests` folder. They can be examined as a reference for the expected output of the program,
as well as examples of how to use the program.