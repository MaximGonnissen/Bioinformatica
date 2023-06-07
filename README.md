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
  - numpy 1.24.3 -- [documentation](https://numpy.org/doc/stable/) `Used for matrix calculations`

## Usage

Run from the root of the project (above the `src` folder):

```shell
python main.py -c <config file> -i <input file> -o <output file>
```

Where `<config file>` is the path to the config file, `<input file>` is the path to the input file, and `<output file>`
is the path to the output file - all relative to the root of the project.

### Example usage

```shell
python ./src/main.py -c ./data/config/config.json -i ./data/input/cs_assignment.fasta -o ./data/output/output.txt msa needleman_wunsch
```

## Testing

Tests are provided in the `tests` folder. They can be examined as a reference for the expected output of the program,
as well as examples of how to use the source as a library.

## Packages

### [src](src)

This is the top level package for the project. It contains the [main](src/main.py) module, which is a CLI interface
for the project. It also contains the [utils](src/utils.py) module, which contains two simple utility methods.
One to create a directory if it does not exist, and one to print alignments in a pretty format.

### [fasta_parser](src/fasta_parser)

This package contains the [fasta_parser](src/fasta_parser/fasta_parser.py) module, which is used to parse FASTA files.
This module allows parsing from both a file and a string, and also provides iterator variants of these methods in case
the FASTA file is too large to fit in memory.

### [msa](src/msa)

This package contains everything related to multiple sequence alignment.

An abstract base class for msa solvers is provided in the [msa_solver](src/msa/msa_solver.py) module. This class
contains a number of base methods which are reused in the concrete implementations of the msa solvers, as well as a
number of abstract methods which must be implemented by these implementations. The main method of this class is the
`solve` method, which chains the different abstract methods together. It also provides a pre- and post-solve abstract
method, which can be used to perform actions before and after the actual solving of the msa. These are not used in the
concrete implementations, but are provided in case they are useful for any future implementations.

The actual implementations of the msa solvers are provided in the [needleman_wunsch](src/msa/needleman_wunsch.py) and
[smith_waterman](src/msa/smith_waterman.py) modules, which contain the NeedlemanWunschMSASolver and
SmithWatermanMSASolver classes respectively. The NeedlemanWunschMSASolver class inherits from the
SmithWatermanMSASolver class, which inherits from the MSASolver class. This is because there is a large amount of
code reuse between the two implementations.

#### [scoring_matrix](src/msa/scoring_matrix)

This subpackage of `msa` contains the [scoring_matrix](src/msa/scoring_matrix/scoring_matrix.py) module. This module
contains the ScoringMatrixEntry and ScoringMatrix classes, which (together) are used to represent a scoring matrix.

The ScoringMatrixEntry class represents a single entry in a scoring matrix. It contains a score and traceback.
A number of methods are implemented to work nicely with the ScoringMatrix class, as well as feel like a general
python type. For example, it implements methods such as `__eq__`, `__lt__`, `__gt__`, `__add__`, `__sub__`, `__mul__`.

The ScoringMatrix class represents a scoring matrix. It contains an n-dimensional numpy array of ScoringMatrixEntry.
Dimensions are chosen based on the number of sequences to align. A plethora of methods are implemented to access and
manipulate the scoring matrix. For example, it implements methods such as `__getitem__`, `__setitem__`, `__iter__`.
In addition, it also implements methods to quickly grab msa relevant information from the scoring matrix, such as
the index with the highest score, whether an index is part of a zero plane, iterating over zero plane indices, etc...

### [psa](src/psa)

This package contains everything related to pairwise sequence alignment.

Technically, this module is not necessary since the msa solvers can also be used to perform pairwise sequence.
However, during development, I found it useful to first implement pairwise sequence alignment as an exercise to
better understand the algorithms, before implementing the msa solvers.

Similarly to the MSA package, an abstract base class for psa solvers is provided in the
[psa_solver](src/psa/psa_solver.py) module. This class likewise contains a number of base and abstract methods, which
are reused in the concrete implementations of the psa solvers. The main method of this class is once again the `solve`
method, and it works analogously to the `solve` method of the MSASolver class.

The actual implementations of the psa solvers are provided in the [needleman_wunsch](src/psa/needleman_wunsch.py) and
[smith_waterman](src/psa/smith_waterman.py) modules, which contain the NeedlemanWunschPSASolver and
SmithWatermanPSASolver classes respectively. The NeedlemanWunschPSASolver class inherits from the
SmithWatermanPSASolver class, which inherits from the PSASolver class.