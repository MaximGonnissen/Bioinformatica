import os
import unittest
from typing import List

from main import main


class TestCLI(unittest.TestCase):
    """
    Test the CLI by running it through os.system().
    """
    input_dir_path = "./test_inputs/"
    output_dir_path = "./test_outputs/"
    config_dir_path = "./test_configs/"
    main_path = "../src/main.py"

    # Activate the virtual environment
    base_command_string = "run ../venv/Scripts/activate.bat && "

    def setUp(self) -> None:
        self.input_file_path = os.path.join(self.input_dir_path, 'input.fasta')
        self.output_file_path = os.path.join(self.output_dir_path, 'output.txt')
        self.config_file_path = os.path.join(self.config_dir_path, 'config.json')

    def get_command_string(self) -> str:
        return self.base_command_string + f"python {self.main_path} -c {self.config_file_path} -i {self.input_file_path} -o {self.output_file_path} "

    def verify_output(self, correct_output: List[str]) -> None:
        """
        Verify the output of the program is correct.
        :param correct_output: List
        """
        with open(self.output_file_path, 'r') as output_file:
            for line in output_file:
                self.assertEqual(line.strip(), correct_output.pop(0).strip())

    def test_cli_psa(self):
        """
        Test the CLI with pairwise sequence alignment.
        """
        command_string = self.get_command_string() + "psa smith_waterman"

        os.system(command_string)

        correct_output_lines = [
            "unknown_J_region_1: FGSGTRL",
            "unknown_J_region_2: FGQGTRL"
        ]

        self.verify_output(correct_output_lines)

    def test_cli_psa2(self):
        """
        Test the CLI with pairwise sequence alignment.
        """
        command_string = self.get_command_string() + "psa needleman_wunsch"

        os.system(command_string)

        correct_output_lines = [
            "unknown_J_region_1: GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2: N-TEA---FFGQGTRL-TVV"
        ]

        self.verify_output(correct_output_lines)

    def test_cli_msa(self):
        """
        Test the CLI with multiple sequence alignment.
        """
        command_string = self.get_command_string() + "msa smith_waterman"

        os.system(command_string)

        correct_output_lines = [
            "unknown_J_region_1: FGSGTRLSIR",
            "unknown_J_region_2: FGQGTRLTVV",
            "unknown_J_region_3: FGSGTRLTVV"
        ]

        self.verify_output(correct_output_lines)

    def test_cli_msa2(self):
        """
        Test the CLI with multiple sequence alignment.
        """
        command_string = self.get_command_string() + "msa needleman_wunsch"

        os.system(command_string)

        correct_output_lines = [
            "unknown_J_region_1: GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2: ----NTEAFFGQGTRL-TVV",
            "unknown_J_region_3: ----NYGYTFGSGTRL-TVV"
        ]

        self.verify_output(correct_output_lines)


class TestCLIDirect(TestCLI):
    """
    Test the CLI by importing the main module.
    """

    def test_cli_psa(self):
        """
        Test the CLI with pairwise sequence alignment.
        """
        args = [
            "-c", self.config_file_path,
            "-i", self.input_file_path,
            "-o", self.output_file_path,
            "psa", "smith_waterman"
        ]

        main(args=args)

        correct_output_lines = [
            "unknown_J_region_1: FGSGTRL",
            "unknown_J_region_2: FGQGTRL"
        ]

        self.verify_output(correct_output_lines)

    def test_cli_psa2(self):
        """
        Test the CLI with pairwise sequence alignment.
        """
        args = [
            "-c", self.config_file_path,
            "-i", self.input_file_path,
            "-o", self.output_file_path,
            "psa", "needleman_wunsch"
        ]

        main(args=args)

        correct_output_lines = [
            "unknown_J_region_1: GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2: N-TEA---FFGQGTRL-TVV"
        ]

        self.verify_output(correct_output_lines)

    def test_cli_msa(self):
        """
        Test the CLI with multiple sequence alignment.
        """
        args = [
            "-c", self.config_file_path,
            "-i", self.input_file_path,
            "-o", self.output_file_path,
            "msa", "smith_waterman"
        ]

        main(args=args)

        correct_output_lines = [
            "unknown_J_region_1: FGSGTRLSIR",
            "unknown_J_region_2: FGQGTRLTVV",
            "unknown_J_region_3: FGSGTRLTVV"
        ]

        self.verify_output(correct_output_lines)

    def test_cli_msa2(self):
        """
        Test the CLI with multiple sequence alignment.
        """
        args = [
            "-c", self.config_file_path,
            "-i", self.input_file_path,
            "-o", self.output_file_path,
            "msa", "needleman_wunsch"
        ]

        main(args=args)

        correct_output_lines = [
            "unknown_J_region_1: GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2: ----NTEAFFGQGTRL-TVV",
            "unknown_J_region_3: ----NYGYTFGSGTRL-TVV"
        ]

        self.verify_output(correct_output_lines)
