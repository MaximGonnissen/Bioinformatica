import unittest

from src.fasta_parser.fasta_parser import *


class TestFastaParser(unittest.TestCase):
    """
    Tests for the fasta_parser module.
    """

    def test_parse_valid_small(self):
        """
        Test that the parser correctly parses a small valid fasta file.
        """
        fasta_file = "fasta_parser_tests/test_inputs/valid_small.fasta"
        expected = {
            "unknown_J_region_1": "GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2": "NTEAFFGQGTRLTVV",
            "unknown_J_region_3": "NYGYTFGSGTRLTVV"
        }

        results = parse(fasta_file)

        self.assertEqual(expected, results)

    def test_parse_generator_valid_small(self):
        """
        Test that the parser correctly parses a small valid fasta file when using the generator method.
        """
        fasta_file = "fasta_parser_tests/test_inputs/valid_small.fasta"
        expected = {
            "unknown_J_region_1": "GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2": "NTEAFFGQGTRLTVV",
            "unknown_J_region_3": "NYGYTFGSGTRLTVV"
        }

        results = dict(parse_generator(fasta_file))

        self.assertEqual(expected, results)

    def test_parse_generator_valid_small_alt(self):
        """
        Test that the parser correctly parses a small valid fasta file when using the generator method.

        Test that the generator yielding works correctly.
        """
        fasta_file = "fasta_parser_tests/test_inputs/valid_small.fasta"
        expected = [
            ("unknown_J_region_1", "GYSSASKIIFGSGTRLSIRP"),
            ("unknown_J_region_2", "NTEAFFGQGTRLTVV"),
            ("unknown_J_region_3", "NYGYTFGSGTRLTVV")
        ]

        for i, (seq_id, seq) in enumerate(parse_generator(fasta_file)):
            self.assertEqual(expected[i][0], seq_id)
            self.assertEqual(expected[i][1], seq)

    def test_parse_valid_large(self):
        """
        Test that the parser can parse a large valid fasta file.

        This test is not exhaustive, but it is a good indicator that the parser can handle large files.
        Only a small partial check is done.
        """
        fasta_file = "fasta_parser_tests/test_inputs/valid_large.fasta"
        expected_partial = {
            "unknown_J_region_1_0": "GYSSASKIIFGSGTRLSIRP" * 100,
            "unknown_J_region_2_1": "NTEAFFGQGTRLTVV" * 100,
            "unknown_J_region_3_2": "NYGYTFGSGTRLTVV" * 100
        }

        results = parse(fasta_file)

        for key in expected_partial:
            self.assertEqual(expected_partial[key], results[key])

    def test_parse_invalid_small(self):
        """
        Test that the parser correctly raises an exception when parsing a small invalid fasta file.
        """
        fasta_file = "fasta_parser_tests/test_inputs/invalid_small.fasta"

        with self.assertRaises(ValueError):
            parse(fasta_file)

    def test_parse_generator_invalid_small(self):
        """
        Test that the parser correctly raises an exception when parsing a small invalid fasta file when using the generator method.
        """
        fasta_file = "fasta_parser_tests/test_inputs/invalid_small.fasta"

        with self.assertRaises(ValueError):
            dict(parse_generator(fasta_file))

    def test_parse_invalid_small2(self):
        """
        Test that the parser correctly raises an exception when parsing a small invalid fasta file.
        """
        fasta_file = "fasta_parser_tests/test_inputs/invalid_small2.fasta"

        with self.assertRaises(ValueError):
            parse(fasta_file)

    def test_parse_generator_invalid_small2(self):
        """
        Test that the parser correctly raises an exception when parsing a small invalid fasta file when using the generator method.
        """
        fasta_file = "fasta_parser_tests/test_inputs/invalid_small2.fasta"

        with self.assertRaises(ValueError):
            dict(parse_generator(fasta_file))

    def test_parse_invalid_small3(self):
        """
        Test that the parser correctly raises an exception when parsing a small invalid fasta file.
        """
        fasta_file = "fasta_parser_tests/test_inputs/invalid_small3.fasta"

        with self.assertRaises(ValueError):
            parse(fasta_file)

    def test_parse_generator_invalid_small3(self):
        """
        Test that the generator variant of the parser does not raise an exception if sequence ids are not unique.

        This is because the generator variant does not check for unique sequence ids, compared to the normal parse method.
        """
        fasta_file = "fasta_parser_tests/test_inputs/invalid_small3.fasta"

        try:
            dict(parse_generator(fasta_file))
        except ValueError:
            self.fail("parse_generator raised ValueError unexpectedly!")

    def test_parse_string(self):
        """
        Test that the parser correctly parses a small valid fasta string when using the string parse method.
        """
        fasta_string = ">unknown_J_region_1\n" \
                       "GYSSASKIIFGSGTRLSIRP\n" \
                       ">unknown_J_region_2\n" \
                       "NTEAFFGQGTRLTVV\n" \
                       ">unknown_J_region_3\n" \
                       "NYGYTFGSGTRLTVV\n"

        expected = {
            "unknown_J_region_1": "GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2": "NTEAFFGQGTRLTVV",
            "unknown_J_region_3": "NYGYTFGSGTRLTVV"
        }

        results = parse_str(fasta_string)

        self.assertEqual(expected, results)

    def test_parse_string_generator(self):
        """
        Test that the parser correctly parses a small valid fasta string when using the string parse generator method.
        """
        fasta_string = ">unknown_J_region_1\n" \
                       "GYSSASKIIFGSGTRLSIRP\n" \
                       ">unknown_J_region_2\n" \
                       "NTEAFFGQGTRLTVV\n" \
                       ">unknown_J_region_3\n" \
                       "NYGYTFGSGTRLTVV\n"

        expected = {
            "unknown_J_region_1": "GYSSASKIIFGSGTRLSIRP",
            "unknown_J_region_2": "NTEAFFGQGTRLTVV",
            "unknown_J_region_3": "NYGYTFGSGTRLTVV"
        }

        results = dict(parse_str_generator(fasta_string))

        self.assertEqual(expected, results)
