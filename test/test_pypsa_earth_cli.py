# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import importlib
import sys
import unittest
from contextlib import ExitStack
from unittest import mock


# The CLI reconfigures stdout/stderr on Windows at import; make that a no-op here.
class _StreamWithoutReconfigure:
    def __init__(self, stream):
        self._stream = stream

    def reconfigure(self, *args, **kwargs):
        pass

    def __getattr__(self, name):
        return getattr(self._stream, name)


with ExitStack() as stack:
    if sys.platform == "win32":
        for name in ("stdout", "stderr"):
            stream = _StreamWithoutReconfigure(getattr(sys, name))
            stack.enter_context(mock.patch.object(sys, name, stream))
    _cli = importlib.import_module("scripts.non_workflow.pypsa_earth_cli")

coerce_value = _cli.coerce_value
prompt_config_value = _cli.prompt_config_value


class TestCoerceValue(unittest.TestCase):
    def test_bool(self):
        cases = (
            ("true", True),
            ("TRUE", True),
            ("yes", True),
            (" yes ", True),
            ("false", False),
            ("FALSE", False),
            ("no", False),
            (" no ", False),
        )
        for value, expected in cases:
            with self.subTest(value=value):
                self.assertIs(coerce_value(value, expected), expected)
        for value in ("1", "maybe", ""):
            with self.subTest(value=value), self.assertRaises(ValueError):
                coerce_value(value, True)

    def test_int(self):
        for value, expected in (("42", 42), ("-3", -3), (" 7 ", 7)):
            with self.subTest(value=value):
                self.assertEqual(coerce_value(value, 0), expected)
        for value in ("3.5", "abc", ""):
            with self.subTest(value=value), self.assertRaises(ValueError):
                coerce_value(value, 0)

    def test_float(self):
        for value, expected in (("3.5", 3.5), ("5", 5.0), (" -3 ", -3.0)):
            with self.subTest(value=value):
                self.assertEqual(coerce_value(value, 0.0), expected)
        for value in ("abc", ""):
            with self.subTest(value=value), self.assertRaises(ValueError):
                coerce_value(value, 0.0)

    def test_str(self):
        for value in ("hello", "", "  padded  ", "false"):
            with self.subTest(value=value):
                self.assertEqual(coerce_value(value, "existing"), value)

    def test_list(self):
        cases = (
            ("10,20", [10], [10, 20]),
            ("10, 20", [10], [10, 20]),
            ("2,3.5", [1.5], [2.0, 3.5]),
            ("true,no", [True], [True, False]),
            ("NG, BJ", ["NG"], ["NG", "BJ"]),
            ("a,,b", [], ["a", "", "b"]),
            ("a,b", [], ["a", "b"]),
        )
        for value, expected, result in cases:
            with self.subTest(value=value):
                self.assertEqual(coerce_value(value, expected), result)
        invalid = (("10, x", [10]), ("2, x", [1.5]), ("maybe", [True]), ("1,", [10]))
        for value, expected in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                coerce_value(value, expected)
        for expected in ([], [""], [10], ["NG"], [1.5], [True]):
            with self.subTest(expected=expected):
                self.assertEqual(coerce_value("", expected), [])
        expected = [10]
        coerce_value("10,20", expected)
        self.assertEqual(expected, [10])

    def test_none(self):
        for value in ("anything", "", "true", "   "):
            with self.subTest(value=value):
                self.assertEqual(coerce_value(value, None), value)

    def test_unsupported_expected_type(self):
        self.assertEqual(coerce_value("value", {}), "value")


class TestPromptConfigValue(unittest.TestCase):
    def test_valid_value(self):
        with mock.patch.object(_cli, "ask", return_value="42"):
            self.assertEqual(prompt_config_value("scenario.clusters", 10), 42)

    def test_retry_then_success(self):
        with (
            mock.patch.object(_cli, "ask", side_effect=["abc", "42"]) as ask,
            mock.patch.object(_cli.console, "print") as output,
        ):
            self.assertEqual(prompt_config_value("scenario.clusters", 10), 42)
        self.assertEqual(ask.call_count, 2)
        self.assertEqual(output.call_count, 1)

    def test_error_message_names_key_and_bad_value(self):
        with (
            mock.patch.object(_cli, "ask", side_effect=["abc", "42"]),
            mock.patch.object(_cli.console, "print") as output,
        ):
            prompt_config_value("scenario.clusters", 10)
        message = output.call_args.args[0]
        self.assertIn("scenario.clusters", message)
        self.assertIn("abc", message)

    def test_markup_input_retries(self):
        with (
            mock.patch.object(_cli, "ask", side_effect=["[/foo]", "42"]),
            mock.patch.object(_cli.console, "print") as output,
        ):
            self.assertEqual(prompt_config_value("scenario.clusters", 10), 42)
        self.assertIn("\\[/foo]", output.call_args.args[0])

    def test_simpl_empty_input(self):
        with mock.patch.object(_cli, "ask", return_value=""):
            self.assertEqual(prompt_config_value("scenario.simpl", ["a"]), [""])

    def test_polymorphic_false(self):
        cases = (
            ("false", False),
            ("FALSE", False),
            (" false ", False),
            ("standard", "standard"),
            ("api", "api"),
        )
        for key in (
            "build_shape_options.worldpop_method",
            "build_shape_options.gdp_method",
        ):
            for value, expected in cases:
                with (
                    self.subTest(key=key, value=value),
                    mock.patch.object(_cli, "ask", return_value=value),
                ):
                    self.assertEqual(prompt_config_value(key, "standard"), expected)


if __name__ == "__main__":
    unittest.main()
