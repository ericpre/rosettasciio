#!/usr/bin/env python3
import argparse
import re
import os
import subprocess

import toml


def run_towncrier(tag):
    cmd = ("towncrier", "build", "--version", tag)

    return subprocess.call(cmd)


def update_fallback_version_in_pyproject(tag, fname="pyproject.toml"):
    with open(fname, "r") as file:
        lines = file.readlines()

    pattern = "fallback_version"
    # Iterate through the lines and find the pattern
    for i, line in enumerate(lines):
        if re.search(pattern, line):
            lines[i] = f"{pattern} = {tag}.dev0\n"
            break

    # Write the updated content back to the file
    with open(fname, "w") as file:
        file.writelines(lines)


if __name__ == "__main__":
    # Get tag argument
    parser = argparse.ArgumentParser()
    parser.add_argument("tag")
    args = parser.parse_args()
    tag = args.tag

    # Update release notes
    run_towncrier(tag)

    # Update fallback version for setuptools_scm
    update_fallback_version_in_pyproject(tag)
