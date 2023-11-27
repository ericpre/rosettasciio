#!/usr/bin/env python3
import argparse
import os
import subprocess

import toml


def run_towncrier(tag):
    cmd = ("towncrier", "build", "--version", tag)

    return subprocess.call(cmd)


def update_fallback_version_in_pyproject(tag, fname="pyproject.toml"):
    pyproject = toml.load(fname)
    version = tag.strip("v").split(".")
    # Default to +1 on minor version
    major, minor = version[0], int(version[1]) + 1

    pyproject["tool"]["setuptools_scm"]["fallback_version"] = f"{major}.{minor}.dev0"

    with open(fname, "w") as f:
        toml.dump(pyproject, f)


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
