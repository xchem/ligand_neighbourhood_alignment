import subprocess
import sys

from ligand_neighbourhood_alignment import __version__


def test_cli_version():
    cmd = [sys.executable, "-m", "ligand_neighbourhood_alignment", "--version"]
    assert subprocess.check_output(cmd).decode().strip() == __version__
