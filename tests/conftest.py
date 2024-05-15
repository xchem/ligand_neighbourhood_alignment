import os
from pathlib import Path
import shutil

import yaml
import pytest

from ligand_neighbourhood_alignment import dt

@pytest.fixture(scope="session")
def constants():
    class Constants:
        TEST_DATA_DIR = "data"
        ASSEMBLIES_FILE = "data/assemblies.yaml"
        TEST_OUTPUT_DIR = "tests/output"
        DATA_PATHS = [
            "data/5rgs.pdb",
            "data/Mpro-i0130.pdb",
            "data/refine_6.split.bound-state.pdb",
            "data/refine_7.split.bound-state.pdb",
            "data/refine_8.split.bound-state.pdb",
            "data/refine_16.split.bound-state.pdb",
        ]

    return Constants()


@pytest.fixture(scope="session")
def test_data_dir(constants):
    return Path(constants.TEST_DATA_DIR)


@pytest.fixture(scope="session")
def test_output_dir(constants):
    path = Path(constants.TEST_OUTPUT_DIR)
    if path.exists():
        shutil.rmtree(path)

    return path


@pytest.fixture(scope="session")
def assemblies_file(
        constants,
):
    path = Path(constants.ASSEMBLIES_FILE)
    return path


@pytest.fixture(scope="session")
def pdb_paths(constants):
    pdb_paths = [Path(path) for path in constants.DATA_PATHS]
    return pdb_paths

@pytest.fixture(scope="session")
def assemblies(constants, assemblies_file):
    with open(assemblies_file, "r") as f:
        dic = yaml.safe_load(f)

    for assembly_id, assembly_info in dic["assemblies"].items():
        assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    return assemblies