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
        DATA_PATHS = {
            "5rgs": "data/5rgs.pdb",
            "8e1y": "data/8e1y.pdb",
            "7ql8": "data/7ql8-pdb-bundle1.pdb",
            "Mpro-i0130": "data/Mpro-i0130.pdb",
            "Mpro-IBM0078":"data/refine_6.split.bound-state.pdb",
            "Mpro-IBM0058":"data/refine_7.split.bound-state.pdb",
            "Mpro-x0107":"data/refine_8.split.bound-state.pdb",
            "Mpro-IBM0045": "data/refine_16.split.bound-state.pdb",
        }

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
    pdb_paths = {key: Path(path) for key, path in constants.DATA_PATHS.items()}
    return pdb_paths


@pytest.fixture(scope="session")
def assemblies(constants, assemblies_file):
    _assemblies = {}
    with open(assemblies_file, "r") as f:
        dic = yaml.safe_load(f)

    for assembly_id, assembly_info in dic["assemblies"].items():
        _assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    return _assemblies


@pytest.fixture(scope="session")
def xtalforms(constants, assemblies_file):
    _xtalforms = {}
    with open(assemblies_file, "r") as f:
        dic = yaml.safe_load(f)["crystalforms"]

    for xtalform_id, xtalform_info in dic.items():
        _xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return _xtalforms
