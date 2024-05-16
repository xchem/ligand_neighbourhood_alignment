from pathlib import Path

import pytest
import yaml
from rich import print as rprint
import gemmi

from ligand_neighbourhood_alignment import alignment_heirarchy


def test_derive_alignment_heirarchy(
        constants,
        assemblies
):
    # Run _derive_alignment_heirarchy
    reference_assemblies = alignment_heirarchy._derive_alignment_heirarchy(assemblies, debug=True)

    # Print important info
    # rprint(reference_assemblies)
    ...


def test_chain_to_biochain(
        constants,
        assemblies,
        xtalforms
):
    rprint(xtalforms)
    for xtalform_name, xtalform in xtalforms.items():
        for xtalform_assembly_name, xtalform_assembly in xtalform.assemblies.items():
            for chain in xtalform_assembly.chains:
                biochain = alignment_heirarchy._chain_to_biochain(chain, xtalform, assemblies)
                rprint(
                    f'Chain to Biochain: {xtalform_name} : {xtalform_assembly_name} : {chain} -> {biochain}'
                )


def test_structure_to_landmarks(
        pdb_paths
):
    st = gemmi.read_structure(str(pdb_paths['Mpro-IBM0045']))
    landmarks = alignment_heirarchy.structure_to_landmarks(st)
    rprint(landmarks)


def test_calculate_assembly_transform(
        pdb_paths,
        assemblies
):
    # Generate assembly structure from reference
    as1_ref = gemmi.read_structure(str(pdb_paths['Mpro-IBM0045']))
    as1 = alignment_heirarchy._get_assembly_st(assemblies['dimer'], as1_ref)

    as2_ref = gemmi.read_structure(str(pdb_paths['7ql8']))
    as2 = alignment_heirarchy._get_assembly_st(assemblies['monomer'], as2_ref)

    # Generate assembly landmarks from assembly structure
    as1_lm = alignment_heirarchy.structure_to_landmarks(as1)
    as2_lm = alignment_heirarchy.structure_to_landmarks(as2)

    # Generate alignment hierarchy
    hierarchy = alignment_heirarchy._derive_alignment_heirarchy(assemblies)

    # Determine transform
    transform = alignment_heirarchy._calculate_assembly_transform(
        ref=as1_lm,
        mov=as2_lm,
        chain=hierarchy['monomer'][1],
        debug=True
    )
    rprint(transform)

def test_calculate_assembly_transform_sequence(assemblies):

    # Generate alignment hierarchy
    hierarchy = alignment_heirarchy._derive_alignment_heirarchy(assemblies)

    #
    transform_sequence = alignment_heirarchy._calculate_assembly_transform_sequence(hierarchy, 'fake_tetramer')
    rprint(transform_sequence)


# @pytest.mark.order(after="test_collator_upload_1")
# def test_aligner_upload_1(constants, assemblies_file, upload_1_dir):
#     a = Aligner(upload_1_dir, constants.METADATA_FILE, assemblies_file)
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         print(a.logger.errors)
#         exit(1)
#     else:
#         a.run()
#
#
# @pytest.mark.order(after="test_aligner_upload_1")
# def test_collator_upload_2(
#     constants,
#     test_dir,
#     upload_2_dir,
#     upload_2_data_dir,
#     config_2_file,
# ):
#     c = Collator(
#         str(config_2_file),
#     )
#     logger = c.logger
#
#     meta = c.validate()
#
#     if meta is None or len(logger.errors) > 0:
#         print("There are errors, cannot continue")
#         print(logger.errors)
#         exit(1)
#     else:
#         c.run(meta)
#     assert (
#         len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {}))
#         != 0
#     )
#
#
# @pytest.mark.order(after="test_collator_upload_2")
# def test_aligner_upload_2(constants, assemblies_file, upload_2_dir):
#     a = Aligner(upload_2_dir, constants.METADATA_FILE, assemblies_file)
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         print(a.logger.errors)
#         exit(1)
#     else:
#         a.run()
#     assert "Mpro-i0130" in [x.name for x in (Path(upload_2_dir) / "aligned_files").glob("*")]


# @pytest.mark.order(after="test_aligner_upload_2")
# def test_collator_upload_3(
#     constants,
#     test_dir,
#     upload_3_dir,
#     upload_3_data_dir,
#     config_3_file,
# ):
#     c = Collator(str(config_3_file), )
#
#     meta, num_errors, num_warnings = c.validate()
#
#     if meta is None or num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         c.run(meta)
#
#     with open(Path(upload_3_dir) / "meta_collator.yaml", 'r') as f:
#         new_meta = yaml.safe_load(f)
#     assert len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})) != 0
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "crystallographic_files").glob("*")]
#
#
# @pytest.mark.order(after="test_collator_upload_3")
# def test_aligner_upload_3(constants, xtalforms_file, upload_3_dir):
#     a = Aligner(upload_3_dir, constants.METADATA_FILE, xtalforms_file, )
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         a.run()
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "aligned_files").glob("*")]
