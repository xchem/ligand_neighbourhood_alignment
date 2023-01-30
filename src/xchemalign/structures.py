import gemmi

from xchemalign.data import LigandNeighbourhood, ResidueID, SystemData
from xchemalign.matching import match_atom


def get_structures(system_data: SystemData):
    structures = {}
    for dataset in system_data.datasets:
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    return structures


def _get_transform():
    ...


def get_transform_from_residues(rs: list[ResidueID], srs, ssrs):
    # Transform from ssrs to srs
    acs = []
    for resid in rs:
        chain, num = resid.chain, resid.residue
        srsca = srs[chain][num][0]["CA"][0]
        ssrsca = ssrs[chain][num][0]["CA"][0]
        acs.append((srsca, ssrsca))

    sup = gemmi.superpose_positions([x[0] for x in acs], [x[1] for x in acs])

    return sup.transform


def get_transforms(
    reference_neighbourhood: LigandNeighbourhood,
    neighbourhood: LigandNeighbourhood,
):

    alignable_cas = {}
    for (
        ligand_1_atom_id,
        ligand_1_atom,
    ) in zip(reference_neighbourhood.atom_ids, reference_neighbourhood.atoms):
        for (
            ligand_2_atom_id,
            ligand_2_atom,
        ) in zip(neighbourhood.atom_ids, neighbourhood.atoms):
            if ligand_1_atom_id.atom == "CA":
                if match_atom(ligand_1_atom, ligand_2_atom, ignore_chain=True):
                    alignable_cas[ligand_1_atom_id] = (
                        gemmi.Position(
                            ligand_1_atom.x,
                            ligand_1_atom.y,
                            ligand_1_atom.z,
                        ),
                        gemmi.Position(
                            ligand_2_atom.x,
                            ligand_2_atom.y,
                            ligand_2_atom.z,
                        ),
                    )

    if len(alignable_cas) < 3:
        return gemmi.Transform(), []

    sup = gemmi.superpose_positions(
        [alignable_ca[0] for alignable_ca in alignable_cas.values()],
        [alignable_ca[1] for alignable_ca in alignable_cas.values()],
    )
    # logger.debug(f"Superposition: rmsd {sup.rmsd} n {len(alignable_cas)}")

    return sup.transform, [ligand_id for ligand_id in alignable_cas.keys()]