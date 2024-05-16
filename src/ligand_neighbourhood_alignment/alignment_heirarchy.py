from rich import print as rprint
import gemmi

from ligand_neighbourhood_alignment import dt, constants

AlignmentHeirarchy = dict[str, tuple[str, str]]


def _derive_alignment_heirarchy(assemblies: dict[str, dt.Assembly], debug=False) -> AlignmentHeirarchy:
    # The Alignment hierarchy is the graph of alignments one must perform in order to get from
    # a ligand canonical site to the Reference Assembly Frame

    # In order to calculate the assembly the following steps are performed:
    # 1. Determine the Assembly priority
    # 2. Determine the Chain priority
    # 3. Find each assembly's reference
    # 4. Check per-chain RMSDs and warn if any are high

    # 1. Determine the Assembly priority
    assembly_priority = {_assembly_name: _j for _j, _assembly_name in enumerate(assemblies)}

    # 2. Determine the Chain priority and map assembly names to chains
    chain_priority = {}
    assembly_chains = {}
    chain_priority_count = 0
    for _assembly_name, _j in assembly_priority.items():
        assembly = assemblies[_assembly_name]
        assembly_chains[_assembly_name] = []
        for _generator in assembly.generators:
            _biological_chain_name = _generator.biomol
            assembly_chains[_assembly_name].append(_biological_chain_name)
            if _biological_chain_name not in chain_priority:
                chain_priority[_biological_chain_name] = chain_priority_count
                chain_priority_count += 1

    if debug:
        rprint(f'Assembly priority')
        rprint(assembly_priority)
        rprint(f'Chain priority')
        rprint(chain_priority)

    # 3. Find each assembly's reference
    reference_assemblies = {}
    for _assembly_name, _assembly in assemblies.items():
        # Get the highest priority chain
        reference_chain = min(
            [_generator.biomol for _generator in _assembly.generators], key=lambda _x: chain_priority[_x]
        )

        # Get the highest priority assembly in which it occurs
        reference_assembly = min(
            [
                _assembly_name
                for _assembly_name in assembly_chains
                if reference_chain in assembly_chains[_assembly_name]
            ],
            key=lambda _x: assembly_priority[_x],
        )
        reference_assemblies[_assembly_name] = (reference_assembly, reference_chain)

    if debug:
        rprint(f'Reference Assemblies')
        rprint(reference_assemblies)

    # 4. Check per-chain RMSDs and warn if any are high
    # TODO

    return reference_assemblies


def _chain_to_biochain(chain_name, xtalform: dt.XtalForm, assemblies: dict[str, dt.Assembly]) -> str:
    for _xtal_assembly_name, _xtal_assembly in xtalform.assemblies.items():
        for _j, _chain_name in enumerate(_xtal_assembly.chains):
            if chain_name == _chain_name:
                return assemblies[_xtal_assembly.assembly].generators[_j].biomol


StructureLandmarks = dict[tuple[str, str, str], tuple[float, float, float]]


def structure_to_landmarks(st):
    landmarks = {}
    for model in st:
        for chain in model:
            for residue in chain:
                if residue.name not in constants.RESIDUE_NAMES:
                    continue
                for atom in residue:
                    pos = atom.pos
                    landmarks[(chain.name, (residue.name, str(residue.seqid.num)), atom.name)] = (pos.x, pos.y, pos.z)

    return landmarks


def _get_assembly_st(as1, as1_ref):
    # Setup new structure to add biochains to
    new_st = gemmi.Structure()
    new_model = gemmi.Model("0")
    new_st.add_model(new_model)

    # Iterate over chain, biochain, transform tuples in the assembly
    for generator in as1.generators:
        # Generate the symmetry operation
        op = gemmi.Op(generator.triplet)

        # Create a clone of base chain to transform
        new_chain = as1_ref[0][generator.chain].clone()

        # Transform the residues
        for residue in new_chain:
            for atom in residue:
                atom_frac = as1_ref.cell.fractionalize(atom.pos)
                new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                new_pos_orth = as1_ref.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))
                atom.pos = gemmi.Position(*new_pos_orth)
        new_chain.name = generator.biomol
        new_st[0].add_chain(new_chain)

    new_st.add_model(new_model)
    return new_st


def _landmark_to_structure(lm):
    st = gemmi.Structure()
    model = gemmi.Model("0")
    st.add_model(model)

    used_chains = []
    used_ress = []
    for (chain, (res_name, seqid), atom), (x, y, z) in lm.items():
        if chain not in used_chains:
            st[0].add_chain(gemmi.Chain(chain))
            used_chains.append(chain)

        if (chain, seqid) not in used_ress:
            new_residue = gemmi.Residue()
            new_residue.name = res_name
            new_residue.seqid = gemmi.SeqId(str(seqid))
            st[0][chain].add_residue(new_residue)
            used_ress.append((chain, seqid))

        new_atom = gemmi.Atom()
        new_atom.name = atom
        pos = gemmi.Position(x, y, z)
        new_atom.pos = pos
        st[0][chain][seqid][0].add_atom(new_atom)

    st.setup_entities()
    return st

    ...

def _get_atoms(st):
    atoms = {}
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    atoms[(chain, res.name, res.seqid.num, atom.name)] = atom

    return atoms

def _calculate_assembly_transform(
        ref=None,
        mov=None,
        chain=None,
        debug=False
):
    # Convert to gemmi structures to use superposition algorithm there
    sup = gemmi.superpose_positions(
        [gemmi.Position(x, y, z) for atom_id, (x, y, z ) in ref.items() if (atom_id[0] == chain) & (atom_id in mov) & (atom_id[2] == 'CA')],
        [gemmi.Position(x, y, z) for atom_id, (x, y, z) in mov.items() if (atom_id[0] == chain) & (atom_id in ref)& (atom_id[2] == 'CA')]
    )
    transform = sup.transform

    # transform to interchangable format and return
    return {
        'vec': transform.vec.tolist(),
        'mat': transform.mat.tolist(),
        'rmsd': sup.rmsd,
        "count": sup.count
    }
