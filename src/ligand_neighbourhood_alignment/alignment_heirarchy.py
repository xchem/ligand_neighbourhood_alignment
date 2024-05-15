from rich import print as rprint

from ligand_neighbourhood_alignment import dt

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
    assembly_priority = {_j: _assembly_name for _j, _assembly_name in enumerate(assemblies)}

    # 2. Determine the Chain priority and map assembly names to chains
    chain_priority = {}
    assembly_chains = {}
    chain_priority_count = 0
    for _j, _assembly_name in assembly_priority.items():
        assembly = assemblies[_assembly_name]
        assembly_chains[_assembly_name] = []
        for _generator in assembly.generators:
            _biological_chain_name = _generator.biomol
            assembly_chains[_assembly_name].append(_biological_chain_name)
            if _biological_chain_name not in chain_priority.values():
                chain_priority[_biological_chain_name] = chain_priority_count
                chain_priority_count += 1

    if debug:
        rprint(chain_priority)


    # 3. Find each assembly's reference
    reference_assemblies = {}
    for _assembly_name, _assembly in assemblies.items():
        # Get the highest priority chain
        reference_chain = min(
            [_generator.chain for _generator in _assembly.generators], key=lambda _x: chain_priority[_x]
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

    # 4. Check per-chain RMSDs and warn if any are high
    # TODO

    return reference_assemblies


def _chain_to_biochain(chain_name, xtalform: dt.XtalForm, assemblies: dict[str, dt.Assembly]) -> str:
    for _xtal_assembly_name, _xtal_assembly in xtalform.assemblies.items():
        for _j, _chain_name in enumerate(_xtal_assembly.chains):
            if chain_name == _chain_name:
                return assemblies[_xtal_assembly.assembly].generators[_j].biomol


StructureLandmarks = dict[tuple[str, str, str], tuple[float, float, float]]


def _calculate_assembly_transform(
    assembly_name: str, alignment_heirarchy: AlignmentHeirarchy, assembly_landmarks: dict[str, StructureLandmarks]
):
    # Get the chain to align to
    ...
