import gemmi
from ligand_neighbourhood_alignment import alignment_heirarchy, dt


def main():
    dtag = 'Mpro-x0107'
    file = 'data/refine_8.split.bound-state.pdb'
    st = gemmi.read_structure(file)
    generators = [
        dt.Generator('A', 'A', 'x,y,z'),
        dt.Generator('B', 'A', '-x,y,-z'),
    ]

    # Setup new structure to add biochains to
    new_st = gemmi.Structure()
    new_st.spacegroup_hm = 'P 1'
    new_st.cell = st.cell
    new_model = gemmi.Model("0")
    new_st.add_model(new_model)

    # Iterate over chain, biochain, transform tuples in the assembly
    for generator in generators:
        # Generate the symmetry operation
        op = gemmi.Op(generator.triplet)

        # Create a clone of base chain to transform
        new_chain = st[0][generator.chain].clone()

        # Transform the residues
        for residue in new_chain:
            for atom in residue:
                atom_frac = st.cell.fractionalize(atom.pos)
                new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                new_pos_orth = st.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))
                atom.pos = gemmi.Position(*new_pos_orth)
        new_chain.name = generator.biomol
        new_st[0].add_chain(new_chain)

    new_st.setup_entities()
    print(new_st.cell)
    print(new_st.spacegroup_hm)
    new_st.write_pdb('Mpro-x0107_fake_P1.pdb')


if __name__ == "__main__":
    main()
