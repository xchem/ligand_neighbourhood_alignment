from pathlib import Path
from typing import Generator

import gemmi
import networkx as nx
from pydantic import BaseModel, validator

from xchemalign import constants

Structure = gemmi.Structure


class DatasetID(BaseModel):
    dtag: str


class LigandID(BaseModel):
    dtag: str
    chain: str
    id: int

    def __eq__(self, other) -> bool:

        try:

            if self.dtag == other.dtag:
                if self.id == other.id:
                    if self.chain == other.chain:
                        return True
            return False
        except Exception:
            return False

    def __hash__(self):
        return hash((self.dtag, self.chain, self.id))

    def to_string(
        self,
    ):
        return f"{self.dtag}~{self.chain}~{self.id}"

    @classmethod
    def from_string(cls, string):
        dtag, chain, id = string.split("~")

        return LigandID(dtag=dtag, chain=chain, id=int(id))


class SymOp(BaseModel):
    operation: str
    image: tuple[int, int, int]


class AtomID(BaseModel):
    chain: str
    residue: int
    atom: str

    def __eq__(self, other) -> bool:
        if self.atom == other.atom:
            if self.residue == other.residue:
                if self.chain == other.chain:
                    return True

        return False

    def __hash__(self):
        return hash((self.chain, self.residue, self.atom))


class Transform(BaseModel):
    vec: list[float]
    mat: list[list[float]]


class Transforms(BaseModel):
    ligand_ids: list[tuple[LigandID, LigandID]]
    transforms: list[Transform]

    def get_transform(self, transform_id: tuple[LigandID, LigandID]):
        for ligand_id, _transform in zip(self.ligand_ids, self.transforms):
            if transform_id == ligand_id:

                transform = gemmi.Transform()
                transform.vec.fromlist(_transform.vec)
                transform.mat.fromlist(_transform.mat)

                return transform

        raise Exception(f"Transform {transform_id} not in transforms!")


class Atom(BaseModel):
    element: str
    atom_id: AtomID
    x: float
    y: float
    z: float
    image: Transform


class AlignableSite(BaseModel):
    id: int
    name: str
    ligand_ids: list[LigandID]
    reference: LigandID


class CanonicalSite(BaseModel):
    id: int
    name: str
    refpdb: str
    atoms: dict[AtomID, Atom]
    literatureref: str
    members: list[LigandID]


class XtalForm(BaseModel):
    id: int
    space_group: int
    unit_cell: tuple[float, float, float, float, float, float]
    members: list[DatasetID]
    transforms: list[Transform]


class XtalForms(BaseModel):
    xtalforms: list[XtalForm]


class XtalFormSite(BaseModel):
    id: int
    canon_site_id: int
    xtal_form_id: int
    code: str
    refpdb: str
    atoms: dict[AtomID, Atom]
    artefact_atoms: dict[AtomID, Atom]
    members: list[LigandID]


class SiteObservation(BaseModel):
    id: int
    ligand_id: LigandID
    xtal_form_site_id: int
    fragalysis_label: str
    description: str
    code: str
    dataset: str
    # compound: int


class LigandBindingEvent(BaseModel):
    id: int
    chain: str
    residue: int
    xmap: str


class LigandBindingEvents(BaseModel):
    ligand_ids: list[LigandID]
    ligand_binding_events: list[LigandBindingEvent]

    def __getitem__(self, lid: LigandID) -> LigandBindingEvent:
        for _lid, lbe in zip(self.ligand_ids, self.ligand_binding_events):
            if lid == _lid:
                return lbe

        raise Exception(f"{lid} : {self.ligand_ids}")


class Dataset(BaseModel):
    dtag: str
    pdb: str
    ligand_binding_events: LigandBindingEvents
    # mtz: str


class SystemData(BaseModel):
    dataset_ids: list[DatasetID]
    datasets: list[Dataset]

    def get_dataset(self, did: str | DatasetID) -> Dataset:
        for _did, dataset in zip(self.dataset_ids, self.datasets):
            if _did == did:
                return dataset

        raise Exception(f"{did} : {self.dataset_ids}")


class LigandNeighbourhood(BaseModel):
    atom_ids: list[AtomID]
    atoms: list[Atom]
    artefact_atom_ids: list[AtomID]
    artefact_atoms: list[Atom]


class LigandNeighbourhoods(BaseModel):
    ligand_ids: list[LigandID]
    ligand_neighbourhoods: list[LigandNeighbourhood]

    def get_neighbourhood(self, ligand_id: LigandID):
        for _ligand_id, _neighbourhood in zip(
            self.ligand_ids, self.ligand_neighbourhoods
        ):
            if _ligand_id == ligand_id:
                return _neighbourhood


class SystemSites(BaseModel):
    canonical_site: dict[int, CanonicalSite]
    xtal_form_site: dict[int, XtalFormSite]
    ligand_ids: list[LigandID]
    site_observation: dict[int, SiteObservation]

    @validator("canonical_site")
    def check_canonical_site_ids(cls, v: dict[LigandID, SiteObservation]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.id

    @validator("canonical_site")
    def check_canonical_site_ids_sequential(
        cls,
        v: dict[LigandID, SiteObservation],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums

    @validator("xtal_form_site")
    def check_xtal_form_site_ids(cls, v: dict[int, XtalFormSite]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.id

    @validator("xtal_form_site")
    def check_xtal_form_site_ids_sequential(
        cls,
        v: dict[int, XtalFormSite],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums

    @validator("site_observation")
    def check_site_observation_ids(cls, v: dict[LigandID, SiteObservation]):
        if not v:
            return
        for site_id, site in v.items():
            assert site_id == site.ligand_id

    @validator("site_observation")
    def check_site_observation_ids_sequential(
        cls,
        v: dict[LigandID, SiteObservation],
    ):
        if not v:
            return
        num_sites: int = len(v)
        site_nums = [site.id for site in v.values()]
        for site_num in range(num_sites):
            assert site_num in site_nums


# class Site(BaseModel):
#     members: list[LigandID]


class Block(BaseModel):
    xi: int
    yi: int
    zi: int
    dx: int
    dy: int
    dz: int
    transform: Transform


class ResidueID(BaseModel):
    chain: str
    residue: int

    def __hash__(self):
        return hash((self.chain, self.residue))


class SubSite(BaseModel):
    id: int
    name: str
    residues: list[ResidueID]
    members: list[LigandID]


class SubSites(BaseModel):
    subsites: list[SubSite]

    def iter(self) -> Generator[tuple[int, SubSite], None, None]:
        for subsite in self.subsites:
            yield subsite.id, subsite


class Site(BaseModel):
    id: int
    subsite_ids: list[int]
    subsites: list[SubSite]
    members: list[LigandID]
    residues: list[ResidueID]

    def iter(self) -> Generator[tuple[int, SubSite], None, None]:
        for subsite_id, subsite in zip(self.subsite_ids, self.subsites):
            yield subsite_id, subsite


class Sites(BaseModel):
    site_ids: list[int]
    sites: list[Site]

    def iter(self) -> Generator[tuple[int, Site], None, None]:
        for site_id, site in zip(self.site_ids, self.sites):
            yield site_id, site


def read_xmap(path: Path):
    m = gemmi.read_ccp4_map(str(path), setup=True)
    return m.grid


def transform_to_gemmi(transform: Transform):
    transform_gemmi = gemmi.Transform()
    transform_gemmi.vec.fromlist(transform.vec)
    transform_gemmi.mat.fromlist(transform.mat)

    return transform_gemmi


def gemmi_to_transform(transform):
    return Transform(vec=transform.vec.tolist(), mat=transform.mat.tolist())


def get_box(neighbourhood: LigandNeighbourhood, xmap, transform: Transform):

    # transform_gemmi = transform_to_gemmi(transform)
    transform_gemmi = transform

    box = gemmi.FractionalBox()
    for atom in neighbourhood.atoms:
        box.extend(
            xmap.cell.fractionalize(
                transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
            )
        )

    for atom in neighbourhood.artefact_atoms:
        box.extend(
            xmap.cell.fractionalize(
                transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
            )
        )

    return box


def write_xmap(
    xmap, path: Path, neighbourhood: LigandNeighbourhood, transform: Transform
):

    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = xmap
    box = get_box(neighbourhood, xmap, transform)
    ccp4.set_extent(box)
    ccp4.update_ccp4_header()

    ccp4.write_ccp4_map(str(path))


def read_graph(path: Path):
    g = nx.read_gml(
        str(path / constants.ALIGNABILITY_GRAPH_FILE_NAME),
        destringizer=lambda x: LigandID.from_string(x),
    )

    return g


def read_neighbourhoods(path: Path):
    neighbourhoods = LigandNeighbourhoods.parse_file(
        str(path / constants.NEIGHBOURHOODS_FILE_NAME)
    )
    return neighbourhoods


def read_sites(path: Path):
    sites = Sites.parse_file(str(path / constants.SITES_FILE_NAME))

    return sites


def read_transforms(path: Path):
    transforms = Transforms.parse_file(
        str(path / constants.TRANSFORMS_FILE_NAME)
    )
    return transforms


def read_structures(system_data: SystemData):
    structures = {}
    for dataset in system_data.datasets:
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    return structures


def read_system_data(path: Path):
    return SystemData.parse_file(str(path / constants.DATA_JSON_PATH))


class SiteTransforms(BaseModel):
    site_transform_ids: list[tuple[int, int]]
    site_transforms: list[Transform]
    subsite_transform_ids: list[tuple[int, int, int]]
    subsite_transforms: list[Transform]

    def get_subsite_transform(self, site_id, subsite_id):
        for subsite_transform_id, subsite_transform in zip(
            self.subsite_transform_ids, self.subsite_transforms
        ):
            if subsite_transform_id[0] == site_id:
                if subsite_transform_id[2] == subsite_id:
                    return subsite_transform

        raise Exception()

    def get_site_transform(
        self,
        site_id,
    ):
        for site_transform_id, subsite_transform in zip(
            self.site_transform_ids, self.site_transforms
        ):
            if site_transform_id[1] == site_id:
                return subsite_transform

        raise Exception()


def save_site_transforms(site_transforms: SiteTransforms, path: Path):
    with open(path / constants.SITES_TRANSFORMS_FILE_NAME, "w") as f:
        f.write(site_transforms.json())


def read_site_transforms(path: Path):
    return SiteTransforms.parse_file(
        str(path / constants.SITES_TRANSFORMS_FILE_NAME)
    )
