import io
import numpy as np
import pathlib
from urllib.parse import urlparse
from numbers import Number

class MiomNetwork:
    """A minimal class to store a metabolic network.

    Attributes:
        S (numpy.ndarray): stoichiometric matrix
        R (numpy.ndarray): Structured array with the reactions. The fields are:
            - id (int): reaction ID
            - lb (float): lower bound
            - ub (float): upper bound
            - subsystem (str): subsystem
            - gpr (str): gene-protein-reaction rule
        M (numpy.ndarray): Structured array with the metabolites. The fields are:
            - id (int): metabolite ID
            - name (str): metabolite name
            - formula (str): metabolite formula
    """
    def __init__(self, S, R, M):
        if S.shape[1] != R.shape[0]:
            raise ValueError("Num of reactions in the stoichiometric matrix "
                             "are different from the number of reaction in R")
        self.S = S
        self.R = R
        self.M = M

    @property
    def num_reactions(self):
        """Number of reactions in the network.

        Returns:
            int: Number of reactions
        """
        return self.R.shape[0]

    @staticmethod
    def _find_reaction(rxn_id, R):
        if isinstance(rxn_id, Number):
            return rxn_id, R[rxn_id]
        for i, r in enumerate(R):
            if r['id'] == rxn_id:
                return i, r
        raise ValueError(f"Cannot find reaction {rxn_id}")

    def find_reaction(self, rxn_id):
        """Find a particular reaction in the metabolic network.

        Args:
            rxn_id (str): Name of the reaction

        Returns:
            numpy.ndarray: Structured array with the information of the reaction.
        """
        return MiomNetwork._find_reaction(rxn_id, self.R)

    def get_reaction_id(self, rxn):
        i, r = self.find_reaction(rxn)
        return i

    def find_reactions(self, rxn_ids):
        return [self.find_reaction(rxn_id) for rxn_id in rxn_ids]

    def find_reactions_from_pathway(self, pathway_name):
        return np.array([1 if pathway_name.lower() in rxn['subsystem'].lower() else 0 for rxn in self.R])

    @property
    def object_size(self):
        bytes = self.S.__sizeof__() + self.R.__sizeof__() + self.M.__sizeof__()
        return bytes / 1024**2




def _is_url(url):
    """
    Determine if the provided string is a valid url
    :param url: string
    :return: True if the string is a URL
    """
    if isinstance(url, str):
        try:
            result = urlparse(url)
            return all([result.scheme, result.netloc])
        except ValueError:
            return False
    return False


def _download(url_file):
    if not _is_url(url_file):
        raise ValueError("Invalid url")
    import tempfile
    import os
    from urllib.request import urlopen
    ext = pathlib.Path(url_file).suffix
    path = os.path.join(tempfile.mkdtemp(), "file" + ext)
    with urlopen(url_file) as rsp, open(path, 'wb') as output:
        output.write(rsp.read())
    return path


def load_gem(model_or_path):
    """Load a metabolic network from a file or URL.

    The method supports any format supported by cobrapy (.xml, .yml, .json, .mat)
    or a miom compressed model (.miom) from a url or a local file path. For the cobra
    supported formats, you need the cobrapy package installed, and for .mat files, you
    need both cobrapy and scipy installed.
    
    Args:
        model_or_path (str): Path to a local file or URL pointing to the metabolic network

    Returns:
        MiomNetwork: A [MiomNetwork][miom.mio] instance with the minimal information 
            of the metabolic network required for simulations. It includes the stoichiometric
            matrix, the list of reactions with the lower and upper bounds, the associated
            genes and GPR rules, and the list of metabolites.
    """
    if isinstance(model_or_path, str):
        file = model_or_path
        if _is_url(model_or_path):
            file = _download(model_or_path)
        ext = pathlib.Path(file).suffix
        if ext == '.miom' or ext == '.xz' or ext == '.npz':
            return _load_compressed_model(file)
        else:
            return _cobra_to_miom(_read_cobra_model(file))
    else:
        return _cobra_to_miom(model_or_path)


def _read_cobra_model(filepath):
    ext = pathlib.Path(filepath).suffix
    if ext == '.mat':
        from cobra.io import load_matlab_model
        return load_matlab_model(filepath)
    elif ext == '.xml':
        from cobra.io import read_sbml_model
        return read_sbml_model(filepath)
    elif ext == '.json':
        from cobra.io import load_json_model
        return load_json_model(filepath)
    elif ext == '.yml':
        from cobra.io import load_yaml_model
        return load_yaml_model(filepath)
    else:
        raise ValueError("Unsupported file format")


def export_gem(miom_network, path_to_exported_file):
    """Export a miom network to a file in the miom format.

    Args:
        miom_network (MiomNetwork): an instance of a MiomNetwork
        path_to_exported_file (str): Path to the exported file (e.g. /path/to/file.miom)
    """
    import lzma
    with io.BytesIO() as npz:
        np.savez_compressed(npz,
                            S=miom_network.S,
                            reactions=miom_network.R,
                            metabolites=miom_network.M)
        compressed = lzma.compress(npz.getbuffer())
        # Store to file
        with open(path_to_exported_file, 'wb') as f_out:
            f_out.write(compressed)


def _cobra_to_miom(model):
    try:
        from cobra.util.array import create_stoichiometric_matrix
    except ImportError as e:
        raise ImportError("Cobrapy package is not installed, "
                          "but required to read and import metabolic networks", e)
    S = create_stoichiometric_matrix(model)
    rxn_data = [(rxn.id, rxn.lower_bound, rxn.upper_bound, rxn.subsystem, rxn.gene_reaction_rule)
                for rxn in model.reactions]
    met_data = [(met.id, met.name, met.formula) for met in model.metabolites]
    id_max_length = max(len(rxn.id) for rxn in model.reactions)
    subsyst_max_length = max((len(rxn.subsystem)) for rxn in model.reactions)
    metid_max_length = max((len(met.id)) for met in model.metabolites)
    metname_max_length = max((len(met.name)) for met in model.metabolites)
    metform_max_length = max((len(met.formula) if met.formula is not None else 0) for met in model.metabolites)
    gpr_max_length = max(len(rxn.gene_reaction_rule) for rxn in model.reactions)
    R = np.array(rxn_data, dtype=[
        ('id', f"<U{id_max_length}"),
        ('lb', 'float'),
        ('ub', 'float'),
        ('subsystem', f"<U{subsyst_max_length}"),
        ('gpr', f"<U{gpr_max_length}")
    ])
    M = np.array(met_data, dtype=[
        ('id', f"<U{metid_max_length}"),
        ('name', f"<U{metname_max_length}"),
        ('formula', f"<U{metform_max_length}")
    ])
    return MiomNetwork(S, R, M)


def _export_cobra_model(model_or_path, path_to_exported_file):
    m = load_gem(model_or_path)
    export_gem(m, path_to_exported_file)


def _load_compressed_model(url_or_filepath):
    file = url_or_filepath
    if _is_url(url_or_filepath):
        file = _download(url_or_filepath)
    ext = pathlib.Path(file).suffix
    if ext == '.xz' or ext == '.miom':
        import lzma
        from io import BytesIO
        with lzma.open(file, 'rb') as f_in:
            M = np.load(BytesIO(f_in.read()))
    else:
        M = np.load(file)
    return MiomNetwork(M['S'], M['reactions'], M['metabolites'])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Export cobrapy models')
    parser.add_argument('in_file', type=str, help='URL or input file (a .xml, .json, .yml, .mat model)')
    parser.add_argument('out_file', type=str, help='Output file (a npz.xz compressed file)')
    args = parser.parse_args()
    print("Exporting", args.in_file, "to", args.out_file)
    _export_cobra_model(args.in_file, args.out_file)
    print("Verifying...")
    m = _load_compressed_model(args.out_file)
    print(m.S.shape, m.R.shape, m.M.shape)
    print("Done")
