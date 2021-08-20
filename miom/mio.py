import io
import numpy as np
import pathlib
from urllib.parse import urlparse
from numbers import Number

# Default repository for loading miom models
DEFAULT_REPOSITORY = "https://github.com/pablormier/miom-gems/raw/main/gems/"


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


def set_repository(url):
    if not _is_url(url):
        raise ValueError("The default repository has to be a valid URL")
    global DEFAULT_REPOSITORY
    if not url.endswith('/'):
        url += '/'
    DEFAULT_REPOSITORY = url


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

    def subnet(self, idxs):
        S_sub = self.S[:, idxs]
        R_sub = self.R[idxs]
        act_met = np.sum(np.abs(S_sub), axis=1) > 0
        M_sub = self.M[act_met]
        S_sub = S_sub[act_met, :]
        return MiomNetwork(S_sub, R_sub, M_sub)

    @property
    def object_size(self):
        bytes = self.S.__sizeof__() + self.R.__sizeof__() + self.M.__sizeof__()
        return bytes / 1024**2


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
        model_or_path (str): Path to a local file or URL pointing to the metabolic network.
            If the string starts with '@', the file will be loaded from the default github
            repository.

    Returns:
        MiomNetwork: A [MiomNetwork][miom.mio] instance with the minimal information 
            of the metabolic network required for simulations. It includes the stoichiometric
            matrix, the list of reactions with the lower and upper bounds, the associated
            genes and GPR rules, and the list of metabolites.
    """
    extensions = ['.xml', '.yml', '.json', '.mat', '.miom']
    if isinstance(model_or_path, str):
        file = model_or_path
        if model_or_path.startswith("@"):
            # Check if the file has a valid file extension
            if not any(file.endswith(ext) for ext in extensions):
                # Assume is a relative url pointing to a version
                file = _download(DEFAULT_REPOSITORY + model_or_path[1:] + "/default.miom")
            else:
                file = _download(DEFAULT_REPOSITORY + model_or_path[1:])
        elif _is_url(model_or_path):
            file = _download(model_or_path)
        ext = pathlib.Path(file).suffix
        if ext == '.miom' or ext == '.xz' or ext == '.npz':
            return _load_compressed_model(file)
        else:
            return cobra_to_miom(_read_cobra_model(file))
    else:
        return cobra_to_miom(model_or_path)


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


def cobra_to_miom(model):
    try:
        from cobra.util.array import create_stoichiometric_matrix
    except ImportError as e:
        raise ImportError("Cobrapy package is not installed, "
                          "but required to read and import metabolic networks", e)
    S = create_stoichiometric_matrix(model)
    subsystems = []
    for rxn in model.reactions:
        subsys = rxn.subsystem
        list_subsystem_rxn = []
        # For .mat models, the subsystem can be loaded as a string repr of a numpy array
        if isinstance(subsys, str) and (subsys.startswith("array(") or subsys.startswith("[array(")):
            from numpy import array
            try:
                subsys = eval(subsys.strip())
            except:
                # Try to create a list
                import re
                subsys = re.findall('\[\'(.*?)\'\]', subsys)
                if len(subsys) == 0:
                    subsys = rxn.subsystem
            # A list containing a numpy array?
            for s in subsys:
                if "tolist" in dir(s):
                    list_subsystem_rxn.extend(s.tolist())
                else:
                    list_subsystem_rxn.append(s)
            if len(list_subsystem_rxn) == 1:
                list_subsystem_rxn = list_subsystem_rxn[0]
            subsystems.append(list_subsystem_rxn)
        elif "tolist" in dir(rxn.subsystem):
            subsystems.append(rxn.subsystem.tolist())
        else:
            subsystems.append(rxn.subsystem)
    rxn_data = [(rxn.id, rxn.name, rxn.lower_bound, rxn.upper_bound, subsystem, rxn.gene_reaction_rule)
                for rxn, subsystem in zip(model.reactions, subsystems)]
    met_data = [(met.id, met.name, met.formula) for met in model.metabolites]
    R = np.array(rxn_data, dtype=[
        ('id', 'object'),
        ('name', 'object'),
        ('lb', 'float'),
        ('ub', 'float'),
        ('subsystem', 'object'),
        ('gpr', 'object')
    ])
    M = np.array(met_data, dtype=[
        ('id', 'object'),
        ('name', 'object'),
        ('formula', 'object')
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
            M = np.load(BytesIO(f_in.read()), allow_pickle=True)
    else:
        M = np.load(file)
    return MiomNetwork(M['S'], M['reactions'], M['metabolites'])

