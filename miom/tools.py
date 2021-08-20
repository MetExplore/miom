import numpy as np
import picos as pc


def irev_blocked_lp(S, R):
    """
    Linear Program to detect all blocked irreversible reactions.
    Code adapted from the blocked.m method from Mojtaba Tefagh, Stephen P. Boyd, 2019 (Stanford University)
    Original code: https://github.com/mtefagh/swiftcore/blob/master/src/blocked.m

    :param S: Stoichiometric matrix
    :param R: Matrix of size Num. R x 3 with the reactions
    :return: tuple LP, V (LP = picos problem, V = optimized variables)
    """
    m, n = S.shape
    # Find irreversible reactions
    irev_rxn = np.array([True if rxn['lb'] >= 0 else False for rxn in R])
    irev = m + np.where(irev_rxn)[0]
    A = np.concatenate((S.T, -np.eye(n)), axis=1)
    lb = np.concatenate(([-np.inf] * m, [0] * n))
    lb[irev] = -1
    ub = np.concatenate(([np.inf] * m, [0] * n))
    V = pc.RealVariable("V", (n + m, 1), lower=lb, upper=ub)
    LP = pc.Problem()
    LP.add_constraint(A[irev_rxn, :] * V <= 0)
    LP.add_constraint(A[np.logical_not(irev_rxn), :] * V == 0)
    obj = np.zeros((m + n, 1))
    obj[irev] = 1
    C = pc.Constant("C", value=obj)
    LP.set_objective('min', C.T * V)
    return LP, V


def swiftcc(S, R, solver='glpk', verbosity=0):
    """
    SWIFTCC method for finding the largest flux consistent subnetwork of
    the original metabolic network. The implementation is based on the
    matlab code from Mojtaba Tefagh, Stephen P. Boyd, 2019 (Stanford University).
    https://github.com/mtefagh/swiftcore/blob/master/src/swiftcc.m

    :param S: Stoichiometric matrix
    :param R: Matrix of size Num. R x 3 with the reactions
    :param solver: picos-compatible linear solver (gplk, cplex, gurobi...)
    :param verbosity: 1 to show solver specific information, 0 to disable output
    :return: indexes of the consistent reactions
    """
    from scipy.linalg import qr
    m, n = S.shape
    rev = np.array([True if rxn['lb'] < 0 else False for rxn in R])
    LP, V = irev_blocked_lp(S, R)
    LP.options.verbosity = verbosity
    LP.options.solver = solver
    LP.solve()
    consistent = np.array([True] * n)
    consistent[np.array(V[m:].value).squeeze() < -0.5] = False
    # Estimate the rank of the consistent matrix using qr
    # (equivalent to np.linalg.rank_matrix)
    tol = np.linalg.norm(S[:, consistent]) * np.finfo(S.dtype).eps
    Q, R, _ = qr(S[:, consistent].T, pivoting=True)
    srank = np.sum(np.abs(np.diag(R)) > tol)
    Z = Q[rev[consistent] == True, srank:]
    consistent[np.logical_and(consistent, rev == True)] = np.diag(np.dot(Z, Z.T)) > tol ** 2
    return np.where(consistent == True)[0]

def consistent_subnetwork(network):
    """Finds the largest consistent subnetwork of the original network.

    Args:
        network (MiomNetwork): A MiomNetwork instance.

    Returns:
        MiomNetwork: Flux consistent subnetwork.
    """
    return network.subnet(swiftcc(network.S, network.R))
