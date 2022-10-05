#!/usr/bin/env python
# -*- tab-width:4;indent-tabs-mode:f;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=4 et sw=4:

import csv
import numpy as np
import sys
import os
from scipy.sparse.linalg import eigsh
import traceback

from fpcore.fingerprint import get_fingerprint_generator
from libcore.zodb_library.library import Library
from libcore.util.util import get_library_cutoff
from seacore.run.scoring import tc_matrix
from seashell.cli.fputil import FingerprintType


"""
Collect target set correlation statistics for a SEA library

Usage:

   python compute_target_set_correlation.py <library>.sea <tc_cutoff> [library,background] <output_statistics>.csv

"""



STATISTIC_FIELD_NAMES=[
	'ref_set_id',
	'ref_set_size',
	'query_set_id',
	'query_set_size',

	'cut_set_sum',
	'gini_coeffient',
	'largest_eigenvalue',
	'smallest_eigenvalue',

#	'ref_n_variable_bits',
#	'ref_cor_evals',
#	'ref_mp_max_prob',
#
#	'query_n_variable_bits',
#	'query_cor_evals',
#	'query_mp_max_prob',
]

def compute_gini_coefficient(y):
    """
    from Gini in the ineq package in R with the finite sample correction:
    G <- (2 * sum(x * 1L:n) / sum(x) - (n+1))/(n-1)
    """
    x = np.sort(y, axis=None)
    n = x.shape[0]
    z = range(1,n+1) / x.sum()
    return (2 * np.dot(x, z) - (n + 1)) / (n - 1)


def correlation_matrix_eigen_values(
	fp_type,
	target_fingerprints):

    #get the fingerprints as an N by p array
    fp_generator = get_fingerprint_generator(fp_type.name)
    fp_array = fp_generator.to_array(target_fingerprints)

    # remove invariant columns
    fp_variate_mask = ~ np.all(fp_array == fp_array[0,:], axis=0)
    fp_array = fp_array[:, fp_variate_mask]
    N, p = fp_array.shape
    # N number of substances in set
    # p number of feature columns that are not invariant

    # compute eigenspectrum of the correlation matrix
    fp_array = (fp_array - np.mean(fp_array, axis=0)) / np.std(fp_array, axis=0)
    fp_cor = np.cov(fp_array.T)
    evals, evecs = np.linalg.eigh(fp_cor)
	return N, p, evals, evecs

def mp_distribution_prob(
	eigen_value,
	p,      # length of fingerprint (only non-invariant bits)
	N):     # N ligands

	gamma = float(p)/float(N)

	try:
		discriminant = \
            ((1.0 + np.sqrt(gamma))**2 - eigen_value) * \
		    (eigen_value - (1.0-np.sqrt(gamma))**2)
		if discriminant < 0:
			return float("-infinity")
		else:
			prob = np.log(np.sqrt(discriminant) / (2.0 * np.pi * gamma * eigen_value))
	except:
		prob = float("-infinity")
	return prob


def get_target_comparison_statistics(
	fp_type,
    tc_cutoff,
    ref_key,
    ref_fps,
    query_key,
    query_fps):

    ref_set_size = len(ref_fps)
    query_set_size = len(query_fps)

    tcs = tc_matrix(ref_fps, query_fps, fp_type.mode)
    tcs = np.asarray(tcs).reshape(ref_set_size, query_set_size)
    tcs[ tcs < tc_cutoff ] = 0

    cut_set_sum = tcs.sum()
    gini_coefficient = compute_gini_coefficient(tcs)
    try:
        largest_eigenvalue = eigsh(tcs, 1, which='LM')[0][0]
    except:
        largest_eigenvalue = None

    try:
        smallest_eigenvalue = eigsh(tcs, 1, sigma=0, which='LM')[0][0]
    except:
        smallest_eigenvalue = None

#	try:
#		ref_N, ref_p, ref_cor_evals, ref_cor_evecs = correlation_matrix_eigen_values(fp_type, ref_fps)
#		ref_mp_max_prob = max([mp_distribution_prob(e, ref_N, ref_p) for e in ref_cor_evals])
#		ref_cor_evals = ";".join(map(str, ref_cor_evals))
#	except:
#		ref_p = None
#		ref_cor_evals = None
#		ref_mp_max_prob = None
#
#	if query_key == ref_key:
#		query_p = ref_p
#		query_cor_evals = ref_cor_evals
#		query_mp_max_prob = ref_mp_max_prob
#	else:
#    	try:
#    		query_N, query_p, query_cor_evals, query_cor_evecs = \
#                correlation_matrix_eigen_values(fp_type, query_fps)
#    		query_mp_max_prob = \
#                max([mp_distribution_max_prob(e, query_n, query_p) for e in query_cor_evals])
#    		query_cor_evals = ";".join(map(str, query_cor_evals))
#    	except:
#    		query_p = None
#    		query_cor_evals = None
#    		query_mp_max_prob = None

    statistics = [
        ref_key,
        ref_set_size,
        query_key,
        query_set_size,
        cut_set_sum,
        gini_coefficient,
        largest_eigenvalue,
        smallest_eigenvalue,
#		ref_p,
#		ref_cor_evals,
#		ref_mp_max_prob,
#		query_p,
#		query_cor_evals,
#		query_mp_max_prob,
    ]
    return statistics

def library_target_statistics(
    library,
    tc_cutoff):

	fp_type = FingerprintType(library.sea_meta)

    statistics = []
    for i, target_key in enumerate(library.sets.keys()):
        target_fps = [
            (library.raw_fingerprints[cid], cid) for cid in library.sets[target_key].cids]
        target_statistics = get_target_comparison_statistics(
			fp_type=fp_type,
            tc_cutoff=tc_cutoff,
            ref_key=target_key.tid,
            ref_fps=target_fps,
            query_key=target_key.tid,
            query_fps=target_fps)

        statistics.append(target_statistics)

    return statistics


def background_target_statistics(
    library,
    tc_cutoff,
    set_size_min,
    set_size_max,
    n_samples):

    fp_format = FingerprintType(library.sea_meta).mode
    setA, n_setA = seacore.background.fingerprint_set(library)
    setB, n_setB = setA, n_setB

    bins = seacore.background.init()
    lower = bins.min**2
    upper = bins.max**2
    pop_size = n_setA * n_setB
    if upper > pop_size:
        upper  = pop_size
        set_limit = min(n_setA, n_setB)
        bins = BinInfo(bins.min, set_limit, bins.count, bins.samples)
    step = math.log(upper - lower , 2)/bins.count
    shuffled = range(bins.count)
    random.shuffle(shuffled)
    finished = set()
    for i in shuffled:
        target = int(2**(i*step)) # log space
        if target in finished:
            # skip integer collisions in low log-space
            continue
        finished.add(target)
        sample_sets(target, bins, setA, n_setA, setB, n_setB, scoring)


    statistics = []
    for set_size in np.random.random_integers(set_size_min, set_size_max, n_samples):
        target_fps = np.random.choice(fps, set_size)
        target_statistics = get_target_comparison_statistics(
            fp_format=fp_format,
            tc_cutoff=tc_cutoff,
            ref_key=target_key.tid,
            ref_fps=target_fps,
            query_key=target_key.tid,
            query_fps=target_fps)
        statistics.append(target_statistics)
    return statistics

def write_statistics(
    statistics, fname):

    with open(fname, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(STATISTIC_FIELD_NAMES)
        [writer.writerow(s) for s in statistics]


def main():
    library_fname = sys.argv[1]
    tc_cutoff = float(sys.argv[2])
    library_or_background = sys.argv[3]
    out_fname = sys.argv[4]

    library = Library(library_fname)
    try:
        if library_or_background == 'library':
            s = library_target_statistics(library, tc_cutoff)
        elif library_or_background == 'background':
            s = background_target_statistics(
                library=library,
                tc_cutoff=tc_cutoff,
                set_size_min=1,
                set_size_max=300,
                n_samples=4000)
        else:
            print "The third argument should be either 'library' or background"
            return 1

        write_statistics(s, out_fname)
    except Exception as e:
        print e.message
        print traceback.format_exc()
    finally:
        if library:
            library.close()

if __name__ == '__main__':
    main()





