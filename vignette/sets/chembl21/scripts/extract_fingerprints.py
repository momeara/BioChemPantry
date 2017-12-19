#!/env/python

"""
For a given sea library, extract bit vectors

"""

import sys
from libcore.zodb_library.library import Library
from seashell.cli.fputil import FingerprintType
from fpcore.fingerprint import (
    get_fingerprint_generator,
)
import numpy as np


library_fname = sys.argv[1]

library = Library(library_fname)
fp_generator = get_fingerprint_generator(FingerprintType(library.sea_meta).name)
for target_key, target_set in library.sets.iteritems():
    fps = []
    cids = []
    for cid in target_set.cids:
        fp = library.raw_fingerprints[cid]
        fps.append((fp,cid))
        cids.append(cid)
    import pdb
    pdb.set_trace()

    fp_array = fp_generator.to_array(fps)

    # remove invariant columns
    fp_variate_mask = ~ np.all(fp_array == fp_array[0,:], axis=0)
    fp_array = fp_array[:, fp_variate_mask]
    N, p = fp_array.shape

    # standarize
    fp_array = (fp_array - np.mean(fp_array, axis=0)) / np.std(fp_array, axis=0)

    # correlation matrix
    # note the correlation matrix is the covariance matrix is of a standardized data array
    fp_cor = np.cov(fp_array.T)

    evals, evecs = np.linalg.eigh(fp_cor)
    
    projections = np.dot(fp_array, evecs)
    
    MP_gamma = p/N
    MP_bound = np.square(1+np.sqrt(MP_gamma))

    exemplars  = cids[np.argsort(projections, axis=0)]

    

    import pdb
    pdb.set_trace()
    pass

                         
library.close()
