"""
============
ICE with CNV
============
"""

import argparse
import sys
import numpy as np
from scipy import sparse
import pandas as pd
import iced
from iced import io


def load_lengths(filename, return_base=False):
    """
    Fast loading of the bed files

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be a bed file

    Returns
    -------
    lengths : the lengths of each chromosomes
    """
    data = pd.read_csv(filename, sep="\t", comment="#", header=None)
    data = data.as_matrix()
    _, idx, lengths = np.unique(data[:, 0], return_counts=True,
                                return_index=True)
    if return_base:
        return lengths[idx.argsort()], data[0, 3]
    else:
        return lengths[idx.argsort()]


def load_counts(filename, lengths=None, base=None):
    """
    Fast loading of a raw interaction counts file

    Parameters
    ----------
    filename : str
        path to the file to load. The file should be of the following format:
        i, j, counts

    lengths : ndarray
        lengths of each chromosomes

    Returns
    --------
    X : the interaction counts file
    """
    n = None
    if lengths is not None:
        n = lengths.sum()
        shape = (n, n)
    else:
        shape = None
    # This is the interaction count files
    dataframe = pd.read_csv(filename, sep="\t", comment="#", header=None)
    row, col, data = dataframe.as_matrix().T

    # If there are NAs remove them
    mask = np.isnan(data)
    if np.any(mask):
        warnings.warn(
            "NAs detected in %s. "
            "Removing NAs and replacing with 0." % filename)
        row = row[np.invert(mask)]
        col = col[np.invert(mask)]
        data = data[np.invert(mask)]

    # XXX We need to deal with the fact that we should not duplicate entries
    # for the diagonal.
    # XXX what if n doesn't exist?
    if base is not None:
        if base not in [0, 1]:
            raise ValueError("indices should start either at 0 or 1")
        col -= base
        row -= base
    else:
        warnings.warn(
            "Attempting to guess whether counts are 0 or 1 based")

        if (col.min() >= 1 and row.min() >= 1) and \
           ((n is None) or (col.max() == n)):
            # This is a hack to deal with the fact that sometimes, the files
            # are indexed at 1 and not 0

            col -= 1
            row -= 1

    if shape is None:
        n = max(col.max(), row.max()) + 1
        shape = (n, n)

    data = data.astype(float)
    counts = sparse.coo_matrix((data, (row, col)), shape=shape)
    return counts


parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--remove-all-zeros-loci", default=False,
                    action="store_true")
parser.add_argument("--filter_low_counts_perc", default=0.02,
                    type=float)
parser.add_argument("--classic", action="store_true", default=False)
parser.add_argument("--bed-file", "-b")
parser.add_argument("--seg-file", "-s")
parser.add_argument("--max_iter", "-m", default=100, type=int,
                    help="Maximum number of iterations")
parser.add_argument("--remove-cnv", "-r", default=False, action="store_true")
parser.add_argument("--outfile", "-o", default=None)
parser.add_argument("--type", "-t", default="independant", type=str,
                    help="Can be 'independant', 'joint', or 'chrjoint'")
parser.add_argument("--dense", default=False, action="store_true")
args = parser.parse_args()

print("Loading files...")
if args.seg_file is not None:
    lengths = load_lengths(args.seg_file)
    base = 1
elif args.bed_file is not None:
    lengths, base = load_lengths(args.bed_file, return_base=True)
else:
    lengths = load_lengths(args.filename.replace(".matrix", ".bed"),
                              return_base=True)

counts = load_counts(args.filename, lengths=lengths, base=base)

if args.type not in ["independant", "joint", "chrjoint"]:
    raise ValueError("Unknown type of normalization")

if args.dense:
    counts = counts.toarray()
    counts = counts.T + counts

# Filter counts
counts = iced.filter.filter_low_counts(
    counts,
    percentage=args.filter_low_counts_perc,
    remove_all_zeros_loci=args.remove_all_zeros_loci,
    copy=False, sparsity=False,
    verbose=100)

# Remove diagonal
if sparse.issparse(counts):
    counts.setdiag(0)
    counts = counts.tocsr()
    counts.eliminate_zeros()
    #counts = counts.tocoo()
else:
    counts[np.diag_indices_from(counts)] = 0


if not args.classic:
    if args.seg_file is not None:
        cnv = np.loadtxt(args.seg_file, usecols=(3, ))
    else:
        cnv = np.loadtxt(args.bed_file, usecols=(3, ))
    breakpoints = np.where((cnv[1:] - cnv[:-1]).astype(bool))[0] + 1
    breakpoints = np.array(list(breakpoints) + list(lengths.cumsum()))
    breakpoints = np.unique(breakpoints)
    breakpoints.sort()
    print(
        "Found %d breakpoints, with %d CNV profiles" % (
            len(breakpoints), len(np.unique(cnv))))
    print("Estimating CNV bias from data...")
    begin, end = 0, 0
    filtered = (np.array(counts.sum(axis=0)).flatten() +
                np.array(counts.sum(axis=1)).flatten()) == 0
    bed = cnv.copy()
    for c in np.unique(bed):
        mask = bed == c
        length = mask.sum() - filtered[mask].sum()
        if length == 0:
            cnv[mask] = 0
            begin = end
            continue
        diag = counts[mask][:, mask]
        cnv[mask] = (counts[mask].sum() + counts[:, mask].sum() -
                     diag.sum())
        cnv[mask] /= length
        begin = end

        # Remove cnv which are 0. They correspond to filtered out regions
        cnv[cnv == 0] = 1
else:
    cnv = None

counts = iced.normalization.ICE_normalization(
    counts, eps=0.1, copy=False,
    verbose=100, max_iter=args.max_iter, counts_profile=cnv)

if args.outfile is not None:
    outfile = args.outfile
elif args.remove_cnv:
    outfile = args.filename + "_removed_cnv_bias.matrix"
elif args.classic:
    outfile = args.filename + "_iced.matrix"
else:
    outfile = args.filename + "cnv_iced.matrix"

if not args.remove_cnv or args.classic:
    print("Saving...")
    counts.row += 1
    counts.col += 1
    io.write_counts(outfile, counts)
    sys.exit(0)

# Remove CNV bias
cnv = bed
print("Converting to dense")
counts = counts.toarray()
counts = counts.T + counts
sum_c = counts.sum(axis=0) == 0
counts[sum_c] = np.nan
counts[:, sum_c] = np.nan

print("")
print("Estimating CNV-effects.")
max_iter = 20
min_iter = 5
bias_es = iced.normalization._ca_utils.estimate_block_biases(counts, lengths,
                                                             cnv, verbose=True)

if args.outfile is not None:
    outfile = args.outfile
else:
    outfile = args.filename + "_removed_cnv_bias.matrix"

if sparse.issparse(counts):
    counts.data /= bias_es
else:
    counts /= bias_es
    counts[np.isnan(counts)] = 0
    counts = sparse.coo_matrix(np.triu(counts))
    bias_es = bias_es[counts.row, counts.col]

counts.row += 1
counts.col += 1

print("Saving...")
io.write_counts(outfile, counts)
counts.data = bias_es

if outfile.endswith(".matrix"):
    outfile = outfile.replace(".matrix", "_cnv_effect.matrix")
else:
    outfile = outfile + "_cnv_effect"

io.write_counts(outfile, counts)
