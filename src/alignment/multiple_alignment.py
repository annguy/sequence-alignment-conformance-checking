"""multiple_alignment

This script is for multiple sequence alignment.
Expend the pairwise alignment for sequence, alignment aligning,
and alignment pair aligning. Using distance matrix to compute 
guid tree.

This file can also be imported as a module and contains the following
main functions:

    * global_realign - returns alignment for Needleman_Wunsch algorithm
    * align_sequences - returns alignment for pair of sequences
    * align_aligns - returns alignment for sequence and alignment or pair 
      of alignments
    * multi_sequence_alignment - returns multiple sequence alignments
    
"""


from functools import partial
import collections
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from itertools import combinations
import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('seaborn')
from more_itertools import locate


_traceback_encoding = {'match': 1, 'vertical-gap': 2, 'horizontal-gap': 3,
                       'uninitialized': -1, 'alignment-end': 0}
def iflatten(iterable):
    for elt in iterable:
        if isinstance(elt, collections.Iterable) and not isinstance(elt, str):
            for sub in flatten(elt):
                yield sub
        else:
            yield elt

def flatten(iterable):
    '''merge iterable'''
    return list(iflatten(iterable))

def merge(*sequences):
    '''merge alignments'''
    return list(map(flatten, zip(*sequences)))
 

def backtrace(traceback_matrix, sequence_a, sequence_b):
    """transcribing traceback_matrix

    Parameters
    ----------
    traceback_matrix : matrix
    sequence_a : iterable
        list or pair of activities
    sequence_b : iterable
        list or pair of activities
   
    Returns
    -------
    align1 : iterable
        aligned activities for the first input
    align2 : iterable
        aligned activities for the second input
    """
    i, j = len(sequence_a), len(sequence_b)
    align1, align2 = [], []
    fill_a, fill_b = '-', '-'
    if any(isinstance(e, (tuple, list)) for e in sequence_a):
        fill_a = ('-',) * len(sequence_a[0])
    if any(isinstance(e, (tuple, list)) for e in sequence_b):
        fill_b = ('-',) * len(sequence_b[0])
    while True:
        p = traceback_matrix[i, j]
        if p == _traceback_encoding['alignment-end']:
            break
        if p == _traceback_encoding['match']:
            align1.append(sequence_a[i - 1])
            align2.append(sequence_b[j - 1])
            i, j = i - 1, j - 1
        elif p == _traceback_encoding['horizontal-gap']:
            align1.append(fill_a)
            align2.append(sequence_b[j - 1])
            j -= 1
        elif p == _traceback_encoding['vertical-gap']:
            align1.append(sequence_a[i - 1])
            align2.append(fill_b)
            i -= 1
        else:
            raise ValueError("Something went terribly wrong.")
    align1 = align1[::-1]
    align2 = align2[::-1]
    return align1, align2


def global_realign(sequence_a, sequence_b, score_dict, gap_penalty=1, scale=1.0):
    """returns alignment for Needleman_Wunsch algorithm
    
    Parameters
    ----------
    sequence_a : iterable
        list or pair of activities
    sequence_b : iterable
        list or pair of activities
    score_dict : diction
        a dictionary holding the scores between all pairwise items 
        in sequence_a and sequence_b.
    gap_penalty : int
        the gap opening penalty used in the analysis.
    scale : float
        the factor by which gap_penalty should be decreased.
        
    Returns
    -------  
    align1 : iterable
        aligned activities for the first input
    align2 : iterable
        aligned activities for the second input
    distance : float    
        distance between those two sequences         
    """
    len1, len2 = len(sequence_a), len(sequence_b)
    traceback_matrix = np.zeros((len1 + 1, len2 + 1), dtype='i')
    score_matrix = np.zeros((len1 + 1, len2 + 1), dtype='f')
    length = np.zeros((len1 + 1, len2 + 1), dtype='f')
    traceback_matrix[0, 0] = _traceback_encoding['alignment-end']
    traceback_matrix[0, 1:] = _traceback_encoding['horizontal-gap']
    traceback_matrix[1:, 0] = _traceback_encoding['vertical-gap']
    for i in range(1, len1 + 1):
        score_matrix[i, 0] = score_matrix[i - 1, 0] + gap_penalty * scale
        length[i, 0] = length[i - 1, 0] + gap_penalty * scale
    for j in range(1, len2 + 1):
        score_matrix[0, j] = score_matrix[0, j - 1] + gap_penalty * scale
        length[0, j] = length[0, j - 1] + gap_penalty * scale
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            gap_a = score_matrix[i - 1, j] + (gap_penalty * scale if traceback_matrix[i - 1, j] == _traceback_encoding['vertical-gap'] else gap_penalty)
            gap_b = score_matrix[i, j - 1] + (gap_penalty * scale if traceback_matrix[i, j - 1] == _traceback_encoding['horizontal-gap'] else gap_penalty)
            match = score_matrix[i - 1, j - 1] + score_dict[i - 1, j - 1]
            if gap_a < match and gap_a <= gap_b:
                score_matrix[i, j] = gap_a
                traceback_matrix[i, j] = _traceback_encoding['vertical-gap']
            elif match <= gap_b:
                score_matrix[i, j] = match
                traceback_matrix[i, j] = _traceback_encoding['match']
            else:
                score_matrix[i, j] = gap_b
                traceback_matrix[i, j] = _traceback_encoding['horizontal-gap']
            p = traceback_matrix[i, j]
            l_gap_a = length[i - 1, j] + (gap_penalty * scale if p == _traceback_encoding['vertical-gap'] else 0)
            l_gap_b = length[i, j - 1] + (gap_penalty * scale if p == _traceback_encoding['horizontal-gap'] else 0)
            l_match = length[i - 1, j - 1] + (score_dict[i - 1, j - 1] if p == _traceback_encoding['match'] else 0)
            length[i, j] = max(l_gap_a, l_gap_b, l_match)
    # normalize the distance
    distance = score_matrix[len1, len2] / length[len1, len2]
    align1, align2 = backtrace(traceback_matrix, sequence_a, sequence_b)
    
    return align1, align2, distance
    


def align_sequences(sequence_a, sequence_b, scoring_fn=None, gap_penalty=1, scale=1.0):
  
    """Align two sequences using the Needleman-Wunsch algorithm
    
    Parameters
    ----------
    sequence_a : iterable
        list or pair of activities
    sequence_b : iterable
        list or pair of activities
    scoring_fn : function
        a distance function
    gap_penalty : int
        the gap opening penalty used in the analysis.
    scale : float
        the factor by which gap_penalty should be decreased.
        
    Returns
    -------  
    align1 : iterable
        aligned activities for the first input
    align2 : iterable
        aligned activities for the second input
    distance : float    
        distance between those two sequences
    """
    
    if scoring_fn is None:
        scoring_fn = lambda a, b: 0.0 if a == b else 2.0
    scores = {(i, j): scoring_fn(sequence_a[i], sequence_b[j])
              for i in range(len(sequence_a)) for j in range(len(sequence_b))}
              
    align1, align2, distance = global_realign(sequence_a, sequence_b, scores, gap_penalty, scale)          
              
    return align1, align2, distance

def align_aligns(sequence_a, sequence_b, scoring_fn=None,
                   gap_penalty=1, scale=1.0, gap_weight=0.5):
                   
    """Align two aligns or align and sequence using the Needleman-Wunsch algorithm
    
    Parameters
    ----------
    sequence_a : iterable
        list or pair of aligns
    sequence_b : iterable
        list or pair of aligns
    scoring_fn : function
        a distance function
    gap_penalty : int
        the gap opening penalty used in the analysis.
    scale : float
        the factor by which gap_penalty should be decreased.
    gap_weight : float
        weight gap between elements in aligns
        
    Returns
    -------  
    align1 : iterable
        aligned element for the first input
    align2 : iterable
        aligned element for the second input
    distance : float    
        distance between those two aligns
    """               
                             
    scores = {}
    for i in range(len(sequence_a)):
        for j in range(len(sequence_b)):
            dist = 0.0
            count = 0.0
            for k in range(len(sequence_a[i])):
                for l in range(len(sequence_b[j])):
                    mi = sequence_a[i][k]
                    mj = sequence_b[j][l]
                    if mi == '-' or mj == '-':
                        dist += gap_weight
                    else:
                        dist += 0.0 if scoring_fn(mi, mj) < 1 else 1.0
                    count += 1.0
                scores[i, j] = dist / count
    
    align1, align2, distance = global_realign (sequence_a, sequence_b, scores, gap_penalty, scale)
    return align1, align2, distance


def pairwise_distances(sequences, fn):
    """compute distance matrix
    
    Parameters
    ----------
    sequences : list
        list of sequences
    scoring_fn : function
        a distance function
    
    Returns
    -------  
    distances : matrix    
        distance matrix for a set of sequences
    """
    distances = np.zeros((len(sequences), len(sequences)))
    for i in range(len(sequences)):
        for j in range(i):
            _, _, distance = fn(sequences[i], sequences[j])
            distances[i, j] = distance
            distances[j, i] = distances[i, j]
    return distances

def multi_sequence_alignment(sequences, scoring_fn=None, linking='single',
                            gap_penalty=1, scale=1.0, gap_weight=1.0):
    """returns multiple sequence alignments
    
    Parameters
    ----------
    sequences : list
        list of sequences
    scoring_fn : function
        a distance function
    linking : str
        linking method
    scale : float
        the factor by which gap_penalty should be decreased.
    gap_penalty : int
        the gap opening penalty used in the analysis.
    scale : float
        the factor by which gap_penalty should be decreased.
    gap_weight : float
        weight gap between elements in aligns
        
    Returns
    -------  
    realignments : list
        aligned multiple sequences
    """                         
                            

    if scoring_fn is None:
        scoring_fn = lambda a, b: 0.0 if a == b else 2.0
    # compute all pairwise distances
    matrix = pairwise_distances(sequences, partial(align_sequences, scoring_fn=scoring_fn,
                                                   gap_penalty=gap_penalty, scale=scale))
    # compute the guiding tree to do the progressive alignment
    Z = linkage(squareform(matrix), method='single') 
    # perform the alignment by iterating through the clusters
    alignments = {}
    n_seqs = len(sequences)
    for cluster_id, (node1, node2, _, _) in enumerate(Z, n_seqs):
        node1, node2 = int(node1), int(node2)
        if node1 < n_seqs and node2 < n_seqs:
            align1, align2, _ = align_sequences(sequences[node1], sequences[node2],
                                                scoring_fn, gap_penalty, scale)
        else:
            if node1 < n_seqs:
                sequence_a, sequence_b = [[elt] for elt in sequences[node1]], alignments[node2]
            elif node2 < n_seqs:
                sequence_a, sequence_b = alignments[node1], [[elt] for elt in sequences[node2]]
            else:
                sequence_a, sequence_b = alignments[node1], alignments[node2]
            align1, align2, _ = align_aligns(sequence_a, sequence_b, scoring_fn, gap_penalty, scale, gap_weight)
        alignments[cluster_id] = merge(align1, align2)
        realignments = list(zip(*map(flatten, alignments[max(alignments)])))
        realignments = [list(i) for i in realignments]
        
    realignments = sort_realignments(realignments, sequences)
        
    return realignments


def sort_realignments(realignments, sequences):
    """returns multiple sequence alignments
    
    Parameters
    ----------
    realignments : list
        aligned multiple sequences
    sequences : list
        list of sequences with fixed order
    
    Returns
    -------  
    realignments_sort : list
        sorted sequences with the same oder as input
    """            

    ## find name for realigned seqs
    sequences = [''.join(acts) for acts in sequences]
    realignments = [''.join(acts) for acts in realignments]
    realigns = [i.replace('-', '') for i in realignments]
    indexes = []
    locate_list = [list(locate(sequences, lambda a: a == i)) for i in realigns] 
    for ps in locate_list:
        for p in ps:
            if p not in indexes:
                indexes.append(p)
                break
            else:
                continue
    realignments_sort = [x for _,x in sorted(zip(indexes, realignments))]
    return realignments_sort

def identity(realigns):
    """returns identities. All other seqs map to the first one
    
    Parameters
    ----------
    realigns : list
        aligned multiple sequences
    Returns
    -------  
    align_ident : list
        list of identities
    """ 
    def count(seq1,seq2):
        c = 0
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                c += 1 
        return c
    base = list(realigns[0])
    length = len(base)
    align_ident = [round(count(base,list(realigns[i]))/length,2) for i in range(1,len(realigns))]
    
    return align_ident