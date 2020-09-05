"""pairwise_alignment

This script is a rewrite of source code of the scikit-bio and provides
algorithms for computing pairwise alignment and distances. 

Don't need the substitution_matrix. The type of input sequence change
from Protein or DNA str to activities/event log

This file can also be imported as a module and contains the following
main functions:

    * Needleman_Wunsch - returns alignment for Needleman_Wunsch algorithm
    * Mith_Waterman - returns alignment for Mith_Waterman algorithm
    * LCS - returns alignment for longest common subsequence algorithm
    * Levenshtein - returns Levenshtein distance between two sequences
    * Damerau_Levenshtein - returns Damerau_Levenshtein distance between 
    two sequences
    
"""


#rewrite the source code of the scikit-bio. Don't need the substitution_matrix
#and 
import numpy as np


_traceback_encoding = {'match': 1, 'vertical-gap': 2, 'horizontal-gap': 3,
                       'uninitialized': -1, 'alignment-end': 0}
					   
def _init_matrices_nw(aln1, aln2, gap_open_penalty, gap_extend_penalty):
    """initialize score matrix and traceback matrix for global alignment

    Parameters
    ----------
    aln1 : list
        list of activities, which is the first sequence to be aligned
    aln2 : list
        list of activities, which is the second sequence to be aligned
    gap_open_penalty : int
    gap_extend_penalty : int
    
    Returns
    -------
    score_matrix: matrix
    traceback_matrix: matrix
    """

    shape = (len(aln2)+1, len(aln1)+1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=np.int)
    traceback_matrix += _traceback_encoding['uninitialized']
    traceback_matrix[0, 0] = _traceback_encoding['alignment-end']

    # cache some values for quicker access
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']

    for i in range(1, shape[0]):
        score_matrix[i, 0] = +gap_open_penalty + ((i-1) * gap_extend_penalty)
        traceback_matrix[i, 0] = vgap

    for i in range(1, shape[1]):
        score_matrix[0, i] = +gap_open_penalty + ((i-1) * gap_extend_penalty)
        traceback_matrix[0, i] = hgap

    return score_matrix, traceback_matrix	

def _init_matrices_sw(aln1, aln2, gap_open_penalty, gap_extend_penalty):

    """initialize score matrix and traceback matrix for local alignment

    Parameters
    ----------
    aln1 : list
        list of activities, which is the first sequence to be aligned
    aln2 : list
        list of activities, which is the second sequence to be aligned
    gap_open_penalty : int
    gap_extend_penalty : int
    
    Returns
    -------
    score_matrix : matrix
    traceback_matrix : matrix
    """
    shape = (len(aln2)+1, len(aln1)+1)
    score_matrix = np.zeros(shape)
    traceback_matrix = np.zeros(shape, dtype=np.int)
    traceback_matrix += _traceback_encoding['uninitialized']
    traceback_matrix[0, :] = _traceback_encoding['alignment-end']
    traceback_matrix[:, 0] = _traceback_encoding['alignment-end']
    return score_matrix, traceback_matrix


	
def _first_largest(scores):
    """ Similar to max, but returns the first element achieving the high score
    If max receives a tuple, it will break a tie for the highest value
    of entry[i] with entry[i+1]. We don't want that here - to better match
    with the results of other tools, we want to be able to define which
    entry is returned in the case of a tie.
    
    """
    result = scores[0]
    for score, direction in scores[1:]:
        if score > result[0]:
            result = (score, direction)
    return result
	
def _compute_score_and_traceback_matrices(
        aln1, 
		aln2, 
		gap_open_penalty=-2, 
		gap_extend_penalty=-2, 
		match_reward=4,
		mismatch_penalty=-2,
        new_alignment_score=0,# for global: -np.inf ; for local : 0.0, 
		init_matrices_f=_init_matrices_sw,
        penalize_terminal_gaps=True):
        
    """Return dynamic programming (score) and traceback matrices.

    A note on the ``penalize_terminal_gaps`` parameter. 
    
    Parameters
    ----------
    aln1 : list
        list of activities, which is the first sequence to be aligned
    aln2 : list
        list of activities, which is the second sequence to be aligned
    gap_open_penalty : int
        the cost for starting a gap
    gap_extend_penalty : int
        the penalty for extending a gap by one residue
    match_reward : int 
        a reward if a pair of elements match
    mismatch_penalty : int 
        a penalty value if a pair of elements doesn't match
    new_alignment_score : const number
        for local alignment is 0 and for global is -inf
    init_matrices_f : function
        function for initialize score_matrix and traceback_matrix
    penalize_terminal_gaps : bool
        When this value is ``False``, this function is no longer truE
        Smith-Waterman/Needleman-Wunsch scoring, but when ``True`` it 
        can result in biologically irrelevant artifacts in Needleman-Wunsch 
        (global) alignments. Specifically, if one sequence is longer than 
        the other (e.g., if aligning a primer sequence to an amplification product,
        or searching for a gene in a genome) the shorter sequence will have 
        a long gap inserted.
        
    Returns
    -------
    score_matrix : matrix
    traceback_matrix : matrix
  
    """
    aln1_length = len(aln1)
    aln2_length = len(aln2)
    # cache some values for quicker/simpler access
    aend = _traceback_encoding['alignment-end']
    match = _traceback_encoding['match']
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']

    new_alignment_score = (new_alignment_score, aend)

    # Initialize a matrix to use for scoring the alignment and for tracing
    # back the best alignment
    score_matrix, traceback_matrix = init_matrices_f(
        aln1, aln2, gap_open_penalty, gap_extend_penalty)

    # Iterate over the characters in aln2 (which corresponds to the vertical
    # sequence in the matrix)
    for aln2_pos, aln2_chars in enumerate(aln2, 1):
        aln2_chars = str(aln2_chars)

        # Iterate over the characters in aln1 (which corresponds to the
        # horizontal sequence in the matrix)
        for aln1_pos, aln1_chars in enumerate(aln1, 1):
            aln1_chars = str(aln1_chars)
            substitution_score = 0
            # compute the score for a match/mismatch
            if aln1_chars==aln2_chars:
                substitution_score = match_reward
            else:
                substitution_score = mismatch_penalty
                #print(mismatch_penalty)

            diag_score = \
                (score_matrix[aln2_pos-1, aln1_pos-1] + substitution_score,
                 match)

            # compute the score for adding a gap in aln2 (vertical)
            if not penalize_terminal_gaps and (aln1_pos == aln1_length):
                # we've reached the end of aln1, so adding vertical gaps
                # (which become gaps in aln1) should no longer
                # be penalized (if penalize_terminal_gaps == False)
                up_score = (score_matrix[aln2_pos-1, aln1_pos], vgap)
            elif traceback_matrix[aln2_pos-1, aln1_pos] == vgap:
                # gap extend, because the cell above was also a gap
                up_score = \
                    (score_matrix[aln2_pos-1, aln1_pos] + gap_extend_penalty,
                     vgap)
            else:
                # gap open, because the cell above was not a gap
                up_score = \
                    (score_matrix[aln2_pos-1, aln1_pos] + gap_open_penalty,
                     vgap)

            # compute the score for adding a gap in aln1 (horizontal)
            if not penalize_terminal_gaps and (aln2_pos == aln2_length):
                # we've reached the end of aln2, so adding horizontal gaps
                # (which become gaps in aln2) should no longer
                # be penalized (if penalize_terminal_gaps == False)
                left_score = (score_matrix[aln2_pos, aln1_pos-1], hgap)
            elif traceback_matrix[aln2_pos, aln1_pos-1] == hgap:
                # gap extend, because the cell to the left was also a gap
                left_score = \
                    (score_matrix[aln2_pos, aln1_pos-1] + gap_extend_penalty,
                     hgap)
            else:
                # gap open, because the cell to the left was not a gap
                left_score = \
                    (score_matrix[aln2_pos, aln1_pos-1] + gap_open_penalty,
                     hgap)

            # identify the largest score, and use that information to populate
            # the score and traceback matrices
            best_score = _first_largest([new_alignment_score, left_score,
                                         diag_score, up_score])
            score_matrix[aln2_pos, aln1_pos] = best_score[0]
            traceback_matrix[aln2_pos, aln1_pos] = best_score[1]
    #print(score_matrix)
    return score_matrix, traceback_matrix
	
	
	
	
#Outputs are aligned_seqs1, aligned_seqs2, best_score,	
	
def _traceback(traceback_matrix, score_matrix, aln1, aln2, start_row, start_col):

    """return alignments by transcribing traceback_matrix and score_matrix

    Parameters
    ----------
    traceback_matrix : matrix
        get from _compute_score_and_traceback_matrices
    score_matrix : matrix
        get from _compute_score_and_traceback_matrices
    aln1 : list
        list of activities, which is the first sequence to be aligned
    aln2 : list
        list of activities, which is the second sequence to be aligned
    start_row : int
        row index for triscribing
    start_col : int
        column index for triscribing
    Returns
    -------
    aligned_seqs1 : list
        aligned activities for the first input
    aligned_seqs2 : list
        aligned activities for the second input
    best_score : float
        similarity score between two sequence
    identity : float
        identity between two sequence
    """

    # cache some values for simpler reference
    
    aend = _traceback_encoding['alignment-end']
    match = _traceback_encoding['match']
    vgap = _traceback_encoding['vertical-gap']
    hgap = _traceback_encoding['horizontal-gap']
    gap_character = '-'

    # initialize the result alignments
    #aln1_sequence_count = 1
    #aligned_seqs1 = [[] for e in range(aln1_sequence_count)]
    aligned_seqs1 = []
    #aln2_sequence_count = 1
    #aligned_seqs2 = [[] for e in range(aln2_sequence_count)]
    aligned_seqs2 = []	
    current_row = start_row
    current_col = start_col
    
    #print(traceback_matrix.shape)
    best_score = score_matrix[current_row, current_col]
    current_value = None
    identity = 0
	
    while current_value != aend:
        current_value = traceback_matrix[current_row, current_col]

        if current_value == match:
            
            #for aligned_seq, input_seq in zip(aligned_seqs1, aln1):
            #    aligned_seq.append(str(input_seq[current_col-1]))
   
            aligned_seqs1.append(aln1[current_col-1])
            #for aligned_seq, input_seq in zip(aligned_seqs2, aln2):
            #    aligned_seq.append(str(input_seq[current_row-1]))
            aligned_seqs2.append(aln2[current_row-1])
            
            if aln1[current_col-1] == aln2[current_row-1]:
                identity += 1
            current_row -= 1
            current_col -= 1
			
        elif current_value == vgap:
            #for aligned_seq in aligned_seqs1:
                #aligned_seq.append(gap_character)
            aligned_seqs1.append(gap_character)
            #for aligned_seq, input_seq in zip(aligned_seqs2, aln2):
                #aligned_seq.append(str(input_seq[current_row-1]))
            aligned_seqs2.append(aln2[current_row-1])
            
            current_row -= 1
			
        elif current_value == hgap:
            #for aligned_seq, input_seq in zip(aligned_seqs1, aln1):
                #print(input_seq)
                #aligned_seq.append(str(input_seq[current_col-1]))
            aligned_seqs1.append(aln1[current_col-1])
            #for aligned_seq in aligned_seqs2:
                #aligned_seq.append(gap_character)
            aligned_seqs2.append(gap_character)
            current_col -= 1
        elif current_value == aend:
            continue
        else:
            raise ValueError(
                "Invalid value in traceback matrix: %s" % current_value)
    aligned_seqs1.reverse()
    aligned_seqs2.reverse()
    if identity:
	
        identity =round( float(identity) / len(aligned_seqs1),2)
    else:
        identity = 0
    return aligned_seqs1, aligned_seqs2, best_score, identity
	

def Needleman_Wunsch(seq1, seq2, gap_open_penalty=-2, gap_extend_penalty=-2, match_reward=1, mismatch_penalty=-2):
	
    """return Needleman_Wunsch alignment

    Parameters
    ----------
    seq1 : list
        list of activities, which is the first sequence to be aligned
    seq2 : list
        list of activities, which is the second sequence to be aligned
    gap_open_penalty : int
        the cost for starting a gap
    gap_extend_penalty : int
        the penalty for extending a gap by one residue
    match_reward : int 
        a reward if a pair of elements match
    mismatch_penalty : int 
        a penalty value if a pair of elements doesn't match
    Returns
    -------
    aligned_seqs1 : list
        aligned activities for the first input
    aligned_seqs2 : list
        aligned activities for the second input
    best_score : float
        similarity score between two sequence
    identity : float
        identity between two sequence
    """					 
    score_matrix, traceback_matrix = _compute_score_and_traceback_matrices(
        seq1, 
        seq2, 
        gap_open_penalty,                         
        gap_extend_penalty,match_reward,mismatch_penalty,
        new_alignment_score=-np.inf,# for global: -np.inf ; for local : 0.0, 
        init_matrices_f=_init_matrices_nw,
        penalize_terminal_gaps=True
        )
		
    current_row = traceback_matrix.shape[0] - 1
    current_col = traceback_matrix.shape[1] - 1   
    aligned_seqs1, aligned_seqs2, best_score, identity =_traceback(traceback_matrix, score_matrix, seq1, seq2, current_row, current_col)
	
    return aligned_seqs1, aligned_seqs2, best_score, identity
	
def Mith_Waterman(seq1, seq2, gap_open_penalty=-2, gap_extend_penalty=-2, match_reward=4, mismatch_penalty=-2):
    """return Mith_Waterman alignment

    Parameters
    ----------
    seq1 : list
        list of activities, which is the first sequence to be aligned
    seq2 : list
        list of activities, which is the second sequence to be aligned
    gap_open_penalty : int
        the cost for starting a gap
    gap_extend_penalty : int
        the penalty for extending a gap by one residue
    match_reward : int 
        a reward if a pair of elements match
    mismatch_penalty : int 
        a penalty value if a pair of elements doesn't match
    Returns
    -------
    aligned_seqs1 : list
        aligned activities for the first input
    aligned_seqs2 : list
        aligned activities for the second input
    best_score : float
        similarity score between two sequence
    identity : float
        identity between two sequence
    """		
    score_matrix, traceback_matrix = _compute_score_and_traceback_matrices(
        seq1, 
        seq2, 
        gap_open_penalty,gap_extend_penalty,match_reward,mismatch_penalty,
        new_alignment_score=0.0,# for global: -np.inf ; for local : 0.0, 
        init_matrices_f=_init_matrices_sw,
        penalize_terminal_gaps=True
        )
    end_row_position, end_col_position = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
      
    aligned_seqs1, aligned_seqs2, best_score, identity =_traceback(traceback_matrix, score_matrix, seq1, seq2, end_row_position, end_col_position)
	
    return aligned_seqs1, aligned_seqs2, best_score, identity

def LCS(seq1, seq2, gap_open_penalty=0, gap_extend_penalty=0, match_reward=1, mismatch_penalty=-np.inf):
    """return Longest common subsequence alignment

    Parameters
    ----------
    seq1 : list
        list of activities, which is the first sequence to be aligned
    seq2 : list
        list of activities, which is the second sequence to be aligned
    gap_open_penalty : int
        the cost for starting a gap
    gap_extend_penalty : int
        the penalty for extending a gap by one residue
    match_reward : int 
        a reward if a pair of elements match
    mismatch_penalty : int 
        a penalty value if a pair of elements doesn't match
    Returns
    -------
    aligned_seqs1 : list
        aligned activities for the first input
    aligned_seqs2 : list
        aligned activities for the second input
    best_score : float
        similarity score between two sequence
    identity : float
        identity between two sequence
    """		    
    score_matrix, traceback_matrix = _compute_score_and_traceback_matrices(
        seq1, 
        seq2, 
        gap_open_penalty,gap_extend_penalty,match_reward,mismatch_penalty,
        new_alignment_score=0.0,# for global: -np.inf ; for local : 0.0, 
        init_matrices_f=_init_matrices_sw,
        penalize_terminal_gaps=True
        )
    end_row_position, end_col_position = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
      
    aligned_seqs1, aligned_seqs2, best_score, identity =_traceback(traceback_matrix, score_matrix, seq1, seq2, end_row_position, end_col_position)
	
    return aligned_seqs1, aligned_seqs2, best_score, identity
	
    

def Damerau_Levenshtein(s1,s2):
    """return Damerau–Levenshtein distance

    Parameters
    ----------
    s1 : list
        list of activities, which is the first sequence to be aligned
    s2 : list
        list of activities, which is the second sequence to be aligned
    
    Returns
    -------
    score : float
        Distance score between two sequence
    
    """		    
    
    maxdist = len(s1) + len(s2)
    shape = (len(s1)+1, len(s2)+1)
    #print(shape)
    score_matrix = np.zeros(shape)
    #score_matrix[0,0] = maxdist
    for i in range(1, shape[0]):
        score_matrix[i, 0] = i
        score_matrix[i, -1] = maxdist   
    for i in range(1, shape[1]):
        score_matrix[0, i] = i
        score_matrix[-1, i] = maxdist
    last_row = {}    
    for row in range(1, shape[0]):
            # Current character in a
        ch_a = s1[row-1]

        # Column of last match on this row: DB in pseudocode
        last_match_col = 0

        for col in range(1, shape[1]):
            # Current character in b
            ch_b = s2[col-1]

            # Last row with matching character
            last_matching_row = last_row.get(ch_b, 0)

            # Cost of substitution
            cost = 0 if ch_a == ch_b else 1

            # Compute substring distance
            score_matrix[row][col] = min(
                score_matrix[row-1][col-1] + cost, # Substitution
                score_matrix[row-1][col] + 1,  # Addition
                score_matrix[row][col-1] + 1,  # Deletion               
                score_matrix[last_matching_row-1][last_match_col-1]+ (row - last_matching_row - 1) + (col - last_match_col - 1)+1)# Transposition
            
            # If there was a match, update last_match_col
            if cost == 0:
                last_match_col = col

        # Update last row for current character
        last_row[ch_a] = row
		
    return score_matrix[-1,-1]

def Levenshtein(s1,s2):
    """return Levenshtein distance

    Parameters
    ----------
    s1 : list
        list of activities, which is the first sequence to be aligned
    s2 : list
        list of activities, which is the second sequence to be aligned
    
    Returns
    -------
    score : float
        Distance score between two sequence
    
    """	
 
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return float(distances[-1])


    
