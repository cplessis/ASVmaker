from Bio import Align
from Bio.Seq import Seq

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5
aligner.mismatch_score = 0
# Algortithm by default for local alignement with the above parameters : "Gotoh"
# These parameters are used in order to return very high identity score when no gaps.


"""This module can be used independantly from the main function of the package. 
"""


def __extract_position(primer_type, primer, alignment):
    """Find the primer start (if primer type == 'rv') or end position (if primer type == 'fw') 
    on the sequence with a given alignment. The function also return the score of the alignment.

    Args:
        primer_type (string): 'fw' if forward, 'rv' if reverse
        primer (string): Primer sequence
        alignment (Bio.Align.PairwiseAlignments): Results from a pairwaise alignment

    Returns:
        list: Informations about the alignment [position, score]
    """    
    seq_pos = alignment.aligned[0][0][0]
    if primer_type == "fw": 
        prim_pos = alignment.aligned[1][0][0]
        position = [seq_pos+len(primer)-prim_pos]
    else:     
        prim_pos = alignment.aligned[1][0][0]
        position = [seq_pos-prim_pos]
    position.append(alignment.score)
    return position

def __primer_infos(primer_type, primer, sequence):
    """Get the the position and score of the best alignment score between 'primer' and 'sequence'.
    End position if 'fw' and start one if 'rv'. 

    Args:
        primer_type (str): 'fw' if forward, 'rv' if reverse 
        primer (str): Primer sequence
        sequence (str): Sequence

    Returns:
        list: The best position and score [postion, score] 
    """    
    al = aligner.align(sequence, primer)
    if len(al) == 1: pos_infos = __extract_position(primer_type, primer, al[0])
    elif len(al) > 1:
        highest, iterator = 0, 0
        for alignement in al:
            if alignement.score > highest:
                highest = alignement.score
                index = iterator
            iterator += 1
        pos_infos = __extract_position(primer_type, primer, al[index])
    else: pos_infos = [0,0]
    return pos_infos

def primer_condition(sequence, forward_primer, reverse_primer, \
    fw_mismatch_tol, rv_mismatch_tol, fw_position, rv_position):
    """Verify if the primers alignments respect the wanted conditions 
    (mismatch tolerance and identical 3 last 3'-end nucleotides). The 
    primers must have at leat three last nuclotide as identical as in 
    the sequence. 

    Args:
        sequence (str): Sequence
        forward_primer (str): forward primer sequence
        reverse_primer (str): reverse primer sequence
        rv_mismatch_tol (int): reverse primer mismatch tolerance
        fw_mismatch_tol (int): forward primer mismatch tolerance
        fw_position (int): primer end position on sequence 
        rv_position (int): primer start position on sequence

    Returns:
        boolean: True if conditions are respected, else False.
    """    
    bool = ((len(forward_primer) - fw_mismatch_tol) <= fw_position[1])
    bool += ((len(reverse_primer) - rv_mismatch_tol) <= rv_position[1])
    bool += (sequence[fw_position[0]-3:fw_position[0]] == forward_primer[-3:])
    bool += (sequence[rv_position[0]:rv_position[0]+3] == reverse_primer[:3])
    return bool == 4

def amplify(sequence, forward_primer,reverse_primer, fw_mismatch_tol, rv_mismatch_tol, trim_prim):
    """Find the best positions of the inputed primers on the sequence in order to return an amplicon.
    The amplicon primers sites are trimmed.
    
    Args:
        sequence (string): Sequence to amplify
        forward_primer (str): Forward primer sequence
        reverse_primer (str): Reverse primer sequence
        fw_mismatch_tol (int): Forward primer mismatch tolerance
        rv_mismatch_tol (int): Reverse primer mismatch tolerance
        trim_prim(bool): True to trim primers on amplicons, else False

    Returns:
        string: The amplicon sequence
    """    
    reverse_primer = Seq(reverse_primer).reverse_complement()
    fw_position = __primer_infos("fw", forward_primer, sequence)
    rv_position = __primer_infos("rv", reverse_primer, sequence)
    if primer_condition(sequence, forward_primer, reverse_primer, \
        fw_mismatch_tol, rv_mismatch_tol, fw_position, rv_position): 
        if trim_prim: amplicon = str(sequence[fw_position[0]:rv_position[0]])
        else: amplicon = str(sequence[fw_position[0]-len(forward_primer)\
            :rv_position[0]+len(reverse_primer)])
    else:
        sequence = Seq(sequence).reverse_complement()
        fw_position = __primer_infos("fw", forward_primer, sequence)
        rv_position = __primer_infos("rv", reverse_primer, sequence)
        if primer_condition(sequence, forward_primer, reverse_primer, \
        fw_mismatch_tol, rv_mismatch_tol, fw_position, rv_position):
            if trim_prim: amplicon = str(sequence[fw_position[0]:rv_position[0]])
            else: amplicon = str(sequence[fw_position[0]-len(forward_primer)\
                :rv_position[0]+len(reverse_primer)])
        else: amplicon = "NA"
    return amplicon
