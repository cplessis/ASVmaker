#######################################################################################################################
#####################        SEQUENCES         ########################################################################
#######################################################################################################################

def max_seq_len(Database):
    """Get maximum sequence length from the data.
    Returns:
        integer: The maximum sequence length
    """        
    max_len = 0
    for access_nb in Database.access_dict:
        if len(Database.access_dict[access_nb]["sequence"]) > max_len:
            max_len = len(Database.access_dict[access_nb]["sequence"])
    return max_len

def min_seq_len(Database):
    """Get the minimum sequence length from the data.
    Returns:
        integer: The minimum sequence length
    """        
    min_len = 100000
    for access_nb in Database.access_dict:
        if len(Database.access_dict[access_nb]["sequence"]) < min_len:
            min_len = len(Database.access_dict[access_nb]["sequence"])
    return min_len

def mean_seq_len(Database):
    """Get the mean sequence length from the data. 
    Returns:
        integer: The mean sequence length
    """        
    len_sum = 0
    for access_nb in Database.access_dict:
        len_sum += len(Database.access_dict[access_nb]["sequence"])
    try: mean_len = len_sum/len(Database.access_dict)
    except ZeroDivisionError: mean_len = 0
    return mean_len

#######################################################################################################################
#####################        TAXON         ############################################################################
#######################################################################################################################

def max_nb_seq_tax(Database):
    """Get the maximum number of sequences per taxon.
    Returns:
        integer: max number of sequences per taxon
    """        
    max_seq = 0
    for taxon in Database.taxon_dict:
        if len(Database.taxon_dict[taxon]) > max_seq:
            max_seq = len(Database.taxon_dict[taxon])
    return max_seq

def min_nb_seq_tax(Database):
    """Get the minimum number of sequences per taxon;
    Returns:
        integer: min number of sequence per taxon
    """        
    min_seq = 100000
    for taxon in Database.taxon_dict:
        if len(Database.taxon_dict[taxon]) < min_seq:
            min_seq = len(Database.taxon_dict[taxon])
    return min_seq

def mean_nb_seq_tax(Database):
    """Get the mean number of sequence per taxon;
    Returns:
        float: mean number of seq per taxon
    """        
    seq_sum = 0
    for taxon in Database.taxon_dict:
        seq_sum += len(Database.taxon_dict[taxon])
    try: mean_seq = seq_sum/len(Database.taxon_dict)
    except ZeroDivisionError: mean_seq = 0
    return mean_seq

#######################################################################################################################
#####################        AMPLICON        ##########################################################################
#######################################################################################################################
    
def get_amplicon_nb(self):
    """Get the total number of amplicons.
    Returns:
        int: the number of amplicons.
    """
    ampl_nb = 0
    for access_nb in self.access_dict:
        if self.access_dict[access_nb]["amplicon"] != "NA":
            ampl_nb += 1
    return ampl_nb

def max_amplicon_len(Database):
    """Get the maximum aplicon length.
    Returns:
        int: max amplicon length
    """        
    max_seq = 0
    for access_nb in Database.access_dict:
        if len(Database.get_amplicon(access_nb)) > max_seq:
            max_seq = len(Database.get_amplicon(access_nb))
    return max_seq

def min_amplicon_len(Database):
    """Get the minimum amplicon length.
    Returns:
        int: min amplicon length
    """        
    min_seq = 100000
    for access_nb in Database.access_dict:
        if (len(Database.get_amplicon(access_nb)) < min_seq) & \
            (Database.get_amplicon(access_nb) != "NA"):
            min_seq = len(Database.get_amplicon(access_nb))
    return min_seq

def mean_amplicon_len(Database):
    """Get the mean amplicon length.
    Returns:
        float: mean amplicon length
    """        
    seq_sum = 0
    na_amplicon = 0
    for access_nb in Database.access_dict:
        if Database.get_amplicon(access_nb) != "NA":
            seq_sum += len(Database.get_amplicon(access_nb))
        else: na_amplicon += 1
    try: mean_seq = seq_sum/(len(Database.access_dict)-na_amplicon)
    except ZeroDivisionError: mean_seq = 0
    return mean_seq

#######################################################################################################################
#####################        COMPLEX       ############################################################################
#######################################################################################################################

def mean_tax_complex(Database):
    """Get the mean number of taxon per complex.
    Returns:
        float: mean number of taxon per complex
    """        
    tax_nb = 0
    for complex_name in Database.complex_dict:
        tax_nb += len(Database.complex_dict[complex_name]["taxon"])
    try: mean_tax = tax_nb/len(Database.complex_dict)    
    except ZeroDivisionError: mean_tax = 0
    return mean_tax

def max_tax_complex(Database):
    """Get the maximum number of taxon per complex;
    Returns:
        integer: max number of taxon per complex
    """        
    max_tax = 0
    for complex_name in Database.complex_dict:
        if len(Database.complex_dict[complex_name]["taxon"]) > max_tax:
            max_tax = len(Database.complex_dict[complex_name]["taxon"])
    return max_tax

def min_tax_complex(Database):
    """Get the minimum number of taxon per complex.
    Returns:
        integer: min number of taxon per complex
    """        
    min_tax = 100000
    for complex_name in Database.complex_dict:
        if len(Database.complex_dict[complex_name]["taxon"]) < min_tax:
            min_tax = len(Database.complex_dict[complex_name]["taxon"])
    if min_tax == 100000: min_tax = 0
    return min_tax

def mean_seq_complex(Database):
    """Get the mean number of sequences per complex.
    Returns:
        float: mean number of sequences per complex
    """        
    seq_nb = 0
    for complex_name in Database.complex_dict:
        seq_nb += len(Database.complex_dict[complex_name]["access"])
    try: mean_seq = seq_nb/len(Database.complex_dict)    
    except ZeroDivisionError: mean_seq = 0
    return mean_seq

def max_seq_complex(Database):
    """Get the maximum number of sequences per complex.
    Returns:
        integer: max number of sequences per complex
    """        
    max_seq = 0
    for complex_name in Database.complex_dict:
        if len(Database.complex_dict[complex_name]["access"]) > max_seq:
            max_seq = len(Database.complex_dict[complex_name]["access"])
    return max_seq

def min_seq_complex(Database):
    """Get the minimum number of sequences per complex.
    Returns:
        integer: min number of sequence per complex
    """        
    min_seq = 100000
    for complex_name in Database.complex_dict:
        if len(Database.complex_dict[complex_name]["access"]) < min_seq:
            min_seq = len(Database.complex_dict[complex_name]["access"])
    if min_seq == 100000: min_tax = 0
    return min_seq