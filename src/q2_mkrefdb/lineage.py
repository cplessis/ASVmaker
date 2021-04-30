import argparse
from os import write

sequences = {}
taxonomy = {}
groups = {}

t_levels = {
        "k": [0, 1],"p": [1, 2],"c": [2, 3],"o": [3, 4],"f": [4, 5],"g": [5, 6],"s": [6, 7],
        'kp': [0, 2],'kc': [0, 3],'ko': [0, 4],'kf': [0, 5],'kg': [0, 6],'ks': [0, 7],
        'pc': [1, 3],'po': [1, 4],'pf': [1, 5],'pg': [1, 6],'ps': [1, 7],
        'co': [2, 4],'cf': [2, 5],'cg': [2, 6],'cs': [2, 7],
        'of': [3, 5],'og': [3, 6],'os': [3, 7],
        'fg': [4, 6],'fs': [4, 7],
        'gs': [5, 7],
    }

def custom_group(access_dict, groups_file, taxonomic_level):
    tl1 = t_levels[taxonomic_level][0]
    tl2 = t_levels[taxonomic_level][1]
    with open(groups_file) as g_file:
        for line in g_file.readlines():
            lineage = line.split()[1:]
            groups[" ".join(lineage)] = line.split()[0]
        new_access_dict = {}
        for access_nb in access_dict:
            new_access_dict[access_nb] = access_dict[access_nb]
            full_lineage = (access_dict[access_nb]["lineage"]+access_dict[access_nb]["taxon"])
            if full_lineage in groups:
                lineage_lvls = full_lineage.split("; ")
                lineage_lvls[tl1:tl2] = groups[full_lineage].split("; ") 
                lineage = "; ".join(lineage_lvls[:-1])+"; "
                taxon = lineage_lvls[-1]+"_"+"group"
                new_access_dict[access_nb]["lineage"] = lineage
                new_access_dict[access_nb]["taxon"] = taxon
    return new_access_dict


if __name__ == "__main__":
    def get_arguments():
        parser = argparse.ArgumentParser(description=__doc__)  
        #---------------------------------------------------
        #                  Init
        #---------------------------------------------------
        parser.add_argument('-is',
                              '--input_sequences',
                              help="File path to the sequences FASTA file to treat.",
                              required=True)
        parser.add_argument('-it',
                              '--input_taxonomy',
                              help="File path to the taxonomy TXT file to treat.",
                              required=True)
        parser.add_argument('-ig',
                              '--input_groups',
                              help="File path to the group TXT file to treat.",
                              required=True)
        parser.add_argument('-tl',
                              '--taxonomic_level',
                              help="Taxonomic level to modify. This can be from kingdom (k) to genus (g).",
                              required=True)
        parser.add_argument('-os',
                              '--output_sequences',
                              help="Name of the sequences FASTA file.",
                              required=True)
        parser.add_argument('-ot',
                              '--output_taxonomy',
                              help="Name of output taxonomy TXT file.",
                              required=True)
        args = parser.parse_args()
        return args

    args = get_arguments()

    input_sequences = args.input_sequences
    input_taxonomy = args.input_taxonomy
    input_groups = args.input_groups
    taxonomic_level = args.taxonomic_level
    output_sequences = args.output_sequences
    output_taxonomy = args.output_taxonomy

    tl1 = t_levels[taxonomic_level][0]
    tl2 = t_levels[taxonomic_level][1]

    # Getting infos from files
    with open(input_sequences) as s_file, open(input_taxonomy) as t_file, \
        open(input_groups) as g_file:
        #Taxonomy
        for line in t_file.readlines():
            seq_id = line.split()[0]
            taxonomy[seq_id] = line.split("\t")[1].strip("\n")
        #Sequence
        for line in s_file.readlines():
            if line[0] == ">": seq_id = line[1:].strip("\n")
            else: sequences[seq_id] = line
        #Groups
        for line in g_file.readlines():
            lineage = line.split()[1:]
            groups[" ".join(lineage)] = line.split()[0]

    with open(output_sequences, "w") as s_file, open(output_taxonomy, "w") as t_file:
        for seq_id in sequences:
            lineage = taxonomy[seq_id]
            if lineage in groups:
                lineage_lvls = lineage.split("; ")
                lineage_lvls[tl1:tl2] =  groups[lineage].split("; ")
                lineage = "; ".join(lineage_lvls)
            s_file.write(">"+seq_id+"\n"+sequences[seq_id]+"\n")
            t_file.write(seq_id+"\t"+lineage+"\n")
