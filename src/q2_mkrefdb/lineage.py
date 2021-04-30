import argparse

sequences = {}
taxonomy = {}
groups = {}

def custom_group(Database, groups_file):
    print("LINEAGE")
    with open(groups_file) as g_file:
        for line in g_file.readlines():
            taxon = line.split()[1]
            groups[taxon] = line.split()[0]
        new_access_dict = {}
        for access_nb in Database.access_dict:
            new_access_dict[access_nb] = Database.access_dict[access_nb]
            access_taxon = Database.access_dict[access_nb]["taxon"]
            if access_taxon in groups:
                new_tax = groups[access_taxon]+"_"+"group"
                new_access_dict[access_nb]["taxon"] = new_tax
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
                lineage_lvls[-1] =  groups[lineage].split("; ")
                lineage = "; ".join(lineage_lvls)
            s_file.write(">"+seq_id+"\n"+sequences[seq_id]+"\n")
            t_file.write(seq_id+"\t"+lineage+"\n")
