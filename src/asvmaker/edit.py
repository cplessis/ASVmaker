def remove_by_id(Database, id_list_csv):
    """Remove all the sequences with IDs on id_list_csv.
        The id_list_csv FILE must be as one line = 'seq_id'. 

        Args:
            id_list_csv (str): ID list CSV or TXT file
        """ 
    seq_dict = {}
    id_set = set()
    with open(id_list_csv) as file:
        for line in file.readlines(): id_set.add(line.strip("\n"))
        for access_nb in Database.access_dict:
            if access_nb not in id_set:
                seq_dict[Database.get_name(access_nb)] = Database.get_sequence(access_nb)
    return seq_dict

def remove_by_taxon(Database, tax_list_csv):
    """Remove all the sequences with IDs on id_list_csv.
        The tax_list_csv FILE must be as one line = 'Genus_specie'. 

        Args:
            tax_list_csv (str): taxon list CSV or TXT file
    """ 
    tax_set = set()
    with open(tax_list_csv) as file:
        for line in file.readlines(): tax_set.add(line.strip("\n"))
        for taxon in tax_set:
            if taxon in Database.taxon_dict:
                for access_nb in Database.taxon_dict[taxon]:
                    Database.seq_dict.pop(Database.get_name(access_nb))
    return Database.seq_dict

def rename_by_id(Database, id_list_csv):
    """Rename all the sequences with IDs on id_list_csv.
        The id_list_csv FILE must be as one line = 'seq_id   new_taxon_name'. 

        Args:
            id_list_csv (str): ID list CSV or TXT file
        """
    id_set = set()
    with open(id_list_csv) as file:
        for line in file.readlines(): id_set.add(line.strip("\n"))
        for seq_id in id_set:
            Database.access_dict[seq_id.split()[0]]["taxon"] = seq_id.split()[1]
    return Database.access_dict

def group_by_id(Database, shared_ext_csv):
    """Group all the sequences with IDs on shared_ext_csv on a commune taxon name.
        The taxon name will become the SA_id of the group. Only one amplicon will represent
        the group after this command is runned. The shared_ext_csv FILE must be the one
        generated with EXPORT shared ampl.

        Args:
            shared_ext_csv (str): shared_ext FILE generated with EXPORT shared ampl.
        """
    id_set = set()
    remove_list = []
    seq_dict = {}
    with open(shared_ext_csv) as file:
        for line in file.readlines(): id_set.add(line.strip("\n"))
        for sa_seq in id_set:
            access_nb = sa_seq.split()[1].split("|")[0]
            genus = Database.get_taxon(access_nb).split("_")[0]
            #except IndexError: genus = Database.get_lineage(access_nb).split("_")[0]
            Database.access_dict[access_nb]["taxon"] = genus+"_"+sa_seq.split()[0]
            remove_list += [i.split("|")[0] for i in sa_seq.split()[2:]]
        for access_nb in Database.access_dict:
            if access_nb not in remove_list:
                seq_dict[Database.get_name(access_nb)] = Database.get_sequence(access_nb)
    return Database.access_dict, seq_dict