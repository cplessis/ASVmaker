import time
from progress.bar import FillingSquaresBar
from pydna.dseqrecord import Dseqrecord

def del_redund_seq_tax(Database):
    """Remove all the redundant sequence by taxon. All the sequences are compared two by two in each group of taxon only, 
    if they are 100% identical one of them is removed.
    """        
    t0 = time.time()
    with FillingSquaresBar('Deleting redundant sequence variant by taxon ', max= len(Database.seq_dict)) as bar:
        ref_seq_dict = {}
        for taxon in Database.taxon_dict:
            seq_set = set()
            for access_nb in Database.taxon_dict[taxon]:
                d_seq = Dseqrecord(Database.get_sequence(access_nb))
                if all(seq not in seq_set for seq in [d_seq.seq.crick, d_seq.seq.watson]):
                    seq_set.add(d_seq.seq.crick)
                    seq_set.add(d_seq.seq.watson)
                    ref_seq_dict[access_nb] = Database.get_sequence(access_nb)
                bar.next()
        seq_dict = {}
        for access_nb, sequence1 in ref_seq_dict.items():
            seq_dict[Database.get_name(access_nb)] = sequence1
        t1 = time.time()
        print("\n   ==> Deletion done in %f seconds."%(t1 - t0))
        return seq_dict

def del_redund_seq_glob(Database):
    """Remove all the redundant sequence globaly. All the sequences are compared two by two, 
    if they are 100% identical one of them is removed.
    """        
    t0 = time.time()
    with FillingSquaresBar('Deleting redundant sequence global ', max= len(Database.seq_dict)) as bar:
        temp_seq_set = set()
        temp_seq_dict = {}
        for seq_name in Database.seq_dict:
            d_seq = Dseqrecord(Database.seq_dict[seq_name])
            if all(seq not in temp_seq_set for seq in [d_seq.seq.crick, d_seq.seq.watson]):
                temp_seq_set.add(d_seq.seq.crick)
                temp_seq_set.add(d_seq.seq.watson)
                temp_seq_dict[seq_name] = Database.seq_dict[seq_name]
            bar.next()
        t1 = time.time()
        print("\n   ==> Deletion done in %f seconds."%(t1 - t0))
        return temp_seq_dict

def del_na_amplicons(Database):
    """Remove all the sequences which have no amplicon ('NA').
    """        
    t0 = time.time()
    with FillingSquaresBar('Deleting NA amplicons ', max= len(Database.access_dict)) as bar:
        temp_seq_dict = {}
        for access_nb in Database.access_dict:
            if Database.get_amplicon(access_nb) != "NA":
                temp_seq_dict[Database.get_name(access_nb)] = Database.get_sequence(access_nb)
            bar.next()
        t1 = time.time()
        print("\n   ==> Deletion done in %f seconds."%(t1 - t0))
        if len(temp_seq_dict) == 0: 
            print("\n!! WARNING !! No amplicons have been generated from the database.\n\
            Please check your PRIMERS parameters and sequences and CREATE a NEW database.\n\n\
                ==> STOPPED the PROGRAM before END.")
            exit()
        return temp_seq_dict

def del_redund_ampli_glob(Database):
    """Remove all the redundant amplicons globaly. All the amplicons are compared two by two, 
    if they are 100% identical one of them is removed.
    """        
    t0 = time.time()
    with FillingSquaresBar('Deleting redundant amplicons global ', max= len(Database.access_dict)) as bar:
        ampli_set = set()
        seq_dict = {}
        for access_nb in Database.access_dict:
            if Database.get_amplicon(access_nb) not in ampli_set:
                ampli_set.add(Database.get_amplicon(access_nb))
                seq_dict[Database.get_name(access_nb)] = Database.get_sequence(access_nb)
            bar.next()
        t1 = time.time()
        print("\n   ==> Deletion done in %f seconds."%(t1 - t0))
        return seq_dict
    
def del_redund_ampli_tax(Database):
    """Remove all the redundant amplicons by taxon. All the amplicons are compared two by two in each group of taxon only, 
    if they are 100% identical one of them is removed.
    """        
    t0 = time.time()
    with FillingSquaresBar('Deleting redundant amplicons by taxon ', max= len(Database.seq_dict)) as bar:
        ref_access_list = []
        for taxon in Database.taxon_dict:
            seq_set = set()
            for access_nb in Database.taxon_dict[taxon]:
                if Database.get_amplicon(access_nb) not in seq_set:
                    seq_set.add(Database.get_amplicon(access_nb))
                    ref_access_list.append(access_nb)
                bar.next()
        seq_dict = {}
        for access_nb in ref_access_list:
            seq_dict[Database.get_name(access_nb)] = Database.get_sequence(access_nb)
        t1 = time.time()
        print("\n   ==> Deletion done in %f seconds."%(t1 - t0))
        return seq_dict

def clean_dataset(Database, wanted_genus1, wanted_genus2, unverified_bool):
        """Keep in the dataset only the specified genus. The unverified bool remove all
        the species which are ended by '.' such as 'Fusarium sp.' if False is specified. 

        Args:
            wanted_genus1 (string): genus to keep 1
            wanted_genus2 (string): genus to keep 2
            unverified_bool (bool): False to to remove unverfied ('.'), else True.  
        """        
        t0 = time.time()
        with FillingSquaresBar('Cleaning dataset ', max= len(Database.seq_dict)) as bar:
            temp_dict = {}
            for seq_name in Database.seq_dict:
                seq_name_list = seq_name.split(" ")
                if seq_name_list[1] == wanted_genus1:
                    if (unverified_bool == False) & (seq_name_list[2][-1] != "."):
                        temp_dict[seq_name] = Database.seq_dict[seq_name]
                    elif unverified_bool == True:
                        temp_dict[seq_name] = Database.seq_dict[seq_name]
                elif (unverified_bool == True) & (seq_name_list[2] == wanted_genus1):
                    temp_dict[seq_name] = Database.seq_dict[seq_name]
                elif seq_name_list[1] == wanted_genus2:
                    if (unverified_bool == False) & (seq_name_list[2][-1] != "."):                        
                        temp_dict[seq_name] = Database.seq_dict[seq_name]
                    elif unverified_bool == True:
                        temp_dict[seq_name] = Database.seq_dict[seq_name]
                elif (unverified_bool == True) & (seq_name_list[2] == wanted_genus2):
                    temp_dict[seq_name] = Database.seq_dict[seq_name]
                bar.next()            
            t1 = time.time()
            print("\n   ==> Dataset cleaned in %f seconds."%(t1 - t0))
            if len(temp_dict) == 0: 
                print("\n!! WARNING !! After cleaning the database there is no sequence left.\n\
                Please check your GENUS parameter and run the command again.\n\n\
                    ==> STOPPED the PROGRAM before END.")
                exit()
            return temp_dict

def group_by_complex(self):
    """Change the lineage and taxon of all species which are in a species complex.
    The taxon name become the complex name and the lineage stop at this level. 
    By this way, species are gathered in complex.  
    """        
    t0 = time.time()
    with FillingSquaresBar('Assembling taxon by complex ', max= len(self.seq_dict)) as bar:
        for access_nb in self.access_dict:
            lineage = self.get_lineage(access_nb)
            if lineage == "NA": continue
            if lineage.split("; ")[-2].split(" ")[-1] in {"complex","group"}:
                taxon = ""
                if lineage.split("; ")[-2].split(" ")[0] == "unclassified":
                    taxon = lineage.split("; ")[-2].split(" ")[1:3]
                    taxon.append("complex")
                    taxon = "_".join(taxon)
                    lineage = "; ".join(lineage.split("; ")[:-3])+"; "
                else:
                    taxon = lineage.split("; ")[-2].split(" ")[:2]
                    taxon.append("complex")
                    taxon = "_".join(taxon)
                    lineage = "; ".join(lineage.split("; ")[:-2])+"; "
                self.access_dict[access_nb]["taxon"] = taxon
                self.access_dict[access_nb]["lineage"] = lineage
            bar.next()
        t1 = time.time()
        print("\n   ==> Taxon assembled by complex in %f seconds."%(t1 - t0))  
        return self.access_dict