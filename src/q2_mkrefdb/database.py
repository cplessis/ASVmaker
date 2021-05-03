import re, time, random
from . import utils, filter
from .pcr import amplify
from .lineage import custom_group
from progress.bar import FillingSquaresBar
from progress.spinner import PixelSpinner
from pydna.dseqrecord import Dseqrecord
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

###########################################################################################################################################
#################### CLASS FastaFile ######################################################################################################
###########################################################################################################################################

class Database:
    """The class Database own all the function and attributes to filter data from a FASTA file. 
    Especially daata from NCBI 'nt', EBI 'ena', DDBJ 'arsa'. 
    In order to use it you will need a FASTA file from one of these database and 2 other one in which
    you need to give the primer for the amplification.

    The Database Class can be used directly but it is more adapted to use it throw the q2_mkrefdb modules (command lines). 
    """    
    def create(self, fasta_file, origin_database, forward_primer_fasta, reverse_primer_fasta, \
        fw_mismatch_tol = 3, rv_mismatch_tol = 3, trim_primers = True):
        """Initiate the Database Class.

        Args:
            fasta_file (string): path to the FASTA file to treat
            origin_database (string): 'ncbi', 'ebi' or 'ddbj'
            forward_primer_fasta (string): path to the forward primer FASTA file
            reverse_primer_fasta (string): path to the reverse primer FASTA file
            fw_mismatch_tol (int, optional): forward primer mismatch tolerance for annealing. Defaults to 3.
            rv_mismatch_tol (int, optional): reverse primer mismatch tolerance for annealing. Defaults to 3.
            trim_primers = True to trim primers on amplicons, else False.
        """              
        self.file_name = fasta_file
        self.origin_database = origin_database
        self.fw_mismatch_tol = fw_mismatch_tol
        self.rv_mismatch_tol = rv_mismatch_tol
        self.trim_primers = trim_primers
        self.forward_primer_list = self.__get_seq_from_fasta(forward_primer_fasta, False)
        self.reverse_primer_list = self.__get_seq_from_fasta(reverse_primer_fasta, False)
        self.seq_dict = self.__get_seq_from_fasta(fasta_file, True)
        self.access_dict = self.__make_access_dict()
        self.taxon_dict = self.__make_taxon_dict()
        self.__make_amplicon_dict()
        self.__make_lineage_dict()
        self.complex_dict = self.__make_complex_dict()
        return

    def import_db(self, database_json):
        db = utils.open_from_json(database_json)
        self.access_dict = db
        self.seq_dict = {}
        for access_nb in self.access_dict:
            self.seq_dict[self.get_name(access_nb)] = self.get_sequence(access_nb)
        self.taxon_dict = self.__make_taxon_dict()
        self.complex_dict = self.__make_complex_dict()


    def __get_seq_from_fasta(self, fasta_file, bool_multiple_seq):
        """Extract the sequences from the FASTA file and initiate the data.

        Args:
            fasta_file (string): the path to the FASTA file
            bool_multiple_seq (boolean): True if the file is the sequence file, else False for primer file.

        Returns:
            list or dict: list if bool_multiple_seq == True, else dict with sequences names as key and sequences as value.
        """                
        seq_dict = {}
        seq_ref = ""
        seq_list = []
        with open(fasta_file) as file:
            for line in file.readlines():
                line = line.rstrip("\n")
                if line == "": pass
                elif (line[0] == ">") & (bool_multiple_seq == True):
                    if line.split(" ")[1] not in {"UNVERIFIED:", "Uncultured"}:
                        seq_dict[line] = ""
                        seq_ref = line
                    else: seq_ref = "pass"
                elif (line[0] == ">") & (bool_multiple_seq == False):
                    seq_list.append(line)
                    seq_ref = line
                elif seq_ref != "pass":
                    if bool_multiple_seq == True: seq_dict[seq_ref] += line
                    else: seq_list.append(line)
        if bool_multiple_seq == True: answer = seq_dict
        else: answer = seq_list
        return answer

    def __make_access_dict(self):
        """Create the access_dict dictionnary in which all the informations will be stored.
        Add the sequences to the dictionnary only if letters [ATCGatcg] in sequence.

        Returns:
            dictionnary: The access_dict dictionnary
        """        
        access_dict = {}
        new_seq_dict = dict(self.seq_dict)
        for seq_name in self.seq_dict:
            if re.search("[^ATCGatcg]", self.seq_dict[seq_name]) == None:
                access_dict[self.__get_access_num(seq_name)] = \
                    {"name":seq_name, "sequence":self.seq_dict[seq_name], \
                        "taxon":seq_name.split()[1]+"_"+seq_name.split()[2]}
            else: new_seq_dict.pop(seq_name)
        self.seq_dict = new_seq_dict    
        return access_dict

    def __update_data(self):
        """Update the information of the class.
        """        
        new_access_dict = {}
        for seq_name in self.seq_dict:
            access_nb = self.__get_access_num(seq_name)
            new_access_dict[access_nb] = self.access_dict[access_nb]
        self.access_dict = new_access_dict
        self.taxon_dict = self.__make_taxon_dict()
        self.complex_dict = self.__make_complex_dict()

    def __max_seq_len(self):
        """Get maximum sequence length from the data.

        Returns:
            integer: The maximum sequence length
        """        
        max_len = 0
        final = ""
        for access_nb in self.access_dict:
            if len(self.access_dict[access_nb]["sequence"]) > max_len:
                max_len = len(self.access_dict[access_nb]["sequence"])
                final = access_nb
        return max_len
    
    def __min_seq_len(self):
        """Get the minimum sequence length from the data.

        Returns:
            integer: The minimum sequence length
        """        
        min_len = 100000
        for access_nb in self.access_dict:
            if len(self.access_dict[access_nb]["sequence"]) < min_len:
                min_len = len(self.access_dict[access_nb]["sequence"])
        return min_len

    def __mean_seq_len(self):
        """Get the mean sequence length from the data. 

        Returns:
            integer: The mean sequence length
        """        
        len_sum = 0
        for access_nb in self.access_dict:
            len_sum += len(self.access_dict[access_nb]["sequence"])
        try: mean_len = len_sum/len(self.access_dict)
        except ZeroDivisionError: mean_len = 0
        return mean_len

    def __make_taxon_dict(self):
        """Create the taxon dictionnary.

        Returns:
            dictionnary: The taxon dictionnary.
        """        
        taxon_dict = {}
        for access_nb in self.access_dict:
            try:
                taxon_dict[self.get_taxon(access_nb)].add(access_nb)
            except KeyError:
                taxon_dict[self.get_taxon(access_nb)] = {access_nb}
        return taxon_dict

    def __make_lineage_dict(self):
        """Add the lineage key to the access_dict. This function use requires an internet connection.
        This function is made of 2 other ones.
        """        
        with FillingSquaresBar('Getting species lineages from EBI taxonomy online ressources ', \
            max= len(self.taxon_dict)) as bar:
                t0 = time.time()  
                self.__modified_taxon = set()
                for taxon in self.taxon_dict:
                    name_list = taxon.split("_")
                    self.__get_lineage_norm(taxon, name_list)
                    bar.next()
                if len(self.__na_tax_str) > 0: print("\nThe following taxons have NO lineage : \n", self.__na_tax_str, "\n")
                self.__update_data()
                t1 = time.time()
                print("\n   ==> Got species lineages in %f seconds."%(t1 - t0))
 
    def __get_lineage_norm(self, taxon, name_list):
        """Used in __make_lineage_dict function. Find the lineage of a specie from EBI. 

        Args:
            taxon (string): The taxon name as 'Genus_specie'
            name_list (list): The taxon name as [Genus, specie]
        """        
        self.__na_tax_str = ""
        url = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"+"%20".join(name_list)
        try: lineage = utils.get_var_from_url(url, "json")[0]['lineage']
        except: lineage = "NA"
        if (taxon[-1] == ".") & (lineage == "NA"): 
            for access_nb in self.taxon_dict[taxon]:
                lineage = self.__complete_lineage(access_nb, url)
                if str(taxon+" => "+self.get_taxon(access_nb)) not in self.__modified_taxon:
                    self.__modified_taxon.add(taxon+" => "+self.get_taxon(access_nb))
                if lineage == "NA": self.__na_tax_str += (taxon+"\n")
                self.access_dict[access_nb]["lineage"] = lineage    
        else:
            for access_nb in self.taxon_dict[taxon]:
                if lineage == "NA": 
                    lineage = self.__complete_lineage(access_nb, url)
                    if str(taxon+" => "+self.get_taxon(access_nb)) not in self.__modified_taxon:
                        self.__modified_taxon.add(taxon+" => "+self.get_taxon(access_nb))
                if lineage == "NA": self.__na_tax_str += (taxon+"\n")
                self.access_dict[access_nb]["lineage"] = lineage

    def __complete_lineage(self, access_nb, url):
        """Third part of the _make_lineage_dict function, used if the specie is not found on the EBI.
        Research the specie on the NCBI database  with Entrez (longer than EBI) then search again the
        lineage on EBI. 

        Args:
            access_nb (string): Sequence accession number related to the specie.
            url (string): url to EBI website

        Returns:
            string: lineage if found, else 'NA'
        """        
        Entrez.email = 'someuser@mail.com'
        handle = Entrez.efetch(db="nucleotide", id=access_nb, rettype="gb", retmode="text")
        result=handle.read().split('\n')
        for line in result:
            if 'ORGANISM' in line :
                self.access_dict[access_nb]["taxon"] = "_".join(line.split()[1:])
                url = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"+"%20".join(line.split()[1:])
        try: return utils.get_var_from_url(url, "json")[0]['lineage']
        except: return "NA"

    def __make_complex_dict(self):
        """Create the complex dictionnary from the lineage informations.

        Returns:
            dictionnary: complex dictionnary
        """        
        complex_dict = {}
        for access_nb in self.access_dict:
            if self.get_lineage(access_nb) != "NA":
                complex_name = self.get_lineage(access_nb).split("; ")[-2]
                if self.get_taxon(access_nb).split("_")[-1] in {"complex", "group"}:
                    complex_name = self.access_dict[access_nb]["taxon"]
                    taxon = "_".join(self.get_name(access_nb).split()[1:3])
                    if  complex_name not in complex_dict:
                        complex_dict[complex_name] = {"access":{access_nb}, "taxon":{taxon}}
                    else:
                        complex_dict[complex_name]["access"].add(access_nb)
                        if taxon not in complex_dict[complex_name]["taxon"]:
                            complex_dict[complex_name]["taxon"].add(taxon)
                elif complex_name.split()[-1] in {"complex", "group"}:
                    taxon = self.get_taxon(access_nb)
                    if complex_name not in complex_dict:
                        complex_dict[complex_name] = {"access":{access_nb}, "taxon":{taxon}}
                    else:
                        complex_dict[complex_name]["access"].add(access_nb)
                        if taxon not in complex_dict[complex_name]["taxon"]:
                            complex_dict[complex_name]["taxon"].add(taxon)
        return complex_dict

    def custom_complex(self, groups_file):
        self.access_dict = custom_group(self, groups_file)
        self.__update_data()

    def __mean_tax_complex(self):
        """Get the mean number of taxon per complex.

        Returns:
            float: mean number of taxon per complex
        """        
        tax_nb = 0
        for complex_name in self.complex_dict:
            tax_nb += len(self.complex_dict[complex_name]["taxon"])
        try: mean_tax = tax_nb/len(self.complex_dict)    
        except ZeroDivisionError: mean_tax = 0
        return mean_tax

    def __max_tax_complex(self):
        """Get the maximum number of taxon per complex;

        Returns:
            integer: max number of taxon per complex
        """        
        max_tax = 0
        for complex_name in self.complex_dict:
            if len(self.complex_dict[complex_name]["taxon"]) > max_tax:
                max_tax = len(self.complex_dict[complex_name]["taxon"])
        return max_tax

    def __min_tax_complex(self):
        """Get the minimum number of taxon per complex.

        Returns:
            integer: min number of taxon per complex
        """        
        min_tax = 100000
        for complex_name in self.complex_dict:
            if len(self.complex_dict[complex_name]["taxon"]) < min_tax:
                min_tax = len(self.complex_dict[complex_name]["taxon"])
        if min_tax == 100000: min_tax = 0
        return min_tax

    def __mean_seq_complex(self):
        """Get the mean number of sequences per complex.

        Returns:
            float: mean number of sequences per complex
        """        
        seq_nb = 0
        for complex_name in self.complex_dict:
            seq_nb += len(self.complex_dict[complex_name]["access"])
        try: mean_seq = seq_nb/len(self.complex_dict)    
        except ZeroDivisionError: mean_seq = 0
        return mean_seq

    def __max_seq_complex(self):
        """Get the maximum number of sequences per complex.

        Returns:
            integer: max number of sequences per complex
        """        
        max_seq = 0
        for complex_name in self.complex_dict:
            if len(self.complex_dict[complex_name]["access"]) > max_seq:
                max_seq = len(self.complex_dict[complex_name]["access"])
        return max_seq

    def __min_seq_complex(self):
        """Get the minimum number of sequences per complex.

        Returns:
            integer: min number of sequence per complex
        """        
        min_seq = 100000
        for complex_name in self.complex_dict:
            if len(self.complex_dict[complex_name]["access"]) < min_seq:
                min_seq = len(self.complex_dict[complex_name]["access"])
        if min_seq == 100000: min_tax = 0
        return min_seq

    def __get_access_num(self, seq_name):
        """Get the sequence accession number from the seqeunce name.

        Args:
            seq_name (string): raw sequence name from seq_dict

        Returns:
            string: accession number
        """        
        access_num = ""
        if self.origin_database == "ncbi":
            access_num = seq_name.split(" ")[0].strip(">")[:-2]
        if self.origin_database == "embl":
            access_num = seq_name.split("|")[1]
        if self.origin_database == "ddbj":
            access_num = seq_name.split("|")[0].strip(">")
        return access_num

    def __max_nb_seq_tax(self):
        """Get the maximum number of sequences per taxon.

        Returns:
            integer: max number of sequences per taxon
        """        
        max_seq = 0
        for taxon in self.taxon_dict:
            if len(self.taxon_dict[taxon]) > max_seq:
                max_seq = len(self.taxon_dict[taxon])
        return max_seq

    def __min_nb_seq_tax(self):
        """Get the minimum number of sequences per taxon;

        Returns:
            integer: min number of sequence per taxon
        """        
        min_seq = 100000
        for taxon in self.taxon_dict:
            if len(self.taxon_dict[taxon]) < min_seq:
                min_seq = len(self.taxon_dict[taxon])
        return min_seq

    def __mean_nb_seq_tax(self):
        """Get the mean number of sequence per taxon;

        Returns:
            float: mean number of seq per taxon
        """        
        seq_sum = 0
        for taxon in self.taxon_dict:
            seq_sum += len(self.taxon_dict[taxon])
        try: mean_seq = seq_sum/len(self.taxon_dict)
        except ZeroDivisionError: mean_seq = 0
        return mean_seq

    def __make_amplicon_dict(self):
        """Add the amplicons keys to the access_dict. Each amplicon is created thanks to the pcr module.
        The default values of mismatch tolerance is 3. 
        """        
        spinner = PixelSpinner("Loading data ")
        for seq_name, sequence in self.seq_dict.items():
            try:
                amplicon_seq = amplify(sequence, self.forward_primer_list[1], self.reverse_primer_list[1], \
                    self.fw_mismatch_tol, self.rv_mismatch_tol, self.trim_primers)
                self.access_dict[self.__get_access_num(seq_name)]["amplicon"] = amplicon_seq
            except ValueError:
                self.access_dict[self.__get_access_num(seq_name)]["amplicon"] = "NA"
            spinner.next()

    def __max_amplicon_len(self):
        """Get the maximum aplicon length.

        Returns:
            int: max amplicon length
        """        
        max_seq = 0
        for access_nb in self.access_dict:
            if len(self.get_amplicon(access_nb)) > max_seq:
                max_seq = len(self.get_amplicon(access_nb))
        return max_seq

    def __min_amplicon_len(self):
        """Get the minimum amplicon length.

        Returns:
            int: min amplicon length
        """        
        min_seq = 100000
        for access_nb in self.access_dict:
            if (len(self.get_amplicon(access_nb)) < min_seq) & \
                (self.get_amplicon(access_nb) != "NA"):
                min_seq = len(self.get_amplicon(access_nb))
        return min_seq

    def __mean_amplicon_len(self):
        """Get the mean amplicon length.

        Returns:
            float: mean amplicon length
        """        
        seq_sum = 0
        na_amplicon = 0
        for access_nb in self.access_dict:
            if self.get_amplicon(access_nb) != "NA":
                seq_sum += len(self.get_amplicon(access_nb))
            else: na_amplicon += 1
        try: mean_seq = seq_sum/(len(self.access_dict)-na_amplicon)
        except ZeroDivisionError: mean_seq = 0
        return mean_seq

#########   INFOS  ######### 

    def get_info(self, title):
        """Return a string to print to see all the actual informations about the data. 

        Args:
            title (string): Title of the output string

        Returns:
            string: description of the actual data informations
        """        
        message = \
        "\n****************************************\n"+\
        "===> "+title.upper()+"\n"+\
        "File name : "+self.file_name+"\n\n"+\
        "--------------SEQUENCES---------------------\n"+\
        "Number of sequences : "+str(len(self.seq_dict))+"\n"+ \
        "Maximum sequences length : "+str(self.__max_seq_len())+"\n"+\
        "Minimum sequences length : "+str(self.__min_seq_len())+"\n"+\
        "Mean sequences length : "+str(self.__mean_seq_len())+"\n\n"+\
        "---------------TAXONS-----------------------\n"+\
        "Number of different taxons : "+str(len(self.taxon_dict))+"\n"+\
        "Maximum number of sequences per taxon : "+str(self.__max_nb_seq_tax())+"\n"+\
        "Minimum number of sequences per taxon : "+str(self.__min_nb_seq_tax())+"\n"+\
        "Mean number of sequences per taxon : "+str(self.__mean_nb_seq_tax())+"\n\n"+\
        "---------------COMPLEXS---------------------\n"+\
        "Number of differents complex : "+str(len(self.complex_dict))+"\n"+\
        "Maximum number of taxons per complex : "+str(self.__max_tax_complex())+"\n"+\
        "Minimum number of taxons per complex : "+str(self.__min_tax_complex())+"\n"+\
        "Mean number of taxons per complex : "+str(self.__mean_tax_complex())+"\n"+\
        "Maximum number of sequences per complex : "+str(self.__max_seq_complex())+"\n"+\
        "Minimum number of sequences per complex : "+str(self.__min_seq_complex())+"\n"+\
        "Mean number of sequences per complex : "+str(self.__mean_seq_complex())+"\n\n"+\
        "--------------AMPLICONS---------------------\n"+\
        "Number of amplicons : "+str(self.__get_amplicon_nb())+"\n"+\
        "Maximum amplicons length : "+str(self.__max_amplicon_len())+"\n"+\
        "Minimum amplicons length : "+str(self.__min_amplicon_len())+"\n"+\
        "Mean amplicons length : "+str(self.__mean_amplicon_len())+"\n"+\
        "Number of redundant amplicons : "+str(len(self.get_shared_amplicons()))+"\n"+\
        "****************************************\n"
        return message

    def get_sequence(self, access_number):
        """Get sequence from the accession number.

        Args:
            access_number (string): accesssion number

        Returns:
            string: sequence
        """        
        return self.access_dict[access_number]["sequence"]

    def get_amplicon(self, access_number):
        """Get the amplicon from the accesssion number.

        Args:
            access_number (string): accession number

        Returns:
            string: amplicon
        """        
        return self.access_dict[access_number]["amplicon"]

    def get_name(self, access_number):
        """Get the sequence raw name (description) from the accession number.

        Args:
            access_number (string): accession number

        Returns:
            string: sequence name (description).
        """        
        return self.access_dict[access_number]["name"]

    def get_taxon(self, access_number):
        """Get the taxon from the accession number.

        Args:
            access_number (string): accession number

        Returns:
            string: taxon
        """        
        return self.access_dict[access_number]["taxon"]

    def get_all_taxons(self):
        """Get the set of all taxons.

        Returns:
            set: all taxons
        """        
        return {taxon for taxon in self.taxon_dict}

    def get_genus(self, access_number):
        """Get the genus from accession number.

        Args:
            access_number (string): accession number

        Returns:
            string: genus related to accession
        """        
        return self.get_taxon(access_number).split("_")[0]
    
    def get_specie(self, access_number):
        """Get specie related to accession number.

        Args:
            access_number (string): accession number

        Returns:
            string: specie related to accession number
        """        
        return self.get_taxon(access_number).split("_")[1]
        
    def get_all_genus(self):
        """Get set of all genus.

        Returns:
            set: all genus
        """        
        return {taxon.split("_", 1)[0] for taxon in self.taxon_dict}

    def get_na_amplicons(self):
        """Get a dictionnary with accession number and all 'NA' amplicons.

        Returns:
            dict: all accessions of 'NA' amplicons as {accession : taxon}
        """        
        na_amplicon_dict = {}
        for access_nb in self.access_dict:
            if self.get_amplicon(access_nb) == "NA":
                na_amplicon_dict[access_nb] = self.get_taxon(access_nb)
        return na_amplicon_dict

    def get_shared_amplicons(self):
        """Get a list of all the sequences which have an identical amplicon.

        Returns:
            list: list of strings which describes sequences with shared amplicon
        """        
        shared_ampli_list = []
        checked_sequences = {}
        for access_nb in self.access_dict:
            amplicon = self.get_amplicon(access_nb)
            if  amplicon not in checked_sequences:
                checked_sequences[amplicon] = {access_nb}
            else:
                checked_sequences[amplicon].add(access_nb)
        for amplicon in checked_sequences:
            if len(checked_sequences[amplicon]) > 1:
                shared_amplicon = ""
                for access_nb in checked_sequences[amplicon]:
                    shared_amplicon += access_nb+"_"+self.get_taxon(access_nb)+"  <==>  "
                shared_ampli_list.append(shared_amplicon)
        return shared_ampli_list
    
    def __get_amplicon_nb(self):
        """Get the total number of amplicons.

        Returns:
            int: the number of amplicons.
        """
        ampl_nb = 0
        for access_nb in self.access_dict:
            if self.access_dict[access_nb]["amplicon"] != "NA":
                ampl_nb += 1
        return ampl_nb

    def get_lineage(self, access_number):
        """Get the lineage related to the accession number.

        Args:
            access_number (string): accession number

        Returns:
            string: lineage
        """        
        return self.access_dict[access_number]["lineage"]

#########   FILTER  ######### 

    def clean_dataset(self, wanted_genus1, wanted_genus2, unverified_bool):
        """Keep in the dataset only the specified genus. The unverified bool remove all
        the species which are ended by '.' such as 'Fusarium sp.' if False is specified. 

        Args:
            wanted_genus1 (string): genus to keep 1
            wanted_genus2 (string): genus to keep 2
            unverified_bool (bool): False to to remove unverfied ('.'), else True.  
        """        
        t0 = time.time()
        with FillingSquaresBar('Cleaning dataset ', max= len(self.seq_dict)) as bar:
            temp_dict = {}
            for seq_name in self.seq_dict:
                seq_name_list = seq_name.split(" ")
                if seq_name_list[1] == wanted_genus1:
                    if (unverified_bool == False) & (seq_name_list[2][-1] != "."):
                        temp_dict[seq_name] = self.seq_dict[seq_name]
                    elif unverified_bool == True:
                        temp_dict[seq_name] = self.seq_dict[seq_name]
                elif (unverified_bool == True) & (seq_name_list[2] == wanted_genus1):
                    temp_dict[seq_name] = self.seq_dict[seq_name]
                elif seq_name_list[1] == wanted_genus2:
                    if (unverified_bool == False) & (seq_name_list[2][-1] != "."):                        
                        temp_dict[seq_name] = self.seq_dict[seq_name]
                    elif unverified_bool == True:
                        temp_dict[seq_name] = self.seq_dict[seq_name]
                elif (unverified_bool == True) & (seq_name_list[2] == wanted_genus2):
                    temp_dict[seq_name] = self.seq_dict[seq_name]
                bar.next()            
            self.seq_dict = temp_dict
            self.__update_data()
            t1 = time.time()
            print("\n   ==> Dataset cleaned in %f seconds."%(t1 - t0))              

    def del_redund_seq_glob(self):
        """Remove all the redundant sequence globaly. All the sequences are compared two by two, 
        if they are 100% identical one of them is removed.
        """
        self.seq_dict = filter.del_redund_seq_glob(self)
        self.__update_data()

    def del_redund_seq_tax(self):
        """Remove all the redundant sequence by taxon. All the sequences are compared two by two in each group of taxon only, 
        if they are 100% identical one of them is removed.
        """        
        self.seq_dict = filter.del_redund_seq_tax(self)
        self.__update_data()

    def del_na_amplicons(self):
        """Remove all the sequences which have no amplicon ('NA').
        """        
        self.seq_dict = filter.del_na_amplicons(self)
        self.__update_data()

    def del_redund_ampli_glob(self):
        """Remove all the redundant amplicons globaly. All the amplicons are compared two by two, 
        if they are 100% identical one of them is removed.
        """        
        self.seq_dict = filter.del_redund_ampli_glob(self)
        self.__update_data()

    def del_redund_ampli_tax(self):
        """Remove all the redundant amplicons by taxon. All the amplicons are compared two by two in each group of taxon only, 
        if they are 100% identical one of them is removed.
        """        
        self.seq_dict = filter.del_redund_ampli_tax(self)
        self.__update_data()

    def group_by_complex(self):
        """Change the lineage and taxon of all species which are in a species complex.
        The taxon name become the complex name and the lineage stop at this level. 
        By this way, species are gathered in complex.  
        """        
        t0 = time.time()
        with FillingSquaresBar('Assembling taxon by complex ', max= len(self.seq_dict)) as bar:
            for access_nb in self.access_dict:
                lineage = self.get_lineage(access_nb)
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
            self.__update_data()
            t1 = time.time()
            print("\n   ==> Taxon assembled by complex in %f seconds."%(t1 - t0))  

#########   EXPORT  ######### 

    def export_tax_id_txt(self, output_file_name, separator):
        """Export the Taxonomic Table for Qiime2. Col1 with accession number and 
        Col2 with Lineage.

        Args:
            output_file_name (string): path to output file
            separator (string): column separator 
        """        
        with open(output_file_name, "w") as taxon_file:
            for access_nb in self.access_dict:
                lineage = self.get_lineage(access_nb)
                if lineage == "NA": lineage = ""
                taxon_file.write(access_nb+separator+lineage+self.get_taxon(access_nb)+"\n")
        print("   ==> Taxonomic table successfully exported.")

    def export_seq_fasta(self, output_file_name, fasta_type):
        """Export the sequences FASTA file. Each sequence is written on a single line 
        if 'qiime' is specified or multiple lines (70bp) if 'ncbi' is specified.
        The description is the accession number of the sequence. 

        Args:
            output_file_name (string): path to output file
            fasta_type (string): 'qiime' or 'ncbi'
        """        
        if fasta_type == "qiime": self.__seq_qiime_fasta(output_file_name)
        elif fasta_type == "ncbi": self.__seq_ncbi_fasta(output_file_name)
        else: print("Please reload this function with 'qiime' or 'ncbi' argument.")
  
    def __seq_ncbi_fasta(self, output_file_name):
        """Called by export_seq_fasta if 'ncbi' is specified.
        This function makes the file.

        Args:
            output_file_name (string): path to output file
        """        
        with open(output_file_name, "w") as fasta_file:
            for seq_name in self.seq_dict:
                fasta_file.write(seq_name+"\n")
                seq_list = re.findall('.{1,70}', self.seq_dict[seq_name])
                for seq in seq_list:
                    fasta_file.write(seq+"\n")
        print("   ==> Sequences variants FASTA file susccessfully exported.")

    def __seq_qiime_fasta(self, output_file_name):
        """Called by export_seq_fasta if 'qiime' is specified.
        This function makes the file.

        Args:
            output_file_name (string): path to output file
        """        
        with open(output_file_name, "w") as fasta_file:
            for access_nb in self.access_dict:
                fasta_file.write(">"+access_nb+"\n"+self.get_sequence(access_nb)+"\n")
        print("   ==> Sequences variants FASTA file susccessfully exported.")
        
    def export_taxon_list(self, output_file_name):
        """Export the taxon list with their number of ASV. Col1 = taxon, col2 = nb of ASV. 

        Args:
            output_file_name (string): path to the output file
        """        
        taxon_list = []
        for taxon in self.taxon_dict:
            taxon_list.append(taxon) 
        utils.export_list_csv(sorted(taxon_list), output_file_name)
        print("   ==> Taxon list susccessfully exported.")

    def export_ampli_fasta(self, output_file_name, fasta_type):
        """Export the amplicon FASTA file. Each amplicon is written on a single line 
        if 'qiime' is specified or multiple lines (70bp) if 'ncbi' is specified.
        The description is the accession number of the amplicon. 

        Args:
            output_file_name (string): path to output file
            fasta_type (string): 'qiime' or 'ncbi'
        """        
        if fasta_type == "qiime": self.__ampli_qiime_fasta(output_file_name)
        elif fasta_type == "ncbi": self.__ampli_ncbi_fasta(output_file_name)
        else: print("Please reload this function with 'qiime' or 'ncbi' argument.")

    def __ampli_ncbi_fasta(self, output_file_name):
        """Called by export_ampli_fasta if 'ncbi' is specified.
        This function makes the file.

        Args:
            output_file_name (string): path to the output file
        """        
        with open(output_file_name, "w") as amplicon_file:
            for access_nb in self.access_dict:
                amplicon_file.write(">"+access_nb+"_"+self.get_taxon(access_nb)+"\n")
                seq_list = re.findall('.{1,70}', self.get_amplicon(access_nb))
                for seq in seq_list:
                    amplicon_file.write(seq+"\n")
        print("   ==> Amplicon FASTA file susccessfully exported.")

    def __ampli_qiime_fasta(self, output_file_name):
        """Called by export_ampli_fasta if 'ncbi' is specified.
        This function makes the file.

        Args:
            output_file_name (string): path to output file
        """        
        with open(output_file_name, "w") as amplicon_file:
            for access_nb in self.access_dict:
                amplicon_file.write(">"+access_nb+"\n"+self.get_amplicon(access_nb)+"\n")
        print("   ==> Amplicon FASTA file susccessfully exported.")

    def export_shared_ampli(self, output_file_name):
        """Export a file with all the accession number which have an identical amplicon.
        Each interaction is represented as 'accession_specie1 <=> accession_specie2'.
        The symbol '<=>' means there a common amplicon between the two sequences. 

        Args:
            output_file_name (string): path to the output file
        """        
        utils.export_list_csv(self.get_shared_amplicons(), output_file_name)
        print("   ==> Shared amplicons list successfully exported.")
    
    def export_complex_dict(self, outputfile_name, separator):
        """Export a file with all the complex line per line. On each line, 
        all the accession number related to this complex are written and also
        all the original species in the complex (all the species which have now the name 
        of the complex instead of there original name).
        Exemple: Fusarium Fujikori complex  {'access': [AJ567689, AK787990], 'taxon': [Fusarium_fujikoroi] }

        Args:
            outputfile_name (string): path to the output file
            separator (string): seprator between the complex name and the complex infos
        """        
        utils.export_dict_csv(self.complex_dict, outputfile_name, separator)
        print("   ==> Complex dictionnary successfully exported.")

    def export_modified_tax(self, output_file_name):
        """Export the list of all taxon which have been modified during the
        lineage construction. For example: 'Gibberella_sp. => Fusarium_oxysporum'.

        Args:
            output_file_name (string): path to output file
        """        
        utils.export_list_csv(self.__modified_taxon, output_file_name)
        print("   ==> Modified taxon list successfully exported.")

    def export_homologous_dict(self, output_file_name):
        """WARNING : This function is not well working. This require a lot ot time to 
        to run and require a very stable internet connexion to avoid any trouble. 

        Export a list of accession number from NCBI 'nt' which have 100% of similarity with
        the final amplicons in our database and are not the same species. 
        This function can be use for database validation. 

        Args:
            output_file_name (string): path to output file 
        """        
        t0 = time.time()
        sequences = random.sample(list(self.access_dict), 30)
        print(sequences)
        homologous = {"GENUS":"ACCESSESSION: HomologousAccession_HomologousGenus"}
        with FillingSquaresBar('Getting Blast result from NCBI servers ', max= len(sequences)) as bar: 
            for access_nb in sequences:
                print(self.access_dict[access_nb]["amplicon"])
                result_handle = NCBIWWW.qblast("blastn", "nt", self.access_dict[access_nb]["amplicon"])
                blast_record = NCBIXML.read(result_handle)
                genus = self.get_genus(access_nb)
                homologous[genus] = {}
                for blast in blast_record.descriptions:
                    if genus not in set(blast.title.split()):
                        try: homologous[genus][access_nb].append(blast.accession+"_"+blast.title.split()[1])
                        except:
                            homologous[genus][access_nb] = blast.accession+"_"+blast.title.split()[1]
                bar.next
        t1 = time.time()
        print("\n   ==> Blast the sequences in %f seconds."%(t1 - t0))
        utils.export_dict_csv(homologous, output_file_name, "  ==>  ")
        print("   ==> Interspecies out of selection amplicons dict successfully exported.")

    def export_access_dict(self, output_file_name):
        """Export the access_dict with all the informations of the data to JSON format.
        WARNING : You must write the file extension as '.json' to avoir Errors.

        Args:
            output_file_name (string): path to output file
        """        
        utils.save_as_json(self.access_dict, output_file_name)
        print("   ==> Access_dict successfully exported.")