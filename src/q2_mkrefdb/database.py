import re, time
from . import utils, filter, lineage, stats, edit, export
from . import merge as mg
from .pcr import amplify
from progress.bar import FillingSquaresBar
from progress.spinner import PixelSpinner
from Bio import Entrez

#######################################################################################################################
#################     CLASS Database     ##############################################################################
#######################################################################################################################

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
            origin_database (string): 'ncbi', 'ebi', 'ddbj' or rnaCentral
            forward_primer_fasta (string): path to the forward primer FASTA file
            reverse_primer_fasta (string): path to the reverse primer FASTA file
            fw_mismatch_tol (int, optional): forward primer mismatch tolerance for annealing. Defaults to 3.
            rv_mismatch_tol (int, optional): reverse primer mismatch tolerance for annealing. Defaults to 3.
            trim_primers = True to trim primers on amplicons, else False.
        """              
        self.file_name = fasta_file.rstrip(".fasta")
        self.origin_database = origin_database
        self.fw_mismatch_tol = fw_mismatch_tol
        self.rv_mismatch_tol = rv_mismatch_tol
        self.trim_primers = trim_primers
        self.forward_primer_list = self.__get_seq_from_fasta(forward_primer_fasta, False)
        self.reverse_primer_list = self.__get_seq_from_fasta(reverse_primer_fasta, False)
        self.seq_dict = self.__get_seq_from_fasta(fasta_file, True)
        if origin_database in {"silva", "unite"}: self.access_dict = self.__make_access_dict2()
        else: self.access_dict = self.__make_access_dict()
        self.taxon_dict = self.__make_taxon_dict()
        self.modified_taxon = set()
        self.__make_amplicon_dict()
        self.__make_lineage_dict()
        self.complex_dict = self.__make_complex_dict()
        self.sa_threshold = 1000
        return

    def import_db(self, database_json):
        db = utils.open_from_json(database_json)
        self.file_name = db["infos"]["file_name"]
        self.origin_database = db["infos"]["origin_db"]
        self.fw_mismatch_tol = db["infos"]["fw_mis_tol"]
        self.rv_mismatch_tol = db["infos"]["rv_mis_tol"]
        self.trim_primers = db["infos"]["trim_prim"]
        self.forward_primer_list = db["infos"]["fw_prim"]
        self.reverse_primer_list = db["infos"]["rv_prim"]
        self.sa_threshold = db["infos"]["sa_threshold"]
        db.pop("infos")
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
                    try:
                        if line.split(" ")[1] not in {"UNVERIFIED:", "Uncultured"}:
                            seq_dict[line] = ""
                            seq_ref = line
                        else: seq_ref = "pass"
                    except IndexError:
                        seq_dict[line] = ""
                        seq_ref = line
                elif (line[0] == ">") & (bool_multiple_seq == False):
                    seq_list.append(line)
                    seq_ref = line
                elif seq_ref != "pass":
                    if bool_multiple_seq == True: seq_dict[seq_ref] += line.replace("U", "T").upper()
                    else: seq_list.append(line.replace("U", "T").upper())
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
                access_dict[self.get_access_from_des(seq_name)] = \
                    {"name":seq_name, "sequence":self.seq_dict[seq_name], \
                        "taxon":seq_name.split()[1]+"_"+seq_name.split()[2]}
            else: new_seq_dict.pop(seq_name)
        self.seq_dict = new_seq_dict    
        return access_dict

    def __make_access_dict2(self):
        access_dict = {}
        new_seq_dict = dict(self.seq_dict)
        for seq_name in self.seq_dict:
            taxon =""
            if self.origin_database == "silva":
                taxon = "_".join(seq_name.split(";")[-1].split()[:2])
                if (taxon not in {"unidentified", "metagenome", "uncultured"}) \
                    & taxon[0].islower():
                    taxon = seq_name.split(";")[-2]+"_"+taxon.split("_")[0]
            if self.origin_database == "unite":
                taxon = "_".join(seq_name.split("|")[1].split(";")[-1].split("_")[2:4])
            if taxon.split("_")[0] not in {"unidentified", "metagenome", "uncultured"}:
                access_dict[self.get_access_from_des(seq_name)] = \
                    {"name":seq_name, "sequence":self.seq_dict[seq_name], "taxon":taxon}
            else: new_seq_dict.pop(seq_name)
        self.seq_dict = new_seq_dict
        return access_dict

    def update_data(self):
        """Update the information of the class.
        """        
        new_access_dict = {}
        for seq_name in self.seq_dict:
            access_nb = self.get_access_from_des(seq_name)
            new_access_dict[access_nb] = self.access_dict[access_nb]
        self.access_dict = new_access_dict
        self.taxon_dict = self.__make_taxon_dict()
        self.complex_dict = self.__make_complex_dict()

    def update_data2(self):
        """Update the information of the class.
        """
        self.taxon_dict = self.__make_taxon_dict()
        self.complex_dict = self.__make_complex_dict()

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
                for taxon in self.taxon_dict:
                    name_list = taxon.split("_")
                    self.__get_lineage_norm(taxon, name_list)
                    bar.next()
                if len(self.__na_tax_str) > 0: print("\nThe following taxons have NO lineage : \n", self.__na_tax_str, "\n")
                self.update_data()
                t1 = time.time()
                print("\n   ==> Got species lineages in %f seconds."%(t1 - t0))
    
    # def __make_lineage_dict2(self):
    #     with FillingSquaresBar('Getting species lineages from EBI taxonomy online ressources ', \
    #         max= len(self.taxon_dict)) as bar:
    #             t0 = time.time()  
    #             for taxon in self.taxon_dict:
    #                 for access_nb in self.taxon_dict[taxon]:
    #                     if self.origin_database == "silva":
    #                         lineage = self.get_name(access_nb).split()[1].split(";")[:-1]
    #                     if self.origin_database == "unite":
    #                         lineage = self.get_name(access_nb).split("|")[1].split(";")[:-1]
    #                     self.access_dict[access_nb]["lineage"] = "; ".join(lineage)+"; "
    #                 bar.next()
    #             self.update_data()
    #             t1 = time.time()
    #             print("\n   ==> Got species lineages in %f seconds."%(t1 - t0))
 
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
        if (taxon[-1] == ".") & (lineage == "NA"): # If taxon end by a dot, then we look for lineage of each accession nb
            n = 0
            for access_nb in self.taxon_dict[taxon]:
                print("OK --> ", access_nb)
                n += 1
                if self.origin_database == "rnaCentral":
                    lineage = self.__complete_lineage2(access_nb)
                elif access_nb[:3] != "UDB": lineage = self.__complete_lineage(access_nb)
                else: lineage = "NA"
                if str(taxon+" => "+self.get_taxon(access_nb)) not in self.modified_taxon:
                    self.modified_taxon.add(taxon+" => "+self.get_taxon(access_nb))
                if lineage == "NA": self.__na_tax_str += (taxon+"\n")
                self.access_dict[access_nb]["lineage"] = lineage    
        else:
            for access_nb in self.taxon_dict[taxon]:
                print(" . or NA --> ", access_nb)
                if lineage == "NA":
                    if self.origin_database == "rnaCentral":
                        lineage = self.__complete_lineage2(access_nb)
                    elif access_nb[:3] != "UDB": lineage = self.__complete_lineage(access_nb)
                    else: lineage = "NA"
                    if str(taxon+" => "+self.get_taxon(access_nb)) not in self.modified_taxon:
                        self.modified_taxon.add(taxon+" => "+self.get_taxon(access_nb))
                if lineage == "NA": self.__na_tax_str += (taxon+"\n")
                self.access_dict[access_nb]["lineage"] = lineage

    def __complete_lineage(self, access_nb):
        """Third part of the _make_lineage_dict function, used if the specie is not found on the EBI.
        Research the specie on the NCBI database  with Entrez (longer than EBI) then search again the
        lineage on EBI. 

        Args:
            access_nb (string): Sequence accession number related to the specie.
            url (string): url to EBI website

        Returns:
            string: lineage if found, else 'NA'
        """
        xml_url = "https://www.ebi.ac.uk/ena/browser/api/text/%s?lineLimit=25"%access_nb
        try: 
            taxon = ""
            lineage = ""
            response = utils.get_var_from_url(xml_url, "str")
            for line in response.split("\n"):
                if line[:2] == "OS": taxon = "_".join(line.split()[1:3])
                if line[:2] == "OC": lineage += " ".join(line.split()[1:])+" "
            self.access_dict[access_nb]["taxon"] = taxon
            return lineage.rstrip(" .")+"; "
        except UnicodeDecodeError:
            Entrez.email = 'someuser@mail.com'
            handle = Entrez.efetch(db="nucleotide", id=access_nb, rettype="gb", retmode="text")
            result=handle.read().split('\n')
            for line in result:
                if 'ORGANISM' in line :
                    self.access_dict[access_nb]["taxon"] = "_".join(line.split()[1:3])
                    url = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"+"%20".join(line.split()[1:3])
                    try: return utils.get_var_from_url(url, "json")[0]['lineage']
                    except: return "NA" 

    def __complete_lineage2(self, access_nb):
        xml_url = "https://www.ebi.ac.uk/ena/browser/api/xml/"+access_nb.split("_")[1]
        response = utils.get_var_from_url(xml_url, "str")
        start = response.find("taxon scientificName=") + len("taxon scientificName=")
        end = response.find("<lineage>")
        substring = response[start:end]
        taxon = substring.split('"')[1].split()
        if "UTF-8" not in taxon: 
            self.access_dict[access_nb]["taxon"] = "_".join(taxon)
            url = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"+"%20".join(taxon)
            try:
                return utils.get_var_from_url(url, "json")[0]['lineage']
            except: 
                return "NA"
        else: return "NA"

    def __make_complex_dict(self):
        """Create the complex dictionnary from the lineage informations.

        Returns:
            dictionnary: complex dictionnary
        """        
        complex_dict = {}
        for access_nb in self.access_dict:
            if self.get_lineage(access_nb) not in {"NA", None, "; "}:
                complex_name = self.get_lineage(access_nb).split("; ")[-2]
                if self.get_taxon(access_nb).split("_")[-1] in {"complex", "group", "clade"}:
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
        self.access_dict = lineage.custom_group(self, groups_file)
        self.update_data()

    def get_access_from_des(self, seq_name):
        """Get the sequence accession number from the seqeunce name.

        Args:
            seq_name (string): raw sequence name from seq_dict

        Returns:
            string: accession number
        """        
        access_num = ""
        if self.origin_database == "ncbi":
            access_num = seq_name.split(" ")[0].strip(">")[:-2]
        if self.origin_database == "ebi":
            access_num = seq_name.split("|")[1]
        if self.origin_database == "ddbj":
            access_num = seq_name.split("|")[0].strip(">")
        if self.origin_database == "rnaCentral":
            access_num = seq_name.split()[0].strip(">")
        if self.origin_database == "unite":
            access_num = seq_name.split("|")[0].strip(">")
        if self.origin_database == "silva":
            access_num = seq_name.split(".")[0].strip(">")
        return access_num

    def __make_amplicon_dict(self):
        """Add the amplicons keys to the access_dict. Each amplicon is created thanks to the pcr module.
        The default values of mismatch tolerance is 3. 
        """        
        spinner = PixelSpinner("Loading data ")
        for seq_name, sequence in self.seq_dict.items():
            try:
                amplicon_seq = amplify(sequence, self.forward_primer_list[1], self.reverse_primer_list[1], \
                    self.fw_mismatch_tol, self.rv_mismatch_tol, self.trim_primers)
                self.access_dict[self.get_access_from_des(seq_name)]["amplicon"] = amplicon_seq
            except ValueError:
                self.access_dict[self.get_access_from_des(seq_name)]["amplicon"] = "NA"
            spinner.next()

#######################################################################################################################
#####################        INFOS       ##############################################################################
#######################################################################################################################

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
        "S Number of sequences : "+str(len(self.seq_dict))+"\n"+ \
        "S Maximum sequences length : "+str(stats.max_seq_len(self))+"\n"+\
        "S Minimum sequences length : "+str(stats.min_seq_len(self))+"\n"+\
        "S Mean sequences length : "+str(stats.mean_seq_len(self))+"\n\n"+\
        "---------------TAXONS-----------------------\n"+\
        "T Number of different taxons : "+str(len(self.taxon_dict))+"\n"+\
        "T Maximum number of sequences per taxon : "+str(stats.max_nb_seq_tax(self))+"\n"+\
        "T Minimum number of sequences per taxon : "+str(stats.min_nb_seq_tax(self))+"\n"+\
        "T Mean number of sequences per taxon : "+str(stats.mean_nb_seq_tax(self))+"\n"+\
        "T Number of SAs taxons : "+str(len(self.get_SAs()))+"\n\n"+\
        "---------------COMPLEXS---------------------\n"+\
        "C Number of differents complex : "+str(len(self.complex_dict))+"\n"+\
        "C Maximum number of taxons per complex : "+str(stats.max_tax_complex(self))+"\n"+\
        "C Minimum number of taxons per complex : "+str(stats.min_tax_complex(self))+"\n"+\
        "C Mean number of taxons per complex : "+str(stats.mean_tax_complex(self))+"\n"+\
        "C Maximum number of sequences per complex : "+str(stats.max_seq_complex(self))+"\n"+\
        "C Minimum number of sequences per complex : "+str(stats.min_seq_complex(self))+"\n"+\
        "C Mean number of sequences per complex : "+str(stats.mean_seq_complex(self))+"\n\n"+\
        "--------------AMPLICONS---------------------\n"+\
        "A Number of amplicons : "+str(stats.get_amplicon_nb(self))+"\n"+\
        "A Maximum amplicons length : "+str(stats.max_amplicon_len(self))+"\n"+\
        "A Minimum amplicons length : "+str(stats.min_amplicon_len(self))+"\n"+\
        "A Mean amplicons length : "+str(stats.mean_amplicon_len(self))+"\n"+\
        "A Number of redundant amplicons : "+str(len(self.get_shared_amplicons(self.sa_threshold)[0]))+"\n"+\
        "****************************************\n"
        return message

    def get_infos_dict(self):
        return {
            "nb seq" : len(self.seq_dict),
            "max seq len" : stats.max_seq_len(self),
            "min seq len" : stats.min_seq_len(self),
            "mean seq len" : stats.mean_seq_len(self),
            "nb tax" : len(self.taxon_dict),
            "max seq nb / tax" : stats.max_nb_seq_tax(self),
            "min seq nb / tax": stats.min_nb_seq_tax(self),
            "mean seq nb / tax" : stats.mean_nb_seq_tax(self),
            "nb SA" : len(self.get_SAs()),
            "nb complex" : len(self.complex_dict),
            "max tax nb / complex" : stats.max_tax_complex(self),
            "min tax nb /complex" : stats.min_tax_complex(self),
            "mean tax nb /complex" : stats.mean_tax_complex(self),
            "max seq nb / complex" : stats.max_seq_complex(self),
            "min seq nb /complex" : stats.min_seq_complex(self),
            "mean seq nb / complex" : stats.mean_seq_complex(self),
            "nb ampl" : stats.get_amplicon_nb(self),
            "max ampl len" : stats.max_amplicon_len(self),
            "min ampl len" : stats.min_amplicon_len(self),
            "mean ampl len" : stats.mean_amplicon_len(self),
            "nb redund ampl" : len(self.get_shared_amplicons(self.sa_threshold)[0]),
        }

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
    
    def get_SAs(self):
        """Get the set of all SA taxons.

        Returns:
            set: all SA taxons
        """
        SAs = set()
        for taxon in self.taxon_dict:
            if taxon.split("_")[-1][:2] == "SA":
                SAs.add(taxon)
        return SAs

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

    def get_shared_amplicons(self, threshold):
        """Get a list of all the sequences which have an identical amplicon.

        Args:
            threshold (int): Threshold number od sequences with the same amplicon.

        Returns:
            list: list of strings which describes sequences with shared amplicon
        """        
        shared_ampli_list = []
        checked_sequences = {}
        sa_number = 0
        sa_dict = {}
        self.sa_threshold = threshold
        for access_nb in self.access_dict:
            amplicon = self.get_amplicon(access_nb)
            if  amplicon not in checked_sequences:
                checked_sequences[amplicon] = {access_nb+"|"+self.get_taxon(access_nb)}
            else:
                checked_sequences[amplicon].add(access_nb+"|"+self.get_taxon(access_nb))
        for amplicon in checked_sequences:
            if threshold > len(checked_sequences[amplicon]) > 1:
                shared_amplicon = ""
                for access_nb in checked_sequences[amplicon]:
                    shared_amplicon += access_nb+"  <--->  "
                shared_ampli_list.append(shared_amplicon)
            elif len(checked_sequences[amplicon]) > threshold:
                sa_number += 1 
                sa_id = "SA"+str(sa_number)
                shared_ampli_list.append(sa_id)
                sa_dict[sa_id] = checked_sequences[amplicon]
        return shared_ampli_list, sa_dict

    def get_lineage(self, access_number):
        """Get the lineage related to the accession number.

        Args:
            access_number (string): accession number

        Returns:
            string: lineage
        """        
        return self.access_dict[access_number]["lineage"]

#######################################################################################################################
#################         FILTER         ##############################################################################
#######################################################################################################################

    def clean_dataset(self, wanted_genus1, wanted_genus2, unverified_bool):
        self.seq_dict = filter.clean_dataset(self, wanted_genus1, wanted_genus2, unverified_bool)
        self.update_data()
    clean_dataset.__doc__ = filter.clean_dataset.__doc__

    def del_redund_seq_glob(self):
        self.seq_dict = filter.del_redund_seq_glob(self)
        self.update_data()
    del_redund_seq_glob.__doc__ = filter.del_redund_seq_glob.__doc__

    def del_redund_seq_tax(self):
        self.seq_dict = filter.del_redund_seq_tax(self)
        self.update_data()
    del_redund_seq_tax.__doc__ = filter.del_redund_seq_tax.__doc__

    def del_na_amplicons(self):
        self.seq_dict = filter.del_na_amplicons(self)
        self.update_data()
    del_na_amplicons.__doc__ = filter.del_na_amplicons.__doc__

    def del_redund_ampli_glob(self):
        self.seq_dict = filter.del_redund_ampli_glob(self)
        self.update_data()
    del_redund_ampli_glob.__doc__ = filter.del_redund_ampli_glob.__doc__

    def del_redund_ampli_tax(self):
        self.seq_dict = filter.del_redund_ampli_tax(self)
        self.update_data()
    del_redund_ampli_tax.__doc__ = filter.del_redund_ampli_tax.__doc__

    def group_by_complex(self):
        self.access_dict = filter.group_by_complex(self)
        self.update_data()
    group_by_complex.__doc__ = filter.group_by_complex.__doc__

#######################################################################################################################
#################         EDIT           ##############################################################################
#######################################################################################################################

    def remove_by_id(self, id_list_csv):        
        self.seq_dict = edit.remove_by_id(self, id_list_csv)
        self.update_data()
    remove_by_id.__doc__ = edit.remove_by_id.__doc__

    def remove_by_taxon(self, tax_list_csv):
        self.seq_dict = edit.remove_by_taxon(self, tax_list_csv)
        self.update_data()
    remove_by_taxon.__doc__ = edit.remove_by_taxon.__doc__

    def rename_by_id(self, id_list_csv):
        self.access_dict = edit.rename_by_id(self, id_list_csv)
        self.update_data()
    rename_by_id.__doc__ = edit.rename_by_id.__doc__

    def group_by_id(self, shared_ext_csv):
        group = edit.group_by_id(self, shared_ext_csv)
        self.access_dict = group[0]
        self.seq_dict = group[1]
        self.update_data()
    group_by_id.__doc__ = edit.group_by_id.__doc__

#######################################################################################################################
#################         MERGE          ##############################################################################
#######################################################################################################################

    def merge(self, Database2, sa_file1, sa_file2, output_Database):
        mg.merge(self,Database2, sa_file1, sa_file2, output_Database)
    merge.__doc__ =  mg.merge.__doc__

    def merge_all(dir_path, s_or_g, outputDB):
        mg.merge_all(dir_path, s_or_g, outputDB)
    merge_all.__doc__ =  mg.merge_all.__doc__

#######################################################################################################################
#################         EXPORT         ##############################################################################
#######################################################################################################################

    def export_tax_id_txt(self, output_file_name, separator):
        export.export_tax_id_txt(self, output_file_name, separator)
    export_tax_id_txt.__doc__ = export.export_tax_id_txt.__doc__

    def export_seq_fasta(self, output_file_name, fasta_type):
        export.export_seq_fasta(self, output_file_name, fasta_type)
    export_seq_fasta.__doc__ = export.export_seq_fasta.__doc__

    def export_taxon_list(self, output_file_name):
        export.export_taxon_list(self, output_file_name)
    export_taxon_list.__doc__ = export.export_taxon_list.__doc__

    def export_ampli_fasta(self, output_file_name, fasta_type):
        export.export_ampli_fasta(self, output_file_name, fasta_type)
    export_ampli_fasta.__doc__ = export.export_ampli_fasta.__doc__

    def export_shared_ampli(self, output_file_name, threshold):
        export.export_shared_ampli(self, output_file_name, threshold)
    export_shared_ampli.__doc__ = export.export_shared_ampli.__doc__
    
    def export_complex_dict(self, output_file_name, separator):
        export.export_complex_dict(self, output_file_name, separator)
    export_complex_dict.__doc__ = export.export_complex_dict.__doc__

    def export_modified_tax(self, output_file_name):
        export.export_modified_tax(self, output_file_name)
    export_modified_tax.__doc__ = export.export_modified_tax.__doc__

    def export_access_dict(self, output_file_name):
        export.export_access_dict(self, output_file_name)
    export_access_dict.__doc__ = export.export_access_dict.__doc__
    
    def export_phylo_annot(self, output_file):
        export.phylo_annotation(self, output_file)
    export_phylo_annot.__doc__ = export.phylo_annotation.__doc__