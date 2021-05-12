import re
from . import utils

def phylo_annotation(Database, output_file):
    """Export an annotation file that can be used in iTOL tool at
    https://itol.embl.de/ after creating the tree with https://ngphylogeny.fr/.

    Args:
        Database (Database): the Database object
        output_file (file.txt): the annotation file for iTOL
    """    
    with open(output_file, "w") as file:
        file.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for taxon in Database.taxon_dict:
            if taxon.split("_")[1][0] != "S": 
                color = utils.random_color()
                for access_nb in Database.taxon_dict[taxon]:
                    file.write(access_nb+"_"+taxon+"\tlabel\trgb"+str(color)+"\tnormal\t1\n")
                    file.write(access_nb+"_"+taxon+"\tbranch\trgb"+str(color)+"\tnormal\t1\n")

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

def export_taxon_list(self, output_file_name):
        """Export the taxon list with their number of ASV. Col1 = taxon, col2 = nb of ASV. 

        Args:
            output_file_name (string): path to the output file
        """        
        taxon_list = []
        for taxon in self.taxon_dict:
            taxon_list.append(taxon+"\t"+str(len(self.taxon_dict[taxon]))) 
        utils.export_list_csv(sorted(taxon_list), output_file_name)
        print("   ==> Taxon list susccessfully exported.")

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

def __seq_phylo_fasta(self, output_file_name):
    """Called by export_seq_fasta if 'phylo' is specified.
    This function makes the file.
    Args:
        output_file_name (string): path to output file
    """        
    with open(output_file_name, "w") as fasta_file:
        for access_nb in self.access_dict:
            fasta_file.write(">"+access_nb+"_"+self.get_taxon(access_nb)+"\n"+\
                self.get_sequence(access_nb)+"\n")
    print("   ==> Sequences variants FASTA file susccessfully exported.")

def export_seq_fasta(self, output_file_name, fasta_type):
    """Export the sequences FASTA file. Each sequence is written on a single line 
    if 'qiime' is specified or multiple lines (70bp) if 'ncbi' is specified.
    The description is the accession number of the sequence.

    Args:
        output_file_name (string): path to output file
        fasta_type (string): 'qiime' or 'ncbi'
    """        
    if fasta_type == "qiime": __seq_qiime_fasta(output_file_name)
    elif fasta_type == "ncbi": __seq_ncbi_fasta(output_file_name)
    elif fasta_type == "phylo": __seq_phylo_fasta(output_file_name)
    else: print("Please reload this function with 'qiime', 'ncbi' or 'phylo' argument.")

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


def __ampli_phylo_fasta(self, output_file_name):
    """Called by export_ampli_fasta if 'phylo' is specified.
    This function makes the file.
    
    Args:
        output_file_name (string): path to output file
    """        
    with open(output_file_name, "w") as amplicon_file:
        for access_nb in self.access_dict:
            amplicon_file.write(">"+access_nb+"_"+self.get_taxon(access_nb)+"\n"+\
                self.get_amplicon(access_nb)+"\n")
    print("   ==> Amplicon FASTA file susccessfully exported.")

def export_ampli_fasta(self, output_file_name, fasta_type):
    """Export the amplicon FASTA file. Each amplicon is written on a single line 
    if 'qiime' is specified or multiple lines (70bp) if 'ncbi' is specified.
    The description is the accession number of the amplicon.

    Args:
        output_file_name (string): path to output file
        fasta_type (string): 'qiime' or 'ncbi'
    """        
    if fasta_type == "qiime": __ampli_qiime_fasta(self, output_file_name)
    elif fasta_type == "ncbi": __ampli_ncbi_fasta(self, output_file_name)
    elif fasta_type == "phylo": __ampli_phylo_fasta(self, output_file_name)
    else: print("Please reload this function with 'qiime' or 'ncbi' argument.")

def export_shared_ampli(self, output_file_name, threshold):
        """Export a file with all the accession number which have an identical amplicon.
        Each interaction is represented as 'accession_specie1 <=> accession_specie2'.
        The symbol '<=>' means there a common amplicon between the two sequences. 

        Args:
            output_file_name (string): path to the output file
        """
        sa_ampli = self.get_shared_amplicons(threshold)
        utils.export_list_csv(sa_ampli[0], output_file_name)
        utils.export_dict_csv2(sa_ampli[1], "_ext.".join(output_file_name.rsplit(".", 1)), "\t")        
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
        utils.export_list_csv(self.modified_taxon, output_file_name)
        print("   ==> Modified taxon list successfully exported.")

def export_access_dict(self, output_file_name):
        """Export the access_dict with all the informations of the data to JSON format.
        WARNING : You must write the file extension as '.json' to avoir Errors.

        Args:
            output_file_name (string): path to output file
        """
        self.access_dict["infos"] = {
            "file_name":self.file_name,
            "origin_db":self.origin_database,
            "fw_mis_tol":self.fw_mismatch_tol,
            "rv_mis_tol":self.rv_mismatch_tol,
            "trim_prim":self.trim_primers,
            "fw_prim":self.forward_primer_list,
            "rv_prim":self.reverse_primer_list,
            "sa_threshold":self.sa_threshold}     
        utils.save_as_json(self.access_dict, output_file_name)
        print("   ==> Access_dict successfully exported.")