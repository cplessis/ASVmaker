import argparse, sys, getpass, datetime, platform
from . import fastafilter as ff
from . import utils

__doc__ =  """This package is released and maintained by the Research and Development Institute 
for Agri-Environnement (IRDA - Quebec). More informations at https://github.com/cplessis/q2_mkrefdb."""

def get_arguments():
    parser = argparse.ArgumentParser(description=__doc__)  
    #---------------------------------------------------
    #                  Init
    #---------------------------------------------------
    parser.add_argument('-i',
                          '--sequences_input',
                          help="File path to the raw sequences FASTA file to treat.",
                          required=True)
    parser.add_argument('-db',
                          '--source_database',
                          choices=['ncbi', 'ebi', 'ddbj'],
                          help="Name of the database from which the FASTA file was obtained.",
                          required=True)                          
    parser.add_argument('-fp',
                          '--forward_primer',
                          help="File path to the forward primer FASTA file.",
                          required=True)
    parser.add_argument('-fmt',
                          '--fw_mismatch_tol',
                          help="Number of mismatch which are accepted for the forward primer annealing.",
                          default=3)
    parser.add_argument('-rp',
                          '--reverse_primer',
                          help="File path to the reverse primer FASTA file.",
                          required=True)
    parser.add_argument('-rmt',
                          '--rv_mismatch_tol',
                          help="Number of mismatch which are accepted for the reverse primer annealing.",
                          default=3)
    parser.add_argument('-inf',
                          '--infos_file',
                          help="Informations file in which all the files treatments are recorded.",
                          default="none")
    parser.add_argument('-dit',
                          '--displ_inf_terminal',
                          help="Display information file in terminal if arg is specified.",
                          action = 'store_true',
                          default=False)
    
    #---------------------------------------------------
    #                Cleaning
    #---------------------------------------------------    
    # By default the cleaning does not occur if args are not specified. 
    # This means all the genus in the fasta file are conserved. 
    parser.add_argument('-g1',
                          '--genus1',
                          help="Name of the first genus to keep.\
                          Ex: Fusarium",
                          default="none")
    parser.add_argument('-g2',
                          '--genus2',
                          help="Name of the second genus to keep.\
                          Ex: Gibberella",
                          default="none")
    parser.add_argument('-u',
                          '--unverified',
                          action = 'store_true',
                          help="Keep the unverified genus if arg is specified.",
                          default=False)
    parser.add_argument('-tbc',
                          '--taxonomy_by_complex',
                          help="Do not group species with same complex lineage if arg is specified.",
                          action = 'store_false',
                          default=True)

    #---------------------------------------------------
    #                Filtering
    #---------------------------------------------------    
    parser.add_argument('-f',
                          '--filtering_type',
                          choices=['taxon', 'global'],
                          help="",
                          default='none')
    parser.add_argument('-rsv',
                          '--redund_seq_variants',
                          action = 'store_false',
                          help="Keep the redundant sequences variants if arg is specified.",
                          default=True)
    parser.add_argument('-ra',
                          '--redundant_amplicon',
                          action = 'store_false',
                          help="Keep the redundant amplicons if arg is specified.",
                          default=True)
    parser.add_argument('-nasv',
                          '--seq_without_ampl',
                          action = 'store_false',
                          help="Keep the sequences without amplicons if arg is specified.",
                          default=True)
    
    #---------------------------------------------------
    #                Export
    #---------------------------------------------------    
    parser.add_argument('-svo',
                          '--seq_variants_output',
                          help="File name of the sequences variants FASTA output.",
                          default=None)
    parser.add_argument('-svop',
                          '--seq_output_phylo',
                          help="File name of the sequences variants FASTA output in NCBI format (for phylogeny).",
                          default=None)
    parser.add_argument('-ao',
                          '--amplicons_output',
                          help="File name of the amplicons FASTA output.",
                          default=None)
    parser.add_argument('-aop',
                          '--ampl_output_phylo',
                          help="File name of the amplicon variants FASTA output in NCBI format (for phylogeny).",
                          default=None)
    parser.add_argument('-to',
                          '--taxonomy_output',
                          help="File name of the taxonomy TXT output.",
                          default=None)
    parser.add_argument('-tlo',
                          '--taxon_list_output',
                          help="File name of the taxon list TXT output.",
                          default=None)
    parser.add_argument('-sao',
                          '--shared_ampl_output',
                          help="File name of the shared amplicons TXT output.",
                          default=None)
    parser.add_argument('-cdo',
                          '--complex_dict_output',
                          help="File name of the complex dictionnary TXT output.",
                          default=None)
    parser.add_argument('-mto',
                          '--modified_tax_output',
                          help="File name of the modified taxon list TXT output.",
                          default=None)
    parser.add_argument('-ado',
                          '--access_dict_output',
                          help="File name of the access dictionnary JSON output.",
                          default=None)
    args = parser.parse_args()
    return args

args = get_arguments()

# INIT
sequences_input = args.sequences_input
source_database = args.source_database
forward_primer = args.forward_primer
reverse_primer = args.reverse_primer
infos_file = args.infos_file
displ_inf_terminal = args.displ_inf_terminal
fw_mismatch_tol = args.fw_mismatch_tol
rv_mismatch_tol = args.rv_mismatch_tol

print("\n\n\n\
=============================================================================================\n\
=========                               Q2_MkRefDb                                  =========\n\
=============================================================================================\n")

# Init input data: OUTPUT FILE INFO
if infos_file == "none":
    output_saver = utils.Logger("temp_output_saver.txt", displ_inf_terminal)
else:
    output_saver = utils.Logger(infos_file, displ_inf_terminal)

# Init input data: OUTPUT FILE INFO
output_saver.write("USERNAME : "+getpass.getuser()+"\n"+\
    "DATE  : "+str(datetime.date.today().strftime("%d/%m/%Y"))+"\n"+\
        "TIME  : "+str(datetime.datetime.now().strftime("%H:%M"))+"\n"+\
            "SYSTEM  : "+str(platform.system())+"\n"+\
                "PLATEFORM  : "+str(platform.platform())+"\n"+\
                    "PYTHON VERSION : "+str(sys.version)+"\n")

output_saver.write("\n\n\n\
=============================================================================================\n\
=========                               Init Data                                   =========\n\
=============================================================================================\n")

# Init input data: FASTA FILE from NCBI, EBI or DDBJ
output_saver.write("Fasta file path : "+sequences_input+"\n")
output_saver.write("Original database of the fasta file : "+source_database+"\n")
output_saver.write("Forward primer path : "+forward_primer+"\n")
output_saver.write("Reverse primer path : "+reverse_primer+"\n")

# Init input data: DATABASE creation
data = ff.Database(sequences_input, source_database, forward_primer, reverse_primer, \
    fw_mismatch_tol, rv_mismatch_tol)
output_saver.write(data.get_info("init data")+"\n")

output_saver.write("\n\n\n\
=============================================================================================\n\
=========                  PART 1/3                Cleaning                         =========\n\
=============================================================================================\n")
# CLEAN
genus1 = args.genus1
genus2 = args.genus2
unverified = args.unverified
taxonomy_by_complex = args.taxonomy_by_complex
seq_without_ampl = args.seq_without_ampl

if genus1 != "none":
    output_saver.write("* You have decided to 'CLEAN UP' your dataset.\n")
    if taxonomy_by_complex == True:
        data.group_by_complex()
        output_saver.write("* You have decided to 'GROUP TAXON by COMPLEX' in your dataset.\n")
        output_saver.write(data.get_info("After grouping by taxa \n"))
        output_saver.write("\n\
=============================================================================================\n")
    genus_set = data.get_all_genus()
    if genus1 in genus_set:
        output_saver.write(genus1+" is in your dataset."+"\n")
    else:
        output_saver.write("WARNING : There is NO genus named "+genus1+" in your dataset."+"\n")
    if (genus2 != "none") & (genus2 in genus_set):
        output_saver.write(genus2+" is in your dataset."+"\n")
    elif genus2 != "none":
        output_saver.write("WARNING : There is NO genus named "+genus2+" in your dataset."+"\n")
    data.clean_dataset(genus1, genus2, unverified)
    if unverified == True: 
        output_saver.write("* You have decided to 'KEEP' the 'UNVERIFIED' and 'UNCULTURED' genus.\n")
    else: 
        output_saver.write("* 'UNVERIFIED' and 'UNCULTURED' genus have been removed.\n")
    output_saver.write(data.get_info("After Cleaning with : "+genus1+" -- "+genus2+" -- UNVERIFIED="+str(unverified))+"\n")
    if seq_without_ampl == True:
        output_saver.write("\n\
=============================================================================================\n")
        data.del_na_amplicons()
        output_saver.write("** 'SEQUENCES VARIANTS' have been removed if 'NO AMPLICONS'.\n")
        output_saver.write(data.get_info("after deleting sequences variants without amplicons")+"\n")
    else:
        output_saver.write("** You decided to 'KEEP' the 'SEQUENCES VARIANTS' with 'NO AMPLICONS'.\n")
else:
    output_saver.write("* You have decided 'NOT' to 'CLEAN UP' your dataset.\n")

output_saver.write("\n\n\n\
=============================================================================================\n\
=========                  PART 2/3                Filtering                        =========\n\
=============================================================================================\n")

# FILTER
filtering_type = args.filtering_type
redund_seq_variants = args.redund_seq_variants
redundant_amplicon = args.redundant_amplicon

if filtering_type != "none":
    output_saver.write("* You have decided to 'FILTER' your dataset.\n")
    if redund_seq_variants == True:
        if filtering_type == "taxon":
            data.del_redund_seq_tax()
            output_saver.write("** SEQUENCES VARIANTS have been removed by TAXON.\n")
        else: 
            data.del_redund_seq_glob()
            output_saver.write("** SEQUENCES VARIANTS have been removed GLOBAL.\n")
    output_saver.write(data.get_info("after deleting redundant sequence variant")+"\n")
    if redundant_amplicon == True:
        output_saver.write("\n\
=============================================================================================\n")
        if filtering_type == "taxon":
            data.del_redund_ampli_tax()
            output_saver.write("* 'REDUNDANT AMPLICONS' have been removed by 'TAXON'.\n")
        else: 
            data.del_redund_ampli_glob()
            output_saver.write("* 'REDUNDANT AMPLICONS' have been removed 'GLOBAL'.\n")
        output_saver.write(data.get_info("after deleting redundant amplicons")+"\n")
else:
    output_saver.write("* You have decided 'NOT' to 'FILTER' your dataset.\n")

output_saver.write("\n\n\n\
=============================================================================================\n\
=========                  PART 3/3                Exporting                        =========\n\
=============================================================================================\n")

# EXPORT
seq_variants_output = args.seq_variants_output
amplicons_output = args.amplicons_output
taxonomy_output = args.taxonomy_output
taxon_list_output = args.taxon_list_output
shared_ampl_output = args.shared_ampl_output
complex_dict_output = args.complex_dict_output
modified_tax_output = args.modified_tax_output
access_dict_output = args.access_dict_output
ampl_output_phylo = args.ampl_output_phylo
seq_output_phylo = args.seq_output_phylo

flag = False
if seq_variants_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your qiime format SEQUENCES VARIANTS as %s.\n"%seq_variants_output)
    data.export_seq_fasta(seq_variants_output, "qiime")

if seq_output_phylo != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your ncbi format SEQUENCES VARIANTS as %s.\n"%seq_output_phylo)
    data.export_seq_fasta(seq_output_phylo, "ncbi")

if amplicons_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your qiime format AMPLICONS as %s.\n"%amplicons_output)
    data.export_ampli_fasta(amplicons_output, "qiime")

if ampl_output_phylo != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your ncbi format AMPLICONS as %s.\n"%ampl_output_phylo)
    data.export_ampli_fasta(ampl_output_phylo, "ncbi")

if taxonomy_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your TAXONOMY file as '%s'.\n"%taxonomy_output)
    data.export_tax_id_txt(taxonomy_output, "\t")

if taxon_list_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your TAXON LIST file as '%s'.\n"%taxon_list_output)
    data.export_taxon_list(taxon_list_output)

if shared_ampl_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your SHARED AMPLICONS file as '%s'.\n"%shared_ampl_output)
    data.export_shared_ampli(shared_ampl_output)

if complex_dict_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your COMPLEX DICT file as '%s'.\n"%complex_dict_output)
    data.export_complex_dict(complex_dict_output, "\t")

if modified_tax_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your MODIFIED TAX LIST file as '%s'.\n"%modified_tax_output)
    data.export_modified_tax(modified_tax_output)

if access_dict_output != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your ACCESS DICT JSON file as '%s'.\n"%modified_tax_output)
    data.export_access_dict(access_dict_output)

if infos_file != None:
    flag = True
    output_saver.write("* You have decided to 'EXPORT' your INFOS file as '%s'.\n"%infos_file)
    output_saver.flush()

if flag == False:
    output_saver.write("* You have decided 'NOT' to 'EXPORT' your data.\n")
    print("ALL DONE ! NO OUTPUT FILES GENERATED.")
else:
    print("\nALL DONE ! YOUR FILE HAVE BEEN SUCCESSFULLY GENERATED. \n")
if infos_file == "none":
    output_saver.flush()
    output_saver.delete()
print("Closing the program... \n")




