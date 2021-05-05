import argparse, sys, getpass, datetime, platform
from . import database as db
from . import utils

__doc__ =  """This package is released and maintained by the Research and Development Institute 
for Agri-Environnement (IRDA - Quebec). More informations at https://github.com/cplessis/q2_mkrefdb."""

def get_arguments():
    parser = argparse.ArgumentParser(description=__doc__)  
    subparsers = parser.add_subparsers()
    subparser_create = subparsers.add_parser("create")
    subparser_filter = subparsers.add_parser("filter")
    subparser_export = subparsers.add_parser("export")
    
    #---------------------------------------------------
    #                  CREATE
    #---------------------------------------------------
    subparser_create.add_argument('-i',
                          '--sequences_input',
                          help="File path to the raw sequences FASTA file to treat.",
                          required=True)
    subparser_create.add_argument('-db',
                          '--source_database',
                          choices=['ncbi', 'ebi', 'ddbj'],
                          help="Name of the database from which the FASTA file was obtained.",
                          required=True)                          
    subparser_create.add_argument('-fp',
                          '--forward_primer',
                          help="File path to the forward primer FASTA file.",
                          required=True)
    subparser_create.add_argument('-fmt',
                          '--fw_mismatch_tol',
                          help="Number of mismatch which are accepted for the forward primer annealing.",
                          type=int,
                          default=3)
    subparser_create.add_argument('-rp',
                          '--reverse_primer',
                          help="File path to the reverse primer FASTA file.",
                          required=True)
    subparser_create.add_argument('-rmt',
                          '--rv_mismatch_tol',
                          help="Number of mismatch which are accepted for the reverse primer annealing.",
                          type=int,
                          default=3)
    subparser_create.add_argument('-tp',
                          '--trim_primers',
                          help="""This option keep the primers on amplicons if specified.
                          By default the primers are trimmed.""",
                          default= True)
    subparser_create.add_argument('-mto',
                          '--modified_tax_output',
                          help="File name of the modified taxon list TXT output.",
                          default=None)
    subparser_create.add_argument('-o',
                          '--output_database_json',
                          help="File name of the access dictionnary JSON output.",
                          required=True)
    subparser_create.add_argument('-inf',
                          '--infos_file',
                          help="Informations file in which all the files treatments are recorded.",
                          default=None)
    subparser_create.add_argument('-dit',
                          '--displ_inf_terminal',
                          help="Display information file in terminal if arg is specified.",
                          action = 'store_true',
                          default=False)
    #---------------------------------------------------
    #                  FILTER
    #---------------------------------------------------
    subparser_filter.add_argument('-i',
                          '--database_json',
                          help="Path to the database JSON file to treat.",
                          required=True)
    subparser_filter.add_argument('-cc',
                          '--custom_complex',
                          help="File path to the custom complex CSV file.",
                          default=None)
    subparser_filter.add_argument('-g1',
                          '--genus1',
                          help="Name of the first genus to keep.\
                          Ex: Fusarium",
                          default="None")
    subparser_filter.add_argument('-g2',
                          '--genus2',
                          help="Name of the second genus to keep.\
                          Ex: Gibberella",
                          default="None")
    subparser_filter.add_argument('-u',
                          '--unverified',
                          action = 'store_true',
                          help="Keep the unverified genus if arg is specified.",
                          default=False)
    subparser_filter.add_argument('-tbc',
                          '--taxonomy_by_complex',
                          help="Do not group species with same complex lineage if arg is specified.",
                          action = 'store_false',
                          default=True)
    subparser_filter.add_argument('-rsv',
                          '--redund_seq_variants',
                          action = 'store_false',
                          help="Keep the redundant sequences variants if arg is specified.",
                          default=True)
    subparser_filter.add_argument('-ra',
                          '--redundant_amplicon',
                          action = 'store_false',
                          help="Keep the redundant amplicons if arg is specified.",
                          default=True)
    subparser_filter.add_argument('-nasv',
                          '--seq_without_ampl',
                          action = 'store_false',
                          help="Keep the sequences without amplicons if arg is specified.",
                          default=True)
    subparser_filter.add_argument('-o',
                          '--output_database_json',
                          help="File name of the access dictionnary JSON output.",
                          default="database.json")
    subparser_filter.add_argument('-inf',
                          '--infos_file',
                          help="Informations file in which all the files treatments are recorded.",
                          default=None)
    subparser_filter.add_argument('-dit',
                          '--displ_inf_terminal',
                          help="Display information file in terminal if arg is specified.",
                          action = 'store_true',
                          default=False)
    #---------------------------------------------------
    #                Export
    #---------------------------------------------------
    subparser_export.add_argument('-i',
                          '--database_json',
                          help="Path to the database JSON file to treat.",
                          required=True) 
    subparser_export.add_argument('-svo',
                          '--seq_variants_output',
                          help="File name of the sequences variants FASTA output.",
                          default=None)
    subparser_export.add_argument('-svop',
                          '--seq_output_phylo',
                          help="File name of the sequences variants FASTA output in NCBI format (for phylogeny).",
                          default=None)
    subparser_export.add_argument('-ao',
                          '--amplicons_output',
                          help="File name of the amplicons FASTA output.",
                          default=None)
    subparser_export.add_argument('-aop',
                          '--ampl_output_phylo',
                          help="File name of the amplicon variants FASTA output in NCBI format (for phylogeny).",
                          default=None)
    subparser_export.add_argument('-to',
                          '--taxonomy_output',
                          help="File name of the taxonomy TXT output.",
                          default=None)
    subparser_export.add_argument('-tlo',
                          '--taxon_list_output',
                          help="File name of the taxon list TXT output.",
                          default=None)
    subparser_export.add_argument('-sao',
                          '--shared_ampl_output',
                          help="File name of the shared amplicons TXT output.",
                          default=None)
    subparser_export.add_argument('-cdo',
                          '--complex_dict_output',
                          help="File name of the complex dictionnary TXT output.",
                          default=None)
    subparser_export.add_argument('-inf',
                          '--infos_file',
                          help="Informations file in which all the files treatments are recorded.",
                          default=None)
    subparser_export.add_argument('-dit',
                          '--displ_inf_terminal',
                          help="Display information file in terminal if arg is specified.",
                          action = 'store_true',
                          default=False)
    
    args = parser.parse_args()
    return args

args = get_arguments()
args_dict = args.__dict__

# Action type detection
action_type = ""
if "source_database" in args_dict: action_type = "create"
if "genus1" in args_dict: action_type = "filter"
if "seq_variants_output" in args_dict: action_type = "export"

print("\n\n\n\
=============================================================================================\n\
=========                               Q2_MkRefDb                                  =========\n\
=============================================================================================\n")

# INIT
infos_file = args.infos_file
displ_inf_terminal = args.displ_inf_terminal

# Init input data: OUTPUT FILE INFO
if infos_file: output_saver = utils.Logger(infos_file, displ_inf_terminal)
else: output_saver = utils.Logger("temp_output_saver.txt", displ_inf_terminal)

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

data = db.Database()

def write_db_infos(database):
    output_saver.write("Fasta file name : "+database.file_name+"\n")
    output_saver.write("Original database of the fasta file : "+database.origin_database+"\n")
    output_saver.write("Forward primer : "+str(database.forward_primer_list) +"\n")
    output_saver.write("Reverse primer : "+str(database.reverse_primer_list)+"\n")
    output_saver.write("Forward mismatch tolerance : "+str(database.fw_mismatch_tol)+"\n")
    output_saver.write("Reverse mismatch tolerance : "+str(database.rv_mismatch_tol)+"\n")
    output_saver.write("Trimmed primer : "+str(database.trim_primers)+"\n")

# Initiate if the CREATE subparser is called
if action_type == "create":
    print("\n\n\n\
------------------------------------\n\
                Create              \n\
------------------------------------\n")
    # Init input data: FASTA FILE from NCBI, EBI or DDBJ
    data.create(args.sequences_input, args.source_database, args.forward_primer, args.reverse_primer, \
        args.fw_mismatch_tol, args.rv_mismatch_tol, args.trim_primers)
    write_db_infos(data)
    output_saver.write(data.get_info("init data")+"\n")

# Initiate if the FILTER or EXPORT subparser is called
if action_type in ["filter", "export"]:
    data.import_db(args.database_json)
    write_db_infos(data)
    output_saver.write(data.get_info("init data")+"\n")


output_saver.write("\n\n\n\
=============================================================================================\n\
=========                            Filter Informations                            =========\n\
=============================================================================================\n")
if action_type == "filter":
    print("\n\n\n\
------------------------------------\n\
                FILTER              \n\
------------------------------------\n")
    genus1 = args.genus1
    genus2 = args.genus2
    unverified = args.unverified
    taxonomy_by_complex = args.taxonomy_by_complex
    seq_without_ampl = args.seq_without_ampl
    custom_complex = args.custom_complex
    redund_seq_variants = args.redund_seq_variants
    redundant_amplicon = args.redundant_amplicon

    if taxonomy_by_complex == True:
        data.group_by_complex()
        output_saver.write("* You have decided to 'GROUP TAXON by COMPLEX' in your dataset.\n")
        output_saver.write(data.get_info("After grouping by complex \n"))
        output_saver.write("\n\
=============================================================================================\n")
    if custom_complex: 
        data.custom_complex(custom_complex)
        output_saver.write("* You have decided to use a 'CUSTOM TAXONOMY'.\n")
        output_saver.write(data.get_info("After grouping by complex with custom taxonomy \n"))
        output_saver.write("\n\
=============================================================================================\n")
    if genus1 != "None":
        output_saver.write("* You have decided to 'CLEAN UP' your dataset.\n")
        genus_set = data.get_all_genus()
        if genus1 in genus_set:
            output_saver.write(genus1+" is in your dataset."+"\n")
        else:
            output_saver.write("WARNING : There is NO genus named "+genus1+" in your dataset."+"\n")
        if (genus2 != "None") & (genus2 in genus_set):
            output_saver.write(genus2+" is in your dataset."+"\n")
        elif genus2 != "None":
            output_saver.write("WARNING : There is NO genus named "+genus2+" in your dataset."+"\n")
        data.clean_dataset(genus1, genus2, unverified)
        if unverified: 
            output_saver.write("* You have decided to 'KEEP' the 'UNVERIFIED' and 'UNCULTURED' genus.\n")
        else: 
            output_saver.write("* 'UNVERIFIED' and 'UNCULTURED' genus have been removed.\n")
        output_saver.write(data.get_info("After Cleaning with : "+genus1+" -- "+genus2+" -- UNVERIFIED="+str(unverified))+"\n")

    if seq_without_ampl:
        output_saver.write("\n\
=============================================================================================\n")
        data.del_na_amplicons()
        output_saver.write("** 'SEQUENCES VARIANTS' have been removed if 'NO AMPLICONS'.\n")
        output_saver.write(data.get_info("after deleting sequences variants without amplicons")+"\n")
    else:
        output_saver.write("** You decided to 'KEEP' the 'SEQUENCES VARIANTS' with 'NO AMPLICONS'.\n")

    if redund_seq_variants:
        output_saver.write("\n\
=============================================================================================\n")
        data.del_redund_seq_tax()
        output_saver.write("** SEQUENCES VARIANTS have been removed by TAXON.\n")
        
    output_saver.write(data.get_info("after deleting redundant sequence variant")+"\n")
    if redundant_amplicon:
        output_saver.write("\n\
=============================================================================================\n")
        data.del_redund_ampli_tax()
        output_saver.write("* 'REDUNDANT AMPLICONS' have been removed by 'TAXON'.\n")
        output_saver.write(data.get_info("after deleting redundant amplicons")+"\n")





output_saver.write("\n\n\n\
=============================================================================================\n\
=========                            Export Informations                            =========\n\
=============================================================================================\n")

if action_type == "export":
    print("\n\n\n\
------------------------------------\n\
                EXPORT              \n\
------------------------------------\n")
    # EXPORT
    seq_variants_output = args.seq_variants_output
    amplicons_output = args.amplicons_output
    taxonomy_output = args.taxonomy_output
    taxon_list_output = args.taxon_list_output
    shared_ampl_output = args.shared_ampl_output
    complex_dict_output = args.complex_dict_output
    ampl_output_phylo = args.ampl_output_phylo
    seq_output_phylo = args.seq_output_phylo

    if seq_variants_output:
        output_saver.write("* You have decided to 'EXPORT' your qiime format SEQUENCES VARIANTS \
            as %s.\n"%seq_variants_output)
        data.export_seq_fasta(seq_variants_output, "qiime")

    if seq_output_phylo:
        output_saver.write("* You have decided to 'EXPORT' your ncbi format SEQUENCES VARIANTS \
            as %s.\n"%seq_output_phylo)
        data.export_seq_fasta(seq_output_phylo, "ncbi")

    if amplicons_output:
        output_saver.write("* You have decided to 'EXPORT' your qiime format AMPLICONS \
            as %s.\n"%amplicons_output)
        data.export_ampli_fasta(amplicons_output, "qiime")

    if ampl_output_phylo:
        output_saver.write("* You have decided to 'EXPORT' your ncbi format AMPLICONS \
            as %s.\n"%ampl_output_phylo)
        data.export_ampli_fasta(ampl_output_phylo, "ncbi")

    if taxonomy_output:
        output_saver.write("* You have decided to 'EXPORT' your TAXONOMY file \
            as '%s'.\n"%taxonomy_output)
        data.export_tax_id_txt(taxonomy_output, "\t")

    if taxon_list_output:
        output_saver.write("* You have decided to 'EXPORT' your TAXON LIST file \
            as '%s'.\n"%taxon_list_output)
        data.export_taxon_list(taxon_list_output)

    if shared_ampl_output:
        output_saver.write("* You have decided to 'EXPORT' your SHARED AMPLICONS file \
            as '%s'.\n"%shared_ampl_output)
        data.export_shared_ampli(shared_ampl_output)

    if complex_dict_output:
        output_saver.write("* You have decided to 'EXPORT' your COMPLEX DICT file \
            as '%s'.\n"%complex_dict_output)
        data.export_complex_dict(complex_dict_output, "\t")

if action_type == "create":
    modified_tax_output = args.modified_tax_output
    if modified_tax_output:
        output_saver.write("* You have decided to 'EXPORT' your MODIFIED TAX LIST file \
            as '%s'.\n"%modified_tax_output)
        data.export_modified_tax(modified_tax_output)

if action_type in ["create", "filter"]:
    output_database_json = args.output_database_json
    output_saver.write("* DATABASE JSON file exported \
        as '%s'.\n"%output_database_json)
    data.export_access_dict(output_database_json)

if infos_file:
    output_saver.write("* INFOS TXT file exported as '%s'.\n"%infos_file)
    output_saver.flush()
else:
    output_saver.flush()
    output_saver.delete()

print("\nALL DONE ! %s FILES HAVE BEEN SUCCESSFULLY GENERATED. \n"%action_type.upper())