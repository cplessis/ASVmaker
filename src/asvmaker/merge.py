import os
from . import database as db

def __check_primers(db1, db2):
    if db1.trim_primers != db2.trim_primers:
        print("!! WARNING !!   Please, before merging make sure primers are kept or trimmed for both database.")
        exit()
    return

def is_SA(taxon):
    if taxon.split("_")[1][:2] == "SA":
        return True
    else : return False

def any_SA(taxon_list):
    for i in taxon_list:
        if is_SA(i):
            return True
    else : return False

def open_sa_file(sa_file):
    sa_dict = {}
    with open(sa_file, "r") as file:
        for line in file.readlines():
            line = line.split()
            sa = line[0] 
            for i in line[1:]:
                access = i.split("|")[0]
                taxon = i.split("|")[1]
                sa_dict[access] = dict(taxon=taxon, sa=sa)
    return sa_dict 

def merge(db1, db2, sa_file1, sa_file2, output_database):
    __check_primers(db1, db2)

    sa1_dict = open_sa_file(sa_file1)
    sa2_dict = open_sa_file(sa_file2)

    max_sa = len(set([sa1_dict[i]["sa"] for  i in sa1_dict]))

    z = lambda x : [sa1_dict[i]["taxon"] for i in sa1_dict if sa1_dict[i]["sa"] == f"SA{str(x)}"]
    sa_out_dict = {f"SA{str(i)}":z(i) for i in range(1, max_sa+1)}

    db1_ampl_dict = {
        db1.access_dict[access]["amplicon"]:{
            "taxon":db1.access_dict[access]["taxon"],
            "access":access
        } for access in db1.access_dict}

    for access in db2.access_dict:
        amplicon = db2.get_amplicon(access)
        
        #if amplicon not in DB1
        if amplicon not in db1_ampl_dict:
            db1.seq_dict[db2.get_name(access)] = db2.get_sequence(access)
            db1.access_dict[access] = db2.access_dict[access]
            if access in sa2_dict:
                # Si amplicon de db2 est un SA
                max_sa += 1
                sa1 = "SA"+str(max_sa)
                sa2 = sa2_dict[access]["sa"]
                db1.access_dict[access]["taxon"] = db2.get_genus(access)+"_"+sa1
                sa_out_dict[sa1] = [sa2_dict[i]["taxon"] for i in sa2_dict if sa2_dict[i]["sa"] == sa2]

        else:
            tax1 = db1_ampl_dict[amplicon]["taxon"]
            tax2 = db2.get_taxon(access)
            sa1 = tax1.split('_')[-1]
            sa2 = tax2.split('_')[-1]
                        
            # print(tax1, "----" ,tax2)
            if tax1 == tax2:
                if is_SA(tax1):
                    # ecrire sa2 avec les valeurs de sa1 dans SA out
                    sa1_access = db1_ampl_dict[amplicon]["access"]
                    sa1_taxons = [sa1_dict[i]["taxon"] for i in sa1_dict if sa1_dict[i]["sa"] == sa1]
                    sa2_taxons = [sa2_dict[i]["taxon"] for i in sa2_dict if sa2_dict[i]["sa"] == sa1]
                    sa_out_dict[sa1] = set(sa1_taxons+sa2_taxons)
                else:
                    # Donc si les taxon sont ont le meme nom : Pas d'action
                    pass
            elif any_SA([tax1, tax2]):
                if is_SA(tax1) and is_SA(tax2):
                    # Prendre valeur des 2 SA et ajouter au SA out
                    sa1_taxons = [sa1_dict[i]["taxon"] for i in sa1_dict if sa1_dict[i]["sa"] == sa1]
                    sa2_taxons = [sa2_dict[i]["taxon"] for i in sa2_dict if sa2_dict[i]["sa"] == sa2]
                    sa_out_dict[sa1] = set(sa1_taxons+sa2_taxons)
                elif is_SA(tax1):
                    # Prendre valeur non SA et ajoute au SA db1
                    sa1_taxons = [sa1_dict[i]["taxon"] for i in sa1_dict if sa1_dict[i]["sa"] == sa1]
                    sa_out_dict[sa1] = set(sa1_taxons+[tax2])
                elif is_SA(tax2):
                    # Prendre valeur non SA et valeurs de SA et crer nouvel SA
                    max_sa += 1
                    sa1 = "SA"+str(max_sa)
                    sa2_taxons = [sa2_dict[i]["taxon"] for i in sa2_dict if sa2_dict[i]["sa"] == sa2]
                    sa_out_dict[sa1] = set([tax1]+sa2_taxons)
                    # MAJ de la db
                    access_1 = db1_ampl_dict[amplicon]["access"]
                    db1.access_dict[access_1]["taxon"] = db1.get_genus(access_1)+"_"+"SA"+str(max_sa)
            else :
                # Donc si tax1 != tax2
                # Creer un nouveau SA 
                max_sa += 1
                sa1 = "SA"+str(max_sa)
                sa_out_dict[sa1] = set([tax1,tax2])
                # MAJ de la db
                access_1 = db1_ampl_dict[amplicon]["access"]
                db1.access_dict[access_1]["taxon"] = db1.get_genus(access_1)+"_"+sa1

    with open(output_database.replace(".json","SA_ext.txt"), "w") as sa_file3:
        for sa in sorted(sa_out_dict, key=lambda x: int(x[2:])):
            sa_file3.write(sa+'\t'+'\t'.join(sa_out_dict[sa])+"\n")


def merge_s():
    pass

def merge_g(db1, db2, sa_dict):
    pass

# def merge_g(db1, db2, sa_dict):
#     db1_ampl_dict = {db1.get_amplicon(access):{"taxon":db1.get_taxon(access), 
#         "access":access} for access in db1.access_dict}
#     max_sa = len(sa_dict)
#     for access in db2.access_dict:
#         ampl2 = db2.get_amplicon(access)
#         if ampl2 in db1_ampl_dict:
#             if ampl2 in sa_dict:
#                 sa_dict[ampl2]["access"].append(access+"|"+db2.get_taxon(access))
#             else:
#                 max_sa += 1
#                 db1acccess = db1_ampl_dict[ampl2]["access"]
#                 sa_dict[ampl2] = {"access":"", "SA":""}
#                 sa_dict[ampl2]["SA"] = max_sa
#                 sa_dict[ampl2]["access"] = [db1acccess+"|"+db1.get_taxon(db1acccess), access+"|"+db2.get_taxon(access)]
#                 db1.access_dict[db1acccess]["taxon"] = "Genus_SA"+str(max_sa)
#                 lin = db1.access_dict[db1acccess]["lineage"].split("; ")
#                 lin[-2] = "Genus_SA"
#                 db1.access_dict[db1acccess]["lineage"] = "; ".join(lin)
#         else:
#             db1.seq_dict[db2.get_name(access)] = db2.get_sequence(access)
#             db1.access_dict[access] = db2.access_dict[access]

def merge_all(dir_path, s_or_g, outputDB):
    pass

# def merge_all(dir_path, s_or_g, outputDB):
#     all_dbs = db.Database()
#     all_dbs.import_db(dir_path+os.listdir(dir_path)[0])
#     sa_dict = {}
#     for file in os.listdir(dir_path)[1:]:
#         curr_db = db.Database()
#         curr_db.import_db(dir_path+file)
#         if s_or_g == "g":
#             merge_g(all_dbs, curr_db, sa_dict)
    
#     all_dbs.update_data2()    
#     all_dbs.export_ampli_fasta(outputDB.rstrip(".json")+"-ampli.fasta", "qiime")
#     all_dbs.export_access_dict(outputDB)
#     with open(outputDB.rstrip(".json")+"-SA_ext.txt", "w") as sa_file3:
#         for sa in sa_dict:
#             sa_file3.write("SA"+str(sa_dict[sa]["SA"])+"\t"+"\t".join(sa_dict[sa]["access"])+"\n")
