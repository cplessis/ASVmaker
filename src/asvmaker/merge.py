import os
from . import database as db

def __check_primers(db1, db2):
    if db1.trim_primers != db2.trim_primers:
        print("!! WARNING !!   Please, before merging make sure primers are kept or trimmed for both database.")
        exit()
    return

def open_sa_file(sa_file):
    sa_dict = {}
    with open(sa_file) as file:
        for line in file.readlines():
            line = line.split()
            sa_dict[line[1]] = {"SA":line[0], "access":line[2:]}
    return sa_dict 

def merge(db1, db2, sa_file1, sa_file2, output_Database):
    __check_primers(db1, db2)

    sa1_dict = open_sa_file(sa_file1)
    sa2_dict = open_sa_file(sa_file2)

    max_sa = len(sa1_dict)
    db1_ampl_set = {access["amplicon"]:access["taxon"] for access in db1.access_dict.values()}

    for access in db2.access_dict:
        if db2.get_amplicon(access) not in db1_ampl_set:
            db1.seq_dict[db2.get_name(access)] = db2.get_sequence(access)
            db1.access_dict[access] = db2.access_dict[access]
            if access in sa2_dict:
                max_sa+=1
                db1.access_dict[access]["taxon"] = db2.get_genus(access)+"_"+"SA"+str(max_sa)
                sa1_dict[access] = {"SA":"SA"+str(max_sa), "access":sa2_dict[access]["access"]}

    db1.update_data2()
    
    with open(output_Database.rstrip(".json")+"SA_ext.txt", "w") as sa_file3:
        for access in sa1_dict:
            sa_file3.write(sa1_dict[access]["SA"]+"\t"+access+"\t"+"\t".join(sa1_dict[access]["access"])+"\n")

def merge_s():
    pass

def merge_g(db1, db2, sa_dict):
    db1_ampl_dict = {db1.get_amplicon(access):{"taxon":db1.get_taxon(access), 
        "access":access} for access in db1.access_dict}
    max_sa = len(sa_dict)
    for access in db2.access_dict:
        ampl2 = db2.get_amplicon(access)
        if ampl2 in db1_ampl_dict:
            if ampl2 in sa_dict:
                sa_dict[ampl2]["access"].append(access+"|"+db2.get_taxon(access))
            else:
                max_sa += 1
                db1acccess = db1_ampl_dict[ampl2]["access"]
                sa_dict[ampl2] = {"access":"", "SA":""}
                sa_dict[ampl2]["SA"] = max_sa
                sa_dict[ampl2]["access"] = [db1acccess+"|"+db1.get_taxon(db1acccess), access+"|"+db2.get_taxon(access)]
                db1.access_dict[db1acccess]["taxon"] = "Genus_SA"+str(max_sa)
                lin = db1.access_dict[db1acccess]["lineage"].split("; ")
                lin[-2] = "Genus_SA"
                db1.access_dict[db1acccess]["lineage"] = "; ".join(lin)
        else:
            db1.seq_dict[db2.get_name(access)] = db2.get_sequence(access)
            db1.access_dict[access] = db2.access_dict[access]

def merge_all(dir_path, s_or_g, outputDB):
    all_dbs = db.Database()
    all_dbs.import_db(dir_path+os.listdir(dir_path)[0])
    sa_dict = {}
    for file in os.listdir(dir_path)[1:]:
        curr_db = db.Database()
        curr_db.import_db(dir_path+file)
        if s_or_g == "g":
            merge_g(all_dbs, curr_db, sa_dict)
    
    all_dbs.update_data2()    
    all_dbs.export_ampli_fasta(outputDB.rstrip(".json")+"-ampli.fasta", "qiime")
    all_dbs.export_access_dict(outputDB)
    with open(outputDB.rstrip(".json")+"-SA_ext.txt", "w") as sa_file3:
        for sa in sa_dict:
            sa_file3.write("SA"+str(sa_dict[sa]["SA"])+"\t"+"\t".join(sa_dict[sa]["access"])+"\n")
