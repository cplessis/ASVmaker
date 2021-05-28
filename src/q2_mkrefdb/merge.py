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

def merge(Database1, Database2, sa_file1, sa_file2, output_Database):
    db1 = Database1
    db2 = Database2
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