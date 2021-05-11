from . import utils

def phylo_annotation(Database, output_file):
    with open(output_file, "w") as file:
        file.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for taxon in Database.taxon_dict:
            if taxon.split("_")[1][0] != "S": 
                color = utils.random_color()
                for access_nb in Database.taxon_dict[taxon]:
                    file.write(access_nb+"_"+taxon+"\tlabel\trgb"+str(color)+"\tnormal\t1\n")
                    file.write(access_nb+"_"+taxon+"\tbranch\trgb"+str(color)+"\tnormal\t1\n")
