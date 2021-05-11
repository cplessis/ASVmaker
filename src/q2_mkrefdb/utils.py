import json, sys, os, urllib.request, random
from pathlib import Path

def save_as_json(variable, file_json):
    with open(file_json, 'w') as saved_file:
    	json.dump(variable, saved_file)

def open_from_json(file_json):
	with open(file_json) as saved_file:
		saved_variable = json.load(saved_file)
	return saved_variable

def export_dict_csv(wanted_dict, output_file_name, separator):
    """Dictionnary as :   'key : x' """
    with open(output_file_name, "w") as output_file:
        for key in wanted_dict:
            output_file.write(key+separator+str(wanted_dict[key])+"\n")

def export_dict_csv2(wanted_dict, output_file_name, separator):
    """Dictionnary as :   'key : {x, y, z}' """
    with open(output_file_name, "w") as output_file:
        for key in wanted_dict:
            line = key
            for subkey in wanted_dict[key]:
                line += separator+subkey
            output_file.write(line+"\n")

def export_list_csv(wanted_list, output_file_name):
    with open(output_file_name, "w") as output_file:
        for key in wanted_list:
            output_file.write(key+"\n")

# To try if a filepath exists or not. 
def try_path(file_to_test):
    flag = False
    while flag == False:
        try:
            Path(file_to_test).resolve(strict=True)
            flag = True
        except FileNotFoundError:
            file_to_test = input("ERROR : Please enter a VALID file path for '%s' !\n>" %file_to_test)
    return file_to_test         

# To print both in terminal and in an output file. 
class Logger(object):
    def __init__(self, log_file, display):
        self.log_file = log_file
        self.terminal = sys.stdout
        self.log = open(log_file, "w")
        self.display = display

    def write(self, message):
        if self.display == True: 
            self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        self.log.close()
    
    def delete(self):
        os.remove(self.log_file)

# Variable and Binary object
def binary_to_variable(the_binary, var_type):
    if var_type == "json": return json.loads(the_binary.decode('utf-8'))

# ULR to Variable
def get_var_from_url(url, var_type):
    with urllib.request.urlopen(url) as response:
        var = binary_to_variable(response.read(), var_type)
        return var

# COLORS
def random_color():
    """Return a random tuple RGB (r, g, b) as (10, 235, 53) for example.
    """
    return tuple([random.choice(range(256)), random.choice(range(256)), random.choice(range(256))])