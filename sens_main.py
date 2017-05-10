import os, sys, csv
from random import randint
import itertools
import thermo_get_properties as thermo_reader
import sampling as sam

###################################
# @author Shimaa_Gamil
# KAUST-CCRC
###################################

# STEPS:
# 1. Read Thermo file
# 2. Caclulate thermo preperties of thermo file (Hf, S, CPs)
# 3. Define the Sampling Method (# of samples, normal distribution)
# 4. Do Sampling
# 5. Generate many files with differnt properties based on sampling values
# 6. Generate Thermo files of the above different thermo properties files


###
#  a function to call the module that calculates thermo properties from therm file
###
def readThermoFile (input_file): 
        path = os.path.dirname(input_file)        
        file_name_only = os.path.basename(input_file).split('.')[0]
        output_file = path + '\\' + file_name_only + '_props.csv'
        thermo_reader.readThermoFile (input_file, output_file)
        return output_file


###
# a function to calculate thermo data for the passed properties file,
# then call a function to rewrite the source thermo file and update thermo data of the butanol species
###
def writeNewThermoFile (butanol_species_list, props_file, source_thermo_file, output_thermo_file):

        butanol_species_props_list = []
        
        with open(props_file, newline='') as props_f:
                
            reader = csv.reader(props_f)
            
            for row in reader:
                    
                    species = row [0]
                    if species in butanol_species_list:
                            butanol_species_props_list.append ( row )                            
        
        
        thermo_reader.reWriteThermoFile (butanol_species_list, butanol_species_props_list, source_thermo_file, output_thermo_file)

###
# a function to calculate thermo data for the passed properties file,
# then call a function to rewrite the source thermo file and update thermo data of the butanol species
###
def getbutanolSpeciesThermo (butanol_species_list, props_file, output_file):

        butanol_species_props_list = []

        butanol_props_file = open(output_file, "wt", newline='')
        writer = csv.writer(butanol_props_file)
        
        with open(props_file, newline='') as props_f:
                
            reader = csv.reader(props_f)
            
            for row in reader:
                    
                    species = row [0]
                    if species in butanol_species_list:
                            writer.writerow (row)                            
        

        butanol_props_file.close () 
                

### 
# main code that runs first
###

# You should replace these file paths with your therm file, and the dic file that
# include your target species for doing sensitivity
thermo_file_path = 'C:\\Shimaa_Work\\Thermo_UQ\\PECS_therm.therm'
species_type_file_path = 'C:\\Shimaa_Work\\Thermo_UQ\\butanol_species_dict.csv'

# Read thermo file and calculate properties
properties_file = readThermoFile (thermo_file_path)
if not os.path.exists(properties_file): # handle the case if the props file is not created for any reason
        print ("An error has occured while calculating Thermo properties")
        sys.exit(0)


# Read the important species to do UQ on their parameters
species_type_dict = {}
butanol_species_list = []

with open(species_type_file_path, newline='') as type_f:

        reader = csv.reader(type_f)

        for row in reader:
                
                butanol_species_list.append (row[0])
                species_type = row[1]
                
                if 'alcohol like' in species_type:
                        species_type_dict[ row[0] ] = sam.ALCOHOL_LIKE_TYPE
                        
                elif 'hydroxy alkyl radical' in  species_type:    
                        species_type_dict[ row[0] ] = sam.HYDROXY_ALKYL_RADICAL_TYPE
                        
                elif 'peroxy radical' in  species_type:    
                        species_type_dict[ row[0] ] = sam.PEROXY_RADICAL_TYPE
                        
                elif 'hydroperoxide' in  species_type:     
                        species_type_dict[ row[0] ] = sam.HYDROPEROXIDE_TYPE
                        
                elif 'alkoxy radical like propoxy' in  species_type:     
                        species_type_dict[ row[0] ] = sam.ALKOXY_RADICAL_LIKE_PROPOXY_TYPE

butanol_species_flag_list = [0] * len (butanol_species_list)

#print (species_type_dict)
#print (butanol_species_list)

output_file = os.path.dirname(properties_file) + '\\butanol_species_thermo_props.csv'
getbutanolSpeciesThermo (butanol_species_list, properties_file, output_file)


# create a directory with a unique key to write output files
unique_key = randint(100, 100000)
path = os.path.dirname(properties_file)
directory = path + '\\' + 'sampling_output_' + str(unique_key) 
if not os.path.exists(directory):
    os.makedirs(directory)

Hf_directory = directory + '\\Hf'
S_directory = directory + '\\S'

if not os.path.exists(Hf_directory):
    os.makedirs(Hf_directory)

if not os.path.exists(S_directory):
    os.makedirs(S_directory)

for species_name in butanol_species_list:
        
        with open(properties_file, newline='') as props_f:
                
                reader = csv.reader(props_f)
                
                file_name_only = os.path.basename(properties_file).split('.')[0]
                Hf_output_props_maxSD_file = Hf_directory + '\\' + file_name_only + '_' + species_name + '_Hf_maxSD' +  '.csv'
                Hf_output_props_minSD_file = Hf_directory + '\\' + file_name_only + '_' + species_name + '_Hf_minSD' +  '.csv'
                Hf_output_thermo_maxSD_file = Hf_directory + '\\thermo_' + species_name + '_Hf_maxSD' +  '.dat'
                Hf_output_thermo_minSD_file = Hf_directory + '\\thermo_' + species_name + '_Hf_minSD' +  '.dat'
                
                S_output_props_maxSD_file = S_directory + '\\' + file_name_only + '_' + species_name + '_S_maxSD' +  '.csv'
                S_output_props_minSD_file = S_directory + '\\' + file_name_only + '_' + species_name + '_S_minSD' +  '.csv'
                S_output_thermo_maxSD_file = S_directory + '\\thermo_' + species_name + '_S_maxSD' +  '.dat'
                S_output_thermo_minSD_file = S_directory + '\\thermo_' + species_name + '_S_minSD' +  '.dat'

                Hf_output_max = open(Hf_output_props_maxSD_file, "wt", newline='')
                Hf_writer_max = csv.writer(Hf_output_max)

                Hf_output_min = open(Hf_output_props_minSD_file, "wt", newline='')
                Hf_writer_min = csv.writer(Hf_output_min)

                S_output_max = open(S_output_props_maxSD_file, "wt", newline='')
                S_writer_max = csv.writer(S_output_max)

                S_output_min = open(S_output_props_minSD_file, "wt", newline='')
                S_writer_min = csv.writer(S_output_min)
                
                for row in reader:

                        species_name_file = row [0]
                        
                        if species_name == species_name_file:
                                Hf = row [1]
                                S = row [2]
                                newHf_max = Hf
                                newS_max = S
                                newHf_min = Hf
                                newS_min = S
                                
                                species_type = species_type_dict [species_name]
                                
                                if species_type == sam.ALCOHOL_LIKE_TYPE:
                                        newHf_max = float(Hf) + sam.ALCOHOL_LIKE_HF_SD
                                        newS_max = float(S) + sam.ALCOHOL_LIKE_S_SD
                                        
                                        newHf_min = float(Hf) - sam.ALCOHOL_LIKE_HF_SD
                                        newS_min = float(S) - sam.ALCOHOL_LIKE_S_SD
                                        
                                elif species_type == sam.HYDROXY_ALKYL_RADICAL_TYPE:
                                        newHf_max = float(Hf) + sam.HYDROXY_ALKYL_RADICAL_HF_SD
                                        newS_max = float(S) + sam.HYDROXY_ALKYL_RADICAL_S_SD

                                        newHf_min = float(Hf) - sam.HYDROXY_ALKYL_RADICAL_HF_SD
                                        newS_min = float(S) - sam.HYDROXY_ALKYL_RADICAL_S_SD
                                        
                                elif species_type == sam.PEROXY_RADICAL_TYPE:
                                        newHf_max = float(Hf) + sam.PEROXY_RADICAL_HF_SD
                                        newS_max = float(S) + sam.PEROXY_RADICAL_S_SD

                                        newHf_min = float(Hf) - sam.PEROXY_RADICAL_HF_SD
                                        newS_min = float(S) - sam.PEROXY_RADICAL_S_SD
                                        
                                elif species_type == sam.HYDROPEROXIDE_TYPE:
                                        newHf_max = float(Hf) + sam.HYDROPEROXIDE_HF_SD
                                        newS_max = float(S) + sam.HYDROPEROXIDE_S_SD

                                        newHf_min = float(Hf) - sam.HYDROPEROXIDE_HF_SD
                                        newS_min = float(S) - sam.HYDROPEROXIDE_S_SD                                      
                                        
                                elif species_type == sam.ALKOXY_RADICAL_LIKE_PROPOXY_TYPE:
                                        newHf_max = float(Hf) + sam.ALKOXY_RADICAL_LIKE_PROPOXY_HF_SD
                                        newS_max = float(S) + sam.ALKOXY_RADICAL_LIKE_PROPOXY_S_SD

                                        newHf_min = float(Hf) - sam.ALKOXY_RADICAL_LIKE_PROPOXY_HF_SD
                                        newS_min = float(S) - sam.ALKOXY_RADICAL_LIKE_PROPOXY_S_SD


                                row_max_Hf = row [:]
                                row_min_Hf = row [:]

                                row_max_S = row [:]
                                row_min_S = row [:]

                                row_max_Hf [1] = newHf_max
                                row_min_Hf [1] = newHf_min
                                
                                row_max_S [2] = newS_max
                                row_min_S [2] = newS_min

                                Hf_writer_max.writerow (row_max_Hf)
                                Hf_writer_min.writerow (row_min_Hf)

                                S_writer_max.writerow (row_max_S)
                                S_writer_min.writerow (row_min_S)

                                print (row)
                                print (row_max_Hf)
                                print (row_min_Hf)
                                print (row_max_S)
                                print (row_min_S)
                                
                                #print ( 'Species Name', species_name , 'Species Type = ' , species_type, 'Hf = ', Hf, 'new Hf = ', newHf)
                                #print ('S = ', S, 'new S = ', newS)

                        else:

                                Hf_writer_max.writerow (row)
                                Hf_writer_min.writerow (row)

                                S_writer_max.writerow (row)
                                S_writer_min.writerow (row)

                Hf_output_max.close ()                        
                Hf_output_min.close ()
                S_output_max.close ()
                S_output_min.close ()

                
                ## Recalculate and write therrmo file based on the new generated properties file
                writeNewThermoFile (butanol_species_list, Hf_output_props_maxSD_file, thermo_file_path, Hf_output_thermo_maxSD_file)
                writeNewThermoFile (butanol_species_list, Hf_output_props_minSD_file, thermo_file_path, Hf_output_thermo_minSD_file)
                writeNewThermoFile (butanol_species_list, S_output_props_maxSD_file, thermo_file_path, S_output_thermo_maxSD_file)
                writeNewThermoFile (butanol_species_list, S_output_props_minSD_file, thermo_file_path, S_output_thermo_minSD_file)
              
                              
#print (len(butanol_species_list))


