import os, shutil
import sys
   
def clean_directory(path_str):
	for filename in os.listdir(path_str):
	    	if os.path.isfile(os.path.join(path_str, filename)):
    			os.remove(os.path.join(path_str, filename))
        
if (len(sys.argv) != 4):
	print("python exec <P1> <P2>")
	print("P1: number of layers (if negative, do not encode layers again)")
	print("P2: layer for the results")
	print("P3: model_type (0, 1, 2)")
	exit()

os.system("preproc images 1.3 filtered_images")

def copy_makerset(models):
	for model in models:
		shutil.copyfile('all_markers/' + model, 'markers/' + model)	

from os import listdir
from os.path import isfile, join
import re

#Computes the final accuracy for detection and delineation in each layer [1-4] for every single model (architecture with a single model)
#Previously a pair of markers has been done in each image. Those markers are saved in the folde 'all_markers'
# The result for delineation is in: final_delineation_results.txt
# The result for detection is in: final_detection_results.txt	
def every_single_model_arch():
	
	mypath = "all_markers"
	onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
	onlyfiles.sort()
	
	detection_results = []
	delineation_results = []
	
	for file in onlyfiles:
		models = [file]
		#print("Models: " + ", ".join(models))
		#print(models)
		clean_directory("markers")
		copy_makerset(models)
		results = run_all()
		model_number = getModelNumber(file)
		detection_results.append("Accuracy for model " + model_number + " :"  + results[0])
		delineation_results.append("Accuracy for model " + model_number + " :"  + results[1])
	
	results = run_all()
	detection_results.append("Accuracy for current detection model:"  + results[0])
	delineation_results.append("Accuracy for current delineation model:"  + results[1])
		
	file = open("final_detection_results.txt", 'r+')
	file.truncate(0)
	file.writelines(detection_results)
	file.close()
	
	file = open("final_delineation_results.txt", 'r+')
	file.truncate(0)
	file.writelines(delineation_results)
	file.close()

#Gives the final accuracy for detection and delineation in each layer [1-4] for a fixed model
# The result for delineation is in: final_delineation_results.txt
# The result for detection is in: final_detection_results.txt	
def fixed_model_arch():
	
	detection_results = []
	delineation_results = []

	results = run_all()
	detection_results.append("Accuracy for current detection model:"  + results[0])
	delineation_results.append("Accuracy for current delineation model:"  + results[1])
	
	print(detection_results)
	print(delineation_results)
	
	file = open("final_detection_results.txt", 'r+')
	file.truncate(0)
	file.writelines(detection_results)
	file.close()
	
	file = open("final_delineation_results.txt", 'r+')
	file.truncate(0)
	file.writelines(delineation_results)
	file.close()
		
def getModelNumber(model_name):

	model_number = "".join(re.split("\D", model_name)).lstrip("0")
	return model_number
	

    		    			
def run_all():
	
	clean_directory("bag")
	clean_directory("layer0")
	clean_directory("layer1")
	clean_directory("layer2")
	clean_directory("layer3")
	clean_directory("salie")
	clean_directory("flim_models")
	clean_directory("boxes")
	clean_directory("objs")
	 
	nlayers      = int(sys.argv[1])
	target_layer = int(sys.argv[2])
	model_type   = int(sys.argv[3])
	
	npts_per_marker = 3
	line = "bag_of_feature_points filtered_images markers {} bag".format(npts_per_marker)
	os.system(line)

	delineation_layers = []
	detection_layers = []
	
	for layer in range(1,nlayers+1):
		line = "create_layer_model bag arch.json {} flim_models".format(layer)
		os.system(line)
	
		print("merge_layer_models process\n")
		line = "merge_layer_models arch.json {} flim_models".format(layer)
		os.system(line)
		
		print("encode_merged_layer process\n")
		line = "encode_merged_layer arch.json {} flim_models".format(layer)
		os.system(line)
		
		print("decode_layer process\n")
		line = "decode_layer {} arch.json flim_models {} salie".format(layer, model_type)
		os.system(line)
		
		print("detection process\n")
		line = "detection salie {} boxes".format(layer)
		os.system(line)
	
		print("delineation process\n")
		line = "delineation salie {} objs".format(layer)
		os.system(line)
	
		file = open("delineation_result.txt", "r")
		for f in file:
			delineation_layers.append(f)
		
		file = open("detection_result.txt", "r")
		for f in file:
			detection_layers.append(f)
		
	return [",".join(detection_layers), ",".join(delineation_layers)]
	#print("encode_layer process\n")
	#line = "encode_layer arch.json {} flim_models".format(layer)
	#os.system(line)
	
fixed_model_arch()
