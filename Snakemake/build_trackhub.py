# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Peak Annotate
# snakemake --snakefile build_trackhub.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile build_trackhub.py --configfile Encode.json --rulegraph | dot -Tsvg > build_trackhub.svg
# ################################### IMPORT ##################################### #


import os
import sys
import glob
sys.path.append(os.path.abspath("./library/"))
import utility
# ################################### FUNCTIONS ################################## #


def build_design_Dict(metadata_Dict):
	#
	design_Dict = {}
	for sample, sample_Dict in metadata_Dict.items():
		#
		if sample_Dict["Design"] not in design_Dict:
			#
			design_Dict[sample_Dict["Design"]] = {}
			design_Dict[sample_Dict["Design"]]["Case"] = []
			design_Dict[sample_Dict["Design"]]["Control"] = []
			if sample_Dict["Type"] == "CASE":
				#
				design_Dict[sample_Dict["Design"]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				design_Dict[sample_Dict["Design"]]["Control"].append(sample)
			else:
				pass
		elif sample_Dict["Design"] in design_Dict:
			#
			if sample_Dict["Type"] == "CASE":
				#
				design_Dict[sample_Dict["Design"]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				design_Dict[sample_Dict["Design"]]["Control"].append(sample)
			else:
				pass
	return design_Dict
# ################################# CONFIGURATION ############################# #

# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
TITLE = config_general_Dict["TITLE"]
INFOLINK = config_general_Dict["INFOLINK"]
HPC_DATASHARE = config_general_Dict["HPC_DATASHARE"]
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
LAYOUT = config_data_Dict["LAYOUT"].lower()
GENOME = config_data_Dict["GENOME"].lower()
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#CONDA
config_conda_Dict = config["CONDA"]
CONDA_INIT = config_conda_Dict["CONDA_INIT"]
CONDA_PY2 = config_conda_Dict["ATAC_Seq_py2"]
CONDA_PY3 = config_conda_Dict["ATAC_Seq_py3"]
ACTIVATE_CONDA_PY2 = config_conda_Dict["ACTIVATE_PY2"]
ACTIVATE_CONDA_PY3 = config_conda_Dict["ACTIVATE_PY3"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#CLUSTER
config_cluster_Dict = config["CLSUTER_CONFIG"]
PROCESSORS = config_cluster_Dict["PROCESSORS"]
MEMORY = config_cluster_Dict["MEMORY"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#METADATA
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
metadata_Dict = utility.build_metadata_dict(METADATA_FILE)
design_Dict = build_design_Dict(metadata_Dict)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#UTILITIES
config_utilities_Dict = config["UTILITIES"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
config_qualimap_Dict = config["QUALIMAP"][GENOME]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#PEAK_CALLING
config_peak_calling_Dict = config["PEAK_CALLING"][GENOME]
MACS2_NARROW_PARAMETERS = config_peak_calling_Dict["MACS2_NARROW"]
MACS2_BROAD_PARAMETERS = config_peak_calling_Dict["MACS2_BROAD"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
config_qualimap_Dict = config["QUALIMAP"][GENOME]
config_reference_Dict = config["REFERENCE"][GENOME]
EFFECTIVE_GENOME_SIZE = config_reference_Dict["EFFECTIVE_GENOME_SIZE"]
# ------------------------------------
# ################################### WILDCARDS ################################ #

alignment_List = []
post_alignment_List = []
peak_calling_List = []
peak_analysis_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	#ALIGNMENT
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	#PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))


for design in design_Dict:
	#
	#ALIGNMENT
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	#
	#PEAK_ANALYSIS
	peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped}.narrowPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped}.broadPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"])))
	peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR}.narrowPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))
	peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR}.broadPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"])))

	for case in design_Dict[design]["Case"]:
		for control in design_Dict[design]["Control"]:
			#
			#PEAK_CALLING
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
	##
	for control in design_Dict[design]["Control"]:
		#
		#PEAK_CALLING
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}_VS_{control}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		#
		#PEAK_ANALYSIS
		peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{overlapped}_VS_{control}.narrowPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{overlapped}_VS_{control}.broadPeak.gz".format(design=design, overlapped="_OVERLAPPED_".join(design_Dict[design]["Case"]), control=control))
		peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{IDR}_VS_{control}.narrowPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
		peak_analysis_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{IDR}_VS_{control}.broadPeak.gz".format(design=design, IDR="_IDR_".join(design_Dict[design]["Case"]), control=control))
#
trackhub_List = [WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/trackhub/igv/{design}_igv_trackhub.xml".format(design=design)]
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		trackhub_List
# ################################### PIPELINE RULES ########################## #


rule IGV_Trackhub:
	"""
	"""
	input:
		alignment_List
		+ post_alignment_List
		+ peak_calling_List
		+ peak_analysis_List
	output:
		igv_trackhub_xml = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/trackhub/igv/{design}_igv_trackhub.xml",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "IGV_Trackhub: {wildcards.design}"
	run:
		IGV_String = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
		<Global name="''' + EXPERIMENT + '''"  infolink="''' + INFOLINK + '''" version="1">
		'''
		for design in design_Dict:
			#
			IGV_String += '''
				<Category name="''' + design + '''">
			'''
			##BAM
			IGV_String += '''
					<Category name="Bam">
			'''
			
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + design + "/alignment/*.bam", recursive=True):
				#
				accessible_each_track = HPC_DATASHARE + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			##PROCESSED_BAM
			IGV_String += '''
					<Category name="Processed_Bam">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + design + "/post_alignment/*.processed.bam", recursive=True):
				#
				accessible_each_track = HPC_DATASHARE + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			##PEAK_CALLING_NARROW
			IGV_String += '''
					<Category name="NarrowPeak_Bed">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/narrowpeak/*.narrowPeak.gz", recursive=True):
				#
				accessible_each_track = HPC_DATASHARE + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			IGV_String += '''
					<Category name="NarrowPeak_Bigwig">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/narrowpeak/*.narrowPeak.bigwig", recursive=True):
				#
				accessible_each_track = HPC_DATASHARE + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			##PEAK_CALLING_BROAD
			IGV_String += '''
					<Category name="BroadPeak_Bed">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/broadpeak/*.broadPeak.gz", recursive=True):
				#
				accessible_each_track = HPC_DATASHARE + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			IGV_String += '''
					<Category name="BroadPeak_Bigwig">
			'''
			for each_track in glob.glob(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/" + design + "/peak_calling/broadpeak/*.broadPeak.bigwig", recursive=True):
				#
				accessible_each_track = HPC_DATASHARE + each_track.split("RTB/datashare")[1]
				sample_name = accessible_each_track.split("/")[-1]
				IGV_String += '''
							<Resource name="''' + design + '_' + sample_name + '''" path="''' + accessible_each_track + '''"/>
				'''
			IGV_String += '''
					</Category>
			'''
			IGV_String += '''
				</Category>
			'''
		IGV_String += '''
		</Global>
		'''
		f = open(output.igv_trackhub_xml, "w")
		f.write(IGV_String)
		f.close()



		
















