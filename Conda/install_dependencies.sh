#!/bin/bash
## conda environment name

if [ $# -eq 0 ]; then
	echo "WORKDIR not provided"
	exit 1
fi

args=("$@")
WORKDIR=${args[0]}
mkdir -p $WORKDIR
mkdir -p ${WORKDIR}/EXECDIR
mkdir -p ${WORKDIR}/REFDIR
mkdir -p ${WORKDIR}/TMPDIR
EXECDIR=${WORKDIR}/EXECDIR
REFDIR=${WORKDIR}/REFDIR
TMPDIR=${WORKDIR}/TMPDIR

echo "EXECDIR: $EXECDIR"
echo "REFDIR: $REFDIR"

ENV_NAME_PY2=CHIP_Seq_py2
ENV_NAME_PY3=CHIP_Seq_py3


if [[ ! -f ${EXECDIR}/.conda_env ]]; then
	touch ${EXECDIR}/.conda_env
	chmod 766 ${EXECDIR}/.conda_env
fi

source ${EXECDIR}/.conda_env

if [[ ! -f $(which conda) ]]; then
	echo "Installing Conda ...."
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${TMPDIR}/Miniconda3-latest-Linux-x86_64.sh
	bash ${TMPDIR}/Miniconda3-latest-Linux-x86_64.sh -bfp ${EXECDIR}/Miniconda3
	echo "export PATH='${EXECDIR}/Miniconda3/bin:$PATH'" > ${EXECDIR}/.conda_env
	source ${EXECDIR}/.conda_env

fi

echo "Conda is Available "
conda --version

if [[ $(conda env list | grep $ENV_NAME_PY2 | wc -l) == "0" ]]; then
	echo "Installing $ENV_NAME_PY2 ...."
	conda create -n $ENV_NAME_PY2 python=2.7.15 -y -c defaults
	echo "$ENV_NAME_PY2 Installed!!!"
else
	echo "$ENV_NAME_PY2 is Available"
fi

if [[ $(conda env list | grep $ENV_NAME_PY3 | wc -l) == "0" ]]; then
	echo "Installing $ENV_NAME_PY3 ...."
	conda create -n $ENV_NAME_PY3 python=3.7.0 -y -c defaults
	echo "$ENV_NAME_PY3 Installed!!!"
else
	echo "$ENV_NAME_PY3 is Available"
fi


if [[ $(conda env list | grep $ENV_NAME_PY2 | wc -l) != "0" ]]; then
	source activate $ENV_NAME_PY2
	conda install -c defaults -c conda-forge bc=1.06.0 cython=0.28.5 boost=1.68.0 blas=1.1.0 openblas=0.3.3 ncurses=6.1.0 pandas=0.23.4 perl=5.26.2.0 libgfortran=3.0.0 pigz=2.3.4 zlib=1.2.11 graphviz=2.38.0 python-levenshtein=0.12.0 jinja2=2.10 gsl=2.4.0 imagemagick=7.0 matplotlib=2.2.0 -y
	conda install -c defaults -c r r=3.2 r-bitops=1.0 r-catools=1.17 -y
	conda install -c defaults -c bioconda r-snow=0.3 r-snowfall=1.84 bioconductor-rsamtools bioconductor-deseq2 r-spp=1.13 ghostscript=9.18.0 trim-galore=0.5.0 fastqc=0.11.7 afterqc=0.9.7 qualimap=2.2.2a multiqc=1.6.0 pysam=0.15.1 pybedtools=0.7.10 pyfaidx=0.5.3 -y
	conda install -c defaults -c anaconda wget=1.19.5 nomkl=3.0.0 numpy=1.14 numpy-base=1.14 scipy=1.1.0 six=1.11.0 python-dateutil=2.7.3 -y
	conda install -c defaults -c bioconda ucsc-bedgraphtobigwig=366 ucsc-fetchchromsizes=366 ucsc-wigtobigwig=366 ucsc-bigwiginfo=366 ucsc-bedclip=366 ucsc-bedtobigbed=366 ucsc-bigbedtobed=366 ucsc-bigwigtowig=366 ucsc-twobittofa=366 ucsc-fatotwobit=366 ucsc-bigwiginfo=366 ucsc-bigwigtobedgraph=366  ucsc-liftover=366 ucsc-bigWigAverageOverBed=366 -y
	conda install -c defaults -c bioconda ucsc-bigwigmerge=366 ucsc-bigbedinfo=366 ucsc-bigbedsummary=366 -y
	conda install -c defaults -c bioconda bedtools=2.27.1 macs2=2.1.1.20160309 cutadapt=1.18.0 bedops=2.4.35 bowtie=1.2.2 bowtie2=2.3.4.2 bwa=0.7.17 samtools=1.9.0 sambamba=0.6.6 picard=2.18.7 metaseq=0.5 preseq=2.0.3 -y
	
	if [[ ! -f ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/NucleoATAC/bin/nucleoatac ]]; then
		echo "Installing NUCLEOATAC..."
		cd ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}
		git clone https://github.com/GreenleafLab/NucleoATAC.git
		cd NucleoATAC
		pip install --upgrade pip
		pip install mmtf-python
		pip install msgpack
		pip install .
		chmod -R 766 ./bin/*
		echo "export PATH='${PWD}/bin:$PATH'" >> ${EXECDIR}/.conda_env
	else
		echo "NucleoATAC is Available"
	fi


	if [[ ! -f ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/script/run_spp_nodups.R ]]; then
		mkdir -p ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/script
		wget https://raw.githubusercontent.com/shenlab-sinai/chip-seq_preprocess/master/bin/run_spp_nodups.R -O ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/script/run_spp_nodups.R
		wget https://raw.githubusercontent.com/amirshams84/ATAC_SEQ/master/script/narrowpeak2hammock.py -O ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/script/narrowpeak2hammock.py
		wget https://raw.githubusercontent.com/amirshams84/ATAC_SEQ/master/script/broadpeak2hammock.py -O ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/script/broadpeak2hammock.py
		chmod -R 766 ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/script/*
		echo "export PATH='${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY2}/script/bin:$PATH'" >> ${EXECDIR}/.conda_env
	else
		echo "Customized script is available"
	fi
	source deactivate

else
	echo "${ENV_NAME_PY2} packages is Available"
fi




if [[ $(conda env list | grep $ENV_NAME_PY3 | wc -l) != "0" ]]; then
	source activate $ENV_NAME_PY3
	conda install -c defaults -c conda-forge bc=1.06.0 cython=0.28.5 boost=1.68.0 blas=1.1.0 openblas=0.3.3 ncurses=6.1.0 pandas=0.23.4 perl=5.26.2.0 libgfortran=3.0.0 pigz=2.3.4 zlib=1.2.11 graphviz=2.38.0 python-levenshtein=0.12.0 jinja2=2.10 gsl=2.4.0 imagemagick=7.0 matplotlib=2.2.0 -y
	conda install -c defaults -c anaconda wget=1.19.5 nomkl=3.0.0 numpy=1.14 numpy-base=1.14 scipy=1.1.0 six=1.11.0 python-dateutil=2.7.3 -y
	conda install -c defaults -c bioconda pysam=0.15.1 pybedtools=0.7.10 pyfaidx=0.5.3 -y
	conda install -c defaults -c bioconda perl-threaded=5.22.0 homer=4.9.1-6 meme=5.0.2 bedtools=2.27.1 idr=2.0.4.2 -y


	if [[ ! -f ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/data/genomes/hg19/chrom.sizes ]]; then
		perl ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/configureHomer.pl -install hg19
	fi
	if [[ ! -f ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/data/genomes/hg38/chrom.sizes ]]; then
		perl ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/configureHomer.pl -install hg38
	fi
	if [[ ! -f ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/data/genomes/mm9/chrom.sizes ]]; then
		perl ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/configureHomer.pl -install mm9
	fi
	if [[ ! -f ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/data/genomes/mm10/chrom.sizes ]]; then
		perl ${EXECDIR}/Miniconda3/envs/${ENV_NAME_PY3}/share/homer-4.9.1-6/configureHomer.pl -install mm10
	fi
	
	
	source deactivate
else
	echo "$ENV_NAME_PY3 installed successfully"
fi

echo "Installing dependencies has been successfully done."


