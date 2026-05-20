
## Additional documentation: command-line code for running a untargeted viromics analysis. 

This document outlines a step-by-step workflow for predicting viruses in metagenomic samples using a combination of binning, phage prediction, and quality assessment tools. Conda/Mamba installation instructions are provided for each tool.

### **Setting variables and creating directories**
```markdown
ASSEMBLY="~/workshop_materials/untargetedViromics/uv_data/assembly.fasta"
mkdir ~/workshop_materials/untargetedViromics/runningTools
PHAGEPRED_WD="~/workshop_materials/untargetedViromics/runningTools/phagePrediction"
QUALPRED_WD="~/workshop_materials/untargetedViromics/runningTools/Qual"
CPUS_PER_TASK=4
cd ~/workshop_materials/untargetedViromics/runningTools
```
### Binning using Reneo
**Do not run**

#### **Description**
Reneo is a binning tool that groups contigs into bins (putative genomes) based on coverage, sequence composition, and a flow decomposition algorithm that allows us to get complete genomes from our contigs. This step is important for reducing complexity and isolating viral contigs.

#### **Conda Installation**
```bash
conda create -n reneo -c bioconda reneo
conda activate reneo
```

#### **Command**
```{.bash}
reneo run --input "${ASSEMBLY}" \
          --reads "${READS}" \
          --minlength 1000 \
          --output "${RENEO_OUT}" \
          --threads ${CPUS_PER_TASK}
echo "reneo done"
```
### Phage Prediction using Jaeger, geNomad, and DeepVirFinder

#### **Description**
#### Jaeger

**Methodology:** Jaeger uses a machine learning model based on convolutional neural networks (CNNs) to predict phage sequences. It analyzes genomic features such as k-mer frequencies, GC content, and sequence length.
**Strengths:**
High sensitivity for detecting novel phages due to its ability to learn complex patterns in viral sequences.
Can process both reads and contigs, making it versatile for different types of metagenomic data.
Performs well on low-abundance sequences and fragmented genomes.

#### geNomad
**Methodology:** geNomad combines homology-based searches (e.g., against protein domain databases like Pfam) with neural network models to classify sequences as viral, plasmid, or chromosomal.
The homology component identifies conserved viral and plasmid proteins.
The machine learning component uses k-mer frequencies and genomic features to refine predictions.
**Strengths:**
High accuracy in distinguishing between viruses, plasmids, and host sequences.
Effective for novel sequences due to its hybrid approach.
Provides detailed annotations, including viral taxonomy and functional genes.

#### DeepVirFinder
**Methodology:** DeepVirFinder employs a deep learning model based on convolutional neural networks (CNNs) to predict viral sequences. It uses k-mer frequencies and sequence composition as input features.
**Strengths:**
High sensitivity for detecting novel viruses, even in the absence of close homologs in reference databases.
Works well on short reads and fragmented contigs.
Can be applied to both DNA and RNA viruses.

#### VirSorter2
**Methodology:** VirSorter2 combines homology-based searches (using curated viral protein databases) with machine learning models to identify phage sequences. It uses genomic features such as gene density, strand shifts, and viral hallmark genes.
**Strengths:**
High accuracy in detecting phages and viral contigs in metagenomic assemblies.
Can classify sequences into lytic, lysogenic, or eukaryotic viruses.
Includes a curated database of viral proteins for improved homology-based detection.

#### **Conda Installation**
```{.bash}
## Jaeger
mamba create -n jaeger -c bioconda jaeger
conda activate jaeger

## geNomad **Do not run**
mamba create -n genomad -c bioconda genomad
conda activate genomad

## virsorter **Do not run**
mamba create -n vs2 -c conda-forge -c bioconda virsorter=2
virsorter setup -d db -j 4

## DeepVirFinder
mamba create --name dvf python=3.6 numpy theano=1.0.3 keras=2.2.4 scikit-learn Biopython h5py=2.10.0
git clone https://github.com/jessieren/DeepVirFinder
cd DeepVirFinder
```

#### **Commands**
```{.bash}
mkdir -p "${PHAGEPRED_WD}"
## DeepVirFinder
conda activate dvf
python dvf.py -i "${ASSEMBLY}" \
              -o "${PHAGEPRED_WD}/deepvirfinder" \
              -l 1000 \
              -c ${CPUS_PER_TASK}

echo "deepvirfinder done"

## Jaeger
conda activate jaeger
Jaeger -i "${ASSEMBLY}"  \
       -o "${PHAGEPRED_WD}/jaeger" \
       -s 2.5 \
       --fsize 1000 \
       --stride 1000

echo "jaeger done"

## geNomad **Do not run**
conda activate genomad
genomad end-to-end --min-score 0.6 \
       --cleanup \
       --threads ${CPUS_PER_TASK} \
       "${ASSEMBLY}" \
       "${PHAGEPRED_WD}/geNomad" \
       "${GENOMAD_DB}"

echo "genomad done"

## virsorter **Do not run**
virsorter run -w "${PHAGEPRED_WD}/virsorter" \
-i "${ASSEMBLY}" \
--include-groups "dsDNAphage,ssDNA, RNA" \
-j 4 \
--min-length 1000 \
--min-score 0.8 \
--provirus-off  \
all
```
### Integrative workflows for phage prediction, taxonomic classification, lifestyle prediction, and host prediction
#### **Conda Installation**
```{.bash}
mamba create -n phabox2 phabox=2.1.10 -c conda-forge -c bioconda -y

## Downloading the database using wget
cd ~/workshop_materials/untargetedViromics/runningTools
wget https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2.zip
unzip phabox_db_v2.zip > /dev/null
```
#### **Command**
```bash
conda activate phabox2
phabox2 --task end_to_end --dbdir phabox_db_v2 \
        --outpth  ${PHAGEPRED_WD}/Phabox_OUT \
        --contigs ${ASSEMBLY} \
        --len 1000 \
        --threads ${CPUS_PER_TASK}
echo "phabox2 done"
```
### **Step 4: Quality Assessment using CheckV**

#### **Description**
CheckV is a tool designed to assess the quality and completeness of viral genomes recovered from metagenomes. It employs a multi-step approach:
**Completeness Estimation:**
Uses a database of high-quality reference viral genomes to identify conserved, single-copy genes (e.g., capsid proteins, terminases, integrases) that serve as markers for viral completeness.
Estimates completeness by detecting the presence/absence of these hallmark genes and their collinearity (gene order conservation).
**Contamination Detection:**
Identifies host contamination (e.g., bacterial or eukaryotic genes) using a database of non-viral sequences.
Flags sequences with atypical GC content, codon usage, or gene content for further scrutiny.

#### **Conda Installation**

```{.bash}
mamba create -n checkv -c conda-forge -c bioconda checkv
conda activate checkv
checkv download_database ./
```

#### **Command**
```{.bash}
mkdir -p "${QUALPRED_WD}"

checkv end_to_end "${ASSEMBLY}" \
       "${QUALPRED_WD}" \
       -d checkv-db-v1.5 \
       -t ${CPUS_PER_TASK}

echo "checkV done"
```
### Phold for Phage Annotation
**Do not Run**

#### **Description**
Phold is a protein structure prediction and annotation tool designed specifically for bacteriophages. It combines homology-based searches with protein language models to annotate phage genomes. Annotation is performed by aligning the predicted structures to a database of structures usinf foldseek. 

#### **Conda Installation**

```{.bash}
mamba create -n pholdENV -c bioconda phold
conda activate pholdENV
phold install #Installs the databases
```

#### **Command**
```{.bash}
phold run -i "${ASSEMBLY}" \
       -o "${PHAGEPRED_WD}/phold" \
       -d "${PHOLD_DB}" \
       -t ${CPUS_PER_TASK} --cpu
echo "phold done"
``