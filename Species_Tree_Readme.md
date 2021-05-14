# Running the species tree pipeline

## Create a new environment
```
conda create --name phylo python=3.7
conda activate phylo
```
## Install packages
Most packages can be installed with conda
```
conda install -c bioconda exonerate iqtree amas newick_utils seqtk jalview

conda install -c smirarab treeshrink
```
Install upp, astral and clone the PhylogenomicsPipelines
```
mkdir ~/pasta-code; cd ~/pasta-code; git clone https://github.com/smirarab/pasta.git; git clone https://github.com/smirarab/sate-tools-linux.git; cd pasta; python setup.py develop --user
git clone https://github.com/smirarab/sepp.git; cd sepp; python setup.py config; python setup.py install; python setup.py upp

cd ~; wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.7.zip; unzip Astral.5.7.7.zip; export PATH=~/Astral/astral.5.7.7.jar:$PATH

cd ~; git clone https://github.com/baileyp1/PhylogenomicsPipelines; export PATH=~/PhylogenomicsPipelines:$PATH
```

## Run example
* 1. Activate environment and enter a working directory

```
conda activate phylo
cd /mnt/shared/scratch/kleempoe/paftol/phylotree/test_2
PATH=~/PhylogenomicsPipelines:$PATH
export ASTRAL=/home/kleempoe/Astral/astral.5.7.7.jar
PATH=/home/kleempoe/Astral/astral.5.7.7.jar:$PATH
```

* 2. Add all samples (.fasta files) in the in_fasta/ folder
* 3. the file in_fasta/sample_tree_tip_info_file.txt contains labels to rename tree tips once the species tree has been computed. It should be a space delimited file in the format "sample label" (example: 120340 Malpighiales_Euphorbiaceae_Euphorbia_canariensis)
* 4. Run the following command

```
make_species_trees_pipeline.sh \
-a \
-D 'dna' \
-A upp \
-t in_fasta/sample_tree_tip_info_file.txt \
-u \
-g ../Angiosperms353_targetSequences_20_geneIds_ONLY.txt \
-F '60 10' \
-m 0.7 \
-O 30 \
-K 0.003 \
-q iqtree2-B1000-nm1000 \
-T \
-L 30 \
-c 26 \
-C 4 \
-Q long \
-Y himem \
-R 5500 \
-U 700000 \
-V 0 \
-W 0 \
-H 1 \
-p run1 \
-X 5921,5596:6:40000  \
in_fasta/*.fasta \
> make_species_trees_pipeline.log 2>&1 &
```

* 5. If it runs correctly (check with squeue), The number of parallel jobs can be increased

```
scontrol update arraytaskthrottle=19 job=
```

The resulting species tree will be located at after_treeshrink_USE_THIS_dna/run1.dna.species_tree.astral_USE_THIS.nwk