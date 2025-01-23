# QuickVariants

### Code for paper "Fast and Accurate Variant Identification Tool for Sequencing-Based Studies"

QuickVariants is a fast and accurate variant identification tool, designed to summarize allele information from read alignments without discarding or filtering the data.

### Install
Requirement: [Java](https://www.java.com/en/download/help/download_options.html) \
Please download the latest QuickVariants here: https://github.com/caozhichongchong/QuickVariants/releases/download/1.1.0/quick-variants-1.1.0.jar

You may install java by `conda install conda-forge::openjdk`\
You may install java and QuickVariants by `conda install caozhichongchong::quick-variants`\
QuickVariants can be found at `$Conda_env_location/bin/quick-variants-VERSION.jar`
### Usage

```
java -Xms10g -Xmx10g -jar quick-variants-VERSION.jar [--out-vcf <out.vcf>] [--out-mutations <out.txt>] --reference <ref.fasta> --in-sam <input.sam> --num-threads num_threads [options]
```
This command converts a SAM file to other formats, most notably .vcf.

**Input**
- `--reference <file>`: The reference to use. Should be in .fasta/.fa/.fna or .fastq/.fq format or a .gz of one of those formats.

**Output formats**\
Summary by reference position, mutation, genome, and raw output are possible.
- Options for output files, mutation counts, exclusion of non-mutations, and more.

**Summary by reference position**

- `--out-vcf <file>` output file to generate containing a description of mutation counts by position
- `--vcf-exclude-non-mutations` if set, the output vcf file will exclude positions where no mutations were detected
- `--distinguish-query-ends <fraction>` (default 0.1) In the output vcf file, we separately display which queries aligned at each position with <fraction> of the end of the query and which didn't.

**Summary by mutation**

- `--out-mutations <file>` output file to generate listing the mutations of the queries compared to the reference genome

- `--distinguish-query-ends <fraction>` (default 0.1) When detecting indels, only consider the middle <fraction> of each query

- `--snp-threshold <min total depth> <min supporting depth fraction>` (default 5, 0.9)
    The minimum total depth and minimum supporting depth fraction required at a position to report it as a point mutation

- `--indel-start-threshold <min total depth> <min supporting depth fraction>` (default 1, 0.8)
    The minimum total (middle) depth and minimum supporting depth fraction required at a position to report it as the start of an insertion or deletion

- `--indel-continue-threshold <min total depth> <min supporting depth fraction>` (default 1, 0.7)
    The minimum total (middle) depth and minimum supporting depth fraction required at a position to report it as a continuation of an insertion or deletion
- `--indel-threshold <min total depth> <min supporting depth fraction>`
    Alias for --indel-start-threshold <min total depth> <min supporting depth frequency> and --indel-continue-threshold <min total depth> <min supporting depth frequency>

**Summary by genome**

- `--out-refs-map-count <file>` the output file to summarize the number of reads mapped to each combination of references

**Raw output**

- `--out-sam <file>` the output file in SAM format

- `--out-unaligned <file>` output file containing unaligned reads. Must have a .fasta or .fastq extension

- `--no-output` if no output is requested, skip writing output rather than throwing an error

**Debugging**

- `-v, --verbose` output diagnostic information

- `-vv` more diagnostic information

- `--verbosity-auto` set verbosity flags dynamically and automatically based on the data and alignment

Multiple output formats may be specified during a single run; for example:

- `--out-sam` out.sam --out-unaligned out.fastq

**Others**

- `--num-threads <count>` number of threads to use at once for processing. Higher values will run more quickly on a system that has that many CPUs available.


### Test
```
java -jar quick-variants-1.1.0.jar --out-vcf Fig4Example1.vcf --reference examples/Fig4/reference.fasta --in-sam examples/Fig4/Example1/90.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig4Example2.vcf --reference examples/Fig4/reference.fasta --in-sam examples/Fig4/Example2/86.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig5Example1.vcf --reference examples/Fig5/reference.fasta --in-sam examples/Fig5/Example1/99.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig5Example2.vcf --reference examples/Fig5/reference.fasta --in-sam examples/Fig5/Example2/102.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig5Example3.vcf --reference examples/Fig5/reference.fasta --in-sam examples/Fig5/Example3/9.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig5Example4.vcf --reference examples/Fig5/reference.fasta --in-sam examples/Fig5/Example4/47.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig6Example1.vcf --reference examples/Fig6/reference.fasta --in-sam examples/Fig6/Example1/88.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig6Example2.vcf --reference examples/Fig6/reference.fasta --in-sam examples/Fig6/Example2/90.sam
java -jar quick-variants-1.1.0.jar --out-vcf Fig6Example3.vcf --reference examples/Fig6/reference.fasta --in-sam examples/Fig6/Example3/9.sam
```
- You can compare your results to pre-generated VCF files located in example/test_results/.
### filtering point mutations and indels
```
python QuickVariants_Pointmutationfilter.py -i Your_QuickVariants_Output_Folder
python QuickVariants_Indelfilter.py -i Your_QuickVariants_Output_Folder
```
- `af`, `allele-frequency` minimum allele frequency
- `rd`, `read-depth` minimum read depth
- Testing
```
unzip test_vcf.zip
python QuickVariants_Pointmutationfilter.py
python QuickVariants_Indelfilter.py
```

### Additional scripts and models
The `benchmark_scripts` folder contains code used to construct the benchmark dataset, filter SNPs and indels in VCF files, and analyze VCF files.

- `SNP_model.py`: Insert in silico mutations and indels randomly into reference genomes, and generate alignment code.
- `SNPfilter.py`: Filter SNPs in VCF files generated by QuickVariants and bcftools.
- `Indelfiter.py`: Filter indels in VCF files generated by QuickVariants and bcftools, and summarize the number of false positives and false negatives for each dataset.
- `SNP_model_compare.py`: Summarize the number of false positive and false negative SNPs for each dataset.
- `VCFAnalysis.ipynb`: Code to analyze point mutations and indels detected from benchmark datasets.
- `COVID_MG.ipynb`: Code to analyze point mutations and indels detected from SARS-COV-2 sewage MG real data.

**Analyzing benchmark datasets used in this study**\
Please download benchmark datasets [here](https://doi.org/10.6084/m9.figshare.25437217)\
Requirements: [bowtie2](https://anaconda.org/bioconda/bowtie2), [bwa](https://anaconda.org/bioconda/bwa), [minimap2](https://anaconda.org/bioconda/minimap2),
[samtools](https://www.htslib.org/download/), [bcftools](https://samtools.github.io/bcftools/howtos/install.html), [python3](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-python.html), [jupyter notebook](https://jupyter.org/install)\
-Gut microbiome WGS data with in silico mutations
```
python SNP_model_covid.py -i Gut_microbiome_benchmark/original_data -o Gut_microbiome_benchmark/
python SNPfilter.py -i Gut_microbiome_benchmark/
python Indelfilter.py -i Gut_microbiome_benchmark/
python SNP_model_compare.py -i Gut_microbiome_benchmark/
```

Point mutations detected: Gut_microbiome_benchmark/SNP_model/merge/*final.txt and Gut_microbiome_benchmark/SNP_model/merge/model.sum.txt\
Indels detected: Gut_microbiome_benchmark/SNP_model/merge/*indel.vcf.filtered and Gut_microbiome_benchmark/SNP_model/merge/modelindelsum.txt

-SARS-COV-2 WGS data with in silico mutations
```
python SNP_model_covid.py -i COVID_benchmark/original_data -fa .fasta -fq _1.fastq -o COVID_benchmark/
python SNPfilter.py -i COVID_benchmark/
python Indelfilter.py -i COVID_benchmark/
python SNP_model_compare.py -i COVID_benchmark/
```

-WGS data simulated with sequencing errors
```
python SNP_model_covid.py -i WGS_simulation_sequencingerror/original_data -fa .fasta -fq _1.fq -o WGS_simulation_sequencingerror/
python SNPfilter.py -i WGS_simulation_sequencingerror/
python Indelfilter.py -i WGS_simulation_sequencingerror/
python SNP_model_compare.py -i WGS_simulation_sequencingerror/
```

-MG data simulated with sequencing errors (20X)
```
python SNP_model_covid.py -i MG_simulation_sequencingerror/original_data -fa .fasta -fq _1.fq -o MG_simulation_sequencingerror/
python SNPfilter.py -i MG_simulation_sequencingerror/
python Indelfilter.py -i MG_simulation_sequencingerror/
python SNP_model_compare.py -i MG_simulation_sequencingerror/
```

-MG data simulated with sequencing errors (100X)
```
python SNP_model_covid.py -i MGBIG_simulation_sequencingerror/original_data -fa .fasta -fq _1.fq -o MGBIG_simulation_sequencingerror/
python SNPfilter.py -i MGBIG_simulation_sequencingerror/
python Indelfilter.py -i MGBIG_simulation_sequencingerror/
python SNP_model_compare.py -i MGBIG_simulation_sequencingerror/
```

-SARS-COV-2 sewage MG real data
```
python SNP_model_covid.py -i COVID_MGSW/original_data -fa .fasta -fq _1.fq -o COVID_MGSW/
python SNPfilter.py -i COVID_MGSW/
python Indelfilter.py -i COVID_MGSW/
python SNP_model_compare.py -i COVID_MGSW/
```
