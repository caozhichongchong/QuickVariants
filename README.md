# QuickVariants

## Fast and accurate genetic variant identification

QuickVariants summarizes allele information from read alignments without discarding or filtering the data.

If you're also interested in sequence alignment, you might be interested in [Mapper](https://github.com/mathjeff/mapper), which first aligns sequences and then identifies variants.

Read more about QuickVariants in [the paper](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-01891-4).

### Install
Requirement: [Java](https://www.java.com/en/download/help/download_options.html) \
Please download the latest QuickVariants here: https://github.com/caozhichongchong/QuickVariants/releases/download/1.1.0/quick-variants-1.1.0.jar

You may install java by `conda install conda-forge::openjdk`\
You may install java and QuickVariants by `conda install caozhichongchong::quick-variants`\
QuickVariants can be found at `$Conda_env_location/bin/quick-variants-VERSION.jar`
### Usage

```
java -Xms10g -Xmx10g -jar quick-variants-VERSION.jar [--out-vcf <out.vcf>] [--out-mutations <out.txt>] --reference <ref.fasta> --in-sam <input1.sam> [--in-sam <input2.sam> ...] --num-threads num_threads [options]
```
This command converts one or more SAM files to other formats, most notably .vcf.

**Input**
- `--reference <file>`: The reference to use. Should be in .fasta/.fa/.fna or .fastq/.fq format or a .gz of one of those formats.

**Output formats**\
Summary by reference position, mutation, genome, and raw output are possible.
- Options for output files, mutation counts, exclusion of non-mutations, and more.

**Summary by reference position**

- `--out-vcf <file>` output file to generate containing a description of mutation counts by position
-  --vcf-omit-support-reads By default, the vcf file has a column showing one or more supporting reads for each variant. If set, the output vcf file will hide the supporting reads for each variant.
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

See [TESTING.md](TESTING.md)

## If you're working on a bioinformatics project and would be interested in some consulting help, check out our website at https://omicode.info !
