# QuickVariants: Fast and accurate genetic variant identification

If you're also interested in sequence alignment, you might be interested in [Mapper](https://github.com/mathjeff/mapper), which first aligns sequences and then identifies variants.

Read more about QuickVariants in [the paper](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-01891-4).

Download the latest release version here: https://github.com/caozhichongchong/QuickVariants/releases/download/1.2.4/quick-variants-1.2.4.jar

Also available as `quick-variants` in Bioconda - see https://bioconda.github.io/recipes/quick-variants/README.html

Contact:\
 Dr. Anni Zhang, caozhichongchong at gmail dot com

## Usage:
  java -jar quick-variants.jar [--out-vcf <out.vcf>] [--out-mutations <out.txt>] [--out-sam <out.sam>] [--out-refs-map-count <counts.txt>] --reference <ref.fasta> --in-sam <input.sam> [options]

    Converts a sam file to other formats, most notably .vcf

  INPUT:

    --reference <file> the reference to use. Should be in .fasta/.fa/.fna or .fastq/.fq format or a .gz of one of those formats.

    --in-ordered-sam <file.sam> the input alignments. If there are Illumina-style paired-end reads, alignments for each mate should be adjacent in the sam file
    --in-unordered-sam <file.sam> same as --in-ordered-sam but doesn't require alignments for matest to be listed in adjacent lines in the file
      This argument is slower than --in-unordered-sam and requires more memory
    --in-sam <file.sam> alias for --in-unordered-sam <file.sam>

  OUTPUT FORMATS:

    Summary by reference position

      --out-vcf <file> output file to generate containing a description of mutation counts by position
        Details about the file format are included in the top of the file
      --vcf-exclude-non-mutations if set, the output vcf file will exclude positions where no mutations were detected
      --vcf-omit-support-reads By default, the vcf file has a column showing one or more supporting reads for each variant. If set, the output vcf file will hide the supporting reads for each variant.
      --distinguish-query-ends <fraction> (default 0.1) In the output vcf file, we separately display which queries aligned at each position with <fraction> of the end of the query and which didn't.

      --snp-threshold <min total depth> <min supporting depth fraction> (default 0, 0)
        The minimum total depth and minimum supporting depth fraction required at a position to report the support for the mutation

      --indel-start-threshold <min total depth> <min supporting depth fraction> (default 0, 0)
        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report support for the start of an insertion or deletion

      --indel-continue-threshold <min total depth> <min supporting depth fraction> (default 0, 0)
        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report support for a continuation of an insertion or deletion
      --indel-threshold <min total depth> <min supporting depth fraction>
        Alias for --indel-start-threshold <min total depth> <min supporting depth frequency> and --indel-continue-threshold <min total depth> <min supporting depth frequency>

    Summary by mutation

      --out-mutations <file> output file to generate listing the mutations of the queries compared to the reference genome
        Details about the file format are included in the top of the file

      --distinguish-query-ends <fraction> (default 0.1) When detecting indels, only consider the middle <fraction> of each query

      --snp-threshold <min total depth> <min supporting depth fraction> (default 5, 0.9)
        The minimum total depth and minimum supporting depth fraction required at a position to report it as a point mutation

      --indel-start-threshold <min total depth> <min supporting depth fraction> (default 1, 0.8)
        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report it as the start of an insertion or deletion

      --indel-continue-threshold <min total depth> <min supporting depth fraction> (default 1, 0.7)
        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report it as a continuation of an insertion or deletion
      --indel-threshold <min total depth> <min supporting depth fraction>
        Alias for --indel-start-threshold <min total depth> <min supporting depth frequency> and --indel-continue-threshold <min total depth> <min supporting depth frequency>

    Summary by genome

      --out-refs-map-count <file> the output file to summarize the number of reads mapped to each combination of references

    Raw output

      --out-sam <file> the output file in SAM format

    --no-output if no output is requested, skip writing output rather than throwing an error

    Debugging

      -v, --verbose output diagnostic information

      -vv more diagnostic information

      --verbosity-auto set verbosity flags dynamically and automatically based on the data and alignment

  Multiple output formats may be specified during a single run; for example:

    --out-mutations out.mutations --out-vcf out.vcf

  OTHER:

    --num-threads <count> number of threads to use at once for processing. Higher values will run more quickly on a system that has that many CPUs available.

    --help output this help message
      If no other arguments are given, exit instead of attempting to summarize alignments

    --version output the version of QuickVariants
      If no other arguments are given, exit instead of attempting to summarize alignments


### Test

See [TESTING.md](TESTING.md)

## If you're working on a bioinformatics project and would be interested in some consulting help, check out our website at https://genomiverse.net/ !
