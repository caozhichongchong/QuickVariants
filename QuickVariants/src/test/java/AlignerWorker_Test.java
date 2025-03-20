package mapper;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;

public class AlignerWorker_Test {

  @Test
  public void simpleTest() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t4M\t*\t*\t4\tACGT\t*";

    String ref  = "ACGTAAAAACGTAAAA";

    String vcf = buildVcf(sam1, ">contig1\n" + ref);

    String expectedVcf = "contig1	1	A		1	1,0	0,0	\n" +
       "contig1	2	C		1	1,0	0,0	\n" +
       "contig1	3	G		1	1,0	0,0	\n" +
       "contig1	4	T		1	1,0	0,0	\n";

    checkVcf(vcf, expectedVcf);
  }


  @Test
  public void pairedEndAlignment() {
    String sam1 = "name1\t32\tcontig1\t1\t255\t10M\tcontig1\t20\t10\tAACCGGTTAT\t*";
    String sam2 = "name1\t16\tcontig1\t21\t255\t10M\tcontig1\t1\t10\tACGTACGTAT\t*";
    String ref  = "AACCGGTTATAAAAAAAAAAACGTACGTATAAAAAAAAAA";

    String vcf = buildVcf(sam1 + "\n" + sam2, ">contig1\n" + ref);

    String expectedVcf = "contig1	1	A		1	1,0	0,0	\n" +
        "contig1	2	A		1	1,0	0,0	\n" +
        "contig1	3	C		1	1,0	0,0	\n" +
        "contig1	4	C		1	1,0	0,0	\n" +
        "contig1	5	G		1	1,0	0,0	\n" +
        "contig1	6	G		1	1,0	0,0	\n" +
        "contig1	7	T		1	1,0	0,0	\n" +
        "contig1	8	T		1	1,0	0,0	\n" +
        "contig1	9	A		1	1,0	0,0	\n" +
        "contig1	10	T		1	1,0	0,0	\n" +
        "contig1	21	A		1	0,1	0,0	\n" +
        "contig1	22	C		1	0,1	0,0	\n" +
        "contig1	23	G		1	0,1	0,0	\n" +
        "contig1	24	T		1	0,1	0,0	\n" +
        "contig1	25	A		1	0,1	0,0	\n" +
        "contig1	26	C		1	0,1	0,0	\n" +
        "contig1	27	G		1	0,1	0,0	\n" +
        "contig1	28	T		1	0,1	0,0	\n" +
        "contig1	29	A		1	0,1	0,0	\n" +
        "contig1	30	T		1	0,1	0,0	\n";

    checkVcf(vcf, expectedVcf);
  }

  @Test
  public void oneReadWithMultipleAlignments() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t4M\t*\t*\t4\tACGT\t*";
    String sam2 = "name1\t0\tcontig1\t9\t255\t4M\t*\t*\t4\tACGT\t*";

    String ref  = "ACGTAAAAACGTAAAA";

    String vcf = buildVcf(sam1 + "\n" + sam2, ">contig1\n" + ref);

    String expectedVcf = "contig1	1	A		0.5	0.5,0	0,0	\n" +
        "contig1	2	C		0.5	0.5,0	0,0	\n" +
        "contig1	3	G		0.5	0.5,0	0,0	\n" +
        "contig1	4	T		0.5	0.5,0	0,0	\n" +
        "contig1	9	A		0.5	0.5,0	0,0	\n" +
        "contig1	10	C		0.5	0.5,0	0,0	\n" +
        "contig1	11	G		0.5	0.5,0	0,0	\n" +
        "contig1	12	T		0.5	0.5,0	0,0	\n" +
        "";

    checkVcf(vcf, expectedVcf);
  }

  @Test
  public void pairedEndReadWithMultipleAlignments() {
    String samA1 = "name1\t32\tcontig1\t1\t255\t4M\tcontig1\t9\t4\tACGT\t*";
    String samA2 = "name1\t16\tcontig1\t9\t255\t4M\tcontig1\t1\t4\tCCCC\t*";
    String samB1 = "name1\t32\tcontig1\t17\t255\t4M\tcontig1\t25\t4\tACGT\t*";
    String samB2 = "name1\t16\tcontig1\t25\t255\t4M\tcontig1\t17\t4\tCCCC\t*";

    String ref  = "ACGTAAAACCCCTTTTACGTAAAACCCC";

    String vcf = buildVcf(samA1 + "\n" + samA2 + "\n" + samB1 + "\n" + samB2, ">contig1\n" + ref);

    String expectedVcf = "contig1	1	A		0.5	0.5,0	0,0	\n" +
        "contig1	2	C		0.5	0.5,0	0,0	\n" +
        "contig1	3	G		0.5	0.5,0	0,0	\n" +
        "contig1	4	T		0.5	0.5,0	0,0	\n" +
        "contig1	9	C		0.5	0,0.5	0,0	\n" +
        "contig1	10	C		0.5	0,0.5	0,0	\n" +
        "contig1	11	C		0.5	0,0.5	0,0	\n" +
        "contig1	12	C		0.5	0,0.5	0,0	\n" +
        "contig1	17	A		0.5	0.5,0	0,0	\n" +
        "contig1	18	C		0.5	0.5,0	0,0	\n" +
        "contig1	19	G		0.5	0.5,0	0,0	\n" +
        "contig1	20	T		0.5	0.5,0	0,0	\n" +
        "contig1	25	C		0.5	0,0.5	0,0	\n" +
        "contig1	26	C		0.5	0,0.5	0,0	\n" +
        "contig1	27	C		0.5	0,0.5	0,0	\n" +
        "contig1	28	C		0.5	0,0.5	0,0	\n" +
        "";

    checkVcf(vcf, expectedVcf);
  }

  @Test
  public void testMissingMate() {
    String samA1 = "name1\t32\tcontig1\t1\t255\t4M\tcontig1\t9\t4\tACGT\t*";
    String samA2 = "name1\t16\tcontig1\t9\t255\t4M\tcontig1\t1\t4\tCCCC\t*";
    String samB1 = "name2\t32\tcontig1\t17\t255\t4M\tcontig1\t25\t4\tACGT\t*";
    String samB2 = "name2\t16\tcontig1\t25\t255\t4M\tcontig1\t17\t4\tCCCC\t*";
    String samC1 = "name3\t32\tcontig1\t33\t255\t4M\tcontig1\t41\t4\tACGT\t*";
    String samC2 = "name3\t16\tcontig1\t41\t255\t4M\tcontig1\t33\t4\tCCCC\t*";
    String goodSam = samA1 + "\n" + samA2;
    String badSam = samA1 + "\n" + samB1 + "\n" + samB2 + "\n" + samC1 + "\n" + samC2;
    String ref = "ACGTAAAACCCCTTTTACGTAAAACCCCTTTTACGTAAAACCCCTTTTACGTAAAACCCCTTTT";

    int goodSamNumErrors = getNumErrors(goodSam);
    if (goodSamNumErrors != 0) {
      fail("Found " + goodSamNumErrors + " errors parsing '" + goodSam + "'");
    }
    int badSamNumErrors = getNumErrors(badSam);
    int badSamExpectedNumErrors = 1;
    if (badSamNumErrors != badSamExpectedNumErrors) {
      fail("Found " + badSamNumErrors + " errors instead of " + badSamExpectedNumErrors + " parsing:\n'" + badSam + "'");
    }
    String vcfFromBadSam = buildVcf(badSam, ">contig1\n" + ref);
    String expectedVcf =
        "contig1	1	A		1	1,0	0,0	\n" +
        "contig1	2	C		1	1,0	0,0	\n" +
        "contig1	3	G		1	1,0	0,0	\n" +
        "contig1	4	T		1	1,0	0,0	\n" +
        "contig1	17	A		1	1,0	0,0	\n" +
        "contig1	18	C		1	1,0	0,0	\n" +
        "contig1	19	G		1	1,0	0,0	\n" +
        "contig1	20	T		1	1,0	0,0	\n" +
        "contig1	25	C		1	0,1	0,0	\n" +
        "contig1	26	C		1	0,1	0,0	\n" +
        "contig1	27	C		1	0,1	0,0	\n" +
        "contig1	28	C		1	0,1	0,0	\n" +
        "contig1	33	A		1	1,0	0,0	\n" +
        "contig1	34	C		1	1,0	0,0	\n" +
        "contig1	35	G		1	1,0	0,0	\n" +
        "contig1	36	T		1	1,0	0,0	\n" +
        "contig1	41	C		1	0,1	0,0	\n" +
        "contig1	42	C		1	0,1	0,0	\n" +
        "contig1	43	C		1	0,1	0,0	\n" +
        "contig1	44	C		1	0,1	0,0	\n" +
        "";

    checkVcf(vcfFromBadSam, expectedVcf);
  }

  @Test
  public void testSupplementaryAlignment() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t4M\t*\t*\t*\tACGT\t*\tSA:Z:contig1,11,+,4M,AAAA,0;";
    String sam2 = "name1\t0\tcontig1\t11\t255\t4M\t*\t*\t*\tAAAA\t*\tSA:Z:contig1,1,+,4M,ACGT,0;";

    String ref  = "ACGTGGGGCCAAAACCCC";

    String vcf = buildVcf(sam1 + "\n" + sam2, ">contig1\n" + ref);

    String expectedVcf =
        "contig1	1	A		1	1,0	0,0	\n" +
        "contig1	2	C		1	1,0	0,0	\n" +
        "contig1	3	G		1	1,0	0,0	\n" +
        "contig1	4	T		1	1,0	0,0	\n" +
        "contig1	11	A		1	1,0	0,0	\n" +
        "contig1	12	A		1	1,0	0,0	\n" +
        "contig1	13	A		1	1,0	0,0	\n" +
        "contig1	14	A		1	1,0	0,0	\n" +
        "";

    checkVcf(vcf, expectedVcf);
  }

  private int getNumErrors(String samLines) {
    GroupedAlignment_Provider samParser = newSamParser(samLines);
    while (samParser.getNextGroup() != null) {
    }
    return samParser.getNumErrors();
  }

  private void checkVcf(String actual, String expected) {
    actual = withoutVcfMetadataLines(actual);
    if (!(expected.equals(actual))) {
      String actualAsCode = "\"" + actual.replace("\n", "\\n\" +\n        \"") + "\";";
      fail("Difference in generated .vcf file.\nactual vcf:\n" + actual + "\nexpected vcf:\n" + expected + "\n code for actual vcf:\n" + actualAsCode);
    }
  }

  private String withoutVcfMetadataLines(String vcf) {
    String[] lines = vcf.split("\n");
    StringBuilder resultBuilder = new StringBuilder();
    for (String line: lines) {
      if (!line.startsWith("#") && line.length() > 0) {
        resultBuilder.append(line);
        resultBuilder.append("\n");
      }
    }
    return resultBuilder.toString();
  }

  private GroupedAlignment_Provider newSamParser(String samRecords) {
    // parse alignments
    StringReader alignmentStringReader = new StringReader(samRecords);
    SamReader alignmentReader = new SamReader(new BufferedReader(alignmentStringReader), "alignments.sam");
    SamAlignment_Provider queryBuilder = new SamAlignment_Provider(alignmentReader, "alignments.sam", false);
    GroupedAlignment_Provider samParser = new Simple_GroupedAlignment_Provider(queryBuilder);
    return samParser;
  }

  private String buildVcf(String samRecords, String referenceGenome) {
    // make alignment listener
    MatchDatabase matchDatabase = new MatchDatabase(0);
    List<AlignmentListener> listeners = new ArrayList<AlignmentListener>();
    listeners.add(matchDatabase);

    // parse alignments
    GroupedAlignment_Provider samParser = newSamParser(samRecords);
    List<List<SamAlignment_Builder>> alignments = getQueries(samParser);

    // parse reference
    BufferedReader referenceReader = new BufferedReader(new StringReader(referenceGenome));
    SequenceProvider referenceParser = new FastaParser(referenceReader, "reference.fasta");
    Sequence referenceContig = referenceParser.getNextSequence().build();
    SequenceDatabase reference = new SequenceDatabase(referenceContig);

    // logger
    Logger logger = new Logger(new PrintWriter());
    // compute alignments
 AlignmentCache resultsCache = new AlignmentCache();
    AlignerWorker worker = new AlignerWorker(reference, 0, listeners, resultsCache, new ArrayDeque<AlignerWorker>());
    worker.process(alignments, logger);

    // format vcf
    ByteArrayOutputStream vcfStream = new ByteArrayOutputStream();
    VcfWriter vcfWriter = new VcfWriter(vcfStream, true);
    try {
      vcfWriter.write(matchDatabase.groupByPosition(), 1);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }

    // return results
    return vcfStream.toString();
  }

  private List<List<SamAlignment_Builder>> getQueries(GroupedAlignment_Provider queryProvider) {
    List<List<SamAlignment_Builder>> queries = new ArrayList<List<SamAlignment_Builder>>();
    while (true) {
      List<SamAlignment_Builder> newQuery = queryProvider.getNextGroup();
      if (newQuery == null)
        return queries;
      queries.add(newQuery);
    }
  }

  private void fail(String message) {
    Assert.fail(message);
  }
}
