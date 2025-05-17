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

public class MutationsWriter_Test {

  @Test
  public void testNoMutations() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t4M\t*\t*\t4\tACGT\t*";

    String ref  = "ACGTAAAAACGTAAAA";

    String mutations = buildMutations(sam1, ">contig1\n" + ref);

    String expectedMutations = "";

    checkMutations(mutations, expectedMutations);
  }

  @Test
  public void testOneMutation() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t5M\t*\t*\t5\tACGTT\t*";

    String ref  = "ACGTAAAAA";

    String mutations = buildMutations(sam1, ">contig1\n" + ref);

    String expectedMutations =
        "contig1	5	A	T	1	1\n" +
        "";

    checkMutations(mutations, expectedMutations);
  }

  @Test
  public void testConsecutiveMutations() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t9M\t*\t*\t5\tACGTTTAAA\t*";

    String ref  = "ACGTAAAAA";

    String mutations = buildMutations(sam1, ">contig1\n" + ref);

    String expectedMutations =
        "contig1	5	A	T	1	1\n" +
        "contig1	6	A	T	1	1\n" +
        "";

    checkMutations(mutations, expectedMutations);
  }

  @Test
  public void testInsertion() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t3M2I3M\t*\t*\t12\tACGGACTT\t*";

    String ref  = "ACGCTT";

    String mutations = buildMutations(sam1, ">contig1\n" + ref);

    String expectedMutations =
        "contig1	3	--	GA	1	1\n" +
        "";

    checkMutations(mutations, expectedMutations);
  }

  @Test
  public void testDeletion() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t6M2D6M\t*\t*\t12\tACGTAACCGGTT\t*";

    String ref  = "ACGTAAAACCGGTT";

    String mutations = buildMutations(sam1, ">contig1\n" + ref);

    String expectedMutations =
        "contig1	7	AA	--	1	1\n" +
        "";

    checkMutations(mutations, expectedMutations);
  }

  @Test
  public void testIgnoringMutationWithLowDepth() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t8M\t*\t*\t12\tACGTAACT";

    String ref  = "ACGTACGT";

    MutationDetectionParameters filter = new MutationDetectionParameters();
    filter.minSNPTotalDepth = 2;

    String filteredMutations = buildMutations(sam1, ">contig1\n" + ref, filter, 0);
    String expectedFilteredMutations =
        "";
    checkMutations(filteredMutations, expectedFilteredMutations);

    String unfilteredMutations = buildMutations(sam1, ">contig1\n" + ref);
    String expectedUnfilteredMutations =
        "contig1	6	C	A	1	1\n" +
        "contig1	7	G	C	1	1\n" +
        "";
    checkMutations(unfilteredMutations, expectedUnfilteredMutations);
  }

  @Test
  public void testIgnoringIndelNearQueryEnd() {
    String sam1 = "name1\t0\tcontig1\t1\t255\t1M1D6M\t*\t*\t12\tAGTAACT\t*";

    String ref  = "ACGTAACT";

    MutationDetectionParameters filter = new MutationDetectionParameters();
    filter.minIndelTotalStartDepth = 1;
    double distinguishQueryEnds = 0.5;

    String filteredMutations = buildMutations(sam1, ">contig1\n" + ref, filter, distinguishQueryEnds);
    String expectedFilteredMutations =
        "";
    checkMutations(filteredMutations, expectedFilteredMutations);

    String unfilteredMutations = buildMutations(sam1, ">contig1\n" + ref);
    String expectedUnfilteredMutations =
        "contig1	2	C	-	1	1\n" +
        "";
    checkMutations(unfilteredMutations, expectedUnfilteredMutations);
  }

  private void checkMutations(String actual, String expected) {
    actual = withoutMetadataLines(actual);
    if (!(expected.equals(actual))) {
      String actualAsCode = "\"" + actual.replace("\n", "\\n\" +\n        \"") + "\";";
      fail("Difference in generated mutations file.\nactual:\n" + actual + "\nexpected:\n" + expected + "\n code for actual:\n" + actualAsCode);
    }
  }

  private String withoutMetadataLines(String original) {
    String[] lines = original.split("\n");
    StringBuilder resultBuilder = new StringBuilder();
    for (String line: lines) {
      if (!line.startsWith("#") && !line.startsWith("CHR") && line.length() > 0) {
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

  private String buildMutations(String samRecords, String referenceGenome) {
    MutationDetectionParameters parameters = MutationDetectionParameters.emptyFilter();
    return buildMutations(samRecords, referenceGenome, parameters, 0);
  }

  private String buildMutations(String samRecords, String referenceGenome, MutationDetectionParameters parameters, double distinguishQueryEnds) {
    // make alignment listener
    MatchDatabase matchDatabase = new MatchDatabase(distinguishQueryEnds);
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

    // format mutations
    ByteArrayOutputStream mutationsStream = new ByteArrayOutputStream();
    MutationsWriter mutationsWriter = new MutationsWriter(mutationsStream, parameters);
    try {
      mutationsWriter.write(matchDatabase.groupByPosition(), 1);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }

    // return results
    return mutationsStream.toString();
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
