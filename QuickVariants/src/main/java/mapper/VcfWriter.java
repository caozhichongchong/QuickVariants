package mapper;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Map;
import java.util.List;
import java.util.TreeMap;

// A VcfWriter writes .vcf files
public class VcfWriter {
  public VcfWriter(String path, boolean includeNonMutations, MutationDetectionParameters mutationsFilter) throws FileNotFoundException, IOException {
    this.initialize(new FileOutputStream(new File(path)));
    this.includeNonMutations = includeNonMutations;
    this.mutationsFilter = mutationsFilter;
  }

  public VcfWriter(OutputStream destination, boolean includeNonMutations) {
    this.initialize(destination);
    this.includeNonMutations = includeNonMutations;
    this.mutationsFilter = MutationDetectionParameters.emptyFilter();
  }

  private void initialize(OutputStream destination) {
    this.fileStream = destination;
    this.bufferedStream = new BufferedOutputStream(destination);
  }

  public void write(Map<Sequence, Alignments> alignments, int numParallelJobs) throws IOException {
    this.writeText("##fileType=\"Vcf summary of variants\"\n");
    this.writeText("\n");
    this.writeText("##commandLine=\"" + MapperMetadata.guessCommandLine() + "\"\n");
    this.writeText("\n");
    this.writeText("#CHROM\tPOS\tREF\tALT\tDP\tDETAILS-MIDDLE\tDETAILS-ENDS\tSUPPORT\n");

    List<VcfFormatRequest> jobs = this.splitJobs(alignments, numParallelJobs);
    int waitIndex = 0; // index of next worker to wait for
    List<VcfFormatterWorker> workers = new ArrayList<VcfFormatterWorker>();

    while (waitIndex < jobs.size()) {
      boolean hasCapacityToLaunchJob = (workers.size() < waitIndex + numParallelJobs);
      boolean hasJobToLaunch = workers.size() < jobs.size();
      if (hasCapacityToLaunchJob && hasJobToLaunch) {
        workers.add(requestFormat(jobs.get(workers.size())));
      } else {
        // cannot launch a new job; wait for one to complete instead
        VcfFormatterWorker worker = workers.get(waitIndex);
        try {
          worker.join();
        } catch (InterruptedException e) {
        }
        // write results
        String result = worker.getResults();
        this.writeText(result);
        this.numReferencePositionsMatched += worker.getNumReferencePositionsMatched();
        // clear worker
        workers.set(waitIndex, null);
        waitIndex++;
      }
    }

    bufferedStream.close();
    fileStream.close();
  }

  private List<VcfFormatRequest> splitJobs(Map<Sequence, Alignments> alignments, int numParallelJobs) {
    List<VcfFormatRequest> jobs = new ArrayList<VcfFormatRequest>();

    TreeMap<String, Sequence> sortedSequences = new TreeMap<String, Sequence>();
    for (Map.Entry<Sequence, Alignments> entry : alignments.entrySet()) {
      sortedSequences.put(entry.getKey().getName(), entry.getKey());
    }
    int maxJobSize = 8192;

    for (Map.Entry<String, Sequence> entry : sortedSequences.entrySet()) {
      Sequence sequence = entry.getValue();
      Alignments alignmentsHere = alignments.get(sequence);
      FilteredAlignments filteredAlignments = new FilteredAlignments(alignmentsHere, this.mutationsFilter);
      int startIndex = 0;
      while (startIndex < sequence.getLength()) {
        int jobSize = Math.min(maxJobSize, Math.max(1, (jobs.size() + 1) / numParallelJobs * maxJobSize));

        int endIndex = Math.min(sequence.getLength(), startIndex + jobSize);
        int length = endIndex - startIndex;
        int jobId = jobs.size();
        jobs.add(new VcfFormatRequest(sequence, startIndex, length, filteredAlignments, jobId));
        startIndex = endIndex;
      }
      this.numReferencePositions += sequence.getLength();
    }
    return jobs;
  }

  private VcfFormatterWorker requestFormat(VcfFormatRequest formatRequest) {
    VcfFormatterWorker worker = new VcfFormatterWorker(this.includeNonMutations);
    worker.request(formatRequest);
    worker.start();
    return worker;
  }

  public long getNumReferencePositions() {
    return numReferencePositions;
  }
  public long getNumReferencePositionsMatched() {
    return numReferencePositionsMatched;
  }

  private void writeText(String text) throws IOException {
    bufferedStream.write(text.getBytes());
  }
  BufferedOutputStream bufferedStream;
  OutputStream fileStream;

  long numReferencePositionsMatched;
  long numReferencePositions;
  boolean includeNonMutations;
  MutationDetectionParameters mutationsFilter;
}
