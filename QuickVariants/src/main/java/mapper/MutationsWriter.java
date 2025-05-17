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

// A MutationsWriter writes mutations observed in the queries compared to the reference
public class MutationsWriter {
  public MutationsWriter(String path, MutationDetectionParameters parameters) throws FileNotFoundException, IOException {
    this.initialize(new FileOutputStream(new File(path)));
    this.parameters = parameters;
  }

  public MutationsWriter(OutputStream destination, MutationDetectionParameters parameters) {
    this.initialize(destination);
    this.parameters = parameters;
  }

  private void initialize(OutputStream destination) {
    this.fileStream = destination;
    this.bufferedStream = new BufferedOutputStream(destination);
  }

  public void write(Map<Sequence, Alignments> alignments, int numParallelJobs) throws IOException {
    this.writeText("# Version 1\n");
    this.writeText("CHROM\tPOS\tREF\tALT\tDEPTH\tTOTAL_DEPTH\n");

    List<MutationsFormatRequest> jobs = this.splitJobs(alignments, numParallelJobs);
    int waitIndex = 0; // index of next worker to wait for
    List<MutationsFormatterWorker> workers = new ArrayList<MutationsFormatterWorker>();

    while (waitIndex < jobs.size()) {
      boolean hasCapacityToLaunchJob = (workers.size() < waitIndex + numParallelJobs);
      boolean hasJobToLaunch = workers.size() < jobs.size();
      if (hasCapacityToLaunchJob && hasJobToLaunch) {
        workers.add(requestFormat(jobs.get(workers.size())));
      } else {
        // cannot launch a new job; wait for one to complete instead
        MutationsFormatterWorker worker = workers.get(waitIndex);
        try {
          worker.join();
        } catch (InterruptedException e) {
        }
        // write results
        String result = worker.getResults();
        this.writeText(result);
        // clear worker
        workers.set(waitIndex, null);
        waitIndex++;
      }
    }

    bufferedStream.close();
    fileStream.close();
  }

  private List<MutationsFormatRequest> splitJobs(Map<Sequence, Alignments> alignments, int numParallelJobs) {
    List<MutationsFormatRequest> jobs = new ArrayList<MutationsFormatRequest>();

    TreeMap<String, Sequence> sortedSequences = new TreeMap<String, Sequence>();
    for (Map.Entry<Sequence, Alignments> entry : alignments.entrySet()) {
      sortedSequences.put(entry.getKey().getName(), entry.getKey());
    }
    int maxJobSize = 8192;

    for (Map.Entry<String, Sequence> entry : sortedSequences.entrySet()) {
      Sequence sequence = entry.getValue();
      Alignments alignmentsHere = alignments.get(sequence);
      FilteredAlignments filteredAlignments = new FilteredAlignments(alignmentsHere, this.parameters);
      int startIndex = 0;
      while (startIndex < sequence.getLength()) {
        int jobSize = Math.min(maxJobSize, Math.max(1, (jobs.size() + 1) / numParallelJobs * maxJobSize));

        int endIndex = Math.min(sequence.getLength(), startIndex + jobSize);
        int length = endIndex - startIndex;
        int jobId = jobs.size();
        jobs.add(new MutationsFormatRequest(sequence, startIndex, length, filteredAlignments, jobId));
        startIndex = endIndex;
      }
    }
    return jobs;
  }

  private MutationsFormatterWorker requestFormat(MutationsFormatRequest formatRequest) {
    MutationsFormatterWorker worker = new MutationsFormatterWorker();
    worker.request(formatRequest);
    worker.start();
    return worker;
  }

  private void writeText(String text) throws IOException {
    bufferedStream.write(text.getBytes());
  }
  BufferedOutputStream bufferedStream;
  OutputStream fileStream;
  MutationDetectionParameters parameters;
}
