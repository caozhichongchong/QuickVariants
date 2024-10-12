package mapper;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

// An UnalignedQuery_Writer writes queries that did not align
public class UnalignedQuery_Writer implements AlignmentListener {
  public UnalignedQuery_Writer(String path, boolean allReadsContainQualityInformation) throws FileNotFoundException {
    if (path.endsWith(".fastq")) {
      if (!allReadsContainQualityInformation) {
        throw new IllegalArgumentException("Cannot write a .fastq file (" + path + ") when not all input reads are in .fastq format");
      }
      this.initialize(new FastqWriter(path));
    } else {
      if (path.endsWith(".fasta")) {
        this.initialize(new FastaWriter(path));
      } else {
        throw new IllegalArgumentException("Unsupported output type (must be .fastq or .fasta): " + path);
      }
    }
  }

  private void initialize(SequenceWriter writer) {
    this.writer = writer;
  }

  public void addAlignments(List<QueryAlignments> alignments) {
    synchronized(this) {
      for (QueryAlignments queryAlignments: alignments) {
        for (int i = 0; i < queryAlignments.getNumComponents(); i++) {
          List<QueryAlignment> subAlignments = queryAlignments.getAlignments(i);
          if (subAlignments.size() < 1) {
            // this sequence is unaligned
            if (queryAlignments.getNumQueries() == 1) {
              // all sequences are unaligned so write all of them
              for (Sequence sequence: queryAlignments.getSequences()) {
                this.writer.write(sequence);
              }
            } else {
              // only this sequence is unaligned so just write this one
              Sequence sequence = queryAlignments.getSequence(i);
              if (sequence != null) {
                this.writer.write(sequence);
              }
            }
          }
        }
      }
    }
  }

  public void close() {
    this.writer.close();
  }

  SequenceWriter writer;
}
