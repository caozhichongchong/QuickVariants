package mapper;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

// A SamWriter writes .sam files
// See https://samtools.github.io/hts-specs/SAMv1.pdf for more information
public class SamWriter implements AlignmentListener {
  public SamWriter(SequenceDatabase sequenceDatabase, String path) throws FileNotFoundException {
    File file = new File(path);
    this.fileStream = new FileOutputStream(file);
    this.bufferedStream = new BufferedOutputStream(fileStream);
    // write header
    this.write("@HD\tVN:1.6\tSO:unsorted\n");
    // write reference sequence names
    this.writeReferenceSequenceNames(sequenceDatabase);
  }

  private void writeReferenceSequenceNames(SequenceDatabase sequenceDatabase) {
    int count = sequenceDatabase.getNumSequences();
    for (int i = 0; i < count; i++) {
      Sequence contig = sequenceDatabase.getSequence(i);
      if (contig.getComplementedFrom() == null)
        this.write("@SQ\tSN:" + contig.getName() + "\tLN:" + contig.getLength() + "\n");
    }
  }

  public void addAlignments(List<List<QueryAlignment>> alignments) {
    this.write(this.format(alignments));
  }

  public void addUnaligned(List<SamAlignment> unalignedQueries) {
  }

  public void close() {
    try {
      this.bufferedStream.close();
    } catch (IOException e) {
    }
    try {
      this.fileStream.close();
    } catch (IOException e) {
    }
  }

  private String format(List<List<QueryAlignment>> alignments) {
    StringBuilder builder = new StringBuilder();
    for (List<QueryAlignment> queryAlignments: alignments) {
      for (QueryAlignment queryAlignment: queryAlignments) {
         for (SequenceAlignment alignment: queryAlignment.getComponents()) {
           Sequence query = alignment.getSequenceA();
           Sequence ref = alignment.getSequenceB();
           // QNAME
           builder.append(query.getSourceName());
           builder.append('\t');
           // FLAG
           builder.append("" + getSamFlags(alignment));
           builder.append('\t');
           // RNAME
           builder.append(ref.getName());
           builder.append('\t');
           // POS
           builder.append(getSamReferencePosition(alignment));
           builder.append('\t');
           // MAPQ
           builder.append("255\t");
           // CIGAR flags
           int queryLengthConsumed = 0;
           for (AlignedBlock block : alignment.getSections()) {
             if (block.getStartIndexA() != queryLengthConsumed) {
               // left side of the query fell off the reference
               builder.append("" + block.getStartIndexA() + "S");
               queryLengthConsumed = block.getStartIndexA();
             }
             if (block.getLengthA() == block.getLengthB()) {
               builder.append("" + block.getLengthA() + "M");
             } else {
               if (block.getLengthA() > block.getLengthB()) {
                 builder.append("" + block.getLengthA() + "I");
               } else {
                 builder.append("" + block.getLengthB() + "D");
               }
             }
             queryLengthConsumed += block.getLengthA();
           }
           if (queryLengthConsumed < query.getLength()) {
             // right side of the query fell off the reference
             builder.append("" + (query.getLength() - queryLengthConsumed) + "S");
           }
           builder.append('\t');
           SequenceAlignment pairedSequenceAlignment = getPaired(queryAlignment, alignment);
           if (pairedSequenceAlignment != null) {
             // RNEXT:
             builder.append(pairedSequenceAlignment.getSequenceB().getName());
             builder.append('\t');
             // PNEXT:
             builder.append(getSamReferencePosition(pairedSequenceAlignment));
             builder.append('\t');
           } else {
             // RNEXT:
             builder.append("*\t");
             // PNEXT:
             builder.append("0\t");
           }
           // TLEN
           builder.append("" + query.getLength());
           builder.append('\t');
           // SEQ
           builder.append(query.getText());
           builder.append('\t');
           // QUAL
           builder.append("*\t");
           // alignment score
           builder.append(formatPenalty(alignment));
           builder.append("\n");
        }
      }
    }
    return builder.toString();
  }

  private String formatPenalty(SequenceAlignment alignment) {
    float score = 0; // not keeping track
    float roundedScore = Math.round(score);
    if (Math.abs(score - roundedScore) < Math.abs(score) * 0.000001) {
      // score is essentially an integer
      return "AS:i:" + (int)roundedScore;
    } else {
      // score is a float
      return "AS:f:" + score;
    }
  }

  private int getSamFlags(SequenceAlignment alignment) {
    int flags = 0;
    if (alignment.isReferenceReversed()) {
      flags += 16;
    }
    return flags;
  }

  // gets the other/paired SequenceAlignment in a QueryAlignment, or null if there is no other
  private SequenceAlignment getPaired(QueryAlignment queryAlignment, SequenceAlignment sequenceAlignment) {
    if (queryAlignment.getNumSequences() != 2)
      return null;
    List<SequenceAlignment> sequenceAlignments = queryAlignment.getComponents();
    if (sequenceAlignments.get(0) == sequenceAlignment)
      return sequenceAlignments.get(1);
    return sequenceAlignments.get(0);
  }

  private int getSamReferencePosition(SequenceAlignment sequenceAlignment) {
    return sequenceAlignment.getSection(0).getStartIndexB() + 1;
  }

  private void write(String text) {
    synchronized(this.bufferedStream) {
      try {
        this.bufferedStream.write(text.getBytes());
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }
  }

  FileOutputStream fileStream;
  BufferedOutputStream bufferedStream;
}
