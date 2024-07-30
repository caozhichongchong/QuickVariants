package mapper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

// A DataLoader reads files (such as .fasta files) and returns the corresponding representation of them
public class DataLoader {
  public static SequenceProvider LoadFrom(List<String> paths, boolean keepQualityData) throws IllegalArgumentException, IOException, FileNotFoundException {
    List<SequenceProvider> providers = new ArrayList<SequenceProvider>();
    for (String path : paths) {
      providers.add(readSequencesFrom(path, keepQualityData, false, null));
    }
    return new SequencesIterator(providers);
  }

  public static SequenceProvider LoadFrom(String path, boolean keepQualityData) throws IllegalArgumentException, IOException, FileNotFoundException {
    return readSequencesFrom(path, keepQualityData, false, null);
  }

  public static SequenceProvider LoadFrom(String path, boolean keepQualityData, boolean groupLinesInSamFiles) throws IllegalArgumentException, IOException, FileNotFoundException {
    return readSequencesFrom(path, keepQualityData, groupLinesInSamFiles, null);
  }

  public static SequenceProvider LoadSam(String path, boolean keepQualityData, boolean groupLinesInSamFiles) throws IllegalArgumentException, IOException, FileNotFoundException {
    return readSequencesFrom(path, keepQualityData, groupLinesInSamFiles, ".sam");
  }

  public static GroupedQuery_Provider ParseSamAlignments(String path, boolean keepQualityData, boolean groupLinesInSamFiles) throws IllegalArgumentException, IOException, FileNotFoundException {
    SequenceProvider sequenceProvider = readSequencesFrom(path, keepQualityData, groupLinesInSamFiles, ".sam");
    QueryProvider queryBuilder = new SimpleQueryProvider(sequenceProvider);
    GroupedQuery_Provider groupProvider = new GroupedQuery_Provider(queryBuilder);
    return groupProvider;
  }

  private static SequenceProvider readSequencesFrom(String path, boolean keepQualityData, boolean groupLinesInSamFiles, String overrideFileType) throws IllegalArgumentException, IOException, FileNotFoundException {
    String effectivePath = path;

    InputStream inputStream;
    if ("-".equals(path)) {
      inputStream = System.in;
      path = "stdin";
      effectivePath = path;
    } else {
      inputStream = new FileInputStream(path);
      String gzipSuffix = ".gz";
      if (path.endsWith(gzipSuffix)) {
        inputStream = new GZIPInputStream(inputStream);
        effectivePath = effectivePath.substring(0, effectivePath.length() - gzipSuffix.length());
      }
    }

    BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

    String fileType = effectivePath;
    if (overrideFileType != null) {
      fileType = overrideFileType;
    }

    if (fileType.endsWith(".fasta") || fileType.endsWith(".fa") || fileType.endsWith(".fna")) {
      return new FastaParser(reader, path);
    }
    if (fileType.endsWith(".fastq") || fileType.endsWith(".fq") || fileType.endsWith(".ca")) {
      return new FastqParser(reader, path, keepQualityData);
    }
    if (fileType.endsWith(".sam")) {
      SamProvider samProvider = new SamReader(reader, path);
      if (groupLinesInSamFiles)
        samProvider = new SamGrouper(samProvider);
      return new SamParser(samProvider, path);
    }

    throw new IllegalArgumentException("Unrecognized file extension: " + effectivePath + ", not .sam or .fasta/.fa or .fastq/.fq/.ca");
  }
}
