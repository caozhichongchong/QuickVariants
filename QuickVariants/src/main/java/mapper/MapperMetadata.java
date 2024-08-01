package mapper;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.util.List;
import java.util.Properties;
import java.lang.management.ManagementFactory;

// Returns metadata about QuickVariants
public class MapperMetadata {

  // The version of QuickVariants that is running
  public static String getVersion() {
    Properties properties = new Properties();
    try {
      properties.load(MapperMetadata.class.getResourceAsStream("/mapper.properties"));
    } catch (IOException e) {
      throw new RuntimeException("Failed to get QuickVariants version", e);
    }
    String version = properties.getProperty("mapper.version", "unknown");
    return version;
  }

  // A guess of the command line used to run this program
  // This might not get the correct filepath for java
  // This might not get exactly the correct filepath of the QuickVariants jar
  public static String guessCommandLine() {
    String javaArgumentsString = String.join(" ", getJavaArguments());
    Path quickVariantsJarPath = getQuickVariantsPath();
    Path workingDirPath = new File(".").toPath();

    Path simplifiedQuickVariantsPath = simplifyPath(quickVariantsJarPath, workingDirPath);

    String[] mainArguments = getMainArguments();
    if (mainArguments == null) {
      return "unknown"; // probably running in a test
    }
    
    String mainArgumentsString = String.join(" ", getMainArguments());
    return "java " + javaArgumentsString + " -jar " + simplifiedQuickVariantsPath + " " + mainArgumentsString;
  }

  private static Path getQuickVariantsPath() {
    try {
      File quickVariantsJar = new File(MapperMetadata.class.getProtectionDomain().getCodeSource().getLocation().toURI());
      return quickVariantsJar.toPath();
    } catch (URISyntaxException e) {
      throw new RuntimeException(e);
    }
  }

  // Returns a simpler path for destPath given a working directory of workingDir
  public static Path simplifyPath(Path destPath, Path workingDir) {
    Path relativePath = workingDir.toAbsolutePath().relativize(destPath.toAbsolutePath());
    // Return relative path if it's shorter
    if (relativePath.toString().length() <= destPath.toString().length()) {
      return relativePath;
    }
    return destPath;
  }

  public static List<String> getJavaArguments() {
    return ManagementFactory.getRuntimeMXBean().getInputArguments();
  }

  public static String[] getMainArguments() {
    return mainArguments;
  }
  public static void setMainArguments(String[] arguments) {
    mainArguments = arguments;
  }

  private static String[] mainArguments;

}
