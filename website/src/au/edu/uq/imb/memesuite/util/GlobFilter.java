package au.edu.uq.imb.memesuite.util;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

import static java.nio.file.FileVisitResult.*;

/**
 */
public class GlobFilter {

  private static class GlobVisitor extends SimpleFileVisitor<Path> {
    private final Path start;
    private List<PathMatcher> matchers;
    private List<File> files;

    public GlobVisitor(Path start, String pattern) {
      this.start = start;
      matchers = new ArrayList<PathMatcher>();
      // split the pattern on unescaped spaces
      boolean space = true;
      boolean escape = false;
      int left = -1;
      for (int i = 0; i < pattern.length(); i++) {
        char c = pattern.charAt(i);
        if (space) {
          if (Character.isWhitespace(c)) continue;
          left = i;
          space = false;
        }
        if (escape) {
          escape = false;
        } else if (c == '\\') {
          escape = true;
        } else if (Character.isWhitespace(c)) {
          matchers.add(FileSystems.getDefault().getPathMatcher("glob:" + pattern.substring(left, i)));
          left = -1;
          space = true;
        }
      }
      if (left >= 0) matchers.add(FileSystems.getDefault().getPathMatcher("glob:" + pattern.substring(left)));
      files = new ArrayList<File>();
    }

    public File[] getFiles() {
      return files.toArray(new File[files.size()]);
    }

    // Invoke the pattern matching
    // method on each file.
    @Override
    public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
      Path relFile = start.relativize(file);
      for (PathMatcher matcher : matchers) {
        if (matcher.matches(relFile)) {
          files.add(file.toFile());
          break;
        }
      }
      return CONTINUE;
    }
  }

  public static File[] find(File dir, String pattern) throws IOException{
    Path start = dir.toPath();
    GlobVisitor globVisitor = new GlobVisitor(start, pattern);
    Files.walkFileTree(start, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, globVisitor);
    return globVisitor.getFiles();
  }

}
