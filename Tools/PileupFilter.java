import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PileupFilter {
  
  private static final String BASE_TABLE = "ACGTURYKMSWBDHVN><*"; // TODO why so many?
  private static final int[] baseIds = new int[256];
  private static int baseCount = 0;
  private static final Map<String, String> parameters = new HashMap<String, String>();
  
  private static String inputFile = null;
  private static String outputFile = null;
  private static int minBaseQuality = 0;
  private static int minMapQuality = 0;
  private static int asciiBase = 33;
  private static boolean filtered = false;
  
  public static void main(String[] args) throws Exception {
    if (!init(args)) {
      return;
    }

    InputStream in = null;
    if (inputFile != null && !inputFile.isEmpty()) {
      in = new FileInputStream(inputFile);
    } else {
      in = System.in;
    }

    OutputStream out = null;
    if (outputFile != null && !outputFile.isEmpty()) {
      out = new FileOutputStream(outputFile);
    } else {
      out = System.out;
    }
    
    BufferedReader r = new BufferedReader(new InputStreamReader(in));
    BufferedWriter w = new BufferedWriter(new OutputStreamWriter(out));
    
    try {
      for (;;) {
        String l = r.readLine();
        if (l == null) {
          break;
        }
        
        String ret = null;
        try {
          ret = process(l);
        } catch (Exception e) {
          throw new Exception(l, e);
        }
        
        if (ret != null) {
          w.write(ret);
        }
      }
    } finally {
      r.close();
      w.close();
    }
  }
  
  private static boolean init(String[] args) {
    // build base id mapping, from char('ACGT...') to number
    for (int i = 0; i < 256; ++i) {
      baseIds[i] = -1;
    }
    baseCount = 0;
    for (char c : BASE_TABLE.toCharArray()) {
      baseIds[c] = baseCount;
      baseCount++;
      char lc = lc(c);
      if (lc != c) {
        baseIds[lc] = baseCount;
        baseCount++;
      }
    }

    // parse parameters
    Pattern argPattern = Pattern.compile("^--([a-zA-Z0-9_]+)(=(.*))?$");
    for (String arg : args) {
      Matcher matcher = argPattern.matcher(arg.trim());
      if (matcher.find()) {
        parameters.put(matcher.group(1).toLowerCase(), matcher.group(3));
      }
    }
    
    if (parameters.isEmpty() || parameters.containsKey("help")) {
      help();
      return false;
    }
    inputFile = parameters.get("inputfile");
    outputFile = parameters.get("outputfile");
    if (parameters.containsKey("minbasequal")) {
      minBaseQuality = Integer.parseInt(parameters.get("minbasequal"));
    }
    if (parameters.containsKey("minmapqual")) {
      minMapQuality = Integer.parseInt(parameters.get("minmapqual"));
    }
    if (parameters.containsKey("asciibase")) {
      asciiBase = Integer.parseInt(parameters.get("asciibase"));
    }
    if (parameters.containsKey("filtered")) {
      filtered = "1".equals(parameters.get("filtered")) || 
                 "true".equals(parameters.get("filtered"));
    }
    return true;
  }
  
  private static void help() {
    System.out.println("\n" +
        "Usage: javac PileupFilter.java\n" + 
        "       java PileupFilter \n" +
        "         --inputfile    input file name or stdin if mssing\n" +
        "         --outputfile    output file name or stdout if mssing\n" +
        "         --minbasequal    Minimum base quality [0]\n" +
        "         --minmapqual     Minimum mapping quality [0]\n" +
        "         --asciibase      For Illumina 1.8+ or Sanger, set to 33; for Illumina 1.3+, set to 64 [33]\n" +
        "         --filtered       Whether to output filtered base quality for each base and strand, instead of original base quality [0]\n" + 
        "         --help           Display help info\n\n" +
        "Format of output file:\n" +
        "       \\$1     chr\n" +
        "       \\$2     position\n" +
        "       \\$3     ref nt\n" +
        "       \\$4     total depth\n" +
        "       \\$5     stack of nt\n" +
        "       \\$6     stack of base quality\n" +
        "       \\$7     stack of mapping quality\n" +
        "       \\$8     depth of ref nt in + strand\n" +
        "       \\$9     depth of ref nt in - strand\n" +
        "       \\$10    depth of alt1 nt in + strand\n" +
        "       \\$11    depth of alt1 nt in - strand\n" +
        "       \\$12    alt1 nt\n" +
        "       \\$13    depth of alt2 nt in + strand\n" +
        "       \\$14    depth of alt2 nt in - strand\n" +
        "       \\$15    alt2 nt\n" +
        "       \\$16    stack of base quality for ref nt in + strand\n" +
        "       \\$17    stack of base quality for ref nt in - strand\n" +
        "       \\$18    stack of base quality for alt1 nt in + strand\n" +
        "       \\$19    stack of base quality for alt1 nt in - strand\n" +
        "       \\$20    stack of base quality for alt2 nt in + strand\n" +
        "       \\$21    stack of base quality for alt2 nt in - strand\n\n");
  }
  
  private static String process(String input) {
    String[] tokens = split(input, '\t', 7);
    
    String chr = tokens[0];
    String coor = tokens[1];
    char ref = tokens[2].charAt(0);
    int num = Integer.parseInt(tokens[3]);
    String readList = tokens[4];
    String baseQualityList = tokens[5];
    String mapQualityList = tokens[6];
    
    int filteredNum = 0;
    StringBuilder filteredRead = new StringBuilder();
    StringBuilder filteredBaseQuality = new StringBuilder();
    StringBuilder filteredMapQuality = new StringBuilder();
    
    Read[] reads = parseReadList(readList, num);
    int[] alle = new int[baseCount];
    StringBuilder[] baseQuality = new StringBuilder[baseCount];
    
    // filter based on quality
    for (int i = 0; i < num; ++i) {
      char ch = reads[i].base;
      if (ch == '.') {
        ch = ref;
      } else if (ch == ',') {
        ch = lc(ref);
      }
      int id = baseIds[ch];
      if (baseQuality[id] == null) {
        baseQuality[id] = new StringBuilder();
      }
      if (isGoodQuality(baseQualityList.charAt(i), minBaseQuality) &&
          isGoodQuality(mapQualityList.charAt(i), minMapQuality)) {
        filteredNum++;
        filteredRead.append(reads[i].value);
        filteredBaseQuality.append(baseQualityList.charAt(i));
        filteredMapQuality.append(mapQualityList.charAt(i));
        alle[id]++;
        baseQuality[id].append(baseQualityList.charAt(i));
      } else if (!filtered) {
        baseQuality[id].append(baseQualityList.charAt(i));
      }
    }
    if (filteredNum == 0) {
      return null;
    }
    
    // sort
    char[] chs = "TCGA".toCharArray();
    char tmp;
    for (int i = 1; i < 4; ++i) {
      char uc1 = chs[i];
      char lc1 = lc(uc1);
      int n1 = alle[baseIds[uc1]] + alle[baseIds[lc1]];
      int k = -1;
      for (int j = 0; j < i; ++j) {
        char uc2 = chs[j];
        char lc2 = lc(uc2);
        int n2 = alle[baseIds[uc2]] + alle[baseIds[lc2]];
        if (uc2 == ref || (uc1 != ref && n1 > n2)) {
          k = j;
          break;
        }
      }
      if (k >= 0) {
        tmp = chs[i];
        for (int j = k; j < i; ++j) {
          chs[j + 1] = chs[j];
        }
        chs[k] = tmp;
      }
    }
    
    // build results
    StringBuilder ret = new StringBuilder();
    ret.append(chr).append('\t').
        append(coor).append('\t').
        append(ref).append('\t').
        append(filteredNum).append('\t').
        append(filteredRead).append('\t').
        append(filteredBaseQuality).append('\t').
        append(filteredMapQuality).append('\t');
    
    ret.append(alle[baseIds[ref]]).append('\t').
        append(alle[baseIds[lc(ref)]]).append('\t');
    
    int[] count = new int[2];
    for (int i = 0; i < 2; ++i) {
      int countUc = alle[baseIds[chs[i]]];
      int countLc = alle[baseIds[lc(chs[i])]];
      count[i] = countUc + countLc;
      char ch =count[i] == 0 ? '.' : chs[i];
      ret.append(countUc).append('\t').
          append(countLc).append('\t').
          append(ch).append('\t');
    }
    
    appendQuality(ret, baseQuality[baseIds[ref]], '\t', false);
    appendQuality(ret, baseQuality[baseIds[lc(ref)]], '\t', false);
    appendQuality(ret, baseQuality[baseIds[chs[0]]], '\t', count[0] == 0);
    appendQuality(ret, baseQuality[baseIds[lc(chs[0])]], '\t', count[0] == 0);
    appendQuality(ret, baseQuality[baseIds[chs[1]]], '\t', count[1] == 0);
    appendQuality(ret, baseQuality[baseIds[lc(chs[1])]], '\n', count[1] == 0);
    
    return ret.toString();
  }
  
  private static void appendQuality(
      StringBuilder sb, StringBuilder quality, char delimiter, boolean ignore) {
    if (quality == null || quality.length() == 0 || ignore) {
      sb.append('|').append(delimiter);
    } else {
      sb.append(quality).append(delimiter);
    }
  }
  
  private static char lc(char ch) {
    return Character.toLowerCase(ch);
  }
  
  private static boolean isGoodQuality(char ch, int minQuality) {
    return ch - asciiBase >= minQuality;
  }
  
  private static Read[] parseReadList(String s, int n) {
    Read[] reads = new Read[n];
    char[] chars = s.toCharArray();
    int p0 = 0;
    for (int i = 0; i < n; ++i) {
      char base = 0;
      boolean startMark = false;
      char mappingQuality = 0;
      boolean endMark = false;
      String insertion = null;
      String deletion = null;

      int p1 = p0;
      if (chars[p1] == '^') {
        startMark = true;
        mappingQuality = chars[p1 + 1];
        p1 += 2;
      }
      base = chars[p1];
      p1++;
      if (p1 < chars.length && (chars[p1] == '+' || chars[p1] == '-')) {
        char op = chars[p1];
        p1++;
        int m = 0;
        while (chars[p1] >= '0' && chars[p1] <= '9') {
          m = m * 10 + (chars[p1] - '0');
          p1++;
        }
        if (op == '+') {
          insertion = new String(chars, p1, m);
        } else {
          deletion = new String(chars, p1, m);;
        }
        p1 += m;
      }
      if (p1 < chars.length && chars[p1] == '$') {
        endMark = true;
        p1++;
      }
      reads[i] = new Read(base, startMark, endMark, mappingQuality, insertion, deletion);
      p0 = p1;
    }
    return reads;
  }
  
  private static String[] split(String s, char delimiter, int n) {
    String[] tokens = new String[n];
    int i = 0;
    int p0 = 0; 
    char[] chars = s.toCharArray();
    for (int p1 = 0; p1 <= chars.length; ++p1) {
      if (p1 == chars.length || chars[p1] == delimiter) {
        tokens[i] = new String(chars, p0, p1 - p0);
        p0 = p1 + 1;
        i++;
        if (i >= n) {
          break;
        }
      }
    }
    return tokens;
  }
}

class Read {
  public final char base;
  public final boolean startMark;
  public final char mappingQuality;
  public final boolean endMark;
  public final String insertion;
  public final String deletion;
  public final String value;
  
  public Read(char base) {
    this(base, false, false, (char) 0, null, null);
  }
  
  public Read(char base, boolean startMark, boolean endMark, 
              char mappingQuality, String insertion, String deletion) {
    this.base = base;
    this.startMark = startMark;
    this.endMark = endMark;
    this.mappingQuality = mappingQuality;
    this.insertion = insertion;
    this.deletion = deletion;
    
    StringBuilder sb = new StringBuilder();
    if (startMark) {
      sb.append('^').append(mappingQuality);
    }
    sb.append(base);
    if (insertion != null) {
      sb.append('+').append(insertion.length()).append(insertion);
    }
    if (deletion != null) {
      sb.append('-').append(deletion.length()).append(deletion);
    }
    if (endMark) {
      sb.append('$');
    }
    value = sb.toString();
  }
  
  @Override
  public String toString() {
    return value;
  }
}
