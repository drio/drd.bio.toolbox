import java.util.*;
import java.io.*;

public class SnpSites {
  private Hashtable<String, Vector<Integer>> sites;
  private FileInputStream in;
  private BufferedReader br;
  private final String NEW_INFO = "##INFO=<ID=RDP,Number=.,Type=String,Description=\"Number of reads in site for sample\">";
  private String vcfFn;
  private Vector<String> bamPaths = new Vector<String>();
  private int chunk_size;
  private int n_threads;
 
  public SnpSites(String fn, int nt, int cs) {
    sites = new Hashtable<String, Vector<Integer>>(); 
    vcfFn = fn;
    n_threads = nt;
    chunk_size = cs;
  }
 
  public void doWork() {
    openVCF();
    processHeader();
    computeSnps();
    System.err.println("# of bams found in header: " + bamPaths.size() );
    closeVCF();
  }
 
  // Iterate over the SNPs and add the new RDP INFO field
  private void addRDP(Vector <String> lines) {
    int colNumberInfo = 7;
    String chrm, coor;
   
    try {
      for (String line : lines) {
        String[] s = line.split("\t");
        chrm       = s[0];
        coor       = s[1];
        for (int i=0; i<s.length; i++) { // i: index in column
          if (i == colNumberInfo) { // append the RDP and print
            StringBuilder rdps = new StringBuilder(";RDP=");
            for (int j=0; j<bamPaths.size(); j++) {
              Vector <Integer> vCov = sites.get(chrm + "_" + coor); // Get coverage at site
              rdps.append(Integer.toString(vCov.get(j)) + ",");
            } 
            rdps.setCharAt(rdps.length()-1, '\t');
            System.out.print(s[i] + rdps);
          } else { 
            if (i == s.length - 1)
              System.out.print(s[i] + "\n");
            else
              System.out.print(s[i] + "\t");
          }
        }
      }
    } catch(Exception e) {
      bailOut(e);
    }
  }
 
  private void computeChunk(String refName, int start, int end, Vector<String> lines) {
    // We have a chunk of snp locations loaded. Let's find the raw coverage for all the bams
    BamCoverage bc = new BamCoverage(bamPaths, this, n_threads);
    bc.findCoverage(refName, start, end);
    addRDP(lines);
  }

  // Iterate over the vcf snps and save the coordinate for a bunch of them. When we have
  // a bunch we can compute the raw coverage and dump the new vcf lines with the new RDP.
  public void computeSnps() {
    try {
      Vector<String> chunkLines = new Vector<String>();
      int nLines = 0, start = 0, end = 0;
      String line = "", refName = "";
      String[] s = new String[]{};
      while ((line = br.readLine()) != null && line.length() != 0) {
        chunkLines.add(line);
        nLines++;
        s = line.split("\t");
        add(s[0], s[1]);
        
        if (nLines == 1)
          start = Integer.parseInt(s[1]);
        
        if (nLines == chunk_size-1)
          refName = s[0]; end = Integer.parseInt(s[1]);
        
        if (nLines == chunk_size) { // We have enough snps, compute coverage
          computeChunk(refName, start, end, chunkLines);
          // Get ready for a new chunk
          nLines = 0;
          sites = new Hashtable<String, Vector<Integer>>();
          chunkLines = new Vector<String>();
        }
      }
      
      if (sites.size() > 0) {
        refName = s[0]; end = Integer.parseInt(s[1]);
        computeChunk(refName, start, end, chunkLines);
      }
    } catch(Exception e) {
      bailOut(e);
    }
  }

  // Read the metadata and right before the header insert the NEW INFO, then 
  // save the bam paths and print the actual header line
  public void processHeader() {
    int col_first_path = 9; // path to first merged bam
    try {
      String line;
      while ((line = br.readLine()) != null && line.length() != 0) {
        char[] a_line = line.toCharArray();
        if (a_line[0] == '#' && a_line[1] != '#') { // line is header
          String[] s = line.split("\t");
          for (int i=col_first_path; i<s.length; i++) // Save bam paths
            bamPaths.add(s[i]);
          // Print the NEW INFO line and the actually header and we are done
          System.out.println(NEW_INFO);
          System.out.println(line);
          break;
        } else {
          System.out.println(line); // print the metadata lines
        }
      }
    } catch(Exception e) {
      bailOut(e);
    }
  }
  
  public void add(String chrm, String coor) {
    Vector<Integer> v = new Vector<Integer>(bamPaths.size());
    for (int i=0; i<bamPaths.size(); i++)
      v.add(0);
    sites.put(chrm + "_" + coor, v);
  }
 
  public void plusOne(int bamNumber, String k) {
    Vector <Integer> v = sites.get(k);
    if (v.size() == 0)
      bailOut(new Exception("Trying to access an unexpected site: " + k));
    int old = v.elementAt(bamNumber);
    v.removeElementAt(bamNumber);
    v.add(bamNumber, old + 1); 
  }
 
  public boolean containsKey(String k) {
    return sites.containsKey(k);
  }
 
  public int size() {
    return sites.size();
  }
 
  private void openVCF() {
    try {
      in = new FileInputStream(vcfFn);
      br = new BufferedReader(new InputStreamReader(in));
    } catch(Exception e) {
      bailOut(e);
    }   
  }
 
  private void closeVCF() {
    try {
      in.close();
      br.close();
    } catch(Exception e) {
      bailOut(e);
    }
  }
  
  private void bailOut(Exception e) {
    System.err.println(e.getMessage());
    e.printStackTrace();
    System.exit(1);
  }
}
