import java.io.*;
import java.util.*;
import net.sf.samtools.*;

public class BamCoverage {
  private int n_threads;
  private Vector<String> bams;
  private SnpSites sites;

  // Create the SAMreaders for all the bams
  public BamCoverage(Vector<String> bPaths, SnpSites s, int n_threads) {
    sites = s;
    bams = bPaths;
    this.n_threads = n_threads;
  }

  // find the coverage in all the sites for all the bams concurrently
  public void findCoverage(String rn, int start, int end) {
    ArrayList<Thread> running = new ArrayList<Thread>();
    ArrayList<Thread> waiting = new ArrayList<Thread>();

    System.err.println("Starting Threads ... [" + Integer.toString(n_threads) + "]");
    for (int i=0; i<bams.size(); i++) {
      if (i < n_threads) {
        Thread t = new Worker(i, rn, start, end);
        t.start();
        running.add(t);
      }
      else {
        waiting.add(new Worker(i, rn, start, end));
      }
    }

    // Poll the status of the running threads and start new ones as
    // necessary
    while (running.size() > 0) {
      try { Thread.sleep(1000); }
      catch (InterruptedException ignore) {}
      for (int i=0; i<running.size(); i++) {
        if (running.get(i).getState() == java.lang.Thread.State.TERMINATED) {
          running.remove(i);
          if (waiting.size() > 0) {
            System.err.println("Starting new thread");
            Thread t = waiting.remove(0);
            t.start();
            running.add(t);
          }
        }
      }
    }
  } 
 
  // Threads that will traverse the bams and save the coverage for the sites
  class Worker extends Thread {
    private int bamNumber, start, end;
    private String refName;
    public Worker(int i, String refName, int start, int end) { 
      this.bamNumber = i;
      this.refName = refName;
      this.start = start;
      this.end = end;
    }

    // Get all the raw coverage for all the samples for all the SNP sites
    public void run() {
      System.err.println("Working on bam: " + bams.get(bamNumber));
      File inputFile       = new File(bams.get(bamNumber));
      SAMFileReader reader = new SAMFileReader(inputFile);
      SAMRecord r;
      for (SAMRecordIterator it = reader.query(refName, start, end, false); it.hasNext();) {
        r = it.next();
        if (!r.getReadUnmappedFlag() &&
             r.getAlignmentStart() >= start && r.getAlignmentEnd() <= end) { // We may have bases from a read outside
          for (int j=r.getAlignmentStart(); j<=r.getAlignmentEnd(); j++) {
            String k = refName + "_" + Integer.toString(j);
            if (sites.containsKey(k))
              sites.plusOne(bamNumber, k);
          }
        }
      }
      reader.close();
      System.err.println("Done with: " + bams.get(bamNumber));
    }
  }
}
