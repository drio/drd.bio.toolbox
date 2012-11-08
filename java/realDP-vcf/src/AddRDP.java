import net.sf.picard.cmdline.*;

public class AddRDP extends CommandLineProgram {
  @Usage
  public String USAGE = getStandardUsagePreamble() +
    "Adds read coverage at site in vcf merged files.\r\n";

  @Option(shortName="VCF", doc="vcf file")
  public String vcfFileName;
  
  @Option(shortName="NT", doc="Number of threads to use.")
  public Integer N_THREADS;
  
  @Option(shortName="CS", doc="Number of vcf entries per chunk.")
  public Integer CHUNK_SIZE;
 
 
  public static void main(String args[]) {
    new AddRDP().instanceMainWithExit(args);
  }
 
  @Override
  protected int doWork() {
    // System.out.println("Working Directory = " + System.getProperty("user.dir"));
    SnpSites ss = new SnpSites(vcfFileName, N_THREADS, CHUNK_SIZE);
    ss.doWork();
    return 0;
  }
}
