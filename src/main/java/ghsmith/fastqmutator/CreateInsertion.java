package ghsmith.fastqmutator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author ghsmith
 */
public class CreateInsertion {
    
    public static List<Insert> inserts = new ArrayList<>();
    static {
        inserts.add(new Insert("chr7", 140534507, "ATCTTTTTT"));
    }
    
    public static void main(String[] args) throws IOException {

        List<SamRecord> samRecords = new ArrayList<>();
        
        {
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            String samLine;
            while((samLine = reader.readLine()) != null) {
                SamRecord samRecord = new SamRecord(samLine, inserts);
                if(!samRecord.isProperlyAligned()) {
                    System.out.println(String.format("%s : rejected (%s)", samRecord.qname, "!isProperlyAligned"));
                    continue;
                }
                if(!samRecord.cigar.matches("[0-9]*M")) {
                    System.out.println(String.format("%s : rejected (%s)", samRecord.qname, "CIGAR = " + samRecord.cigar));
                    continue;
                }
                System.out.println(String.format("%s : accepted for insert %s", samRecord.qname, samRecord.insert));
                samRecords.add(samRecord);
            }
        }
        
        String fastqFileName = args[0];
        
        {
            BufferedReader fastqReader = new BufferedReader(new FileReader(fastqFileName));
            PrintWriter fastqWriter = new PrintWriter(new FileWriter(new File(fastqFileName.replaceAll(".fastq$", ".mutated.fastq"))));
            String[] fastqRecord = new String[4];
            while((fastqRecord[0] = fastqReader.readLine()) != null) {
                fastqRecord[1] = fastqReader.readLine();
                fastqRecord[2] = fastqReader.readLine();
                fastqRecord[3] = fastqReader.readLine();
                for(SamRecord samRecord : samRecords) {
                    if(fastqRecord[0].substring(1).startsWith(samRecord.qname)) { // substring accounts for "@" at beginning of line
                        // the SAM record is ambiguous w.r.t. R1/R2
                        if(!samRecord.seq.equals(fastqRecord[1])) {
                            continue;
                        }
                        System.out.println(String.format("%s : match in %s", samRecord.qname, fastqFileName));
                    }
                }
            }
            fastqReader.close();
            fastqWriter.close();
        }
        
    }
    
}
