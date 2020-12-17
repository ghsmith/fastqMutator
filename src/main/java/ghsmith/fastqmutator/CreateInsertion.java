package ghsmith.fastqmutator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
                if(samRecord.insert == null) {
                    System.out.println(String.format("%s : rejected (%s)", samRecord.qname, "no applicable insert found"));
                    continue;
                }
                assert samRecord.simpleLength() == samRecord.seq.length();
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
                        assert fastqRecord[1].length() == samRecord.simpleLength();
                        String newSequence = String.format("%s%s%s",
                            fastqRecord[1].substring(0, samRecord.insert.position - samRecord.pos),
                            samRecord.isReverseComplement() ? reverseComplement(samRecord.insert.sequence) : samRecord.insert.sequence,
                            fastqRecord[1].substring(samRecord.insert.position)
                        );
                        while(newSequence.length() > samRecord.simpleLength()) {
                            if(Math.random() < 0.5) {
                                newSequence = newSequence.substring(1);
                            }
                            else {
                                newSequence = newSequence.substring(0, newSequence.length() - 1);
                            }
                        }
                        System.out.println(String.format("...old: %s", fastqRecord[1]));
                        System.out.println(String.format("...new: %s", newSequence));
                    }
                }
            }
            fastqReader.close();
            fastqWriter.close();
        }
        
    }

    public static Map<Character, Character> complementMap = new HashMap<>();
    static {
        complementMap.put('A', 'T');
        complementMap.put('C', 'G');
        complementMap.put('G', 'C');
        complementMap.put('T', 'A');
    }
    public static String reverseComplement(String sequence) {
        StringBuilder reverseComplementSequence = new StringBuilder();
        for(int x = sequence.length() - 1; x >= 0; x++) {
            assert complementMap.keySet().contains(sequence.charAt(x));
            reverseComplementSequence.append(complementMap.get(sequence.charAt(x)));
        }
        return reverseComplementSequence.toString();
    }
    
}
