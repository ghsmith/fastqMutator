package ghsmith.fastqmutator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
        inserts.add(new Insert("chr7", 140534508, "ATCTTTTTT"));
    }
    
    public static void main(final String[] args) throws IOException {

        List<SamRecord> samRecords = new ArrayList<>();
        
        {
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            String samLine;
            while((samLine = reader.readLine()) != null) {
                SamRecord samRecord = new SamRecord(samLine, inserts);
                if(!samRecord.isProperlyAligned()) {
                    System.out.println(String.format("read %s : rejected (%s)", samRecord.qname, "!isProperlyAligned"));
                    continue;
                }
                if(!samRecord.cigar.matches("[0-9]*M")) {
                    System.out.println(String.format("read %s : rejected (%s)", samRecord.qname, "CIGAR = " + samRecord.cigar));
                    continue;
                }
                if(samRecord.insert == null) {
                    System.out.println(String.format("read %s : rejected (%s)", samRecord.qname, "no applicable insert found"));
                    continue;
                }
                assert samRecord.simpleLength() == samRecord.seq.length();
                System.out.println(String.format("read %s : accepted for insert %s", samRecord.qname, samRecord.insert));
                samRecords.add(samRecord);
            }
            System.out.println(String.format("%d reads accepted", samRecords.size()));
        }
        
        final int[] readCountTotal = new int[] { samRecords.size() };
        final int[] readCountMatched = new int[] { 0 };
        List<Thread> threads = new ArrayList<>();
        
        for(int x = 0; x < args.length; x++) {
            final String final_fastqFileName = args[x];
            Thread t = new Thread() {
                public String fastqFileName = final_fastqFileName;
                @Override
                public void run() {
                    try {
                        Boolean read1File = fastqFileName.contains("_R1_");
                        Boolean read2File = fastqFileName.contains("_R2_");
                        BufferedReader fastqReader = new BufferedReader(new FileReader(fastqFileName), 536870912);
                        BufferedWriter fastqWriter = new BufferedWriter(new FileWriter(new File("mutated-" + fastqFileName)), 536870912);
                        String[] fastqRecord = new String[4];
                        while((fastqRecord[0] = fastqReader.readLine()) != null) {
                            fastqRecord[1] = fastqReader.readLine();
                            fastqRecord[2] = fastqReader.readLine();
                            fastqRecord[3] = fastqReader.readLine();
                            for(SamRecord samRecord : samRecords) {
                                if((samRecord.isFirstSegment() && read1File) || (samRecord.isLastSegment() && read2File)) {
                                    if(fastqRecord[0].substring(1, fastqRecord[0].indexOf(' ')).equals(samRecord.qname)) {
                                        assert samRecord.seq.equals(samRecord.isReverseComplement() ? reverseComplement(fastqRecord[1]) : fastqRecord[1]);
                                        String newSequence = String.format("%s%s%s",
                                            samRecord.seq.substring(0, samRecord.insert.position - samRecord.pos),
                                            samRecord.insert.sequence,
                                            samRecord.seq.substring(samRecord.insert.position - samRecord.pos)
                                        );
                                        newSequence = samRecord.isReverseComplement() ? reverseComplement(newSequence) : newSequence;
                                        String newSequenceForPrinting = String.format("%s%s%s",
                                            samRecord.seq.substring(0, samRecord.insert.position - samRecord.pos),
                                            samRecord.insert.sequence,
                                            samRecord.seq.substring(samRecord.insert.position - samRecord.pos)
                                        );
                                        newSequenceForPrinting = samRecord.isReverseComplement() ? reverseComplement(newSequenceForPrinting) : newSequenceForPrinting;
                                        String oldSequenceForPrinting = String.format("%s%" + samRecord.insert.sequence.length() + "s%s",
                                            samRecord.seq.substring(0, samRecord.insert.position - samRecord.pos),
                                            " ",
                                            samRecord.seq.substring(samRecord.insert.position - samRecord.pos)
                                        );
                                        oldSequenceForPrinting = samRecord.isReverseComplement() ? reverseComplement(oldSequenceForPrinting) : oldSequenceForPrinting;
                                        Integer trimLeft = 0;
                                        Integer trimRight = 0;
                                        while(newSequence.length() > samRecord.simpleLength()) {
                                            if(Math.random() < 0.5) {
                                                newSequence = newSequence.substring(1);
                                                trimLeft++;
                                            }
                                            else {
                                                newSequence = newSequence.substring(0, newSequence.length() - 1);
                                                trimRight++;
                                            }
                                        }
                                        synchronized(readCountMatched) {
                                            System.out.println(String.format("[%5d/%5d] %s : %s", ++readCountMatched[0], readCountTotal[0], fastqFileName, samRecord.qname));
                                            System.out.println(String.format("...rname   : %s", samRecord.rname));
                                            System.out.println(String.format("...position: %s", samRecord.pos));
                                            System.out.println(String.format("...rev comp: %s", samRecord.isReverseComplement()));
                                            System.out.println(String.format("...first   : %s", samRecord.isFirstSegment()));
                                            System.out.println(String.format("...last    : %s", samRecord.isLastSegment()));
                                            System.out.println(String.format("...old seq : %s", oldSequenceForPrinting));
                                            System.out.println(String.format("...new seq : %s", newSequenceForPrinting));
                                            System.out.print  (              "...trimming: ");
                                            for(int x = 0; x < trimLeft; x++) { System.out.print("-"); }
                                            for(int x = 0; x < samRecord.simpleLength(); x++) { System.out.print("+"); }
                                            for(int x = 0; x < trimRight; x++) { System.out.print("-"); }
                                            System.out.println();
                                        }
                                        fastqRecord[1] = newSequence;
                                    }
                                }
                            }
                            fastqWriter.write(fastqRecord[0]); fastqWriter.newLine();
                            fastqWriter.write(fastqRecord[1]); fastqWriter.newLine();
                            fastqWriter.write(fastqRecord[2]); fastqWriter.newLine();
                            fastqWriter.write(fastqRecord[3]); fastqWriter.newLine();
                        }
                        fastqReader.close();
                        fastqWriter.close();
                    }
                    catch(Exception e) {
                        throw new RuntimeException(e);
                    }
                }
            };
            t.start();
        }
        
    }

    public static Map<Character, Character> complementMap = new HashMap<>();
    static {
        complementMap.put('A', 'T');
        complementMap.put('C', 'G');
        complementMap.put('G', 'C');
        complementMap.put('T', 'A');
        complementMap.put(' ', ' ');
    }
    public static String complement(String sequence) {
        StringBuilder complementSequence = new StringBuilder();
        for(int x = 0; x < sequence.length(); x++) {
            assert complementMap.keySet().contains(sequence.charAt(x));
            complementSequence.append(complementMap.get(sequence.charAt(x)));
        }
        return complementSequence.toString();
    }
    public static String reverseComplement(String sequence) {
        return new StringBuilder(complement(sequence)).reverse().toString();
    }
    
}
