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
        
        String fastqFileName = args[0];
        
        {
            int matchedReadCount = 0;
            BufferedReader fastqReader = new BufferedReader(new FileReader(fastqFileName));
            PrintWriter fastqWriter = new PrintWriter(new FileWriter(new File(fastqFileName.replaceAll(".fastq$", ".mutated.fastq"))));
            String[] fastqRecord = new String[4];
            while((fastqRecord[0] = fastqReader.readLine()) != null) {
                fastqRecord[1] = fastqReader.readLine();
                fastqRecord[2] = fastqReader.readLine();
                fastqRecord[3] = fastqReader.readLine();
                for(SamRecord samRecord : samRecords) {
                    if(fastqRecord[0].split(" ")[0].equals("@" + samRecord.qname)) {
                        if(!samRecord.seq.equals(samRecord.isReverseComplement() ? reverseComplement(fastqRecord[1]) : fastqRecord[1])) {
                            continue;
                        }
                        System.out.println(String.format("[%5d] %s : %s", ++matchedReadCount, fastqFileName, samRecord.qname));
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
                        System.out.println(String.format("...SAM seq : %s", samRecord.seq));
                        System.out.println(String.format("...FQ seq  : %s", fastqRecord[1]));
                        System.out.println(String.format("...position: %s", samRecord.pos));
                        System.out.println(String.format("...rev comp: %s", samRecord.isReverseComplement()));
                        System.out.println(String.format("...first   : %s", samRecord.isFirstSegment()));
                        System.out.println(String.format("...last    : %s", samRecord.isLastSegment()));
                        System.out.println(String.format("...old seq : %s", oldSequenceForPrinting));
                        System.out.println(String.format("...new new : %s", newSequenceForPrinting));
                        System.out.print  (              "...trimming: ");
                        for(int x = 0; x < trimLeft; x++) { System.out.print("-"); }
                        for(int x = 0; x < samRecord.simpleLength(); x++) { System.out.print("+"); }
                        for(int x = 0; x < trimRight; x++) { System.out.print("-"); }
                        System.out.println();
                        fastqWriter.println(fastqRecord[0]);
                        fastqWriter.println(newSequence);
                        fastqWriter.println(fastqRecord[2]);
                        fastqWriter.println(fastqRecord[3]);
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
        complementMap.put(' ', ' ');
    }
    public static String reverseComplement(String sequence) {
        StringBuilder reverseComplementSequence = new StringBuilder();
        for(int x = sequence.length() - 1; x >= 0; x--) {
            assert complementMap.keySet().contains(sequence.charAt(x));
            reverseComplementSequence.append(complementMap.get(sequence.charAt(x)));
        }
        return reverseComplementSequence.toString();
    }
    
}
