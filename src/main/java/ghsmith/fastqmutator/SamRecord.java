package ghsmith.fastqmutator;

/**
 *
 * @author ghsmith
 */
public class SamRecord {

    public String qname;
    public Integer flag;
    public String rname;
    public Integer pos;
    public String cigar;
    public String seq;
    
    public SamRecord(String samLine) {
        String[] fields = samLine.split("\t");
        qname = fields[0];
        flag = Integer.valueOf(fields[1]);
        rname = fields[2];
        pos = Integer.valueOf(fields[3]);
        cigar = fields[5];
        seq = fields[9];
    }
    
    public boolean isProperlyAligned() {
        return (flag & 0x2) > 0;
    }

    public boolean isReverseComplement() {
        return (flag & 0x10) > 0;
    }
    
}
