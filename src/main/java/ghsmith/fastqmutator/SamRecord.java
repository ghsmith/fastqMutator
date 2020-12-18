package ghsmith.fastqmutator;

import java.util.List;

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

    public String lane;
    
    public Insert insert;
    
    public SamRecord(String samLine, List<Insert> inserts) {
        
        String[] fields = samLine.split("\t");
        qname = fields[0];
        flag = Integer.valueOf(fields[1]);
        rname = fields[2];
        pos = Integer.valueOf(fields[3]);
        cigar = fields[5];
        seq = fields[9];
        
        lane = qname.split(":")[3];
        
        for(Insert insert : inserts) {
            if(
                this.simpleLength() != null
                && insert.chromosome.equals(this.rname)
                && insert.position >= this.pos
                && insert.position < this.pos + this.simpleLength()
            ) {
                this.insert = insert;
            }
        }
        
    }
    
    public boolean isProperlyAligned() {
        return (flag & 0x2) > 0;
    }

    public boolean isReverseComplement() {
        return (flag & 0x10) > 0;
    }

    public boolean isFirstSegment() {
        return (flag & 0x40) > 0;
    }

    public boolean isLastSegment() {
        return (flag & 0x80) > 0;
    }
    
    public Integer simpleLength() {
        if(cigar.matches("[0-9]*M")) {
            return Integer.valueOf(cigar.substring(0, cigar.length() - 1));
        }
        return null;
    }
    
}
