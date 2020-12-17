package ghsmith.fastqmutator;

/**
 *
 * @author ghsmith
 */
public class Insert {
    
    public String chromosome;
    public Integer position;
    public String sequence;

    public Insert(String chromosome, Integer position, String sequence) {
        this.chromosome = chromosome;
        this.position = position;
        this.sequence = sequence;
    }

    @Override
    public String toString() {
        return "{" + chromosome + ", " + position + ", " + sequence + '}';
    }
    
}
