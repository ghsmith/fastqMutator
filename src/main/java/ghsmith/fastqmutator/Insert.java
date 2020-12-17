package ghsmith.fastqmutator;

/**
 *
 * @author ghsmith
 */
public class Insert {
    
    public String chromosome;
    public Integer position;
    public String insertionSequence;

    public Insert(String chromosome, Integer position, String insertionSequence) {
        this.chromosome = chromosome;
        this.position = position;
        this.insertionSequence = insertionSequence;
    }

    @Override
    public String toString() {
        return "{" + chromosome + ", " + position + ", " + insertionSequence + '}';
    }
    
}
