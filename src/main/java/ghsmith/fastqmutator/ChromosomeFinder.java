package ghsmith.fastqmutator;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author ghsmith
 */
public class ChromosomeFinder {
    
    public String genomeFileName;
    public Map<String, Chromosome> chromosomeMap = new HashMap<>();

    public ChromosomeFinder(String genomeFileName) {
        this.genomeFileName = genomeFileName;
    }
    
    public Map<String, Chromosome> getChromosomeMap() throws FileNotFoundException, IOException {
        if(chromosomeMap == null) {
            BufferedReader faReader = new BufferedReader(new FileReader(genomeFileName), 52428800);
            Chromosome chromosome = null;
            StringBuilder sequence = null;
            String faLine;
            while((faLine = faReader.readLine()) != null) {
                if(faLine.charAt(0) == '>') {
                    if(chromosome != null) {
                        chromosome.sequence = sequence.toString().toCharArray();
                    }
                    chromosome = new Chromosome(faLine.substring(1));
                    chromosomeMap.put(chromosome.name, chromosome);
                    continue;
                }
                sequence.append(faLine);
            }
            chromosome.sequence = sequence.toString().toCharArray();
            faReader.close();
        }
        return chromosomeMap;
    }
    
}
