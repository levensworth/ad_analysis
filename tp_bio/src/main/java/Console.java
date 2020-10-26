import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.GenbankReaderHelper;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.List;

public class Console {

    public static final String FILE_PATH = "/Users/levensworth/IdeaProjects/tp_bio/ale.gb.txt";


    public static void main(String[] args) {
        /*
        * read file genbank
        * parse orf from file
        * get orf into a record
        * write record as fasta formatted file
        *
         * */



        try {
            File genbankFile = new File(FILE_PATH);

            LinkedHashMap<String, DNASequence> sequences =
                    GenbankReaderHelper.readGenbankDNASequence( genbankFile );


            sequences.forEach((String name, DNASequence seq) -> System.out.println(seq.getSequenceAsString()));

        } catch (Exception e) {
            e.printStackTrace();
        }


    }
}
