import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Calin Popa 20158726
 * Ming-Xia Delvas 20104038
 */
public class ReadFastq {

    // pénalitées
    static int match = 4;
    static int mismatch = 4;
    static int indel = 8;

    // initialisation de la matrice de scores
    static int[][] matScore;

    // initialisation des indices pour les scores maximums
    static int max_i;
    static int max_j;

    static int longueurChev; 

    // intialisation d une liste contenant les sequences du fichier reads
    static ArrayList<String> sequence = new ArrayList<String>();

    // fonction qui retourne le score
    public static int scoreMax(String s1, String s2) {
        matScore = new int[s1.length() + 1][s2.length() + 1];

        // initialization des premieres lignes et colonnes de la matrice des scores
        for (int i = 1; i <= s1.length(); i++) {
            matScore[i][0] = 0;
        }
        for (int j = 1; j <= s2.length(); j++) {
            matScore[0][j] = matScore[0][j - 1] - indel;
        }
        // remplissage de la matrice de scores
        for (int i = 1; i <= s1.length(); i++) {
            for (int j = 1; j <= s2.length(); j++) {
                int scoreDiag = matScore[i - 1][j - 1] +
                        (s1.charAt(i - 1) == s2.charAt(j - 1) ? match : // meme symbol
                                -mismatch); // symbol different
                int scoreLeft = matScore[i][j - 1] - indel; // insertion
                int scoreUp = matScore[i - 1][j] - indel; // deletion
                // on prend le max
                matScore[i][j] = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
            }
        }

        int score = 0;
        max_i = s1.length();
        max_j = 0;
        for (int i = 0; i <= s2.length(); i++) {
            if (score < matScore[max_i][i]) {
                score = matScore[max_i][i];
                max_j = i;
            }
        }

        return score;
    }

    // fonction qui retourne une paire de sequence
    public static String alignerSeq(String s1, String s2) {
        scoreMax(s1, s2);

        // initialization de l'alignement
        String alignS1 = "";
        String alignS2 = "";

        // indices pour la boucle while
        int i = max_i;
        int j = max_j;

        // Boucle pour retourner a la case de départ
        // avec un alignement au final.
        while (j > 0) {
            if (i > 0) {
                int scoreDiag = matScore[i - 1][j - 1] +
                        (s1.charAt(i - 1) == s2.charAt(j - 1) ? match : // meme symbol
                                -mismatch); // symbol different
                int scoreLeft = matScore[i][j - 1] - indel; // insertion
                int scoreUp = matScore[i - 1][j] - indel; // deletion

                if (Math.max(Math.max(scoreDiag, scoreLeft), scoreUp) == scoreDiag) {
                    alignS1 = s1.charAt(i - 1) + alignS1;
                    alignS2 = s2.charAt(j - 1) + alignS2;
                    i--;
                    j--;
                }
                if (Math.max(Math.max(scoreDiag, scoreLeft), scoreUp) == scoreUp) {
                    alignS1 = s1.charAt(i - 1) + alignS1;
                    alignS2 = "-" + alignS2;
                    i--;
                }
                if (Math.max(Math.max(scoreDiag, scoreLeft), scoreUp) == scoreLeft) {
                    alignS1 = "-" + alignS1;
                    alignS2 = s2.charAt(j - 1) + alignS2;
                    j--;
                }
            } else {
                alignS1 = "-" + alignS1;
                alignS2 = s2.charAt(j - 1) + alignS2;
                j--;
            }
        }

        // Calcul de la longueur de chevauchement
        longueurChev = alignS1.length();

        // resultat des alignements que la fonction retourne
        String alignement = "alignement:\n" + alignS1 + "\n" + alignS2;
        return alignement;

    }

    public static String printResult(String s1, String s2) {
        // resultat que la fonction retourne
        String resultat = "score: " + scoreMax(s1, s2) + "\n" + alignerSeq(s1, s2)
                + "\n"
                + "longueur du chevauchement: " + longueurChev + "\n";
        System.out.println(resultat);
        return resultat;
    }

    public static ArrayList<String> readFastq(String file) throws IOException {
        int lineNumber = 1;

        try {
            // lire le fichier fastq
            BufferedReader bufferedReader = new BufferedReader(new FileReader(file));
            String line = null;
            // ajoute chaque sequence du fichier read dans un arrayList
            while ((line = bufferedReader.readLine()) != null) {
                if (lineNumber % 4 == 2) {
                    sequence.add(line);
                }
                lineNumber++;
            }

            for (var i = 0; i < 20; i++) {
                for (var j = 0; j < 20; j++) {
                   printResult(sequence.get(i), sequence.get(j));
                }
            }

            bufferedReader.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        return sequence;
    }

    public static Integer[][] assemblageFrag() {
        Integer matrix[][] = new Integer[20][20];
 
        // calcule la matrice des score maximales des chevauchements
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 20; j++) {
                 if (i != j) {
                    matrix[i][j] = scoreMax(sequence.get(i), sequence.get(j));
                } else {
                    matrix[i][j] = null;
                }
                System.out.print(matrix[i][j] +  "\t");
            }
            System.out.println();
        }

        return matrix;
    }

    public static String cadresAvecAA() throws IOException{
        File file=new File("./sequence.fasta");
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String st;
            String cadre1="";
            String AA1="";
            String cadre2;
            String AA2="";
            String cadre3;
            String AA3="";
            int bool = 0;
            while ((st = br.readLine()) != null){
                if(bool>0){
                    cadre1=cadre1 + st;
                    
                }
                bool=1;
            }
            
            cadre2=cadre1.substring(1, cadre1.length()-3);
            cadre3=cadre1.substring(2,cadre1.length()-2);
            for(int i=0; i<cadre1.length()/3;i=i+3){
                String acide=String.valueOf(cadre1.charAt(i))+String.valueOf(cadre1.charAt(i+1))+String.valueOf(cadre1.charAt(i+2));
                
                //System.out.println(acide.equals("AGG"));
                if(acide.equals("TTT")||acide.equals("TTC")){
                    AA1=AA1+"F";
                }
                if(acide.equals("TTA")||acide.equals("TTG")||acide.equals("CTT")||acide.equals("CTC")||acide.equals("CTA")||acide.equals("CTG")){
                    AA1=AA1+"L";
                }
                if(acide.equals("ATT")||acide.equals("ATC")||acide.equals("ATA")){
                    AA1=AA1+"I";
                }
                if(acide.equals("ATG")){
                    AA1=AA1+"M";
                }
                if(acide.equals("GTT")||acide.equals("GTC")||acide.equals("GTA")||acide.equals("GTG")){
                    AA1=AA1+"V";
                }
                if(acide.equals("TCT")||acide.equals("TCC")||acide.equals("TCA")||acide.equals("TCG")){
                    AA1=AA1+"S";
                }
                if(acide.equals("CCT")||acide.equals("CCC")||acide.equals("CCA")||acide.equals("CCG")){
                    AA1=AA1+"P";
                }
                if(acide.equals("ACT")||acide.equals("ACC")||acide.equals("ACA")||acide.equals("ACG")){
                    AA1=AA1+"T";
                }
                if(acide.equals("GCT")||acide.equals("GCC")||acide.equals("GCA")||acide.equals("GCG")){
                    AA1=AA1+"A";
                }
                if(acide.equals("TAT")||acide.equals("TAC")){
                    AA1=AA1+"Y";
                }
                if(acide.equals("TAA")||acide.equals("TAG")||acide.equals("TGA")){
                    AA1=AA1+"*";
                }
                if(acide.equals("CAT")||acide.equals("CAC")){
                    AA1=AA1+"H";
                }
                if(acide.equals("CAA")||acide.equals("CAG")){
                    AA1=AA1+"Q";
                }
                if(acide.equals("AAT")||acide.equals("AAC")){
                    AA1=AA1+"N";
                }
                if(acide.equals("AAA")||acide.equals("AAG")){
                    AA1=AA1+"K";
                }
                if(acide.equals("GAT")||acide.equals("GAC")){
                    AA1=AA1+"D";
                }
                if(acide.equals("GAA")||acide.equals("GAG")){
                    AA1=AA1+"E";
                }
                if(acide.equals("TGT")||acide.equals("TGC")){
                    AA1=AA1+"C";
                }
                if(acide.equals("TGG")){
                    AA1=AA1+"W";
                }
                if(acide.equals("CGT")||acide.equals("CGC")||acide.equals("CGA")||acide.equals("CGG")||acide.equals("AGA")||acide.equals("AGG")){
                    AA1=AA1+"R";
                }
                if(acide.equals("AGT")||acide.equals("AGC")){
                    AA1=AA1+"S";
                }
                if(acide.equals("GGT")||acide.equals("GGC")||acide.equals("GGA")||acide.equals("GGG")){
                    AA1=AA1+"G";
                }
                
            }
            
            for(int i=0; i<cadre2.length()/3;i=i+3){
                String acide=""+cadre2.charAt(i)+cadre2.charAt(i+1)+cadre2.charAt(i+2);
                if(acide.equals("TTT")||acide.equals("TTC")){
                    AA2=AA2+"F";
                }
                if(acide.equals("TTA")||acide.equals("TTG")||acide.equals("CTT")||acide.equals("CTC")||acide.equals("CTA")||acide.equals("CTG")){
                    AA2=AA2+"L";
                }
                if(acide.equals("ATT")||acide.equals("ATC")||acide.equals("ATA")){
                    AA2=AA2+"I";
                }
                if(acide.equals("ATG")){
                    AA2=AA2+"M";
                }
                if(acide.equals("GTT")||acide.equals("GTC")||acide.equals("GTA")||acide.equals("GTG")){
                    AA2=AA2+"V";
                }
                if(acide.equals("TCT")||acide.equals("TCC")||acide.equals("TCA")||acide.equals("TCG")){
                    AA2=AA2+"S";
                }
                if(acide.equals("CCT")||acide.equals("CCC")||acide.equals("CCA")||acide.equals("CCG")){
                    AA2=AA2+"P";
                }
                if(acide.equals("ACT")||acide.equals("ACC")||acide.equals("ACA")||acide.equals("ACG")){
                    AA2=AA2+"T";
                }
                if(acide.equals("GCT")||acide.equals("GCC")||acide.equals("GCA")||acide.equals("GCG")){
                    AA2=AA2+"A";
                }
                if(acide.equals("TAT")||acide.equals("TAC")){
                    AA2=AA2+"Y";
                }
                if(acide.equals("TAA")||acide.equals("TAG")||acide.equals("TGA")){
                    AA2=AA2+"*";
                }
                if(acide.equals("CAT")||acide.equals("CAC")){
                    AA2=AA2+"H";
                }
                if(acide.equals("CAA")||acide.equals("CAG")){
                    AA2=AA2+"Q";
                }
                if(acide.equals("AAT")||acide.equals("AAC")){
                    AA2=AA2+"N";
                }
                if(acide.equals("AAA")||acide.equals("AAG")){
                    AA2=AA2+"K";
                }
                if(acide.equals("GAT")||acide.equals("GAC")){
                    AA2=AA2+"D";
                }
                if(acide.equals("GAA")||acide.equals("GAG")){
                    AA2=AA2+"E";
                }
                if(acide.equals("TGT")||acide.equals("TGC")){
                    AA2=AA2+"C";
                }
                if(acide.equals("TGG")){
                    AA2=AA2+"W";
                }
                if(acide.equals("CGT")||acide.equals("CGC")||acide.equals("CGA")||acide.equals("CGG")||acide.equals("AGA")||acide.equals("AGG")){
                    AA2=AA2+"R";
                }
                if(acide.equals("AGT")||acide.equals("AGC")){
                    AA2=AA2+"S";
                }
                if(acide.equals("GGT")||acide.equals("GGC")||acide.equals("GGA")||acide.equals("GGG")){
                    AA2=AA2+"G";
                }
                
            }
            for(int i=0; i<cadre3.length()/3;i=i+3){
                String acide=""+cadre3.charAt(i)+cadre3.charAt(i+1)+cadre3.charAt(i+2);
                if(acide.equals("TTT")||acide.equals("TTC")){
                    AA3=AA3+"F";
                }
                if(acide.equals("TTA")||acide.equals("TTG")||acide.equals("CTT")||acide.equals("CTC")||acide.equals("CTA")||acide.equals("CTG")){
                    AA3=AA3+"L";
                }
                if(acide.equals("ATT")||acide.equals("ATC")||acide.equals("ATA")){
                    AA3=AA3+"I";
                }
                if(acide.equals("ATG")){
                    AA3=AA3+"M";
                }
                if(acide.equals("GTT")||acide.equals("GTC")||acide.equals("GTA")||acide.equals("GTG")){
                    AA3=AA3+"V";
                }
                if(acide.equals("TCT")||acide.equals("TCC")||acide.equals("TCA")||acide.equals("TCG")){
                    AA3=AA3+"S";
                }
                if(acide.equals("CCT")||acide.equals("CCC")||acide.equals("CCA")||acide.equals("CCG")){
                    AA3=AA3+"P";
                }
                if(acide.equals("ACT")||acide.equals("ACC")||acide.equals("ACA")||acide.equals("ACG")){
                    AA3=AA3+"T";
                }
                if(acide.equals("GCT")||acide.equals("GCC")||acide.equals("GCA")||acide.equals("GCG")){
                    AA3=AA3+"A";
                }
                if(acide.equals("TAT")||acide.equals("TAC")){
                    AA3=AA3+"Y";
                }
                if(acide.equals("TAA")||acide.equals("TAG")||acide.equals("TGA")){
                    AA3=AA3+"*";
                }
                if(acide.equals("CAT")||acide.equals("CAC")){
                    AA3=AA3+"H";
                }
                if(acide.equals("CAA")||acide.equals("CAG")){
                    AA3=AA3+"Q";
                }
                if(acide.equals("AAT")||acide.equals("AAC")){
                    AA3=AA3+"N";
                }
                if(acide.equals("AAA")||acide.equals("AAG")){
                    AA3=AA3+"K";
                }
                if(acide.equals("GAT")||acide.equals("GAC")){
                    AA3=AA3+"D";
                }
                if(acide.equals("GAA")||acide.equals("GAG")){
                    AA3=AA3+"E";
                }
                if(acide.equals("TGT")||acide.equals("TGC")){
                    AA3=AA3+"C";
                }
                if(acide.equals("TGG")){
                    AA3=AA3+"W";
                }
                if(acide.equals("CGT")||acide.equals("CGC")||acide.equals("CGA")||acide.equals("CGG")||acide.equals("AGA")||acide.equals("AGG")){
                    AA3=AA3+"R";
                }
                if(acide.equals("AGT")||acide.equals("AGC")){
                    AA3=AA3+"S";
                }
                if(acide.equals("GGT")||acide.equals("GGC")||acide.equals("GGA")||acide.equals("GGG")){
                    AA3=AA3+"G";
                }
                
            }
            String rep=cadre1+"\n\n"+AA1+"\n\n"+cadre2+"\n\n"+AA2+"\n\n"+cadre3+"\n\n"+AA3 +"\n\n"+"Le codon start (M) se trouve dans chaque cadre de lecture";
            System.out.println(rep);
            return rep;
        }

    }

    public static void main(String[] args) throws IOException {
        alignerSeq("CATCCTTCT", "CCTTTCACC");
        try {
            readFastq("./reads.fq");
            assemblageFrag();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        cadresAvecAA();

    }

}
