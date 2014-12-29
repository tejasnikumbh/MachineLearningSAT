package ca.ubc.cs.beta.models.fastrf;

import java.util.*;
import java.io.*;
import ca.ubc.cs.beta.models.fastrf.utils.Utils;

public class PairwiseClassificationForests implements java.io.Serializable {
    private static final long serialVersionUID = 1593817496029385921L;
    
    public int numForests;
    private RandomForest matrix[][];
    private boolean reverse[][];
    
    // Transforms these forests were built with. Set these manually if you want to use them.
    public int[] keptColumns;
    public double[] scale;
    public double[] bias;
    
    public PairwiseClassificationForests(int numForests) {
        if (numForests <= 0) {
            throw new RuntimeException("Invalid number of forests in PairwiseClassificationForests: " + numForests);
        }
        this.numForests = numForests;
        matrix = new RandomForest[numForests][];
        reverse = new boolean[numForests][];
        for (int i=0; i < numForests; i++) {
            matrix[i] = new RandomForest[numForests-1-i];
            reverse[i] = new boolean[numForests-1-i];
            for (int j=0; j < numForests-1-i; j++) {
                matrix[i][j] = null;
            }
        }
    }
    
    /**
     * Adds a new model between a pair a and b.
     * a and b should be in [0, numForests).
     */
    public void add(RandomForest forest, int a, int b) {
        if (a == b) {
            throw new RuntimeException("a == b: " + a);
        }
        boolean reversed = false;
        if (a > b) {
            int t = a;
            a = b;
            b = t;
            reversed = true;
        }
        b = b - 1 - a;
        if (matrix[a][b] != null) {
            throw new RuntimeException("Adding a pair that already exists: " + a + ", " + b+a+1);
        }
        matrix[a][b] = forest;
        reverse[a][b] = reversed;
    }
    
    /**
     * Checks to see if all forests have been added, and throws error if not.
     */
    public void checkInputs() {
        for (int i=0; i < numForests; i++) {
            for (int j=0; j < numForests-1-i; j++) {
                if (matrix[i][j] == null) {
                    throw new RuntimeException("Missing pair: " + i + ", " + j+i+1);
                }
            }
        }
    }
    
    /**
     * Classifies the given instantiations of features
     * Assumes that the classes for the forests are either 1 or 2,
     * Where 1 means b wins, and 2 means a wins, where a and b are from
     * the call to add(forest, a, b).
     *
     * @returns a matrix of size X.length where index i contains the 
     * most popular response for X[i]
     */
    public static double[] classify(PairwiseClassificationForests PCF, double[][] X) {
        double[][] votes = new double[X.length][PCF.numForests*(PCF.numForests-1)/2];
        for (int i=0, count=0; i < PCF.numForests; i++) {
            for (int j=0; j < PCF.numForests-1-i; j++) {
                double[] res = RandomForest.classify(PCF.matrix[i][j], X);
                for (int k=0; k < res.length; k++) {
                    votes[k][count] = (res[k] == 1 + (PCF.reverse[i][j] ? 1 : 0)) ? j+i+1 : i;
                }
                count++;
            }
        }
        double[] retn = new double[X.length];
        for (int i=0; i < X.length; i++) {
            double[] best = Utils.mode(votes[i]);
            if (best.length == 1) {
                retn[i] = best[0];
            } else {
                // Multiple bests, reclassify only between them
                double[] faceoff = new double[best.length*(best.length-1)/2];
                double[][] input = new double[1][];
                input[0] = X[i];
                for (int j=0, count=0; j < best.length; j++) {
                    for (int k=j+1; k < best.length; k++) {
                        int a = (int)best[j], b = (int)best[k];
                        if (a > b) {
                            int t = a;
                            a = b;
                            b = t;
                        }
                        b = b - 1 - a;
                        int res = (int)(RandomForest.classify(PCF.matrix[a][b], input)[0]);
                        faceoff[count++] = (res == 1 + (PCF.reverse[a][b] ? 1 : 0)) ? a+b+1 : a;
                    }
                }
                best = Utils.mode(faceoff);
                retn[i] = best[(int)(Math.random()*best.length)];
            }
        }
        return retn;
    }
     
    /**
     * Ranks the forests based on the given instantiations of features
     * Assumes that the classes for the forests are either 1 or 2,
     * Where 1 means b wins, and 2 means a wins, where a and b are from
     * the call to add(forest, a, b).
     *
     * If freqout is not null, then the frequencies of the forests for each row of X is given there.
     *
     * @returns a matrix of size X.length * PCF.numForests where index i contains the 
     * response for X[i] ranked from most popular to least popular.
     */
    public static double[][] rank(PairwiseClassificationForests PCF, double[][] X, ArrayList<HashMap<Double, Integer> > freqout) {
        double[][] votes = new double[X.length][PCF.numForests*(PCF.numForests-1)/2];
        for (int i=0, count=0; i < PCF.numForests; i++) {
            for (int j=0; j < PCF.numForests-1-i; j++) {
                double[] res = RandomForest.classify(PCF.matrix[i][j], X);
                for (int k=0; k < res.length; k++) {
                    votes[k][count] = (res[k] == 1 + (PCF.reverse[i][j] ? 1 : 0)) ? j+i+1 : i;
                }
                count++;
            }
        }
        double[][] retn = new double[X.length][];
        for (int i=0; i < X.length; i++) {
            HashMap<Double, Integer> freq = Utils.getFreq(votes[i]);
            if (freqout != null) {
                freqout.add(freq);
            }
            Double[] d = new Double[PCF.numForests];
            for (int j=0; j < PCF.numForests; j++) d[j] = new Double(j);
            TreeMap<Integer, ArrayList<Double> > freqinv = getInvFreq(freq, d); 

            double[] reti = new double[PCF.numForests];
            int numProcessed = 0;

            for (ArrayList<Double> v : freqinv.values()) {
                if (v.size() == 1) {
                    reti[numProcessed++] = v.get(0);
                } else {
                    // Multiple forests with the same frequency. reclassify only between them
                    double[] faceoff = new double[v.size()*(v.size()-1)/2];
                    double[][] input = new double[1][];
                    input[0] = X[i];
                    for (int j=0, count=0; j < v.size(); j++) {
                        for (int k=j+1; k < v.size(); k++) {
                            int a = v.get(j).intValue(), b = v.get(k).intValue();
                            if (a > b) {
                                int t = a;
                                a = b;
                                b = t;
                            }
                            b = b - 1 - a;
                            int res = (int)(RandomForest.classify(PCF.matrix[a][b], input)[0]);
                            faceoff[count++] = (res == 1 + (PCF.reverse[a][b] ? 1 : 0)) ? a+b+1 : a;
                        }
                    }
               
                    TreeMap<Integer, ArrayList<Double> > subinvfreq = getInvFreq(Utils.getFreq(faceoff), v.toArray(new Double[0]));
                    for (ArrayList<Double> subv : subinvfreq.values()) {
                        Collections.shuffle(subv);
                        for (int j=0; j < subv.size(); j++) {
                            reti[numProcessed++] = subv.get(j);
                        }
                    }
                }
            }
            retn[i] = reti;
        }
        return retn;
    }

    public static double[][] rank(PairwiseClassificationForests PCF, double[][] X) {
        return rank(PCF, X, null);
    }

    private static <K> TreeMap<Integer, ArrayList<K> > getInvFreq(HashMap<K,Integer> freq, K[] values) {
        TreeMap<Integer, ArrayList<K> > freqinv = new TreeMap<Integer, ArrayList<K> >();
        for (int j=0; j < values.length; j++) {
            K k = values[j];
            if (freq.get(k) == null) freq.put(k, 0);
            Integer f = -freq.get(k);
            if (freqinv.get(f) == null) freqinv.put(f, new ArrayList<K>());
            freqinv.get(f).add(k);
        }
        return freqinv;
    }

    /** 
     * Prepares the PCF for classification
     * dropysub specifies whether the ysub array should be cleared after preprocessing
     * @see RandomForest.preprocessForestForClassification
     */
    public static void preprocess(PairwiseClassificationForests PCF, boolean dropysub) {
        for (int i=0; i < PCF.numForests; i++) {
            for (int j=0; j < PCF.numForests-1-i; j++) {
                RandomForest.preprocessForestForClassification(PCF.matrix[i][j], dropysub);
            }
        }
    }
    
    public static void serialize(PairwiseClassificationForests PCF, String filename) {   
	    try {
			FileOutputStream fileOut = new FileOutputStream(filename);
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(PCF);
			out.close();
			fileOut.close();
		} catch (IOException i) {
			i.printStackTrace();
		}
    }
    
    public static PairwiseClassificationForests deserialize(String filename) { 
        PairwiseClassificationForests PCF = null;
	    try {
	        FileInputStream fileIn =
	                      new FileInputStream(filename);
	        ObjectInputStream in = new ObjectInputStream(fileIn);
	        try {
				PCF = (PairwiseClassificationForests) in.readObject();
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			}
	        in.close();
	        fileIn.close();
	    } catch(IOException i) {
	        i.printStackTrace();
	    }
        return PCF;
    }
}
