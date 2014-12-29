package ca.ubc.cs.beta.models.fastrf;

import java.util.*;
import java.io.*;
import ca.ubc.cs.beta.models.fastrf.utils.Utils;

public class RandomForest implements java.io.Serializable {
    private static final long serialVersionUID = 5204746081208095706L;
    
    public int numTrees;
    public Regtree[] Trees;
    
    public int logModel;
    
    // Transforms this forest was built with. Set these manually if you want to use them.
    public int[] keptColumns;
    public double[] scale;
    public double[] bias;
    
    public RandomForest(int numtrees, int logModel) {
        if (numtrees <= 0) {
            throw new RuntimeException("Invalid number of regression trees in forest: " + numtrees);
        }
        this.logModel = logModel;
        numTrees = numtrees;
        Trees = new Regtree[numtrees];
    }
    
    public static RandomForest learnModel(int numTrees, double[][] allTheta, double[][] allX, int[][] theta_inst_idxs, double[] y, RegtreeBuildParams params) {
        Random r = params.random;
        if (r == null) {
            r = new Random();
            if (params.seed != -1) {
                r.setSeed(params.seed);
            }
        }        
        
        int N = y.length;
        
        // Do bootstrap sampling for data for each tree.
        int[][] dataIdxs = new int[numTrees][N];
        for (int i = 0; i < numTrees; i++) {
            for (int j = 0; j < N; j++) {
                dataIdxs[i][j] = r.nextInt(N);
            }
        }
        return learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, dataIdxs, params);
    }
    /**
     * Learns a random forest putting the specified data points into each tree.
     * @see RegtreeFit.fit
     */
    public static RandomForest learnModel(int numTrees, double[][] allTheta, double[][] allX, int[][] theta_inst_idxs, double[] y, int[][] dataIdxs, RegtreeBuildParams params) {       
        if (dataIdxs.length != numTrees) {
            throw new RuntimeException("length(dataIdxs) must be equal to numtrees.");
        }
        
        RandomForest rf = new RandomForest(numTrees, params.logModel);
        for (int i = 0; i < numTrees; i++) {
            int N = dataIdxs[i].length;
            int[][] this_theta_inst_idxs = new int[N][];
            double[] thisy = new double[N];
            for (int j=0; j<N; j++) {
                int idx = dataIdxs[i][j];
                this_theta_inst_idxs[j] = theta_inst_idxs[idx];
                thisy[j] = y[idx];
            }
            rf.Trees[i] = RegtreeFit.fit(allTheta, allX, this_theta_inst_idxs, thisy, params);
        }
        return rf;
    }
    
    /**
     * Propogates data points down the regtree and returns a numtrees*X.length vector of node #s
     * specifying which node each data point falls into.
     * @params X a numdatapoints*numvars matrix
     */
    public static int[][] fwd(RandomForest forest, double[][] X) {
        int[][] retn = new int[forest.numTrees][X.length];
        for (int i=0; i < forest.numTrees; i++) {
            int[] result = RegtreeFwd.fwd(forest.Trees[i], X);
            System.arraycopy(result, 0, retn[i], 0, result.length);
        }
        return retn;
    }
    
    /**
     * Gets a prediction for the given instantiations of features
     * @returns a matrix of size X.length*2 where index (i,0) is the prediction for X[i] 
     * and (i,1) is the variance of that prediction. See Matlab code for how var is calculated.
     */
    public static double[][] apply(RandomForest forest, double[][] X) {
		double[][] retn = new double[X.length][2]; // mean, var
        for (int i=0; i < forest.numTrees; i++) {
            int[] result = RegtreeFwd.fwd(forest.Trees[i], X);
            for (int j=0; j < X.length; j++) {
				double pred = forest.Trees[i].nodepred[result[j]];
                if (forest.logModel>0) {
                    pred = Math.log10(pred);
                }
                retn[j][0] += pred;
                retn[j][1] += forest.Trees[i].nodevar[result[j]]+pred*pred;
            }
        }
        for (int i=0; i < X.length; i++) {
            retn[i][0] /= forest.numTrees;
            retn[i][1] /= forest.numTrees;
            retn[i][1] -= retn[i][0]*retn[i][0];
            retn[i][1] = retn[i][1] * ((forest.numTrees+0.0)/Math.max(1, forest.numTrees-1));
        }
        return retn;
    }
    
    /**
     * Classifies the given instantiations of features
     * @returns a matrix of size X.length where index i contains the 
     * most popular response for X[i]
     */
    public static double[] classify(RandomForest forest, double[][] X) {
        // This currently uses numTrees*X.length memory in order to
        // take advantage of batch forwarding of Xs. 
        // There might be some way of reducing memory but I haven't thought about it yet.
        double[][] votes = new double[X.length][forest.numTrees];
        for (int i=0; i < forest.numTrees; i++) {
            double[] res = Regtree.classify(forest.Trees[i], X);
            for (int j=0; j < res.length; j++) {
                votes[j][i] = res[j];
            }
        }
        
        double[] retn = new double[X.length];
        for (int i=0; i < X.length; i++) {
            double[] best = Utils.mode(votes[i]);
            retn[i] = best[(int)(Math.random()*best.length)];
        }
        return retn;
    }
    
    public static double[][] applyMarginal(RandomForest forest, int[] tree_idxs_used, double[][] Theta) {
        return applyMarginal(forest, tree_idxs_used, Theta, null);
    }
    
    /**
     * Gets a prediction for the given configurations and instances
     * @returns a matrix of size Theta.length*2 where index (i,0) is the prediction for Theta[i] 
     * and (i,1) is the variance of that prediction. See Matlab code for how var is calculated.
     * @see RegtreeFwd.marginalFwd
     */
    public static double[][] applyMarginal(RandomForest forest, int[] tree_idxs_used, double[][] Theta, double[][] X) {
        int nTheta = Theta.length, nTrees = tree_idxs_used.length;
		double[][] retn = new double[nTheta][2]; // mean, var
        
        for (int i=0; i < nTrees; i++) {
            Object[] result = RegtreeFwd.marginalFwd(forest.Trees[tree_idxs_used[i]], Theta, X);
            double[] preds = (double[])result[0];
            double[] vars = (double[])result[1];

            for (int j=0; j < nTheta; j++) {
                double pred = preds[j];
                if (forest.logModel>0) {
                    pred = Math.log10(pred);
                }
                retn[j][0] += pred;
                retn[j][1] += vars[j]+pred*pred;
            }
        }
        
        for (int i=0; i < nTheta; i++) {
            retn[i][0] /= nTrees;
            retn[i][1] /= nTrees;
			retn[i][1] -= retn[i][0]*retn[i][0];
            retn[i][1] = retn[i][1] * ((forest.numTrees+0.0)/Math.max(1, forest.numTrees-1));
        }
        return retn;
    }
    
    /** 
     * Prepares the random forest for marginal predictions.
     * @see RegtreeFwd.preprocess_inst_splits
     */
    public static RandomForest preprocessForest(RandomForest forest, double[][] X) {
        RandomForest prepared = new RandomForest(forest.numTrees, forest.logModel);
        for (int i=0; i < forest.numTrees; i++) {
            prepared.Trees[i] = RegtreeFwd.preprocess_inst_splits(forest.Trees[i], X);
        }
        return prepared;
    }
    
    /** 
     * Prepares the random forest for classification
     * dropysub specifies whether the ysub array should be cleared after preprocessing
     * @see RegtreeFwd.preprocess_for_classification
     */
    public static void preprocessForestForClassification(RandomForest forest, boolean dropysub) {
        for (int i=0; i < forest.numTrees; i++) {
            RegtreeFwd.preprocess_for_classification(forest.Trees[i], dropysub);
        }
    }
    
    public static void serialize(RandomForest forest, String filename) {   
	    try {
			FileOutputStream fileOut = new FileOutputStream(filename);
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(forest);
			out.close();
			fileOut.close();
		} catch (IOException i) {
			i.printStackTrace();
		}
    }
    
    public static RandomForest deserialize(String filename) { 
        RandomForest forest = null;
	    try {
	        FileInputStream fileIn =
	                      new FileInputStream(filename);
	        ObjectInputStream in = new ObjectInputStream(fileIn);
	        try {
				forest = (RandomForest) in.readObject();
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			}
	        in.close();
	        fileIn.close();
	    } catch(IOException i) {
	        i.printStackTrace();
	    }
        return forest;
    }
}
