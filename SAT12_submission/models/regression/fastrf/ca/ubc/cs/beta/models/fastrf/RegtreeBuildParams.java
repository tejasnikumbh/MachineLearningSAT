package ca.ubc.cs.beta.models.fastrf;

import java.util.Random;

/**
 * @param catDomainSizes a vector of size X[0].length indicating the size of the domain for the corresponding categorical feature, or 0 if the feature is continuous.
 * @param condParents the conditional parents of the given variable (column index into X). Empty if no parents. NOT SORTED
 * @param condParentVals the OK values of all the conditional parents of the given variable (each element is a matrix). Indexed wrt condParents.
 * @param splitMin the minimum number of data points in each node.
 * @param ratioFeatures the percentage (0, 1] of the features to use in deciding the split at each node.
 * @param logModel whether to build a model in log space
 * @param storeResponses whether to store the responses themselves in the leaves or just to store some statistic (sum, sumofsquares, and leaf size)
 * @param seed -1 means don't use a seed (i.e. create a new Random but don't call setSeed). 
 * @param random will be used instead of seed if it's not null.
 */
public class RegtreeBuildParams implements java.io.Serializable {    
	private static final long serialVersionUID = -6803645785543626390L;
	public int[] catDomainSizes;
    public int[][] condParents = null;
    public int[][][] condParentVals = null;
    public int[] topologicalParameterOrder = null;
    
    public int splitMin;
    public double ratioFeatures;
    public int logModel;
    public boolean storeResponses;
    
    public long seed = -1;
    public Random random = null;
    
    /**
     * Matlab cell arrays don't carry over to Java very well so create the conditional arrays from the cell arrays
     * Called from matlab instead of setting the condParents and condParentVals arrays
     */
    public void conditionalsFromMatlab(int[] cond, int[] condParent, Object[] condParentValsObj, int nvars) {
        condParents = new int[nvars][];
        condParentVals = new int[nvars][][];
        
        for (int i=0; i < nvars; i++) {
            int count = 0;
            for (int j=0; j < cond.length; j++) {
                if (cond[j]-1 == i) {
                    count++;
                }
            }
            condParents[i] = new int[count];
            condParentVals[i] = new int[count][];
            
            count = 0;
            for (int j=0; j < cond.length; j++) {
                if (cond[j]-1 == i) {
                    condParents[i][count] = condParent[j]-1;
                    if(condParentValsObj[j] instanceof int[]) {
                        condParentVals[i][count] = (int[])condParentValsObj[j];
                    } else {
                        condParentVals[i][count] = new int[1];
                        condParentVals[i][count][0] = (Integer)condParentValsObj[j];
                    }
                    count++;
                }
            }
        }
    }
}