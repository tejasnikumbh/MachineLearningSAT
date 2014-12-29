package ca.ubc.cs.beta.models.fastrf;

import java.util.*;
import ca.ubc.cs.beta.models.fastrf.utils.Utils;

public class RegtreeFwd {
    /**
     * Propogates data points down the regtree and returns a 1*X.length vector of node #s
     * specifying which node each data point falls into.
     * @param tree the regtree to use
     * @param X a numdatapoints*numvars matrix
     */
    public static int[] fwd(Regtree tree, double[][] X) {  
        int numdata = X.length;
        int numnodes = tree.var.length;
        if (tree.cut.length != numnodes) {
            throw new RuntimeException("cut must be Nx1 vector.");
        }
        if (tree.children.length != numnodes) {
            throw new RuntimeException("children must be Nx2 matrix.");
        }
        
        int[] result = new int[numdata];
        
        for (int i=0; i < numdata; i++) {
            int thisnode = 0;
            while(true) {
                int splitvar = tree.var[thisnode];
                if (splitvar == 0) {
                    // This node not split, store results.
                    result[i] = thisnode;
                    break;
                }
                double cutoff = tree.cut[thisnode];
                int left_kid = tree.children[thisnode][0];
                int right_kid = tree.children[thisnode][1];
                // Determine if the points goes left or goes right
                if(Double.isNaN(X[i][Math.abs(splitvar)-1])){
                    throw new RuntimeException("In fwd, trying to split on variable " + splitvar + " (1-based, negative means categorical), but data point number " + i + " is NaN for that.");
                }
                if (splitvar > 0) { 
                    // continuous variable
                    thisnode = (X[i][splitvar-1] <= cutoff ? left_kid : right_kid);
                } else { 
                    // categorical variable
                    int x = (int)X[i][-splitvar-1];
                    if (x<=0){
                        throw new RuntimeException("Input error in Regtree.fwd: categoricals have to be integers >= 1");
                    }

                    int split = tree.catsplit[(int)cutoff][x-1];
                    if (split == 0) thisnode = left_kid;
                    else if (split == 1) thisnode = right_kid;
                    else throw new RuntimeException("Missing value -- not allowed in this implementation.");
                }
            }
        }
        return result;
    }
    
    /**
     * Propogates configurations(Theta) and instances(X) down the tree, and returns a 1*Theta.length vector of 
     * marginal prediction for each configuration (summed across each of the specified instances).
     * @param tree the regtree to use
     * @param Theta a vector of configuration parameters
     * @param X a vector of instance parameters. If the tree has already been preprocessed, this argument is ignored.
     */
    public static Object[] marginalFwd(Regtree tree, double[][] Theta, double[][] X) {
        if (Theta == null || Theta.length == 0) {
            throw new RuntimeException("Theta must not be empty");
        }
        int thetarows = Theta.length;
        int thetacols = Theta[0].length;
        int numnodes = tree.node.length;
        if (numnodes == 0) {
            throw new RuntimeException("Tree must exist.");
        }
        if (tree.cut.length != numnodes) {
            throw new RuntimeException("cut must be Nx1 vector.");
        }
        if (tree.parent.length != numnodes) {
            throw new RuntimeException("parent must be Nx1 matrix.");
        }
        if (tree.children.length != numnodes) {
            throw new RuntimeException("children must be Nx2 matrix.");
        }
        
        if (!tree.preprocessed) {
            tree = preprocess_inst_splits(tree, X);
        }
        
        double[] result = new double[thetarows];
        double[] vars = new double[thetarows];
        
        LinkedList<Integer> queue = new LinkedList<Integer>();
        
        for (int i=0; i < thetarows; i++) {
            vars[i] = 0;
            queue.add(0);
            while(!queue.isEmpty()) {
                int thisnode = queue.poll();
                while(true) {
                    int splitvar = tree.var[thisnode];
                    double cutoff = tree.cut[thisnode];
                    int left_kid = tree.children[thisnode][0];
                    int right_kid = tree.children[thisnode][1];

                    if (splitvar == 0) {
                        // We are in leaf node. store results.
                        result[i] += tree.weightedpred[thisnode];
                        break;
                    } else if (Math.abs(splitvar) > thetacols) {
                        // Splitting on instance - pass this instance down both children
                        queue.add(right_kid);
                        thisnode = left_kid;
                    } else {
                        if(Double.isNaN(Theta[i][Math.abs(splitvar)-1])){
                            throw new RuntimeException("In marginalFwd, trying to split on variable " + splitvar + " (1-based, negative means categorical), but data point number " + i + " is NaN for that.");
                        }
                        if (splitvar > 0) { // continuous
                            thisnode = (Theta[i][splitvar-1] <= cutoff ? left_kid : right_kid);
                        } else { // categorical
                            int x = (int) Theta[i][-splitvar-1];
                            if (x<=0){
                                throw new RuntimeException("Input error in in Regtree.marginalFwd: categoricals have to be integers >= 1");
                            }
                            int split = tree.catsplit[(int)cutoff][x-1];
                            if (split == 0) thisnode = left_kid;
                            else if (split == 1) thisnode = right_kid;
                            else throw new RuntimeException("Missing value -- not allowed in this implementation.");
                        }
                    }
                }
            }
        }
        Object[] retn = new Object[2];
        retn[0] = result;
        retn[1] = vars;
        return retn;
    }
    
    /**
     * Preprocesses the regtree for marginal predictions using the specified instances. 
     * A call to preprocess_inst_splits(tree,X) followed by marginalFwd(tree,Theta,null) 
     * is equivalent to the call marginalFwd(tree,Theta,X)
     * except that marginalFwd preprocesses the tree internally.
     */
    public static Regtree preprocess_inst_splits(Regtree tree, double[][] X) {
        tree = new Regtree(tree);
        
        int numnodes = tree.node.length;
        if (numnodes == 0) {
            throw new RuntimeException("Tree must exist.");
        }
        if (tree.cut.length != numnodes) {
            throw new RuntimeException("cut must be Nx1 vector.");
        }
        if (tree.nodepred.length != numnodes) {
            throw new RuntimeException("nodepred must be Nx1 vector.");
        }
        if (tree.parent.length != numnodes) {
            throw new RuntimeException("parent must be Nx1 matrix.");
        }
        if (tree.children.length != numnodes) {
            throw new RuntimeException("children must be Nx2 matrix.");
        }
        
        tree.weights = new double[numnodes];
        tree.weightedpred = new double[numnodes];
        tree.weightedvar = new double[numnodes];
        for (int i=0; i < numnodes; i++) {
            tree.weights[i] = 0;
			if (tree.var[i] != 0) continue;
            
            tree.weightedpred[i] = tree.nodepred[i];
            tree.weightedvar[i] = tree.nodevar[i];

        }
        
        if (X == null) {
            tree.preprocessed = true;
            return tree;
        }
       
        int numinsts = X.length;       
        int thetacols = tree.npred - X[0].length;

        LinkedList<Integer> queue = new LinkedList<Integer>();

        for (int i=0; i < numinsts; i++) {
            queue.add(0);
            while(!queue.isEmpty()) {
                int thisnode = queue.poll();
                while(true) {
                    int splitvar = tree.var[thisnode];
                    double cutoff = tree.cut[thisnode];
                    int left_kid = tree.children[thisnode][0];
                    int right_kid = tree.children[thisnode][1];

                    if (splitvar == 0) {
                        // We are in leaf node. 
                        tree.weights[thisnode]++;
                        break;
                    } else if (Math.abs(splitvar) <= thetacols) {
                        // Splitting on Theta - pass this instance down both children
                        queue.add(right_kid);
                        thisnode = left_kid;
                    } else {
                        if(Double.isNaN(X[i][Math.abs(splitvar)-1-thetacols])){
                            throw new RuntimeException("In preprocess_inst_splits, trying to split on variable " + splitvar + " (1-based, negative means categorical), but data point number " + i + " is NaN for that.");
                        }
                        if (splitvar > 0) { // continuous
                            thisnode = (X[i][splitvar-1-thetacols] <= cutoff ? left_kid : right_kid);
                        } else { // categorical
                            int x = (int)X[i][-splitvar-1-thetacols];
                            int split = tree.catsplit[(int)cutoff][x-1];
                            if (split == 0) thisnode = left_kid;
                            else if (split == 1) thisnode = right_kid;
                            else throw new RuntimeException("Missing value -- not allowed in this implementation.");
                        }
                    }
                }
            }
        }

		for (int i=0; i < numnodes; i++) {
            tree.weights[i] /= numinsts;
            tree.weightedpred[i] *= tree.weights[i];
            tree.weightedvar[i] *= tree.weights[i] * tree.weights[i];
        }

        // cut leaf splits on instances. 
        cut_instance_leaf_split_helper(tree, thetacols, 0);
        
        tree.preprocessed = true;
        return tree;
    }
    
    // Returns whether thisnode is a leaf at the end of the function
    private static int cut_instance_leaf_split_helper(Regtree tree, int thetacols, int thisnode) {
		if (tree.var[thisnode] == 0) return 1;

        int left_kid = tree.children[thisnode][0], right_kid = tree.children[thisnode][1];
        int ret = 0;
        
		if (
               cut_instance_leaf_split_helper(tree, thetacols, left_kid)
                   + 
               cut_instance_leaf_split_helper(tree, thetacols, right_kid) == 2 // This is so we don't short-circuit
               && Math.abs(tree.var[thisnode]) > thetacols
            ) {
                // both children are leaves, and this is a split on an instance
                make_into_leaf(tree, thisnode);
                ret = 1;
        }
        
        tree.weights[thisnode] = tree.weights[left_kid] + tree.weights[right_kid];
        if (ret == 0 && tree.weights[thisnode] == 0) {
            make_into_leaf(tree, thisnode);
            ret = 1;
        }
        return ret;
    }
    
    private static void make_into_leaf(Regtree tree, int thisnode) {
        int left_kid = tree.children[thisnode][0], right_kid = tree.children[thisnode][1];
        tree.children[thisnode][0] = 0;
        tree.children[thisnode][1] = 0;
		tree.var[thisnode] = 0;

		tree.weightedpred[thisnode] = tree.weightedpred[left_kid] + tree.weightedpred[right_kid];
        tree.weightedvar[thisnode] = tree.weightedvar[left_kid] + tree.weightedvar[right_kid];

		/* for debugging only
        if (tree.resultsStoredInLeaves) {
            tree.ysub[thisnode] = new double[tree.nodesize[thisnode]];
            System.arraycopy(tree.ysub[left_kid], 0, tree.ysub[thisnode], 0, tree.nodesize[left_kid]);
            System.arraycopy(tree.ysub[right_kid], 0, tree.ysub[thisnode], tree.nodesize[left_kid], tree.nodesize[right_kid]);
        } else {
            tree.ysub[thisnode][0] = tree.ysub[left_kid][0] + tree.ysub[right_kid][0]; // sum
            tree.ysub[thisnode][1] = tree.ysub[left_kid][1] + tree.ysub[right_kid][1]; // sumsq
        }
        tree.recalculateStats(thisnode);
		*/
    }
    
    /**
     * Preprocesses the regtree for classification. Stores an array of the most popular responses for each leaf.
     * dropysub specifies whether the ysub array should be cleared after preprocessing
     */
    public static void preprocess_for_classification(Regtree tree, boolean dropysub) {
        if (!tree.resultsStoredInLeaves) {
            throw new RuntimeException("Classification can only be done if the tree was built with the resultsStoredInLeaves flag on.");
        }
        tree.bestClasses = new double[tree.numNodes][];
        for (int i=0; i < tree.numNodes; i++) {
            if (tree.var[i] != 0) // not a leaf
                continue;
            tree.bestClasses[i] = Utils.mode(tree.ysub[i]);
        }
        tree.preprocessed_for_classification = true;
        if (dropysub) {
            tree.ysub = null;
            tree.resultsStoredInLeaves = false;
        }
    }
}
