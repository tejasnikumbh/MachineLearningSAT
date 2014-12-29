package ca.ubc.cs.beta.models.fastrf.utils;

import java.util.*;

public class Utils {   
    public static final double l10e = Math.log10(Math.E);

    public static int sum(int[] arr) {
        if (arr == null || arr.length == 0) return 0;
        int l = arr.length;
        int res = 0;
        for (int i=0; i < l; i++) {
            res += arr[i];
        }
        return res;
    }
    
    public static double sum(double[] arr) {
        if (arr == null || arr.length == 0) return 0;
        int l = arr.length;
        double res = 0;
        for (int i=0; i < l; i++) {
            res += arr[i];
        }
        return res;
    }
    
    public static double mean(double[] arr) {
        if (arr == null || arr.length == 0) return 0;
        return sum(arr)/arr.length;
    }
    
    public static double var(double[] arr) {
        if (arr == null || arr.length == 0) return 0;
        int l = arr.length;
        
        double sum = 0, sumSq = 0;
        for (int i=0; i < l; i++) {
            sum += arr[i];
            sumSq += arr[i] * arr[i];
        }
        return (sumSq - sum*sum/l)/Math.max(l-1, 1);
    }
    
    public static double median(double[] arr) {
        if (arr == null || arr.length == 0) return Double.NaN;
        int l = arr.length;
        Arrays.sort(arr);
        return arr[(int)Math.floor(l/2.0)] / 2 + arr[(int)Math.ceil(l/2.0)] / 2;
    }
    
    public static double prod(double[] arr, int start, int end) {
        double result = 1;
        for (int i=start; i < end; i++) {
            result *= arr[i];
        }
        return result;
    }

    public static HashMap<Double, Integer> getFreq(double[] arr) {
        HashMap<Double, Integer> freq = new HashMap<Double, Integer>();
        for (int i=0; i < arr.length; i++) {
            Integer cur = freq.get(arr[i]);
            if (cur == null) cur = 0;
            freq.put(arr[i], cur+1);
        }
        return freq;
    }
    
    public static double[] mode(double[] arr) {
        HashMap<Double, Integer> freq = new HashMap<Double, Integer>();
        int highestFreq = 0;
        ArrayList<Double> best = new ArrayList<Double>();
        for (int i=0; i < arr.length; i++) {
            Integer cur = freq.get(arr[i]);
            if (cur == null) cur = 0;
            cur++;
            if (cur > highestFreq) {
                highestFreq = cur;
                best.clear();
            }
            if (cur == highestFreq) best.add(arr[i]);
            freq.put(arr[i], cur);
        }
        double[] retn = new double[best.size()];
        for (int i=0; i < best.size(); i++) {
            retn[i] = best.get(i);
        }
        return retn;
    }
    
    /**
     * Transforms the given inputs. Assumes keptColumns is 0-indexed, sorted in increasing order,
     * and has the same length as scale and bias.
     */
    public static double[][] transform(double[][] X, int[] keptColumns, double[] scale, double[] bias) {
        if (X == null) throw new RuntimeException("X is null in transform.");
        if (X.length == 0 || keptColumns == null) return X;
        if (scale == null || bias == null || keptColumns.length != scale.length || keptColumns.length != bias.length) {
            throw new RuntimeException("Badly formatted inputs to transform.");
        }
        
        double[][] ret = new double[X.length][keptColumns.length];
        for (int i=0; i < X.length; i++) {
            for (int j=0, k=0; j < X[i].length && k < keptColumns.length; j++) {
                if (j == keptColumns[k]) {
                    ret[i][k] = (X[i][j] - bias[k]) / scale[k];
                    k++;
                }
            }
        }
        return ret;
    }
}
