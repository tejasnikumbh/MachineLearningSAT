package ca.ubc.cs.beta.models.fastrf.utils;

import java.io.Serializable;
import ca.ubc.cs.beta.models.fastrf.*;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;

public final class DebugTreeBuildInputs implements Serializable {
	private static final long serialVersionUID = -3226427544555346375L;
	public double[][] allTheta; 
	public double[][] allX;
	public int[][] dataIdxs;
	public double[] y;
	public RegtreeBuildParams params;
    
    public static void serializeDataForDebug(double[][] allTheta, double[][] allX, int[][] dataIdxs, double[] y, RegtreeBuildParams params, String filename) {   
    	DebugTreeBuildInputs t = new DebugTreeBuildInputs();
    	t.allTheta = allTheta;
    	t.allX = allX;
    	t.dataIdxs = dataIdxs;
    	t.y = y;
    	t.params = params;

	    try {
			FileOutputStream fileOut = new FileOutputStream(filename + ".ser");
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(t);
			out.close();
			fileOut.close();
		} catch (IOException i) {
			i.printStackTrace();
		}
    }
}
