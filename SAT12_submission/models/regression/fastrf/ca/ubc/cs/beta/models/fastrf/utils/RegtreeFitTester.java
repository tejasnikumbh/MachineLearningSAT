package ca.ubc.cs.beta.models.fastrf.utils;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import ca.ubc.cs.beta.models.fastrf.*;

public class RegtreeFitTester {
	public static void main(String [] args) {
		 DebugTreeBuildInputs t = null;
	     try {
	        FileInputStream fileIn =
	                      new FileInputStream("ca/ubc/cs/beta/models/fastrf/cplex12-CORLAT-Matrix-30percent.ser");
	        ObjectInputStream in = new ObjectInputStream(fileIn);
	        try {
				t = (DebugTreeBuildInputs) in.readObject();
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			}
	        in.close();
	        fileIn.close();
	    }catch(IOException i) {
	        i.printStackTrace();
	        return;
	    }
	    
	    RegtreeFit.fit(t.allTheta, t.allX, t.dataIdxs, t.y, t.params);
	}
}
