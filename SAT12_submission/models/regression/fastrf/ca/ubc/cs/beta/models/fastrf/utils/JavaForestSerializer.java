package ca.ubc.cs.beta.models.fastrf.utils;

import ca.ubc.cs.beta.models.fastrf.*;
import java.io.*;

public final class JavaForestSerializer {
    
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
