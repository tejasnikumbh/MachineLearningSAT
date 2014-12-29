import java.io.*;
import java.util.*;
import ca.ubc.cs.beta.models.fastrf.*;
import ca.ubc.cs.beta.models.fastrf.utils.*;

public class SingleClassificationTester {
    public static void main(String[] args) throws Exception {
         if (args.length < 2) {
             System.out.println("Error: Please specify both the model file and a csv file containing the features.");
             return;
         }     
         
         System.out.println("Reading in the Forest");
         RandomForest JavaForest = RandomForest.deserialize(args[0]);
 
         System.out.println("Reading in features...");
         Vector<double[]> v = new Vector<double[]>();
         BufferedReader r = new BufferedReader(new FileReader(args[1]));
         String line = r.readLine();
         while(line != null) {
             String[] split = line.split(",");
             double[] data = new double[split.length];
             for (int i=0; i < split.length; i++) {
                 data[i] = Double.parseDouble(split[i]);
             }
             v.add(data);
             line = r.readLine();
         }
         r.close();
         double[][] X = v.toArray(new double[1][0]);

         System.out.println("Transforming features...");
         X = Utils.transform(X, JavaForest.keptColumns, JavaForest.scale, JavaForest.bias);
         
         System.out.println("Classifying...");
         double[] results = RandomForest.classify(JavaForest, X);

         System.out.println("Results:");
         for (int i=0; i < results.length; i++) {
            System.out.println((int)results[i]);
         }
    }
}
