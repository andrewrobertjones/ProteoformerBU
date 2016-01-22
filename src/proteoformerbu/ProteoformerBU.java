/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package proteoformerbu;

import java.io.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.rosuda.JRI.*;

/**
 *
 * @author jonesar
 */
public class ProteoformerBU {

    private String[] fileColumns;
    private ArrayList<Integer> intensityColumns = null;
    private String sep = "\t";                          //separator
    private HashMap<String,ArrayList<BigInteger>> allPepQuantData = new HashMap();    //map from pepID to allQuant value
    private HashMap<String,BigInteger[]> allProtQuantData = new HashMap();    //map from prot accession to total abundance for each column
    private HashMap<String,BigInteger[]> allNormProtQuantData = new HashMap();  //map from prot accession to normalised total abundance for each column
    private HashMap<String,String> pepToProtAcc = new HashMap();
    private HashMap<String,ArrayList<String>> protToPeps = new HashMap();
    private HashMap<String,String[][]> protAccToDistanceMatrix = new HashMap();
    private String outputFolder = "E:/Work/ProteomicsSoftware/PBU/outputs/";
    private ArrayList<String> matrixFilesToProcess = new ArrayList();
    private ArrayList<String> matrixRawFilesToProcess = new ArrayList();
    private int MAX_PROT_ACC_LENGTH = 60;   //cut-off for long lists of protein groups
    //Test comment in here
    
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        File maxQuantFile = new File("E:/Work/DropBox/Dropbox/DocStore/GrantApps/2016_ERC/Data/QuantCorrelations/PXD001550/modificationSpecificPeptides.txt");
        ProteoformerBU pbu = new ProteoformerBU();
        pbu.readPepData(maxQuantFile);
        pbu.computeProtQuantValues();
        //pbu.printProtQuantValue();
        
        //Not currently used
        //pbu.normaliseProtQuantValues();
        
        pbu.createPeptideCorrelations();
        pbu.processAllMatricesForR();
        
    }
    
    /*
    * Method to read data from MaxQuant modificationSpecificPeptides.txt to retrieve quant values
    */
    private void readPepData(File mqFile){
        
        int rowCounter=0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(mqFile));
            String line;
            while ((line = br.readLine()) != null) {
                
               String[] cells = line.split(sep);
               if(rowCounter==0){
                   fileColumns = cells;              
                   intensityColumns = getIntensityColumns(); //Get columns that contain quant data                   
               }
               else{
                    String protAcc = cells[5];
                   
                    if(protAcc!=null && !protAcc.equals("")){ //Ignore any peps without an assigned protein
                   
                        String pepModString = cells[0] + "__mods:" + cells[1];
                        pepModString = pepModString.replaceAll(" ", "_").replaceAll("\\(", "").replaceAll("\\)", "").replaceAll(",",";");

                        ArrayList<BigInteger> quantValues = new ArrayList();
                        for(Integer col : intensityColumns){
                            quantValues.add(new BigInteger(cells[col.intValue()]));        //Grab values from cells and put into quantValues array
                        }
                        pepToProtAcc.put(pepModString, protAcc);     //add map from pepID to protID

                        ArrayList<String> peps = null;
                        if(protToPeps.containsKey(protAcc)){
                            peps = protToPeps.get(protAcc);
                        }
                        else{
                            peps = new ArrayList();
                        }                   
                        peps.add(pepModString);
                        protToPeps.put(protAcc, peps);
                        allPepQuantData.put(pepModString, quantValues);
                        //System.out.println(protAcc+ "---> " + pepModString + "->" + quantValues.toString());
                   }
               }
               rowCounter++;
            }
        }
        catch(IOException e){
            e.printStackTrace();
        }        
    }
    
    private ArrayList<Integer> getIntensityColumns(){
        ArrayList<Integer> intColumns = new ArrayList();        
        for(int i =0; i<fileColumns.length;i++){
            String header = fileColumns[i];
            if(!header.equals("Intensity")){
                if(header.contains("Intensity")){
                    //System.out.print("header: " + header + " ");
                    intColumns.add(new Integer(i));
                }
            }
        }        
        return intColumns;
    }
    
    private void computeProtQuantValues(){
        
        int totalColumns = intensityColumns.size();
        for(String protAcc : protToPeps.keySet()){
            ArrayList<String> peps = protToPeps.get(protAcc);
            BigInteger[] protQuantData = new BigInteger[totalColumns]; 
            for(String pep : peps){ //loop on rows (peps)
               ArrayList<BigInteger> pepIntValues = allPepQuantData.get(pep);
               
               for(int i=0;i<totalColumns;i++){   //loop on columns (different condition)
                   BigInteger protIntensity = null;
                   if(protQuantData[i]==null){
                       protIntensity=new BigInteger("0");                       
                   }
                   else{
                       protIntensity = protQuantData[i];
                   }
                   protIntensity= protIntensity.add(pepIntValues.get(i));
                   protQuantData[i] = protIntensity;
               }
            }
            allProtQuantData.put(protAcc,protQuantData);
        }        
    }
    
    /*
    Algorithm:
    Sum protein quant in each column to get total protein abundance - find median value
    For each value in array (except median), then mutliply by 1/total protein abundance
    */
    private void normaliseProtQuantValues(){
        
        BigInteger[] totalProteinAbundances = new BigInteger[intensityColumns.size()];
                
        //Initialize all values to zero
        for(int i=0; i<totalProteinAbundances.length;i++){
            totalProteinAbundances[i]=new BigInteger("0");
        }
        
        for(String protAcc : allProtQuantData.keySet()){
            BigInteger[] protQuantValues = allProtQuantData.get(protAcc);
            
            for(int i=0;i<protQuantValues.length;i++){  
                totalProteinAbundances[i] = totalProteinAbundances[i].add(protQuantValues[i]);
            }
        }
        
        /*
        for(BigInteger total : totalProteinAbundances){
            System.out.print(total + "\t");
       }
        System.out.println("");
        
        */
        
        int medianPos = medianPositionInArray(totalProteinAbundances);
        //System.out.println("middlepos:" + medianPos + " value:" + totalProteinAbundances[medianPos].intValue());
        
        double[] normalisationFactors = new double [intensityColumns.size()];
        for(int i=0;i<normalisationFactors.length;i++){
            BigInteger medianAbundance = totalProteinAbundances[medianPos];
            normalisationFactors[i]=totalProteinAbundances[i].doubleValue()/totalProteinAbundances[medianPos].doubleValue();            
        }
        //System.out.println(Arrays.toString(normalisationFactors));
        

        //Second loop to create normalised values
        for(String protAcc : allProtQuantData.keySet()){
            BigInteger[] protQuantValues = allProtQuantData.get(protAcc);
            BigInteger[] normProtQuantValues = new BigInteger[protQuantValues.length];
            for(int i=0;i<protQuantValues.length;i++){                  
                BigInteger normalisedProtQuant = new BigDecimal(protQuantValues[i].doubleValue() * normalisationFactors[i]).toBigInteger();
                normProtQuantValues[i]=normalisedProtQuant;
            }
            allNormProtQuantData.put(protAcc,normProtQuantValues);
        }      
    }
    
    
    // quick and dirty class to get approx middle value in aray
    public static int medianPositionInArray(BigInteger[] values) {
        double[] tempArray = new double[values.length];
         
        int i=0;
        for(BigInteger value : values){
            tempArray[i] = value.doubleValue();
        }
        
        Arrays.sort(tempArray);
        return tempArray.length/2;
    }
    
    private void printProtQuantValue(){
        
        
        for(String protAcc : allProtQuantData.keySet()){
            System.out.print(protAcc + "\t");
            for(BigInteger protQuant : allProtQuantData.get(protAcc)){
                System.out.print(protQuant + "\t");
            }
            System.out.print("\n");
        }
        
    }
    
    /*
    * Within each protein:
    *
    * Perform n * n-1 correlations on all peptides within the set
    * Put into a matrix structure
    */
    private void createPeptideCorrelations(){
        
        /* This code does regular doubles, alternative code computes pepQuants as proportion of protein abundance
        HashMap<String,double[]> pepQuantsAsDoubles = new HashMap();
        for(String pep : allPepQuantData.keySet()){
            ArrayList<BigInteger> pepValues = allPepQuantData.get(pep);
            double[] doubleValues = new double[pepValues.size()];
            for(int i=0; i<pepValues.size();i++){
                doubleValues[i] = pepValues.get(i).doubleValue();
            }
            pepQuantsAsDoubles.put(pep, doubleValues);
        }
        */
        
        HashMap<String,double[]> pepQuantsAsDoubles = new HashMap();
        HashMap<String,double[]> pepRawQuantsAsDoubles = new HashMap();
        for(String pep : allPepQuantData.keySet()){
            String protAcc = pepToProtAcc.get(pep);            
            ArrayList<BigInteger> pepValues = allPepQuantData.get(pep);
            double[] doubleValues = new double[pepValues.size()];
            double[] doubleRawValues = new double[pepValues.size()];
            for(int i=0; i<pepValues.size();i++){
                double totalProtAbundance = allProtQuantData.get(protAcc)[i].doubleValue();
                doubleValues[i] = pepValues.get(i).doubleValue()/totalProtAbundance;
                doubleRawValues[i] = pepValues.get(i).doubleValue();
                
            }
            pepQuantsAsDoubles.put(pep, doubleValues);
            pepRawQuantsAsDoubles.put(pep, doubleRawValues);
        }
        
        
        for(String protAcc : protToPeps.keySet()){
            //System.out.print(protAcc + " ");
            ArrayList<String> peps = protToPeps.get(protAcc);
            
            String headerRow = "";
            String[][] matrix = new String[peps.size()][peps.size()];
            String[][] matrixRaw = new String[peps.size()][peps.size()];
            
            /*
            Prepare String matrix for R as follows:
            ASGQAFELILSPR_Phospho_STY	ASGQAFELILSPR_2Phospho_STY	RASGQAFELILSPR_2Phospho_STY	RASGQAFELILSPR_Phospho_STY	RKSHEAEVLK_Phospho_STY	SKESVPEFPLSPPK_Phospho_STY	DLSLEEIQK_Unmodified	DLSLEEIQK_Phospho_STY
            0	1.369	0.969	0.771	0.881	1.365	0.921	0.916
            1.369	0	1.416	1.199	0.958	1.491	0.936	0.577
            0.969	1.416	0	1.064	0.918	1.412	0.612	0.934
            0.771	1.199	1.064	0	1.092	0.917	1.043	1.264
            0.881	0.958	0.918	1.092	0	1.426	0.806	0.952
            1.365	1.491	1.412	0.917	1.426	0	1.446	1.101
            0.921	0.936	0.612	1.043	0.806	1.446	0	1.465
            0.916	0.577	0.934	1.264	0.952	1.101	1.465	0
            */
            HashMap<String, double[]> allPepQuants = new HashMap();
            HashMap<String, double[]> allRawPepQuants = new HashMap();
            
            for(int i=0;i<peps.size();i++){   //outer loop                
                String outerPep = peps.get(i);
                headerRow+=outerPep+",";
                double[] outerPepVals = pepQuantsAsDoubles.get(outerPep);
                double[] outerPepRawVals = pepRawQuantsAsDoubles.get(outerPep);                
                
                allPepQuants.put(outerPep,outerPepVals);
                allRawPepQuants.put(outerPep,outerPepRawVals);
                 
                for(int j=0;j<peps.size();j++){
                    String innerPep = peps.get(j);
                    double[] innerPepVals = pepQuantsAsDoubles.get(innerPep);
                    double[] innerPepRawVals = pepRawQuantsAsDoubles.get(innerPep);
                    PearsonsCorrelation pCorr = new PearsonsCorrelation();
                    double corr = pCorr.correlation(outerPepVals, innerPepVals);
                    double distance = 1-corr;
                    matrix[i][j] = ""+distance;  
                    
                    //Now run again for raw values
                    corr = pCorr.correlation(outerPepRawVals, innerPepRawVals);
                    distance = 1-corr;
                    matrixRaw[i][j] = ""+distance;
                 }
            }
            
            if(protAcc.length()<=MAX_PROT_ACC_LENGTH && matrix.length >=3){
                //protAccToDistanceMatrix.put(protAcc, matrix); //not currently used
                headerRow = headerRow.substring(0,headerRow.length()-1); //Remove final comma
                String outFile = printMatrix(protAcc,headerRow,matrix, false);
                printPepQuantValues(protAcc,allPepQuants,false);
                printPepQuantValues(protAcc,allRawPepQuants,true);

                if(!outFile.equals("ERROR")){   //Internal error code if file could not be written
                    matrixFilesToProcess.add(outFile);
                }
                
                headerRow = headerRow.substring(0,headerRow.length()-1); //Remove final comma
                outFile = printMatrix(protAcc,headerRow,matrixRaw,true);
                if(!outFile.equals("ERROR")){   //Internal error code if file could not be written
                    matrixRawFilesToProcess.add(outFile);
                }
            }
            //System.out.println(protAcc + "--> " + matrix.toString());
        }
    }
    
    private void printPepQuantValues(String protAcc, HashMap<String,double[]> pepQuants, Boolean isRaw){
        String outFile = "";
        
        if(isRaw){
            outFile = protAcc + "_quantMatrix.csv";
        }
        else{
            outFile = protAcc + "_quantMatrixRaw.csv";
        }
        
        try{
            Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFolder + outFile), "utf-8")); 
            
            String matrix="";
            for (String pep : pepQuants.keySet()){
                double[] quants = pepQuants.get(pep);
                matrix+=pep+",";
                for (int j = 0; j<quants.length; j++){
                   matrix+= "" + quants[j]+",";
                }
                matrix = matrix.substring(0, matrix.length()-1);//Remove extra comma at end
                matrix += "\n";
           }
            writer.write(matrix);
            writer.close();
        }
        catch(IOException e){
            e.printStackTrace();
            outFile = "ERROR";
        }
    }
    
    
    /*
    * Prints out 2D matrix of distance values to csv for processing with R
    */
    private String printMatrix(String protAcc, String header, String[][] matrix, Boolean isRaw){        

        protAcc = protAcc.replaceAll(";", "_").replaceAll(":", "_").replaceAll(",", "_");
        String outFile = "";
        
        if(isRaw){
            outFile = protAcc + "_CorrMatrix.csv";
        }
        else{
            outFile = protAcc + "_CorrMatrixRaw.csv";
        }
        
        
        try{
            Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFolder + outFile), "utf-8")); 
            writer.write(header+"\n");
            String matrixAsString = "";
            for (int i = 0; i<matrix.length; i++){
                for (int j = 0; j<matrix[i].length; j++){
                   matrixAsString += matrix[i][j] + ",";
                } 
                matrixAsString = matrixAsString.substring(0, matrixAsString.length()-1);//Remove extra comma at end
                matrixAsString += "\n";
           }
            writer.write(matrixAsString+"\n");
            writer.close();
        }
        catch(IOException e){
            e.printStackTrace();
            outFile = "ERROR";
        }
        return outFile;        
    }
    
    
    /*
    * Helper method to loop through all files and call appropriate R routines
    */
    private void processAllMatricesForR(){
        
        String batch = "";
        for(String matrixFile : matrixFilesToProcess){
            String rScriptFile = createDistanceTreeInR(matrixFile);
            batch += "RScript.exe "+rScriptFile+ "\n";
        }   
        
        //Run again for matrices based on raw data i.e. without normlisation based on summed (protein) peptide abundance
        for(String matrixFile : matrixRawFilesToProcess){
            String rScriptFile = createDistanceTreeInR(matrixFile);
            batch += "RScript.exe "+rScriptFile+ "\n";           
        }  
        
        String rBatchFile = "RunAllRComms.bat";
        
        try{
            Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFolder + rBatchFile), "utf-8")); 
            writer.write(batch);           
            writer.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
    }
    
    
    private String createDistanceTreeInR(String matrixFile){
        /*
        
        jpeg(filename = "P55201_CorrMatrix.jpeg",width = 480, height = 480, units = "px", pointsize = 12,quality = 75)
plot(treee, "u")
title("Distance tree of P55201 protein group")
dev.off()
        */
        String treePDFFile = matrixFile.substring(0,matrixFile.lastIndexOf("csv"))+"pdf";
        String treeJPG = matrixFile.substring(0,matrixFile.lastIndexOf("csv"))+"jpg";
        String rScriptFile = matrixFile.substring(0,matrixFile.lastIndexOf("csv"))+"r";
        String protAcc = matrixFile.substring(0,matrixFile.lastIndexOf("_CorrMatrix"));
               
        String rCode = "library(ape)\n"
                + "setwd(\""+outputFolder+"\")\n"
                + "pep_data = read.csv(\""+matrixFile+"\", header = TRUE)\n"
                + "mat <- data.matrix(pep_data)\n"
                + "rownames(mat) <- colnames(mat)\n"
                + "treee <- nj(mat)\n"
                + "pdf(\""+treePDFFile+"\",width=14,height=10)\n"
                + "par(cex=0.8)\n"
                + "plot(treee, \"u\")\n"
                + "title(\"Distance tree of " + protAcc + " protein group\")\n"                
                + "dev.off()\n"
                + "jpeg(filename = \""+treeJPG +"\",width = 480, height = 480, units = \"px\", pointsize = 12,quality = 75)"
                + "plot(treee, \"u\")\n"
                + "title(\"Distance tree of " + protAcc + " protein group\")\n"                
                + "dev.off()\n"
                ;
        
        
        //Write a temp file for R to process
    
        try{
            Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFolder + rScriptFile), "utf-8")); 
            writer.write(rCode+"\n");           
            writer.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        
        
        //Rengine engine = createRengine();
        //engine.eval(rCode);
        //Rengine re = new Rengine(new String[]{"--vanilla"}, false, null);
        /*
        try{
            Runtime rt = Runtime.getRuntime();
           // String command = "RScript.exe "+outputFolder + "tempTree.R";
            String command = "R CMD BATCH "+outputFolder + "tempTree.R";
            Process pr = rt.exec(command);
            //System.out.println("Comm:" + command);
            InputStream error = pr.getErrorStream();
            
            pr.waitFor();        
            for (int i = 0; i < error.available(); i++) {
                System.out.println("Error1:" + error.read());                
            }         
            pr.destroyForcibly();
            
        }
        catch(IOException e){           
            e.printStackTrace();
        }
        catch(InterruptedException ie){
           ie.printStackTrace();
        }
        */
        return rScriptFile;
    }
    
    /**
    * Method to generate a R thread inside the java application.
    * 
    * @return Generated R engine.
    */
   private Rengine createRengine() {
		// ensure, that the right versions of R and Java are available
		if (!Rengine.versionCheck()) {
				System.err
				    .println("** JRI R-Engine: Version mismatch - Java files don't match library version.");
				System.exit(1);
		}
		System.out.println("\n------------------------------");
		System.out.println("Creating JRI R-Engine");
		// arguments which should be passed to R
		String[] args = new String[3];
		args[0] = "--quiet"; // Don't print startup message
		args[1] = "--no-restore"; // Don't restore anything
		args[2] = "--no-save";// Don't save workspace at the end of the session
		// generate new R engine
		Rengine re = new Rengine(args, false, null);
		System.out.println("JRI R-Engine created, waiting for R...");
		// wait until thread to create R is ready
		if (!re.waitForR()) {
				System.out.println("Cannot load R");
				return null;
		}
		// print R engine arguments
		System.out.print("JRI R-Engine call: ");
		for (int i = 0; i < args.length; i++) {
				System.out.print(args[i] + " ");
		}
		System.out.println("...done!");
		System.out.println("------------------------------\n");
		// return the R engine
		return re;
    }
}





