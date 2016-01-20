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
    
    //Test comment in here
    
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        File maxQuantFile = new File("E:\\Work\\DropBox\\Dropbox\\DocStore\\GrantApps\\2016_ERC\\Data\\QuantCorrelations\\PXD001550\\modificationSpecificPeptides.txt");
        ProteoformerBU pbu = new ProteoformerBU();
        pbu.readPepData(maxQuantFile);
        pbu.computeProtQuantValues();
        //pbu.printProtQuantValue();
        
        //Not currently used
        //pbu.normaliseProtQuantValues();
        
        pbu.createPeptideCorrelations();
        
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
                        pepModString = pepModString.replaceAll(" ", "_").replaceAll("\\(", "").replaceAll("\\)", "");

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
        for(String pep : allPepQuantData.keySet()){
            String protAcc = pepToProtAcc.get(pep);            
            ArrayList<BigInteger> pepValues = allPepQuantData.get(pep);
            double[] doubleValues = new double[pepValues.size()];
            for(int i=0; i<pepValues.size();i++){
                double totalProtAbundance = allProtQuantData.get(protAcc)[i].doubleValue();
                doubleValues[i] = pepValues.get(i).doubleValue()/totalProtAbundance;
                
            }
            pepQuantsAsDoubles.put(pep, doubleValues);
        }
        
        
        for(String protAcc : protToPeps.keySet()){
            System.out.print(protAcc + " ");
            ArrayList<String> peps = protToPeps.get(protAcc);
            for(int i=0;i<peps.size()-1;i++){   //outer loop
                String outerPep = peps.get(i);
                double[] outerPepVals = pepQuantsAsDoubles.get(outerPep);
                 
                for(int j=i+1;j<peps.size();j++){
                    String innerPep = peps.get(j);
                    double[] innerPepVals = pepQuantsAsDoubles.get(innerPep);
                    PearsonsCorrelation pCorr = new PearsonsCorrelation();
                    double corr = pCorr.correlation(outerPepVals, innerPepVals);
                    System.out.print("Outer:" + outerPep + " inner: " + innerPep + " corr: " + corr + "   ");                     
                 }
                System.out.println();
                
            }
            
        }
    }
    
    
    
}
