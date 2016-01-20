/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package proteoformerbu;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.math.BigInteger;
import java.util.Arrays;

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
        pbu.normaliseProtQuantValues();
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
                totalProteinAbundances[i].add(protQuantValues[i]);
            }
            
        }
        
        for(BigInteger total : totalProteinAbundances){
            System.out.print(total + "\t");
       }
        System.out.println("");
        
        int medianPos = medianPositionInArray(totalProteinAbundances);
        System.out.println("middlepos:" + medianPos + " value:" + totalProteinAbundances[medianPos].intValue());
        
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
    
    
    
}
