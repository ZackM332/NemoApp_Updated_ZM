package edu.uwb.nemolib;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Subsystem to store relative frequencies for each subgraph pattern / label and
 * calculate the p-values and z-scores for those labels. Performs calculations
 * upon initialization.
 */
public final class RelativeFrequencyAnalyzer {
	
	// labels mapped to a list containing relativeFrequencies
	private Map<String, Double> targetLabelToRelativeFrequencies;
	private Map<String, Double> randomLabelToMeanRelativeFrequencies;
	private Map<String, Double> directLabelToMeanRelativeFrequencies;
	private Map<String, Double> zScores;
	private Map<String, Double> pValues;
        private Map<String, Double> zDScores; // z^ scores from the direct method
        private Map<String, List<String>> rep; // subgraph representation for each graph
        private int motifSize;
        // debug
        //private StringBuilder debug;
        
	// do not allow default constructor
	private RelativeFrequencyAnalyzer() {
		throw new AssertionError();
	}
	
	/**
	 * Constructor for a Statistical Analysis object.
	 *
	 * @param randGraphRelFreqs   Labels paired with lists of relative
	 *                            frequencies of subgraph patterns found in
	 *                            a random graph pool
	 * @param targetGraphRelFreqs Labels paired with relative frequencies as
	 *                            found in the target network.
         * 
         * @param direct              Labels paired with relative frequencies
         *                            from the direct method
	 */
	public RelativeFrequencyAnalyzer(Map<String, List<Double>> randGraphRelFreqs,
	                                 Map<String, Double> targetGraphRelFreqs, 
                                         Map<String, List<Double>> direct,
                                         int motifSize) {
                
		this.targetLabelToRelativeFrequencies = targetGraphRelFreqs;
		
		this.zScores = new HashMap<>();
		calculateZScores(randGraphRelFreqs, targetGraphRelFreqs);
		this.pValues = new HashMap<>();
		calculatePValues(randGraphRelFreqs, targetGraphRelFreqs);
                
                this.zDScores = new HashMap<>();
                calculateZDScores(direct, targetGraphRelFreqs);
                
		this.randomLabelToMeanRelativeFrequencies = calcRandMeans(randGraphRelFreqs);
		this.directLabelToMeanRelativeFrequencies = calcRandMeans(direct);
                
                HashSet<String> alllabels = new HashSet<String>();
                alllabels.addAll(randGraphRelFreqs.keySet());
                alllabels.addAll(targetGraphRelFreqs.keySet());
                alllabels.addAll(direct.keySet());
                this.motifSize = motifSize;
                this.rep = getRepresentation(alllabels);
	}
	
        private Map<String, List<String>> getRepresentation(HashSet<String> alllabels)
        {
            Map<String, List<String>> res = new HashMap<>();
            
            List<String> g6s = new ArrayList<>();
            
            // preparing input file:
            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter("subgraphs2.txt"));
                for (String key : alllabels) {
                    writer.write(key);
                    g6s.add(key);
                    writer.write('\n');
                }
                writer.close();
            }catch (IOException e){
                System.err.println(e.getMessage());
                //return "Error: subgraphs.txt write failed";
                return res;
            }


            // processing input file
            Runtime rt = Runtime.getRuntime();
            Process proc;
            try {
                proc = rt.exec("./showg -a subgraphs2.txt");
                
            }catch (IOException e){
                System.err.println(e.getMessage());
                //return "Error: showg.exe did not work or path to showg.exe is wrong";
                return res;
            }
            
            try {
                proc.waitFor();
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
            }

            BufferedReader stdInput = new BufferedReader(new
                    InputStreamReader(proc.getInputStream()));

            // Read the output from the command
            List<String> adjMatrixs = new ArrayList<>();
            String s = null;
            try {
                while ((s = stdInput.readLine()) != null) {
                    if(!s.equals("")){
                        adjMatrixs.add(s);
                    }
                }
            }catch (IOException e){
                System.err.println(e.getMessage());
                //return "Error: read failed";
                return res;
            }
            
            int k = 0;
            for(String g6: g6s){
                List<String> adj = new ArrayList<>();
                for(int i = 0; i <= this.motifSize; i++){
                    if(adjMatrixs.get(k * (this.motifSize + 1) + i).startsWith("Graph")){
                        continue;
                    }
                    adj.add(adjMatrixs.get(k * (this.motifSize + 1) + i));
                }
                res.put(g6, adj);
                k++;
            }
            
            return res;
        }
        
	public Map<String, Double> getRandMeans() {
		return randomLabelToMeanRelativeFrequencies;
	}
	
	
	private Map<String, Double> calcRandMeans(Map<String, List<Double>> randGraphRelFreqs) {
		Map<String, Double> result = new HashMap<>();
		for (Map.Entry<String, List<Double>> entry : randGraphRelFreqs.entrySet()) {
			double mean = 0.0;
			for (double freq : entry.getValue()) {
				mean += freq;
			}
			mean /= entry.getValue().size();
			result.put(entry.getKey(), mean);
		}
		return result;
	}
	
	/**
	 * Get the z-scores for this RelativeFrequencyAnalyzer object.
	 *
	 * @return a map containing labels and corresponding z-scores
	 */
	public Map<String, Double> getZScores() {
		return zScores;
	}
	
	/**
	 * Get the z-score for a specified label for the data sets contained in
	 * this RelativeFrequencyAnalyzer object
	 *
	 * @param label the label for which to get the z-score
	 * @return the z-score for the given label
	 */
	private double getZScore(String label) {
		return zScores.getOrDefault(label, 0.0);
	}
	
        /**
	 * Get the z-score for a specified label for the data sets contained in
	 * this RelativeFrequencyAnalyzer object
	 *
	 * @param label the label for which to get the z-score
	 * @return the z-score for the given label
	 */
	private double getZDScore(String label) {
		return zDScores.getOrDefault(label, 0.0);
	}
        
        
	/**
	 * Get the p-values for this RelativeFrequencyAnalyzer object.
	 *
	 * @return a map containing labels and corresponding p-values
	 */
	public Map<String, Double> getPValues() {
		return pValues;
	}
	
	/**
	 * Get the p-value for a specified label for the data sets contained in
	 * this RelativeFrequencyAnalyzer object
	 *
	 * @param label the label for which to get the z-score
	 * @return the z-score for the given label
	 */
	private double getPValue(String label) {
		return pValues.getOrDefault(label, 0.0);
	}
	
	// TODO Don't recalculate randMean
	// calculates z-scores for this Statistical Analysis object
	private void calculateZScores(Map<String, List<Double>> randGraphRelFreqs,
	                              Map<String, Double> targetGraphRelFreqs) {
		
		//System.out.println("randGraphRelFreqs..." + randGraphRelFreqs);
		
		HashSet<String> alllabels = new HashSet<String>();
		
		alllabels.addAll(targetGraphRelFreqs.keySet());
		alllabels.addAll(randGraphRelFreqs.keySet());
		for (String label : alllabels) {
			double zScore = 0.0;
			if (!randGraphRelFreqs.containsKey(label))
				zScore = Double.POSITIVE_INFINITY;
			else {
				
				List<Double> freqs = randGraphRelFreqs.get(label);
				double targetGraphFreq = targetGraphRelFreqs.getOrDefault(label, 0.0);
				
				zScore = calcZscore(freqs, targetGraphFreq);
			}
			zScores.put(label, zScore);
			
		}
  
		/*for (Map.Entry<String, List<Double>> labelFreqPair :
				randGraphRelFreqs.entrySet()) {
			String label = labelFreqPair.getKey();
			List<Double> freqs = labelFreqPair.getValue();
			double randMean = calcMean(freqs);
			double randStdDev = calcStdDev(randMean, freqs);
			double targetGraphFreq =
					targetGraphRelFreqs.getOrDefault(label, 0.0);
			double zScore = 0.0;
			if (randStdDev != 0) {
				zScore = (targetGraphFreq - randMean) / randStdDev;
			}
			zScores.put(label, zScore);
		}
            */
	}
        
	private double calcZscore(List<Double> freqs, double targetGraphFreq) {
		double randMean = calcMean(freqs);
		double randStdDev = calcStdDev(randMean, freqs);
		double numerator = targetGraphFreq - randMean;
		
                // debug
                /*
		NumberFormat nf = new DecimalFormat("0.000");
                this.debug.append("Calculating z-score").append(String.format("%n"));
                this.debug.append("randMean").append(nf.format(randMean)).append(String.format("%n"));
                this.debug.append("randStdDev").append(nf.format(randStdDev)).append(String.format("%n"));
                this.debug.append("numerator").append(nf.format(numerator)).append(String.format("%n"));
                */
                
		double zScore = 0.0;
		if (randStdDev != 0) {
			zScore = numerator / randStdDev;
		} else {
			if (numerator > 0)
				zScore = Double.POSITIVE_INFINITY;
			if (numerator == 0)
				zScore = 0.0;
			else
				zScore = Double.NEGATIVE_INFINITY;
		}
		return zScore;
	}
	
        // TODO Don't recalculate randMean
	// calculates z^-scores for this Statistical Analysis object
	private void calculateZDScores(Map<String, List<Double>> directGraphRelFreqs,
	                              Map<String, Double> targetGraphRelFreqs) {
		
		//System.out.println("randGraphRelFreqs..." + randGraphRelFreqs);
		
		HashSet<String> alllabels = new HashSet<String>();
		
		alllabels.addAll(targetGraphRelFreqs.keySet());
		alllabels.addAll(directGraphRelFreqs.keySet());
		for (String label : alllabels) {
			double zScore = 0.0;
			if (!directGraphRelFreqs.containsKey(label))
				zScore = Double.POSITIVE_INFINITY;
			else {
				
				List<Double> freqs = directGraphRelFreqs.get(label);
				double targetGraphFreq = targetGraphRelFreqs.getOrDefault(label, 0.0);
				
				zScore = calcZDscore(freqs, targetGraphFreq);
			}
			zDScores.put(label, zScore);
			
		}
	}
        
        private double calcZDscore(List<Double> freqs, double targetGraphFreq) {
		double randMean = calcMean(freqs);
		double directStdDev = calcStdDev(randMean, freqs);
		double numerator = targetGraphFreq - randMean;
		
                // debug
                /*
		NumberFormat nf = new DecimalFormat("0.000");
                this.debug.append("Calculating z^-score").append(String.format("%n"));
                this.debug.append("randMean").append(nf.format(randMean)).append(String.format("%n"));
                this.debug.append("directStdDev").append(nf.format(directStdDev)).append(String.format("%n"));
                this.debug.append("numerator").append(nf.format(numerator)).append(String.format("%n"));
                */
                
		double zScore = 0.0;
		if (directStdDev != 0) {
			zScore = numerator / directStdDev;
		} else {
			if (numerator > 0)
				zScore = Double.POSITIVE_INFINITY;
			if (numerator == 0)
				zScore = 0.0;
			else
				zScore = Double.NEGATIVE_INFINITY;
		}
		return zScore;
	}
	
	// calculates the standard deviation for a random
	private double calcStdDev(double randMean, List<Double> values) {
		double variance = 0.0;
		for (double value : values) {
			double distance = value - randMean;
			double squaredDistance = Math.pow(distance, 2);
			variance += squaredDistance;
		}
		return Math.sqrt(variance / (values.size() - 1));
	}
	
	// calculates the p-values for this object
	private void calculatePValues(Map<String, List<Double>> randGraphRelFreqs,
	                              Map<String, Double> targetGraphRelFreqs) {
		HashSet<String> alllabels = new HashSet<String>();
		
		alllabels.addAll(targetGraphRelFreqs.keySet());
		alllabels.addAll(randGraphRelFreqs.keySet());
		
		for (String label : alllabels) {
			double pValue = calcPValue(label, randGraphRelFreqs,
					targetGraphRelFreqs);
			pValues.put(label, pValue);
		}
		/*for(String label : randGraphRelFreqs.keySet())
		{
			double pValue = calcPValue(label, randGraphRelFreqs,
					targetGraphRelFreqs);
			pValues.put(label, pValue);
		}*/
	}
	
	// calculates a p-value for an specified label
	private double calcPValue(String label,
	                          Map<String, List<Double>> randomGraphRelFreqs,
	                          Map<String, Double> targetGraphRelFreqs) {
		// if a label appears in the target graph that didn't show up in any
		// random graphs, clearly it's a network motif. This scenario shouldn't
		// happen for a reasonable number of random graphs
		if (!randomGraphRelFreqs.containsKey(label)) {
			return 0;
		}
		
		// This shouldn't happen with the current design, but in case the design
		// changes to include functionality to display all labels found in the
		// random graphs instead of just those in the target graph, this will
		// ensure those labels are not identified as network motifs
		if (!targetGraphRelFreqs.containsKey(label)) {
			return 1;
		}
		
		int prePValue = 0;
		List<Double> randFreqs = randomGraphRelFreqs.get(label);
		double targetFreq = targetGraphRelFreqs.get(label);
		for (double randFreq : randFreqs) {
			if (randFreq >= targetFreq) {
				prePValue++;
			}
		}
		double randFreqCount = randFreqs.size();
		return (double) prePValue / randFreqCount;
	}
	
	//
	private double calcMean(List<Double> values) {
		double total = 0.0;
		for (double randFreq : values) {
			total += randFreq;
		}
		return total / values.size();
	}
	
	@Override
	public String toString() {
		NumberFormat nf = new DecimalFormat("0.000");
		StringBuilder sb = new StringBuilder();
                
                // debug
                //sb.append(this.debug.toString());
                
		sb.append("Label\tRelFreq\t\tRandMeanFreq\t\tDirectMeanFreq\tZ-Score\tP-Value");
                
                sb.append("\tZ^-Value");
                
		sb.append(String.format("%n"));
		for (String label : zScores.keySet()) {
			sb.append(label).append("\t\t");
			if (targetLabelToRelativeFrequencies.containsKey(label)) {
				double targetGraphRelFreqPerc =
						targetLabelToRelativeFrequencies.get(label) * 100.0;
				//System.out.println("label="+label+",Freq="+targetGraphRelFreqPerc);
				sb.append(nf.format(targetGraphRelFreqPerc));
			} else {
				sb.append(nf.format(0.0));
			}
			sb.append("%\t\t");
			if (randomLabelToMeanRelativeFrequencies.containsKey(label)) {
				double randomGraphRelFreqPerc =
						randomLabelToMeanRelativeFrequencies.get(label) * 100.0;
				sb.append(nf.format(randomGraphRelFreqPerc));
			} else {
				sb.append(nf.format(0.0));
			}
                        
                        // direct
                        sb.append("%\t\t");
			if (directLabelToMeanRelativeFrequencies.containsKey(label)) {
				double directGraphRelFreqPerc =
						directLabelToMeanRelativeFrequencies.get(label) * 100.0;
				sb.append(nf.format(directGraphRelFreqPerc));
			} else {
				sb.append(nf.format(0.0));
			}
                        
			sb.append("%\t\t\t");
			double zScore = getZScore(label);
			sb.append(nf.format(zScore)).append("\t");
			double pValue = getPValue(label);
			sb.append(nf.format(pValue));
                        
                        sb.append("\t");
                        // direct method results:
                        double zDScore = getZDScore(label);
			sb.append(nf.format(zDScore));
                        
			sb.append(String.format("%n"));
		}
                
                // subgraph representation
                sb.append("Subgraph representation:").append(String.format("%n"));
                for (String label : rep.keySet()) {
                    sb.append(label).append(":").append(String.format("%n"));
                    
                    sb.append("0|| ");
                    for (int i = 0; i < this.motifSize; i++){
                        sb.append(i + 1);
                    }
                    
                    sb.append(String.format("%n")).append("==");
                    for (int i = 0; i < this.motifSize; i++){
                        sb.append("=");
                    }
                    sb.append(String.format("%n"));
                    
                    int j = 1;
                    for (String line : rep.get(label)){
                        sb.append(j).append("|| ");
                        sb.append(line).append(String.format("%n"));
                        j++;
                    }
                }
                
		return sb.toString();
	}
}