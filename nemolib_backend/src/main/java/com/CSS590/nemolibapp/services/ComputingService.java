package com.CSS590.nemolibapp.services;

import com.CSS590.nemolibapp.configure.FileStorageProperties;
import com.CSS590.nemolibapp.model.FileResponse;
import com.CSS590.nemolibapp.model.ResponseBean;
import edu.uwb.nemolib.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.logging.Level;

/**
 * @author Yangxiao on 3/5/2019.
 */
@Service
public class ComputingService {
	
	private final Logger logger = LogManager.getLogger(ComputingService.class);
	private final Path dirPath;
	private final String dirPathSep;
        
        // debug
        //private StringBuilder debug;
	
	public ComputingService(FileStorageProperties fileStorageProperties) {
		this.dirPath = Paths.get(fileStorageProperties.getWorkDir()).toAbsolutePath().normalize();
		this.dirPathSep = this.dirPath.toString() + File.separator;
		logger.debug("Current workdir: " + this.dirPath);
                
                //this.debug = new StringBuilder();
                //this.debug.append(" - DEBUG - COMPUTINGSERVICE - ").append(String.format("%n"));
                
	}
	
	/**
	 * Provide three output options from here.
	 * Subgraph results can be provided into three options.
	 * 1. SubgraphCount: provide frequency for each pattern
	 * 2. SubgraphProfile: Provide frequency as well as the pattern's concentration on each vertex
	 * 3. SubgraphCollection: Provide frequency as well as the instances of each pattern written in the filename given
	 * <p>
	 * The first option is "SubgraphCount", which provides the frequencies of each pattern.
	 * If the graph or motif size is big, this method is recommended.
	 * To go different option, just comment out all the method from this line until encounter
	 */
	public boolean CalculateNetworkMotif(final String fileName, final int motifSize, final int randGraphCount,
	                                     final boolean directed, final List<Double> prob, final ResponseBean responseBean) {
		final long time = System.currentTimeMillis();
		logger.info("Start CalculateNetworkMotif");
		
		if (motifSize < 3) {
			System.err.println("Motif size must be 3 or larger");
			responseBean.setResults("Motif size must be 3 or larger");
			return false;
		}
		logger.debug("Parsing target graph...");
		final Graph targetGraph;
		try {
			targetGraph = GraphParser.parse(fileName, directed);
		} catch (IOException e) {
			System.err.println("Could not process " + fileName);
			System.err.println(e.getMessage());
			responseBean.setResults("Could not process " + fileName);
			return false;
		}
		
		SubgraphCount subgraphCount = new SubgraphCount();
		SubgraphEnumerator targetGraphESU = new ESU(prob);
		
		TargetGraphAnalyzer trgtGraphAnalyzer = new TargetGraphAnalyzer(targetGraphESU, subgraphCount);
		
		Map<String, Double> tgtFreqMap = trgtGraphAnalyzer.analyze(targetGraph, motifSize);
		
		logger.debug("Target label to relative frequency: " + tgtFreqMap);
		logger.debug("Step 2: Generating " + randGraphCount + " random graph...");
		SubgraphEnumerator randESU = new RandESU(prob);
		logger.debug("Random graph analyze...");
		RandomGraphAnalyzer randGraphAnalyzer = new RandomGraphAnalyzer(randESU, randGraphCount);
		Map<String, List<Double>> randFreqMap = randGraphAnalyzer.analyze(targetGraph, motifSize);
		
                logger.debug("Obtain direct method results");
                // String direct = getDM(directed, targetGraph, motifSize, prob, tgtFreqMap, randFreqMap);
                Map<String, List<Double>> direct = getDM(directed, targetGraph, motifSize, prob, tgtFreqMap, randFreqMap, randGraphCount);
                
		logger.debug("Step 3: Determine network motifs through statistical analysis...");
		RelativeFrequencyAnalyzer relativeFreqAnalyzer = new RelativeFrequencyAnalyzer(randFreqMap, tgtFreqMap, direct, motifSize);
		setRes(responseBean, time, relativeFreqAnalyzer.toString(), targetGraph, subgraphCount.getSubgraphCount());
		return true;
	}
	
	/**
	 * The second option is "SubgraphProfile", which maps each vertex concentration to each pattern.
	 * The frequency of each pattern is also provided
	 * If the graph or motif size is big, this method is recommended.
	 * To go different option, just comment out all the method from this line until encounter
	 */
	public String CalculateNemoProfile(String uuid, String fileName, int motifSize, int randGraphCount,
	                                   boolean directed, List<Double> prob, FileResponse responseBean) {
		logger.info("Start CalculateNemoProfile");
		final long time = System.currentTimeMillis();
		
		if (motifSize < 3) {
			System.err.println("Motif size must be 3 or larger");
			responseBean.setResults("Motif size must be 3 or larger");
			return null;
		}
		
		logger.debug("Parsing target graph...");
		Graph targetGraph;
		try {
			targetGraph = GraphParser.parse(fileName, directed);
		} catch (IOException e) {
			System.err.println("Could not process " + fileName);
			System.err.println(e.getMessage());
			responseBean.setResults("Could not process " + fileName);
			return null;
		}
		
		SubgraphProfile subgraphCount = new SubgraphProfile();
		SubgraphEnumerator targetGraphESU = new ESU(prob);
		
		TargetGraphAnalyzer trgtGraphAnalyzer = new TargetGraphAnalyzer(targetGraphESU, subgraphCount);
		Map<String, Double> tgtFreqMap = trgtGraphAnalyzer.analyze(targetGraph, motifSize);
		logger.debug("Target label to relative frequency: " + tgtFreqMap);
		logger.debug("Step 2: Generating " + randGraphCount + " random graph...");
		SubgraphEnumerator randESU = new RandESU(prob);
		logger.debug("Random graph analyze...");
		RandomGraphAnalyzer randGraphAnalyzer = new RandomGraphAnalyzer(randESU, randGraphCount);
		
		Map<String, List<Double>> randFreqMap = randGraphAnalyzer.analyze(targetGraph, motifSize);
                
                logger.debug("Obtain direct method results");
                Map<String, List<Double>> direct = getDM(directed, targetGraph, motifSize, prob, tgtFreqMap, randFreqMap, randGraphCount);
                // String direct = getDM(directed, targetGraph, motifSize, prob, tgtFreqMap, randFreqMap);
                
		// logger.debug("random Label To Relative Frequencies=" + randFreqMap);
		logger.debug("Step 3: Determine network motifs through statistical analysis...");
		RelativeFrequencyAnalyzer relativeFreqAnalyzer = new RelativeFrequencyAnalyzer(randFreqMap, tgtFreqMap, direct, motifSize);
		
		String resFileName = uuid + "_subProfile.txt";
		logger.debug("Building results based on pvalue < 0.05 " + resFileName);
		// NemoProfileBuilder.buildwithPvalue(subgraphCount, relativeFreqAnalyzer,
		// 		0.05, this.dirPathSep + resFileName, targetGraph.getNameToIndexMap());
		NemoProfileBuilder.buildwithPvalue(subgraphCount, relativeFreqAnalyzer,
				1.0, this.dirPathSep + resFileName, targetGraph.getNameToIndexMap());
		Map<String, Integer> freqVector = NemoProfileBuilder.getNemoFrequencyVector(subgraphCount, relativeFreqAnalyzer,
				1.0, targetGraph.getNameToIndexMap(), motifSize);
		
                
		logger.trace("SubgraphProfile Compete");
		setResWithMap(responseBean, time, relativeFreqAnalyzer.toString(), freqVector, targetGraph, subgraphCount.getSubgraphCount());
		return setFileRes(resFileName, responseBean);
	}
	
	private void setResWithMap(FileResponse responseBean, long time, String relaFreqAna, Map<String, Integer> freqVector, 
                    Graph targetGraph, int subGraphs) {
		StringBuilder freqV = new StringBuilder();
		for (Map.Entry<String, Integer> entry : freqVector.entrySet()) {
			freqV.append(entry.getKey()).append(": ").append(entry.getValue()).append("\t");
		}
		responseBean.setResults("Running time = " + (System.currentTimeMillis() - time) + "ms\n" + 
                        "Nodes: " + targetGraph.getSize() + "\n" + 
                        "Edges: " + targetGraph.getEdgeCount() + "\n" +
                        "Subgraph Count: " + subGraphs + "\n" + relaFreqAna + "\n" + freqV.toString());
                        //+ "\n" + this.debug.toString());
                
                //this.debug = new StringBuilder();
                //this.debug.append(" - DEBUG - COMPUTINGSERVICE - ").append(String.format("%n"));
	}
	
	/**
	 * The last option is "SubgraphCollection", which write all instances of each pattern, and frequency of each
	 * pattern
	 * It is recommended to use for moderate graph size or motif size.
	 * To go different option, just comment out all the method from this line until encounter 33333333333333333333333333333333333333333
	 */
	public String CalculateNemoCollection(String uuid, String fileName, int motifSize, int randGraphCount,
	                                      boolean directed, List<Double> prob, FileResponse responseBean) {
		logger.info("Start CalculateNemoCollection");
		final long time = System.currentTimeMillis();
		
		if (motifSize < 3) {
			System.err.println("Motif size must be 3 or larger");
			responseBean.setResults("Motif size must be 3 or larger");
			return null;
		}
		
		logger.debug("Parsing target graph...");
		Graph targetGraph;
		try {
			targetGraph = GraphParser.parse(fileName, directed);
		} catch (IOException e) {
			System.err.println("Could not process " + fileName);
			System.err.println(e.getMessage());
			responseBean.setResults("Could not process " + fileName);
			return null;
		}
		
		// If want to save the name to index map to the file
		// targetGraph.write_nametoIndex("Name_Index.txt");
		
		// If want to provide collections with instances written "Results.txt" file.
		// SubgraphCollection subgraphCount = new SubgraphCollection("Results.txt");
		
		// Default file name is "SubgraphCollectionResult.txt"
		
		// Create subgraphCount instance which will collect results in SubgraphCollectionG6.txt
		SubgraphCollection subgraphCount = new SubgraphCollection(this.dirPathSep + uuid + "_subG6.txt");
		
		SubgraphEnumerator targetGraphESU = new ESU(prob);
		TargetGraphAnalyzer trgtGraphAnalyzer = new TargetGraphAnalyzer(targetGraphESU, subgraphCount);
		Map<String, Double> tgtFreqMap = trgtGraphAnalyzer.analyze(targetGraph, motifSize);
		logger.debug("Target label to relative frequency: " + tgtFreqMap);
		logger.debug("Step 2: Generating " + randGraphCount + " random graph...");
		SubgraphEnumerator randESU = new RandESU(prob);
		RandomGraphAnalyzer randGraphAnalyzer = new RandomGraphAnalyzer(randESU, randGraphCount);
		Map<String, List<Double>> randFreqMap = randGraphAnalyzer.analyze(targetGraph, motifSize);
		
                logger.debug("Obtain direct method results");
                // String direct = getDM(directed, targetGraph, motifSize, prob, tgtFreqMap, randFreqMap);
                Map<String, List<Double>> direct = getDM(directed, targetGraph, motifSize, prob, tgtFreqMap, randFreqMap, randGraphCount);
                
		// logger.debug("random Label To Relative Frequencies=" + randFreqMap);
		logger.debug("Step 3: Determine network motifs through statistical analysis...");
		RelativeFrequencyAnalyzer relativeFreqAnalyzer = new RelativeFrequencyAnalyzer(randFreqMap, tgtFreqMap, direct, motifSize);
		
		// This is optional, if the user want to collect all subgraphs with canonical label in a file
		// Write the nemocollection result based on zscore thresh (anything with >=2 is collected) .
		// logger.trace("Write the nemocollection result based on zscore thresh (anything with >=2 is collected) .");
		// NemoCollectionBuilder.buildwithZScore(subgraphCount, relativeFreqAnalyzer,
		// 		2.0, "NemoCollectionZscore.txt", targetGraph.getNameToIndexMap());
		
		// logger.trace("Write the nemocollection result based on pvalue thresh (anything with <0.05 is collected).");
		// NemoCollectionBuilder.buildwithPvalue(subgraphCount, relativeFreqAnalyzer,
		// 		0.05, "NemoCollectionPValue.txt", targetGraph.getNameToIndexMap());
		
		String resFileName = uuid + "_subCol.txt";
		logger.trace("Write the subgraph collection to " + resFileName);
		NemoCollectionBuilder.buildwithPvalue(subgraphCount, relativeFreqAnalyzer,
				1.0, this.dirPathSep + resFileName, targetGraph.getNameToIndexMap());
		
		logger.trace("NemoCollection Compete");
		
                
		setRes(responseBean, time, relativeFreqAnalyzer.toString(), targetGraph, subgraphCount.getSubgraphCount());
		return setFileRes(resFileName, responseBean);
	}
	
	private void setRes(ResponseBean responseBean, long time, String relaFreqAna, Graph targetGraph, int subGraphs) {
		responseBean.setResults("Running time = " + (System.currentTimeMillis() - time) + "ms\n" + 
                        "Nodes: " + targetGraph.getSize() + "\n" + 
                        "Edges: " + targetGraph.getEdgeCount() + "\n" +
                        "Subgraph Count: " + subGraphs + "\n" + relaFreqAna);
                        //+ "\n" + this.debug.toString());
                
                //this.debug = new StringBuilder();
                //this.debug.append(" - DEBUG - COMPUTINGSERVICE - ").append(String.format("%n"));
	}
	
	private String setFileRes(String resFileName, FileResponse responseBean) {
		File file = new File(this.dirPathSep + resFileName);
		if (file.exists()) {
			responseBean.setFilename(resFileName);
			responseBean.setSize(file.length() / 1024);
			return resFileName;
		} else {
			return "no";
		}
	}
        
        // get the results of the direct method
        public Map<String, List<Double>> getDM(boolean directed, Graph inputGraph, int motifSize, List<Double> prob, 
                Map<String, Double> tgtFreqMap, Map<String, List<Double>> randFreqMap, int graphCount)
        {
            // copied from DirectMethodBuilderTest.java
            // modified for use with nemolib app
            
            //NumberFormat nf = new DecimalFormat("0.000");
            
            List<String> g6s = new ArrayList<>();
            List<Long> motifs = new ArrayList<>();
            Map<String, Long> g6_to_number = new HashMap<>();
            Map<Long, String> number_to_g6 = new HashMap<>();

            Map<String, List<Double>> g6_to_res = new HashMap<>();

            HashSet<String> alllabels = new HashSet<String>();
            alllabels.addAll(tgtFreqMap.keySet());
            alllabels.addAll(randFreqMap.keySet());
            
            logger.debug("Writing to subgraphs.txt");
            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter("subgraphs.txt"));
                for (String key : alllabels) {
                    writer.write(key);
                    g6s.add(key);
                    writer.write('\n');
                }
                writer.close();
            }catch (IOException e){
                System.err.println(e.getMessage());
                //return "Error: subgraphs.txt write failed";
                return g6_to_res;
            }


            logger.debug("Running showg with subgraphs.txt");
            Runtime rt = Runtime.getRuntime();
            Process proc;
            try {
                proc = rt.exec("./showg -a subgraphs.txt");
            }catch (IOException e){
                System.err.println(e.getMessage());
                //return "Error: showg.exe did not work or path to showg.exe is wrong";
                return g6_to_res;
            }
            
            logger.debug("Waiting for showg to finish");
            try {
                proc.waitFor();
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
            }
            
            BufferedReader stdInput = new BufferedReader(new
                    InputStreamReader(proc.getInputStream()));

            logger.debug("Reading output from showg");
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
                return g6_to_res;
            }

            logger.debug("Analyzing output from showg");
            int k = 0;
            for(String g6: g6s){
                logger.debug(g6);
                List<String> adj = new ArrayList<>();
                for(int i = 0; i <= motifSize; i++){
                    logger.debug("i = " + i);
                    logger.debug("checking for \"Graph\"");
                    if(adjMatrixs.get(k * (motifSize + 1) + i).startsWith("Graph")){
                        logger.debug(adjMatrixs.get(k * (motifSize + 1) + i));
                        continue;
                    }
                    logger.debug("adding to adj");
                    logger.debug(adjMatrixs.get(k * (motifSize + 1) + i));
                    adj.add(adjMatrixs.get(k * (motifSize + 1) + i));
                }
                logger.debug("loop done, k = " + k);
                k++;
                logger.debug("Setting number_code with adjMatrixToNumber(adj)");
                long number_code = adjMatrixToNumber(adj);
                logger.debug("number_code set");
                
                g6_to_number.put(g6, number_code);
                number_to_g6.put(number_code, g6);
                motifs.add(number_code);
            }
            
            logger.debug("Using direct method builder to analyze results");
            
            long[] degree_r = DirectMethodBuilder.getDegrees(inputGraph);
            long[] degree_c = DirectMethodBuilder.getDegrees(inputGraph);

            //int num_samples = 100000;
            int num_samples = 100;
            
            DirectMethodBuilder directMethodBuilder;
            
            if(inputGraph.getDir()){
                directMethodBuilder = new DirectMethodBuilder(inputGraph.getSize(), degree_r, degree_c, new Random());
            }else{
                directMethodBuilder = new DirectMethodBuilder(inputGraph.getSize(), degree_r, new Random());
            }
            
            double[][] res = new double[motifs.size()][graphCount]; //result
            
            for(int i = 0; i != motifs.size(); i++){
                
                for (int j = 0; j != graphCount; j++){
                    logger.debug("getting motif sample log i = " + i + " j = " + j);
                    res[i][j] = directMethodBuilder.motif_sample_log(motifs.get(i), (short) motifSize,  num_samples);
                    
                }
            }
            
            logger.debug("testing for NaN");
            double max = -50000000.0;
            for (int i = 0; i!= motifs.size(); ++i) {
                for (int j = 0; j != graphCount; j++){
                    if (res[i][j] > max && !(res[i][j] != res[i][j]))  //test for NaN!
                    { max = res[i][j];	}
                }
            }
            
            logger.debug("setting up sum");
            // declare and initialize sum
            double[] sum = new double[graphCount];
            for (int i = 0; i != graphCount; i++){
                sum[i] = 0.0;
            }
            
            for (int i = 0; i!= motifs.size(); ++i) {
                for (int j = 0; j != graphCount; j++){
                    if (!(res[i][j] != res[i][j]))  //test for NaN!
                    {
                        //cout << res[i];
                        res[i][j] = Math.exp(res[i][j] - max);
                        sum[j] += res[i][j];
                    }
                    else
                        res[i][j] = 0.0;
                }
            }
            
            logger.debug("finishing relative frequencies");
            for (int i = 0; i!= motifs.size(); ++i) {
                String mot = number_to_g6.get(motifs.get(i));
                
                List<Double> temp = new ArrayList<>();
                
                
                for (int j = 0; j != graphCount; j++){
                    temp.add(res[i][j] / sum[j]);
                }
                g6_to_res.put(mot, temp);
                
            }
            
            logger.debug("Direct method finished");
            
            // return the direct method results
            return g6_to_res;
        }
        
        public static long adjMatrixToNumber(List<String> adj){
            StringBuilder num_code = new StringBuilder();
            for(String line: adj){
                num_code.append(line);
            }
            return Long.parseLong(num_code.toString(), 2);
        }
}
