package com.CSS590.nemolibapp.model;

/**
 * @author Yangxiao on 3/13/2019.
 */
public class NetworkMotifResponse implements ResponseBean {
	
	private String message;
	private String results;
	private String optional;
	
	public NetworkMotifResponse(int motifSize, int randSize, boolean directed, String fileName) {
		this.message = "Motif size: " + motifSize + "\n" +
				"Random graph size: " + randSize + "\n" +
				"Directed: " + directed + "\n" +
				"File: " + fileName + "\n";
	}
	
	public NetworkMotifResponse() {
	}
	
	public static NetworkMotifResponse initWithMessage(String message) {
		NetworkMotifResponse res = new NetworkMotifResponse();
		res.setMessage(message);
		return res;
	}
	
	@Override
	public void setMessage(String message) {
		this.message = message;
	}
	
	@Override
	public String getMessage() {
		return message;
	}
	
	@Override
	public void setResults(String results) {
		this.results = results;
	}
	
	@Override
	public String getResults() {
		return results;
	}
	
	@Override
	public void setOptional(String optional) {
		this.optional = optional;
	}
	
	@Override
	public String getOptional() {
		return optional;
	}
}
