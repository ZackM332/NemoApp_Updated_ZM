package com.CSS590.nemolibapp.model;

/**
 * @author Yangxiao Wang on 7/28/2019
 */
public class FileResponse extends NetworkMotifResponse {
	
	private String message;
	private String results;
	private String filename;
	private String url;
	private long size;
	
	
	public FileResponse(int motifSize, int randSize, boolean directed, String fileName) {
		super(motifSize, randSize, directed, fileName);
	}
	
	public FileResponse() {
	}
	
	public String getFilename() {
		return filename;
	}
	
	public void setFilename(String filename) {
		this.filename = filename;
	}
	
	public String getUrl() {
		return url;
	}
	
	public void setUrl(String url) {
		this.url = url;
	}
	
	public long getSize() {
		return size;
	}
	
	public void setSize(long size) {
		this.size = size;
	}
}
