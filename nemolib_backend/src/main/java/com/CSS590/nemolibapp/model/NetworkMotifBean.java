package com.CSS590.nemolibapp.model;

import org.springframework.web.multipart.MultipartFile;

/**
 * @author Yangxiao on 3/13/2019.
 */
public class NetworkMotifBean {
	
	private int motifSize;
	private int randSize;
	private MultipartFile file;
	
	public int getMotifSize() {
		return motifSize;
	}
	
	public void setMotifSize(int motifSize) {
		this.motifSize = motifSize;
	}
	
	public int getRandSize() {
		return randSize;
	}
	
	public void setRandSize(int randSize) {
		this.randSize = randSize;
	}
	
	public MultipartFile getFile() {
		return file;
	}

	public void setFile(MultipartFile file) {
		this.file = file;
	}
}
