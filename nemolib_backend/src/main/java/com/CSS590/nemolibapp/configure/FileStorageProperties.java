package com.CSS590.nemolibapp.configure;

import org.springframework.boot.context.properties.ConfigurationProperties;

/**
 * @author Yangxiao on 3/6/2019.
 */

@ConfigurationProperties(prefix = "file")
public class FileStorageProperties {
	private String uploadDir;
	private String workDir;
	public String getUploadDir() {
		return uploadDir;
	}
	
	public void setUploadDir(String uploadDir) {
		this.uploadDir = uploadDir;
	}
	
	public String getWorkDir() {
		return workDir;
	}
	
	public void setWorkDir(String workDir) {
		this.workDir = workDir;
	}
	
}
