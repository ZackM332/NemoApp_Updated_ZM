package com.CSS590.nemolibapp.model;

/**
 * @author Yangxiao on 3/14/2019.
 */

public interface ResponseBean {
	
	void setMessage(String message);
	
	String getMessage();
	
	void setResults(String results);
	
	String getResults();
	
	void setOptional(String optional);
	
	String getOptional();
}
