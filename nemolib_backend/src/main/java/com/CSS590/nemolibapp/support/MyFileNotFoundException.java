package com.CSS590.nemolibapp.support;

import org.springframework.http.HttpStatus;
import org.springframework.web.bind.annotation.ResponseStatus;

/**
 * @author Yangxiao Wang on 7/27/2019
 */
@ResponseStatus(HttpStatus.NOT_FOUND)
public class MyFileNotFoundException extends RuntimeException {
	public MyFileNotFoundException(String message) {
		super(message);
	}
	
	public MyFileNotFoundException(String message, Throwable cause) {
		super(message, cause);
	}
}
