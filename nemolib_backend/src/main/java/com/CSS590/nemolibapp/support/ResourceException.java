package com.CSS590.nemolibapp.support;

import org.springframework.http.HttpStatus;

/**
 * @author Yangxiao Wang on 7/27/2019
 */
public class ResourceException extends RuntimeException {
	
	private final HttpStatus httpStatus;
	private final String message;
	
	/**
	 * Constructs a new runtime exception with the specified HttpStatus code and detail message.
	 * The cause is not initialized, and may subsequently be initialized by a call to {@link #initCause}.
	 *
	 * @param httpStatus the http status.  The detail message is saved for later retrieval by the {@link
	 *                   #getHttpStatus()} method.
	 * @param message    the detail message. The detail message is saved for later retrieval by the {@link
	 *                   #getMessage()} method.
	 * @see HttpStatus
	 */
	public ResourceException(HttpStatus httpStatus, String message) {
		this.message = message;
		this.httpStatus = httpStatus;
	}
	
	/**
	 * Gets the HTTP status code to be returned to the calling system.
	 *
	 * @return http status code.  Defaults to HttpStatus.INTERNAL_SERVER_ERROR (500).
	 * @see HttpStatus
	 */
	public HttpStatus getHttpStatus() {
		return httpStatus;
	}
	
	@Override
	public String getMessage() {
		return message;
	}
}
