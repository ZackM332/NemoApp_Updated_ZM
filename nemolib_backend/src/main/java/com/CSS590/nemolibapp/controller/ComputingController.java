package com.CSS590.nemolibapp.controller;

import com.CSS590.nemolibapp.model.FileResponse;
import com.CSS590.nemolibapp.model.NetworkMotifResponse;
import com.CSS590.nemolibapp.model.ResponseBean;
import com.CSS590.nemolibapp.services.ComputingService;
import com.CSS590.nemolibapp.services.StorageService;
import com.CSS590.nemolibapp.support.ResourceException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.io.Resource;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpStatus;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import org.springframework.web.servlet.support.ServletUriComponentsBuilder;

import javax.servlet.http.HttpServletRequest;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;


/**
 * @author Yangxiao on 3/5/2019.
 */
@RestController
@RequestMapping("compute")
public class ComputingController {
	
	private final Logger logger = LogManager.getLogger(ComputingController.class);
	
	@Autowired
	private StorageService storageService;
	
	@Autowired
	private ComputingService cService;
	
	@CrossOrigin()
	@ExceptionHandler(ResourceException.class)
	@RequestMapping(value = "networkmotif", method = RequestMethod.POST)
	public ResponseBean getNetworkMotif(@RequestParam(name = "motifSize") int motifSize,
	                                    @RequestParam(name = "randSize") int randGraph,
	                                    @RequestParam(name = "directed") int directed,
	                                    @RequestParam(name = "file") MultipartFile file,
	                                    @RequestParam(value = "prob[]") Double[] prob) {
		logger.debug("get network motif");
		
		if (file == null || file.isEmpty()) {
			throw new ResourceException(HttpStatus.BAD_REQUEST, "Error! File is empty");
		}
		logger.trace("Start getNetworkMotif, file name " + file.getOriginalFilename());
		
		ResponseBean responseBean = new NetworkMotifResponse(motifSize, randGraph,
				directed == 1, file.getOriginalFilename());
		Path filePath = storageService.storeFile(file);
		
		if (filePath == null) {
			responseBean.setResults("Error, cannot upload file: " + file.getName());
			throw new ResourceException(HttpStatus.BAD_REQUEST,
					"Error, cannot upload file" + file.getOriginalFilename());
		}
		
		List<Double> probs = Arrays.asList(prob);
		boolean success = cService.CalculateNetworkMotif(filePath.toString(), motifSize, randGraph,
				directed == 1, probs, responseBean);
		if (!success) {
			throw new ResourceException(HttpStatus.BAD_REQUEST, "Unknown error");
		}
		logger.trace("Return result");
		logger.info(responseBean.getMessage() + '\n' + responseBean.getResults());
		return responseBean;
	}
	
	@CrossOrigin()
	@ExceptionHandler(ResourceException.class)
	@RequestMapping(value = "nemoprofile", method = RequestMethod.POST)
	public FileResponse getNemoProfile(@RequestParam(name = "uuid") String uuid,
	                                   @RequestParam(name = "motifSize") int motifSize,
	                                   @RequestParam(name = "randSize") int randGraph,
	                                   @RequestParam(name = "directed") int directed,
	                                   @RequestParam(name = "file") MultipartFile file,
	                                   @RequestParam(value = "prob[]") Double[] prob) {
		logger.debug("get NemoProfile with id: " + uuid);
		
		if (file == null || file.isEmpty()) {
			throw new ResourceException(HttpStatus.BAD_REQUEST, "Error! File is empty");
		}
		FileResponse responseBean = new FileResponse(motifSize, randGraph,
				directed == 1, file.getOriginalFilename());
		Path filePath = storageService.storeFile(file);
		
		if (filePath == null) {
			responseBean.setResults("Error, cannot upload file: " + file.getName());
			throw new ResourceException(HttpStatus.BAD_REQUEST,
					"Error, cannot upload file" + file.getOriginalFilename());
		}
		
		List<Double> probs = Arrays.asList(prob);
		String filename = cService.CalculateNemoProfile(uuid, filePath.toString(), motifSize, randGraph,
				directed == 1, probs, responseBean);
		logger.info(responseBean.getMessage() + '\n' + responseBean.getResults());
		return processResults(filename, responseBean);
	}
	
	@CrossOrigin()
	@ExceptionHandler(ResourceException.class)
	@RequestMapping(value = "nemocollect", method = RequestMethod.POST)
	public FileResponse getNemoCollection(@RequestParam(name = "uuid") String uuid,
	                                      @RequestParam(name = "motifSize") int motifSize,
	                                      @RequestParam(name = "randSize") int randGraph,
	                                      @RequestParam(name = "directed") int directed,
	                                      @RequestParam(name = "file") MultipartFile file,
	                                      @RequestParam(value = "prob[]") Double[] prob) {
		logger.debug("get NemoCollection with id: " + uuid);
		
		if (file == null || file.isEmpty()) {
			throw new ResourceException(HttpStatus.BAD_REQUEST, "Error! File is empty");
		}
		FileResponse responseBean = new FileResponse(motifSize, randGraph,
				directed == 1, file.getOriginalFilename());
		Path filePath = storageService.storeFile(file);
		
		if (filePath == null) {
			responseBean.setResults("Error, cannot upload file: " + file.getName());
			throw new ResourceException(HttpStatus.BAD_REQUEST,
					"Error, cannot upload file" + file.getOriginalFilename());
		}
		
		List<Double> probs = Arrays.asList(prob);
		String filename = cService.CalculateNemoCollection(uuid, filePath.toString(), motifSize, randGraph,
				directed == 1, probs, responseBean);
		logger.info(responseBean.getMessage() + '\n' + responseBean.getResults());
		return processResults(filename, responseBean);
	}
	
	@CrossOrigin()
	@RequestMapping(value = "/downloadFile/{fileName}", method = RequestMethod.GET)
	public ResponseEntity<Resource> downloadFile(@PathVariable("fileName") String fileName,
	                                             HttpServletRequest request) {
		logger.debug("Start downloadFile " + fileName);
		
		// Load file as Resource
		Resource resource = storageService.loadAsResource(fileName);
		
		// Try to determine file's content type
		String contentType = null;
		try {
			contentType = request.getServletContext().getMimeType(resource.getFile().getAbsolutePath());
		} catch (IOException ex) {
			logger.info("Could not determine file type.");
		}
		
		// Fallback to the default content type if type could not be determined
		if (contentType == null) {
			contentType = "octet-stream";
		}
		
		return ResponseEntity.ok()
				.contentType(MediaType.parseMediaType(contentType))
				.header(HttpHeaders.CONTENT_DISPOSITION, "attachment; filename=\"" + resource.getFilename() + "\"")
				.body(resource);
	}
	
	private FileResponse processResults(String filename, FileResponse responseBean) {
		if (filename == null) {
			throw new ResourceException(HttpStatus.BAD_REQUEST, "Unknown error");
		} else if (filename.equals("no")) {
			responseBean.setOptional("No file generated");
			return responseBean;
		}
		setFileDownloadUrl(responseBean, filename);
		logger.trace("Return result");
		return responseBean;
	}
	
	private void setFileDownloadUrl(FileResponse responseBean, String name) {
		String fileDownloadUri = ServletUriComponentsBuilder.fromCurrentContextPath()
				.path("/compute/downloadFile/")
				.path(name)
				.toUriString();
		logger.debug(fileDownloadUri);
		responseBean.setUrl(fileDownloadUri);
	}
}
