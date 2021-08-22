package com.CSS590.nemolibapp.services;

import com.CSS590.nemolibapp.configure.FileStorageProperties;
import com.CSS590.nemolibapp.support.MyFileNotFoundException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.core.io.UrlResource;
import org.springframework.stereotype.Service;
import org.springframework.util.FileSystemUtils;
import org.springframework.util.StringUtils;
import org.springframework.web.multipart.MultipartFile;
import org.springframework.core.io.Resource;

import java.io.IOException;
import java.net.MalformedURLException;
import java.nio.file.*;
import java.util.stream.Stream;

/**
 * @author Yangxiao on 3/5/2019.
 */
@Service
public class StorageService {
	private final Path dirPath;
	private final Path downLoadPath;
	private final Path tempPath;
	Logger logger = LogManager.getLogger(StorageService.class);
	
	public StorageService(FileStorageProperties fileStorageProperties) {
		
		this.tempPath = Paths.get("/tmp/nemoback/upload_tmp").toAbsolutePath().normalize();
		this.dirPath = Paths.get(fileStorageProperties.getUploadDir()).toAbsolutePath().normalize();
		this.downLoadPath = Paths.get(fileStorageProperties.getWorkDir()).toAbsolutePath().normalize();
		
		try {
			Files.createDirectories(this.tempPath);
			Files.createDirectories(this.dirPath);
			Files.createDirectories(this.downLoadPath);
		} catch (Exception ex) {
			logger.error("Could not create the directory : " + ex.getMessage());
		}
		
		// try {
		// 	Files.createDirectories(this.downLoadPath);
		// } catch (Exception ex) {
		// 	logger.error("Could not create the directory : " + ex.getMessage());
		// }
		//
		// try {
		// 	Files.createDirectories(this.dirPath);
		// } catch (Exception ex) {
		// 	logger.error("Could not create the directory : " + ex.getMessage());
		// }
	}
	
	public Path storeFile(MultipartFile file) {
		// Normalize file name
		String fileName = StringUtils.cleanPath(file.getOriginalFilename());
		
		try {
			// Check if the file's name contains invalid characters
			if (fileName.contains("..")) {
				System.out.println("Sorry! Filename contains invalid path sequence " + fileName);
			}
			// Copy file to the target location (Replacing existing file with the same name)
			Path targetLocation = this.dirPath.resolve(fileName);
			Files.copy(file.getInputStream(), targetLocation, StandardCopyOption.REPLACE_EXISTING);
			return targetLocation;
		} catch (IOException ex) {
			logger.error("Could not store file: " + fileName + ". Please try again!\n" + ex.getMessage());
			return null;
		}
	}
	
	Stream<Path> loadAll() {
		return null;
	}
	
	Path loadFile(String filename) {
		return null;
	}
	
	public Resource loadAsResource(String fileName) {
		try {
			Path file = downLoadPath.resolve(fileName);
			Resource resource = new UrlResource(file.toUri());
			if (resource.exists() || resource.isReadable()) {
				return resource;
			} else {
				logger.error("Cannot load file " + fileName + ", file does not exist or not readable");
				throw new MyFileNotFoundException("File not found " + fileName);
			}
		} catch (MalformedURLException ex) {
			logger.error("Cannot load file " + fileName);
			throw new MyFileNotFoundException("File not found " + fileName, ex);
		}
	}
	
	private void deleteAll() {
		FileSystemUtils.deleteRecursively(dirPath.toFile());
	}
	
}
