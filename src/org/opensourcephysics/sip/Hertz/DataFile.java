package org.opensourcephysics.sip.Hertz;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Calendar;

public class DataFile {
	private static Calendar calendar;
	private final static String baseDir = "data";
	private String initialConfigurations;
	private String name;
	private String data;
	private BufferedWriter bw;
	private Path filePath;
	
	public DataFile(String name, String configurations){
		this.name = name;
		if(DataFile.calendar == null){
			calendar = Calendar.getInstance();
		}
		String creationDate = calendar.get(Calendar.DAY_OF_MONTH) + " " + 
							  (calendar.get(Calendar.MONTH) + 1) + " " + // calendar.MONTH starts counting at 0
							  calendar.get(Calendar.YEAR);
		String creationTime = calendar.get(Calendar.HOUR_OF_DAY) + ":" + 
							  calendar.get(Calendar.MINUTE);
		Path dir = Paths.get(DataFile.baseDir + "/" + creationDate);
		
		// create directory to store datafiles, categorized by date if it doesn't exist
		try {
			Files.createDirectories(dir);
		} catch (IOException e) {
			e.printStackTrace();
		}

		filePath = Paths.get(dir + "/" + creationTime + " " + name + ".dat");
		
		try {
			bw = Files.newBufferedWriter(filePath, 
					Charset.forName("US-ASCII"), filePath.toFile().exists()? StandardOpenOption.APPEND : StandardOpenOption.CREATE);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		record(configurations);
	}
	
	/**
	 * Adds the data to the buffer.
	 * @param data Data to be added
	 */
	public void record(String data){
		try {
			bw.write(data);
			bw.newLine();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			System.out.println("record method called without calling constructor!");
		}
	}
	
	/**
	 * Writes the data to the file.
	 */
	public void write(){
		try{
			bw.flush();
		} catch (IOException e){
			e.printStackTrace();
			System.out.println("Error writing to file.");
		}
	}
	
	/**
	 * Closes the stream for writing.
	 */
	public void close(){
		try{
			bw.close();
		} catch(IOException e){
			e.printStackTrace();
			System.out.println("Problem closing the buffer stream.");
		}
	}
}
