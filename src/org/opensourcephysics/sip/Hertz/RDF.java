package org.opensourcephysics.sip.Hertz;

import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.opensourcephysics.numerics.PBC;

/**
 * RDF calculates the radial distribution function from a Monte Carlo simulation
 * 
 * @authors Alan Denton and Matt Urich
 * @version 1.0 24-01-2015
 * 
 */

public class RDF {
    private int N;
    private double [] x, y, z;
    private double [] gr;
    private double side;
    private double radBinWidth;
    private int bins;
    private double maxR;
    private int count = 0;
    private String extension; // for writing the file
    
    public RDF(double x[], double y[], double z[], double side, double radBinWidth, String extension){
	N = x.length;
	this.x = x;
	this.y = y;
	this.z = z;
	this.side = side;
	this.radBinWidth = radBinWidth;
	this.extension = extension; 
	maxR = side/2.0;
	bins = (int) (maxR/radBinWidth) + 1;
	gr = new double [bins];
    }
    
    public void update(){
	count++;
	for(int i = 0; i < N; i++){
	    for(int j = 0; j < N; j++){
		if(i != j){
		    double rx = PBC.separation(Math.abs(x[i]-x[j]), side);
		    double ry = PBC.separation(Math.abs(y[i]-y[j]), side);
		    double rz = PBC.separation(Math.abs(z[i]-z[j]), side);
		    double r2=rx*rx+ry*ry+rz*rz;
		    if(r2 < (maxR*maxR)){
			double r = Math.sqrt(r2);
			int bin = (int) (r/radBinWidth);
			gr[bin]++;
		    }
		}
	    }
	}
    }
    
    public void normalize(){
	double rho = (double) N/(side*side*side); // number density
	double norm = 4*Math.PI*rho*Math.pow(radBinWidth,3)*N*count; // normalization
	
	for(int i = 0; i < gr.length; i++){
	    double r = radBinWidth*i; // approximate value of radius in bin i (rounded to bottom of bin)
	    if(r <= maxR){
		gr[i] = gr[i] / (norm*(i+1)*(i+1));
	    }
	}
    }

    public void writeRDF(){ // write normalized gr data to a space delimited file

	normalize();

	try{
	    File radFile = new File("data/radialdistribution" + extension + ".txt");
	    
	    if (!radFile.exists()) { // if file doesn't exist, create it
		radFile.createNewFile();
	    }
	    
	    FileWriter fwrite = new FileWriter(radFile.getAbsoluteFile());
	    BufferedWriter bwrite = new BufferedWriter(fwrite);
	    for (int i=0;i<gr.length;i++){
		bwrite.write(i + " " + gr[i]);
		bwrite.newLine();
	    }
	    
	    bwrite.close();
	}
	
	catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
