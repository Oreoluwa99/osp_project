package org.opensourcephysics.sip.HertzMixture;

import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.opensourcephysics.numerics.PBC;

/**
 * RDFmixture computes radial distribution functions from a Monte Carlo simulation of a microgel mixture. 
 * 
 * @author Alan Denton and Matthew Urich
 * @version 1.1 09-11-2017
 * 
 */

public class RDFmixture {
    private int N1, N2; // numbers of particles of species 1 and 2
    private double [] x1, y1, z1, x2, y2, z2; // particle coordinates
    private double [] gr11, gr12, gr22; // radial distribution functions
    private double side; // box side length
    private double radBinWidth; // histogram bin width
    private int bins; // number of bins in histogram
    private double maxR; // range of histogram
    private int count = 0;
    private String extension; // for writing the data file
    
    public RDFmixture(double x1[], double y1[], double z1[], double x2[], double y2[], double z2[], double side, double radBinWidth, String extension){
	N1 = x1.length;
	N2 = x2.length;
	this.x1 = x1;
	this.y1 = y1;
	this.z1 = z1;
	this.x2 = x2;
	this.y2 = y2;
	this.z2 = z2;
	this.side = side;
	this.radBinWidth = radBinWidth;
	this.extension = extension; 
	maxR = side/2.; // range = half the box length
	bins = (int) (maxR/radBinWidth) + 1;
	gr11 = new double [bins];
	gr12 = new double [bins];
	gr22 = new double [bins];
    }
    
    public void update(){
	count++;
	for(int i = 0; i < N1; i++){ // histogram for g11(r)
	    for(int j = 0; j < N1; j++){
		if(i != j){
		    double rx = PBC.separation(Math.abs(x1[i]-x1[j]), side);
		    double ry = PBC.separation(Math.abs(y1[i]-y1[j]), side);
		    double rz = PBC.separation(Math.abs(z1[i]-z1[j]), side);
		    double r2=rx*rx+ry*ry+rz*rz;
		    if(r2 < maxR*maxR){
			double r = Math.sqrt(r2);
			int bin = (int) (r/radBinWidth);
			gr11[bin]++;
		    }
		}
	    }
	}

	for(int i = 0; i < N1; i++){ // histogram for g12(r) [same as g21(r)]
	    for(int j = 0; j < N2; j++){
		double rx = PBC.separation(Math.abs(x1[i]-x2[j]), side);
		double ry = PBC.separation(Math.abs(y1[i]-y2[j]), side);
		double rz = PBC.separation(Math.abs(z1[i]-z2[j]), side);
		double r2=rx*rx+ry*ry+rz*rz;
	        if(r2 < maxR*maxR){
		    double r = Math.sqrt(r2);
		    int bin = (int) (r/radBinWidth);
		    gr12[bin]++;
	        }
	    }
	}

	for(int i = 0; i < N2; i++){ // histogram for g22(r)
	    for(int j = 0; j < N2; j++){
		if(i != j){
		    double rx = PBC.separation(Math.abs(x2[i]-x2[j]), side);
		    double ry = PBC.separation(Math.abs(y2[i]-y2[j]), side);
		    double rz = PBC.separation(Math.abs(z2[i]-z2[j]), side);
		    double r2=rx*rx+ry*ry+rz*rz;
		    if(r2 < maxR*maxR){
			double r = Math.sqrt(r2);
			int bin = (int) (r/radBinWidth);
			gr22[bin]++;
		    }
		}
	    }
	}

    }
    
    public void normalize(){
	double rho = (double) N1/(side*side*side); // number density of species 1
	double norm = 4*Math.PI*rho*Math.pow(radBinWidth,3)*N1*count; // normalization
	
	for(int i = 0; i < gr11.length; i++){
	    double r = radBinWidth*i; // approximate value of radius in bin i (rounded to bottom of bin)
	    if(r <= maxR){
		gr11[i] = gr11[i] / (norm*(i+1)*(i+1));
		gr12[i] = gr12[i] / (norm*(i+1)*(i+1));
		gr22[i] = gr22[i] / (norm*(i+1)*(i+1));
	    }
	}
    }

    public void writeRDF(){ // write normalized g(r) data to a space delimited file

	normalize();

	try{
	    File radFile = new File("radialdistribution" + extension + ".txt");
	    
	    // if file doesn't exist, create it
	    if (!radFile.exists()) {
		radFile.createNewFile();
	    }
	    
	    FileWriter fwrite = new FileWriter(radFile.getAbsoluteFile());
	    BufferedWriter bwrite = new BufferedWriter(fwrite);
	    for (int i=0;i<gr11.length;i++){
		bwrite.write(i + " " + gr11[i] + " " + gr12[i] + " " + gr22[i]);
		bwrite.newLine();
	    }
	    
	    bwrite.close();
	}
	
	catch (IOException e) {
	    e.printStackTrace();
	}
    }
    
}
