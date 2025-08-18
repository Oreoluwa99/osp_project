package org.opensourcephysics.sip.HertzMixture;

import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.opensourcephysics.numerics.PBC;

/**
 * SSFmixture computes the static structure factor from a Monte Carlo simulation of a microgel mixture
 * 
 * @author Alan Denton and Matt Urich
 * @version 1.0 24-01-2015
 * 
 */

public class SSFmixture {
    private int N1, N2;
    private double [] x1, y1, z1, x2, y2, z2; // particle coordinates
    private double [] sk, sk11, sk12, sk22; // partial static structure factors
    private double side, d;
    private double deltaK; // wavenumber increment
    private int bins; // number of bins
    private double rnn; // nearest-neighbor distance
    private int count = 0;
    private String extension; // for writing the data file

    public SSFmixture(double x1[], double y1[], double z1[], double x2[], double y2[], double z2[], double side, double d, double deltaK, String extension){
    
	N1 = x1.length; // number of particles of type 1
	N2 = x2.length; // number of particles of type 2
	this.x1 = x1;
	this.y1 = y1;
	this.z1 = z1;
	this.x2 = x2;
	this.y2 = y2;
	this.z2 = z2;
	this.side = side;
	this.d = d;
	this.deltaK = deltaK;
	this.extension = extension;

	rnn = d/Math.sqrt(2); // only for an fcc lattice
	bins = (int) (6*2*Math.PI/rnn/deltaK) + 1; // 6 multiples of 2*pi/(nearest-neighbor distance)
	sk = new double[bins];
	sk11 = new double[bins];
	sk12 = new double[bins];
	sk22 = new double[bins];

    }

    public void update(){

	double rx, ry, rz, r;

	count++;

	// ignore short-k (long-wavelength) limit, where finite-size corrections would be needed 
	for (int k=40;k<bins;k++){ 

	    for (int i=0;i<N1;i++){ // S11
		for (int j=i+1;j<N1;j++){
		    rx = x1[i]-x1[j];
		    ry = y1[i]-y1[j];
		    rz = z1[i]-z1[j];
		    r = Math.sqrt(rx*rx+ry*ry+rz*rz);
		    sk11[k] += Math.sin(k*deltaK*r)/(k*deltaK*r);
		}
	    }

	    for (int i=0;i<N1;i++){ // S12 (same as S21)
		for (int j=0;j<N2;j++){
		    rx = x1[i]-x2[j];
		    ry = y1[i]-y2[j];
		    rz = z1[i]-z2[j];
		    r = Math.sqrt(rx*rx+ry*ry+rz*rz);
		    sk12[k] += Math.sin(k*deltaK*r)/(k*deltaK*r);
		}
	    }

	    for (int i=0;i<N2;i++){ // S22
		for (int j=i+1;j<N2;j++){
		    rx = x2[i]-x2[j];
		    ry = y2[i]-y2[j];
		    rz = z2[i]-z2[j];
		    r = Math.sqrt(rx*rx+ry*ry+rz*rz);
		    sk22[k] += Math.sin(k*deltaK*r)/(k*deltaK*r);
		}
	    }

	}
    }

    public void normalize(){
	for (int k=0;k<bins;k++){
            // total structure factor
            sk[k] = (2*sk11[k] + sk12[k] + 2*sk22[k])/(N1+N2)/count + 1; 
            // partial structure factors
	    sk11[k] = 2*sk11[k]/(N1+N2)/count + (double)N1/(double)(N1+N2);
	    sk12[k] = sk12[k]/(N1+N2)/count;
	    sk22[k] = 2*sk22[k]/(N1+N2)/count + (double)N2/(double)(N1+N2);
	}
    }

    public void writeSSF(){
	normalize();

	try{
	    File ssfFile = new File("ssf" + extension + ".txt");
	    
	    // if file doesn't exist, create it
	    if (!ssfFile.exists()) {
		ssfFile.createNewFile();
	    }
	    
	    FileWriter fwrite = new FileWriter(ssfFile.getAbsoluteFile());
	    BufferedWriter bwrite = new BufferedWriter(fwrite);
	    for (int k=0;k<sk.length;k++){
		bwrite.write(k*deltaK + "  " + sk[k] + "  " +  sk11[k] + "  " + sk12[k] + "  " + sk22[k]);
		bwrite.newLine();
	    }
	    
	    bwrite.close();
	}
	
	catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
