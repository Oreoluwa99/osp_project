package org.opensourcephysics.sip.Hertz;

import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.opensourcephysics.numerics.PBC;

/**
 * SSF calculates the static structure factor from a Monte Carlo simulation
 * 
 * @authors Alan Denton and Matt Urich
 * @version 1.0 24-01-2015
 * 
 */

public class SSF {
    private int N;
    private double [] x, y, z;
    private double [] sk;
    private double side, d;
    private double deltaK;
    private int bins;
    private double rnn;
    private int count = 0;
    private String extension; //this is just for writing the file

    public SSF(double x[], double y[], double z[], double side, double d, double deltaK, String extension){
	N = x.length;
	this.x = x;
	this.y = y;
	this.z = z;
	this.side = side;
	this.d = d;
	this.deltaK = deltaK;
	this.extension = extension;

	rnn = d/Math.sqrt(2); // this is only for an FCC lattice
	bins = (int) (3*2*Math.PI/rnn/deltaK) + 1; // 6 multiples of 2*pi/nearest neighbor distance
	sk = new double[bins];

    }

    public void update(){

	double rx, ry, rz, r;

	count++;

	for (int k=100;k<bins;k++){ // k starts at 100 to avoid long wavelength limit 
	    for (int i=0;i<N;i++){
		for (int j=i+1;j<N;j++){

		    rx = x[i]-x[j];
		    ry = y[i]-y[j];
		    rz = z[i]-z[j];

		    r = Math.sqrt(rx*rx+ry*ry+rz*rz);

		    sk[k] += Math.sin(k*deltaK*r)/(k*deltaK*r);
		}
	    }
	}
    }

    public void normalize(){
	for (int k=0;k<bins;k++){
	    sk[k] = 2*sk[k]/(N*count)+1;
	}
    }

    public void writeSSF(){
	normalize();

	try{
	    File ssfFile = new File("data/ssf_and_rdf_data/x_link_3.4E-5/ssf" + extension + ".txt");
	    
	    if (!ssfFile.exists()) { // if file doesn't exist, create it
		ssfFile.createNewFile();
	    }
	    
	    FileWriter fwrite = new FileWriter(ssfFile.getAbsoluteFile());
	    BufferedWriter bwrite = new BufferedWriter(fwrite);
	    for (int k=0;k<sk.length;k++){
		bwrite.write((k*deltaK) + " " + sk[k]);
		bwrite.newLine();
	    }
	    
	    bwrite.close();
	}
	
	catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
    
