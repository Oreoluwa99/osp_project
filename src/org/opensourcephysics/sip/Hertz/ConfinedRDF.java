package org.opensourcephysics.sip.Hertz;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.opensourcephysics.numerics.PBC;

/**
 * ConfinedRDF calculates the radial distribution function from a Monte Carlo simulation
 * 
 * @authors Oreoluwa Alade and Alan Denton
 * @version 1.0 31-10-2026
 * 
 */

public class ConfinedRDF {
    private int N;
    private double [] x, y, z;
    private double [] gr;
    private double side;
    private double radBinWidth;
    private int bins;
    private double maxR;
    private int count = 0;
    private String extension; // for writing the file
	private double Lx; // box length in x direction
	private double young; // Hertz calibration, used in the folder name
	private int xBins = 200; // bins for the centre distribution across the slit
	private double xBinWidth;
	private double [] xHist; // counts of particle centres at each x
    
	public ConfinedRDF(double x[], double y[], double z[], double side, double Lx, double radBinWidth, String extension, double young){
		N = x.length;
		this.x = x;
		this.y = y;
		this.z = z;
		this.Lx = Lx; // box length in x direction
		this.side = side;
		this.young = young;
		this.radBinWidth = radBinWidth;
		this.extension = extension; 
		maxR = side/2.0;
		bins = (int) Math.ceil(maxR / radBinWidth);
		gr = new double [bins];
		xBinWidth = Lx / xBins;
		xHist = new double [xBins];
    }
    
    public void update(){
	count++;
	for(int i = 0; i < N; i++){
		int b = (int) (x[i] / xBinWidth); // record where centre i sits across the slit
		if (b >= xBins) b = xBins - 1;
		if (b < 0) b = 0;
		xHist[b]++;

	    for(int j = 0; j < N; j++){
		if(i != j){
			double rx = Math.abs(x[i] - x[j]);
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

	// probability that two independent centres are k x-bins apart (both signs)
	private double[] separationWeights() {
		double total = 0;
		for (int k = 0; k < xBins; k++) total += xHist[k];

		double[] p = new double[xBins];
		for (int k = 0; k < xBins; k++) p[k] = xHist[k] / total;

		double[] w = new double[xBins];
		for (int k = 0; k < xBins; k++) {
			double sum = 0;
			for (int m = 0; m + k < xBins; m++) sum += p[m] * p[m + k];
			w[k] = (k == 0) ? sum : 2.0 * sum;
		}
		return w;
	}

	// fraction of randomly placed partners that land in the shell [rLow, rHigh]
	private double shellFraction(double rLow, double rHigh, double[] w) {
		double fraction = 0;
		for (int k = 0; k < xBins; k++) {
			double xi = k * xBinWidth;
			double outer = Math.max(rHigh*rHigh - xi*xi, 0.0);
			double inner = Math.max(rLow*rLow - xi*xi, 0.0);
			fraction += w[k] * Math.PI * (outer - inner);
		}
		return fraction / (side * side);
	}

	private double normalizedBin(int i, double[] w) {
        if (N < 2 || count == 0) {
            throw new IllegalStateException("No particle pairs were sampled.");
        }

        double rLow = i * radBinWidth;
        double rHigh = Math.min((i + 1) * radBinWidth, maxR);

        double expectedPairs =
            count * (double) N * (N - 1)
            * shellFraction(rLow, rHigh, w);

        if (expectedPairs == 0) return 0.0;
        return gr[i] / expectedPairs;
    }

	public void writeRDF() {
		File directory = new File(
			"data/confined/particle_number_study/cubic_FCC_start/Yc" + young + "/rdf/N" + N
		);

		if (!directory.exists() && !directory.mkdirs()) {
			throw new IllegalStateException(
				"Could not create RDF directory: " + directory
			);
		}

		File radFile = new File(
			directory,
			"confined_rdf_N" + N
			+ "_Lx" + Lx
			+ "_run" + extension + ".txt"
		);

		if (radFile.exists()) {
			throw new IllegalStateException(
				"RDF file already exists: " + radFile
			);
		}

		double[] w = separationWeights();

		try (BufferedWriter bwrite =
				new BufferedWriter(new FileWriter(radFile))) {

			bwrite.write("# Lx=" + Lx + " N=" + N
						+ " samples=" + count + " side=" + side);
			bwrite.newLine();
			bwrite.write("# r g(r)");
			bwrite.newLine();

			for (int i = 0; i < gr.length; i++) {
				double rLow = i * radBinWidth;
				double rHigh = Math.min((i + 1) * radBinWidth, maxR);
				double rMid = (rLow + rHigh) / 2.0;

				bwrite.write(rMid + " " + normalizedBin(i, w));
				bwrite.newLine();
			}

		} catch (IOException e) {
			throw new RuntimeException(
				"Could not write confined RDF: " + radFile, e
			);
		}
	}
}