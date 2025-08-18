package org.opensourcephysics.sip.Hertz;

import java.awt.Color;
import java.text.DecimalFormat;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.PlotFrame;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;


/**
 * Performs a Monte Carlo simulation of volume-charged ionic microgels interacting via a
 * Hertz-Yukawa effective pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @author Based on MC code ch15/LJParticles, modified for microgels by Matt Urich and Alan Denton.
 * @version 2.2 05.05.2020
 * 
 */


public class VolumeChargedMicrogelsApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE,WRITE_RADIAL,WRITE_ALL;};
    VolumeChargedMicrogels particles = new VolumeChargedMicrogels();
    PlotFrame energyData = new PlotFrame("steps", "E/N", "Mean energy");
    PlotFrame pressureData = new PlotFrame("steps", "PV/NkT", "Mean pressure");
    PlotFrame sizeData = new PlotFrame("steps", "alpha", "Mean swelling ratio");
    Display3DFrame display3d = new Display3DFrame("3D Frame");
    int dataPoints;
    int maxDataPoints;
    ElementSphere nanoSphere[];
    boolean added = false;
    boolean structure;
    DataFile [] dataFiles;
    RDF rdf;
    SSF ssf;

         /**
	 * Initializes the model.
	 */
	public void initialize() {
		added = false;
		particles.N = control.getInt("N");
		String configuration = control.getString("Initial Configuration");
                particles.dryR = control.getDouble("Dry Radius [nm]");
		particles.dryVolFrac= control.getDouble("Dry Volume Fraction");
		particles.nMon = control.getDouble("nMon");
		particles.nChains = control.getDouble("nChains");
		particles.chi = control.getDouble("chi");
		particles.B = control.getDouble("B Parameter");
		particles.Z = control.getDouble("Z");
		particles.csalt = control.getDouble("[Salt/microM]");
		particles.lBjerrum= control.getDouble("Bjerrum length [nm]");
		particles.tolerance = control.getDouble("Displacement Tolerance");
		particles.atolerance = control.getDouble("Radius Change Tolerance");
		particles.delay = control.getDouble("Delay");
		particles.snapshotDelay = control.getInt("Snapshot Delay");
		particles.stop = control.getInt("Stop");
		particles.maxSwellingRatio = control.getDouble("Maximum Swelling Ratio");
		particles.sizeBinWidth = control.getDouble("Size Bin Width");
		particles.radBinWidth = control.getDouble("g(r) Bin Width");
		particles.deltaK = control.getDouble("Delta k");
		particles.fileExtension = control.getString("File Extension");
		structure = control.getBoolean("Calculate Structure");

		particles.initialize(configuration);

                System.out.println("reservoirSwellingRatio: " + particles.reservoirR);
                System.out.println("reservoirVolFrac: " + particles.reservoirVolFrac);
		
		if(display3d != null) display3d.dispose(); // closes old simulation frame if present
		display3d = new Display3DFrame("3D Frame");
		display3d.setPreferredMinMax(0, particles.side, 0, particles.side, 0, particles.side);
		display3d.setSquareAspect(true);
		//energyData.setPreferredMinMax(0, 1000, -10, 10);
		//pressureData.setPreferredMinMax(0, 1000, 0, 2);
		
		// add simple3d.Element particles to the arrays 
		if (!added) { // particles can be added only once
			nanoSphere = new ElementSphere[particles.N];

			for (int i = 0; i < particles.N; i++) {
				nanoSphere[i] = new ElementSphere();
				display3d.addElement(nanoSphere[i]);
			}
			added = true; 
		}
		
		// initialize visualization elements for nanoparticles
		for (int i = 0; i < particles.N; i++) {
			nanoSphere[i].setSizeXYZ(particles.a[i]*2, particles.a[i]*2, particles.a[i]*2);
			nanoSphere[i].getStyle().setFillColor(Color.RED);
			nanoSphere[i].setXYZ(particles.x[i], particles.y[i], particles.z[i]);
		}
		if (structure){
		    rdf = new RDF(particles.x,particles.y,particles.z,particles.side,particles.radBinWidth,particles.fileExtension);
		    ssf = new SSF(particles.x,particles.y,particles.z,particles.side,particles.d,particles.deltaK,particles.fileExtension);
		}
	}


	
	/**
	 * Does a simulation step.
	 */
	public void doStep() {
		// logical step in the VolumeChargedMicrogels class	    	    
		particles.step();

		if (particles.steps<=particles.stop){
		    if (particles.steps>particles.delay){
			if ((particles.steps-particles.delay) % particles.snapshotDelay == 0){
			    particles.sizeDistributions();
			    if (structure){
				rdf.update();
				ssf.update();
			    }
			    particles.calculateVolumeFraction();
			}
		    }
		}
		    
		if (particles.steps == particles.stop){
		    writeData();
		}
		
		energyData.append(0, particles.steps, particles.meanEnergy());
		//System.out.println(particles.meanEnergy()); // test
		pressureData.append(1, particles.steps, particles.meanPressure());
		sizeData.append(2, particles.steps, particles.meanSwellingRatio());
		
		// update steps
		display3d.setMessage("Number of steps: " + particles.steps);

		// visualization updates
		if(control.getBoolean("Visualization On")){
			for (int i = 0; i < particles.N; i++) {
			    nanoSphere[i].setSizeXYZ(particles.a[i]*2, particles.a[i]*2, particles.a[i]*2);
			    nanoSphere[i].getStyle().setFillColor(Color.RED);
			    nanoSphere[i].setXYZ(particles.x[i], particles.y[i], particles.z[i]);
			}
		}
	}

	/**
	 * Resets the model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		control.setValue("N", 108);
		control.setValue("Initial Configuration", "fcc");
                control.setValue("Dry Radius [nm]", 10);
		control.setValue("Dry Volume Fraction", 0.01);
		control.setValue("nMon", 200000);
		control.setValue("nChains", 100);
		control.setValue("chi", 0.5);
		control.setValue("B Parameter", 15000);
		control.setValue("Z", 1000);
		control.setValue("[Salt/microM]", 0);
		control.setValue("Bjerrum length [nm]", 0.714); // Bjerrum length 
		control.setValue("Displacement Tolerance", 0.1);
		control.setValue("Radius Change Tolerance", 0.05);
		control.setValue("Delay", 5000);//changed
		control.setValue("Snapshot Delay", 100);
		control.setValue("Stop", 15000);//changed
		control.setValue("Maximum Swelling Ratio", 4); // upper range of size distribution
		control.setValue("Size Bin Width", .001);
		control.setValue("g(r) Bin Width", .005);
		control.setValue("Delta k", .005);
		control.setValue("File Extension", "1");
		control.setValue("Calculate Structure", false);
		control.setAdjustableValue("Visualization On", true);
		initialize();
	}


	public void stop() {
		
    		control.println("Number of MC steps = "+particles.steps);
    		control.println("<E/N> = "+decimalFormat.format(particles.meanEnergy()));
    		control.println("<PV/NkT> = "+decimalFormat.format(particles.meanPressure()));
	}

	public void writeData(){
	    
// normalize the size distribution
	    for (int i=0;i<(particles.maxSwellingRatio/particles.sizeBinWidth);i++){
		particles.sizeDist[i] = particles.sizeDist[i]/((particles.stop-particles.delay)/particles.snapshotDelay); 
	    }
	    
	    
// normalize the volume fraction
	    particles.volFrac = particles.volFrac/((particles.stop-particles.delay)/particles.snapshotDelay); 
	    
	    // first, write the test information to a file
	    try{
		File testInfo = new File("microgeltestinfo" + particles.fileExtension + ".txt");
		
		// if file doesn't exist, create it
		if (!testInfo.exists()){
		    testInfo.createNewFile();
		}
		
		FileWriter fw = new FileWriter(testInfo.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("N: " + particles.N);
		bw.newLine();
		bw.write("n monomers: " + particles.nMon);
		bw.newLine();
		bw.write("n chains: " + particles.nChains);
		bw.newLine();
		bw.write("interaction parameter: " + particles.chi);
		bw.newLine();
		bw.write("B: " + particles.B);
		bw.newLine();
		bw.write("Z: " + particles.Z);
		bw.newLine();
		bw.write("dry radius [nm]: " + particles.dryR);
		bw.newLine();
		bw.write("gamma: " + particles.gamma);
		bw.newLine();
		bw.write("csalt [microM]: " + particles.csalt);
		bw.newLine();
		bw.write("kappa * (dry radius): " + particles.k);
		bw.newLine();
		bw.write("reservoir swelling ratio: " + particles.reservoirR);
		bw.newLine();
		bw.write("mean swelling ratio: " + particles.meanRadius());
		bw.newLine();
		bw.write("dry volume fraction: " + particles.dryVolFrac);
		bw.newLine();
		bw.write("reservoir volume fraction: " + particles.reservoirVolFrac);
		bw.newLine();
		bw.write("actual volume fraction: " + particles.volFrac);
		bw.newLine();
		bw.write("box length / dry radius: " + particles.side);
		bw.newLine();
		bw.write("steps: " + particles.steps);
		bw.newLine();
		bw.write("delay: " + particles.delay);
		bw.newLine();
		bw.write("snapshot delay: " + particles.snapshotDelay);
		bw.newLine();
		bw.write("tolerance: " + particles.tolerance);
		bw.newLine();
		bw.write("radius tolerance: " + particles.atolerance);
		bw.newLine();
		bw.write("size bin width: " + particles.sizeBinWidth);
		bw.newLine();
		bw.write("rad bin width: " + particles.radBinWidth);
		bw.newLine();
		bw.write("mean energy: " + particles.meanEnergy());
		bw.newLine();
		bw.write("mean pressure: " + particles.meanPressure());
		
		bw.close();
	    }
	    
	    catch (IOException e) {
		e.printStackTrace();
	    }
	    
	    // now write size distributions
	    
	    try {
		File sizeFile = new File("microgelsizes" + particles.fileExtension + ".txt");
		
		// if file doesn't exist, create it
		if (!sizeFile.exists()) {
		    sizeFile.createNewFile();
		}
		
		FileWriter fwrite = new FileWriter(sizeFile.getAbsoluteFile());
		BufferedWriter bwrite = new BufferedWriter(fwrite);
		for (int i=0;i<(particles.maxSwellingRatio/particles.sizeBinWidth);i++){
		    bwrite.write(i + " " + particles.sizeDist[i]);
		    bwrite.newLine();
		}
		
		bwrite.close();
		
	    }
	    
	    catch (IOException e) {
		e.printStackTrace();
	    } // writes size distributions and radial distribution function to a file

	    if (structure){
		rdf.writeRDF();
		ssf.writeSSF();
	    }
	    
	}
	
	/**
	 * Start the Java application.
	 * 
	 * @param args
	 *            command line parameters
	 */
	public static void main(String[] args) { // set up animation control
		@SuppressWarnings("unused")
		SimulationControl control = SimulationControl.createApp(new VolumeChargedMicrogelsApp());
	}
}
