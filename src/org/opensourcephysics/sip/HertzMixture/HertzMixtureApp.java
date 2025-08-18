package org.opensourcephysics.sip.HertzMixture;

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
 * HertzMixtureApp performs a Monte Carlo simulation of a binary mixture of microgels interacting 
 * via the Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @author Alan Denton 
 * @version 1.2 07.09.2021
 * 
 */


public class HertzMixtureApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE,WRITE_RADIAL,WRITE_ALL;};
    HertzMixture particles = new HertzMixture();
    PlotFrame energyData = new PlotFrame("MC steps", "<E_pair>/N", "Mean pair energy");
    //PlotFrame pressureData = new PlotFrame("MC steps", "PV/NkT", "Mean Pressure");
    PlotFrame sizeData = new PlotFrame("MC steps", "particle radius [nm]", "Mean Particle Size");
    //PlotFrame sizeData = new PlotFrame("MC steps", "swelling ratio", "Mean Swelling Ratios");
    Display3DFrame display3d = new Display3DFrame("Simulation Snapshot");
    int dataPoints;
    int maxDataPoints;
    ElementSphere nanoSphere1[], nanoSphere2[];
    boolean added = false;
    boolean structure;
    RDFmixture rdf;
    SSFmixture ssf;

	
         /**
	 * Initializes the model.
	 */
	public void initialize() {
		added = false;
		particles.N1 = control.getInt("N1");
		particles.N2 = control.getInt("N2");
		String configuration = control.getString("Initial configuration");
		particles.dryRadius1 = control.getDouble("Dry Radius 1 [nm]");
		particles.dryRadius2 = control.getDouble("Dry Radius 2 [nm]");
		particles.xLinkFrac1 = control.getDouble("x-link Fraction 1");
		particles.xLinkDensityRatio = control.getDouble("x-link Density Ratio (2/1)");
		particles.dryVolFrac1 = control.getDouble("Dry Volume Fraction 1");
		particles.Young = control.getDouble("Young's Calibration"); // 10-1000
		particles.chi1 = control.getDouble("chi1");
		particles.chi2 = control.getDouble("chi2");
                particles.maxSwellingRatio = control.getDouble("Maximum Swelling Ratio");
		particles.tolerance = control.getDouble("Displacement Tolerance");
		particles.atolerance = control.getDouble("Radius Change Tolerance");
		particles.delay = control.getDouble("Interval");
		particles.stop = control.getInt("Stop");
		particles.snapshotInterval = control.getInt("Snapshot Interval");
		particles.sizeBinWidth = control.getDouble("Size Bin Width");
		particles.grBinWidth = control.getDouble("g(r) Bin Width");
		particles.deltaK = control.getDouble("Delta k");
		particles.fileExtension = control.getString("File Extension");
		structure = control.getBoolean("Calculate Structure");
		particles.initialize(configuration);

                System.out.println("nMon1: " + particles.nMon1);
                System.out.println("nChains1: " + particles.nChains1);
                System.out.println("nMon2: " + particles.nMon2);
                System.out.println("nChains2: " + particles.nChains2);
                System.out.println("reservoirSwellingRatio1: " + particles.reservoirSR1);
                System.out.println("reservoirSwellingRatio2: " + particles.reservoirSR2);
                System.out.println("reservoirVolFrac: " + particles.reservoirVolFrac);

		if(display3d != null) display3d.dispose(); // closes old simulation frame if present
		display3d = new Display3DFrame("3D Frame");
		display3d.setPreferredMinMax(0, particles.side, 0, particles.side, 0, particles.side);
		display3d.setSquareAspect(true);
		//energyData.setPreferredMinMax(0, 1000, -10, 10);
		//pressureData.setPreferredMinMax(0, 1000, 0, 2);
		
		// add simple3d.Element particles to the arrays 
		if (!added) { // particles can be added only once
			nanoSphere1 = new ElementSphere[particles.N1];
			nanoSphere2 = new ElementSphere[particles.N2];

			for (int i = 0; i < particles.N1; i++) {
				nanoSphere1[i] = new ElementSphere();
				display3d.addElement(nanoSphere1[i]);
			}
			for (int i = 0; i < particles.N2; i++) {
				nanoSphere2[i] = new ElementSphere();
				display3d.addElement(nanoSphere2[i]);
			}
			added = true; 
		}
		
		// initialize visualization elements for nanoparticles
		for (int i = 0; i < particles.N1; i++) {
			nanoSphere1[i].setSizeXYZ(particles.a1[i]*2, particles.a1[i]*2, particles.a1[i]*2);
			//nanoSphere1[i].getStyle().setFillColor(Color.GREEN);
			nanoSphere1[i].getStyle().setFillColor(Color.RED); // choose favorite color
			//nanoSphere1[i].getStyle().setFillColor(new Color(92, 146, 237, 100));  // light blue with half transparency
			nanoSphere1[i].setXYZ(particles.x1[i], particles.y1[i], particles.z1[i]);
		}
		for (int i = 0; i < particles.N2; i++) {
			nanoSphere2[i].setSizeXYZ(particles.a2[i]*particles.dryRadiusRatio*2, particles.a2[i]*particles.dryRadiusRatio*2, particles.a2[i]*particles.dryRadiusRatio*2);
			//nanoSphere2[i].getStyle().setFillColor(Color.YELLOW);
			//nanoSphere2[i].getStyle().setFillColor(Color.GREEN);
			nanoSphere2[i].getStyle().setFillColor(Color.BLUE);
			nanoSphere2[i].setXYZ(particles.x2[i], particles.y2[i], particles.z2[i]);
		}

		if (structure){
		    rdf = new RDFmixture(particles.x1,particles.y1,particles.z1,particles.x2,particles.y2,particles.z2,particles.side,particles.grBinWidth,particles.fileExtension);
		    ssf = new SSFmixture(particles.x1,particles.y1,particles.z1,particles.x2,particles.y2,particles.z2,particles.side,particles.d,particles.deltaK,particles.fileExtension);
		}

	}


	
	/**
	 * Does a simulation step.
	 */
	public void doStep() { // logical step in HertzMixture class	    	    
		particles.step();

		if (particles.steps > particles.delay/10.){
                    particles.scale = 1.;
                }

		if (particles.steps <= particles.stop){
		    if (particles.steps > particles.delay){
			if ((particles.steps-particles.delay) % particles.snapshotInterval == 0){
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
		
		energyData.append(0, particles.steps, particles.meanPairEnergy());
		//pressureData.append(1, particles.steps, particles.meanPressure());
/* Plot particle radii: */
		sizeData.append(0, particles.steps, particles.meanSwellingRatio1()*particles.dryRadius1);
		sizeData.append(2, particles.steps, particles.meanSwellingRatio2()*particles.dryRadius2);
		//sizeData.append(0, particles.steps, particles.meanSwellingRatio1());
		//sizeData.append(2, particles.steps, particles.meanSwellingRatio2());
		
		// update steps
		display3d.setMessage("Number of steps: " + particles.steps);

		// visualization updates
		if(control.getBoolean("Visualization on")){
			for (int i = 0; i < particles.N1; i++) {
			    nanoSphere1[i].setSizeXYZ(particles.a1[i]*2, particles.a1[i]*2, particles.a1[i]*2);
			    //nanoSphere1[i].getStyle().setFillColor(Color.GREEN);
			    nanoSphere1[i].setXYZ(particles.x1[i], particles.y1[i], particles.z1[i]);
			}
			for (int i = 0; i < particles.N2; i++) {
			    nanoSphere2[i].setSizeXYZ(particles.a2[i]*particles.dryRadiusRatio*2, particles.a2[i]*particles.dryRadiusRatio*2, particles.a2[i]*particles.dryRadiusRatio*2);
			    //nanoSphere2[i].getStyle().setFillColor(Color.YELLOW);
			    nanoSphere2[i].setXYZ(particles.x2[i], particles.y2[i], particles.z2[i]);
			}
		}
	}

	/**
	 * Resets the model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		control.setValue("N1", 250);
		control.setValue("N2", 250);
		control.setValue("Initial configuration", "random-FCC");
		//control.setValue("Initial configuration", "disordered-FCC");
		//control.setValue("Initial configuration", "FCC");
		//control.setValue("Initial configuration", "SC");
		control.setValue("Dry Radius 1 [nm]", 40);
		control.setValue("Dry Radius 2 [nm]", 50);
		control.setValue("x-link Fraction 1", 0.001);
		control.setValue("x-link Density Ratio (2/1)", 1);
		control.setValue("Dry Volume Fraction 1", 0.01);
		control.setValue("Young's Calibration", 10); // 10-1000
		control.setValue("chi1", 0);
		control.setValue("chi2", 0);
		control.setValue("Maximum Swelling Ratio", 6);
		control.setValue("Displacement Tolerance", 0.1);
		control.setValue("Radius Change Tolerance", 0.05);
		control.setValue("Interval", 500);
		control.setValue("Stop", 1500);
		control.setValue("Snapshot Interval", 100);
		//control.setValue("Interval", 50000);
		//control.setValue("Stop", 150000);
		//control.setValue("Snapshot Interval", 100);
		control.setValue("Size Bin Width", .001);
		control.setValue("g(r) Bin Width", .005);
		control.setValue("Delta k", .005);
		control.setValue("File Extension", "1");
		control.setValue("Calculate Structure", false);
		control.setAdjustableValue("Visualization on", true);
		initialize();
	}


	public void stop() {
		
    		control.println("Number of MC steps = "+particles.steps);
    		control.println("<E/N> = "+decimalFormat.format(particles.meanEnergy()));
		control.println("<E_pair>/N = "+decimalFormat.format(particles.meanPairEnergy()));
                control.println("<F>/N = "+decimalFormat.format(particles.meanFreeEnergy()));
                //control.println("PV/NkT = "+decimalFormat.format(particles.meanPressure()));
	}

	public void writeData(){
	    
	    for (int i=0;i<(particles.maxSwellingRatio/particles.sizeBinWidth);i++){
		particles.sizeDist1[i] = particles.sizeDist1[i]/((particles.stop-particles.delay)/particles.snapshotInterval); //normalize size distribution
		particles.sizeDist2[i] = particles.sizeDist2[i]/((particles.stop-particles.delay)/particles.snapshotInterval); //normalize size distribution
	    }
	    
	    
	    particles.volFrac = particles.volFrac/((particles.stop-particles.delay)/particles.snapshotInterval); // fixing the volume fraction
	    
	    // first write test information to a file
	    try{
		File testInfo = new File("systemInfo" + particles.fileExtension + ".txt");
		
		// if file doesn't exist, create it
		if (!testInfo.exists()){
		    testInfo.createNewFile();
		}
		
		FileWriter fw = new FileWriter(testInfo.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("N1 : " + particles.N1);
		bw.newLine();
		bw.write("N2 : " + particles.N2);
		bw.newLine();
		bw.write("Dry Radius 1 : " + particles.dryRadius1);
		bw.newLine();
		bw.write("Dry Radius 2 : " + particles.dryRadius2);
		bw.newLine();
		bw.write("box length : " + particles.side);
		bw.newLine();
		bw.write("N monomers 1 : " + particles.nMon1);
		bw.newLine();
		bw.write("N monomers 2 : " + particles.nMon2);
		bw.newLine();
		bw.write("N chains 1 : " + particles.nChains1);
		bw.newLine();
		bw.write("N chains 2 : " + particles.nChains2);
		bw.newLine();
		bw.write("Flory interaction parameter 1 : " + particles.chi1);
		bw.newLine();
		bw.write("Flory interaction parameter 2 : " + particles.chi2);
		bw.newLine();
		bw.write("Young's calibration: " + particles.Young);
		bw.newLine();
		bw.write("B11 : " + particles.B11);
		bw.newLine();
		bw.write("B22 : " + particles.B22);
		bw.newLine();
		bw.write("B12 : " + particles.B12);
		bw.newLine();
		bw.write("dry radius ratio (2/1) : " + particles.dryRadiusRatio);
		bw.newLine();
		bw.write("x-link fraction 1 : " + particles.xLinkFrac1);
		bw.newLine();
		bw.write("x-link density ratio (2/1) : " + particles.xLinkDensityRatio);
		bw.newLine();
		bw.write("reservoir swelling ratio 1 : " + particles.reservoirSR1);
		bw.newLine();
		bw.write("reservoir swelling ratio 2 : " + particles.reservoirSR2);
		bw.newLine();
		bw.write("mean particle radius 1 : " + particles.meanRadius1());
		bw.newLine();
		bw.write("mean particle radius 2 : " + particles.meanRadius2());
		bw.newLine();
		bw.write("dry volume fraction 1 : " + particles.dryVolFrac1);
		bw.newLine();
		bw.write("dry volume fraction 2 : " + particles.dryVolFrac2);
		bw.newLine();
		bw.write("reservoir volume fraction : " + particles.reservoirVolFrac);
		bw.newLine();
		bw.write("actual volume fraction : " + particles.volFrac);
		bw.newLine();
		bw.write("steps : " + particles.steps);
		bw.newLine();
		bw.write("delay : " + particles.delay);
		bw.newLine();
		bw.write("snapshot delay : " + particles.snapshotInterval);
		bw.newLine();
		bw.write("tolerance " + particles.tolerance);
		bw.newLine();
		bw.write("radius tolerance : " + particles.atolerance);
		bw.newLine();
		bw.write("size bin width : " + particles.sizeBinWidth);
		bw.newLine();
		bw.write("rad bin width: " + particles.grBinWidth);
		bw.newLine();
                bw.write("Mean pair energy <E_pair>/N [kT]: " + particles.meanPairEnergy());
                bw.newLine();
                //bw.write("Mean pressure PV/NkT: " + particles.meanPressure());
                //bw.newLine();
                bw.write("Mean energy <E>/N [kT]: " + particles.meanEnergy());
                bw.newLine();
                bw.write("Mean free energy <F>/N [kT]: " + particles.meanFreeEnergy());
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
		    bwrite.write(i + " " + particles.sizeDist1[i] + " " + particles.sizeDist2[i]);
		    bwrite.newLine();
		}
		
		bwrite.close();
		
	    }
	    
	    catch (IOException e) {
		e.printStackTrace();
	    } // write size distributions and radial distribution functions to a file

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
		SimulationControl control = SimulationControl.createApp(new HertzMixtureApp());
	}
}
