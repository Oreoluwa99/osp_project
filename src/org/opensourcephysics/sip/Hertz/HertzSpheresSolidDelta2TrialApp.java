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

// imported for the purpose of creating an array to store the values of the dryVolFrac and free energies
import java.util.List;
import java.util.ResourceBundle.Control;
import java.util.ArrayList;


/**
 * HertzSpheresApp performs a Monte Carlo simulation of microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @authors Alan Denton and Matt Urich
 * @version 1.2 25-05-2021
 * 
 */

 public class HertzSpheresSolidDelta2TrialApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
	HertzSpheresCopyDelta2 particles = new HertzSpheresCopyDelta2();
//	HertzSpheresSolidPhaseEdited particles = new HertzSpheresSolidPhaseEdited();
//	HertzSpheresSolid particles = new HertzSpheresSolid();
//	HertzSpheresSolidPhaseTrialCode particles = new HertzSpheresSolidPhaseTrialCode();
    PlotFrame energyData = new PlotFrame("MC steps", "<E_pair>/N", "Mean pair energy per particle");
    PlotFrame pressureData = new PlotFrame("MC steps", "PV/NkT", "Mean pressure");
    PlotFrame sizeData = new PlotFrame("MC steps", "alpha", "Mean swelling ratio");
    Display3DFrame display3d = new Display3DFrame("Simulation animation");
    int dataPoints;
    int maxDataPoints;
    ElementSphere nanoSphere[];
    boolean added = false;
    boolean structure;
    RDF rdf;
    SSF ssf;
	double dryVolFracStart, dryVolFracMax, dryVolFrac;
	double totalSpringEnergy, lambda, totalVol, totalSumOfEnergies, initialEnergy;
	boolean setLambda = true;
	boolean incrementDryVolFrac = true;
	double totalEinsteinPotential = 0;
	double totalEinsteinEnergyPerVolume;
	double totalSumOfEnergiesPerVolume, totalEnergy;
	double floryRehnerfreeEnergyPerVolume;
	double einsteinFreeEnergyPerVolume;
	double totalPairEnergyPerVolume;
	double totalPairEnergy;
	boolean isFirstIteration;
	double ePairSum=0;
	double calculatedPressure;
	double delta1;
	double idealEinsteinFreeEnergyPerVolume;
	double boltzmannFactorAccumulatorPerSteps;
	double initialFreeEnergy, latestFreeEnergy, initialVolume, latestVolume, volumeChange, changeInFreeEnergy;
	double helmhotzFreeEnergyPerVolume, helmhotzFreeEnergy, centerOfMassContribution;
	double springEnergySum=0, delta2;
	double uLattice, boltzmannFactorPerN, einsteinFreeEnergy, idealContribution1, idealContribution2;

	List<Double> dryVolFracs = new ArrayList<>();
	List<Double> totalSumOfEnergiesList = new ArrayList<>();
	List<Double> totalPairEnergyPerVolumes = new ArrayList<>();
	List<Double> meanFreeEnergyPerVolumes = new ArrayList<>();
	List<Double> initialEnergies = new ArrayList<>();
	List<Double> meanPressures = new ArrayList<>();
	List<Double> calculatedPressures = new ArrayList<>();
	List<Double> totalVolumes = new ArrayList<>();
	List<Double> totalEnergyList = new ArrayList<>();
	List<Double> helmhotzFreeEnergyPerVolumes = new ArrayList<>();
	List<Double> delta1List = new ArrayList<>();
	List<Double> delta2List = new ArrayList<>();

	DecimalFormat decimalFormat = new DecimalFormat("#.###"); // to round my dryVolFrac values to 3 dp
	
    /**
	 * Initializes the model.
	 */
	public void initialize() {
		
		added = false;
		particles.dlambda = control.getDouble("Lambda increment");//the coupling constant increment
		dryVolFracStart = control.getDouble("DryVolFracStart");
		dryVolFracMax = control.getDouble("DryVolFrac Max");
		particles.dphi = control.getDouble("DryVolFrac increment");
		particles.springConstant = control.getDouble("Spring constant"); //the spring constant
		particles.N = control.getInt("N"); // number of particles
		String configuration = control.getString("Initial configuration");
		particles.initConfig = configuration;
        particles.dryR = control.getDouble("Dry radius [nm]");
        particles.xLinkFrac = control.getDouble("x-link fraction");
        particles.Young = control.getDouble("Young's calibration"); // 10-1000
		particles.chi = control.getDouble("chi"); // Flory-Rehner interaction parameter
		particles.tolerance = control.getDouble("Displacement tolerance");
		particles.atolerance = control.getDouble("Radius change tolerance");
		particles.delay = control.getDouble("Delay");
		particles.snapshotInterval= control.getInt("Snapshot interval");
		particles.stop = control.getInt("Stop");
        particles.maxRadius = control.getDouble("Maximum radial distance");
		particles.sizeBinWidth = control.getDouble("Size bin width");
		particles.grBinWidth = control.getDouble("g(r) bin width");
		particles.deltaK = control.getDouble("Delta k");
		particles.fileExtension = control.getString("File extension");
		structure = control.getBoolean("Calculate structure");

		/*
		 * set the value of dryVolFrac to the initial starting value of the 
		 * dry Volume Fraction after which it becomes false every other times
		 */
		if (incrementDryVolFrac){
			dryVolFrac = dryVolFracStart;
			particles.dryVolFrac = dryVolFrac;
			incrementDryVolFrac = false;
		}

		particles.initialize(configuration);

		/*
		 * This block of code sets the original value of lambda to
		 * particles.dlamba * 0.5 only once and that's the first time only!
		 */

		if (setLambda){
			lambda = particles.dlambda * 0.5;
			particles.lambda = lambda;
			setLambda = false;
		}

		/*
		 * sets the first value of initialFreeEnergy and initialVolume
		 * to first value in the list of totalSumOfEnergiesList and totalVolumes
		 * 
		 * This is used for the calculation of the meanPressure() for comparison
		 * with the one from the varial theorem
		 */

		if (isFirstIteration) {
			initialFreeEnergy = totalSumOfEnergiesList.get(0);
			initialVolume = totalVolumes.get(0);
			isFirstIteration = false;
		}

		

        // write out system parameters
        System.out.println("nMon: " + particles.nMon);
        System.out.println("nChains: " + particles.nChains);
        System.out.println("reservoirSwellingRatio: " + particles.reservoirSR);
        System.out.println("reservoirVolFrac: " + particles.reservoirVolFrac);

		if(display3d != null) display3d.dispose(); // closes old simulation frame if present
		display3d = new Display3DFrame("Simulation animation");
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
		
		// initialize visualization elements for particles
		for (int i = 0; i < particles.N; i++) {
			nanoSphere[i].setSizeXYZ(particles.a[i]*2, particles.a[i]*2, particles.a[i]*2);
			nanoSphere[i].getStyle().setFillColor(Color.RED);
			nanoSphere[i].setXYZ(particles.x[i], particles.y[i], particles.z[i]);
		}
		
		if (structure){
		    rdf = new RDF(particles.x,particles.y,particles.z,particles.side,particles.grBinWidth,particles.fileExtension);
		    ssf = new SSF(particles.x,particles.y,particles.z,particles.side,particles.d,particles.deltaK,particles.fileExtension);
		}
	}

	
	/**
	 * Does a simulation step.
	 */
	public void doStep() { // logical step in the HertzSpheres class	    	    
		particles.step();

        // if initial configuration is random, no particle interactions for delay/10 
        if (particles.steps > particles.delay/10.){
            particles.scale = 1.;
        }

		if (particles.steps <= particles.stop){
		    if (particles.steps > particles.delay){
				if ((particles.steps-particles.delay) % particles.snapshotInterval == 0){
			    	particles.sizeDistribution();
			    	if (structure){ // accumulate statistics for structural properties
						rdf.update(); // g(r)
						ssf.update(); // S(k)
			    	}
			    	particles.calculateVolumeFraction();
				}
		    }
		}
		    
		if (particles.steps == particles.stop){
			System.out.println("DryVolFrac = " + dryVolFrac + ", dryVolFracMax = " + dryVolFracMax);
			if (dryVolFrac < dryVolFracMax){
			//	control.println("Number of MC steps = " + particles.steps);
			//	control.println("Initial Energy U_0 = " + particles.initialEnergy);
				control.println("DryVolFrac = " + dryVolFrac);
				control.println("The coupling constant = " + particles.lambda);
			
				if (lambda==1){
					totalVol = particles.totalVol; // set totalVol to the one in the particles class
					initialEnergy=particles.initialEnergy; //set the intial energy to the one in the particles class

					/* difference between interacting einstein crystal and system of interest */
					delta2 = -(springEnergySum*particles.dlambda)/particles.steps/particles.N;
					
					/* calculate the floryRehner free energy at lambda = 1 */
					floryRehnerfreeEnergyPerVolume = (particles.meanFreeEnergy()*particles.N)/totalVol; 

					/* for computing delta1 */
					uLattice=initialEnergy/particles.N;

					/* update the arrays of the pressures and delta2 */
					meanPressures.add(particles.meanPressure()); // from the varial theorem
					delta2List.add(delta2);

					lambda=0; // set lambda to 0
					particles.lambda = lambda;
					this.initialize();
					return;
				}

				if (lambda==0){
					/* calculates the second part of delta1 */
					boltzmannFactorAccumulatorPerSteps = particles.boltzmannFactorAccumulator/particles.steps;
					boltzmannFactorPerN = Math.log(boltzmannFactorAccumulatorPerSteps)/particles.N;
					delta1 = uLattice-boltzmannFactorPerN;

					/* the center of mass contribution to the free energy */
					centerOfMassContribution = Math.log(particles.N)/particles.N;

					einsteinFreeEnergy=1.5*(1-1/(double)particles.N)*Math.log((particles.springConstant*Math.pow(totalVol, 2/(double)3))/Math.PI); //the einstein free energy contribution
					idealContribution1=(1+2/particles.N)*Math.log(particles.N);
					idealContribution2=Math.log(2*Math.PI)/(2*particles.N);
					idealEinsteinFreeEnergyPerVolume=((0.75/Math.PI)*dryVolFrac*(Math.log(dryVolFrac))-1.0); // the ideal einstein free energy per volume

					/* the excess helmhotz free energy for a crystalline solid */
					helmhotzFreeEnergy = einsteinFreeEnergy-idealContribution1-idealContribution2+1.0+centerOfMassContribution+delta1+delta2;
					helmhotzFreeEnergyPerVolume = (helmhotzFreeEnergy*particles.N)/totalVol;

					/* add up all the free energies together */
					totalSumOfEnergiesPerVolume = idealEinsteinFreeEnergyPerVolume+helmhotzFreeEnergyPerVolume+floryRehnerfreeEnergyPerVolume;
					
					/*
					  Update the array with the latest values of dryVolFrac,
					  totalSumOfEnergies, delta1, delta2 and 
					  helmhotzFreeEnergyPerVolume, meanFreeEnergyPerVolume, 
					  and meanPressure(), and totalVolumes
					*/
					dryVolFracs.add(dryVolFrac);
					meanFreeEnergyPerVolumes.add(floryRehnerfreeEnergyPerVolume);
					totalVolumes.add(totalVol);
					delta1List.add(delta1);
					helmhotzFreeEnergyPerVolumes.add(helmhotzFreeEnergyPerVolume);
					totalSumOfEnergiesList.add(totalSumOfEnergiesPerVolume);

					// set the latestFreeEnergy and latestVolume to the next value in the totalSumOfEnergiesList and totalVolumes
					latestFreeEnergy = totalSumOfEnergiesList.get(totalSumOfEnergiesList.size()-1);
					latestVolume = totalVolumes.get(totalVolumes.size()-1);

					/*calculates the difference in the free energies and volumes*/
					changeInFreeEnergy = latestFreeEnergy - initialFreeEnergy;
					volumeChange = latestVolume - initialVolume;

					//performes the pressure calculation: -F/V - V*d(F/V)/dv
					calculatedPressure = -(totalSumOfEnergiesPerVolume)-totalVol*(changeInFreeEnergy/volumeChange);

					//update the list of the calculatedPressures with the latest value of calculatedPressure
					calculatedPressures.add(calculatedPressure);

					//updates the initialVolume to the latestVolume as well as the firstInitialEnergy to the latestInitialEnergy
					initialVolume = latestVolume;
					initialFreeEnergy = latestFreeEnergy;

					ePairSum = 0; //reset the ePairSum to zero
					dryVolFrac += particles.dphi; //increment the dry volume fraction
					particles.dryVolFrac = dryVolFrac;
					lambda = particles.dlambda*0.5; //reset lambda to be dlambda/2
					particles.lambda = lambda;
					this.initialize();
					return;
				}

				if (lambda<1){
					springEnergySum += particles.springEnergyAccumulator; //accumulate springEnergyAccumulator from 0 to 0.95
					lambda += particles.dlambda; // increment lambda
					particles.lambda = lambda;
				}

				if (lambda>1){
					lambda = 1; //set lambda to 1 for a fully interacting system
					particles.lambda = lambda;
					this.initialize();
					return;
				}

				this.initialize();
				return;

			}
			writeData();
		}
		
        // plot mean energy, pressure, swelling ratio
		energyData.append(0, particles.steps, particles.meanPairEnergy());
		pressureData.append(1, particles.steps, particles.meanPressure());
		sizeData.append(2, particles.steps, particles.meanRadius());
		
		display3d.setMessage("Number of steps: " + particles.steps); // update steps

		if(control.getBoolean("Visualization on")){ // visualization updates
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
		control.setValue("Lambda increment", 0.1);
		control.setValue("DryVolFracStart", 0.02);
		control.setValue("DryVolFrac Max", 0.04);
		control.setValue("DryVolFrac increment", 0.002);
		control.setValue("Initial configuration", "FCC");
		control.setValue("Spring constant", 10000);
		//control.setValue("Initial configuration", "random-FCC");
		control.setValue("N", 108); // number of particles
		//control.setValue("N", 500); for FCC lattice, N/4 should be a perfect cube
        control.setValue("Dry radius [nm]", 50);
        control.setValue("x-link fraction", 0.001);
        //control.setValue("Dry volume fraction", 0.01);
        control.setValue("Young's calibration", 10); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
		//control.setValue("chi", 0.2);
        control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.05);
		control.setValue("Delay", 0); // steps after which statistics collection starts
		control.setValue("Snapshot interval", 100); // steps separating successive samples 
		control.setValue("Stop", 1000); // steps after which statistics collection stops
		control.setValue("Size bin width", .001); // bin width of particle radius histogram 
		control.setValue("g(r) bin width", .005); // bin width of g(r) histogram 
		control.setValue("Delta k", .005); // bin width of S(k) histogram 
		control.setValue("File extension", "1");
		control.setValue("Calculate structure", false); // true means calculate g(r) and S(k)
		control.setAdjustableValue("Visualization on", true);
	}

	public void stop() {
		particles.lambda = lambda;
		control.println("The coupling constant = " + particles.lambda);
    	control.println("Number of MC steps = "+particles.steps);
		control.println("Einstein Free Energy = " + decimalFormat.format(particles.einsteinFreeEnergy()));
        control.println("<E_pair>/N = "+decimalFormat.format(particles.meanPairEnergy()));
        control.println("<F>/N = "+decimalFormat.format(particles.meanFreeEnergy()));
    	control.println("PV/NkT = "+decimalFormat.format(particles.meanPressure()));
	}

	public void writeData(){
	    
        for (int i = 0; i < particles.maxRadius/particles.grBinWidth; i++){
			particles.sizeDist[i] = particles.sizeDist[i]/((particles.stop-particles.delay)/particles.snapshotInterval); // normalize size distribution
	    }
	    
	    particles.volFrac = particles.volFrac/((particles.stop-particles.delay)/particles.snapshotInterval); // average system volume fraction in equilibrium 
	    
	    // write system parameters to a file in the data subdirectory
	    try{
			File systemInfo = new File("data/systemInfo" + particles.fileExtension + ".txt");
		
			if (!systemInfo.exists()){ // if file doesn't exist, create it
		    systemInfo.createNewFile();
			}
		
			FileWriter fw = new FileWriter(systemInfo.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
		
			bw.write("Number of particles: " + particles.N);
			bw.newLine();
            bw.write("Initial configuration: " + particles.initConfig);
            bw.newLine();
            bw.write("Dry microgel radius [nm]: " + particles.dryR);
            bw.newLine();
			bw.write("Box length [units of dry radius]: " + particles.side);
			bw.newLine();
			bw.write("Number of monomers: " + particles.nMon);
			bw.newLine();
			bw.write("Number of chains: " + particles.nChains);
			bw.newLine();
			bw.write("Flory interaction parameter (chi): " + particles.chi);
			bw.newLine();
            bw.write("Young's calibration factor: " + particles.Young);
            bw.newLine();
            bw.write("x-link fraction: " + particles.xLinkFrac);
            bw.newLine();
            bw.write("Reservoir swelling ratio: " + particles.reservoirSR);
            bw.newLine();
			bw.write("Mean swelling ratio: " + particles.meanRadius());
			bw.newLine();
            bw.write("Dry volume fraction: " + particles.dryVolFrac);
            bw.newLine();
			bw.write("Reservoir volume fraction: " + particles.reservoirVolFrac);
			bw.newLine();
			bw.write("Actual volume fraction: " + particles.volFrac);
			bw.newLine();
			bw.write("MC steps: " + particles.steps);
			bw.newLine();
			bw.write("Equilibration steps: " + particles.delay);
			bw.newLine();
			bw.write("Snapshot interval: " + particles.snapshotInterval);
			bw.newLine();
			bw.write("Displacement tolerance: " + particles.tolerance);
			bw.newLine();
			bw.write("Particle radius change tolerance: " + particles.atolerance);
			bw.newLine();
			bw.write("Particle radius bin width: " + particles.sizeBinWidth);
			bw.newLine();
			bw.write("g(r) bin width: " + particles.grBinWidth);
			bw.newLine();
            bw.write("Mean pair energy <E_pair>/N [kT]: " + particles.meanPairEnergy());
            bw.newLine();
			bw.write("Mean pressure PV/NkT: " + particles.meanPressure());
			bw.newLine();
		
			bw.newLine();
            bw.write("Mean free energy per particle <F>/N [kT]: " + particles.meanFreeEnergy());
			bw.close();
	    }
	    
	    catch (IOException e) {
			e.printStackTrace();
	    }
	    
	    // write size distribution to a file in the data subdirectory
	    
	    try {
			File sizeFile = new File("data/microgelSize" + particles.fileExtension + ".txt");
		
			if (!sizeFile.exists()) { // if file doesn't exist, create it
		    	sizeFile.createNewFile();
			}
		
			FileWriter fwrite = new FileWriter(sizeFile.getAbsoluteFile());
			BufferedWriter bwrite = new BufferedWriter(fwrite);
			for (int i=0; i<particles.numberBins; i++){
		    	bwrite.write(i + " " + particles.sizeDist[i]);
		    	bwrite.newLine();
			}
		
			bwrite.close();
		
	    }
	    
	    catch (IOException e) {
			e.printStackTrace();
	    } 

		try {
			File outputFile = new File("data/dryVolFracVsFreeEnergySolidPhase.txt");
		
			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}

			FileWriter fw1 = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);

			// write system parameters to the file in addittion to dryVolFracs and totalFree energies
			bw1.write("Number of particles: " + particles.N);
			bw1.newLine();
			bw1.write("Initial configuration: " + particles.initConfig);
			bw1.newLine();
			bw1.write("Dry microgel radius [nm]: " + particles.dryR);
			bw1.newLine();
			bw1.write("Box length [units of dry radius]: " + particles.side);
			bw1.newLine();
			bw1.write("Number of monomers: " + particles.nMon);
			bw1.newLine();
			bw1.write("Number of chains: " + particles.nChains);
			bw1.newLine();
			bw1.write("Flory interaction parameter (chi): " + particles.chi);
			bw1.newLine();
			bw1.write("Young's calibration factor: " + particles.Young);
			bw1.newLine();
			bw1.write("x-link fraction: " + particles.xLinkFrac);
			bw1.newLine();
			bw1.write("Reservoir swelling ratio: " + particles.reservoirSR);
			bw1.newLine();
			bw1.write("Mean swelling ratio: " + particles.meanRadius());
			bw1.newLine();
			bw1.write("Dry volume fraction: " + particles.dryVolFrac);
			bw1.newLine();
			bw1.write("Reservoir volume fraction: " + particles.reservoirVolFrac);
			bw1.newLine();
			bw1.write("Actual volume fraction: " + particles.volFrac);
			bw1.newLine();
			bw1.write("MC steps: " + particles.steps);
			bw1.newLine();
			bw1.write("Equilibration steps: " + particles.delay);
			bw1.newLine();
			bw1.write("Snapshot interval: " + particles.snapshotInterval);
			bw1.newLine();
			bw1.write("Displacement tolerance: " + particles.tolerance);
			bw1.newLine();
			bw1.write("Particle radius change tolerance: " + particles.atolerance);
			bw1.newLine();
			bw1.write("Particle radius bin width: " + particles.sizeBinWidth);
			bw1.newLine();
			bw1.write("g(r) bin width: " + particles.grBinWidth);
			bw1.newLine();
			bw1.write("Mean pair energy <E_pair>/N [kT]: " + particles.meanPairEnergy());
			bw1.newLine();
			bw1.write("The coupling constant increment - dlambda: " + particles.dlambda);
			bw1.newLine();	
			bw1.write("DryVolFrac, totalSumOfEnergies, calculatedPressures(), meanPressures(), delta1, delta2");
			bw1.newLine();	

			for (int i = 0; i < dryVolFracs.size(); i++) {
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				bw1.write(roundedDryVolFrac + ", " + totalSumOfEnergiesList.get(i) + ", " + calculatedPressures.get(i) + ", " + meanPressures.get(i) + ", " + delta1List.get(i) + ", " + delta2List.get(i));
				bw1.newLine();
			}

			bw1.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}

            // write radial distribution function, static structure factor to files in data subdirectory

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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresSolidDelta2TrialApp());
	}
}
