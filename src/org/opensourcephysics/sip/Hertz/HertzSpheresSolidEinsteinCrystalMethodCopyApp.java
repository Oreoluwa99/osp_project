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

// the edited copy of the HertzSpheresCode for the solid phase

/**
 * HertzSpheresApp performs a Monte Carlo simulation of microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @authors Alan Denton and Matt Urich
 * @version 1.2 25-05-2021
 * 
 */

 public class HertzSpheresSolidEinsteinCrystalMethodCopyApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
	HertzSpheresSolidEinsteinCrystalMethod particles = new HertzSpheresSolidEinsteinCrystalMethod(); // with equilibration
//	HertzSpheresSolidEinsteinCrystalCopy particles = new HertzSpheresSolidEinsteinCrystalCopy(); // without equilibration
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
	double lambda, totalVol;
	boolean setLambda = true;
	boolean incrementDryVolFrac = true;
	double deltaF2Accumulator = 0, deltaF2, deltaF1, deltaF1PerVol, deltaF2PerVol;
	double upperLimit, lowerLimit, xIncrement;
	double f0, f0overN, f0PerVolume;
	double x;
	double uLatticePerN;
	double fPair, fPairPerVol, totalFreeEnergy;
	double floryFEperVol;
	boolean isFirstIteration;
	double initialFreeEnergy, latestFreeEnergy, initialVolume, latestVolume;
	double calculatedPressure;
	double volumeChange, freeEnergyChange;
	double initialFloryFreeEnergy, latestFloryRehnerFreeEnergy;
	double latestFloryFreeEnergy, latestPairFreeEnergy;
	double floryPressure, dFRFE_dPhi, floryPContribution, initialPairFreeEnergy, dFPair_dPhi, pairPressure, pairPContribution;

	/* Lists to store various calculated values */
	List<Double> dryVolFracs = new ArrayList<>();
	List<Double> delta2List = new ArrayList<>();
	List<Double> deltaF1List = new ArrayList<>();
	List<Double> deltaF2PerVolumeList = new ArrayList<>();
	List<Double> f0PerVolumeList = new ArrayList<>();
	List<Double> boltzmannFactorPerNList = new ArrayList<>();
	List<Double> floryFEperVolList = new ArrayList<>();
	List<Double> totalSumOfEnergiesList = new ArrayList<>();
	List<Double> totalVolList = new ArrayList<>();
	List<Double> calculatedPressures = new ArrayList<>();
	List<Double> meanPressures = new ArrayList<>();
	List<Double> pairFreeEnergyList = new ArrayList<>();
	List<Double> floryRehnerPressuresList = new ArrayList<>();
	List<Double> pairPressuresList = new ArrayList<>();
	List<Double> varialPressureList = new ArrayList<>();
	List<Double> reservoirVolFracList = new ArrayList<>();
	List<Double> swellingRatioList = new ArrayList<>();

	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values to 3 dp
	
    /**
	 * Initializes the model.
	 */
	public void initialize() {
		
		added = false;
		particles.dlambda = control.getDouble("Lambda increment");//the coupling constant increment
		dryVolFracStart = control.getDouble("DryVolFracStart");
		dryVolFracMax = control.getDouble("DryVolFrac Max");
		particles.dphi = control.getDouble("DryVolFrac increment");
		xIncrement = control.getDouble("limit increment");//to increment the limits
		particles.springConstant = control.getDouble("Spring constant"); //the spring constant
	//	upperLimit = control.getDouble("upper limit");
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


		/* change of variables for calculating delta2 */
		upperLimit = Math.log(particles.springConstant+Math.exp(3.5));
		lowerLimit = 3.5;
		xIncrement = (upperLimit-lowerLimit)/10.0;

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
		 * set the value of x and lambda from changing the variables in the integration
		 */
	
		if (setLambda){
			x = lowerLimit + 0.5*xIncrement; // set the value of x (for the change of variables)
			lambda = (Math.exp(x)-Math.exp(3.5))/particles.springConstant; 
			particles.lambda = lambda;
			setLambda = false;
		}

		/*
		 sets the value of initialpairEnergyTotalPerVolume to the first 
		 value in the list of pairEnergyTotalPerVolumes
		*/

		/*
		if (isFirstIteration) {
			initialFreeEnergy = totalSumOfEnergiesList.get(0); // the first initial energy in the list
			initialFloryFreeEnergy = floryFEperVolList.get(0); // first flory-rehner free energy in the list
			initialPairFreeEnergy = pairFreeEnergyList.get(0);
			initialVolume = totalVolList.get(0); //the first volume in the list
			isFirstIteration = false;
		}
		*/

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
				control.println("DryVolFrac = " + dryVolFrac);
			//	control.println("Coupling constant = "+lambda);
				control.println("The value of x = "+x);
				if (x<upperLimit){
					// 	accumulate deltaF2 
					deltaF2Accumulator += (particles.meanSpringEnergyAccumulator()*Math.exp(x))/particles.springConstant;
					x += xIncrement; // increment x
					lambda = (Math.exp(x)-Math.exp(3.5))/particles.springConstant;
					particles.lambda = lambda;

					if (x>upperLimit){
						/* calculate delta F2 */
						deltaF2 = -(deltaF2Accumulator*xIncrement);
						deltaF2PerVol = (deltaF2*particles.N)/particles.totalVol;

						/* Calculate a0: the addittion of delta3(change between unconstrained solid and solid with fixed center of mass) and the free energy of the refrence ideal Einstein crystal */
						f0 = -3*(particles.N-1)/2.0*Math.log(Math.PI/particles.springConstant)-Math.log(particles.N)/2.0-Math.log(particles.totalVol);
						f0overN = f0/particles.N;
						f0PerVolume = (f0overN*particles.N)/particles.totalVol;

						lambda = 1; // set the value of lambda to 1
						particles.lambda = lambda;

						/* update the f0 and deltaF2 lists */
						f0PerVolumeList.add(f0PerVolume);
						deltaF2PerVolumeList.add(deltaF2PerVol);

						this.initialize();
						return;
						
					}
					this.initialize();
					return;
				}

					if (lambda == 1){
						/* Calculate the Flory–Rehner free energy */
						floryFEperVol = (particles.meanFreeEnergy()*particles.N)/particles.totalVol;
						lambda = 0; // set lambda to 0
						particles.lambda = lambda;

						/* update the meanFreeEnergy list */
						floryFEperVolList.add(floryFEperVol);

						this.initialize();
						return;
					
					}

					if (lambda == 0){
						uLatticePerN=particles.initialEnergy/particles.N; // the constant lattice energy
						/* the free energy change between an interacting einstein crystal and a non-interacting einstein crystal */
						deltaF1 = uLatticePerN-(Math.log(particles.meanboltzmannFactor())/particles.N); 
						deltaF1PerVol = (deltaF1*particles.N)/particles.totalVol;
						deltaF1List.add(deltaF1PerVol);

						fPair = f0overN+deltaF1+deltaF2; // the addition of all the free energy of a molecular solid
					//	control.println("fPair = " + fPair);
						fPairPerVol = (fPair*particles.N)/particles.totalVol; // pair free energy per volume
						pairFreeEnergyList.add(fPairPerVol);
						control.println("pairEPerVol = "+fPairPerVol);

						
						totalFreeEnergy = fPairPerVol+floryFEperVol;
						control.println("totalF " + totalFreeEnergy);

						dryVolFracs.add(dryVolFrac);

						double currentDryVolFrac = dryVolFrac;
						double previousDryVolFrac = currentDryVolFrac - particles.dphi;
						double nextDryVolFrac = currentDryVolFrac + particles.dphi;

						int currentIndex = dryVolFracs.size() - 1; // Index of the current dryVolFrac
						if (currentIndex >= 0) {
							// Calculate the indices of the previous and next dryVolFrac
							int previousIndex = currentIndex - 1;
							int nextIndex = currentIndex + 1;
					
							control.println("currentIndex = "+currentIndex);
							control.println("previousIndex = "+previousIndex);
							control.println("nextIndex = "+nextIndex);

							if (previousIndex >= 0 && nextIndex < dryVolFracs.size()) {
								// Get the values of pairEperVol at the previous and next indices
								double previousPairEperVol = pairFreeEnergyList.get(previousIndex);
								double nextPairEperVol = pairFreeEnergyList.get(nextIndex);
					
								// Calculate the latestPairEperVol
								latestPairFreeEnergy = (nextPairEperVol - previousPairEperVol) / (2.0 * particles.dphi);
								control.println("latestPairEPerVol = " + latestPairFreeEnergy);
					
								// Get the values of floryFperVol at the previous and next indices
								double previousFREperVol = floryFEperVolList.get(previousIndex);
								double nextFREperVol = floryFEperVolList.get(nextIndex);
					
								// Calculate the latestFREperVol
								latestFloryFreeEnergy = (nextFREperVol - previousFREperVol) / (2.0 * particles.dphi);
								control.println("latestFREperVol = " + latestFloryFreeEnergy);
							}
						
								dFRFE_dPhi = (latestFloryFreeEnergy-initialFloryFreeEnergy)/(2*particles.dphi); // change in free energy by dryVolFrac change
								dFPair_dPhi = (latestPairFreeEnergy - initialPairFreeEnergy)/(2*particles.dphi);

								floryPressure = ((dFRFE_dPhi)*dryVolFrac-floryFEperVol);
								pairPressure = ((dFPair_dPhi)*dryVolFrac-fPairPerVol);
						
								floryPContribution = (floryPressure*particles.totalVol)/particles.N; // the pressure contribution from the FloryRehner energy
								floryRehnerPressuresList.add(floryPContribution); // (PV/NKT)flory-rehner

								pairPContribution = (pairPressure*particles.totalVol)/particles.N;
								pairPressuresList.add(pairPContribution+1.0); // (PV/NKT)pair

								initialFreeEnergy = totalFreeEnergy;
								initialFloryFreeEnergy = latestFloryFreeEnergy; // update the initial flory-rehner energy
								initialPairFreeEnergy = latestPairFreeEnergy; // update the initial pair free energy
							//	initialVolume = latestVolume; // reset initial volume

					}

						/* add up all the free energies contributing to the interaction of the particles */
						
						x = lowerLimit + 0.5*xIncrement; // reset the value of x
						dryVolFrac += particles.dphi; // increment the dry volume fraction
						particles.dryVolFrac = dryVolFrac;
						lambda = (Math.exp(x)-Math.exp(3.5))/particles.springConstant; // reset lambda
						particles.lambda = lambda;
						deltaF2Accumulator = 0; // reset deltaF2Accumulator to 0
						
						this.initialize();
						return;
					}
				

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
	//	control.setValue("upper limit", upperLimit);
		control.setValue("DryVolFrac Max", 0.04);
		control.setValue("DryVolFrac increment", 0.002);
		control.setValue("Initial configuration", "FCC");
		control.setValue("limit increment", xIncrement);
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
			File outputFile = new File("data/dryVolFracVsFreeEnergySolidEinsteinCrystalMethodFile.txt");
		
			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}

			FileWriter fw1 = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);

			// write system parameters to the file 
			bw1.write("Number of particles: " + particles.N);
			bw1.newLine();
			bw1.write("springConstant: " + particles.springConstant);
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
			bw1.write("DryVolFrac,      deltaF1,        deltaF2,      f0,    floryRehnerFreeEnergy,   totalFreeEnergy,   (PV/NKT)pairVarial,	(PV/NKT)pairCalculated, (PV/NKT)floryRehner, PV/NKT(totalCalculated)");
			bw1.newLine();	

			for (int i = 0; i < dryVolFracs.size(); i++) {
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				bw1.write(roundedDryVolFrac + ", " + deltaF1List.get(i) + ", " + deltaF2PerVolumeList.get(i) + ", " + f0PerVolumeList.get(i) + ", " + floryFEperVolList.get(i) + ", " + totalSumOfEnergiesList.get(i) + ", " + meanPressures.get(i) + ", " + pairPressuresList.get(i) + ", " + floryRehnerPressuresList.get(i) + ", " + calculatedPressures.get(i));
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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresSolidEinsteinCrystalMethodCopyApp());
	}
}
