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
 * 
 * This Java program performs a Monte Carlo simulation of soft microgel particles 
 * arranged in a crystalline (FCC) configuration and interacting via the Hertzian 
 * elastic potential. The swelling behavior of each particle is modeled using the 
 * Flory–Rehner theory, and the free energy of the system is computed using the 
 * Frenkel–Ladd method adapted for deformable particles.
 * 
 * This works just fine for both the interpenetration and facet algorithm.
 * 
  * Key features:
 * - Iterative adjustment of dry volume fraction to explore phase behavior
 * - Gauss–Legendre quadrature integration for free energy calculations
 * - Calculation of:
 *    • Pairwise interaction free energy
 *    • Flory–Rehner free energy
 *    • Total Helmholtz free energy per volume
 *    • Lindemann parameter for melting criterion
 *    • Nearest-neighbor distances and structure-based potentials
 * - Outputs radial distribution function (g(r)) and static structure factor (S(k))
 * - Data export to text files for further analysis
 * - Real-time 3D visualization of particle configurations using Open Source Physics
 * 
 * 
 * @authors Alan Denton and Oreoluwa Alade
 * @version 1.2 April 2025
 * 
 */

public class HertzSpheresSolidPhaseFreeEnergyApp extends AbstractSimulation {
	public enum WriteModes {WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
	// HertzSpheresSolidPhase particles = new HertzSpheresSolidPhase(); 	
	//HertzSpheresSolidPhaseCopy particles = new HertzSpheresSolidPhaseCopy();															
	//HertzSpheresNonLocalFacet_AlansEdits particles = new HertzSpheresNonLocalFacet_AlansEdits();

	/* For the interpenetration algorithm */
	HertzInterpenetrationFreeEnergies particles = new HertzInterpenetrationFreeEnergies(); 

	/* For the facet algorithm */
	//HertzSpheresFCCNearestNeighbors particles = new HertzSpheresFCCNearestNeighbors(); // For the facet algorithm
	
	PlotFrame energyData = new PlotFrame("MC steps", "<E_pair>/N", "Mean pair energy per particle");
	PlotFrame pressureData = new PlotFrame("MC steps", "PV/NkT", "Mean pressure");
	PlotFrame sizeData = new PlotFrame("MC steps", "alpha", "Mean swelling ratio");
	Display3DFrame display3d = new Display3DFrame("Simulation animation");
	int dataPoints;
	int maxDataPoints;
	int weightIteration = 0, pointIteration = 0; // set and reset the iterations
	ElementSphere nanoSphere[];
	boolean added = false;
	boolean structure;
	RDF rdf;
	SSF ssf;
	double dryVolFracStart, dryVolFracMax, dryVolFrac;
	double lambda, totalVol;
	boolean setLambda = true;
	boolean incrementDryVolFrac = true;
	double variableChanged;
	double einsteinFreeEnergy;
	double deltaFAccumulator = 0;
	double gaussPoint, gaussWeight;
	double fPair, fPairPerVol, totalFreeEnergy;
	double floryFEperVol;
	double initialFreeEnergy, latestFreeEnergy, initialVolume, latestVolume;
	double UPairMinusUSpringPerVolume;
	double floryPressure, dFRFE_dPhi, floryPContribution, initialPairFreeEnergy, dFPair_dPhi, pairPressure, pairPContribution;
	double chemicalPotential;
	double uPairPerVol;
	double density, nnDistance;
	double newHertzianPotential;
	boolean gaussLegendrePoint = true;
	boolean gaussLegendreWeight = true;
	double deltaFPerVol;
	double springConstant;
	double referenceFR, referenceFRPerN, referenceFRPerVol;
	double lindemannParameter;

	/* Lists to store various calculated values */
	List<Double> dryVolFracs = new ArrayList<>();
	List<Double> dryVolFracsEdited = new ArrayList<>();
	List<Double> floryFEperVolList = new ArrayList<>();
	List<Double> totalSumOfEnergiesList = new ArrayList<>();
	List<Double> totalVolList = new ArrayList<>();
	List<Double> calculatedPressures = new ArrayList<>();
	List<Double> meanPressures = new ArrayList<>();
	List<Double> pairFreeEnergyList = new ArrayList<>();
	List<Double> floryRehnerPressuresList = new ArrayList<>();
	List<Double> pairPressuresList = new ArrayList<>();
	List<Double> floryRehnerPressuresListEdited = new ArrayList<>();
	List<Double> pairPressuresListEdited = new ArrayList<>();
	List<Double> chemicalPotentialList = new ArrayList<>();
	List<Double> reservoirVolFracList = new ArrayList<>();
	List<Double> swellingRatioList = new ArrayList<>();
	List<Double> uPairPerVolList = new ArrayList<>();
	List<Double> newHertzianPotentialList = new ArrayList<>();
	List<Double> gaussLegendreWeights = new ArrayList<>();
	List<Double> gaussLegendrePoints = new ArrayList<>();
	List<Double> referenceFRPerVolList = new ArrayList<>();
	List<Double> springConstantList = new ArrayList<>();
	List<Double> lindemannParameterList = new ArrayList<>();

	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values to 3 dp

	/**
	 * Initializes the model.
	 */
	public void initialize() {

		added = false;
		//particles.dlambda = control.getDouble("Lambda increment");// the coupling constant increment
		dryVolFracStart = control.getDouble("DryVolFracStart");
		dryVolFracMax = control.getDouble("DryVolFrac Max");
		particles.dphi = control.getDouble("DryVolFrac increment");
		//xIncrement = control.getDouble("limit increment");// to increment the limits
		// particles.springConstant = control.getDouble("Spring constant"); // the spring constant
		// upperLimit = control.getDouble("upper limit");
		//particles.xLinkFrac = control.getDouble("x-link fraction");
		//dxLink = control.getDouble("x-Link increment");
		//xLinkFracMax = control.getDouble("x-link fraction max");
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
		particles.snapshotInterval = control.getInt("Snapshot interval");
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
		if (incrementDryVolFrac) {
			dryVolFrac = dryVolFracStart;
			particles.dryVolFrac = dryVolFrac;
			incrementDryVolFrac = false;
		}

		particles.initialize(configuration);

		if (gaussLegendreWeight){	
			// initialize the gaussLegendreWeights list	
			gaussLegendreWeights.add(0.23692689);
			gaussLegendreWeights.add(0.47862867);
			gaussLegendreWeights.add(0.56888889);
			gaussLegendreWeights.add(0.47862867);
			gaussLegendreWeights.add(0.23692689);
			gaussLegendreWeight = false;
		}

		if (gaussLegendrePoint){
			//initialize the gaussLegendrePoints list
			gaussLegendrePoints.add(-0.90617985);
        	gaussLegendrePoints.add(-0.53846931);
        	gaussLegendrePoints.add(0.);
        	gaussLegendrePoints.add(0.53846931);
        	gaussLegendrePoints.add(0.90617985);
			gaussLegendrePoint = false;
		}

		/*
		 * set the value of lambda
		 */
		
		if (setLambda) {
			lambda = 1;
			particles.lambda = lambda;
			setLambda = false;
		}
		
		// write out system parameters
		System.out.println("nMon: " + particles.nMon);
		System.out.println("nChains: " + particles.nChains);
		System.out.println("reservoirSwellingRatio: " + particles.reservoirSR);
		System.out.println("reservoirVolFrac: " + particles.reservoirVolFrac);

		if (display3d != null)
		display3d.dispose(); // closes old simulation frame if present
		display3d = new Display3DFrame("Simulation animation");
		display3d.setPreferredMinMax(0, particles.side, 0, particles.side, 0, particles.side);
		display3d.setSquareAspect(true);
		// energyData.setPreferredMinMax(0, 1000, -10, 10);
		// pressureData.setPreferredMinMax(0, 1000, 0, 2);

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
			nanoSphere[i].setSizeXYZ(particles.a[i] * 2, particles.a[i] * 2, particles.a[i] * 2);
			nanoSphere[i].getStyle().setFillColor(Color.RED);
			nanoSphere[i].setXYZ(particles.x[i], particles.y[i], particles.z[i]);
		}

		if (structure) {
			rdf = new RDF(particles.x, particles.y, particles.z, particles.side, particles.grBinWidth, particles.fileExtension);
			ssf = new SSF(particles.x, particles.y, particles.z, particles.side, particles.d, particles.deltaK, particles.fileExtension);
		}
	}

	/**
	 * Does a simulation step.
	 */
	public void doStep() { // logical step in the HertzSpheres class
		particles.step();

		// if initial configuration is random, no particle interactions for delay/10
		if (particles.steps > particles.delay / 10.) {
			particles.scale = 1.;
		}

		if (particles.steps <= particles.stop) {
			if (particles.steps > particles.delay) {
				if ((particles.steps - particles.delay) % particles.snapshotInterval == 0) {
					particles.sizeDistribution();
					if (structure) { // accumulate statistics for structural properties
						rdf.update(); // g(r)
						ssf.update(); // S(k)
					}
					particles.calculateVolumeFraction();
				}
			}
		}

		if (particles.steps == particles.stop) {
			System.out.println("DryVolFrac = " + dryVolFrac + ", dryVolFracMax = " + dryVolFracMax);
			if (dryVolFrac < dryVolFracMax) {
				control.println("phi0: "+dryVolFrac);
				control.println("xLinkFrac: "+ particles.xLinkFrac);
				if (lambda == 1) { // the real or interacting solid

					/* Compute the spring constant at lambda = 1*/
					springConstant = 3/(2*particles.meanSquareDisplacement()); // the springConstant
					
					control.println("springConstant: "+springConstant);

					particles.springConstant = springConstant;
					springConstantList.add(springConstant); // update the list
					
					/* Calculate the Flory–Rehner free energy for a fully interacting system */
					floryFEperVol = (particles.meanFreeEnergy()*particles.N)/particles.totalVol; // Flory-RehnerFR/vol
					floryFEperVolList.add(floryFEperVol); // update the list

					control.println("floryFEperVol: "+floryFEperVol);

					// the reference free energy (the free energy of the ideal Einstein Crystal)
					referenceFR = -3*(particles.N-1)/2.0*Math.log(Math.PI/particles.springConstant)-Math.log(particles.N)/2.0-Math.log(particles.totalVol);
					referenceFRPerN = referenceFR/particles.N; // reference free energy per particles
					referenceFRPerVol = (referenceFRPerN*particles.N)/particles.totalVol; // reference free energy per volume
					referenceFRPerVolList.add(referenceFRPerVol); // update the list

					//einsteinFreeEnergy = particles.einsteinFreeEnergy()*(particles.N/particles.totalVol); //quantity F_Ein
					//control.println("F_Ein: "+einsteinFreeEnergy);

					meanPressures.add(particles.meanPressure()); // pressure from virial theorem
					density = particles.N/particles.totalVol; 
					nnDistance = (1/Math.sqrt(2))*Math.pow((4/density), 1/3.0); // the nearest neighbour distance
					
					// the lindemannParameter 
					lindemannParameter = Math.sqrt(particles.meanSquareDisplacement())/nnDistance;
					lindemannParameterList.add(lindemannParameter);
					control.println("lindemannParameter: "+lindemannParameter);

					// as recommended by Zacahareli and co
					if (nnDistance>(2*particles.meanRadius())){
						newHertzianPotential = 0;
					}
					else{
						newHertzianPotential = particles.B*Math.pow((1-nnDistance/(2*particles.meanRadius())), 2.5);
					}

					control.println("newHertzianPotential: "+newHertzianPotential);

					newHertzianPotentialList.add(newHertzianPotential);
					
					uPairPerVol = particles.meanPairEnergy()*(particles.N/particles.totalVol); // the pair energy per volume
					
					control.println("uPairPerVol: "+uPairPerVol);

					gaussPoint = gaussLegendrePoints.get(pointIteration); // the Gauss-Legendre point
					gaussWeight = gaussLegendreWeights.get(weightIteration); // the Gauss-Legendre weight
					// term coming from changing the integration limits
					variableChanged = 0.5*(gaussPoint*(Math.log(springConstant+Math.exp(3.5))-3.5)+3.5+Math.log(springConstant+Math.exp(3.5)));
					lambda = (Math.exp(variableChanged)-Math.exp(3.5))/particles.springConstant;
					particles.lambda = lambda;
				
					this.initialize();
					return;
				}
				
				else{
					
					deltaFAccumulator += (0.5*gaussWeight)*(Math.log(particles.springConstant+Math.exp(3.5))-3.5)*(Math.exp(variableChanged))/particles.springConstant*(particles.meanPairEnergy()-particles.meanSpringEnergy()); //<U_pair - U_spring>
					// increment the counters
					pointIteration++;
					weightIteration++;	

					if (pointIteration < gaussLegendrePoints.size() && weightIteration < gaussLegendreWeights.size()) {				
						gaussWeight = gaussLegendreWeights.get(weightIteration);
						gaussPoint = gaussLegendrePoints.get(pointIteration); // set the value of x (for the change of variables)
						// comes from changing the limits for the Gauss-Legendre integration
						variableChanged = 0.5*(gaussPoint*(Math.log(springConstant+Math.exp(3.5))-3.5)+3.5+Math.log(springConstant+Math.exp(3.5)));
						// compute the value of lambda
						lambda = (Math.exp(variableChanged)-Math.exp(3.5))/particles.springConstant;
						particles.lambda = lambda;

						this.initialize();
						return;
					}	

					else{
						deltaFPerVol = deltaFAccumulator*(particles.N/particles.totalVol); //ΔF/V
						//control.println("UpairMinusUspringSum: "+UPairMinusUSpringSum);
						
						fPairPerVol = deltaFPerVol; // the pairEnergy of Interaction

						control.println("fPairPerVol: "+fPairPerVol);

						// the total free energy
						totalFreeEnergy = referenceFRPerVol+deltaFPerVol+floryFEperVol;
						control.println("totalF: " + totalFreeEnergy);

						// update the array lists 
						dryVolFracs.add(dryVolFrac);
						uPairPerVolList.add(uPairPerVol);
						totalVolList.add(particles.totalVol);
						reservoirVolFracList.add(particles.reservoirVolFrac); // the reservoir volume fraction list
						swellingRatioList.add(particles.meanRadius()); // the swelling ratio list
						pairFreeEnergyList.add(fPairPerVol);
						totalSumOfEnergiesList.add(totalFreeEnergy);

						// increment the dry volume fraction
						dryVolFrac += particles.dphi;
						particles.dryVolFrac = dryVolFrac;
						//reset lambda to 1 for the interacting solid
						lambda = 1; 
						particles.lambda = lambda;
						deltaFAccumulator = 0; //reset UPairMinusUSpring
						//reset the iterations
						pointIteration = 0;
						weightIteration = 0;

						this.initialize();
						return;
					}
			    }
				
			}

			if (dryVolFracs.size() > 2){
				// Variables to store free energy values for the current, previous, and next pairs
				double currentPairFreeEnergy = Double.NaN;
				double previousPairFreeEnergy = Double.NaN;
				double nextPairFreeEnergy = Double.NaN;

				// Variables to store free energy values for the current, previous, and next Flory interactions
				double currentFloryFreeEnergy = Double.NaN;
				double previousFloryFreeEnergy = Double.NaN;
				double nextFloryFreeEnergy = Double.NaN;
				int i;

				/* Initialize the lists */
				floryRehnerPressuresListEdited.add(0, 0.0);
				pairPressuresListEdited.add(0, 0.0);
				calculatedPressures.add(0,0.0);
				chemicalPotentialList.add(0,0.0);
				
				// Iterate through dry volume fractions, excluding the first and last elements
				for (i = 1; i < dryVolFracs.size() - 1; i++) {
					// Retrieve the current, previous, and next dry volume fractions
					double currentDryVolFrac = dryVolFracs.get(i);
					double previousDryVolFrac = dryVolFracs.get(i - 1);
					double nextDryVolFrac = (i + 1 < dryVolFracs.size()) ? dryVolFracs.get(i + 1) : Double.NaN;

					// Store the current dry volume fraction in a modified list
					dryVolFracsEdited.add(currentDryVolFrac);
					// Retrieve the last value from the modified list
					double value = dryVolFracsEdited.get(i - 1);
					control.println("value = " + value);

					double currentVolume = totalVolList.get(i); // the volume of the system for that dry volume fraction
					control.println("currentVolumeIndex = " + currentVolume);

					control.println("Current Index = " + currentDryVolFrac);
					control.println("Previous Index = " + previousDryVolFrac);
					control.println("Next Index = " + (Double.isNaN(nextDryVolFrac) ? "N/A" : nextDryVolFrac));

					currentPairFreeEnergy = pairFreeEnergyList.get(i); // pair free energy for that dry volume fraction
					currentFloryFreeEnergy = floryFEperVolList.get(i); // FR free energy for that dry volume fraction
					control.println("Current Pair Free Energy = " + currentPairFreeEnergy);
					control.println("Current Flory Free Energy = " + currentFloryFreeEnergy);

					// Retrieve the initial values for free energy from the previous iteration
					initialFreeEnergy = totalSumOfEnergiesList.get(i - 1); 
					previousPairFreeEnergy = pairFreeEnergyList.get(i - 1); 
					previousFloryFreeEnergy = floryFEperVolList.get(i - 1);
					control.println("Previous Pair Free Energy = " + previousPairFreeEnergy);
					control.println("Previous Flory Free Energy = " + previousFloryFreeEnergy);

					// Check if there is a next iteration available
					if (i + 1 < dryVolFracs.size()) {
						// Retrieve the next pair and Flory free energy values
						nextPairFreeEnergy = pairFreeEnergyList.get(i + 1);
						nextFloryFreeEnergy = floryFEperVolList.get(i + 1);
						// Retrieve the latest free energy value from the next iteration
						latestFreeEnergy = totalSumOfEnergiesList.get(i + 1);
					}
					control.println("Next Pair Free Energy = " + nextPairFreeEnergy);
					control.println("Next Flory Free Energy = " + nextFloryFreeEnergy);

					// Calculate chemical potential using the derivative of the free energy with respect to dry volume fraction
					chemicalPotential = (4.0*Math.PI/3.0)*(latestFreeEnergy-initialFreeEnergy)/(2.0*particles.dphi); // dF/dPhi0
					
					// Calculate derivatives of free energies with respect to dry volume fraction
					dFRFE_dPhi = (nextFloryFreeEnergy - previousFloryFreeEnergy) / (2.0 * particles.dphi); // dF_FR/dPhi0
					dFPair_dPhi = (nextPairFreeEnergy - previousPairFreeEnergy) / (2.0 * particles.dphi); // dF_pair/dPhi0

					control.println("dFRFE_dPhi = " + dFRFE_dPhi);
					control.println("dFPair_dPhi = " + dFPair_dPhi);

					// Calculate Flory and pair pressures using the derivatives and free energy values
					floryPressure = ((dFRFE_dPhi) * currentDryVolFrac - currentFloryFreeEnergy); // phi0*dF_pair/dPhi0 - F_pair
					pairPressure = ((dFPair_dPhi) * currentDryVolFrac - currentPairFreeEnergy); //  phi0*dF_FR/dPhi0 - F_FR

					// Calculate contributions to pressure from Flory and pair interactions
					floryPContribution = (floryPressure * currentVolume) / particles.N; // P_FR*V/N(KT)													
					pairPContribution = (pairPressure * currentVolume) / particles.N; // P_pair*V/N(KT)

					control.println("floryContribution = " + floryPContribution);
					//control.println("pairPContribution = " + pairPContribution);

					/* Update the lists accordingly */
					floryRehnerPressuresListEdited.add(floryPContribution); 
					pairPressuresListEdited.add(pairPContribution); // the F0 (eq. 48 from Vegas et. al.) already incldues the ideal gas free energy
					calculatedPressures.add(floryPContribution + pairPContribution); // the calculated pressure from the derivatives of the free energies
					chemicalPotentialList.add(chemicalPotential);
				}
			}
			writeData();
		}
		

		// plot mean energy, pressure, swelling ratio
		energyData.append(0, particles.steps, particles.meanPairEnergy());
		pressureData.append(1, particles.steps, particles.meanPressure());
		sizeData.append(2, particles.steps, particles.meanRadius());

		display3d.setMessage("Number of steps: " + particles.steps); // update steps

		if (control.getBoolean("Visualization on")) { // visualization updates
			for (int i = 0; i < particles.N; i++) {
				nanoSphere[i].setSizeXYZ(particles.a[i] * 2, particles.a[i] * 2, particles.a[i] * 2);
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
		//control.setValue("Lambda increment", 0.1);

		// the original ones
		// control.setValue("DryVolFracStart", 0.001);
		// // control.setValue("upper limit", upperLimit);
		// control.setValue("DryVolFrac Max", 0.003);
		// control.setValue("DryVolFrac increment", 0.0001);

		control.setValue("DryVolFracStart", 0.012);
		control.setValue("DryVolFrac Max", 0.016);
		control.setValue("DryVolFrac increment", 0.0001);
		control.setValue("Initial configuration", "FCC");
		// control.setValue("Spring constant", 10000); // Spring constant: 2.035 for alpha/KT = 100
		// control.setValue("Initial configuration", "random-FCC");
		control.setValue("N", 32); // number of particles
		control.setValue("x-link fraction", 0.001);
		// control.setValue("N", 500); for FCC lattice, N/4 should be a perfect cube
		control.setValue("Dry radius [nm]", 50);
		// control.setValue("Dry volume fraction", 0.01);
		control.setValue("Young's calibration", 1.0); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
		// control.setValue("chi", 0.2);
		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.0); // initially set to 0.05
		control.setValue("Delay", 10000); // steps after which statistics collection starts
		control.setValue("Snapshot interval", 100); // steps separating successive samples
		control.setValue("Stop", 20000); // steps after which statistics collection stops
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
		control.println("Number of MC steps = " + particles.steps);
		control.println("<E_pair>/N = " + decimalFormat.format(particles.meanPairEnergy()));
		control.println("<F>/N = " + decimalFormat.format(particles.meanFreeEnergy()));
		control.println("PV/NkT = " + decimalFormat.format(particles.meanPressure()));
	}

	public void writeData() {

		for (int i = 0; i < particles.maxRadius / particles.grBinWidth; i++) {
			particles.sizeDist[i] = particles.sizeDist[i]/ ((particles.stop - particles.delay) / particles.snapshotInterval); // normalize size distribution
		}

		particles.volFrac = particles.volFrac / ((particles.stop - particles.delay) / particles.snapshotInterval); // average system volume fraction in equilibrium

		// write system parameters to a file in the data subdirectory
		try {
			File systemInfo = new File("data/systemInfo" + particles.fileExtension + ".txt");

			if (!systemInfo.exists()) { // if file doesn't exist, create it
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
			// bw.write("Dry volume fraction: " + particles.dryVolFrac);
			// bw.newLine();
			// bw.write("Reservoir volume fraction: " + particles.reservoirVolFrac);
			// bw.newLine();
			// bw.write("Actual volume fraction: " + particles.volFrac);
			// bw.newLine();
			bw.write("DryVolFrac increment: " + particles.dphi);
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
			for (int i = 0; i < particles.numberBins; i++) {
				bwrite.write(i + " " + particles.sizeDist[i]);
				bwrite.newLine();
			}

			bwrite.close();

		}

		catch (IOException e) {
			e.printStackTrace();
		}

		try {
			//File outputFile = new File("data/HertzianSpheres"+particles.fileExtension+".txt");
			File outputFile = new File("data/Free_Energy_Data/Interpenetration_Method_No_Size_Changes_Short_Run"+particles.fileExtension+".txt");

			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}

			FileWriter fw1 = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);

			// write system parameters to the file
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
			// bw1.write("Number of chains: " + particles.nChains);
			// bw1.newLine();
			bw1.write("Flory interaction parameter (chi): " + particles.chi);
			bw1.newLine();
			bw1.write("Young's calibration factor: " + particles.Young);
			bw1.newLine();
			bw1.write("x-link fraction: " + particles.xLinkFrac);
			bw1.newLine();
			bw1.write("DryVolFrac increment: " + particles.dphi);
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
			bw1.write("phi0,	mu/KT,		(PV/NkT)_total,		1+(PV/NkT)_virial,		 1+(PV/NkT)_pair/calculated,		(PV/NkT)_FR,      <F_total>/V,      EntropySolid,		Freference/V,           <F_FR>/V,     <alpha>,      zeta, 		newHertzianPotential, 		SpringConstant, 		lindemannParameter");
			bw1.newLine();

			for (int i = 1; i < dryVolFracs.size() - 1; i++) {
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				bw1.write(roundedDryVolFrac + ", " + chemicalPotentialList.get(i) + ", " + calculatedPressures.get(i) + ", " + + meanPressures.get(i) +  ", " + pairPressuresListEdited.get(i) + ", " + floryRehnerPressuresListEdited.get(i) + ", " + totalSumOfEnergiesList.get(i) + ", " + (uPairPerVolList.get(i)-totalSumOfEnergiesList.get(i)+1.50) + ", " + referenceFRPerVolList.get(i) + ", " +  floryFEperVolList.get(i) + ", " + swellingRatioList.get(i) + ", " + reservoirVolFracList.get(i) + ", " + newHertzianPotentialList.get(i) + ", " + springConstantList.get(i) + ", " + lindemannParameterList.get(i));
				bw1.newLine();
			}
			bw1.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// write radial distribution function, static structure factor to files in data
		// subdirectory

		if (structure) {
			rdf.writeRDF();
			ssf.writeSSF();
		}

	}

	/**
	 * Start the Java application.
	 * 
	 * @param args
	 * command line parameters
	 */
	public static void main(String[] args) { // set up animation control
			@SuppressWarnings("unused")
			SimulationControl control = SimulationControl.createApp(new HertzSpheresSolidPhaseFreeEnergyApp());
	}
}
