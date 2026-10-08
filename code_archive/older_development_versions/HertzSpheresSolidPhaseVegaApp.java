package org.opensourcephysics.sip.Hertz;

import java.io.File;
import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;
import java.text.DecimalFormat;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;

// imported for the purpose of creating an array to store the values of the dryVolFrac and free energies
import java.util.List;
import java.util.ResourceBundle.Control;
import java.util.ArrayList;

/**
 * HertzSpheresSolidPhaseDevelopmentApp
 *
 * Purpose:
 *   This program computes the free energy and associated thermodynamic parameters 
 *   of soft microgel systems using both the interpenetration and facet algorithms. 
 *   Monte Carlo simulations are performed for particles interacting via the Hertzian 
 *   elastic potential, with swelling behavior modeled by Flory–Rehner theory. 
 *   Free energy is evaluated using the Frenkel–Ladd method adapted for deformable particles.
 *
 * Algorithms Supported:
 *   - Interpenetration algorithm: accounts for overlap-dependent mixing free energy
 *   - Facet algorithm: accounts for surface contacts and volume caps between particles
 *
 * Key Features:
 *   - Iterative adjustment of dry volume fraction to explore bulk phase behavior
 *   - Gauss–Legendre quadrature integration for Frenkel–Ladd free energy calculation
 *   - Calculation of:
 *       • Pairwise interaction free energy
 *       • Flory–Rehner free energy
 *       • Helmholtz free energy per volume
 *       • Lindemann parameter for melting criterion
 *       • Pressure (Flory, pair, and total contributions)
 *       • Chemical potential
 *   - Structural analysis:
 *       • Radial distribution function g(r)
 *       • Static structure factor S(k)
 *   - Real-time 3D visualization of particle configurations
 *
 * Authors: Alan Denton and Oreoluwa Alade
 * Last Modified: 2026-08-25
 */

public class HertzSpheresSolidPhaseVegaApp extends AbstractSimulation {
	public enum WriteModes {WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};

	/* For the interpenetration algorithm */
	//HertzSpheresInterpenetration particles = new HertzSpheresInterpenetration(); // For the optimized interpenetration algorithm
	/* For the facet algorithm */
	// HertzSpheresFCCFacetFreeEnergy particles = new HertzSpheresFCCFacetFreeEnergy(); // For the facet algorithm with free energies
	
	HertzSpheresFCCFacetFreeEnergyVega particles = new HertzSpheresFCCFacetFreeEnergyVega();
	
	PlotFrame energyData = new PlotFrame("MC steps", "<E_pair>/N", "Mean pair energy per particle");
	PlotFrame pressureData = new PlotFrame("MC steps", "PV/NkT", "Mean pressure");
	PlotFrame sizeData = new PlotFrame("MC steps", "alpha", "Mean swelling ratio");
	Display3DFrame display3d = new Display3DFrame("Simulation animation");
	int weightIteration = 0, pointIteration = 0; // set and reset the iterations
	ElementSphere nanoSphere[];
	boolean added = false;
	boolean structure;
	RDF rdf;
	SSF ssf;
	double dryVolFracStart;
	double dryVolFracMax;
	double dryVolFrac;
	double lambda;
	double totalVol;
	boolean setLambda = true;
	boolean incrementDryVolFrac = true;
	double variableChanged;
	double deltaA2Accumulator = 0;
	double deltaA2PerVol;
	double gaussPoint;
	double gaussWeight;
	double fPair;
	double fPairPerVol;
	double totalFreeEnergy;
	double floryFEperVol;
	double initialFreeEnergy;
	double latestFreeEnergy;
	double floryPressure;
	double dFRFE_dPhi;
	double floryPContribution;
	double dFPair_dPhi;
	double pairPressure;
	double pairPContribution;
	double chemicalPotential;
	double uPairPerVol;
	double density;
	double nnDistance;
	double newHertzianPotential;
	boolean gaussLegendrePoint = true;
	boolean gaussLegendreWeight = true;
	double springConstant;
	double referenceFR;
	double referenceFRPerN;
	double referenceFRPerVol;
	double deltaA1;
	double deltaA1PerVol;
	double lindemannParameter;
	double maxValidPhi = 0.9641;
	boolean deltaA1Done = false;
	String outputDirectory = "data/Solid_Phase_Free_Energy/Vega_Implementation/xlink-3e-5/";

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
	List<Double> volumefractionList = new ArrayList<>();
	List<Double> softnessList = new ArrayList<>();
	List<Double> totalPressureList = new ArrayList<>();
	List<Double> einsteinPressureList = new ArrayList<>();
	List<Double> fPairPerVolList = new ArrayList<>();
	List<Double> virialPlusThermoList = new ArrayList<>();
	List<Double> deltaA1PerVolList = new ArrayList<>();
	List<Double> deltaA2PerVolList = new ArrayList<>();

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
		//particles.springConstant = control.getDouble("Spring constant"); // the spring constant
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
		particles.fileExtension = control.getString("Run number");
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

		if (springConstant > 0.0) {
			particles.springConstant = springConstant;
		}

		if (gaussLegendreWeight){
			// initialize the gaussLegendreWeights list
			gaussLegendreWeights.add(0.0666713443086881);
			gaussLegendreWeights.add(0.149451349150581);
			gaussLegendreWeights.add(0.219086362515982);
			gaussLegendreWeights.add(0.269266719309996);
			gaussLegendreWeights.add(0.295524224714753);
			gaussLegendreWeights.add(0.295524224714753);
			gaussLegendreWeights.add(0.269266719309996);
			gaussLegendreWeights.add(0.219086362515982);
			gaussLegendreWeights.add(0.149451349150581);
			gaussLegendreWeights.add(0.0666713443086881);
			gaussLegendreWeight = false;
		}

		if (gaussLegendrePoint){
			//initialize the gaussLegendrePoints list
			gaussLegendrePoints.add(-0.973906528517172);
			gaussLegendrePoints.add(-0.865063366688985);
			gaussLegendrePoints.add(-0.679409568299024);
			gaussLegendrePoints.add(-0.433395394129247);
			gaussLegendrePoints.add(-0.148874338981631);
			gaussLegendrePoints.add( 0.148874338981631);
			gaussLegendrePoints.add( 0.433395394129247);
			gaussLegendrePoints.add( 0.679409568299024);
			gaussLegendrePoints.add( 0.865063366688985);
			gaussLegendrePoints.add( 0.973906528517172);
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

				if (particles.calculateDeltaA1) {

					// Vega Step 1
					double logMeanBoltzmannFactor =
						particles.logMeanBoltzmannFactor();

					deltaA1 =
						particles.initialEnergy
						- logMeanBoltzmannFactor;

					deltaA1PerVol =
						deltaA1 / particles.totalVol;

					control.println("=== Vega Step 1 ===");
					control.println("U_lattice = " + particles.initialEnergy);
					control.println("H_real_lattice = " + particles.initialRealEnergy);
					control.println("log <exp[-(U_H-U_lattice)]> = "
						+ logMeanBoltzmannFactor);
					control.println("deltaA1 = " + deltaA1);
					control.println("deltaA1PerVol = " + deltaA1PerVol);

					deltaA1Done = true;
					particles.calculateDeltaA1 = false;

					// Now begin Vega Step 2
					gaussPoint = gaussLegendrePoints.get(pointIteration);
					gaussWeight = gaussLegendreWeights.get(weightIteration);

					variableChanged =
						0.5 * (
							gaussPoint
							* (Math.log(springConstant + Math.exp(3.5)) - 3.5)
							+ 3.5
							+ Math.log(springConstant + Math.exp(3.5))
						);

					lambda =
						(Math.exp(variableChanged) - Math.exp(3.5))
						/ particles.springConstant;

					particles.lambda = lambda;

					this.initialize();
					return;
				}

				if (lambda == 1) { // the real or interacting solid

					/* Compute the spring constant at lambda = 1*/
					springConstant = 3/(2*particles.meanSquareDisplacement()); // the springConstant
					control.println("springConstant: "+springConstant);
					control.println("meanSquareDisplacement at lambda=1: " + particles.meanSquareDisplacement());

					control.println("springConstant: "+springConstant);

					particles.springConstant = springConstant;
					
					/* Calculate the Flory–Rehner free energy for a fully interacting system */
					floryFEperVol = (particles.meanFreeEnergy()*particles.N)/particles.totalVol; // Flory-RehnerFR/vol

					control.println("floryFEperVol: "+floryFEperVol);

					control.println("mixFRSR = " + particles.mixFRSR);
					control.println("elasticFRSR = " + particles.elasticFRSR);
					control.println("meanFR before subtraction = " + (particles.freeEnergyAccumulator/particles.N/particles.numberOfConfigurations + particles.totalFRSR));

					control.println("reservoirSR = " + particles.reservoirSR);
    				control.println("totalFRSR per particle = " + (particles.mixFRSR + particles.elasticFRSR));
    				control.println("mean FR per particle before subtraction = " + (particles.freeEnergyAccumulator/particles.N/particles.numberOfConfigurations + particles.mixFRSR + particles.elasticFRSR));

					// the reference free energy (the free energy of the ideal Einstein Crystal)
					referenceFR = -3*(particles.N-1)/2.0*Math.log(Math.PI/particles.springConstant)-Math.log(particles.N)/2.0-Math.log(particles.totalVol);
					
					double A0Spring =
						-3.0*(particles.N-1)/2.0
						* Math.log(Math.PI/particles.springConstant);

					double A0N =
						-0.5*Math.log(particles.N);

					double A0Volume =
						-Math.log(particles.totalVol);

					referenceFRPerN = referenceFR / particles.N;
					referenceFRPerVol = referenceFR / particles.totalVol;

					control.println("=== Vega A0 check ===");
					control.println("A0 spring term = " + A0Spring);
					control.println("A0 N term = " + A0N);
					control.println("A0 volume term = " + A0Volume);
					control.println("A0 total = " + referenceFR);
					control.println("A0/V = " + referenceFRPerVol);

					density = particles.N/particles.totalVol;
					nnDistance = (1/Math.sqrt(2))*Math.pow((4/density), 1/3.0); // the nearest neighbour distance
					
					// the lindemannParameter
					lindemannParameter = Math.sqrt(particles.meanSquareDisplacement())/nnDistance;
					control.println("lindemannParameter: "+lindemannParameter);

					// as recommended by Zacahareli and co
					if (nnDistance>(2*particles.meanRadius())){
						newHertzianPotential = 0;
					}
					else{
						newHertzianPotential = particles.B*Math.pow((1-nnDistance/(2*particles.meanRadius())), 2.5);
					}

					// control.println("newHertzianPotential: "+newHertzianPotential);
					
					uPairPerVol = particles.meanPairEnergy()*(particles.N/particles.totalVol); // the pair energy per volume
					
					control.println("uPairPerVol: "+uPairPerVol);

					// After the physical solid run, perform Vega Step 1
					particles.calculateDeltaA1 = true;

					this.initialize();
					return;
				}
				
				else{
					
					deltaA2Accumulator +=
						(0.5 * gaussWeight)
						* (Math.log(particles.springConstant + Math.exp(3.5)) - 3.5)
						* Math.exp(variableChanged)
						/ particles.springConstant
						* (-particles.meanSpringEnergy());

					control.println(
						"GL point " + (pointIteration + 1)
						+ ": lambda = " + lambda
						+ ", spring = " + particles.meanSpringEnergy()
						+ ", DeltaA2 integrand = "
						+ (-particles.meanSpringEnergy())
					);
										
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
						control.println("--- Gauss-Legendre point " + pointIteration + " ---");
						control.println("lambda = " + lambda);
						control.println("meanPairEnergy = " + particles.meanPairEnergy());
						control.println("meanFreeEnergy = " + particles.meanFreeEnergy());
						control.println("meanSpringEnergy = " + particles.meanSpringEnergy());
						control.println("meanSquareDisplacement = " + particles.meanSquareDisplacement());

						deltaA2PerVol =
							deltaA2Accumulator
							* (particles.N / particles.totalVol); //ΔF/V
						

						// the total free energy
						totalFreeEnergy =
							referenceFRPerVol
							+ deltaA1PerVol
							+ deltaA2PerVol;
						control.println("totalF: " + totalFreeEnergy);

						control.println("=== phi0 = " + dryVolFrac + " ===");
						control.println("initialEnergy = " + particles.initialEnergy);
						control.println("referenceFRPerVol = " + referenceFRPerVol);
						control.println("deltaA2PerVol = " + deltaA2PerVol);
						control.println("floryFEperVol = " + floryFEperVol);
						control.println("totalFreeEnergy solid = " + totalFreeEnergy);

						double phi = particles.meanVolFrac();

						control.println("phi: " + phi);
						control.println("softness: " + (1.0 / particles.B));

						// Stop if the facet model has reached its validity limit
						if (phi >= maxValidPhi) {
							control.println("Facet-model validity limit reached.");
							control.println("Stopping because phi = " + phi
								+ " >= " + maxValidPhi);

							writeData();
							return;
						}

						// Update the array lists only for valid states
						dryVolFracs.add(dryVolFrac);
						uPairPerVolList.add(uPairPerVol);
						totalVolList.add(particles.totalVol);
						reservoirVolFracList.add(particles.reservoirVolFrac);
						swellingRatioList.add(particles.meanRadius());

						referenceFRPerVolList.add(referenceFRPerVol); // A0/V
						deltaA1PerVolList.add(deltaA1PerVol);         // ΔA1/V
						deltaA2PerVolList.add(deltaA2PerVol);         // ΔA2/V
						totalSumOfEnergiesList.add(totalFreeEnergy);  // Asol/V

						volumefractionList.add(phi);
						softnessList.add(1.0 / particles.B);
						springConstantList.add(springConstant);
						lindemannParameterList.add(lindemannParameter);
						newHertzianPotentialList.add(newHertzianPotential);

						// increment the dry volume fraction
						dryVolFrac += particles.dphi;
						particles.dryVolFrac = dryVolFrac;
						//reset lambda to 1 for the interacting solid
						lambda = 1;
						particles.lambda = lambda;
						deltaA2Accumulator = 0; //reset UPairMinusUSpring
						//reset the iterations
						pointIteration = 0;
						weightIteration = 0;
						deltaA1 = 0.0;
						deltaA1PerVol = 0.0;
						deltaA1Done = false;
						particles.calculateDeltaA1 = false;

						this.initialize();
						return;
					}
				}
				
			}

			if (false && dryVolFracs.size() > 2) {
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
				calculatedPressures.add(0, 0.0);
				chemicalPotentialList.add(0, 0.0);
				totalPressureList.add(0, 0.0);  
				einsteinPressureList.add(0, 0.0);   
				fPairPerVolList.add(0, 0.0);   
				virialPlusThermoList.add(0, 0.0); 
								
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

					// correct pressure from total free energy derivative
					double dFtotal_dPhi = (totalSumOfEnergiesList.get(i+1) - totalSumOfEnergiesList.get(i-1)) 
										/ (2.0 * particles.dphi);
					double totalPressureContribution = (dFtotal_dPhi * currentDryVolFrac - totalSumOfEnergiesList.get(i)) 
													* totalVolList.get(i) / particles.N;
					totalPressureList.add(totalPressureContribution); // store it

					double dFref_dPhi = (referenceFRPerVolList.get(i+1) - referenceFRPerVolList.get(i-1))
                  / (2.0 * particles.dphi);
					double einsteinPressure = (dFref_dPhi * currentDryVolFrac - referenceFRPerVolList.get(i))
								* totalVolList.get(i) / particles.N;
					double virialPlusThermo = meanPressures.get(i) + totalPressureContribution;

					einsteinPressureList.add(einsteinPressure);
					fPairPerVolList.add(pairFreeEnergyList.get(i));
					virialPlusThermoList.add(virialPlusThermo);

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

		// Density scan
		control.setValue("DryVolFracStart", 0.0020);
		control.setValue("DryVolFrac Max", 0.0030001);
		control.setValue("DryVolFrac increment", 1.1111111111111112E-4);

		// System
		control.setValue("Initial configuration", "FCC");
		control.setValue("N", 108);

		// Microgel parameters
		control.setValue("Dry radius [nm]", 50);
		control.setValue("x-link fraction", 0.00003);
		control.setValue("Young's calibration", 1.0);
		control.setValue("chi", 0);

		// Monte Carlo
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.0);
		control.setValue("Delay", 7000);
		control.setValue("Snapshot interval", 100);
		control.setValue("Stop", 70000);

		// Structural analysis
		control.setValue("Maximum radial distance", 10);
		control.setValue("Size bin width", 0.001);
		control.setValue("g(r) bin width", 0.005);
		control.setValue("Delta k", 0.005);
		control.setValue("Calculate structure", false);

		// Output
		control.setValue("Run number", "1");

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
		// Normalize distributions (unchanged)
		for (int i = 0; i < particles.maxRadius / particles.grBinWidth; i++) {
			particles.sizeDist[i] = particles.sizeDist[i] / ((particles.stop - particles.delay) / particles.snapshotInterval);
		}
		particles.volFrac = particles.volFrac / ((particles.stop - particles.delay) / particles.snapshotInterval);

		// 1. Write system parameters to data/systemInfo*.txt
		try {
				File systemInfo = new File(
					outputDirectory
					+ "systemInfo_run"
					+ particles.fileExtension
					+ ".txt"
			);

			File systemDir = systemInfo.getParentFile();  // "data/"
			if (systemDir != null && !systemDir.exists()) {
				systemDir.mkdirs();
			}
			if (!systemInfo.exists()) {
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
			bw.write("Mean pressure PV/NkT: " + particles.meanPressure());
			bw.newLine();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// 2. Write size distribution to data/microgelSize*.txt
		try {
			File sizeFile = new File(
				outputDirectory
				+ "microgelSize_run"
				+ particles.fileExtension
				+ ".txt"
			);
			
			File sizeDir = sizeFile.getParentFile();  // "data/"
			if (sizeDir != null && !sizeDir.exists()) {
				sizeDir.mkdirs();
			}
			if (!sizeFile.exists()) {
				sizeFile.createNewFile();
			}

			FileWriter fwrite = new FileWriter(sizeFile.getAbsoluteFile());
			BufferedWriter bwrite = new BufferedWriter(fwrite);
			for (int i = 0; i < particles.numberBins; i++) {
				bwrite.write(i + " " + particles.sizeDist[i]);
				bwrite.newLine();
			}
			bwrite.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Write solid-phase free-energy results
		try {
				File outputFile = new File(
					outputDirectory
					+ "FacetSolid_run"
					+ particles.fileExtension
					+ ".txt"
			);

			System.out.println("Output file path: " + outputFile.getAbsolutePath());  // Debug log

			File outputDir = outputFile.getParentFile();
			if (outputDir != null && !outputDir.exists()) {
				outputDir.mkdirs();
			}
			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}

			FileWriter fw1 = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);

			// System params header
			bw1.write("Starting dry volume fraction: " + dryVolFracStart);
			bw1.newLine();
			bw1.write("Maximum dry volume fraction: " + dryVolFracMax);
			bw1.newLine();
			bw1.write("DryVolFrac increment: " + particles.dphi);
			bw1.newLine();
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

			bw1.write(
				"phi0, phi, A0_per_V, deltaA1_per_V, deltaA2_per_V, "
				+ "Asol_per_V, SpringConstant, LindemannParameter, kT_over_E"
			);
			bw1.newLine();			
			// Data rows (with NaN guard for safety)
			for (int i = 0; i < dryVolFracs.size(); i++) {

				double roundedDryVolFrac =
					Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));

				bw1.write(
					roundedDryVolFrac + ", "
					+ volumefractionList.get(i) + ", "
					+ referenceFRPerVolList.get(i) + ", "
					+ deltaA1PerVolList.get(i) + ", "
					+ deltaA2PerVolList.get(i) + ", "
					+ totalSumOfEnergiesList.get(i) + ", "
					+ springConstantList.get(i) + ", "
					+ lindemannParameterList.get(i) + ", "
					+ softnessList.get(i)
				);

				bw1.newLine();
			}
			bw1.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Write RDF/SSF if enabled (unchanged)
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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresSolidPhaseVegaApp());
	}
}
