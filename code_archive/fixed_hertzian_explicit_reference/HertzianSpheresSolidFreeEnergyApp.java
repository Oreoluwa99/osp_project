package org.opensourcephysics.sip.Hertz;

import java.io.File;
import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;
import java.text.DecimalFormat;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.display3d.simple3d.ElementSphere;

// imported for the purpose of creating an array to store the values of the dryVolFrac and free energies
import java.util.List;
import java.util.ArrayList;
import java.util.Locale;

/**
 * HertzianSpheresSolidFreeEnergyApp.java
 *
 * Purpose:
 *   Computes the FCC Helmholtz free energy of fixed-size Hertzian spheres by
 *   Frenkel–Ladd thermodynamic integration to a centre-of-mass-constrained
 *   Einstein crystal.
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
 * Version: 1.2 (April 2025)
 */

public class HertzianSpheresSolidFreeEnergyApp extends AbstractSimulation {
	public enum WriteModes {WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
	//HertzSpheresSolidPhase particles = new HertzSpheresSolidPhase();

	/* For the interpenetration algorithm */
	//HertzSpheresInterpenetration particles = new HertzSpheresInterpenetration(); // For the optimized interpenetration algorithm

	/* For the facet algorithm */
	// HertzSpheresFCCFacetFreeEnergy particles = new HertzSpheresFCCFacetFreeEnergy(); // For the facet algorithm with free energies
	HertzianSpheresModelSolidPhase particles = new HertzianSpheresModelSolidPhase(); // For the facet algorithm with free energies

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
	double dryVolFracStart;
	double dryVolFracMax;
	double dryVolFrac;
	double lambda;
	double totalVol;
	boolean setLambda = true;
	boolean incrementDryVolFrac = true;
	double variableChanged;
	double einsteinFreeEnergy;
	double deltaFAccumulator = 0;
	double gaussPoint;
	double gaussWeight;
	double fPair;
	double fPairPerVol;
	double totalFreeEnergy;
	double floryFEperVol;
	double initialFreeEnergy;
	double latestFreeEnergy;
	double initialVolume;
	double latestVolume;
	double UPairMinusUSpringPerVolume;
	double floryPressure;
	double dFRFE_dPhi;
	double floryPContribution;
	double initialPairFreeEnergy;
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
	double deltaFPerVol;
	double springConstant;
	double referenceFR;
	double referenceFRPerN;
	double referenceFRPerVol;
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
	List<Double> volumefractionList = new ArrayList<>();
	List<Double> softnessList = new ArrayList<>();
	List<Double> totalPressureList = new ArrayList<>();
	List<Double> einsteinPressureList = new ArrayList<>();
	List<Double> fPairPerVolList = new ArrayList<>();

	/*
	 * Temperatures used for the fluid--FCC part of the Hertzian-sphere phase
	 * diagram.  T* = 1/(beta epsilon).  These density windows are deliberately
	 * broad pilot windows; after locating the crossing, only the neighbouring
	 * densities need a finer follow-up run.
	 */
	final double[] betaEpsilonValues = {180.0, 160.0, 150.0, 140.0, 135.0, 130.0};
	final double[] solidDensityStarts = {1.50, 1.55, 1.60, 1.65, 1.70, 1.75};
	final double[] solidDensityMaxes = {1.801, 1.901, 1.951, 2.001, 2.051, 2.101};
	int temperatureIndex = 0;
	boolean temperatureSeriesComplete = false;

	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values to 3 dp

	/**
	 * Initializes the model.
	 */
	public void initialize() {

		added = false;
		//particles.dlambda = control.getDouble("Lambda increment");// the coupling constant increment
		dryVolFracStart = solidDensityStarts[temperatureIndex];
		dryVolFracMax = solidDensityMaxes[temperatureIndex];
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

		particles.sigma = 1.0;
		particles.betaEpsilon = betaEpsilonValues[temperatureIndex];
		particles.density = dryVolFrac;

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
		if (temperatureSeriesComplete) {
			return;
		}
		particles.step();

		// if initial configuration is random, no particle interactions for delay/10

		if (particles.steps <= particles.stop) {
			if (particles.steps > particles.delay) {
				if ((particles.steps - particles.delay) % particles.snapshotInterval == 0) {
					particles.sizeDistribution();
					if (structure) { // accumulate statistics for structural properties
						rdf.update(); // g(r)
						ssf.update(); // S(k)
					}
					particles.volFracAccumulator += particles.calculateVolumeFraction();
				}
			}
		}

		if (particles.steps == particles.stop) {
			System.out.println("DryVolFrac = " + dryVolFrac + ", dryVolFracMax = " + dryVolFracMax);
			if (dryVolFrac < dryVolFracMax) {
				control.println("rho*: "+dryVolFrac);
				control.println("xLinkFrac: "+ particles.xLinkFrac);
				if (lambda == 1) { // the real or interacting solid

					/* Compute the spring constant at lambda = 1*/
					springConstant = 3/(2*particles.meanSquareDisplacement()); // the springConstant
					control.println("springConstant: "+springConstant);
					control.println("meanSquareDisplacement at lambda=1: " + particles.meanSquareDisplacement());

					control.println("springConstant: "+springConstant);

					particles.springConstant = springConstant;
					springConstantList.add(springConstant); // update the list
					
					/* Flory–Rehner energy is absent from the fixed-size Hertzian benchmark. */
					floryFEperVol = 0.0;
					floryFEperVolList.add(floryFEperVol); // update the list

					control.println("floryFEperVol: "+floryFEperVol);

					// control.println("mixFRSR = " + particles.mixFRSR);
					// control.println("elasticFRSR = " + particles.elasticFRSR);
					// control.println("meanFR before subtraction = " + (particles.freeEnergyAccumulator/particles.N/particles.numberOfConfigurations + particles.totalFRSR));

					// control.println("reservoirSR = " + particles.reservoirSR);
    				// control.println("totalFRSR per particle = " + (particles.mixFRSR + particles.elasticFRSR));
    				// control.println("mean FR per particle before subtraction = " + (particles.freeEnergyAccumulator/particles.N/particles.numberOfConfigurations + particles.mixFRSR + particles.elasticFRSR));

					// the reference free energy (the free energy of the ideal Einstein Crystal)
					referenceFR = -3*(particles.N-1)/2.0*Math.log(Math.PI/particles.springConstant)-Math.log(particles.N)/2.0-Math.log(particles.totalVol);
					referenceFRPerN = referenceFR/particles.N; // reference free energy per particles
					referenceFRPerVol = (referenceFRPerN*particles.N)/particles.totalVol; // reference free energy per volume
					referenceFRPerVolList.add(referenceFRPerVol); // update the list

					meanPressures.add(particles.meanPressure()); // pressure from virial theorem
					density = particles.N/particles.totalVol;
					nnDistance = (1/Math.sqrt(2))*Math.pow((4/density), 1/3.0); // the nearest neighbour distance
					
					// the lindemannParameter
					lindemannParameter = Math.sqrt(particles.meanSquareDisplacement())/nnDistance;
					lindemannParameterList.add(lindemannParameter);
					control.println("lindemannParameter: "+lindemannParameter);

					// as recommended by Zacahareli and co
					if (nnDistance >= particles.sigma) {
						newHertzianPotential = 0.0;
					}
					else {
						newHertzianPotential = particles.betaEpsilon
							* Math.pow(1.0 - nnDistance / particles.sigma, 2.5);
					}

					// control.println("newHertzianPotential: "+newHertzianPotential);

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
						control.println("--- Gauss-Legendre point " + pointIteration + " ---");
						control.println("lambda = " + lambda);
						control.println("meanPairEnergy = " + particles.meanPairEnergy());
						control.println("meanSpringEnergy = " + particles.meanSpringEnergy());
						control.println("meanSquareDisplacement = " + particles.meanSquareDisplacement());
						control.println("meanPairEnergy - meanSpringEnergy = " + (particles.meanPairEnergy() - particles.meanSpringEnergy()));
				
						deltaFPerVol = deltaFAccumulator*(particles.N/particles.totalVol); //ΔF/V
						
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
						volumefractionList.add(particles.meanVolFrac()); // the volume fraction list
						control.println("phi: " + particles.meanVolFrac());
						double reducedTemperature =
							1.0 / particles.betaEpsilon;

						softnessList.add(reducedTemperature);
						control.println(
							"T* = " + reducedTemperature
						);

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
				calculatedPressures.add(0, 0.0);
				chemicalPotentialList.add(0, 0.0);
				totalPressureList.add(0, 0.0);  
				einsteinPressureList.add(0, 0.0);   
				fPairPerVolList.add(0, 0.0);   
								
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

					double densitySpacing = nextDryVolFrac - previousDryVolFrac;

					// For f = beta F sigma^3/V, beta mu = df/drho*.
					chemicalPotential =
						(latestFreeEnergy - initialFreeEnergy)
						/ densitySpacing;

					// Component derivatives with respect to rho*.
					dFRFE_dPhi = (nextFloryFreeEnergy - previousFloryFreeEnergy) / densitySpacing;
					dFPair_dPhi = (nextPairFreeEnergy - previousPairFreeEnergy) / densitySpacing;

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
										/ densitySpacing;
					double totalPressureContribution = (dFtotal_dPhi * currentDryVolFrac - totalSumOfEnergiesList.get(i)) 
													* totalVolList.get(i) / particles.N;
					totalPressureList.add(totalPressureContribution); // store it

					double dFref_dPhi = (referenceFRPerVolList.get(i+1) - referenceFRPerVolList.get(i-1))
								  / densitySpacing;
					double einsteinPressure = (dFref_dPhi * currentDryVolFrac - referenceFRPerVolList.get(i))
								* totalVolList.get(i) / particles.N;
					einsteinPressureList.add(einsteinPressure);
					fPairPerVolList.add(pairFreeEnergyList.get(i));

					/* Update the lists accordingly */
					floryRehnerPressuresListEdited.add(floryPContribution);
					pairPressuresListEdited.add(pairPContribution); // the F0 (eq. 48 from Vegas et. al.) already incldues the ideal gas free energy
					calculatedPressures.add(floryPContribution + pairPContribution); // the calculated pressure from the derivatives of the free energies
					chemicalPotentialList.add(chemicalPotential);
				}
			}
			writeData();
			if (advanceToNextTemperature()) {
				this.initialize();
				return;
			}
			temperatureSeriesComplete = true;
			control.println("Completed all fluid--FCC temperatures.");
			return;
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

	private String currentOutputDirectory() {
		double betaEpsilon = betaEpsilonValues[temperatureIndex];
		String betaLabel = String.format(Locale.US, "%.0f", betaEpsilon);
		String temperatureLabel = String.format(Locale.US, "%.8f", 1.0 / betaEpsilon)
			.replaceFirst("0+$", "")
			.replaceFirst("\\.$", "");
		return "data/Hertzian_Spheres_Benchmark/betaEpsilon_"
			+ betaLabel + "_Tstar_" + temperatureLabel + "/solid/";
	}

	private boolean advanceToNextTemperature() {
		if (temperatureIndex + 1 >= betaEpsilonValues.length) {
			return false;
		}

		temperatureIndex++;
		clearTemperatureData();
		dryVolFracStart = solidDensityStarts[temperatureIndex];
		dryVolFracMax = solidDensityMaxes[temperatureIndex];
		dryVolFrac = dryVolFracStart;
		particles.dryVolFrac = dryVolFrac;
		particles.density = dryVolFrac;
		particles.betaEpsilon = betaEpsilonValues[temperatureIndex];
		lambda = 1.0;
		particles.lambda = lambda;
		incrementDryVolFrac = false;

		control.println("");
		control.println("Starting betaEpsilon = " + particles.betaEpsilon
			+ ", T* = " + (1.0 / particles.betaEpsilon));
		control.println("FCC density window = " + dryVolFracStart
			+ " to " + dryVolFracMax);
		return true;
	}

	private void clearTemperatureData() {
		dryVolFracs.clear();
		dryVolFracsEdited.clear();
		floryFEperVolList.clear();
		totalSumOfEnergiesList.clear();
		totalVolList.clear();
		calculatedPressures.clear();
		meanPressures.clear();
		pairFreeEnergyList.clear();
		floryRehnerPressuresList.clear();
		pairPressuresList.clear();
		floryRehnerPressuresListEdited.clear();
		pairPressuresListEdited.clear();
		chemicalPotentialList.clear();
		reservoirVolFracList.clear();
		swellingRatioList.clear();
		uPairPerVolList.clear();
		newHertzianPotentialList.clear();
		referenceFRPerVolList.clear();
		springConstantList.clear();
		lindemannParameterList.clear();
		volumefractionList.clear();
		softnessList.clear();
		totalPressureList.clear();
		einsteinPressureList.clear();
		fPairPerVolList.clear();
		deltaFAccumulator = 0.0;
		pointIteration = 0;
		weightIteration = 0;
	}

	/**
	 * Resets the model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		temperatureIndex = 0;
		temperatureSeriesComplete = false;
		clearTemperatureData();
		incrementDryVolFrac = true;
		setLambda = true;
		lambda = 1.0;

		control.setValue("DryVolFracStart", solidDensityStarts[0]);
		control.setValue("DryVolFrac Max", solidDensityMaxes[0]);
		control.setValue("DryVolFrac increment", 0.020);
		control.setValue("Initial configuration", "FCC");
		// control.setValue("Spring constant", 10000); // Spring constant: 2.035 for alpha/KT = 100
		control.setValue("N", 108); // number of particles
		control.setValue("x-link fraction", 0.00005);
		// control.setValue("N", 500); for FCC lattice, N/4 should be a perfect cube
		control.setValue("Dry radius [nm]", 50);
		control.setValue("Young's calibration", 1.0); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0);
		control.setValue("Delay", 10000); // steps after which statistics collection starts
		control.setValue("Snapshot interval", 100); // steps separating successive samples
		control.setValue("Stop", 20000); // steps after which statistics collection stops
		control.setValue("Size bin width", .001); // bin width of particle radius histogram
		control.setValue("g(r) bin width", .005); // bin width of g(r) histogram
		control.setValue("Delta k", .005); // bin width of S(k) histogram
		control.setValue("File extension", "_phaseDiagram_run1");
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
		String outputDirectory = currentOutputDirectory();
		File temperatureDirectory = new File(outputDirectory);
		if (!temperatureDirectory.exists()) {
			temperatureDirectory.mkdirs();
		}

		// Normalize distributions (unchanged)
		for (int i = 0; i < particles.maxRadius / particles.grBinWidth; i++) {
			particles.sizeDist[i] = particles.sizeDist[i] / ((particles.stop - particles.delay) / particles.snapshotInterval);
		}
		particles.volFrac = particles.volFrac / ((particles.stop - particles.delay) / particles.snapshotInterval);

		// 1. Write system parameters to data/systemInfo*.txt
		try {
			File systemInfo = new File(
				outputDirectory + "systemInfo"
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
			bw.write("Reduced-density increment: " + particles.dphi);
			bw.newLine();
			bw.write("betaEpsilon: " + particles.betaEpsilon);
			bw.newLine();
			bw.write("Tstar: " + (1.0 / particles.betaEpsilon));
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
				outputDirectory + "microgelSize"
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

		// 3. Write Facet_data to data/APS_2026/Facet/Facet_data*.txt
		try {
			File outputFile = new File(
				outputDirectory + "Free_Energy_SolidPhase"
					+ particles.fileExtension
					+ ".txt"
			);
			System.out.println("Output file path: " + outputFile.getAbsolutePath());  // Debug log

			File outputDir = outputFile.getParentFile();  // "data/APS_2026/Facet/"
			if (outputDir != null && !outputDir.exists()) {
				outputDir.mkdirs();
			}
			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}

			FileWriter fw1 = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);

			// System params header
			bw1.write("Starting reduced density rho*: " + dryVolFracStart);
			bw1.newLine();
			bw1.write("Maximum reduced density rho*: " + dryVolFracMax);
			bw1.newLine();
			bw1.write("Reduced-density increment: " + particles.dphi);
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
			bw1.write("betaEpsilon: " + particles.betaEpsilon);
			bw1.newLine();
			bw1.write("Tstar: " + (1.0 / particles.betaEpsilon));
			bw1.newLine();
			bw1.write("Reduced-density increment: " + particles.dphi);
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

			bw1.write("rho, geometric_volume_fraction, mu_from_dF/kT, Z_from_dF_total, Pstar_from_dF_total, Z_virial_measured, Pstar_virial_measured, Z_einstein_contribution, Z_pair_contribution, Z_FR_contribution, <F_total>/V, <F_pair>/V, <S_solid>/V, <F_reference>/V, <F_FR>/V, radius, zeta, beta_u_Hertz_nn, SpringConstant, LindemannParameter, Tstar");
			bw1.newLine();			
			// Data rows (with NaN guard for safety)
			for (int i = 1; i < dryVolFracs.size() - 1; i++) {
				if (i >= chemicalPotentialList.size() || Double.isNaN(chemicalPotentialList.get(i))) {
					System.err.println("Warning: Skipping row " + i + " due to NaN/missing data.");
					continue;
				}
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				bw1.write(roundedDryVolFrac + ", " + 
					volumefractionList.get(i) + ", " + 
					chemicalPotentialList.get(i) + ", " + 
					totalPressureList.get(i) + ", " +
					(totalPressureList.get(i) * dryVolFracs.get(i)) + ", " +
					meanPressures.get(i) + ", " +
					(meanPressures.get(i) * dryVolFracs.get(i)) + ", " +
					einsteinPressureList.get(i) + ", " +
					pairPressuresListEdited.get(i) + ", " + 
					floryRehnerPressuresListEdited.get(i) + ", " + 
					totalSumOfEnergiesList.get(i) + ", " +
					fPairPerVolList.get(i) + ", " +
					(uPairPerVolList.get(i) - totalSumOfEnergiesList.get(i)
						+ 1.50 * dryVolFracs.get(i)) + ", " + 
					referenceFRPerVolList.get(i) + ", " + 
					floryFEperVolList.get(i) + ", " + 
					swellingRatioList.get(i) + ", " + 
					reservoirVolFracList.get(i) + ", " + 
					newHertzianPotentialList.get(i) + ", " + 
					springConstantList.get(i) + ", " + 
					lindemannParameterList.get(i) + ", " + 
					softnessList.get(i));
				bw1.newLine();
			}
			bw1.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			File rawFile = new File(
				outputDirectory + "FCC_FreeEnergyRaw"
				+ particles.fileExtension
				+ ".txt"
			);

			File rawDir = rawFile.getParentFile();
			if (rawDir != null && !rawDir.exists()) {
				rawDir.mkdirs();
			}

			BufferedWriter rawWriter =
				new BufferedWriter(new FileWriter(rawFile));

			rawWriter.write(
				"rho, F_total_perV, F_total_perN, F_pair_perV, F_reference_perV, "
				+ "Z_virial_measured, springConstant, lindemannParameter, Tstar"
			);
			rawWriter.newLine();

			for (int i = 0; i < dryVolFracs.size(); i++) {
				rawWriter.write(
					dryVolFracs.get(i) + ", "
					+ totalSumOfEnergiesList.get(i) + ", "
					+ (totalSumOfEnergiesList.get(i) / dryVolFracs.get(i)) + ", "
					+ pairFreeEnergyList.get(i) + ", "
					+ referenceFRPerVolList.get(i) + ", "
					+ meanPressures.get(i) + ", "
					+ springConstantList.get(i) + ", "
					+ lindemannParameterList.get(i) + ", "
					+ softnessList.get(i)
				);
				rawWriter.newLine();
			}

			rawWriter.close();
		}
		catch (IOException e) {
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
			SimulationControl control = SimulationControl.createApp(new HertzianSpheresSolidFreeEnergyApp());
	}
}
