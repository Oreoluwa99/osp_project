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
 * HertzSpheresApp performs a Monte Carlo simulation of microgels interacting
 * via the
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free
 * energy.
 * 
 * @authors Alan Denton and Matt Urich
 * @version 1.2 25-05-2021
 * 
 */

public class HertzSpheresSolidPhaseGaussLegendreApp extends AbstractSimulation {
	public enum WriteModes {WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
	HertzSpheresSolidPhaseCopy particles = new HertzSpheresSolidPhaseCopy();																			
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
	double x, y;
	double uLatticePerN;
	double fPair, fPairPerVol, totalFreeEnergy;
	double floryFEperVol;
	boolean isFirstIteration;
	double initialFreeEnergy, latestFreeEnergy, initialVolume, latestVolume;
	double calculatedPressure;
	double volumeChange, freeEnergyChange;
	double latestFloryRehnerFreeEnergy;
	double floryPressure, dFRFE_dPhi, floryPContribution, initialPairFreeEnergy, dFPair_dPhi, pairPressure, pairPContribution;
	double chemicalPotential;
	double uPairPerVol;
	double logBoltzmannFactor;
	double deltaF1MinusUlattice;
	double density, nnDistance, lindemannParameter;
	double newHertzianPotential;
	boolean gaussLegendrePoint = true;
	boolean gaussLegendreWeight = true;
	int weightIteration = 0, pointIteration = 0;
	double exponentialOfX, exponentialOfXAccumulator = 0, sinOfXAccumulator = 0, sinOfX, integralOfOneAccumulator = 0;


	/* Lists to store various calculated values */
	List<Double> dryVolFracs = new ArrayList<>();
	List<Double> delta2List = new ArrayList<>();
	List<Double> deltaF1List = new ArrayList<>();
	List<Double> dryVolFracsEdited = new ArrayList<>();
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
	List<Double> floryRehnerPressuresListEdited = new ArrayList<>();
	List<Double> pairPressuresListEdited = new ArrayList<>();
	List<Double> totalVarList = new ArrayList<>();
	List<Double> chemicalPotentialList = new ArrayList<>();
	List<Double> reservoirVolFracList = new ArrayList<>();
	List<Double> swellingRatioList = new ArrayList<>();
	List<Double> uPairPerVolList = new ArrayList<>();
	List<Double> logBoltzmannFactorList = new ArrayList<>();
	List<Double> deltaF1MinusUlatticeList = new ArrayList<>();
	List<Double> newHertzianPotentialList = new ArrayList<>();
	List<Double> gaussLegendreWeights = new ArrayList<>();
	List<Double> gaussLegendrePoints = new ArrayList<>();
	List<Double> exponentialOfXList = new ArrayList<>();
	List<Double> trigonmetricExamplesList = new ArrayList<>();
	List<Double> lindemannParameterList = new ArrayList<>();
	List<Double> oneList = new ArrayList<>();


	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values to 3 dp

	/**
	 * Initializes the model.
	 */
	public void initialize() {

		added = false;
		particles.dlambda = control.getDouble("Lambda increment");// the coupling constant increment
		dryVolFracStart = control.getDouble("DryVolFracStart");
		dryVolFracMax = control.getDouble("DryVolFrac Max");
		particles.dphi = control.getDouble("DryVolFrac increment");
		xIncrement = control.getDouble("limit increment");// to increment the limits
		particles.springConstant = control.getDouble("Spring constant"); // the spring constant
		// upperLimit = control.getDouble("upper limit");
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
			gaussLegendreWeights.add(0.01761401);
			gaussLegendreWeights.add(0.04060143);
			gaussLegendreWeights.add(0.06267205);
			gaussLegendreWeights.add(0.08327674);
			gaussLegendreWeights.add(0.10193012);
			gaussLegendreWeights.add(0.11819453);
			gaussLegendreWeights.add(0.13168864);
			gaussLegendreWeights.add(0.14209611);
			gaussLegendreWeights.add(0.14917299);
			gaussLegendreWeights.add(0.15275339);
			gaussLegendreWeights.add(0.15275339);
			gaussLegendreWeights.add(0.14917299);
			gaussLegendreWeights.add(0.14209611);
			gaussLegendreWeights.add(0.13168864);
			gaussLegendreWeights.add(0.11819453);
			gaussLegendreWeights.add(0.10193012);
			gaussLegendreWeights.add(0.08327674);
			gaussLegendreWeights.add(0.06267205);
			gaussLegendreWeights.add(0.04060143);
			gaussLegendreWeights.add(0.01761401);
			gaussLegendreWeight = false;
		}

		if (gaussLegendrePoint){
			//initialize the gaussLegendrePoints list
			gaussLegendrePoints.add(-0.9931286);
        	gaussLegendrePoints.add(-0.96397193);
        	gaussLegendrePoints.add(-0.91223443);
        	gaussLegendrePoints.add(-0.83911697);
        	gaussLegendrePoints.add(-0.74633191);
        	gaussLegendrePoints.add(-0.63605368);
        	gaussLegendrePoints.add(-0.510867);
        	gaussLegendrePoints.add(-0.37370609);
        	gaussLegendrePoints.add(-0.22778585);
        	gaussLegendrePoints.add(-0.07652652);
        	gaussLegendrePoints.add(0.07652652);
       	 	gaussLegendrePoints.add(0.22778585);
       		gaussLegendrePoints.add(0.37370609);
        	gaussLegendrePoints.add(0.510867);
        	gaussLegendrePoints.add(0.63605368);
        	gaussLegendrePoints.add(0.74633191);
        	gaussLegendrePoints.add(0.83911697);
        	gaussLegendrePoints.add(0.91223443);
        	gaussLegendrePoints.add(0.96397193);
        	gaussLegendrePoints.add(0.9931286);
			gaussLegendrePoint = false;
		}

		/*
		 * set the value of x and lambda from changing the variables in the integration
		 */
		
		if (setLambda) {
			x = gaussLegendrePoints.get(pointIteration); // set the value of x (for the change of variables)
			y = (Math.log(particles.springConstant+Math.exp(3.5))+3.5)/2 + x*(Math.log(particles.springConstant+Math.exp(3.5))-3.5)/2;
			lambda = (Math.exp(y) - Math.exp(3.5))/particles.springConstant;
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
				if (lambda != 1 && lambda!=0){
					control.println("DryVolFrac = " + dryVolFrac);
					deltaF2Accumulator += -(particles.meanSpringEnergy()*(Math.log(particles.springConstant+Math.exp(3.5))-3.5)/2*(Math.exp(y))/particles.springConstant)*gaussLegendreWeights.get(weightIteration);
					control.println("deltaF2Accumulator = "+deltaF2Accumulator);
					//control.println("The value of x = "+x);
					
					pointIteration++;
					weightIteration++;
					
						if (pointIteration < gaussLegendrePoints.size() && weightIteration < gaussLegendreWeights.size()) {				
							x = gaussLegendrePoints.get(pointIteration); // set the value of x (for the change of variables)
							y = (Math.log(particles.springConstant+Math.exp(3.5))+3.5)/2 + x*(Math.log(particles.springConstant+Math.exp(3.5))-3.5)/2;
							lambda = (Math.exp(y) - Math.exp(3.5))/particles.springConstant;
							particles.lambda = lambda;
						}		

						else{
							/* reset the accumulators */
							deltaF2 = deltaF2Accumulator;
							deltaF2Accumulator = 0;
							control.println("deltaF2 = "+deltaF2);
							deltaF2PerVol = (deltaF2 * particles.N) / particles.totalVol; // ΔF2/vol
							control.println("deltaF2PerVol = "+deltaF2PerVol);
							deltaF2PerVolumeList.add(deltaF2PerVol);

							lambda = 1; // set the value of lambda to 1 for a fully interacting system
							particles.lambda = lambda;

							this.initialize();
							return;
						}				

						this.initialize();
						return;
					
				}
				
				if (lambda == 1) {
					control.println("lambda = "+lambda);
					/* Calculate the Flory–Rehner free energy */
					floryFEperVol = (particles.meanFreeEnergy() * particles.N) / particles.totalVol; // Flory-RehnerFR/vol
					floryFEperVolList.add(floryFEperVol); // update the list
					
					/* 	
						Calculate F0: the addittion of delta3(change between unconstrained solid and
						solid with fixed center of mass and the free energy of the refrence ideal
					 	Einstein crystal
					*/
					f0 = -3*(particles.N-1)/2.0*Math.log(Math.PI/particles.springConstant)-Math.log(particles.N)/2.0-Math.log(particles.totalVol);
					f0overN = f0/particles.N; // per particles
					f0PerVolume = (f0overN*particles.N)/particles.totalVol; // per volume
					f0PerVolumeList.add(f0PerVolume);

					meanPressures.add(particles.meanPressure()); // pressure from virial theorem
					density = particles.N/particles.totalVol;
					nnDistance = (1/Math.sqrt(2)) * Math.pow((4/density), 1/3.0); // the nearest neighbour distance
					
					// Calculate the Lindemann parameter <r>^2/rnn
					lindemannParameter = particles.meanSquareDisplacement()/(nnDistance * nnDistance);
					lindemannParameterList.add(lindemannParameter);

					if (nnDistance>(2*particles.meanRadius())){
						newHertzianPotential = 0;
					}
					else{
						newHertzianPotential = particles.B*Math.pow((1-nnDistance/(2*particles.meanRadius())), 2.5);
					}
					newHertzianPotentialList.add(newHertzianPotential);
					
					uPairPerVol = particles.meanPairEnergy()*(particles.N/particles.totalVol); // the pair energy per volume
					uPairPerVolList.add(uPairPerVol);

					reservoirVolFracList.add(particles.reservoirVolFrac); // the reservoir volume fraction list
					swellingRatioList.add(particles.meanRadius()); // the swelling ratio list

					lambda = 0; // reset lambda to 0
					particles.lambda = lambda;
					
					this.initialize();
					return;

				}

				if (lambda == 0) {
					uLatticePerN = particles.initialEnergy / particles.N; // the constant lattice energy
					/*
					 * F1: The free energy change between an interacting einstein crystal and a non-interacting einstein crystal
					 */
					logBoltzmannFactor = Math.log(particles.meanboltzmannFactor())/particles.N;
					deltaF1 = uLatticePerN - logBoltzmannFactor;
					deltaF1PerVol = (deltaF1 * particles.N) / particles.totalVol;
					deltaF1MinusUlattice = deltaF1-uLatticePerN;

					fPair = f0overN + deltaF1 + deltaF2; // the addition of all the components that make up the pair free energies of a molecular solid (the pair energies)
					fPairPerVol = (fPair * particles.N) / particles.totalVol; // pair free energy per volume
					
					/* the addittion of all the free energies contributing to the interaction of the particles */
					totalFreeEnergy = fPairPerVol + floryFEperVol;
					control.println("fPairPerVol: " + fPairPerVol);
					control.println("floryFEperVol: "+floryFEperVol);
					control.println("totalF: " + totalFreeEnergy);

					/* update the array lists */
					dryVolFracs.add(dryVolFrac);
					deltaF1List.add(deltaF1PerVol);
					totalVolList.add(particles.totalVol);
					pairFreeEnergyList.add(fPairPerVol);
					totalSumOfEnergiesList.add(totalFreeEnergy);
					logBoltzmannFactorList.add(logBoltzmannFactor);
					deltaF1MinusUlatticeList.add(deltaF1MinusUlattice);

					
					dryVolFrac += particles.dphi; // increment the dry volume fraction
					particles.dryVolFrac = dryVolFrac;
					
					// reset the iterations for the points and weights of the Gauss-Legendre polynomial
					pointIteration = 0;
					weightIteration = 0;
					x = gaussLegendrePoints.get(weightIteration); // reset the value of x
					y = (Math.log(particles.springConstant+Math.exp(3.5))+3.5)/2 + x*(Math.log(particles.springConstant+Math.exp(3.5))-3.5)/2;
					lambda = (Math.exp(y) - Math.exp(3.5))/particles.springConstant;	// reset lambda
					particles.lambda = lambda;
					

					this.initialize();
					return;
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
					control.println("pairPContribution = " + pairPContribution);

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
		control.setValue("Lambda increment", 0.1);
		control.setValue("DryVolFracStart", 0.016);
		// control.setValue("upper limit", upperLimit);
		control.setValue("DryVolFrac Max", 0.020);
		control.setValue("DryVolFrac increment", 0.0001);
		control.setValue("Initial configuration", "FCC");
		control.setValue("limit increment", xIncrement);
		control.setValue("Spring constant", 10000); // Spring constant: 2.035 for alpha/KT = 100
		// control.setValue("Initial configuration", "random-FCC");
		control.setValue("N", 108); // number of particles
		// control.setValue("N", 500); for FCC lattice, N/4 should be a perfect cube
		control.setValue("Dry radius [nm]", 50);
		control.setValue("x-link fraction", 0.001);
		// control.setValue("Dry volume fraction", 0.01);
		control.setValue("Young's calibration", 0.3); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
		// control.setValue("chi", 0.2);
		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0); // initially set to 0.05
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
			File outputFile = new File("data/HertzianSpheres"+particles.fileExtension+".txt");

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
			bw1.write("phi0,	mu/KT,		(PV/NkT)_total,		1+(PV/NkT)_virial,		1+(PV/NkT)_pair/calculated,		(PV/NkT)_FR,      <F_total>/V,      EntropySolid,		F0/V,       F1/V,     F2/V,     <F_FR>/V,     <alpha>,      zeta, 		deltaF1MinusUlattice, 		newHertzianPotential, 		lindemannParameter");
			bw1.newLine();

			for (int i = 1; i < dryVolFracs.size() - 1; i++) {
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				bw1.write(roundedDryVolFrac + ", " + chemicalPotentialList.get(i) + ", " + calculatedPressures.get(i) + ", " + + meanPressures.get(i) + ", " + pairPressuresListEdited.get(i) + ", " + floryRehnerPressuresListEdited.get(i) + ", " + totalSumOfEnergiesList.get(i) + ", " + (uPairPerVolList.get(i)-totalSumOfEnergiesList.get(i)+1.50) + ", " + f0PerVolumeList.get(i) + ", " + deltaF1List.get(i) + ", " + deltaF2PerVolumeList.get(i) + ", " +  floryFEperVolList.get(i) + ", " + swellingRatioList.get(i) + ", " + reservoirVolFracList.get(i) + ", " + deltaF1MinusUlatticeList.get(i) + ", " + newHertzianPotentialList.get(i) + ", " + lindemannParameterList.get(i));
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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresSolidPhaseGaussLegendreApp());
	}
}
