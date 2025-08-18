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

// imported for the purpose of creating an array to store the values of the dryVolFrac and totalSums
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

public class HertzSpheresLiquidDensityApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
    HertzSpheresLiquidDensity particles = new HertzSpheresLiquidDensity();
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
	double lambda, dlambda;
	double freeEnergyTotal, pairEnergyTotal;
	double dryVolFrac, dryVolFracStart, dryVolFracMax;
	double totalSum;
	boolean incrementDryVolFrac = true;
	boolean flag = false;
	boolean resetLambda = true;
	double totalVol;
	double floryFperVol, idealFreeEnergy;
	boolean isFirstIteration;
	double pairEnergy = 0;
	double n, pairPressure, calculatedPressure;
	double initialfloryFperVol, initialpairEperVol;
	double stirlingApprox, pairEperVol, totalFreeEnergies;
	boolean firstResidualPart;
	double floryPressure, deltaFpairDividedByDeltaPhi, dFRFE_dPhi, dFPair_dPhi;
	double freeEnergyPressure, floryPContribution, pairPressureContribution;
	double pairPContribution;
	double initialFreeEnergy, latestFreeEnergy;
	double latestFREperVol, latestPairEperVol;
	double currentDryVolFrac, previousDryVolFrac, nextDryVolFrac;
	int currentIndex, previousIndex, nextIndex;
	double currentPairEnergy, previousPairEnergy, nextPairEnergy;
	double currentFRFreeEnergy, previousFRFreeEnergy, nextFRFreeEnergy;
	double chemicalPot;
	double latestFE, initialFE, chemicalPotential;
	double uPairPerVol;
	double density, nnDistance, newHertzianPotential;


	List<Double> dryVolFracs = new ArrayList<>();
    List<Double> totalSums = new ArrayList<>();
	List<Double> chemicalPotList = new ArrayList<>();
	List<Double> freeEnergyTotals = new ArrayList<>();
	List<Double> floryFperVolList = new ArrayList<>();
	List<Double> idealFreeEnergyList = new ArrayList<>();
	List<Double> pairEperVolList = new ArrayList<>();
	List<Double> calculatedPressures = new ArrayList<>();
	List<Double> meanPressures = new ArrayList<>();
	List<Double> uPairPerVolList = new ArrayList<>();
	List<Double> pairPressuresList = new ArrayList<>();
	List<Double> varialPressuresList = new ArrayList<>();
	List<Double> floryRehnerPressuresList = new ArrayList<>();
	List<Double> totalPressuresList = new ArrayList<>();
	List<Double> reservoirVolFracList = new ArrayList<>();
	List<Double> swellingRatioList = new ArrayList<>();
	List<Double> floryRehnerPressuresListEdited = new ArrayList<>();
	List<Double> pairPressuresListEdited = new ArrayList<>();
	List<Double> dryVolFracsEdited = new ArrayList<>();
	List<Double> totalVolList = new ArrayList<>();
	List<Double> newHertzianPotentialList = new ArrayList<>();


	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values

    /**
	* Initializes the model.
	*/
	public void initialize() {
		
		added = false;
		dryVolFracMax = control.getDouble("DryVolFrac Max");
		particles.dphi = control.getDouble("DryVolFrac increment");
		particles.N = control.getInt("N"); // number of particles
		String configuration = control.getString("Initial configuration");
		particles.initConfig = configuration;
        particles.dryR = control.getDouble("Dry radius [nm]");
        particles.xLinkFrac = control.getDouble("x-link fraction");
        //particles.dryVolFrac = control.getDouble("Dry volume fraction");
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
		  set the value of dryVolFrac to the initial starting value of the 
		  dry Volume Fraction after which it becomes false every other times
		*/
		if (incrementDryVolFrac){		
			n = 1;							
			dryVolFrac = particles.dphi/2.0;
		//	dryVolFrac = particles.dphi;
			particles.dryVolFrac = dryVolFrac;
			incrementDryVolFrac = false;
		}
		
		particles.initialize(configuration);

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
			System.out.println("DryVolFrac = " + dryVolFrac + " dryVolFracMax= " + dryVolFracMax);
			if (dryVolFrac < dryVolFracMax){
					// increment the pairEnergy 
					density = (particles.N/particles.totalVol);
					pairEnergy += (particles.meanPressure()-1)*2/n; // n is a double that increases by a factor of 2

					control.println("Number of MC steps = "+particles.steps);
					control.println("DryVolFrac = " + dryVolFrac);
					control.println("pairEnergy = "+pairEnergy);

					// calculate the stirling approximation
					stirlingApprox = (0.75/Math.PI)*dryVolFrac*Math.log(2*Math.PI*particles.N)/(2*particles.N); // stirling approximation to the free energy

					// calculate the ideal free energy contribution
					idealFreeEnergy = (0.75/Math.PI)*dryVolFrac*(Math.log(density)-1.0)+stirlingApprox; // ideal gas free energy contribution
					idealFreeEnergyList.add(idealFreeEnergy); // update the list

					uPairPerVol = particles.meanPairEnergy()*(particles.N/particles.totalVol);
					uPairPerVolList.add(uPairPerVol);

					pairEperVol = (pairEnergy*particles.N)/particles.totalVol; // the total pair energy per volume
					control.println("fPair/Vol = "+pairEperVol);

					floryFperVol = (particles.meanFreeEnergy() * particles.N)/particles.totalVol; // flory-rehner free energy per volume
					control.println("fFlory/Vol = "+floryFperVol);

					nnDistance = (1/Math.sqrt(2)) * Math.pow((4/density), 1/3.0);
					newHertzianPotential = (particles.B * Math.pow((1 - (nnDistance/(2 * particles.meanRadius()))), 2.5));
					control.println("nnDistance = "+nnDistance);
					control.println("(1 - (nnDistance/(2 * particles.meanRadius()))) = " + (1 - (nnDistance/(2 * particles.meanRadius()))));
					newHertzianPotentialList.add(newHertzianPotential);

					/* update the lists accordingly */
					pairEperVolList.add(pairEperVol);
					floryFperVolList.add(floryFperVol);
					
					totalFreeEnergies = floryFperVol+idealFreeEnergy+pairEperVol; // total free energy contribution
					control.println("totalF/V = "+totalFreeEnergies);
	
					/* Update the list */
					totalSums.add(totalFreeEnergies);
					dryVolFracs.add(dryVolFrac); 
					totalVolList.add(particles.totalVol);
					reservoirVolFracList.add(particles.reservoirVolFrac); // the reservoir volume fraction
					swellingRatioList.add(particles.meanRadius()); // the swelling ratio
					meanPressures.add(particles.meanPressure()); // quantity PV/NkT //the 1 is coming from the ideal gas
					dryVolFrac += particles.dphi; // increment the dry volume fraction

					particles.dryVolFrac = dryVolFrac;
					n+=2; // the denominator that arises from changing the variables for the pair energy

					this.initialize();
					return;

			}
			//if (dryVolFracs.size() > 1) { // checks whether the size list is greater than 1
				// NaN is used to indicate that a numeric value is not available
				double currentPairFreeEnergy = Double.NaN;
				double previousPairFreeEnergy = Double.NaN;
				double nextPairFreeEnergy = Double.NaN;
				double currentFloryFreeEnergy = Double.NaN;
				double previousFloryFreeEnergy = Double.NaN;
				double nextFloryFreeEnergy = Double.NaN;
				int i;

				/* Initialize the lists */

				floryRehnerPressuresListEdited.add(0, 0.0);
				pairPressuresListEdited.add(0, 0.0);
				chemicalPotList.add(0,0.0);
				calculatedPressures.add(0,0.0);

				for (i = 1; i < dryVolFracs.size() - 1; i++){
					double currentDryVolFrac = dryVolFracs.get(i);
					double previousDryVolFrac = dryVolFracs.get(i - 1);
					double nextDryVolFrac = (i + 1 < dryVolFracs.size()) ? dryVolFracs.get(i + 1) : Double.NaN;

					dryVolFracsEdited.add(currentDryVolFrac);
					double value = dryVolFracsEdited.get(i - 1);
					control.println("value = " + value);

					double currentfloryFE = floryFperVolList.get(i);
					control.println("currentFloryFEIndex = " + currentfloryFE);

					double currentPairFE = pairEperVolList.get(i);
					control.println("currentPairFEIndex = " + currentPairFE);

					double currentVolume = totalVolList.get(i); // the volume of the system for that dry volume fraction
					control.println("currentVolumeIndex = " + currentVolume);

					currentPairFreeEnergy = pairEperVolList.get(i); // pair free energy for that dry volume fraction
					currentFloryFreeEnergy = floryFperVolList.get(i); // FR free energy for that dry volume fraction
					control.println("Current Pair Free Energy = " + currentPairFreeEnergy);
					control.println("Current Flory Free Energy = " + currentFloryFreeEnergy);

					initialFE = totalSums.get(i-1);
					previousPairFreeEnergy = pairEperVolList.get(i - 1);
					previousFloryFreeEnergy = floryFperVolList.get(i - 1);
					control.println("Previous Pair Free Energy = " + previousPairFreeEnergy);
					control.println("Previous Flory Free Energy = " + previousFloryFreeEnergy);

					if (i + 1 < dryVolFracs.size()) {
						nextFloryFreeEnergy = floryFperVolList.get(i + 1);
						nextPairFreeEnergy = pairEperVolList.get(i + 1);
						latestFE = totalSums.get(i + 1);
					}
					
					control.println("Next Pair Free Energy = " + nextPairFreeEnergy);
					control.println("Next Flory Free Energy = " + nextFloryFreeEnergy);

					chemicalPotential = (4.0*Math.PI/3.0)*(latestFE - initialFE)/(2.0 * particles.dphi); // calculate the chemical potential
					dFRFE_dPhi = (nextFloryFreeEnergy - previousFloryFreeEnergy) / (2.0 * particles.dphi); // dF_FR/dPhi0
					dFPair_dPhi = (nextPairFreeEnergy - previousPairFreeEnergy) / (2.0 * particles.dphi); // dF_pair/dPhi0

					control.println("dFRFE_dPhi = " + dFRFE_dPhi);
					control.println("dFPair_dPhi = " + dFPair_dPhi);

					floryPressure = ((dFRFE_dPhi) * currentDryVolFrac - currentfloryFE); //  phi0*dF_FR/dPhi0 - F_FR
					pairPressure = ((dFPair_dPhi) * currentDryVolFrac - currentPairFE); // phi0*dF_pair/dPhi0 - F_pair
					
					floryPContribution = (floryPressure * currentVolume) / particles.N; // P_FR*V/N(KT)
					pairPContribution = (pairPressure * currentVolume) / particles.N; // P_pair*V/N(KT)

					floryRehnerPressuresListEdited.add(floryPContribution); 
					pairPressuresListEdited.add(pairPContribution + 1.0);
					calculatedPressures.add(floryPContribution + pairPContribution + 1.0);
					chemicalPotList.add(chemicalPotential);

				}
			//}
			
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
		control.setValue("DryVolFrac increment", 0.0001);
		control.setValue("DryVolFrac Max", 0.02);
		control.setValue("Initial configuration", "FCC");
		//control.setValue("Initial configuration", "random-FCC");
		control.setValue("N", 32); // number of particles
		//control.setValue("N", 500); for FCC lattice, N/4 should be a perfect cube
        control.setValue("Dry radius [nm]", 50);
        control.setValue("x-link fraction", 0.001);
        // control.setValue("Dry volume fraction", 0.01);
        control.setValue("Young's calibration", 1.0); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
		//control.setValue("chi", 0.2);
        control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0);
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
    	control.println("Number of MC steps = "+particles.steps);
		control.println("<Epair>/N/lambda= " + decimalFormat.format(particles.meanPairEnergy()/particles.lambda));
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
			File outputFile = new File("data/HertzSpheresLiquid" + particles.fileExtension + ".txt");
		
			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}

			FileWriter fw1 = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);

			// write system parameters to the file in addittion to dryVolFracs and totalFree energies
			bw1.write("Number of particles: " + particles.N);
			bw1.newLine();
			bw1.write("dryVolFracMax: " + dryVolFracMax);
			bw1.newLine();
			bw1.write("dryVolFracIncrement: " + particles.dphi);
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
			bw1.write("DryVolFrac, mu/kT, (PV/NKT)total, (PV/NKT)_Virial, (PV/NKT)_pair/calculated, (PV/NKT)_FR, <F_total>/V, QMelting, <F_id>/V, <F_pair>/V, <F_FR>/V,  <alpha>, Zeta");
			bw1.newLine();		

			for (int i = 1; i < dryVolFracs.size() - 1; i++) {
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				bw1.write(roundedDryVolFrac + ", " + chemicalPotList.get(i) + ", " + (meanPressures.get(i)+floryRehnerPressuresListEdited.get(i)) + ", " + meanPressures.get(i) + ", " + pairPressuresListEdited.get(i) + ", " +  floryRehnerPressuresListEdited.get(i) + ", " + totalSums.get(i) + ", "  + (uPairPerVolList.get(i)-totalSums.get(i) + 1.50) + ", " + idealFreeEnergyList.get(i)  + ", " + pairEperVolList.get(i)+ ", " + floryFperVolList.get(i) + ", " + swellingRatioList.get(i) + ", " + reservoirVolFracList.get(i));
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
	 * command line parameters
	 */
	public static void main(String[] args) { // set up animation control
			@SuppressWarnings("unused")
			SimulationControl control = SimulationControl.createApp(new HertzSpheresLiquidDensityApp());
	}
}
