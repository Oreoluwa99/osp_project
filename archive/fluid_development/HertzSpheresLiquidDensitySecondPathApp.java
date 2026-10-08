package org.opensourcephysics.sip.Hertz;

import java.awt.Color;
import java.text.DecimalFormat;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.PlotFrame;

import java.util.List;
import java.util.ResourceBundle.Control;
import java.util.ArrayList;

/**
 * HertzSpheresApp performs a Monte Carlo simulation of microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * This is specifically for the fluid phase and incorporates the second path to getting the free energy 
 * for the fluid phase via integration over x-link fraction at fixed density and adding that to the 
 * ideal and Flory-Rehner free energy contributions
 * 
 * @authors Alan Denton and Oreoluwa Alade
 * @version 1.2 25-05-2025
 * 
 */

public class HertzSpheresLiquidDensitySecondPathApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
    HertzSpheresLiquidDensitySecondPath particles = new HertzSpheresLiquidDensitySecondPath();
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
	double freeEnergyTotal, pairEnergyTotal;
	double dryVolFrac, dryVolFracStart, dryVolFracMax;
	double totalSum;
	boolean incrementDryVolFrac = true;
	boolean flag = false;
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
	double uPairPerVol=0;
	double density, nnDistance, newHertzianPotential;
	double xLinkFracStart, xLinkFracMax;
	double dxLink;
	boolean incrementXLinkFrac = true;
	double xLinkFrac;
	double xLinkFracMin;

	List<Double> dryVolFracs = new ArrayList<>();
	List<Double> total_free_energy_path2 = new ArrayList<>();
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
	List<Double> secondVirialCoefficientList = new ArrayList<>();
	List<Double> reducedB2List = new ArrayList<>();
	List<Double> hardSphereB2List = new ArrayList<>();
	List<Double> volumefractionList = new ArrayList<>();

	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values

    /**
	* Initializes the model.
	*/
	public void initialize() {
		
		    added = false;
			dryVolFracMax = control.getDouble("DryVolFrac Max");
			particles.dphi = control.getDouble("DryVolFrac increment");
			particles.N = control.getInt("N");
			String configuration = control.getString("Initial configuration");
			particles.initConfig = configuration;
			particles.dryR = control.getDouble("Dry radius [nm]");
			
			// Get x-link parameters
			xLinkFracMin = control.getDouble("x-link fraction min");  // ADD THIS
			xLinkFracMax = control.getDouble("x-link fraction max");
			dxLink = control.getDouble("x-Link increment");
			
			particles.Young = control.getDouble("Young's calibration");
			particles.chi = control.getDouble("chi");
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
		  set the value of dryVolFrac to the initial starting value of the 
		  dry Volume Fraction after which it becomes false every other times
		*/
		if (incrementDryVolFrac) {
			n = 1;
			dryVolFrac = particles.dphi / 2.0;
			particles.dryVolFrac = dryVolFrac;
			incrementDryVolFrac = false;
		}
    
		xLinkFrac = xLinkFracMax;  // Use Th instead of T*
    	particles.xLinkFrac = xLinkFrac;
    
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

				if (dryVolFrac < dryVolFracMax) {
					
					// ================================================
					// PATH 1: Density integration at T_h (xLinkFracMax)
					// ================================================
					
					density = (particles.N / particles.totalVol);
					pairEnergy += (particles.meanPressure() - 1) * 2 / n;
					
					control.println("\n=== PATH 1: DENSITY POINT AT Th ===");
					control.println("DryVolFrac = " + dryVolFrac);
					
					// Calculate free energy components at Th
					stirlingApprox = (0.75/Math.PI) * dryVolFrac * Math.log(2*Math.PI*particles.N) / (2*particles.N);
					idealFreeEnergy = (0.75/Math.PI) * dryVolFrac * (Math.log(density) - 1.0) + stirlingApprox;
					floryFperVol = (particles.meanFreeEnergy() * particles.N) / particles.totalVol;
					
					// Path 1 free energy (at Th)
					pairEperVol = (pairEnergy * particles.N) / particles.totalVol;
					double F_path1_total = idealFreeEnergy + floryFperVol + pairEperVol;

					// B2 calculation
					double currentAlpha = particles.meanRadius();
					double B2 = particles.secondVirialCoefficient(currentAlpha, currentAlpha);
					double sigma = 2.0 * currentAlpha;
					double B2_HS = particles.hardSphereB2(sigma);
					double B2_star = B2 / B2_HS;
					
					// Store Path 1 data
					idealFreeEnergyList.add(idealFreeEnergy);
					floryFperVolList.add(floryFperVol);
					secondVirialCoefficientList.add(B2);
					hardSphereB2List.add(B2_HS);
					reducedB2List.add(B2_star);
					swellingRatioList.add(currentAlpha);
					dryVolFracs.add(dryVolFrac);
					totalVolList.add(particles.totalVol);
					reservoirVolFracList.add(particles.reservoirVolFrac);
					volumefractionList.add(particles.meanVolFrac());
					meanPressures.add(particles.meanPressure());
					pairEperVolList.add(pairEperVol);
					
					// ================================================
					// PATH 2: Temperature integration from Th to T*
					// ================================================
					
					control.println("\n=== PATH 2: TEMPERATURE INTEGRATION ===");
					control.println("Integrating from Th (x=" + xLinkFracMax + ") DOWN to T* (x=" + xLinkFracMin + ")");
					control.println("At FIXED density: phi0 = " + dryVolFrac);
					
					// Store current x-link value
					double originalXlink = particles.xLinkFrac;
					
					double F_temperature_integral = 0.0;
					
					List<Double> xValues = new ArrayList<>();
					List<Double> uPairValues = new ArrayList<>();
					
					// START from Th and go DOWN to T*
					double x_current = xLinkFracMax;
					int xlinkCount = 0;
					
					while (x_current >= xLinkFracMin) {
						
						control.println("\n  x = " + x_current);
						
						particles.xLinkFrac = x_current;
						particles.dryVolFrac = dryVolFrac;  // Fixed density
						
						// Collect pair energy
						double uPair = particles.meanPairEnergy();

						double integrand = uPair / (x_current * particles.N); 

						xValues.add(x_current);
						uPairValues.add(uPair);
						
						// Accumulate temperature integral
						F_temperature_integral += integrand * dxLink;
						
						x_current -= dxLink;  // Go DOWN
						xlinkCount++;
					}
					
					// -----------------------------------------------
					// Combine Path 1 + Temperature integral = Path 2
					// -----------------------------------------------
					
					double F_temp_integral_perVol = (F_temperature_integral * particles.N) / particles.totalVol;
					
					// Path 2 formula: F2 = [F_id(Th) + F_FR(Th) + Virial(Th)] + Temp_integral
					//                    =  Path 1 result              + Temp_integral
					double F_path2_total = F_path1_total + F_temp_integral_perVol;
					total_free_energy_path2.add(F_path2_total);
					
					// ================================================
					// OUTPUT
					// ================================================
					
					control.println("\n=== PATH 2 COMPLETE ===");
					control.println("Path 1 at Th:          " + F_path1_total);
					control.println("Temp integral:         " + F_temp_integral_perVol);
					control.println("Path 2 total:          " + F_path2_total);
					
					// Store temperature integral
					uPairPerVol = F_temp_integral_perVol;
					uPairPerVolList.add(uPairPerVol);
					
					// Restore original x-link
					particles.xLinkFrac = originalXlink;
					
					// ================================================
					// Increment to next density
					// ================================================
					
					dryVolFrac += particles.dphi;
					particles.dryVolFrac = dryVolFrac;
					n += 2;
					
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

					initialFE = total_free_energy_path2.get(i-1);
					previousPairFreeEnergy = pairEperVolList.get(i - 1);
					previousFloryFreeEnergy = floryFperVolList.get(i - 1);
					control.println("Previous Pair Free Energy = " + previousPairFreeEnergy);
					control.println("Previous Flory Free Energy = " + previousFloryFreeEnergy);

					if (i + 1 < dryVolFracs.size()) {
						nextFloryFreeEnergy = floryFperVolList.get(i + 1);
						nextPairFreeEnergy = pairEperVolList.get(i + 1);
						latestFE = total_free_energy_path2.get(i + 1);
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
		control.setValue("DryVolFrac increment", 0.001);
		control.setValue("DryVolFrac Max", 0.02);
		control.setValue("Initial configuration", "FCC");
		control.setValue("N", 108);
		control.setValue("Dry radius [nm]", 50);
		
		// X-link fraction settings for Path 2
		control.setValue("x-link fraction min", 0.001);  // 0.00003
		control.setValue("x-link fraction max", 0.01);  // 0.0001
		control.setValue("x-Link increment", 0.001);  // 0.000003
		
		control.setValue("Young's calibration", 1.0);
		control.setValue("chi", 0);
		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.05);
		control.setValue("Delay", 10000);
		control.setValue("Snapshot interval", 100);
		control.setValue("Stop", 20000);
		control.setValue("Size bin width", .001);
		control.setValue("g(r) bin width", .005);
		control.setValue("Delta k", .005);
		control.setValue("File extension", "1");
		control.setValue("Calculate structure", false);
		control.setAdjustableValue("Visualization on", true);
	}

	public void stop() {
		// particles.lambda = lambda;
    	control.println("Number of MC steps = "+particles.steps);
		// control.println("<Epair>/N/lambda= " + decimalFormat.format(particles.meanPairEnergy()/particles.lambda));
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
			File outputFile = new File("data/APS_2026/Liquid_Phase/HertzSpheresLiquidSecondPath" + particles.fileExtension + ".txt");
		
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
			bw1.write("#Minimum x-link fraction: " + xLinkFracMin);
			bw1.newLine();
			bw1.write("#Maximum x-link fraction: " + xLinkFracMax);
			bw1.newLine();
			bw1.write("#x-link fraction increment: " + dxLink);
			bw1.newLine();
			bw1.write("Reservoir swelling ratio: " + particles.reservoirSR);
			bw1.newLine();
			// bw1.write("Reservoir volume fraction: " + particles.reservoirVolFrac);
			// bw1.newLine();
			// bw1.write("Actual volume fraction: " + particles.volFrac);
			// bw1.newLine();
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
			// bw1.write("Mean pair energy <E_pair>/N [kT]: " + particles.meanPairEnergy());
			// bw1.newLine();
			// bw1.write("The coupling constant increment - dlambda: " + particles.dlambda);
			// bw1.newLine();	
			bw1.write("phi0, phi, mu/kT, (PV/NKT)total, (PV/NKT)_Virial, (PV/NKT)_pair/calculated, (PV/NKT)_FR, <F_total_path2>/V, QMelting, <F_id>/V, <F_pair>/V, <F_FR>/V,  <alpha>, Zeta, B2, B2_HS, B2*");			bw1.newLine();		

			for (int i = 1; i < dryVolFracs.size() - 1; i++) {
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				bw1.write(roundedDryVolFrac + ", " + volumefractionList.get(i) + ", " + chemicalPotList.get(i) + ", " + (meanPressures.get(i)+floryRehnerPressuresListEdited.get(i)) + ", " + meanPressures.get(i) + ", " + pairPressuresListEdited.get(i) + ", " +  floryRehnerPressuresListEdited.get(i) + ", " + total_free_energy_path2.get(i) + ", "  + (uPairPerVolList.get(i)-total_free_energy_path2.get(i) + 1.50) + ", " + idealFreeEnergyList.get(i)  + ", " + pairEperVolList.get(i)+ ", " + floryFperVolList.get(i) + ", " + swellingRatioList.get(i) + ", " + reservoirVolFracList.get(i) + ", " + secondVirialCoefficientList.get(i) + ", " + hardSphereB2List.get(i) + ", " + reducedB2List.get(i));
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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresLiquidDensitySecondPathApp());
	}
}
