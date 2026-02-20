package org.opensourcephysics.sip.Hertz;

import java.awt.Color;
import java.text.DecimalFormat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

// imported for the purpose of creating an array to store the values of the dryVolFrac and totalSums
import java.util.List;
import java.util.ResourceBundle.Control;
import java.util.ArrayList;

/**
 * HertzSpheresApp performs a Monte Carlo simulation of microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @authors Alan Denton and Oreoluwa Alade
 * @version 1.2 25-05-2025
 * 
 */

public class HertzSpheresLiquidVirialCoefficientApp extends AbstractSimulation {
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
	List<Double> fExPerVolList = new ArrayList<>();
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
	List<Double> secondVirialCoefficientList = new ArrayList<>();
	List<Double> reducedB2List = new ArrayList<>();
	List<Double> hardSphereB2List = new ArrayList<>();
	List<Double> volumefractionList = new ArrayList<>();
	List<Double> virialCoefficientList = new ArrayList<>();


	// Virial Coefficient Lists
	List<Double> rhoList = new ArrayList<>();
	List<Double> zList = new ArrayList<>();
	List<Double> yList = new ArrayList<>();      // y = (Z-1)/rho
	List<Double> yMinusB2List = new ArrayList<>(); // optional: y - B2

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
			dryVolFrac = particles.dphi;     // start at phi0_min = dphi0
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


	private double[] fitVirialNoIntercept(int maxPower) {

		/**
		 * Fits higher-order virial coefficients (B3, B4, ..., B_{p+2})
		 * from simulation EOS data using least squares (no intercept).
		 *
		 * Fits:
		 *      (Z - 1)/ρ - B2 = B3 ρ + B4 ρ^2 + ... + B_{p+2} ρ^p
		 *
		 * Input:
		 *      maxPower = highest power of ρ in the fit.
		 *
		 * Returns:
		 *      [B3, B4, ..., B_{p+2}]
		 */

		int nPts = rhoList.size();          // number of density data points collected, e.g., 10 data points
		int p = maxPower;                   // number of fitting parameters,  e.g., 8 (fitting 8 parameters)

		double[] y = new double[nPts];      // target vector: y = (Z - 1)/ρ - B2
		double[][] x = new double[nPts][p]; // design matrix: columns = ρ^1 ... ρ^p (the independent variables)

		for (int i = 0; i < nPts; i++) {

			double rho = rhoList.get(i);     // density at state i, e.g., ρ = 0.0001
			y[i] = yMinusB2List.get(i);      // corresponding y-value for regression

			double rpow = rho;               // initialize with ρ^1

			for (int j = 0; j < p; j++) {
				x[i][j] = rpow;              // fill column j with ρ^(j+1)
				rpow *= rho;                 // update to next power: ρ^(j+2)
			}
		}

		OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression(); // create OLS solver
		ols.setNoIntercept(true);       // enforce zero intercept (no constant term in virial expansion)
		ols.newSampleData(y, x);        // provide regression data (y = Xβ)

		return ols.estimateRegressionParameters(); // returns fitted coefficients [B3, B4, ..., B_{p+2}]
	}


	
	public void doStep() {
		particles.step();

		// if initial configuration is random, no particle interactions for delay/10
		if (particles.steps > particles.delay / 10.0) {
			particles.scale = 1.0;
		}

		if (particles.steps <= particles.stop) {
			if (particles.steps > particles.delay) {
				if ((particles.steps - particles.delay) % particles.snapshotInterval == 0) {
					particles.sizeDistribution();
					if (structure) {
						rdf.update();
						ssf.update();
					}
					particles.calculateVolumeFraction();
				}
			}
		}

		if (particles.steps == particles.stop) {
			System.out.println("DryVolFrac = " + dryVolFrac + " dryVolFracMax= " + dryVolFracMax);

			// ===========================
			// Collect one data point per phi0
			// ===========================
			if (dryVolFrac < dryVolFracMax) {

				double rho   = particles.N / particles.totalVol;
				double alpha = particles.meanRadius();

				// measured EOS
				double Z = particles.meanPressure();

				// Compute B2 only ONCE (at first density point)
				double B2;
					
				// B2 depends on alpha (keep per point)
				B2 = particles.secondVirialCoefficient(alpha, alpha);
				
				// diagnostics: hard sphere normalisation
				double sigma   = 2.0 * alpha;
				double B2_HS   = particles.hardSphereB2(sigma);
				double B2_star = B2 / B2_HS;

				// virial fitting target: yMinusB2 = (Z-1)/rho - B2
				double y = (Z - 1.0) / rho;

				rhoList.add(rho);
				yMinusB2List.add(y - B2);
				secondVirialCoefficientList.add(B2);

				// keep these so writeData() can print them
				hardSphereB2List.add(B2_HS);
				reducedB2List.add(B2_star);

				// store measured quantities used in output
				meanPressures.add(Z);
				totalVolList.add(particles.totalVol);
				reservoirVolFracList.add(particles.reservoirVolFrac);
				volumefractionList.add(particles.meanVolFrac());
				swellingRatioList.add(alpha);
				dryVolFracs.add(dryVolFrac);

				// store FR per volume (state specific)
				floryFperVol = particles.meanFreeEnergy() * (particles.N / particles.totalVol);
				floryFperVolList.add(floryFperVol);

				// store Ideal term 
				stirlingApprox = (0.75 / Math.PI) * dryVolFrac
						* Math.log(2.0 * Math.PI * particles.N) / (2.0 * particles.N);

				idealFreeEnergy = (0.75 / Math.PI) * dryVolFrac * (Math.log(rho) - 1.0) + stirlingApprox;
				idealFreeEnergyList.add(idealFreeEnergy);

				// uPair/V (used this in QMelting)
				uPairPerVol = particles.meanPairEnergy() * (particles.N / particles.totalVol);
				uPairPerVolList.add(uPairPerVol);

				// increment phi0 and rerun
				dryVolFrac += particles.dphi;
				particles.dryVolFrac = dryVolFrac;

				this.initialize();
				return;
			}

			// ===========================
			// Scan complete: fit virial coefficients
			// ===========================
			int maxPower = 8; // fits B3..B7

			double[] virial = fitVirialNoIntercept(maxPower); // returns [B3, B4, ..., B_{maxPower+2}]

			virialCoefficientList.clear(); // clear in case of multiple runs
			for (int k = 0; k < virial.length; k++) {
				virialCoefficientList.add(virial[k]);
			} // Copies the fitted coefficients from the array into an ArrayList

			control.println("===== Virial fit from yMinusB2 =====");
			for (int k = 0; k < virial.length; k++) { // print B3, B4, ..., B_{maxPower+2}
				int Bindex = k + 3;
				control.println("B" + Bindex + " = " + virial[k]);
			}

			// ===========================
			// Build outputs per state i
			// ===========================
			
			int nStates = rhoList.size(); // e.g., 10 density points

			fExPerVolList.clear();                   // will store f_ex/V from virial
			varialPressuresList.clear();             // Z_virial
			floryRehnerPressuresListEdited.clear();  // Z_FR = Z_total - Z_virial
			calculatedPressures.clear();             // reconstructed Z_total = Z_virial + Z_FR
			totalSums.clear();                       // F_total/V
			chemicalPotList.clear();                 // mu/kT

			for (int i = 0; i < nStates; i++) { // initialize with zeros (will be overwritten for i=1..nStates-2)
				fExPerVolList.add(0.0);
				varialPressuresList.add(0.0);
				floryRehnerPressuresListEdited.add(0.0);
				calculatedPressures.add(0.0);
				totalSums.add(0.0);
				chemicalPotList.add(0.0);
			}

			for (int i = 1; i < nStates - 1; i++) {
				// // Skip first and last points for chemical potential calculation (central difference)

				double rho = rhoList.get(i); // density at state i
				double B2  = secondVirialCoefficientList.get(i); // B2 at state i 

				// Z_virial = 1 + B2*rho + B3*rho^2 + B4*rho^3 + ...
				double Zvir = 1.0 + B2 * rho;
				for (int k = 0; k < virial.length; k++) { // Loop over B3, B4, ..., B_{maxPower+2}
					int n = k + 3; // n = 3 for B3, n=4 for B4, etc.
					Zvir += virial[k] * Math.pow(rho, n - 1);
				}
				varialPressuresList.set(i, Zvir); // Replace the 0.0 at index i

				// Calculate F_ex/V (Excess Free Energy) from virial expansion:
				// f_ex/V = B₂ρ² + B₃ρ³/2 + B₄ρ⁴/3 + ...
				double fEx = B2 * rho * rho;
				for (int k = 0; k < virial.length; k++) {
					int n = k + 3;
					fEx += virial[k] * Math.pow(rho, n) / (n - 1.0);
				}
				fExPerVolList.set(i, fEx);

				// total free energy per volume
				double fTotal = floryFperVolList.get(i) + idealFreeEnergyList.get(i) + fEx; // F_total/V = F_FR/V + F_ideal/V + F_ex/V
				totalSums.set(i, fTotal);

				// FR contribution to Z (in Z-units)
				double Ztotal = meanPressures.get(i); // Measured from simulation
				double Zfr = Ztotal - Zvir; // Flory-Rehner contribution

				floryRehnerPressuresListEdited.set(i, Zfr);
				calculatedPressures.set(i, Zvir + Zfr); // should ~ Ztotal (fit error)
			}

			// ===========================
			// Chemical potential from central difference on F_total/V
			// ===========================
			for (int i = 1; i < nStates - 1; i++) {
				double Fnext = totalSums.get(i + 1);
				double Fprev = totalSums.get(i - 1);

				double mu = (4.0 * Math.PI / 3.0) * (Fnext - Fprev) / (2.0 * particles.dphi);
				chemicalPotList.set(i, mu);
			}

			writeData();
		}

		// plot mean energy, pressure, swelling ratio
		energyData.append(0, particles.steps, particles.meanPairEnergy());
		pressureData.append(1, particles.steps, particles.meanPressure());
		sizeData.append(2, particles.steps, particles.meanRadius());

		display3d.setMessage("Number of steps: " + particles.steps);

		if (control.getBoolean("Visualization on")) {
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
		control.setValue("DryVolFrac increment", 0.0001);
		control.setValue("DryVolFrac Max", 0.0022);
		control.setValue("Initial configuration", "FCC");
		//control.setValue("Initial configuration", "random-FCC");
		control.setValue("N", 108); // number of particles
		//control.setValue("N", 500); for FCC lattice, N/4 should be a perfect cube
        control.setValue("Dry radius [nm]", 50);
        control.setValue("x-link fraction", 0.00005); // 0.001
        // control.setValue("Dry volume fraction", 0.01);
        control.setValue("Young's calibration", 1.0); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
        control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.05);
		control.setValue("Delay", 10000); // steps after which statistics collection starts
		control.setValue("Snapshot interval", 100); // steps separating successive samples 
		control.setValue("Stop", 100000); // steps after which statistics collection stops
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
			File outputFile = new File("data/APS_2026/Liquid_Phase/HertzSpheresLiquidVirialCoefficient5e-5" + particles.fileExtension + ".txt");
		
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
			// bw1.write("Mean pair energy <E_pair>/N [kT]: " + particles.meanPairEnergy());
			// bw1.newLine();
			bw1.write("The coupling constant increment - dlambda: " + particles.dlambda);
			bw1.newLine();	
			bw1.write("phi0, phi, rho, mu/kT, (PV/NKT)total, (PV/NKT)_Virial, (PV/NKT)_reconstructed, (PV/NKT)_FR, <F_total>/V, QMelting, <F_id>/V, <F_ex>/V, <F_FR>/V, <alpha>, Zeta, B2, B2_HS, B2*");
			bw1.newLine();

			for (int i = 1; i < dryVolFracs.size() - 1; i++) {

				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));

				double Z_total = meanPressures.get(i);                        // measured PV/NkT
				double rho     = rhoList.get(i);                                  // density at state i
				double Z_vir   = varialPressuresList.get(i);                 // virial from fit
				double Z_fr    = floryRehnerPressuresListEdited.get(i);       // Z_total - Z_vir
				double Z_rec   = calculatedPressures.get(i);                 // Z_vir + Z_fr (diagnostic)

				double FtotV = totalSums.get(i);
				double FidV  = idealFreeEnergyList.get(i);
				double FexV  = fExPerVolList.get(i);                       // excess free energy per volume
				double FfrV  = floryFperVolList.get(i);

				double qMelting = (uPairPerVolList.get(i) - FtotV + 1.50);

				bw1.write(
					roundedDryVolFrac + ", " +
					volumefractionList.get(i) + ", " +
					rho + ", " +
					chemicalPotList.get(i) + ", " +
					Z_total + ", " +
					Z_vir + ", " +
					Z_rec + ", " +
					Z_fr + ", " +
					FtotV + ", " +
					qMelting + ", " +
					FidV + ", " +
					FexV + ", " +
					FfrV + ", " +
					swellingRatioList.get(i) + ", " +
					reservoirVolFracList.get(i) + ", " +
					secondVirialCoefficientList.get(i) + ", " +
					hardSphereB2List.get(i) + ", " +
					reducedB2List.get(i)
				);
				bw1.newLine();
			}


			bw1.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		// Write virial coefficients to a separate file
		try {
			File virialFile = new File("data/APS_2026/Liquid_Phase/VirialCoefficientsXlink4e-5" + particles.fileExtension + ".txt");
			
			// Create parent directory if it doesn't exist
			File virialDir = virialFile.getParentFile();
			if (virialDir != null && !virialDir.exists()) {
				virialDir.mkdirs();
			}
			
			if (!virialFile.exists()) {
				virialFile.createNewFile();
			}

			FileWriter fw2 = new FileWriter(virialFile.getAbsoluteFile());
			BufferedWriter bw2 = new BufferedWriter(fw2);

			// Write header with system parameters
			bw2.write("===== Virial Coefficients for Hertzian Microgels =====");
			bw2.newLine();
			bw2.write("x-link fraction: " + particles.xLinkFrac);
			bw2.newLine();
			bw2.write("Young's calibration: " + particles.Young);
			bw2.newLine();
			bw2.write("chi: " + particles.chi);
			bw2.newLine();
			bw2.write("Dry radius [nm]: " + particles.dryR);
			bw2.newLine();
			bw2.write("Number of density points: " + rhoList.size());
			bw2.newLine();
			bw2.write("Density range: rho_min = " + rhoList.get(0) + ", rho_max = " + rhoList.get(rhoList.size()-1));
			bw2.newLine();
			bw2.newLine();

			// Write B2 values at each density point
			bw2.write("===== B2 at Each Density Point =====");
			bw2.newLine();
			bw2.write("rho, B2, B2_HS, B2*");
			bw2.newLine();
			for (int i = 0; i < rhoList.size(); i++) {
				bw2.write(rhoList.get(i) + ", " + 
						secondVirialCoefficientList.get(i) + ", " + 
						hardSphereB2List.get(i) + ", " + 
						reducedB2List.get(i));
				bw2.newLine();
			}
			bw2.newLine();

			// Write higher-order virial coefficients (B3, B4, ..., B10)
			bw2.write("===== Higher-Order Virial Coefficients (from OLS fit) =====");
			bw2.newLine();
			bw2.write("Coefficient, Value");
			bw2.newLine();
			for (int k = 0; k < virialCoefficientList.size(); k++) {
				int Bindex = k + 3;  // B3, B4, B5, ..., B10
				bw2.write("B" + Bindex + ", " + virialCoefficientList.get(k));
				bw2.newLine();
			}

			bw2.close();
			System.out.println("Virial coefficients written to: " + virialFile.getAbsolutePath());
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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresLiquidVirialCoefficientApp());
	}
}
