package org.opensourcephysics.sip.Hertz;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.PlotFrame;

/**
 * Monte Carlo simulation application for the fluid phase of
 * compressible microgels described by the nonlocal facet model.
 *
 * <p>The application scans dry volume fraction, collects the
 * equation of state, fits higher-order virial coefficients,
 * reconstructs the liquid free-energy density, and estimates
 * the chemical potential by finite differences.</p>
 *
 * @author Oreoluwa Alade
 * @author Alan R. Denton
 */

public class HertzSpheresNonLocalFacetFluidFixedTestApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
    //HertzSpheresLiquidDensity particles = new HertzSpheresLiquidDensity();
	// HertzSpheresNonLocalFacet particles = new HertzSpheresNonLocalFacet(); -- the original code I used for the fluid phase
	HertzSpheresNonLocalFacet_FluidFixed particles = new HertzSpheresNonLocalFacet_FluidFixed();
    PlotFrame energyData = new PlotFrame("MC steps", "<E_pair>/N", "Mean pair energy per particle");
    PlotFrame pressureData = new PlotFrame("MC steps", "PV/NkT", "Mean pressure");
    PlotFrame sizeData = new PlotFrame("MC steps", "alpha", "Mean swelling ratio");
    Display3DFrame display3d = new Display3DFrame("Simulation animation");
    ElementSphere nanoSphere[];
    boolean added = false;
    boolean structure;
    RDF rdf;
    SSF ssf;
	double lambda=1;
	boolean incrementDryVolFrac = true;
	double dryVolFrac;
	double dryVolFracMax;
	double floryFperVol;
	double stirlingApprox;
	double uPairPerVol;
	double fixedB2;
	int maxPower;

	List<Double> dryVolFracs = new ArrayList<>();
    List<Double> totalSums = new ArrayList<>();
	List<Double> chemicalPotList = new ArrayList<>();
	List<Double> floryFperVolList = new ArrayList<>();
	List<Double> idealFreeEnergyList = new ArrayList<>();
	List<Double> fExPerVolList = new ArrayList<>();
	List<Double> meanPressures = new ArrayList<>();
	List<Double> uPairPerVolList = new ArrayList<>();
	List<Double> virialPressuresList = new ArrayList<>();
	List<Double> reservoirVolFracList = new ArrayList<>();
	List<Double> swellingRatioList = new ArrayList<>();
	List<Double> totalVolList = new ArrayList<>();
	List<Double> secondVirialCoefficientList = new ArrayList<>();
	List<Double> reducedB2List = new ArrayList<>();
	List<Double> hardSphereB2List = new ArrayList<>();
	List<Double> volumefractionList = new ArrayList<>();
	List<Double> virialCoefficientList = new ArrayList<>();
	List<Double> residualPressuresList = new ArrayList<>();


	// Virial Coefficient Lists
	List<Double> rhoList = new ArrayList<>();
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
		particles.fileExtension = control.getString("Run number");
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

		// calculate the second virial coefficient so we don't have to recalculate it for every state point
		if (dryVolFracs.isEmpty()) {
			double alpha0 = particles.reservoirSR;
			fixedB2 = particles.secondVirialCoefficient(alpha0, alpha0);

			System.out.println("Fixed B2 = " + fixedB2);
			System.out.println("Reference alpha0 = " + alpha0);
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
			System.out.println("i=" + i + ", rho=" + rho + ", yMinusB2=" + yMinusB2List.get(i));
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

				System.out.println("phi0 = " + dryVolFrac
					+ " | Z_measured = " + Z
					+ " | rho = " + rho
					+ " | alpha = " + alpha);

				// Use the fixed B2 calculated from the reservoir swelling ratio
				double B2 = fixedB2;
									
				control.println("B2 = " + B2);
				
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

				// System.out.println("phi0 = " + dryVolFrac);
				// System.out.println("meanFreeEnergy per particle = " + particles.meanFreeEnergy());
				// System.out.println("N/totalVol = " + (particles.N / particles.totalVol));
				// System.out.println("F_FR per volume = " + floryFperVol);

				floryFperVolList.add(floryFperVol);

				// store Ideal term 
				stirlingApprox = (0.75 / Math.PI) * dryVolFrac
						* Math.log(2.0 * Math.PI * particles.N) / (2.0 * particles.N);

				double idealFreeEnergy = (0.75 / Math.PI) * dryVolFrac * (Math.log(rho) - 1.0) + stirlingApprox;
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
			maxPower = control.getInt("Max virial power");

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
			virialPressuresList.clear();             // Z_virial
			residualPressuresList.clear();           // Z_residual
			totalSums.clear();                       // F_total/V
			chemicalPotList.clear();                 // mu/kT

			for (int i = 0; i < nStates; i++) {
				fExPerVolList.add(0.0);
				virialPressuresList.add(0.0);
				residualPressuresList.add(0.0);
				totalSums.add(0.0);
				chemicalPotList.add(0.0);
			}

			for (int i = 0; i < nStates; i++) {
				// Build free-energy and EOS quantities for every density point

				double rho = rhoList.get(i); // density at state i
				double B2 = secondVirialCoefficientList.get(i); // fixed B2

				// Z_virial = 1 + B2*rho + B3*rho^2 + B4*rho^3 + ...
				double Zvir = 1.0 + B2 * rho;
				for (int k = 0; k < virial.length; k++) { // Loop over B3, B4, ..., B_{maxPower+2}
					int n = k + 3; // n = 3 for B3, n=4 for B4, etc.
					Zvir += virial[k] * Math.pow(rho, n - 1);
				}
				virialPressuresList.set(i, Zvir); // Replace the 0.0 at index i

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

				// Residual between measured EOS and fitted virial EOS
				double Ztotal = meanPressures.get(i); // Measured from simulation
				double Zresidual = Ztotal - Zvir;

				residualPressuresList.set(i, Zresidual);
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

		// Reset scan state
		incrementDryVolFrac = true;
		dryVolFrac = 0.0;
		fixedB2 = 0.0;
		maxPower = 0;

		// Clear data from any previous run
		dryVolFracs.clear();
		totalSums.clear();
		chemicalPotList.clear();
		floryFperVolList.clear();
		idealFreeEnergyList.clear();
		fExPerVolList.clear();
		meanPressures.clear();
		uPairPerVolList.clear();
		virialPressuresList.clear();
		residualPressuresList.clear();

		reservoirVolFracList.clear();
		swellingRatioList.clear();
		totalVolList.clear();
		secondVirialCoefficientList.clear();
		reducedB2List.clear();
		hardSphereB2List.clear();
		volumefractionList.clear();
		virialCoefficientList.clear();

		rhoList.clear();
		yMinusB2List.clear();

		// Default simulation parameters
		control.setValue("DryVolFrac increment", 0.0001);
		control.setValue("DryVolFrac Max", 0.003);
		control.setValue("Initial configuration", "FCC");
		control.setValue("N", 108);
		control.setValue("Dry radius [nm]", 50);
		control.setValue("x-link fraction", 0.00003);
		control.setValue("Young's calibration", 1.0);
		control.setValue("chi", 0);
		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.00);
		control.setValue("Max virial power", 7);
		control.setValue("Delay", 10000);
		control.setValue("Snapshot interval", 100);
		control.setValue("Stop", 100000);
		control.setValue("Size bin width", 0.001);
		control.setValue("g(r) bin width", 0.005);
		control.setValue("Delta k", 0.005);
		control.setValue("Run number", "1");
		control.setValue("Calculate structure", false);
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
			File outputFile = new File(
				"data/Fluid_Phase_Free_Energy/Hertzian_Spheres/"
				+ "Hertzian_Spheres_maxPower"
				+ maxPower
				+ "_run"
				+ particles.fileExtension
				+ ".txt"
			);

			File outputDir = outputFile.getParentFile();

			if (outputDir != null && !outputDir.exists() && !outputDir.mkdirs()) {
				throw new IOException(
					"Could not create output directory: " + outputDir.getAbsolutePath()
				);
			}

			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}

			FileWriter fw1 = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);

			// write system parameters to the file in addittion to dryVolFracs and totalFree energies
			bw1.write("This fits B3 to B" + (maxPower + 2) + " into the virial expansion since maxPower is set to " + maxPower + ".");
			bw1.newLine();
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
			// bw1.write("The coupling constant increment - dlambda: " + particles.dlambda);
			// bw1.newLine();	
			
			// Column headers - simplified pressure columns
			bw1.write("phi0, phi, rho, mu/kT, Z_measured, Z_virial_Fit, Z_residual, <F_total>/V, QMelting, <F_id>/V, <F_ex>/V, <F_FR>/V, <alpha>, Zeta, B2_rho, B2_phi0, B2_HS, B2*");
			bw1.newLine();

			double factor = 0.75 / Math.PI; // Conversion factor

			for (int i = 1; i < dryVolFracs.size() - 1; i++) {

				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));

				double Z_measured = meanPressures.get(i);     // Measured from simulation (virial theorem)
				double Z_virial = virialPressuresList.get(i); // Calculated from virial expansion
				double Z_residual = residualPressuresList.get(i); // Residual between measured and fitted EOS
				double rho = rhoList.get(i);

				double FtotV = totalSums.get(i);
				double FidV = idealFreeEnergyList.get(i);
				double FexV = fExPerVolList.get(i);
				double FfrV = floryFperVolList.get(i);

				double qMelting = (uPairPerVolList.get(i) - FtotV + 1.50);

				double B2_rho = secondVirialCoefficientList.get(i);
				double B2_phi0 = B2_rho * factor; // Convert to phi0-space

				bw1.write(
					roundedDryVolFrac + ", " +
					volumefractionList.get(i) + ", " +
					rho + ", " +
					chemicalPotList.get(i) + ", " +
					Z_measured + ", " +
					Z_virial + ", " +
					Z_residual + ", " +
					FtotV + ", " +
					qMelting + ", " +
					FidV + ", " +
					FexV + ", " +
					FfrV + ", " +
					swellingRatioList.get(i) + ", " +
					reservoirVolFracList.get(i) + ", " +
					B2_rho + ", " +
					B2_phi0 + ", " +
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
			File virialFile = new File(
				"data/Fluid_Phase_Free_Energy/Hertzian_Spheres/"
				+ "VirialCoefficients_maxPower"
				+ maxPower
				+ "_run"
				+ particles.fileExtension
				+ ".txt"
			);

			// Create parent directory if it doesn't exist
			File virialDir = virialFile.getParentFile();
			if (virialDir != null && !virialDir.exists() && !virialDir.mkdirs()) {
				throw new IOException(
					"Could not create virial output directory: "
					+ virialDir.getAbsolutePath()
				);
			}
						
			if (!virialFile.exists() && !virialFile.createNewFile() && !virialFile.exists()) {
				throw new IOException("Could not create virial output file: " + virialFile.getAbsolutePath());
			}

			FileWriter fw2 = new FileWriter(virialFile.getAbsoluteFile());
			BufferedWriter bw2 = new BufferedWriter(fw2);

			// Write header with system parameters
			bw2.write("This fits B3 to B" + (maxPower + 2) + " into the virial expansion since maxPower is set to " + maxPower + ".");
			bw2.newLine();
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
			bw2.write("Number of density points: " + dryVolFracs.size());
			bw2.newLine();
			bw2.write("Dry volume fraction range: phi0_min = " + dryVolFracs.get(0) + ", phi0_max = " + dryVolFracs.get(dryVolFracs.size()-1));
			bw2.newLine();
			bw2.newLine();

			// Conversion factor
			double factor = 0.75 / Math.PI;

			// Write B2 values at each density point (in BOTH rho-space and phi0-space)
			bw2.write("===== B2 at Each Density Point =====");
			bw2.newLine();
			bw2.write("phi0, B2_rho, B2_phi0, B2_HS, B2*");
			bw2.newLine();
			for (int i = 0; i < dryVolFracs.size(); i++) {
				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));
				double B2_rho = secondVirialCoefficientList.get(i);
				double B2_phi0 = B2_rho * factor; // Convert to phi0-space
				
				bw2.write(roundedDryVolFrac + ", " + 
						B2_rho + ", " + 
						B2_phi0 + ", " + 
						hardSphereB2List.get(i) + ", " + 
						reducedB2List.get(i));
				bw2.newLine();
			}
			bw2.newLine();

			// Write higher-order virial coefficients IN RHO-SPACE (B3, B4, ..., B10)
			bw2.write("===== Virial Coefficients in RHO-SPACE (from OLS fit) =====");
			bw2.newLine();
			bw2.write("Coefficient, Value");
			bw2.newLine();
			for (int k = 0; k < virialCoefficientList.size(); k++) {
				int Bindex = k + 3;  // B3, B4, B5, ..., B10
				bw2.write("B" + Bindex + ", " + virialCoefficientList.get(k));
				bw2.newLine();
			}
			bw2.newLine();

			// *** NEW SECTION: Write coefficients in PHI0-SPACE ***
			bw2.write("===== Virial Coefficients in PHI0-SPACE =====");
			bw2.newLine();
			bw2.write("Conversion: rho = (0.75/pi) * phi0");
			bw2.newLine();
			bw2.write("B_n_tilde = B_n * (0.75/pi)^(n-1)");
			bw2.newLine();
			bw2.newLine();
			bw2.write("Coefficient, Value");
			bw2.newLine();
			
			double powerOfFactor = factor * factor; // Start with (0.75/pi)^2 for B3
			
			for (int k = 0; k < virialCoefficientList.size(); k++) {
				int Bindex = k + 3;
				double B_tilde = virialCoefficientList.get(k) * powerOfFactor;
				bw2.write("B" + Bindex + "_tilde, " + B_tilde);
				bw2.newLine();
				powerOfFactor *= factor; // Increment power for next coefficient
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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresNonLocalFacetFluidFixedTestApp());
	}
}
