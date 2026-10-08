/*
 * HertzianSpheresFluidFreeEnergyApp.java
 *
 * Monte Carlo calculation of the equation of state and Helmholtz free energy
 * of a fluid of fixed-size Hertzian spheres.
 *
 * The particles interact through
 *
 *   beta*u(r) = betaEpsilon*(1 - r/sigma)^(5/2),  r < sigma,
 *             = 0,                                r >= sigma.
 *
 * The equation of state is measured from the virial theorem. The excess
 * Helmholtz free energy per particle is obtained by thermodynamic integration:
 *
 *   beta*F_ex/N = integral from 0 to rho of [(Z(rho') - 1)/rho'] d rho',
 *
 * where Z = beta*P/rho is the compressibility factor.
 *
 * Thermodynamic quantities are calculated as
 *
 *   beta*F/N           = ln(rho) - 1 + beta*F_ex/N,
 *   beta*mu            = beta*F/N + Z,
 *   beta*P*sigma^3     = rhoStar*Z.
 *
 * Reduced units:
 *
 *   sigma              = 1,
 *   betaEpsilon        = 180, 160, 150, 140, 135, 130,
 *   TStar              = 1/betaEpsilon,
 *   rhoStar            = rho*sigma^3.
 *
 * This application is intended to reproduce the fixed-size Hertzian-sphere
 * fluid data reported by Pamies, Cacciuto, and Frenkel,
 * "Phase diagram of Hertzian spheres," arXiv:0811.2227.
 *
 * Authors: Oreoluwa Alade and Alan R. Denton
 */

package org.opensourcephysics.sip.Hertz;

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.Random;
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
import java.util.Locale;

/**
 * HertzSpheresApp performs a Monte Carlo simulation of microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @authors Alan Denton and Oreoluwa Alade
 * Last modified: 2026-08-24
 * 
 */

public class HertzianSpheresFluidFreeEnergyApp extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
    //HertzSpheresLiquidDensity particles = new HertzSpheresLiquidDensity();
	// HertzSpheresNonLocalFacet particles = new HertzSpheresNonLocalFacet();
	HertzianSpheresModel particles = new HertzianSpheresModel();
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
	double density;
	double densityMax;
	double densityIncrement;
	boolean firstDensity = true;
	double uPairPerVol;
	int maxPower;
	double dryVolFrac;
	double dryVolFracMax;
	boolean incrementDryVolFrac = true;
	double maxValidPhi = 0.9641;
	String outputDirectory;

	/* Temperatures and broad fluid pilot windows for the F--FCC boundary. */
	final double[] betaEpsilonValues = {180.0, 160.0, 150.0, 140.0, 135.0, 130.0};
	final double[] fluidDensityMaxes = {1.801, 1.901, 1.951, 2.001, 2.051, 2.101};
	int temperatureIndex = 0;
	boolean temperatureSeriesComplete = false;

	List<Double> dryVolFracs = new ArrayList<>();
    List<Double> totalSums = new ArrayList<>();
	List<Double> chemicalPotList = new ArrayList<>();
	List<Double> idealFreeEnergyList = new ArrayList<>();
	List<Double> fExPerVolList = new ArrayList<>();
	// Thermodynamic integration of the directly measured EOS.
	List<Double> fExMeasuredPerParticleList = new ArrayList<>();
	List<Double> fExMeasuredPerVolList = new ArrayList<>();
	List<Double> totalMeasuredPerVolList = new ArrayList<>();
	List<Double> chemicalPotMeasuredList = new ArrayList<>();
	List<Double> measuredReducedPressureList = new ArrayList<>();
	List<Double> calculatedPressures = new ArrayList<>();
	List<Double> meanPressures = new ArrayList<>();
	List<Double> uPairPerVolList = new ArrayList<>();
	List<Double> reservoirVolFracList = new ArrayList<>();
	List<Double> virialPressuresList = new ArrayList<>();
	List<Double> swellingRatioList = new ArrayList<>();
	List<Double> floryRehnerPressuresListEdited = new ArrayList<>();
	List<Double> totalVolList = new ArrayList<>();
	List<Double> secondVirialCoefficientList = new ArrayList<>();
	List<Double> reducedB2List = new ArrayList<>();
	List<Double> hardSphereB2List = new ArrayList<>();
	List<Double> volumefractionList = new ArrayList<>();
	List<Double> virialCoefficientList = new ArrayList<>();
	List<Double> rhoFitCurveList = new ArrayList<>();
	List<Double> zFitCurveList = new ArrayList<>();
	List<Double> fExFitCurveList = new ArrayList<>();
	List<Double> muExFitCurveList = new ArrayList<>();
	List<Double> pressureFitCurveList = new ArrayList<>();
	List<Double> phi0FitCurveList = new ArrayList<>();
	List<Double> phiFitCurveList = new ArrayList<>();
	List<Double> alphaFitCurveList = new ArrayList<>();
	// Virial Coefficient Lists
	List<Double> rhoList = new ArrayList<>();
	List<Double> yMinusB2List = new ArrayList<>(); // optional: y - B2
	// Free energy density contributions: f_id/V, f_FR/V, f_ex/V, f_total/V
	List<Double> fidFitCurveList  = new ArrayList<>();  // ideal gas free energy density
	List<Double> ffrFitCurveList  = new ArrayList<>();  // Flory-Rehner free energy density
	List<Double> fexFitCurveList2 = new ArrayList<>();  // excess free energy density from virial integration
	List<Double> ftotalFitCurveList = new ArrayList<>();  // total free energy density = F_id + F_FR + F_ex
	// Chemical potential contributions: mu = (4pi/3) * dF/dphi0
	List<Double> muIdFitCurveList  = new ArrayList<>(); // ideal gas contribution to mu
	List<Double> muFRFitCurveList  = new ArrayList<>(); // Flory-Rehner contribution to mu
	List<Double> muExFitCurveList2 = new ArrayList<>(); // excess contribution to mu
	// Pressure contributions: PV/NkT = (4pi/3)/phi0 * (phi0*dF/dphi0 - F_fit)
	List<Double> pvIdFitCurveList  = new ArrayList<>(); // ideal gas contribution to PV/NkT
	List<Double> pvFrFitCurveList  = new ArrayList<>(); // Flory-Rehner contribution to PV/NkT
	List<Double> pvExFitCurveList  = new ArrayList<>(); // excess contribution to PV/NkT

	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values

    /**
	* Initializes the model.
	*/
	public void initialize() {
		
		added = false;
		dryVolFracMax = fluidDensityMaxes[temperatureIndex];
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
			// dryVolFrac = particles.dphi;     // start at phi0_min = dphi0
			dryVolFrac = control.getDouble("DryVolFrac Start");
			particles.dryVolFrac = dryVolFrac;
			incrementDryVolFrac = false;
		}
		
		// Fixed-size Hertzian-sphere benchmark.
		particles.sigma = 1.0;
		particles.betaEpsilon = betaEpsilonValues[temperatureIndex];
		outputDirectory = currentOutputDirectory();
		File temperatureDirectory = new File(outputDirectory);
		if (!temperatureDirectory.exists()) {
			temperatureDirectory.mkdirs();
		}

		// Temporarily interpret the legacy dryVolFrac variable as rho*.
		particles.density = dryVolFrac;

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

	/**
	 * Integrates the directly measured compressibility factor:
	 *
	 *   beta F_ex/N = integral_0^rho [(Z(rho') - 1)/rho'] d rho'.
	 *
	 * The rho -> 0 endpoint is supplied by the exact second virial
	 * coefficient, since (Z - 1)/rho -> B2.  The measured state points
	 * are then integrated cumulatively with the trapezoidal rule.
	 */
	private void buildMeasuredEosThermodynamics() {
		fExMeasuredPerParticleList.clear();
		fExMeasuredPerVolList.clear();
		totalMeasuredPerVolList.clear();
		chemicalPotMeasuredList.clear();
		measuredReducedPressureList.clear();

		if (rhoList.isEmpty()) {
			return;
		}

		/*
		 * An absolute thermodynamic integration must begin in the dilute
		 * regime.  A high-density segment cannot be integrated directly from
		 * rho = 0 with one large trapezoid; it must first be joined to the
		 * preceding low-density EOS data in post-processing.
		 */
		if (rhoList.get(0) > 0.1) {
			System.err.println(
				"Measured-EOS free energy not calculated: the first density is "
				+ rhoList.get(0)
				+ ". Join this segment to the low-density EOS before integration."
			);
			for (int i = 0; i < rhoList.size(); i++) {
				fExMeasuredPerParticleList.add(Double.NaN);
				fExMeasuredPerVolList.add(Double.NaN);
				totalMeasuredPerVolList.add(Double.NaN);
				chemicalPotMeasuredList.add(Double.NaN);
				measuredReducedPressureList.add(
					rhoList.get(i) * meanPressures.get(i)
				);
			}
			return;
		}

		double previousRho = 0.0;
		double previousIntegrand = secondVirialCoefficientList.get(0);
		double cumulativeFExPerParticle = 0.0;

		for (int i = 0; i < rhoList.size(); i++) {
			double rho = rhoList.get(i);
			double zMeasured = meanPressures.get(i);

			if (rho <= previousRho) {
				throw new IllegalStateException(
					"Densities must be positive and strictly increasing for "
					+ "measured-EOS thermodynamic integration."
				);
			}

			double currentIntegrand = (zMeasured - 1.0) / rho;
			cumulativeFExPerParticle += 0.5
					* (previousIntegrand + currentIntegrand)
					* (rho - previousRho);

			double fExMeasuredV = rho * cumulativeFExPerParticle;
			double fFrV = (particles.atolerance == 0.0)
					? 0.0
					: floryRehnerFreeEnergy(swellingRatioList.get(i), rho);
			double fTotalMeasuredV = idealFreeEnergyList.get(i)
					+ fFrV
					+ fExMeasuredV;
			double muMeasured = fTotalMeasuredV / rho + zMeasured;

			fExMeasuredPerParticleList.add(cumulativeFExPerParticle);
			fExMeasuredPerVolList.add(fExMeasuredV);
			totalMeasuredPerVolList.add(fTotalMeasuredV);
			chemicalPotMeasuredList.add(muMeasured);
			measuredReducedPressureList.add(rho * zMeasured);

			previousRho = rho;
			previousIntegrand = currentIntegrand;
		}
	}

	private double floryRehnerFreeEnergy(double alpha, double rho) {
		/**
		 * Computes the Flory-Rehner free energy density f_FR/V analytically
		 * from the swelling ratio alpha, referenced to the reservoir state.
		 *
		 * Formula (per particle, then converted to per volume):
		 *   f_mix    = nMon * [(alpha^3 - 1)*ln(1 - 1/alpha^3) + chi*(1 - 1/alpha^3)]
		 *   f_el     = 1.5 * nChains * (alpha^2 - ln(alpha) - 1)
		 *   f_FR     = f_mix + f_el - totalFRSR   (reservoir reference subtracted)
		 *   f_FR/V   = f_FR * rho                 (per volume = per particle * number density)
		 */

		double alpha3 = alpha * alpha * alpha;

		// Mixing contribution (Flory-Huggins)
		double mixF = particles.nMon * (
			(alpha3 - 1.0) * Math.log(1.0 - 1.0 / alpha3)
			+ particles.chi * (1.0 - 1.0 / alpha3)
		);

		// Elastic contribution (Gaussian network)
		double elasticF = 1.5 * particles.nChains * (alpha * alpha - Math.log(alpha) - 1.0);

		// Reservoir state reference (single particle in solvent, infinite dilution)
		double reservoirSr = particles.reservoirSR;
		double SR3 = reservoirSr * reservoirSr * reservoirSr;
		double mixFRSR = particles.nMon * (
			(SR3 - 1.0) * Math.log(1.0 - 1.0 / SR3)
			+ particles.chi * (1.0 - 1.0 / SR3)
		);
		double elasticFRSR = 1.5 * particles.nChains * (reservoirSr * reservoirSr - Math.log(reservoirSr) - 1.0);
		double totalFRSR = mixFRSR + elasticFRSR;

		// Per-particle FR free energy referenced to reservoir
		double fFRperParticle = mixF + elasticF - totalFRSR;

		// Convert to per-volume
		return fFRperParticle * rho;
	}

	private void buildVirialFitCurve(double[] virial, int nFitPoints) {

		/**
		 * Builds the virial fit curve and computes thermodynamic quantities
		 * (free energy, chemical potential, pressure) for the fluid phase.
		 *
		 * Steps:
		 *   Step 0 - fit polynomial to alpha vs phi0 from discrete simulation points
		 *   Step 1 - build fine 500-point grid; compute f_id, f_FR (analytic), f_ex at each point
		 *   Step 2 - fit separate polynomials to each free energy contribution vs phi0
		 *   Step 3 - differentiate each polynomial analytically to get mu and PV/NkT
		 *
		 * Key change from previous version:
		 *   f_FR is now computed analytically from the Flory-Rehner formula using a
		 *   smooth polynomial alpha(phi0), rather than linearly interpolating the
		 *   discrete simulation values of meanFreeEnergy(). This gives a physically
		 *   consistent and smooth f_FR on the fine grid.
		 */

		// Clear all output lists before filling them
		rhoFitCurveList.clear();
		zFitCurveList.clear();
		fExFitCurveList.clear();
		muExFitCurveList.clear();
		pressureFitCurveList.clear();
		phi0FitCurveList.clear();
		phiFitCurveList.clear();
		alphaFitCurveList.clear();
		fidFitCurveList.clear();
		ffrFitCurveList.clear();
		fexFitCurveList2.clear();
		ftotalFitCurveList.clear();
		muIdFitCurveList.clear();
		muFRFitCurveList.clear();
		muExFitCurveList2.clear();
		pvIdFitCurveList.clear();
		pvFrFitCurveList.clear();
		pvExFitCurveList.clear();

		double rhoMin  = rhoList.get(0);
		double rhoMax  = rhoList.get(rhoList.size() - 1);
		// In the fixed-size benchmark the scan variable is rho* itself.
		double factor  = 1.0;

		// ===========================
		// Step 0: fit polynomial to alpha vs phi0 from discrete simulation points
		// This gives a smooth alpha(phi0) to use in the Flory-Rehner formula
		// ===========================
		int nDisc = rhoList.size(); // number of discrete simulation points (~28)

		double[] phi0Disc  = new double[nDisc];
		double[] alphaDisc = new double[nDisc];
		double[] phiDisc   = new double[nDisc];

		for (int i = 0; i < nDisc; i++) {
			phi0Disc[i]  = rhoList.get(i);                  // legacy name; value is rho*
			alphaDisc[i] = swellingRatioList.get(i);        // alpha from simulation
			phiDisc[i]   = volumefractionList.get(i);       // phi from simulation
		}

		// Fit polynomial alpha(phi0) to discrete simulation points
		double[] coeffsAlpha = fitPoly(alphaDisc, phi0Disc, nDisc, maxPower);
		double[] coeffsPhi   = fitPoly(phiDisc,   phi0Disc, nDisc, maxPower);

		// ===========================
		// Step 1: build fine grid arrays
		// Compute f_id, f_FR (analytic), f_ex at each of the 500 grid points
		// ===========================
		double[] phi0Arr = new double[nFitPoints];
		double[] phiArr  = new double[nFitPoints];
		double[] rhoArr  = new double[nFitPoints];
		double[] zvirArr = new double[nFitPoints];
		double[] FidArr  = new double[nFitPoints];
		double[] FfrArr  = new double[nFitPoints];
		double[] FexArr  = new double[nFitPoints];
		double[] FtotArr = new double[nFitPoints];
		double fixedB2 = particles.secondVirialCoefficient();

		for (int i = 0; i < nFitPoints; i++) {

			// Uniformly spaced density grid from rhoMin to rhoMax
			double rho  = rhoMin + (rhoMax - rhoMin) * i / (nFitPoints - 1);
			double phi0 = rho; // legacy output name; value is rho*

			// Smooth alpha from polynomial fit to discrete simulation points
			double alpha = evalPoly(coeffsAlpha, phi0, maxPower);

			// Smooth phi from polynomial fit
			double phi = evalPoly(coeffsPhi, phi0, maxPower);

			// B2 computed analytically from smooth alpha
			double B2 = fixedB2;

			// Z_virial: Z = 1 + B2*rho + B3*rho^2 + ...
			double Zvir = 1.0 + B2 * rho;
			for (int k = 0; k < virial.length; k++) {
				int nn = k + 3;
				Zvir += virial[k] * Math.pow(rho, nn - 1);
			}

			// Excess free energy density from integrating virial EOS:
			// f_ex = B2*rho^2 + B3*rho^3/2 + B4*rho^4/3 + ...
			double FexV = B2 * rho * rho;
			for (int k = 0; k < virial.length; k++) {
				int nn = k + 3;
				FexV += virial[k] * Math.pow(rho, nn) / (nn - 1.0);
			}

			// Ideal-gas free-energy density, including the finite-N correction.
			double stirling = rho
							* Math.log(2.0 * Math.PI * particles.N) / (2.0 * particles.N);
			double FidV = rho * (Math.log(rho) - 1.0) + stirling;

			// Total free energy density
			// double FfrV = floryRehnerFreeEnergy(alpha, rho);
			double FfrV = (particles.atolerance == 0.0) ? 0.0 : floryRehnerFreeEnergy(alpha, rho);
			double FtotV = FidV + FfrV + FexV;

			// Store in arrays for polynomial fitting in Step 2
			phi0Arr[i] = phi0;
			phiArr[i]  = phi;
			rhoArr[i]  = rho;
			zvirArr[i] = Zvir;
			FidArr[i]  = FidV;
			FfrArr[i]  = FfrV;
			FexArr[i]  = FexV;
			FtotArr[i] = FtotV;

			// Store Z and structural quantities for output
			rhoFitCurveList.add(rho);
			zFitCurveList.add(Zvir);
			fExFitCurveList.add(FexV);
			phi0FitCurveList.add(phi0);
			phiFitCurveList.add(phi);
			alphaFitCurveList.add(alpha);
		}

		/*
		 * Evaluate the virial thermodynamics directly.  maxPower controls only
		 * the number of fitted virial coefficients.  It must not also be used
		 * to refit rho*log(rho) and the completed free-energy curve with a
		 * low-degree polynomial.  Any smoothing degree needed for coexistence
		 * is a separate post-processing choice.
		 */
		double B2 = fixedB2;
		double finiteNCorrection =
			Math.log(2.0 * Math.PI * particles.N) / (2.0 * particles.N);

		for (int i = 0; i < nFitPoints; i++) {
			double rho = rhoArr[i];
			double fId = FidArr[i];
			double F_fr = FfrArr[i];
			double F_ex = FexArr[i];
			double F_tot = FtotArr[i];

			double mu_id = Math.log(rho) + finiteNCorrection;
			double mu_fr = 0.0; // Flory-Rehner is absent in this benchmark.
			double mu_ex = 2.0 * B2 * rho;
			double PV_ex = B2 * rho;

			for (int k = 0; k < virial.length; k++) {
				int n = k + 3;
				double rhoPower = Math.pow(rho, n - 1);
				mu_ex += (n / (n - 1.0)) * virial[k] * rhoPower;
				PV_ex += virial[k] * rhoPower;
			}

			double mu_tot = mu_id + mu_fr + mu_ex;
			double PV_id = 1.0;
			double PV_fr = 0.0;
			double PV_tot = PV_id + PV_fr + PV_ex;

			fidFitCurveList.add(fId);
			ffrFitCurveList.add(F_fr);
			fexFitCurveList2.add(F_ex);
			ftotalFitCurveList.add(F_tot);

			muIdFitCurveList.add(mu_id);
			muFRFitCurveList.add(mu_fr);
			muExFitCurveList2.add(mu_ex);
			muExFitCurveList.add(mu_tot);

			pvIdFitCurveList.add(PV_id);
			pvFrFitCurveList.add(PV_fr);
			pvExFitCurveList.add(PV_ex);
			pressureFitCurveList.add(PV_tot);
		}
	}

	private double[] fitPoly(double[] Y, double[] phi0Arr, int nPts, int degree) {
		/**
		 * Fits a polynomial of given degree to Y vs phi0Arr using OLS.
		 *
		 * The polynomial has the form:
		 *      Y(phi0) = c0 + c1*phi0 + c2*phi0^2 + ... + cn*phi0^n
		 *
		 * Input:
		 *      Y        = array of free energy values to fit (e.g. FidArr, FfrArr, FexArr)
		 *      phi0Arr  = array of dry volume fraction values (x-axis)
		 *      nPts     = number of data points (500 in our case)
		 *      degree   = polynomial degree (equals maxPower)
		 *
		 * Returns:
		 *      coeffs = [c0, c1, c2, ..., cn] where c0 is the intercept
		 */

		// Build the design matrix X where column j contains phi0^(j+1)
		double[][] X = new double[nPts][degree];
		for (int i = 0; i < nPts; i++) {
			double p = phi0Arr[i];
			for (int j = 0; j < degree; j++) {
				X[i][j] = Math.pow(p, j + 1); // phi0^1, phi0^2, ..., phi0^degree
			}
		}

		// Fit polynomial using OLS
		OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
		ols.setNoIntercept(false); // include constant term c0
		ols.newSampleData(Y, X);

		// Return fitted coefficients [c0, c1, c2, ..., cn]
		return ols.estimateRegressionParameters();
	}

	private double evalPoly(double[] coeffs, double p, int degree) {
		/**
		 * Evaluates a fitted polynomial at a given phi0 value.
		 *
		 * Computes:
		 *      f(phi0) = c0 + c1*phi0 + c2*phi0^2 + ... + cn*phi0^n
		 *
		 * Input:
		 *      coeffs = polynomial coefficients [c0, c1, ..., cn] from fitPoly()
		 *      p      = phi0 value at which to evaluate
		 *      degree = polynomial degree
		 *
		 * Returns:
		 *      the free energy value at phi0 = p
		 */

		double val = coeffs[0]; // start with intercept c0
		for (int j = 0; j < degree; j++) {
			val += coeffs[j + 1] * Math.pow(p, j + 1); // add c1*p + c2*p^2 + ...
		}
		return val;
	}

	private double evalPolyDeriv(double[] coeffs, double p, int degree) {
		/**
		 * Evaluates the analytical derivative of a fitted polynomial at phi0.
		 *
		 * Computes:
		 *      dF/dphi0 = c1 + 2*c2*phi0 + 3*c3*phi0^2 + ... + n*cn*phi0^(n-1)
		 *
		 * This is used to compute mu and PV/NkT analytically without
		 * any numerical differentiation error.
		 *
		 * Input:
		 *      coeffs = polynomial coefficients [c0, c1, ..., cn] from fitPoly()
		 *      p      = phi0 value at which to evaluate the derivative
		 *      degree = polynomial degree
		 *
		 * Returns:
		 *      dF/dphi0 at phi0 = p
		 */

		double deriv = 0.0; // c0 vanishes under differentiation
		for (int j = 0; j < degree; j++) {
			deriv += coeffs[j + 1] * (j + 1) * Math.pow(p, j); // (j+1)*c_{j+1}*p^j
		}
		return deriv;
	}
		
	public void doStep() {
		if (temperatureSeriesComplete) {
			return;
		}
		particles.step();

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

			double phi = particles.meanVolFrac();

			if (dryVolFrac < dryVolFracMax){

				double rho   = particles.N / particles.totalVol;
				double alpha = particles.meanRadius();

				// measured EOS
				double Z = particles.meanPressure();
				double reducedPressure = rho * Z;
				
				System.out.println(
					"rho* = " + rho
					+ " | Z = " + Z
					+ " | P* = " + reducedPressure
				);

				System.out.println("phi0 = " + dryVolFrac
					+ " | phi = " + phi
					+ " | Z_measured = " + Z
					+ " | rho = " + rho
					+ " | alpha = " + alpha);

				// Compute B2 
				// B2 depends on alpha and is evaluated separately at each state point
				double B2 = particles.secondVirialCoefficient();
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
				volumefractionList.add(phi);
				swellingRatioList.add(alpha);
				dryVolFracs.add(dryVolFrac);

				// store FR per volume (state specific)

				// System.out.println("phi0 = " + dryVolFrac);
				// System.out.println("meanFreeEnergy per particle = " + particles.meanFreeEnergy());
				// System.out.println("N/totalVol = " + (particles.N / particles.totalVol));
				// System.out.println("F_FR per volume = " + floryFperVol);

				// store Ideal term 
				double stirlingApprox =
						rho * Math.log(2.0 * Math.PI * particles.N)
						/ (2.0 * particles.N);

				double idealFreeEnergy =
						rho * (Math.log(rho) - 1.0)
						+ stirlingApprox;
				idealFreeEnergyList.add(idealFreeEnergy);

				// uPair/V (used this in QMelting)
				uPairPerVol = particles.meanPairEnergy() * (particles.N / particles.totalVol);
				uPairPerVolList.add(uPairPerVol);

				// Increase the reduced density and rescale the current configuration.
				dryVolFrac += particles.dphi;

				if (dryVolFrac < dryVolFracMax) {
					particles.changeDensity(dryVolFrac);
					return;
				}

			}

			// ===========================
			// Scan complete: fit virial coefficients
			// ===========================
			maxPower = control.getInt("Max virial power");

			double[] virial = fitVirialNoIntercept(maxPower); // returns [B3, B4, ..., B_{maxPower+2}]

			buildVirialFitCurve(virial, 500);

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
			floryRehnerPressuresListEdited.clear();  // Z_FR = Z_total - Z_virial
			calculatedPressures.clear();             // reconstructed Z_total = Z_virial + Z_FR
			totalSums.clear();                       // F_total/V
			chemicalPotList.clear();                 // mu/kT

			for (int i = 0; i < nStates; i++) {
				fExPerVolList.add(0.0);
				virialPressuresList.add(0.0);
				floryRehnerPressuresListEdited.add(0.0);
				calculatedPressures.add(0.0);
				totalSums.add(0.0);
				chemicalPotList.add(0.0);
			}

			for (int i = 0; i < nStates; i++) {

				double rho = rhoList.get(i); // density at state i
				double B2 = secondVirialCoefficientList.get(i); // exact B2 at state i 

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
				double alphaI        = swellingRatioList.get(i);
				double rhoI        = rhoList.get(i);
				double fFrAnalytic = (particles.atolerance == 0.0) ? 0.0 : floryRehnerFreeEnergy(alphaI, rhoI);
				double fTotal       = fFrAnalytic + idealFreeEnergyList.get(i) + fEx;
				totalSums.set(i, fTotal);

				// FR contribution to Z (in Z-units)
				double zTotal = meanPressures.get(i); // Measured from simulation
				double zFr = zTotal - Zvir; // Flory-Rehner contribution

				floryRehnerPressuresListEdited.set(i, zFr);
				calculatedPressures.set(i, Zvir + zFr); // should ~ Ztotal (fit error)
			}

			// ===========================
			// Chemical potential from fitted virial EOS
			// ideal + excess analytic, FR from derivative of F_FR/V
			// ===========================
			for (int i = 0; i < nStates; i++) {

				double rho = rhoList.get(i);
				double B2 = secondVirialCoefficientList.get(i);

				// Ideal contribution
				double muIdeal = Math.log(rho)
						+ Math.log(2.0 * Math.PI * particles.N) / (2.0 * particles.N);

				// Excess (virial) contribution
				double muEx = 2.0 * B2 * rho;

				for (int k = 0; k < virial.length; k++) {
					int n = k + 3;
					muEx += (n / (n - 1.0)) * virial[k] * Math.pow(rho, n - 1);
				}

				// Flory-Rehner contribution (numerical derivative of F_FR/V)
				// double fFRNext = floryRehnerFreeEnergy(swellingRatioList.get(i + 1), rhoList.get(i + 1));
				// double fFRPrev = floryRehnerFreeEnergy(swellingRatioList.get(i - 1), rhoList.get(i - 1));
				// double muFR = (fFRNext - fFRPrev) / (rhoNext - rhoPrev);

				double muFR;
				if (particles.atolerance == 0.0) {
					muFR = 0.0;
				} else if (i == 0) {
					muFR = (floryRehnerFreeEnergy(swellingRatioList.get(1), rhoList.get(1)) -
							floryRehnerFreeEnergy(swellingRatioList.get(0), rhoList.get(0))) /
							(rhoList.get(1) - rhoList.get(0));
				} else if (i == nStates - 1) {
					muFR = (floryRehnerFreeEnergy(swellingRatioList.get(i), rhoList.get(i)) -
							floryRehnerFreeEnergy(swellingRatioList.get(i - 1), rhoList.get(i - 1))) /
							(rhoList.get(i) - rhoList.get(i - 1));
				} else {
					muFR = (floryRehnerFreeEnergy(swellingRatioList.get(i + 1), rhoList.get(i + 1)) -
							floryRehnerFreeEnergy(swellingRatioList.get(i - 1), rhoList.get(i - 1))) /
							(rhoList.get(i + 1) - rhoList.get(i - 1));
				}

				double muTotal = muIdeal + muEx + muFR;

				chemicalPotList.set(i, muTotal);
			}

			// Independent thermodynamic-integration route based directly on
			// the measured compressibility factors.  This does not alter the
			// virial-fit route above.
			buildMeasuredEosThermodynamics();

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

		display3d.setMessage("Number of steps: " + particles.steps);

		if (control.getBoolean("Visualization on")) {
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
			+ betaLabel + "_Tstar_" + temperatureLabel + "/fluid/";
	}

	private boolean advanceToNextTemperature() {
		if (temperatureIndex + 1 >= betaEpsilonValues.length) {
			return false;
		}

		temperatureIndex++;
		clearTemperatureData();
		dryVolFrac = control.getDouble("DryVolFrac Start");
		dryVolFracMax = fluidDensityMaxes[temperatureIndex];
		particles.dryVolFrac = dryVolFrac;
		particles.density = dryVolFrac;
		particles.betaEpsilon = betaEpsilonValues[temperatureIndex];
		outputDirectory = currentOutputDirectory();
		incrementDryVolFrac = false;

		control.println("");
		control.println("Starting betaEpsilon = " + particles.betaEpsilon
			+ ", T* = " + (1.0 / particles.betaEpsilon));
		control.println("Fluid density window = " + dryVolFrac
			+ " to " + dryVolFracMax);
		return true;
	}

	private void clearTemperatureData() {
		dryVolFracs.clear();
		totalSums.clear();
		chemicalPotList.clear();
		idealFreeEnergyList.clear();
		fExPerVolList.clear();
		fExMeasuredPerParticleList.clear();
		fExMeasuredPerVolList.clear();
		totalMeasuredPerVolList.clear();
		chemicalPotMeasuredList.clear();
		measuredReducedPressureList.clear();
		calculatedPressures.clear();
		meanPressures.clear();
		uPairPerVolList.clear();
		reservoirVolFracList.clear();
		virialPressuresList.clear();
		swellingRatioList.clear();
		floryRehnerPressuresListEdited.clear();
		totalVolList.clear();
		secondVirialCoefficientList.clear();
		reducedB2List.clear();
		hardSphereB2List.clear();
		volumefractionList.clear();
		virialCoefficientList.clear();
		rhoFitCurveList.clear();
		zFitCurveList.clear();
		fExFitCurveList.clear();
		muExFitCurveList.clear();
		pressureFitCurveList.clear();
		phi0FitCurveList.clear();
		phiFitCurveList.clear();
		alphaFitCurveList.clear();
		rhoList.clear();
		yMinusB2List.clear();
		fidFitCurveList.clear();
		ffrFitCurveList.clear();
		fexFitCurveList2.clear();
		ftotalFitCurveList.clear();
		muIdFitCurveList.clear();
		muFRFitCurveList.clear();
		muExFitCurveList2.clear();
		pvIdFitCurveList.clear();
		pvFrFitCurveList.clear();
		pvExFitCurveList.clear();
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

		control.setValue("DryVolFrac Start", 0.05);
		control.setValue("DryVolFrac increment", 0.05);
		control.setValue("DryVolFrac Max", fluidDensityMaxes[0]);

		control.setValue("Initial configuration", "FCC");
		control.setValue("N", 108);

		control.setValue("Dry radius [nm]", 50);
		control.setValue("x-link fraction", 0.00003);
		control.setValue("Young's calibration", 1.0);
		control.setValue("chi", 0);

		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0);

		// maxPower = 1 fits B3 only;
		// maxPower = 2 fits B3 and B4, etc.
		control.setValue("Max virial power", 1);

		control.setValue("Delay", 10000);
		control.setValue("Snapshot interval", 100);
		control.setValue("Stop", 20000);

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
			File systemInfo = new File(
				outputDirectory
				+ "systemInfo_run"
				+ particles.fileExtension
				+ ".txt"
			);
		
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
			File sizeFile = new File(
				outputDirectory
				+ "microgelSize_run"
				+ particles.fileExtension
				+ ".txt"
			);
		
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
					outputDirectory
					+ "Fluid_EOS_maxPower"
					+ maxPower
					+ "_run"
					+ particles.fileExtension
					+ ".txt"
			);

			// Create parent directory if it does not exist
			File outputDir = outputFile.getParentFile();
			if (outputDir != null && !outputDir.exists()) {
				outputDir.mkdirs();
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
			bw1.write("betaEpsilon: " + particles.betaEpsilon);
			bw1.newLine();
			bw1.write("Tstar: " + (1.0 / particles.betaEpsilon));
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
			bw1.write("phi0, phi, rho, mu_virialFit/kT, Z_measured, Z_virial_Fit, <F_total_virialFit>/V, QMelting, <F_id>/V, <F_ex_virialFit>/V, <F_FR>/V, <alpha>, Zeta, B2_rho, B2_phi0, B2_HS, B2*, <F_ex_virialFit>/N, <F_ex_measuredTI>/N, mu_measuredTI/kT, P_virialFit*, P_measured*, <F_ex_measuredTI>/V, <F_total_measuredTI>/V");
			bw1.newLine();

			double factor = 1.0; // In this benchmark the legacy phi0 field equals rho*.

			for (int i = 0; i < dryVolFracs.size(); i++) {

				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));

				double Z_measured = meanPressures.get(i);     // Measured from simulation (virial theorem)
				double rho = rhoList.get(i);
				double Z_virial = virialPressuresList.get(i); // Calculated from virial expansion

				double FtotV = totalSums.get(i);
				double FidV = idealFreeEnergyList.get(i);
				double FexV = fExPerVolList.get(i);
				// double FfrV = floryRehnerFreeEnergy(swellingRatioList.get(i), rhoList.get(i));
				double FfrV = (particles.atolerance == 0.0) ? 0.0 : 
					floryRehnerFreeEnergy(swellingRatioList.get(i), rhoList.get(i));

				double qMelting = (uPairPerVolList.get(i) - FtotV + 1.50);

				double B2_rho = secondVirialCoefficientList.get(i);
				double B2_phi0 = B2_rho * factor; // Convert to phi0-space
				double fExVirialPerParticle = FexV / rho;
				double pressureVirialFit = rho * Z_virial;

				bw1.write(
					roundedDryVolFrac + ", " +
					volumefractionList.get(i) + ", " +
					rho + ", " +
					chemicalPotList.get(i) + ", " +
					Z_measured + ", " +
					Z_virial + ", " +
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
					reducedB2List.get(i) + ", " +
					fExVirialPerParticle + ", " +
					fExMeasuredPerParticleList.get(i) + ", " +
					chemicalPotMeasuredList.get(i) + ", " +
					pressureVirialFit + ", " +
					measuredReducedPressureList.get(i) + ", " +
					fExMeasuredPerVolList.get(i) + ", " +
					totalMeasuredPerVolList.get(i)
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
					outputDirectory
					+ "VirialCoefficients_maxPower"
					+ maxPower
					+ "_run"
					+ particles.fileExtension
					+ ".txt"
			);
			
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

			// In this benchmark the legacy phi0 field equals rho*.
			double factor = 1.0;

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
			
				double powerOfFactor = factor * factor; // Unity for rho* = legacy phi0.
			
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

		try {
				File fitCurveFile = new File(
				outputDirectory
				+ "VirialFitCurve_maxPower"
				+ maxPower
				+ "_run"
				+ particles.fileExtension
				+ ".txt"
			);

			File fitCurveDir = fitCurveFile.getParentFile();
			if (fitCurveDir != null && !fitCurveDir.exists()) {
				fitCurveDir.mkdirs();
			}

			if (!fitCurveFile.exists()) {
				fitCurveFile.createNewFile();
			}

			FileWriter fwFit = new FileWriter(fitCurveFile.getAbsoluteFile());
			BufferedWriter bwFit = new BufferedWriter(fwFit);

			bwFit.write("phi0, phi, rho, Z_fit_curve, " +
            "F_id, F_fr, F_ex, F_total, " +
            "mu_id, mu_fr, mu_ex, mu_kT, " +
            "PV_id, PV_fr, PV_ex, PV_NkT, " +
            "alpha");
			bwFit.newLine();

			for (int i = 0; i < rhoFitCurveList.size(); i++) {
			bwFit.write(
				phi0FitCurveList.get(i)    + ", " + // dry volume fraction
				phiFitCurveList.get(i)     + ", " + // swollen volume fraction
				rhoFitCurveList.get(i)     + ", " + // number density
				zFitCurveList.get(i)       + ", " + // Z from virial EOS
				fidFitCurveList.get(i)     + ", " + // ideal free energy density
				ffrFitCurveList.get(i)     + ", " + // Flory-Rehner free energy density
				fexFitCurveList2.get(i)    + ", " + // excess free energy density
				ftotalFitCurveList.get(i)    + ", " + // total free energy density
				muIdFitCurveList.get(i)    + ", " + // ideal contribution to mu
				muFRFitCurveList.get(i)    + ", " + // FR contribution to mu
				muExFitCurveList2.get(i)   + ", " + // excess contribution to mu
				muExFitCurveList.get(i)    + ", " + // total mu
				pvIdFitCurveList.get(i)    + ", " + // ideal contribution to PV/NkT
				pvFrFitCurveList.get(i)    + ", " + // FR contribution to PV/NkT
				pvExFitCurveList.get(i)    + ", " + // excess contribution to PV/NkT
				pressureFitCurveList.get(i) + ", " + // total PV/NkT
				alphaFitCurveList.get(i)            // swelling ratio
			);
			bwFit.newLine();
		}

			bwFit.close();
		}
		catch (IOException e) {
			e.printStackTrace();
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
			SimulationControl control = SimulationControl.createApp(new HertzianSpheresFluidFreeEnergyApp());
	}
}
