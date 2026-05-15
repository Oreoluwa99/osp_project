/*
	This is the second part of the HertzSpheresLiquidVirialCoefficientApp.java file, which is the main
	application class for performing a Monte Carlo simulation of microgels interacting via the Hertz
	elastic pair potential and swelling according to the Flory-Rehner free energy. This version includes
	methods for fitting higher-order virial coefficients from the simulation data, building the virial
	fit curve, and computing thermodynamic quantities such as the free energy, chemical potential, and
	pressure for the fluid phase. The code also includes methods for polynomial fitting and interpolation
	to analyze the simulation results and extract the virial coefficients.
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

/**
 * HertzSpheresApp performs a Monte Carlo simulation of microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @authors Alan Denton and Oreoluwa Alade
 * @version 1.2 13-05-2026
 * 
 */

public class HertzSpheresLiquidVirialCoefficient_V2App extends AbstractSimulation {
    public enum WriteModes{WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
    //HertzSpheresLiquidDensity particles = new HertzSpheresLiquidDensity();
	HertzSpheresNonLocalFacet particles = new HertzSpheresNonLocalFacet();
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
	double lambda=1;
	double dryVolFrac, dryVolFracMax;
	boolean incrementDryVolFrac = true;
	boolean flag = false;
	double totalVol;
	double floryFperVol;
	double pairEnergy = 0;
	double n;
	double uPairPerVol;
	int maxPower;

	List<Double> dryVolFracs = new ArrayList<>();
    List<Double> totalSums = new ArrayList<>();
	List<Double> chemicalPotList = new ArrayList<>();
	List<Double> floryFperVolList = new ArrayList<>();
	List<Double> idealFreeEnergyList = new ArrayList<>();
	List<Double> fExPerVolList = new ArrayList<>();
	List<Double> calculatedPressures = new ArrayList<>();
	List<Double> meanPressures = new ArrayList<>();
	List<Double> uPairPerVolList = new ArrayList<>();
	List<Double> varialPressuresList = new ArrayList<>();
	List<Double> reservoirVolFracList = new ArrayList<>();
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
	List<Double> FidFitCurveList  = new ArrayList<>();  // ideal gas free energy density
	List<Double> FfrFitCurveList  = new ArrayList<>();  // Flory-Rehner free energy density
	List<Double> FexFitCurveList2 = new ArrayList<>();  // excess free energy density from virial integration
	List<Double> FtotFitCurveList = new ArrayList<>();  // total free energy density = F_id + F_FR + F_ex

	// Chemical potential contributions: mu = (4pi/3) * dF/dphi0
	List<Double> muIdFitCurveList  = new ArrayList<>(); // ideal gas contribution to mu
	List<Double> muFRFitCurveList  = new ArrayList<>(); // Flory-Rehner contribution to mu
	List<Double> muExFitCurveList2 = new ArrayList<>(); // excess contribution to mu

	// Pressure contributions: PV/NkT = (4pi/3)/phi0 * (phi0*dF/dphi0 - F_fit)
	List<Double> PVIdFitCurveList  = new ArrayList<>(); // ideal gas contribution to PV/NkT
	List<Double> PVFRFitCurveList  = new ArrayList<>(); // Flory-Rehner contribution to PV/NkT
	List<Double> PVExFitCurveList  = new ArrayList<>(); // excess contribution to PV/NkT

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

	private double interpolateB2(double rho) {

		/**
		 * Performs linear interpolation to estimate the second virial coefficient B2
		 * at an arbitrary density rho using precomputed discrete data points.
		 *
		 * Mathematical form:
		 *   B2(rho) = B21 + ((rho - r1) / (r2 - r1)) * (B22 - B21)
		 *
		 * where:
		 *   (r1, B21) and (r2, B22) are consecutive data points such that r1 ≤ rho ≤ r2.
		 */

		int n = rhoList.size();  // total number of available (rho, B2) data points

		// If rho is below the smallest sampled density, return the first B2 value
		if (rho <= rhoList.get(0)) {
			return secondVirialCoefficientList.get(0);
		}

		// If rho is above the largest sampled density, return the last B2 value
		if (rho >= rhoList.get(n - 1)) {
			return secondVirialCoefficientList.get(n - 1);
		}

		// Iterate through consecutive pairs (r1, r2) to find where rho lies
		for (int i = 0; i < n - 1; i++) {

			double r1 = rhoList.get(i);       // lower bound density
			double r2 = rhoList.get(i + 1);   // upper bound density

			// Check if rho is within this interval
			if (rho >= r1 && rho <= r2) {

				double B21 = secondVirialCoefficientList.get(i);     // B2 at r1
				double B22 = secondVirialCoefficientList.get(i + 1); // B2 at r2

				// --- Linear interpolation ---
				// Compute B2(rho) assuming a linear variation between r1 and r2
				return B21 + (rho - r1) * (B22 - B21) / (r2 - r1);
			}
		}

		// Return the last value as a safety measure
		return secondVirialCoefficientList.get(n - 1);
	}

	private void buildVirialFitCurve(double[] virial, int nFitPoints) {

		/**
		 * Builds the virial fit curve and computes thermodynamic quantities
		 * (free energy, chemical potential, pressure) for the fluid phase.
		 *
		 * Three steps:
		 *   Step 1 - build fine grid of 500 points and compute f_id, f_FR, f_ex at each point
		 *   Step 2 - fit separate polynomials to each free energy contribution
		 *   Step 3 - differentiate each polynomial analytically to get mu and PV/NkT
		 *
		 * This separation allows us to see exactly which contribution (ideal, FR, excess)
		 * dominates the chemical potential and pressure at each density.
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

		FidFitCurveList.clear();
		FfrFitCurveList.clear();
		FexFitCurveList2.clear();
		FtotFitCurveList.clear();
		muIdFitCurveList.clear();
		muFRFitCurveList.clear();
		muExFitCurveList2.clear();
		PVIdFitCurveList.clear();
		PVFRFitCurveList.clear();
		PVExFitCurveList.clear();

		double rhoMin = rhoList.get(0);
		double rhoMax = rhoList.get(rhoList.size() - 1);
		double factor = 0.75 / Math.PI; // conversion: rho = factor * phi0

		// ===========================
		// Step 1: build fine grid arrays
		// Compute f_id, f_FR, f_ex at each of the 500 grid points
		// ===========================
		double[] phi0Arr = new double[nFitPoints];
		double[] phiArr  = new double[nFitPoints];
		double[] rhoArr  = new double[nFitPoints];
		double[] ZvirArr = new double[nFitPoints];
		double[] FidArr  = new double[nFitPoints]; // ideal gas free energy density at each point
		double[] FfrArr  = new double[nFitPoints]; // Flory-Rehner free energy density at each point
		double[] FexArr  = new double[nFitPoints]; // excess free energy density at each point
		double[] FtotArr = new double[nFitPoints]; // total free energy density at each point

		for (int i = 0; i < nFitPoints; i++) {

			// uniformly spaced density grid from rhoMin to rhoMax
			double rho   = rhoMin + (rhoMax - rhoMin) * i / (nFitPoints - 1);
			double phi0  = rho / factor;                      // dry volume fraction
			double alpha = interpolateAlpha(rho);             // swelling ratio at this rho
			double phi   = interpolatePhi(rho);               // swollen volume fraction at this rho

			// Flory-Rehner free energy density interpolated from simulation data
			// referenced to the reservoir state (single particle in solvent)
			double FfrV = interpolate(rhoList, floryFperVolList, rho);

			// B2 computed analytically from the interpolated swelling ratio alpha
			double B2 = particles.secondVirialCoefficient(alpha, alpha);

			// Z_virial smooth curve: Z = 1 + B2*rho + B3*rho^2 + ...
			double Zvir = 1.0 + B2 * rho;
			for (int k = 0; k < virial.length; k++) {
				int n = k + 3; // n=3 for B3, n=4 for B4, etc.
				Zvir += virial[k] * Math.pow(rho, n - 1);
			}

			// Excess free energy density from integrating (Z-1)/rho over rho:
			// f_ex/V = B2*rho^2 + B3*rho^3/2 + B4*rho^4/3 + ...
			double FexV = B2 * rho * rho;
			for (int k = 0; k < virial.length; k++) {
				int n = k + 3;
				FexV += virial[k] * Math.pow(rho, n) / (n - 1.0);
			}

			// Ideal gas free energy density:
			// f_id/V = (3/4pi) * phi0 * (ln(rho) - 1) + Stirling correction
			double stirling = factor * phi0
							* Math.log(2.0 * Math.PI * particles.N) / (2.0 * particles.N);
			
			double FidV = factor * phi0 * (Math.log(rho) - 1.0) + stirling;

			// Total free energy density: f_total = f_id + f_FR + f_ex
			double FtotV = FidV + FfrV + FexV;

			// Store in arrays for polynomial fitting in Step 2
			phi0Arr[i] = phi0;
			phiArr[i]  = phi;
			rhoArr[i]  = rho;
			ZvirArr[i] = Zvir;
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

		// ===========================
		// Step 2: fit separate polynomials to each free energy contribution
		// Using maxPower as the polynomial degree for consistency with the virial fit
		// ===========================
		int polyDegree = maxPower; // same order as virial EOS fit for consistency

		// Fit one polynomial per contribution
		// coeffs = [c0, c1, c2, ..., cn] where f(phi0) = c0 + c1*phi0 + c2*phi0^2 + ...
		double[] coeffsId  = fitPoly(FidArr,  phi0Arr, nFitPoints, polyDegree); // ideal
		double[] coeffsFR  = fitPoly(FfrArr,  phi0Arr, nFitPoints, polyDegree); // Flory-Rehner
		double[] coeffsEx  = fitPoly(FexArr,  phi0Arr, nFitPoints, polyDegree); // excess
		double[] coeffsTot = fitPoly(FtotArr, phi0Arr, nFitPoints, polyDegree); // total

		// ===========================
		// Step 3: evaluate mu and PV/NkT analytically for each contribution
		//
		// Chemical potential:   mu = (4pi/3) * dF/dphi0
		// Pressure:    PV/NkT = (4pi/3)/phi0 * (phi0 * dF/dphi0 - F_fit)
		//
		// Both computed from the analytical derivative of the polynomial —
		// no numerical differentiation error
		// ===========================
		double prefac = 4.0 * Math.PI / 3.0; // prefactor = 4pi/3 = 1/factor

		for (int i = 0; i < nFitPoints; i++) {
			double p = phi0Arr[i]; // current phi0 value

			// Evaluate each free energy contribution from its polynomial fit
			double F_id  = evalPoly(coeffsId,  p, polyDegree); // f_id at phi0
			double F_fr  = evalPoly(coeffsFR,  p, polyDegree); // f_FR at phi0
			double F_ex  = evalPoly(coeffsEx,  p, polyDegree); // f_ex at phi0
			double F_tot = evalPoly(coeffsTot, p, polyDegree); // f_total at phi0

			// Evaluate analytical derivatives dF/dphi0 for each contribution
			double dFid  = evalPolyDeriv(coeffsId,  p, polyDegree); // df_id/dphi0
			double dFfr  = evalPolyDeriv(coeffsFR,  p, polyDegree); // df_FR/dphi0
			double dFex  = evalPolyDeriv(coeffsEx,  p, polyDegree); // df_ex/dphi0
			double dFtot = evalPolyDeriv(coeffsTot, p, polyDegree); // df_total/dphi0

			// Chemical potential contributions: mu = (4pi/3) * dF/dphi0
			double mu_id  = prefac * dFid;  // ideal gas contribution to mu
			double mu_fr  = prefac * dFfr;  // Flory-Rehner contribution to mu
			double mu_ex  = prefac * dFex;  // excess contribution to mu
			double mu_tot = prefac * dFtot; // total chemical potential

			// Pressure contributions: PV/NkT = (4pi/3)/phi0 * (phi0*dF/dphi0 - F_fit)
			double PV_id  = prefac / p * (p * dFid  - F_id);  // ideal contribution to PV/NkT
			double PV_fr  = prefac / p * (p * dFfr  - F_fr);  // FR contribution to PV/NkT
			double PV_ex  = prefac / p * (p * dFex  - F_ex);  // excess contribution to PV/NkT
			double PV_tot = prefac / p * (p * dFtot - F_tot); // total PV/NkT

			// Store free energy contributions
			FidFitCurveList.add(F_id);
			FfrFitCurveList.add(F_fr);
			FexFitCurveList2.add(F_ex);
			FtotFitCurveList.add(F_tot);

			// Store chemical potential contributions
			muIdFitCurveList.add(mu_id);
			muFRFitCurveList.add(mu_fr);
			muExFitCurveList2.add(mu_ex);
			muExFitCurveList.add(mu_tot);     // existing list reused for total mu

			// Store pressure contributions
			PVIdFitCurveList.add(PV_id);
			PVFRFitCurveList.add(PV_fr);
			PVExFitCurveList.add(PV_ex);
			pressureFitCurveList.add(PV_tot); // existing list reused for total PV/NkT
		}
	}

	private double interpolate(List<Double> xList, List<Double> yList, double x) {

		/**
		 * Performs linear interpolation to estimate y(x) from discrete data points.
		 *
		 * Given a set of ordered pairs (x_i, y_i), this function finds the interval
		 * [x0, x1] such that x0 ≤ x ≤ x1, and estimates y using:
		 *
		 *     y(x) = y0 + ((x - x0) / (x1 - x0)) * (y1 - y0)
		 *
		 */

		// Loop through consecutive pairs (x0, x1)
		for (int i = 0; i < xList.size() - 1; i++) {

			double x0 = xList.get(i);       // lower bound of interval
			double x1 = xList.get(i + 1);   // upper bound of interval

			// Check if x lies within this interval
			if (x >= x0 && x <= x1) {

				// Compute interpolation weight (fraction between x0 and x1)
				double t = (x - x0) / (x1 - x0);

				// Linearly interpolate y between y0 and y1
				return yList.get(i) + t * (yList.get(i + 1) - yList.get(i));
			}
		}

		// Fallback: if x is outside the range, return the last value
		return yList.get(yList.size() - 1);
	}

	private double interpolateZfit(double rho) {

		/**
		 * Performs linear interpolation on the fitted Z(ρ) curve to estimate
		 * the compressibility factor Z at an arbitrary density rho.
		 *
		 * Uses precomputed lists:
		 *   rhoFitCurveList → sampled density values (sorted)
		 *   zFitCurveList   → corresponding Z values from virial fit
		 *
		 * For a given rho, the method finds the interval [r1, r2] such that:
		 *      r1 ≤ rho ≤ r2
		 *
		 * and computes:
		 *
		 *      Z(rho) = z1 + ((rho - r1) / (r2 - r1)) * (z2 - z1)
		 *
		 * This ensures a smooth interpolation of the virial EOS curve.
		 */

		// --- Lower boundary ---
		// If rho is below the minimum fitted density, return the first Z value
		if (rho <= rhoFitCurveList.get(0)) {
			return zFitCurveList.get(0);
		}

		// --- Upper boundary ---
		// If rho is above the maximum fitted density, return the last Z value
		if (rho >= rhoFitCurveList.get(rhoFitCurveList.size() - 1)) {
			return zFitCurveList.get(zFitCurveList.size() - 1);
		}

		// --- Find the interval containing rho ---
		for (int i = 0; i < rhoFitCurveList.size() - 1; i++) {

			double r1 = rhoFitCurveList.get(i);       // lower density bound
			double r2 = rhoFitCurveList.get(i + 1);   // upper density bound

			// Check if rho lies between r1 and r2
			if (rho >= r1 && rho <= r2) {

				double z1 = zFitCurveList.get(i);     // Z at r1
				double z2 = zFitCurveList.get(i + 1); // Z at r2

				// Compute interpolation fraction
				double t = (rho - r1) / (r2 - r1);

				// Linear interpolation of Z
				return z1 + t * (z2 - z1);
			}
		}

		// --- Fallback ---
		// Should not normally be reached if rho is within bounds
		return 0.0;
	}

	private double interpolateAlpha(double rho) {

		/**
		 * Performs linear interpolation to estimate the swelling ratio α (alpha)
		 * at an arbitrary density rho using discrete simulation data.
		 *
		 * Uses:
		 *   rhoList             → sampled densities (sorted)
		 *   swellingRatioList   → corresponding swelling ratios α(ρ)
		 *
		 * For a given rho, the method finds the interval [rho1, rho2] such that:
		 *      rho1 ≤ rho ≤ rho2
		 *
		 * and computes:
		 *
		 *      α(rho) = a1 + ((rho - rho1) / (rho2 - rho1)) * (a2 - a1)
		 *
		 * Assumptions:
		 *   - rhoList is sorted in ascending order
		 *   - swelling ratio varies smoothly with density
		 */

		int n = rhoList.size();  // number of available data points

		// --- Lower boundary ---
		// If rho is below the minimum sampled density, return first alpha value
		if (rho <= rhoList.get(0)) {
			return swellingRatioList.get(0);
		}

		// --- Upper boundary ---
		// If rho is above the maximum sampled density, return last alpha value
		if (rho >= rhoList.get(n - 1)) {
			return swellingRatioList.get(n - 1);
		}

		// --- Find the interval containing rho ---
		for (int i = 0; i < n - 1; i++) {

			double rho1 = rhoList.get(i);       // lower density bound
			double rho2 = rhoList.get(i + 1);   // upper density bound

			// Check if rho lies within this interval
			if (rho >= rho1 && rho <= rho2) {

				double a1 = swellingRatioList.get(i);     // alpha at rho1
				double a2 = swellingRatioList.get(i + 1); // alpha at rho2

				// Linear interpolation of alpha
				return a1 + (rho - rho1) * (a2 - a1) / (rho2 - rho1);
			}
		}

		// --- Fallback ---
		// Return last value as a safety measure (should rarely be reached)
		return swellingRatioList.get(n - 1);
	}

	private double interpolatePhi(double rho) {

		/**
		 * Performs linear interpolation to estimate the swollen volume fraction φ
		 * at an arbitrary density rho using discrete simulation data.
		 *
		 * Uses:
		 *   rhoList              → sampled densities (sorted)
		 *   volumefractionList  → corresponding swollen volume fractions φ(ρ)
		 *
		 * For a given rho, the method finds the interval [rho1, rho2] such that:
		 *      rho1 ≤ rho ≤ rho2
		 *
		 * and computes:
		 *
		 *      φ(rho) = φ1 + ((rho - rho1) / (rho2 - rho1)) * (φ2 - φ1)
		 *
		 * Assumptions:
		 *   - rhoList is sorted in ascending order
		 *   - φ varies smoothly with density
		 */

		int n = rhoList.size();  // number of available data points

		// --- Lower boundary ---
		// If rho is below the minimum sampled density, return the first φ value
		if (rho <= rhoList.get(0)) {
			return volumefractionList.get(0);
		}

		// --- Upper boundary ---
		// If rho is above the maximum sampled density, return the last φ value
		if (rho >= rhoList.get(n - 1)) {
			return volumefractionList.get(n - 1);
		}

		// --- Find the interval containing rho ---
		for (int i = 0; i < n - 1; i++) {

			double rho1 = rhoList.get(i);       // lower density bound
			double rho2 = rhoList.get(i + 1);   // upper density bound

			// Check if rho lies within this interval
			if (rho >= rho1 && rho <= rho2) {

				double phi1 = volumefractionList.get(i);     // φ at rho1
				double phi2 = volumefractionList.get(i + 1); // φ at rho2

				// Linear interpolation of φ
				return phi1 + (rho - rho1) * (phi2 - phi1) / (rho2 - rho1);
			}
		}

		// --- Fallback ---
		// Return last value as a safety measure (should rarely be reached)
		return volumefractionList.get(n - 1);
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

				// Compute B2 
				double B2;
					
				// B2 depends on alpha (keep per point)
				B2 = particles.secondVirialCoefficient(alpha, alpha);
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
				double stirlingApprox = (0.75 / Math.PI) * dryVolFrac
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
				// Skip first and last points for chemical potential calculation (central difference)

				double rho = rhoList.get(i); // density at state i
				double B2 = secondVirialCoefficientList.get(i); // exact B2 at state i 

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
			// Chemical potential from fitted virial EOS
			// ideal + excess analytic, FR from derivative of F_FR/V
			// ===========================
			for (int i = 1; i < nStates - 1; i++) {

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
				double rhoNext = rhoList.get(i + 1);
				double rhoPrev = rhoList.get(i - 1);

				double fFRNext = floryFperVolList.get(i + 1);
				double fFRPrev = floryFperVolList.get(i - 1);

				double muFR = (fFRNext - fFRPrev) / (rhoNext - rhoPrev);

				double muTotal = muIdeal + muEx + muFR;

				chemicalPotList.set(i, muTotal);
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
		control.setValue("DryVolFrac increment", 0.0001);    // 0.0002
		control.setValue("DryVolFrac Max", 0.003);     // 0.0025
		control.setValue("Initial configuration", "FCC");
		//control.setValue("Initial configuration", "random-FCC");
		control.setValue("N", 108); // number of particles
		//control.setValue("N", 500); for FCC lattice, N/4 should be a perfect cube
        control.setValue("Dry radius [nm]", 50);
        control.setValue("x-link fraction", 0.00003); // 0.001
        // control.setValue("Dry volume fraction", 0.01);
        control.setValue("Young's calibration", 1.0); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
        control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.05);
		control.setValue("Max virial power", 1); // number of virial coefficients to fit (B3 to B_{maxPower+2})
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
			File outputFile = new File("data/APS_2026/Liquid_Phase/Facet/Hertzian_Spheres/FacetLiquid_maxPower" + maxPower + "_" + particles.fileExtension + ".txt");
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
			bw1.write("phi0, phi, rho, mu/kT, Z_measured, Z_virial_Fit, <F_total>/V, QMelting, <F_id>/V, <F_ex>/V, <F_FR>/V, <alpha>, Zeta, B2_rho, B2_phi0, B2_HS, B2*");
			bw1.newLine();

			double factor = 0.75 / Math.PI; // Conversion factor

			for (int i = 1; i < dryVolFracs.size() - 1; i++) {

				double roundedDryVolFrac = Double.parseDouble(decimalFormat.format(dryVolFracs.get(i)));

				double Z_measured = meanPressures.get(i);     // Measured from simulation (virial theorem)
				double rho = rhoList.get(i);
				double Z_virial = varialPressuresList.get(i); // Calculated from virial expansion

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
			File virialFile = new File("data/APS_2026/Liquid_Phase/Facet/Hertzian_Spheres/VirialCoefficientsMaxPower" + maxPower + "_" + particles.fileExtension + ".txt");			
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

		try {
			File fitCurveFile = new File("data/APS_2026/Liquid_Phase/Facet/Hertzian_Spheres/VirialFitCurve_maxPower"
					+ maxPower + "_" + particles.fileExtension + ".txt");

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
				FidFitCurveList.get(i)     + ", " + // ideal free energy density
				FfrFitCurveList.get(i)     + ", " + // Flory-Rehner free energy density
				FexFitCurveList2.get(i)    + ", " + // excess free energy density
				FtotFitCurveList.get(i)    + ", " + // total free energy density
				muIdFitCurveList.get(i)    + ", " + // ideal contribution to mu
				muFRFitCurveList.get(i)    + ", " + // FR contribution to mu
				muExFitCurveList2.get(i)   + ", " + // excess contribution to mu
				muExFitCurveList.get(i)    + ", " + // total mu
				PVIdFitCurveList.get(i)    + ", " + // ideal contribution to PV/NkT
				PVFRFitCurveList.get(i)    + ", " + // FR contribution to PV/NkT
				PVExFitCurveList.get(i)    + ", " + // excess contribution to PV/NkT
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
			SimulationControl control = SimulationControl.createApp(new HertzSpheresLiquidVirialCoefficient_V2App());
	}
}
