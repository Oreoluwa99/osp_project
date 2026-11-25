/*
 * ConfinedMicrogelsStabilityApp.java
 *
 * Runs a Monte Carlo simulation of compressible microgels between
 * two fixed, flat walls separated by Lx. The walls stay at x = 0
 * and x = Lx throughout the run.
 *
 * Particle centers move in x between the walls. The y and z
 * directions are periodic; the x direction is not.
 *
 * The particles start from an FCC lattice, with its layers spread
 * evenly across the gap between the walls. At each dry volume
 * fraction, the app rebuilds the lattice, equilibrates the particles,
 * and computes the in-plane Lindemann ratio from the mean squared
 * displacement in y and z. The dry volume fraction is then increased
 * and the run repeated, so the crystal's stability is mapped at
 * fixed Lx.
 *
 * Each run gives the Lindemann ratio, volume fraction, and mean
 * swelling ratios as functions of the dry volume fraction. Repeating
 * the scan over several wall separations maps the stability diagram.
 *
 * Authors: Oreoluwa Alade and Alan Denton
 */

package org.opensourcephysics.sip.Hertz;

import java.io.File;
import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;
import java.text.DecimalFormat;
import java.util.List;
import java.util.ResourceBundle.Control;
import java.util.ArrayList;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.display3d.simple3d.ElementBox;

public class ConfinedMicrogelsStabilityApp extends AbstractSimulation {
	public enum WriteModes {WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
	ConfinedMicrogelsStability particles = new ConfinedMicrogelsStability();
	// Counts sampled particle-center positions between the two walls
	double[] positionDist;
	boolean firstDensityRun = true;
	StringBuilder sizeOutput;
	StringBuilder densityOutput;
	StringBuilder ratioOutput;
	double alphaSum;
	double modifiedAlphaSum;
	StringBuilder swellingOutput;

	PlotFrame energyData = new PlotFrame("MC steps", "<E_pair>/N", "Mean pair energy per particle");
	PlotFrame pressureData = new PlotFrame("MC steps", "PV/NkT", "Mean pressure");
	PlotFrame sizeData = new PlotFrame("MC steps", "alpha", "Mean swelling ratio");
	//PlotFrame capToGelVolData = new PlotFrame("MC steps", "Vc/Vm", "Cap volume to Microgel Volume");
	
	Display3DFrame display3d = new Display3DFrame("Simulation animation");
	
	ElementSphere[] nanoSphere;
	ElementBox rightWall;
	boolean added = false;

	double dryVolFracStart, dryVolFracMax, dryVolFrac;
	double lambda;
	boolean structure;
	ConfinedRDF rdf;
	SSF ssf;
	boolean incrementDryVolFrac = true;
	boolean scanComplete = false;

	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values to 3 dp

	/**
	 * Initializes the model.
	 */
	public void initialize() {

		scanComplete = false;
		alphaSum = 0.0;
		modifiedAlphaSum = 0.0;
		added = false;

		//particles.dlambda = control.getDouble("Lambda increment");// the coupling constant increment
		dryVolFracStart = control.getDouble("DryVolFracStart");
		dryVolFracMax = control.getDouble("DryVolFrac Max");
		particles.dphi = control.getDouble("DryVolFrac increment");
		particles.N = control.getInt("N"); // number of particles
		String configuration = control.getString("Initial configuration");
		particles.initConfig = configuration;

		particles.nLayers = control.getInt("Layers"); // number of layers between the walls

		if (configuration.toUpperCase().equals("LAYERED")) {
			particles.nCols = 7; // particles in each row
			particles.nRows = 8; // rows in each layer; 8/7 is close to 2/sqrt(3), so a triangular layer fits the square box
			particles.N = particles.nLayers * particles.nCols * particles.nRows; // 56 particles per layer
		}

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
		particles.Lx = control.getDouble("Wall separation Lx");
		particles.wallYc = particles.Young;
		structure = control.getBoolean("Calculate structure");
		if (incrementDryVolFrac) {
			dryVolFrac = dryVolFracStart;
			particles.dryVolFrac = dryVolFrac;
			incrementDryVolFrac = false;
		}

		particles.initialize(configuration);
		control.println("Lx = " + particles.Lx + ", Yc = " + particles.Young
				+ ", N = " + particles.N + ", side = " + decimalFormat.format(particles.side));
		positionDist = new double[100]; // Start a new position count for this run; all 100 entries begin at zero
		
		// write out system parameters
		System.out.println("nMon: " + particles.nMon);
		System.out.println("nChains: " + particles.nChains);
		System.out.println("reservoirSwellingRatio: " + particles.reservoirSR);
		System.out.println("reservoirVolFrac: " + particles.reservoirVolFrac);

		if (display3d != null)
		display3d.dispose(); // closes old simulation frame if present
		display3d = new Display3DFrame("Simulation animation");
		display3d.setPreferredMinMax(0, particles.Lx, 0, particles.side, 0, particles.side);
		display3d.setSquareAspect(true);

		double wallThickness = 0.08;

		// Left wall at x = 0
		ElementBox leftWall = new ElementBox(); // describes each wall
		leftWall.setSizeXYZ(wallThickness, particles.side, particles.side); // sets the size of the wall. It is particles.side wide and tall with wallThickness 
		leftWall.setXYZ(0.0, particles.side / 2.0, particles.side / 2.0); // sets the position of the wall
		leftWall.getStyle().setFillColor(new Color(0, 0, 0, 30)); // gives the wall a blue fill
		leftWall.getStyle().setLineColor(Color.BLACK); // makes the outline blue
		display3d.addElement(leftWall); //makes it appear in the animation

		// Right wall at x = Lx
		rightWall = new ElementBox(); // describes each wall
		rightWall.setSizeXYZ(wallThickness, particles.side, particles.side); // sets the size of the wall. It is particles.side wide and tall with wallThickness
		rightWall.setXYZ(particles.Lx, particles.side / 2.0, particles.side / 2.0); // sets the position of the wall
		rightWall.getStyle().setFillColor(new Color(0, 0, 0, 30)); // gives the wall a blue fill
		rightWall.getStyle().setLineColor(Color.BLACK); // makes the outline blue
		display3d.addElement(rightWall); //makes it appear in the animation

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
			// nanoSphere[i].getStyle().setFillColor(Color.RED);
			nanoSphere[i].getStyle().setFillColor(new Color(220, 40, 40, 80));
			// Show y and z inside the periodic display box.
			double displayY = ((particles.y[i] % particles.side) + particles.side) % particles.side; //  Calculate the y position within the display box
			double displayZ = ((particles.z[i] % particles.side) + particles.side) % particles.side; //  Calculate the z position within the display box

			nanoSphere[i].setXYZ(particles.x[i], displayY, displayZ);
		}

		if (structure) {
			rdf = new ConfinedRDF(particles.x, particles.y, particles.z, particles.side, particles.Lx, particles.grBinWidth, particles.fileExtension,particles.Young);
			ssf = new SSF(particles.x, particles.y, particles.z, particles.side, particles.d, particles.deltaK, particles.fileExtension);
		}
	}

	/**
	 * Does a simulation step.
	 */
	public void doStep() { // logical step in the HertzSpheres class
		if (scanComplete) return; // stop if scan is complete

		particles.step();

		// if initial configuration is random, no particle interactions for delay/10
		if (particles.steps > particles.delay / 10.) {
			particles.scale = 1.;
		}

		if (particles.steps <= particles.stop) {
			if (particles.steps > particles.delay) {
				if ((particles.steps - particles.delay) % particles.snapshotInterval == 0) { // decide when to take a snapshot
					particles.sizeDistribution(); // counts the sizes of all 32 particles

					// Record the current accepted swelling ratios at this snapshot.
					for (int i = 0; i < particles.N; i++) {
						double alpha = particles.a[i];
						double sphereVolume = (4.0 / 3.0) * Math.PI
								* alpha * alpha * alpha;

						double remainingFraction = 1.0
								- particles.capVolSumBeforeMoveArray[i] / sphereVolume;

						if (remainingFraction <= 0.0) {
							throw new IllegalStateException(
									"Invalid cap volume for particle " + i);
						}

						double alphaTilde = alpha * Math.cbrt(remainingFraction);

						alphaSum += alpha;
						modifiedAlphaSum += alphaTilde;
					}


					// Count each particle's x position at this sampling time
					double dx = particles.Lx / positionDist.length;

					for (int i = 0; i < particles.N; i++) {
						// Find the entry corresponding to this particle's x position
						int b = (int) (particles.x[i] / dx);

						// Include a particle whose center is exactly at the right wall
						if (b == positionDist.length) b--;

						positionDist[b]++; // Record one particle at this x position
					}

					if (structure) { // accumulate statistics for structural properties
						rdf.update(); // g(r)
						ssf.update(); // S(k)
					}

					particles.calculateVolumeFraction();
				}
			}
		}

		if (particles.steps == particles.stop) {

			// Add this dry volume fraction's results.
			// At the final dry volume fraction, write the complete scan files.
			writeData();

			System.out.println("State point finished: Lx=" + particles.Lx
				+ ", phi0=" + dryVolFrac
				+ ", samples=" + particles.numberOfConfigurations);

			// Move to the next dry volume fraction and rebuild the lattice.
			dryVolFrac += particles.dphi;

			if (dryVolFrac > dryVolFracMax) {
				control.println("Scan finished at Lx = " + particles.Lx + ", files written");
				return;
			}

			control.println("Starting phi0 = " + decimalFormat.format(dryVolFrac));
			particles.dryVolFrac = dryVolFrac;
			initialize();
			return;
		}
		
		// plot mean energy, pressure, swelling ratio, total energy
		energyData.append(0, particles.steps, particles.meanPairEnergy());
		pressureData.append(1, particles.steps, particles.meanPressure());
		
		double meanAlphaNow = 0.0;
		for (int i = 0; i < particles.N; i++) {
			meanAlphaNow += particles.a[i];
		}
		meanAlphaNow /= particles.N;

		sizeData.append(2, particles.steps, meanAlphaNow);

		//capToGelVolData.append(0, particles.steps, particles.meanCapToGelVol());

		display3d.setMessage("Number of steps: " + particles.steps); // update steps

		if (control.getBoolean("Visualization on")) { // visualization updates
			for (int i = 0; i < particles.N; i++) {
				nanoSphere[i].setSizeXYZ(particles.a[i] * 2, particles.a[i] * 2, particles.a[i] * 2);
				// nanoSphere[i].getStyle().setFillColor(Color.RED);
				nanoSphere[i].getStyle().setFillColor(new Color(220, 40, 40, 80));
				// Show y and z inside the periodic display box.
				double displayY = ((particles.y[i] % particles.side) + particles.side) % particles.side;
				double displayZ = ((particles.z[i] % particles.side) + particles.side) % particles.side;

				nanoSphere[i].setXYZ(particles.x[i], displayY, displayZ);
			}
		}

	}

	/**
	 * Resets the model to its default state.
	 */
	public void reset() {

		incrementDryVolFrac = true;
		firstDensityRun = true;
		enableStepsPerDisplay(true);
		//control.setValue("Lambda increment", 0.1);
		control.setValue("DryVolFracStart", 0.005);
		control.setValue("DryVolFrac Max", 0.0205);
		control.setValue("DryVolFrac increment", 0.001);
		control.setValue("Initial configuration", "LAYERED");
		control.setValue("Layers", 5); // number of layers between the walls (used when the start is LAYERED)
		control.setValue("N", 280); // N/4 must be a perfect cube for FCC: 32, 108, 256
		control.setValue("Dry radius [nm]", 50);
		control.setValue("x-link fraction", 0.001);
		control.setValue("Young's calibration", 1); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.05);
		control.setValue("Wall separation Lx", 35.0);
		control.setValue("Delay", 10000); // steps after which statistics collection starts
		control.setValue("Snapshot interval", 50); // steps separating successive samples
		control.setValue("Stop", 100000); // steps after which statistics collection stops
		control.setValue("Size bin width", .001); // bin width of particle radius histogram
		control.setValue("g(r) bin width", 0.1); // bin width of g(r) histogram
		control.setValue("Delta k", .005); // bin width of S(k) histogram
		control.setValue("File extension", "prod");
		// control.setValue("Calculate structure", false); // true means calculate g(r) and S(k)
		control.setValue("Calculate structure", false); // true means calculate g(r) and S(k)
		control.setAdjustableValue("Visualization on", true);
	}

	public void stop() {
		control.println("Number of MC steps = " + particles.steps);
		control.println("xLinkFrac = "+particles.xLinkFrac);
		control.println("<E_pair>/N = " + decimalFormat.format(particles.meanPairEnergy()));
		control.println("<F>/N = " + decimalFormat.format(particles.meanFreeEnergy()));
		control.println("PV/NkT = " + decimalFormat.format(particles.meanPressure()));
	}

	public void writeData() {

		double M = particles.numberOfConfigurations;

		if (M == 0) {
			throw new IllegalStateException("No configurations were sampled.");
		}

		// phi0 is set for the run; phi is averaged over sampled configurations.
		double phi0 = particles.dryVolFrac;
		double phi = particles.meanVolFrac();

		double meanAlpha = alphaSum / (M * particles.N);
		double meanAlphaTilde = modifiedAlphaSum / (M * particles.N);

		double lindemann = Math.sqrt(particles.meanSquareDisplacement()) / particles.nnDistance; // in-plane Lindemann ratio

		double dx = particles.Lx / positionDist.length;
		double meanDensity = particles.N
				/ (particles.Lx * particles.side * particles.side);

		// Set up the four outputs when the first dry volume fraction finishes
		if (firstDensityRun) {
			swellingOutput = new StringBuilder();
			sizeOutput = new StringBuilder();
			densityOutput = new StringBuilder();
			ratioOutput = new StringBuilder();

			String[] parameters = {
				"# Wall separation Lx is listed in each data row",
				"# Ly = Lz changes with phi0",
				"# Number of particles: " + particles.N,
				"# Initial configuration: " + particles.initConfig,
				"# Layers: " + particles.nLayers,
				"# Dry microgel radius [nm]: " + particles.dryR,
				"# Number of monomers: " + particles.nMon,
				"# Number of chains: " + particles.nChains,
				"# Flory chi: " + particles.chi,
				"# Young's calibration: " + particles.Young,
				// the wall details
				"# Wall separation: " + particles.Lx,
				"# Dry volume fraction start: " + dryVolFracStart,
				"# Dry volume fraction max: " + dryVolFracMax,
				"# Dry volume fraction increment: " + particles.dphi,
				"# Wall Hertz calibration: " + particles.wallYc,
				"# Cross-link fraction: " + particles.xLinkFrac,
				"# MC steps: " + particles.steps,
				"# Equilibration steps: " + particles.delay,
				"# Snapshot interval: " + particles.snapshotInterval,
				"# Displacement tolerance: " + particles.tolerance,
				"# Radius change tolerance: " + particles.atolerance,
				"# Size bin width: " + particles.sizeBinWidth,
				"# Position bins: " + positionDist.length
			};

			for (String parameter : parameters) {
				swellingOutput.append(parameter).append('\n');
				sizeOutput.append(parameter).append('\n');
				densityOutput.append(parameter).append('\n');
				ratioOutput.append(parameter).append('\n');
			}

			swellingOutput.append(
					"\n# Lx, phi0, phi, meanAlpha, meanAlphaTilde, L\n");
			sizeOutput.append("\n# Lx, phi0, phi, alpha, P(alpha)\n");
			densityOutput.append("\n# Lx, phi0, phi, x, n(x)\n");
			ratioOutput.append(
					"\n# Lx, phi0, phi, x, n(x)/meanDensity\n");

			firstDensityRun = false;
		}

		// One row of mean swelling ratios for this dry volume fraction
		swellingOutput.append(particles.Lx).append(", ")
				.append(phi0).append(", ")
				.append(phi).append(", ")
				.append(meanAlpha).append(", ")
				.append(meanAlphaTilde).append(", ")
				.append(lindemann).append('\n');

		control.println("Lx = " + particles.Lx + ", Yc = " + particles.Young
				+ ", phi0 = " + decimalFormat.format(phi0)
				+ ", phi = " + decimalFormat.format(phi)
				+ ", alpha = " + decimalFormat.format(meanAlpha)
				+ ", L = " + decimalFormat.format(lindemann));

		// Add this state point's swelling-ratio distribution to memory.
		for (int b = 0; b < particles.numberBins; b++) {
			double alpha = (b + 0.5) * particles.sizeBinWidth;
			double probability = particles.sizeDist[b]
					/ (M * particles.N * particles.sizeBinWidth);

			sizeOutput.append(particles.Lx).append(", ")
					.append(phi0).append(", ")
					.append(phi).append(", ")
					.append(alpha).append(", ")
					.append(probability).append('\n');
		}

		// Add this state point's number-density profile and density ratio to memory.
		for (int b = 0; b < positionDist.length; b++) {
			double x = (b + 0.5) * dx;

			double numberDensity = positionDist[b]
					/ (M * particles.side * particles.side * dx);

			double densityRatio = numberDensity / meanDensity;

			densityOutput.append(particles.Lx).append(", ")
					.append(phi0).append(", ")
					.append(phi).append(", ")
					.append(x).append(", ")
					.append(numberDensity).append('\n');

			ratioOutput.append(particles.Lx).append(", ")
					.append(phi0).append(", ")
					.append(phi).append(", ")
					.append(x).append(", ")
					.append(densityRatio).append('\n');
		}

		// Keep results in memory until the final dry volume fraction finishes.
		if (dryVolFrac + particles.dphi <= dryVolFracMax) {
			return;
		}

		// The full scan is finished: create the folder and four files.
		try {
			File directory = new File(
				"data/confined/stability/Yc" + particles.Young + "/N" + particles.N
			);

			if (!directory.exists() && !directory.mkdirs()) {
				throw new IOException("Could not create " + directory);
			}

			String scanName = "N" + particles.N
				+ "_Yc" + particles.wallYc
				+ "_Lx" + particles.Lx
				+ "_phi0" + dryVolFracStart
				+ "_to_" + dryVolFrac
				+ "_dphi" + particles.dphi
				+ "_run" + particles.fileExtension + ".txt";

			File swellingFile = new File(directory,
					"swelling_means_" + scanName);
			File sizeFile = new File(directory,
					"size_distribution_" + scanName);
			File densityFile = new File(directory,
					"number_density_" + scanName);
			File ratioFile = new File(directory,
					"density_ratio_" + scanName);

			// Protect files from an earlier scan.
			if (swellingFile.exists()
					|| sizeFile.exists()
					|| densityFile.exists()
					|| ratioFile.exists()) {
				throw new IllegalStateException(
						"Scan files already exist in " + directory
						+ ". Choose a new file extension before running.");
			}

			try (BufferedWriter swellingWriter =
						new BufferedWriter(new FileWriter(swellingFile));
				BufferedWriter sizeWriter =
						new BufferedWriter(new FileWriter(sizeFile));
				BufferedWriter densityWriter =
						new BufferedWriter(new FileWriter(densityFile));
				BufferedWriter ratioWriter =
						new BufferedWriter(new FileWriter(ratioFile))) {

				swellingWriter.write(swellingOutput.toString());
				sizeWriter.write(sizeOutput.toString());
				densityWriter.write(densityOutput.toString());
				ratioWriter.write(ratioOutput.toString());
			}

		} catch (IOException e) {
			throw new RuntimeException(
					"Could not write stability data.", e);
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
			SimulationControl control = SimulationControl.createApp(new ConfinedMicrogelsStabilityApp());
	}
}
