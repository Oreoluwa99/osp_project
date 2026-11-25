/*
 * ConfinedHardSpheresApp.java
 *
 * Runs a Monte Carlo simulation of fixed-size hard spheres between
 * two parallel hard walls for comparison with confined hard-sphere
 * results.
 *
 * The left wall stays at x = 0. The right wall starts at x = Lx and
 * can move inward between simulation stages. Spheres cannot overlap
 * each other or cross either wall. The y and z directions are
 * periodic; the x direction is not.
 *
 * At each wall separation, the app equilibrates the spheres and
 * records the number-density profile n(x) and the density ratio
 * n(x) / meanDensity. The number of spheres and lateral box
 * dimensions remain fixed during compression.
 *
 * A wall decrement of zero runs one fixed-wall stage.
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

public class ConfinedHardSpheresApp extends AbstractSimulation {
	public enum WriteModes {WRITE_NONE, WRITE_RADIAL, WRITE_ALL;};
	ConfinedHardSpheres particles = new ConfinedHardSpheres();
	// Counts sampled particle-center positions between the two walls
	double[] positionDist;
	boolean firstDensityRun = true;
	StringBuilder sizeOutput;
	StringBuilder densityOutput;
	StringBuilder ratioOutput;

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
	RDF rdf;
	SSF ssf;
	boolean incrementDryVolFrac = true;
	public double wallDecrement;
	public double LxMin;
	boolean scanComplete = false;

	DecimalFormat decimalFormat = new DecimalFormat("#.#######"); // to round my dryVolFrac values to 3 dp

	/**
	 * Initializes the model.
	 */
	public void initialize() {

		scanComplete = false;
		added = false;

		//particles.dlambda = control.getDouble("Lambda increment");// the coupling constant increment
		dryVolFracStart = control.getDouble("DryVolFracStart");
		dryVolFracMax = control.getDouble("DryVolFrac Max");
		particles.dphi = control.getDouble("DryVolFrac increment");
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
		particles.Lx = control.getDouble("Wall separation Lx");
		wallDecrement = control.getDouble("Wall decrement dLx");
		LxMin = control.getDouble("Minimum wall separation Lx");
		particles.wallYc = control.getDouble("Wall Hertz calibration");
		structure = control.getBoolean("Calculate structure");
		if (incrementDryVolFrac) {
			dryVolFrac = dryVolFracStart;
			particles.dryVolFrac = dryVolFrac;
			incrementDryVolFrac = false;
		}

		particles.initialize(configuration);
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
		leftWall.getStyle().setFillColor(new Color(40, 100, 220, 20)); // gives the wall a blue fill
		leftWall.getStyle().setLineColor(Color.BLUE); // makes the outline blue
		display3d.addElement(leftWall); //makes it appear in the animation

		// Right wall at x = Lx
		rightWall = new ElementBox(); // describes each wall
		rightWall.setSizeXYZ(wallThickness, particles.side, particles.side); // sets the size of the wall. It is particles.side wide and tall with wallThickness
		rightWall.setXYZ(particles.Lx, particles.side / 2.0, particles.side / 2.0); // sets the position of the wall
		rightWall.getStyle().setFillColor(new Color(40, 100, 220, 20)); // gives the wall a blue fill
		rightWall.getStyle().setLineColor(Color.BLUE); // makes the outline blue
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
			nanoSphere[i].getStyle().setFillColor(Color.RED);
			// Show y and z inside the periodic display box.
			double displayY = ((particles.y[i] % particles.side) + particles.side) % particles.side;
			double displayZ = ((particles.z[i] % particles.side) + particles.side) % particles.side;

			nanoSphere[i].setXYZ(particles.x[i], displayY, displayZ);
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

			// Add this wall position's measurements to the output buffers.
			writeData();

			System.out.println("Stage finished: Lx=" + particles.Lx
				+ ", dLx=" + wallDecrement
				+ ", LxMin=" + LxMin
				+ ", samples=" + particles.numberOfConfigurations);

			// No wall movement is needed after the final stage.
			if (wallDecrement == 0.0 || particles.Lx <= LxMin) {
				scanComplete = true;
				stopSimulation(); // this stops the simulation automatically
				return; // exit before another particle move 
			}

			double oldLx = particles.Lx;
			double newLx = Math.max(LxMin, oldLx - wallDecrement);

			// Bring the same particle centers into the narrower gap.
			for (int i = 0; i < particles.N; i++) {
				particles.x[i] *= newLx / oldLx;
			}

			particles.Lx = newLx;
			particles.totalVol = newLx * particles.side * particles.side;
			particles.dryVolFrac =
					(4.0 * Math.PI * particles.N) / (3.0 * particles.totalVol);
			dryVolFrac = particles.dryVolFrac;

			// Clear measurements from the previous wall position.
			particles.steps = 0;
			particles.numberOfConfigurations = 0;
			particles.volFracAccumulator = 0;
			particles.energyAccumulator = 0;
			particles.pairEnergyAccumulator = 0;
			particles.freeEnergyAccumulator = 0;
			particles.virialAccumulator = 0;
			particles.squaredDisplacementAccumulator = 0;
			particles.springEnergyAccumulator = 0;

			for (int b = 0; b < particles.sizeDist.length; b++) {
				particles.sizeDist[b] = 0;
			}
			positionDist = new double[100];

			// Refresh stored pair energies for the new geometry.
			// At steps = 0, this does not record a snapshot.
			particles.calculateTotalEnergy(particles.lambda);

			rightWall.setXYZ(newLx, particles.side / 2.0,
							particles.side / 2.0);
			display3d.setPreferredMinMax(
					0, newLx, 0, particles.side, 0, particles.side);

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
				nanoSphere[i].getStyle().setFillColor(Color.RED);
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
		control.setValue("DryVolFracStart", 0.010);
		control.setValue("DryVolFrac Max", 0.014);
		control.setValue("DryVolFrac increment", 0.001);
		control.setValue("Initial configuration", "FCC");
		control.setValue("N", 32); // number of particles
		control.setValue("Dry radius [nm]", 50);
		control.setValue("x-link fraction", 0.001);
		control.setValue("Young's calibration", 1); // 10-1000
		control.setValue("chi", 0); // Flory interaction parameter
		control.setValue("Maximum radial distance", 10);
		control.setValue("Displacement tolerance", 0.1);
		control.setValue("Radius change tolerance", 0.05);
		control.setValue("Wall Hertz calibration", 1.0);
		control.setValue("Wall separation Lx", 20.0);
		control.setValue("Wall decrement dLx", 1.0);
		control.setValue("Minimum wall separation Lx", 18.0);
		control.setValue("Delay", 10000); // steps after which statistics collection starts
		control.setValue("Snapshot interval", 50); // steps separating successive samples
		control.setValue("Stop", 20000); // steps after which statistics collection stops
		control.setValue("Size bin width", .001); // bin width of particle radius histogram
		control.setValue("g(r) bin width", .002); // bin width of g(r) histogram
		control.setValue("Delta k", .005); // bin width of S(k) histogram
		control.setValue("File extension", "1");
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

		double dx = particles.Lx / positionDist.length;
		double meanDensity = particles.N
				/ (particles.Lx * particles.side * particles.side);

		// Set up the three outputs when the first density run finishes.
		if (firstDensityRun) {
			sizeOutput = new StringBuilder();
			densityOutput = new StringBuilder();
			ratioOutput = new StringBuilder();

			String[] parameters = {
				"# Wall separation Lx is listed in each data row",
				"# Ly = Lz: " + particles.side,
				"# Number of particles: " + particles.N,
				"# Sampled configurations per run: " + M,
				"# Initial configuration: " + particles.initConfig,
				"# Dry microgel radius [nm]: " + particles.dryR,
				"# Number of monomers: " + particles.nMon,
				"# Number of chains: " + particles.nChains,
				"# Flory chi: " + particles.chi,
				"# Young's calibration: " + particles.Young,
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
				sizeOutput.append(parameter).append('\n');
				densityOutput.append(parameter).append('\n');
				ratioOutput.append(parameter).append('\n');
			}

			sizeOutput.append("\n# Lx, phi0, phi, alpha, P(alpha)\n");
			densityOutput.append("\n# Lx, phi0, phi, x, n(x)\n");
			ratioOutput.append("\n# Lx, phi0, phi, x, n(x)/meanDensity\n");

			firstDensityRun = false;
		}

		// Add this run's swelling-ratio distribution to memory.
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

		// Add this run's number-density profile and density ratio to memory.
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

		// Keep the results in memory until the final wall position finishes
		if (wallDecrement > 0.0 && particles.Lx > LxMin) {
			return;
		}

		// The full scan is finished: now create the folder and three files.
		try {
			File directory = new File("data/confined/wall_compression_01");

			if (!directory.exists() && !directory.mkdirs()) {
				throw new IOException("Could not create " + directory);
			}

			String scanName = "Lx" + control.getDouble("Wall separation Lx")
				+ "_to_" + particles.Lx
				+ "_dLx" + wallDecrement
				+ "_" + particles.fileExtension + ".txt";

			File sizeFile = new File(directory,
					"size_distribution_" + scanName);
			File densityFile = new File(directory,
					"number_density_" + scanName);
			File ratioFile = new File(directory,
					"density_ratio_" + scanName);

			// Protect files from an earlier scan.
			if (sizeFile.exists()
					|| densityFile.exists()
					|| ratioFile.exists()) {
				throw new IllegalStateException(
						"Scan files already exist in " + directory
						+ ". Choose a new scan folder before running.");
			}

			try (BufferedWriter sizeWriter =
						new BufferedWriter(new FileWriter(sizeFile));
				BufferedWriter densityWriter =
						new BufferedWriter(new FileWriter(densityFile));
				BufferedWriter ratioWriter =
						new BufferedWriter(new FileWriter(ratioFile))) {

				sizeWriter.write(sizeOutput.toString());
				densityWriter.write(densityOutput.toString());
				ratioWriter.write(ratioOutput.toString());
			}

		} catch (IOException e) {
			throw new RuntimeException("Could not write density-scan data.", e);
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
			SimulationControl control = SimulationControl.createApp(new ConfinedHardSpheresApp());
	}
}
