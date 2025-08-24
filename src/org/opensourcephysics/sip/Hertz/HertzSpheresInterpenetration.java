package org.opensourcephysics.sip.Hertz;

import java.io.File;
import java.util.Random;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;

/*
 * HertzSpheresInterpenetration.java
 *
 * Purpose:
 *   Performs Monte Carlo simulations of compressible microgels using the Hertz interpenetration model and Flory-Rehner theory of polymer swelling.
 * 
 * Features:
 *   - Models full pairwise interactions using the Hertzian elastic potential and Flory–Rehner free energy.
 *   - Computes thermodynamic and structural properties:
 *       * Total and pair interaction energies
 *       * Flory–Rehner free energy contributions
 *       * Virial pressure and volume fraction
 *       * Spring (Einstein) energy and squared displacements
 *   - Tracks size distributions and swelling behavior of particles
 *   - Implements advanced interpenetration logic for overlapping microgels
 *   - Includes Einstein solid reference system for free energy estimation
 *
 * Authors: Alan Denton and Matt Urich
 * Edited by: Oreoluwa Alade
 *
 * Last Confirmed Working: September 2025
 * Framework: Open Source Physics (OSP)
 */

public class HertzSpheresInterpenetration {
   
      public int N; // number of particles
      public String initConfig; // initial configuration of particles
      public int nx; // number of columns and rows in initial crystal lattice
      public double side, totalVol; // side length and volume of cubic simulation box
      // Note: lengths are in units of dry microgel radius, energies in thermal (kT) units
      public double d; // distance between neighboring lattice sites
      public double x[], y[], z[], a[]; // coordinates (x, y, z) and radius (a) of particles
      public double x0[], y0[], z0[];//to store the initial positions
      public double energy[], energy1[]; // total energies of particles (Hertzian plus Flory-Rehner)
      public double pairEnergy[][], newPairEnergy[][]; // pair energies of particles
      public double totalEnergy, totalPairEnergy, totalFreeEnergy; // total energy of system
      public double totalVirial; // total virial of system (for pressure calculation)
      // various accumulators for computing thermodynamic properties
      public double energyAccumulator, pairEnergyAccumulator, freeEnergyAccumulator, virialAccumulator;
      public double tolerance, atolerance; // tolerances for trial p[article displacements, radius changes
      public double monRadius, nMon, nChains, chi, xLinkFrac; // microgel parameters
      public double B, Young; // prefactor of Hertz pair potential, Young's modulus calibration
      // if initial configuration is random, particles do not interact (scale=0) for delay/10 steps
      public double scale; // scale factor (0 or 1) for prefactor of Hertz potential
      public double dryR; // dry radius of particles
      public double reservoirSR; // reservoir swelling ratio (infinite dilution)
      // volume fractions of dry, swollen, and fully swollen (dilute) particles
      public double dryVolFrac, volFrac, reservoirVolFrac, dryVolFracStart, dphi;
      public int steps; // number of Monte Carlo (MC) steps
      public double delay, stop; // MC steps after which statistics are collected and not collected
      public int snapshotInterval; // interval by which successive samples are separated
      public double sizeDist[], sizeBinWidth; // particle radius histogram and bin width
      public int numberBins; // number of histogram bins
      public double grBinWidth, maxRadius, deltaK; // bin widths and range
      public double meanR; // mean particle radius
      public double lambda, dlambda; // coupling constant parameters
      public double mixFRSR, elasticFRSR, totalFRSR;
      public double pairEnergySum; //parameters for the lambda = 1 system
      public double initialEnergy; // initialEnergy when the particles are on their lattice sites
      public double displacement, trialDisplacementDistance; // difference in coordinates
      public double springEnergy, springEnergyAccumulator, springEnergySum, springConstant;//the parameters to calculate the spring energy and accumulate it
      public double totalEinsteinEnergy, totalEinsteinPotential; //the potential energy associated with the eistein solid
      public double dSpringEnergy, springEnergy0; // the spring energy parameters
      public double pairPotentialCorrection; // change in harmonic potential energy
      public double boltzmannFactorAccumulator, boltzmannFactor, squaredDisplacementAccumulator, volFracAccumulator;
      public double dxOverN, dyOverN, dzOverN; // change in displacement per particle trial move
      public String fileExtension; // file extension for output files
      public double numberOfConfigurations; // number of configurations collected for averaging
      public double density, nnDistance;  // density of particles and nearest neighbor distance
      public double squaredDisplacement, squaredDisplacementSum;  // squared displacement of particles
      public double oldFloryFR, newFloryFR;  // old and new Flory-Rehner free energies
      public double volumeOfSolvent;  // volume of solvent in units of dry radius cubed
      public double mixF_FR, mixF_Interpenetration;   // mixing free energy terms

   /**
    * Initialize the model.
    *
    * @param configuration
    * Initial lattice structure
    */
   public void initialize(String configuration) {
      x = new double[N]; // particle coordinates
      y = new double[N];
      z = new double[N];
      
      //the coordinates for the new particles
      x0 = new double[N];
      y0 = new double[N];
      z0 = new double[N];

      a = new double[N]; // particle radii (units of dry radius)

      monRadius = 0.3; // estimate of monomer radius [nm]
      nMon = 0.63*Math.pow(dryR, 3.)/Math.pow(monRadius, 3.); // number of monomers in a microgel
      nChains = xLinkFrac*nMon; // number of chains in a microgel
      reservoirSR = reservoirSwellingRatio(nMon, nChains, chi);
      volumeOfSolvent = ((4*Math.PI)/3.0 * Math.pow(monRadius, 3))/Math.pow(dryR, 3); // in units of dry radius

      for (int i=0; i<N; i++){
	      a[i] = reservoirSR; // initial particle radii (fully swollen)
      }

      energy = new double[N]; // total energies of particles
      pairEnergy = new double[N][N]; // pair energies of particles
      newPairEnergy = new double[N][N]; // for lambda = 1
      steps = 0;
      energyAccumulator = 0;
      pairEnergyAccumulator = 0; //for the lamba = 1 system
      freeEnergyAccumulator = 0;
      virialAccumulator = 0;
      springEnergyAccumulator = 0;
      boltzmannFactorAccumulator = 0;
      numberOfConfigurations = 0;
      squaredDisplacementAccumulator = 0;
      volFracAccumulator = 0;
      dxOverN = 0; // change in displacement per particle trial move
      dyOverN = 0;
      dzOverN = 0;

      side = Math.cbrt(4.*Math.PI*N/dryVolFrac/3.); // side length of cubic simulation box
      totalVol = side*side*side; // box volume [units of dry radius cubed]

      numberBins = (int) (maxRadius/sizeBinWidth);

      sizeDist = new double[(int) (numberBins)];

      for(int i=0; i<numberBins; i++){ // initialize size histogram
         sizeDist[i] = 0;
      }

      // initialize positions
      if (configuration.toUpperCase().equals("SC")) {
         scale = 1; // no scaling of Hertz interactions
         setSCpositions();
      }
      if (configuration.toUpperCase().equals("FCC")) {
         scale = 1; // no scaling of Hertz interactions
         setFCCpositions();
      }
      if (configuration.toUpperCase().equals("BCC")) {
         scale = 1; // no scaling of Hertz interactions
         setBCCpositions();
      }
      if (configuration.toUpperCase().equals("RANDOM-FCC")) {
         scale = 0; // particles initially do not interact to allow positions to randomize
         setFCCrandomPositions();
      }
      if (configuration.toUpperCase().equals("RANDOM-BCC")) {
         scale = 0; // particles initially do not interact to allow positions to randomize
         setBCCrandomPositions();
      }
   }

   /**
    * Place particles on sites of a simple cubic lattice.
    */
   public void setSCpositions() {
      System.out.println("SC");
      int ix, iy, iz;
      double dnx = Math.cbrt(N);
      d = side / dnx; // distance between neighboring lattice sites
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // N is not a perfect cube
      }

      int i = 0;
      for (iy = 0; iy < nx; iy++) { // loop through particles in a column
         for (ix = 0; ix < nx; ix++) { // loop through particles in a row
            for (iz = 0; iz < nx; iz++) {
               if (i < N) { // check for remaining particles
                  x[i] = ix * d;
                  y[i] = iy * d;
                  z[i] = iz * d;
                  i++;
               }
            }
         }
      }
      calculateTotalEnergy(lambda); // initial energy
   }

   /**
    * Place particles on sites of an FCC lattice.
    */
   public void setFCCpositions() {
      System.out.println("FCC");
      int ix, iy, iz;
      double dnx = Math.cbrt(N/4.);
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // N/4 is not a perfect cube
      }

      int i = 0;
      for (ix = 0; ix < 2*nx; ix++) { // loop through particles in a row
         for (iy = 0; iy < 2*nx; iy++) { // loop through particles in a column
            for (iz = 0; iz < 2*nx; iz++) { // loop through particles in a layer
               if (i < N) { // check for remaining particles
                  if ((ix+iy+iz)%2 == 0) { // check for remaining particles
                     //the final equilibrium positions of all the particles
                     x[i] = ix * d/2.;
                     y[i] = iy * d/2.;
                     z[i] = iz * d/2.;

                     // initial displacements
                     x0[i] = x[i];
                     y0[i] = y[i];
                     z0[i] = z[i];
                     i++;

                  }
               }
            }
         }
      }
      calculateTotalEnergy(lambda);
      initialEnergy = totalPairEnergy; // initial energy when lambda = 1
   }

   /**
   * Place particles on sites of a BCC lattice.
   */
   public void setBCCpositions() {
      System.out.println("BCC");
      int ix, iy, iz;
      double dnx = Math.cbrt(N/2.);
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // N/2 is not a perfect cube
      }

      int i = 0;
      for (ix = 0; ix < 2*nx; ix++) { // loop through particles in a row
         for (iy = 0; iy < 2*nx; iy++) { // loop through particles in a column
            for (iz = 0; iz < 2*nx; iz++) { // loop through particles in a layer
               if (i < N) { // check for remaining particles
                  if ((ix*iy*iz)%2 == 1 || (ix%2 == 0 && iy%2 == 0 && iz%2 ==0)) {
                     x[i] = ix * d/2.;
                     y[i] = iy * d/2.;
                     z[i] = iz * d/2.;
                     i++;
                  }
               }
            }
         }
      }
      calculateTotalEnergy(lambda); // initial energy
   }

   /**
    * Place particles at positions randomly displaced from sites of an FCC lattice.
    */
   public void setFCCrandomPositions() {
      System.out.println("random FCC");
      int ix, iy, iz;
      double dnx = Math.cbrt(N/4.);
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // N/4 is not a perfect cube
      }

      int i = 0;
      for (ix = 0; ix < 2*nx; ix++) { // loop through particles in a row
         for (iy = 0; iy < 2*nx; iy++) { // loop through particles in a column
            for (iz = 0; iz < 2*nx; iz++) { // loop through particles in a layer
               if (i < N) { // check for remaining particles
                  if ((ix+iy+iz)%2 == 0) { // check for remaining particles
                     x[i] = PBC.position(ix * d/2. + (Math.random()-0.5) * d, side);
                     y[i] = PBC.position(iy * d/2. + (Math.random()-0.5) * d, side);
                     z[i] = PBC.position(iz * d/2. + (Math.random()-0.5) * d, side);
                     i++;
                  }
               }
            }
         }
      }
      calculateTotalEnergy(lambda); // initial energy
   }

   /**
   * Place particles at positions randomly displaced from sites of a BCC lattice.
   */
   public void setBCCrandomPositions() {
      System.out.println("random BCC");
      int ix, iy, iz;
      double dnx = Math.cbrt(N/2.);
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // N/2 is not a perfect cube
      }

      int i = 0;
      for (ix = 0; ix < 2*nx; ix++) { // loop through particles in a row
         for (iy = 0; iy < 2*nx; iy++) { // loop through particles in a column
            for (iz = 0; iz < 2*nx; iz++) { // loop through particles in a layer
               if (i < N) { // check for remaining particles
                  if ((ix*iy*iz)%2 == 1 || (ix%2 == 0 && iy%2 == 0 && iz%2 ==0)) { // check for remaining particles
                     x[i] = PBC.position(ix * d/2. + (Math.random()-0.5) * d, side);
                     y[i] = PBC.position(iy * d/2. + (Math.random()-0.5) * d, side);
                     z[i] = PBC.position(iz * d/2. + (Math.random()-0.5) * d, side);
                     i++;
                  }
               }
            }
         }
      }
      calculateTotalEnergy(lambda); // initial energy
   }

      /**
    * Do a Monte Carlo simulation step.
    */
   public void step() { // performs a trial move of every particle
      steps++;
      double dxtrial;
      double dytrial;
      double dztrial;
      double datrial;
      double xij;
      double yij;
      double zij;
      double de;
      double r;
      double r2;
      double sigma;
      double mixF;
      double elasticF;
      double totalF;
      double capVolSum;

      for (int i = 0; i < N; i++) { // attempt a trial move (displacement and size change)

         dxtrial = tolerance * 2.0 * (Math.random() - 0.5);
         dytrial = tolerance * 2.0 * (Math.random() - 0.5);
         dztrial = tolerance * 2.0 * (Math.random() - 0.5);
         datrial = atolerance * 2.0 * (Math.random() - 0.5);

         // Euclidean distance of trial displacements
         trialDisplacementDistance = dxtrial*dxtrial + dytrial*dytrial + dztrial*dztrial;

         // displacement relative to initial positions (with COM shift)
         double dxi = x[i] - x0[i] - dxOverN;
         double dyi = y[i] - y0[i] - dyOverN;
         double dzi = z[i] - z0[i] - dzOverN;

         // Flory-Rehner free energy BEFORE size changes
         double ai0 = a[i];
         double ai0_2 = ai0 * ai0;
         double ai0Cubed = ai0_2 * ai0;
         mixF    = nMon * ((ai0Cubed - 1.0) * Math.log(1.0 - 1.0 / ai0Cubed) + chi * (1.0 - 1.0 / ai0Cubed));
         elasticF = 1.5 * nChains * (ai0_2 - Math.log(ai0) - 1.0);
         oldFloryFR = elasticF + mixF;

         // apply trial move
         x[i] += dxtrial;
         y[i] += dytrial;
         z[i] += dztrial;
         a[i] += datrial;

         // spring energy correction
         double dot = dxtrial*dxi + dytrial*dyi + dztrial*dzi;
         pairPotentialCorrection = 2.0 * springConstant * Math.abs(dot)+ (1.0 - 1.0 / N) * trialDisplacementDistance;

         mixF_Interpenetration = 0.0;

         double ai = a[i];
         double ai2 = ai * ai;
         double ai3 = ai2 * ai;
         double phipi = 1.0 / ai3;

         for (int j = 0; j < N; j++) {
               if (j == i) continue; // consider interactions with other particles

               xij = PBC.separation(x[i] - x[j], side);
               yij = PBC.separation(y[i] - y[j], side);
               zij = PBC.separation(z[i] - z[j], side);
               r2  = xij*xij + yij*yij + zij*zij;
               sigma = ai + a[j];
               double sigma2 = sigma * sigma;

               // Define phipi and phipj without Math.pow
               double aj = a[j];
               double aj2 = aj * aj;
               double aj3 = aj2 * aj;
               double phipj = 1.0 / aj3;

               // Rewriting the expressions in terms of phipi and phipj
               double firstMixOverlap   = (1.0 - phipi - phipj) * Math.log(1.0 - phipi - phipj);
               double secondMixOverlap  = chi * (1.0 - phipi - phipj) * (phipi + phipj);
               double firstMixNoOverlap = (1.0 - phipi) * Math.log(1.0 - phipi) + chi * (1.0 - phipi) * phipi;
               double secondMixNoOverlap= (1.0 - phipj) * Math.log(1.0 - phipj) + chi * (1.0 - phipj) * phipj;

               if (r2 < sigma2) {
                  r = Math.sqrt(r2);

                  // prefactor B
                  B = scale * Young * nChains * (sigma * sigma) * Math.sqrt(ai * aj) / (ai3 + aj3);

                  // heights of the caps
                  double hi = (aj - ai + r) * (aj + ai - r) / (2.0 * r);
                  double hj = (ai - aj + r) * (ai + aj - r) / (2.0 * r);

                  // cap volumes
                  double hi2 = hi * hi, hj2 = hj * hj;
                  double capVolForI = (Math.PI * hi2 * (3.0 * ai - hi)) / 3.0;
                  double capVolForJ = (Math.PI * hj2 * (3.0 * aj - hj)) / 3.0;
                  capVolSum = capVolForI + capVolForJ;

                  mixF_FR = (capVolSum / volumeOfSolvent)
                           * (firstMixOverlap + secondMixOverlap - firstMixNoOverlap - secondMixNoOverlap);

                  newPairEnergy[i][j] = mixF_FR;
                  mixF_Interpenetration += mixF_FR;
               } else {
                  mixF_FR = 0.0;
                  newPairEnergy[i][j] = 0.0;
                  
               }
         }

         // Flory-Rehner single-particle free energy AFTER the trial move
         double ai_new = a[i];
         double ai_new2 = ai_new * ai_new;
         double ai_new3 = ai_new2 * ai_new;
         mixF = nMon * ((ai_new3 - 1.0) * Math.log(1.0 - 1.0 / ai_new3) + chi * (1.0 - 1.0 / ai_new3));
         elasticF = 1.5 * nChains * (ai_new2 - Math.log(ai_new) - 1.0);
         newFloryFR = elasticF + mixF;

         // Metropolis criterion
         de = lambda * (mixF_Interpenetration - energy[i]) + (1.0 - lambda) * (pairPotentialCorrection) + (newFloryFR - oldFloryFR);

         if (Math.exp(-de) < Math.random()) { // reject move
               x[i] -= dxtrial;
               y[i] -= dytrial;
               z[i] -= dztrial;
               a[i] -= datrial;
         } else { // accept move and update energies
               dxOverN += dxtrial / (double) N;
               dyOverN += dytrial / (double) N;
               dzOverN += dztrial / (double) N;

               energy[i] = mixF_Interpenetration;

               for (int j = 0; j < N; ++j) {
                  if (j == i) continue;
                  energy[j] += newPairEnergy[i][j] - pairEnergy[i][j];
                  pairEnergy[i][j] = newPairEnergy[i][j];
                  pairEnergy[j][i] = newPairEnergy[i][j];
               }
         }
      }
      calculateTotalEnergy(lambda); // new total energy
   }


   public void calculateTotalEnergy(double lambda) {
      totalEnergy = 0;
      totalVirial = 0;
      totalPairEnergy = 0;
      totalFreeEnergy = 0;
      springEnergySum = 0;
      squaredDisplacementSum = 0;

      for (int i = 0; i < N; ++i) {
         mixF_Interpenetration = 0;
         double virialSum = 0;

         // displacement from initial position
         double dxi = x[i] - x0[i] - dxOverN;
         double dyi = y[i] - y0[i] - dyOverN;
         double dzi = z[i] - z0[i] - dzOverN;

         // squared displacement + spring energy
         squaredDisplacement = dxi*dxi + dyi*dyi + dzi*dzi;
         squaredDisplacementSum += squaredDisplacement;
         springEnergy = springConstant * squaredDisplacement;
         springEnergySum += springEnergy;

         for (int j = 0; j < N; ++j) {
               if (j == i) continue;  // skip self

               double xij = PBC.separation(x[i] - x[j], side);
               double yij = PBC.separation(y[i] - y[j], side);
               double zij = PBC.separation(z[i] - z[j], side);
               double r2  = xij*xij + yij*yij + zij*zij;
               double sigma = a[i] + a[j];

               // avoid repeated Math.pow
               double ai3 = a[i] * a[i] * a[i];
               double aj3 = a[j] * a[j] * a[j];
               double phipi = 1.0 / ai3;
               double phipj = 1.0 / aj3;

               // mixing terms
               double overlapTerm = 1 - phipi - phipj;
               double firstMixOverlap = overlapTerm * Math.log(overlapTerm);
               double secondMixOverlap = chi * overlapTerm * (phipi + phipj);
               double firstMixNoOverlap = (1 - phipi) * Math.log(1 - phipi) + chi * (1 - phipi) * phipi;
               double secondMixNoOverlap = (1 - phipj) * Math.log(1 - phipj) + chi * (1 - phipj) * phipj;

               if (r2 < sigma * sigma) {
                  double r = Math.sqrt(r2);
                  B = scale * Young * nChains * (sigma * sigma) * Math.sqrt(a[i] * a[j]) / (ai3 + aj3);

                  // lens cap volumes
                  double hi = (a[j] - a[i] + r) * (a[j] + a[i] - r) / (2.0 * r);
                  double hj = (a[i] - a[j] + r) * (a[i] + a[j] - r) / (2.0 * r);
                  double capVolForI = (Math.PI * hi * hi * (3.0 * a[i] - hi)) / 3.0;
                  double capVolForJ = (Math.PI * hj * hj * (3.0 * a[j] - hj)) / 3.0;
                  double capVolSum = capVolForI + capVolForJ;

                  mixF_FR = (capVolSum / volumeOfSolvent) *
                           (firstMixOverlap + secondMixOverlap - firstMixNoOverlap - secondMixNoOverlap);
               } else {
                  mixF_FR = 0;
               }

               pairEnergy[i][j] = mixF_FR;
               mixF_Interpenetration += mixF_FR;
         }

         // Flory-Rehner free energy
         double ai2 = a[i] * a[i];
         double ai3 = a[i] * ai2;
         double mixF = nMon * ((ai3 - 1) * Math.log(1 - 1.0 / ai3) + chi * (1 - 1.0 / ai3));
         double elasticF = 1.5 * nChains * (ai2 - Math.log(a[i]) - 1);
         double totalF = elasticF + mixF;

         energy[i] = mixF_Interpenetration;

         // reservoir correction
         double resSR2 = reservoirSR * reservoirSR;
         double resSR3 = resSR2 * reservoirSR;
         mixFRSR = nMon * ((resSR3 - 1) * Math.log(1 - 1.0 / resSR3) + chi * (1 - 1.0 / resSR3));
         elasticFRSR = 1.5 * nChains * (resSR2 - Math.log(reservoirSR) - 1);
         totalFRSR = mixFRSR + elasticFRSR;

         totalF -= totalFRSR;

         totalFreeEnergy += totalF;
         totalPairEnergy += mixF_Interpenetration;
         totalVirial += virialSum;
      }

      if (steps > delay && (steps - delay) % snapshotInterval == 0) {
         numberOfConfigurations++;
         totalVirial *= 0.5; // correct for double counting
         totalEnergy = totalPairEnergy + totalFreeEnergy;

         pairEnergyAccumulator += totalPairEnergy;
         energyAccumulator += totalEnergy;
         freeEnergyAccumulator += totalFreeEnergy;
         virialAccumulator += totalVirial;
         springEnergyAccumulator += springEnergySum;
         squaredDisplacementAccumulator += squaredDisplacementSum;
         volFracAccumulator += calculateVolumeFraction();
      }
   }

   // mean energy per particle [kT units]
   public double meanEnergy() {
      return energyAccumulator/N/numberOfConfigurations; // quantity <E>/N
   }

    // mean microgel volume fraction
   public double meanVolFrac() {
      return volFracAccumulator/numberOfConfigurations; // Total volume fraction over number of configurations
   }

   // mean pair energy per particle [kT units]
   public double meanPairEnergy() {
      return pairEnergyAccumulator/N/numberOfConfigurations; // quantity <E_pair>/N
   }

   // einstein  free energy per particle (KT units)
   public double einsteinFreeEnergy() {
      // Assuming the dimension is 3D
      int d = 3; //the dimension
      return (initialEnergy-(d/2.0)*N*Math.log(Math.PI/springConstant))/N;
   }

   public double meanSquareDisplacement(){
      return squaredDisplacementAccumulator/N/numberOfConfigurations; // quantity <r>^2
   }

   public double meanFreeEnergy() { // mean free energy per particle
      return freeEnergyAccumulator/N/numberOfConfigurations; // quantity <F>/N
   }

   public double meanboltzmannFactor(){ // mean boltzmannFactor
      return boltzmannFactorAccumulator/numberOfConfigurations;
   }

   public double meanSpringEnergy(){ // mean springEnergyAccumulator
      return springEnergyAccumulator/N/numberOfConfigurations;
   }

   public double meanPressure() { // mean pressure (dimensionless)
      double meanVirial;
      meanVirial = virialAccumulator/numberOfConfigurations;
      return 1+(1./3.)*meanVirial/N; // quantity PV/NkT //the 1 is coming from the ideal gas
   }

   public void sizeDistribution() {
      int bin;
      for (int i=0; i<N; i++){
         bin = (int) (a[i]/sizeBinWidth);
         sizeDist[bin]++;
      }
   }

   public double meanRadius() { // mean radius [units of dry radius]
      double sum = 0, sumr = 0;
      if (steps > delay){
         for (int i=0; i<numberBins; i++){
            sum += sizeDist[i];
            sumr += i*sizeBinWidth*sizeDist[i];
         }
         meanR = sumr/sum;
        // System.out.println("mberBins " + numberBins);
      }
      else
      meanR = 0;
        // System.out.println("meanR " + meanR);
      return meanR;
   }

    /* compute mean volume fraction after stopping */
    public double calculateVolumeFraction() { // instantaneous volume fraction
         double microgelVol, overlapVol;
         double xij, yij, zij, r, r2, sigma, da2, amin;

         microgelVol = 0;
         for (int i=0; i<N; i++) {
               microgelVol += (4./3.)*Math.PI*a[i]*a[i]*a[i]; // sum volumes of spheres
               for (int j=i+1; j<N; j++){
                  xij = PBC.separation(x[i]-x[j], side);
                  yij = PBC.separation(y[i]-y[j], side);
                  zij = PBC.separation(z[i]-z[j], side);
                  r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
                  sigma = a[i]+a[j]; // sum of radii of two particles [units of dry radius]
                  da2 = Math.pow(a[i]-a[j], 2); // difference of radii squared
                  if (r2 < sigma*sigma){ // particles are overlapping
                     if (r2 < da2){ // one particle is entirely inside the other
                           amin = Math.min(a[i],a[j]); // minimum of radii
                           overlapVol = (4./3.)*Math.PI*amin*amin*amin; // volume of smaller particle
                     }
                     else { // one particle is NOT entirely inside the other
                           r = Math.sqrt(r2);
                           overlapVol = Math.PI*(sigma-r)*(sigma-r)*(r2+2*r*sigma-3*da2)/12./r; // volume of lens-shaped overlap region
                     }
                     microgelVol -= overlapVol; // subtract overlap volume so we don't double-count
                  }
               }
         }
         return microgelVol/totalVol;
      }

      /* reservoir swelling ratio as root of Flory-Rehner pressure */
      public double reservoirSwellingRatio(double nMon, double nChains, double chi){
         Function f = new FloryRehnerPressure(nMon, nChains, chi);
         double xleft = 1.1;
         double xright = 11.;
         double epsilon = 1.e-06;
         double x = Root.bisection(f, xleft, xright, epsilon);
         return x;
      }

      class FloryRehnerPressure implements Function{
         double nMon, nChains, chi;

         /**
            * Constructs the FloryRehnerPressure function with the given parameters.
            * @param _nMon double
            * @param _nChains double
            * @param _chi double
            */
         FloryRehnerPressure(double _nMon, double _nChains, double _chi) {
            nMon = _nMon;
            nChains = _nChains;
            chi = _chi;
         }

         public double evaluate(double x) {
            double pMixing = nMon*(x*x*x*Math.log(1.-1./x/x/x)+1.+chi/x/x/x);
            double pElastic = nChains*(x*x-0.5);
            double pressure = pMixing + pElastic;
            return pressure;
         }

      }

}
/*
 * Open Source Physics software is free software; you can redistribute
 * it and/or modify it under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation; either version 2 of the License,
 * or(at your option) any later version.

 * Code that uses any portion of the code in the org.opensourcephysics package
 * or any subpackage (subdirectory) of this package must must also be be released
 * under the GNU GPL license.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
 * or view the license online at http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2007  The Open Source Physics project
 *                     http://www.opensourcephysics.org
 */
