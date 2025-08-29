/**
 * HertzSpheresFCCNearestNeighbors.java
 * 
 * This class simulates microgel particles arranged on a Face-Centered Cubic (FCC) lattice.
 * It performs Monte Carlo simulations where each microgel interacts only with its 12 nearest neighbors.
 *
 * Key Features:
 * - FCC lattice configuration with nearest-neighbor precomputation.
 * - Monte Carlo steps include particle displacement, radius adjustments, and energy computations.
 * - Efficient simulation of swelling, overlap, and elastic interactions using the Flory-Rehner theory and Hertzian potentials.
 *
 * Usage:
 *   - Initialize particles using FCC lattice (`setFCCpositions`).
 *   - Perform Monte Carlo simulation steps (`step` method).
 *
 * Authors: Oreoluwa Alade, Alan Denton
 * Version: 2.0 (Updated April 2025)
 */

package org.opensourcephysics.sip.Hertz;

import java.awt.Color;
import java.text.DecimalFormat;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

// import java.util.function.Function;
import org.opensourcephysics.numerics.*;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;

public class HertzSpheresNonLocalFacet_Free_Energies {
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
      public double dxOverN, dyOverN, dzOverN; // center of mass
      public String fileExtension;
      public double numberOfConfigurations;
      public double density;
      public double squaredDisplacement, squaredDisplacementSum;
      public double oldFloryFR, newFloryFR; // Old and new Flory-Rehner free energy
      public Random random = new Random(); // Random number generator for Monte Carlo moves
      public double volOfMicrogel[]; //volume of each microgel
      // Arrays to hold Flory-Rehner free energy for microgels after a move (for J and I)
      public double floryFRJAfterMove[], floryFRJBeforeMove[];
      public double floryFRIBeforeMove, floryFRIAfterMove;
      public double newNMon[], cumulativeCapVol[], capVol[]; //new number of monomers, cap volume and cummulative cap volume for each microgel
      public double floryFRAfterMove[], floryFRBeforeMove[]; //Flory-Rehner free energy for microgels before and after a move
      public double capVolJ[], capVolJBeforeMove[], capVolJAfterMove[];
      public double capVolSumBeforeMoveArray[], capVolSumAfterMoveArray[];
      public double capVolAfterTrialMove[][], capVolBeforeTrialMove[][]; // stores the cap volume for each j microgel before and after a move
      public double totalVolFrac;
      public double nonOverlappingVol, volumeFraction;
      public int[][] nearestNeighbors; // Each microgel has exactly 12 nearest neighbors
      

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
      floryFRAfterMove = new double[N];
      floryFRBeforeMove = new double[N];

      capVolSumBeforeMoveArray = new double[N];
      capVolSumAfterMoveArray = new double[N];
      capVolJ = new double[N];
      capVolJBeforeMove = new double[N];
      capVolJAfterMove = new double[N];
      capVol = new double[N];
      newNMon = new double[N];
      cumulativeCapVol = new double[N];
      capVolAfterTrialMove = new double[N][N];
      capVolBeforeTrialMove = new double[N][N];
      volOfMicrogel = new double[N];
      floryFRJAfterMove = new double[N];
      floryFRJBeforeMove = new double[N];

      a = new double[N]; // particle radii (units of dry radius)
      monRadius = 0.3; // estimate of monomer radius [nm]
      nMon = 0.63*Math.pow(dryR, 3.)/Math.pow(monRadius, 3.); // number of monomers in a microgel
      //nMon = Math.pow(dryR, 3.)/Math.pow(monRadius, 3.); // number of monomers in a microgel
      nChains = xLinkFrac*nMon; // number of chains in a microgel
      reservoirSR = reservoirSwellingRatio(nMon, nChains, chi);

      reservoirVolFrac = dryVolFrac*Math.pow(reservoirSR, 3.);
      // // Initialize radii with the fully swollen reservoirSR
      for (int i = 0; i < N; i++) {
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
      //springEnergyAccumulator = 0;
      boltzmannFactorAccumulator = 0;
      numberOfConfigurations = 0;
      squaredDisplacementAccumulator = 0;
      volFracAccumulator = 0;
      dxOverN = 0; // change in displacement per particle trial move
      dyOverN = 0;
      dzOverN = 0;
      volFrac = 0; // counter for system volume fraction

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
      // Now, precompute neighbors
      computeNearestNeighbors();
      calculateTotalEnergy(lambda);
      initialEnergy = totalPairEnergy; // initial energy when lambda = 1
   }

   /**
    * Precompute the 12 nearest neighbors for each microgel particle on an FCC lattice.
    *
    * For each particle:
    * - Calculate the Euclidean distance squared to every other particle
    * - Sort and identify exactly the 12 nearest particles
    * - Store their indices for efficient computation of interactions during simulation steps
    */
   private void computeNearestNeighbors() {
      nearestNeighbors = new int[N][12]; // FCC lattice has exactly 12 nearest neighbors
      for (int i = 0; i < N; i++) {
         double[] distances = new double[N]; // distances squared from particle i to others
         int[] indices = new int[N]; // particle indices

         // Calculate distance squared to each particle
         for (int j = 0; j < N; j++) {
            if (j == i) {
               distances[j] = Double.MAX_VALUE; // particle cannot be its own neighbor
            } else {
                  double dx = PBC.separation(x[i] - x[j], side);
                  double dy = PBC.separation(y[i] - y[j], side);
                  double dz = PBC.separation(z[i] - z[j], side);
                  distances[j] = dx * dx + dy * dy + dz * dz;
            }
            indices[j] = j;
         }

         // Sort distances and identify the 12 nearest neighbors
         sortNeighbors(distances, indices);

         // Store the indices of the 12 closest neighbors for particle i
         System.arraycopy(indices, 0, nearestNeighbors[i], 0, 12);
      }
   }


   /**
    * Sorts particle indices based on their distance from a reference particle (ascending order).
    * Utilizes a simple selection sort, efficient enough for small N.
    *
    * @param distances array of squared distances from the reference particle
    * @param indices array of particle indices corresponding to distances
    */
   private void sortNeighbors(double[] distances, int[] indices) {
      int n = distances.length;
      for (int i = 0; i < 12; i++) {  // Only sort until the first 12 nearest neighbors are found
         int minIndex = i;
         for (int j = i + 1; j < n; j++) {
            if (distances[j] < distances[minIndex]) {
                  minIndex = j;  // found smaller distance
            }
         }

         // swap distance elements
         double tmpDist = distances[i];
         distances[i] = distances[minIndex];
         distances[minIndex] = tmpDist;

         // swap corresponding indices
         int tmpIdx = indices[i];
         indices[i] = indices[minIndex];
         indices[minIndex] = tmpIdx;
      }
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

      // small constants & locals reused a lot
      final double FOUR_THIRDS_PI = 4.0 / 3.0 * Math.PI;
      
      double dxtrial;
      double dytrial;
      double dztrial;
      double datrial;
      double xij, yij, zij, de, r, r2, sigma;
      double mixF, elasticF;
      double capVolForI;
      double capVolForJ;
      int i, j;
      double dFlory, newSWR;
      double capVolSumAfterMoveI;

      // ---------- prepass: cap volumes & FR energy (BEFORE move) ----------
      for (i = 0; i < N; i++) {
         double ai = a[i];
         double ai2 = ai * ai; // cache powers of ai
         double ai3 = ai2 * ai;

         double capVolSumBeforeMoveI = 0.0;
         volOfMicrogel[i] = FOUR_THIRDS_PI * ai3;   // cache vol(i)

         // loop through exactly 12 nearest neighbours of particle i
         for (int neighborIdx = 0; neighborIdx < 12; neighborIdx++) {
               j = nearestNeighbors[i][neighborIdx];

               double aj = a[j];
               double aj2 = aj * aj;
               double aj3 = aj2 * aj;
               volOfMicrogel[j] = FOUR_THIRDS_PI * aj3;

               xij = PBC.separation(x[i] - x[j], side);
               yij = PBC.separation(y[i] - y[j], side);
               zij = PBC.separation(z[i] - z[j], side);

               r2 = xij * xij + yij * yij + zij * zij;
               sigma = ai + aj;
               double sigma2 = sigma * sigma;

               if (r2 < sigma2) { // overlap: compute cap volumes
                  r = Math.sqrt(r2);
                  // heights of the caps
                  double aj_minus_ai_plus_r = (aj - ai + r);
                  double aj_plus_ai_minus_r = (aj + ai - r);
                  double hi = (aj_minus_ai_plus_r * aj_plus_ai_minus_r) / (2.0 * r);

                  double ai_minus_aj_plus_r = (ai - aj + r);
                  double ai_plus_aj_minus_r = (ai + aj - r);
                  double hj = (ai_minus_aj_plus_r * ai_plus_aj_minus_r) / (2.0 * r);

                  double hi2 = hi * hi, hj2 = hj * hj;

                  capVolForI = (Math.PI * hi2 * (3.0 * ai - hi)) / 3.0; // cap volume of i due to j
                  capVolForJ = (Math.PI * hj2 * (3.0 * aj - hj)) / 3.0; // cap volume of j due to i

                  capVolSumBeforeMoveI += capVolForI;
                  capVolBeforeTrialMove[i][j] = capVolForJ;
                  capVolBeforeTrialMove[j][i] = capVolForI;
               } else {
                  capVolBeforeTrialMove[i][j] = 0.0; // no overlap
                  capVolBeforeTrialMove[j][i] = 0.0;
               }
         }

         capVolSumBeforeMoveArray[i] = capVolSumBeforeMoveI; // store for later use

         // FR free energy BEFORE the move for microgel i (with precomputed powers)
         double nonOverlapFrac = 1.0 - (volOfMicrogel[i] > 0.0 ? (capVolSumBeforeMoveArray[i] / volOfMicrogel[i]) : 0.0);
         // guard against tiny negatives due to fp error
         if (nonOverlapFrac < 0.0) nonOverlapFrac = 0.0;
         newSWR = ai * Math.cbrt(nonOverlapFrac);

         double swr2 = newSWR * newSWR;
         double swr3 = swr2 * newSWR;
         double inv_swr3 = 1.0 / swr3;

         mixF     = nMon * ((swr3 - 1.0) * Math.log(1.0 - inv_swr3) + chi * (1.0 - inv_swr3));
         elasticF = 1.5 * nMon * xLinkFrac * (swr2 - Math.log(newSWR) - 1.0);
         floryFRBeforeMove[i] = mixF + elasticF;
      }

      // ---------- main MC loop: propose & accept/reject ----------
      for (i = 0; i < N; i++) {
         double oldFreeEnergyI = floryFRBeforeMove[i];

         // trial displacements
         dxtrial = tolerance * (2.0 * Math.random() - 1.0);
         dytrial = tolerance * (2.0 * Math.random() - 1.0);
         dztrial = tolerance * (2.0 * Math.random() - 1.0);
         datrial = atolerance * (2.0 * Math.random() - 1.0);

         // squared Euclidean length of trial displacement
         double dx2 = dxtrial * dxtrial, dy2 = dytrial * dytrial, dz2 = dztrial * dztrial;
         trialDisplacementDistance = dx2 + dy2 + dz2;

         // apply trial move
         x[i] += dxtrial;
         y[i] += dytrial;
         z[i] += dztrial;
         a[i] += datrial;

         // displacement relative to initial positions (with COM shift)
         double dxi = x[i] - x0[i] - dxOverN;
         double dyi = y[i] - y0[i] - dyOverN;
         double dzi = z[i] - z0[i] - dzOverN;

         // |dot| and spring correction
         double dot = dxtrial * dxi + dytrial * dyi + dztrial * dzi;
         // correction to the spring energy due to COM shift
         pairPotentialCorrection = 2.0 * springConstant * Math.abs(dot) + (1.0 - 1.0 / (double) N) * trialDisplacementDistance;

         // update vol(i) with new radius a[i]
         double ai = a[i];
         double ai2 = ai * ai;
         double ai3 = ai2 * ai;
         volOfMicrogel[i] = FOUR_THIRDS_PI * ai3;

         pairEnergySum = 0.0;
         capVolSumAfterMoveI = 0.0;

         double dFloryForJ = 0.0;

         for (int neighborIdx = 0; neighborIdx < 12; neighborIdx++) { // loop through exactly 12 nearest neighbours of particle i
               j = nearestNeighbors[i][neighborIdx]; // index of neighbor j

               // NOTE: a[j] did not change in this trial; volOfMicrogel[j] from prepass is valid.
               double aj = a[j];
               double aj2 = aj * aj;
               double aj3 = aj2 * aj;

               xij = PBC.separation(x[i] - x[j], side);
               yij = PBC.separation(y[i] - y[j], side);
               zij = PBC.separation(z[i] - z[j], side);

               r2 = xij * xij + yij * yij + zij * zij;
               sigma = ai + aj;
               double sigma2 = sigma * sigma;

               double oldFreeEnergyJ = floryFRBeforeMove[j]; // old FR free energy for microgel j (before move of i)

               if (r2 < sigma2) {
                  r = Math.sqrt(r2);
                  // Hertz energy
                  double Bamp = scale * Young * nChains * (sigma2) * Math.sqrt(ai * aj) / (ai3 + aj3);
                  double oneMinus = 1.0 - r / sigma;
                  double sqrtOneMinus = Math.sqrt(Math.max(0.0, oneMinus));
                  double derivPart = Bamp * oneMinus * sqrtOneMinus;
                  double HertzEnergy = derivPart * oneMinus;
                  newPairEnergy[i][j] = HertzEnergy;
                  pairEnergySum += HertzEnergy;

                  // cap volumes hi, hj
                  double aj_minus_ai_plus_r = (aj - ai + r);
                  double aj_plus_ai_minus_r = (aj + ai - r);
                  double hi = (aj_minus_ai_plus_r * aj_plus_ai_minus_r) / (2.0 * r);

                  double ai_minus_aj_plus_r = (ai - aj + r);
                  double ai_plus_aj_minus_r = (ai + aj - r);
                  double hj = (ai_minus_aj_plus_r * ai_plus_aj_minus_r) / (2.0 * r);

                  double hi2 = hi * hi, hj2 = hj * hj;
                  capVolForJ = (Math.PI * hj2 * (3.0 * aj - hj)) / 3.0;
                  capVolForI = (Math.PI * hi2 * (3.0 * ai - hi)) / 3.0;

                  capVolSumAfterMoveI += capVolForI;
                  capVolAfterTrialMove[i][j] = capVolForJ;
                  capVolAfterTrialMove[j][i] = capVolForI;

                  // FR energy for j after the move (i moved)
                  double newCapVolJ = capVolSumBeforeMoveArray[j] + capVolAfterTrialMove[i][j] - capVolBeforeTrialMove[i][j];
                  double nonOverlapFracJ = 1.0 - newCapVolJ / volOfMicrogel[j];
                  if (nonOverlapFracJ < 0.0) nonOverlapFracJ = 0.0;
                  double newSWRj = aj * Math.cbrt(nonOverlapFracJ);
                  double swrj2 = newSWRj * newSWRj;
                  double swrj3 = swrj2 * newSWRj;
                  double inv_swrj3 = 1.0 / swrj3;

                  mixF     = nMon * ((swrj3 - 1.0) * Math.log(1.0 - inv_swrj3) + chi * (1.0 - inv_swrj3));
                  elasticF = 1.5 * nMon * xLinkFrac * (swrj2 - Math.log(newSWRj) - 1.0);
                  floryFRAfterMove[j] = mixF + elasticF;

                  dFloryForJ += (floryFRAfterMove[j] - oldFreeEnergyJ);
               } else {
                  capVolAfterTrialMove[i][j] = 0.0;
                  capVolAfterTrialMove[j][i] = 0.0;

                  // FR energy for j after the move (no overlap contribution from i)
                  double newCapVolJ = capVolSumBeforeMoveArray[j]; // unchanged minus/plus zeros
                  double nonOverlapFracJ = 1.0 - newCapVolJ / volOfMicrogel[j];
                  if (nonOverlapFracJ < 0.0) nonOverlapFracJ = 0.0;
                  double newSWRj = aj * Math.cbrt(nonOverlapFracJ);
                  double swrj2 = newSWRj * newSWRj;
                  double swrj3 = swrj2 * newSWRj;
                  double inv_swrj3 = 1.0 / swrj3;

                  mixF     = nMon * ((swrj3 - 1.0) * Math.log(1.0 - inv_swrj3) + chi * (1.0 - inv_swrj3));
                  elasticF = 1.5 * nMon * xLinkFrac * (swrj2 - Math.log(newSWRj) - 1.0);
                  floryFRAfterMove[j] = mixF + elasticF;

                  dFloryForJ += (floryFRAfterMove[j] - oldFreeEnergyJ);
               }
         }

         // FR free energy for microgel i AFTER the move
         double nonOverlapFracI = 1.0 - (capVolSumAfterMoveI / volOfMicrogel[i]);
         if (nonOverlapFracI < 0.0) nonOverlapFracI = 0.0;
         newSWR = ai * Math.cbrt(nonOverlapFracI);
         double swri2 = newSWR * newSWR;
         double swri3 = swri2 * newSWR;
         double inv_swri3 = 1.0 / swri3;

         mixF     = nMon * ((swri3 - 1.0) * Math.log(1.0 - inv_swri3) + chi * (1.0 - inv_swri3));
         elasticF = 1.5 * nMon * xLinkFrac * (swri2 - Math.log(newSWR) - 1.0);
         floryFRAfterMove[i] = mixF + elasticF;

         double dFloryForI = floryFRAfterMove[i] - oldFreeEnergyI; // change in FR energy for i
         dFlory = dFloryForI + dFloryForJ;

         // Metropolis criterion
         de = lambda * (pairEnergySum - energy[i]) + (1.0 - lambda) * (pairPotentialCorrection) + dFlory;

         if (Math.exp(-de) < Math.random()) {
               // reject move: revert everything that changed for i
               x[i] -= dxtrial;
               y[i] -= dytrial;
               z[i] -= dztrial;
               a[i] -= datrial;

         } else {
               // accept move
               dxOverN += dxtrial / (double) N;
               dyOverN += dytrial / (double) N;
               dzOverN += dztrial / (double) N;

               energy[i] = pairEnergySum;

               // update neighbors’ pair energies & cap volumes cache
               for (int neighborIdx = 0; neighborIdx < 12; neighborIdx++) {
                  j = nearestNeighbors[i][neighborIdx];
                  energy[j] += newPairEnergy[i][j] - pairEnergy[i][j];
                  pairEnergy[i][j] = newPairEnergy[i][j];
                  pairEnergy[j][i] = newPairEnergy[i][j];

                  capVolBeforeTrialMove[i][j] = capVolAfterTrialMove[i][j];
                  capVolBeforeTrialMove[j][i] = capVolAfterTrialMove[j][i];
               }

               // update prepass arrays for i
               capVolSumBeforeMoveArray[i] = capVolSumAfterMoveI;
               floryFRBeforeMove[i] = floryFRAfterMove[i];
               for (int neighborIdx = 0; neighborIdx < 12; neighborIdx++) {
                  j = nearestNeighbors[i][neighborIdx];
                  floryFRBeforeMove[j] = floryFRAfterMove[j];
               }
         }
      }

      calculateTotalEnergy(lambda); // new total energy
   }

   public void calculateTotalEnergy(double lambda) {
      // reset totals
      totalEnergy = 0.0;         // when lambda = 1
      totalVirial = 0.0;
      totalPairEnergy = 0.0;     // when lambda = 1
      totalFreeEnergy = 0.0;
      springEnergySum = 0.0;     // total spring energy of the system
      squaredDisplacementSum = 0.0;

      // ---------- constants that don't depend on i/j ----------
      final double twoPointFive = 2.5;
      final double onePointFive_nMon_x = 1.5 * nMon * xLinkFrac;

      // reservoir terms are invariant across i
      final double resSR   = reservoirSR;
      final double resSR2  = resSR * resSR;
      final double resSR3  = resSR2 * resSR;
      final double invRes3 = 1.0 / resSR3;
      final double mixFRSR_local =
      nMon * ((resSR3 - 1.0) * Math.log(1.0 - invRes3) + chi * (1.0 - invRes3));
      final double elasticFRSR_local = 1.5 * nChains * (resSR2 - Math.log(resSR) - 1.0);
      final double totalFRSR_local = mixFRSR_local + elasticFRSR_local;

      double virialSum;

      for (int i = 0; i < N; ++i) {
         pairEnergySum = 0.0;   // for the lambda = 1 system
         virialSum = 0.0;

         // squared displacement relative to initial pos (with COM shift)
         double dxi = x[i] - x0[i] - dxOverN;
         double dyi = y[i] - y0[i] - dyOverN;
         double dzi = z[i] - z0[i] - dzOverN;
         double sd  = dxi*dxi + dyi*dyi + dzi*dzi;

         squaredDisplacement = sd;
         squaredDisplacementSum += sd;

         // single-particle spring energy
         double sE = springConstant * sd;
         springEnergySum += sE;

         // cache radius powers for i
         double ai  = a[i];
         double ai2 = ai * ai;
         double ai3 = ai2 * ai;
         double sqrt_ai = Math.sqrt(ai);

         // pair loop over 12 nn
         for (int neighborIdx = 0; neighborIdx < 12; neighborIdx++) {
               int j = nearestNeighbors[i][neighborIdx];

               double xij = PBC.separation(x[i] - x[j], side);
               double yij = PBC.separation(y[i] - y[j], side);
               double zij = PBC.separation(z[i] - z[j], side);

               double r2 = xij*xij + yij*yij + zij*zij;

               double aj  = a[j];
               double aj2 = aj * aj;
               double aj3 = aj2 * aj;
               double sqrt_aj = Math.sqrt(aj);

               double sigma  = ai + aj;
               double sigma2 = sigma * sigma;

               if (r2 < sigma2) {
                  double r = Math.sqrt(r2);

                  // Hertz pair potential amplitude (lambda = 1 system)
                  double Bamp = scale * Young * nChains * sigma2 * (sqrt_ai * sqrt_aj) / (ai3 + aj3); 

                  // (1 - r/sigma) and its 1.5 power without pow()
                  double oneMinus = 1.0 - r / sigma;
                  if (oneMinus < 0.0) oneMinus = 0.0; // guard for fp noise
                  double sqrtOneMinus = Math.sqrt(oneMinus);
                  double derivPart = Bamp * oneMinus * sqrtOneMinus; // (1-r/sigma)^{1.5}

                  // energy (1-r/sigma)^{2.5}
                  double HertzEnergy = derivPart * oneMinus;
                  pairEnergySum += HertzEnergy;
                  pairEnergy[i][j] = HertzEnergy; // store pair energy

                  // virial: x·f + y·f + z·f = (2.5/sigma)*derivPart/r * r2  = (2.5/sigma)*derivPart*r
                  virialSum += (twoPointFive / sigma) * derivPart * r;
               }
         }

         // Flory–Rehner single-particle free energy
         double inv_ai3 = 1.0 / ai3;
         double mixF = nMon * ((ai3 - 1.0) * Math.log(1.0 - inv_ai3) + chi * (1.0 - inv_ai3));
         double elasticF = onePointFive_nMon_x * (ai2 - Math.log(ai) - 1.0);
         double totalF = (elasticF + mixF) - totalFRSR_local;

         energy[i] = pairEnergySum; // current particle’s pair energy (lambda system)

         totalPairEnergy += pairEnergySum; // when lambda = 1
         totalFreeEnergy += totalF;
         totalVirial += virialSum;
      }

      totalPairEnergy *= 0.5; // correct for double counting pairs

      if (steps > delay) {
         if ((steps - delay) % snapshotInterval == 0) { // sample at intervals
               numberOfConfigurations++;
               totalVirial *= 0.5; // correct for double counting pairs
               totalEnergy = totalPairEnergy + totalFreeEnergy;

               pairEnergyAccumulator += totalPairEnergy; // when lambda = 1
               energyAccumulator += totalEnergy;         // running totals
               freeEnergyAccumulator += totalFreeEnergy;
               virialAccumulator += totalVirial;

               squaredDisplacementAccumulator += squaredDisplacementSum; // ⟨r^2⟩
               springEnergyAccumulator += springEnergySum;               // ⟨U_spring⟩

               // Accumulate volume fraction
               volFracAccumulator += calculateVolumeFraction();
         }
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

   // public double meanboltzmannFactor(){ // mean boltzmannFactor
   //    return boltzmannFactorAccumulator/numberOfConfigurations;
   // }

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
         if (bin >= 0 && bin < sizeDist.length){
            sizeDist[bin]++;
         }
         else{
            System.err.println("Error: bin index out of bounds: " + bin);
         }
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
         //volFrac += microgelVol/totalVol;

         return microgelVol/totalVol;
   }

    /* reservoir swelling ratio as root of Flory-Rehner pressure */
      public double reservoirSwellingRatio(double nMon, double nChains, double chi){
         Function f = new FloryRehnerPressure(nMon, nChains, chi);
         double xleft = 1.1;
         double xright = 15;
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
