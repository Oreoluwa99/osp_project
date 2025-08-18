/*
 * The major reason for this code is to account for the non-local faceting of the microgels, 
 * where faceting is assumed to redistribute polymer to a smaller volume (thus reducing the swelling ratio), 
 * rather than “crushing” the polymer in the cap regions.
 * 
 * ✅ NEW in this version:
 * This version now explicitly accounts for **only the 12 nearest neighbors** 
 * in the **FCC crystal structure** during Monte Carlo updates and energy calculations,
 * improving computational performance and better representing the local environment.
 *
 * The step() method is where the change took place.
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
 import org.opensourcephysics.numerics.Function;
 
 /**
  * 
  * @authors Oreoluwa Alade and Alan Denton
  * @version 4.6 - 08-24
  * @version 4.7 - 30.01.2025
  * @version 4.8 - 03.02.2025 (edited by Alan Denton)
  * 
  */
 
 public class HertzSpheresNonLocalFacetV2 {
    
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
     // public Random random = new Random(); // Random number generator for Monte Carlo moves
     public double volOfMicrogel[]; //volume of each microgel
     // Arrays to hold Flory-Rehner free energy for microgels after a move (for J and I)
     public double floryFRJAfterMove[], floryFRJBeforeMove[];
     public double floryFRIBeforeMove, floryFRIAfterMove;
     public double newNMon[], cumulativeCapVol[], capVol[]; //new number of monomers, cap volume and cummulative cap volume for each microgel
      // Alan: capVol[] is never used!
     public double floryFRAfterMove[], floryFRBeforeMove[]; //Flory-Rehner free energy for microgels before and after a move
     public double capVolJ[], capVolJBeforeMove[], capVolJAfterMove[]; // Alan: these arrays are never used!
     public double capVolSumBeforeMoveArray[], capVolSumAfterMoveArray[]; // Alan: these arrays are never used!
     public double capVolAfterTrialMove[][], capVolBeforeTrialMove[][]; // stores the cap volume for each j microgel before and after a move
     public double capVolumesBeforeMove[];
     public double capVolumesAfterMove[]; // Alan: need this array to store cap volumes after trial move
     public double totalVolFrac;
     public double nonOverlappingVol, volumeFraction;
     public int[][] neighborList;

    
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
       capVolJ = new double[N]; // Alan: never used!
       capVolJBeforeMove = new double[N]; // Alan: never used!
       capVolJAfterMove = new double[N]; // Alan: never used!
       capVol = new double[N]; // Alan: capVol[] is never used!
       newNMon = new double[N];
       cumulativeCapVol = new double[N];
       capVolAfterTrialMove = new double[N][N];
       capVolBeforeTrialMove = new double[N][N];
       volOfMicrogel = new double[N];
       floryFRJAfterMove = new double[N];
       floryFRJBeforeMove = new double[N];
       capVolumesBeforeMove = new double[N];
       capVolumesAfterMove = new double[N]; // Alan: need this array to store cap volumes after trial move
 
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

       neighborList = new int[N][12];

       // For generating the 12 nearest neighbors
       for (int i = 0; i < N; i++) {
            double[] dist = new double[N];
            Integer[] indices = new Integer[N];
            for (int j = 0; j < N; j++) {
               indices[j] = j;
               if (i == j) {
                     dist[j] = Double.MAX_VALUE; // exclude self
               } else {
                     double dx = PBC.separation(x[i] - x[j], side);
                     double dy = PBC.separation(y[i] - y[j], side);
                     double dz = PBC.separation(z[i] - z[j], side);
                     dist[j] = dx*dx + dy*dy + dz*dz;
               }
            }  
            java.util.Arrays.sort(indices, java.util.Comparator.comparingDouble(o -> dist[o]));
            for (int k = 0; k < 12; k++) {
               neighborList[i][k] = indices[k];
            }
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
 
      //  // set the seed of the random number generators
      //  random.setSeed(12345);
 
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
       double dxtrial, dytrial, dztrial, datrial;
       double HertzEnergy;
       double xij, yij, zij, de, r, r2, sigma;
       double mixF, elasticF, totalF;
       double capVolForI;
       double da2;
       double totalFRBeforeMove, totalFRAfterMove;
       int i,j;
       double dFlory;
       double capVolForJ;
       double floryRehnerFR;
       boolean isFirstIteration, isSecondIteration;
       boolean overlap=false;
       double newNMon;
       double newSWR;
       double capVolSumBeforeMoveJ, capVolSumAfterMoveI, capVolSumBeforeMoveI, capVolSumAfterMoveJ;
       // Alan: capVolSumAfterMoveJ is never used!
  
       // Compute total cap volumes and FR free energy for all microgels before the move (the first time through)
       for (i = 0; i < N; i++) {
          capVolSumBeforeMoveI = 0;
          //capVolSumBeforeMoveArray[i] = 0; // Alan: this array is never used!
          volOfMicrogel[i] = (4.0 / 3.0) * Math.PI * Math.pow(a[i], 3);

          // consider interactions with other microgels
          for (int k = 0; k < 12; k++) {
             j = neighborList[i][k];
             if (j == i) continue;
                xij = PBC.separation(x[i] - x[j], side);
                yij = PBC.separation(y[i] - y[j], side);
                zij = PBC.separation(z[i] - z[j], side);

                r2 = xij * xij + yij * yij + zij * zij;
                sigma = a[i] + a[j];
 
                // if the microgels are overlapping
                if (r2 < sigma * sigma) {
                   r = Math.sqrt(r2);
                   double hi = (a[j]-a[i]+r)*(a[j]+a[i]-r)/(2.0*r);
                   double hj = (a[i]-a[j]+r)*(a[i]+a[j]-r)/(2.0*r);
                   //System.out.println("hi: "+hi);
                   // the individual caps of the microgels
                   capVolForI = (Math.PI*Math.pow(hi, 2)*(3.0*a[i]-hi))/3.0;
                   if (capVolForI < -1e-10) { System.out.println("STOP: capVolForI < 0 !!!");} // Alan: check for negative cap volumes
                   //System.out.println("capVolI: "+capVolForI);
                   capVolForJ = (Math.PI*Math.pow(hj, 2)*(3.0*a[j]-hj))/3.0; 
                   capVolSumBeforeMoveI += capVolForI;
                   //capVolJBeforeMove[j] = capVolForJ; // Alan: this array is never used!
                   capVolBeforeTrialMove[i][j] = capVolForJ;
                   capVolBeforeTrialMove[j][i] = capVolForI; 
                }
                else{ //if they are not overlapping
                   capVolBeforeTrialMove[i][j] = 0;
                   capVolBeforeTrialMove[j][i] = 0;
                }

          }

	       //System.out.println("capVolSumBeforeMoveI = " + capVolSumBeforeMoveI);
          capVolSumBeforeMoveArray[i] = capVolSumBeforeMoveI; // Alan: this array is never used!
          capVolumesBeforeMove[i] = capVolSumBeforeMoveI; // dummy array
          // Alan: Why do you call this a "dummy array"?
          if (capVolSumBeforeMoveI < -1e-10) { System.out.println("STOP: capVolSumBeforeMoveI < 0 !!!");} // Alan: check for negative cap volumes
          // new swelling ratio
          newSWR = a[i] * Math.pow(((volOfMicrogel[i] - capVolSumBeforeMoveI) / volOfMicrogel[i]), 1.0 / 3.0); // original swelling ratio expression
          mixF = nMon*((newSWR*newSWR*newSWR-1)*Math.log(1-1/(newSWR*newSWR*newSWR))+chi*(1-1/(newSWR*newSWR*newSWR)));

          elasticF = 1.5 * nMon * xLinkFrac * (newSWR * newSWR - Math.log(newSWR) - 1);
          // the Flory-Rehner free energy before the move for each microgel
          floryFRBeforeMove[i] = mixF + elasticF;

       }

      for (i = 0; i < N; i++) {
          
         // Update the old free energy for the current state before making the trial move
         double oldFreeEnergyI = floryFRBeforeMove[i];
	      //System.out.println("oldFreeEnergyI = " + oldFreeEnergyI);

         dxtrial = tolerance*2.*(Math.random()-0.5);
         dytrial = tolerance*2.*(Math.random()-0.5);
         dztrial = tolerance*2.*(Math.random()-0.5);
         datrial = atolerance*2.*(Math.random()-0.5);
         
         // Euclidean distance of trial displacements
         trialDisplacementDistance = (dxtrial*dxtrial)+(dytrial*dytrial)+(dztrial*dztrial);

         x[i] += dxtrial; // trial displacement 
         y[i] += dytrial;                       
         z[i] += dztrial;
         a[i] += datrial; // trial radius change
  
         //the Euclidean distance between the current position and the initial position of each particle and the accumulated shift in the center of mass
         double dxi = x[i]-x0[i]-dxOverN;
         double dyi = y[i]-y0[i]-dyOverN;
         double dzi = z[i]-z0[i]-dzOverN;

         volOfMicrogel[i] = (4.0 / 3.0) * Math.PI * Math.pow(a[i], 3);
 
         pairEnergySum = 0;
         capVolSumAfterMoveI = 0;
         //capVolSumAfterMoveJ = 0; // this variable is never used!
         double dFloryForJ = 0; // accumulates the change in Flory-Rehner free energy for the j microgels
         double oldFloryForJ = 0; // Alan: need to zero these sums
         double newFloryForJ = 0; // Alan

         pairPotentialCorrection = 2*springConstant*Math.sqrt((dxtrial*dxi+dytrial*dyi+dztrial*dzi)*(dxtrial*dxi+dytrial*dyi+dztrial*dzi))+(1-1/(double)N)*trialDisplacementDistance; // change in spring energy

         for(int k = 0; k < 12; k++) {
            j = neighborList[i][k];
 
             volOfMicrogel[j] = (4.0 / 3.0) * Math.PI * Math.pow(a[j], 3);
 
             if (j == i) continue;
                // Save current free energy as old free energy for microgel j
                double oldFreeEnergyJ = floryFRBeforeMove[j];
                
                xij = PBC.separation(x[i] - x[j], side);
                yij = PBC.separation(y[i] - y[j], side);
                zij = PBC.separation(z[i] - z[j], side);
  
                r2 = xij * xij + yij * yij + zij * zij;
                sigma = a[i] + a[j];
                
                // if the microgels are overlapping
                if (r2 < sigma * sigma) {
                   r = Math.sqrt(r2);
                   B = scale * Young * nChains * Math.pow(sigma, 2.0) * Math.sqrt(a[i] * a[j]) / (Math.pow(a[i], 3.0) + Math.pow(a[j], 3.0));
                   HertzEnergy = B * Math.pow(1 - r / sigma, 2.5);
                   newPairEnergy[i][j] = HertzEnergy;
                   pairEnergySum += HertzEnergy;

                   //System.out.println("x: " + x[i] + ", y: " + y[i] + ", z: " + z[i] + ", ai: " + a[i] + ", aj: " + a[j] + "d: "+r);
  
                   double hi = (a[j]-a[i]+r)*(a[j]+a[i]-r)/(2.0*r);

                   //System.out.println("hi: "+hi);

                   double hj = (a[i]-a[j]+r)*(a[i]+a[j]-r)/(2.0*r);
                   // compute cap volume
                   capVolForJ = (Math.PI * Math.pow(hj, 2) * (3.0 * a[j] - hj)) / 3.0;
                   capVolForI = (Math.PI * Math.pow(hi, 2) * (3.0 * a[i] - hi)) / 3.0;

                   //System.out.println("capVolI: "+ capVolForI); // output the cap volume for microgel i
                   capVolSumAfterMoveI += capVolForI;

                   capVolAfterTrialMove[i][j] = capVolForJ;
                   capVolAfterTrialMove[j][i] = capVolForI; 
 
                   // Calculate FR free energy for microgel j after the move
                   //double newCapVolJ = capVolSumBeforeMoveArray[j] + capVolAfterTrialMove[i][j] - capVolBeforeTrialMove[i][j];
                   double newCapVolJ = capVolumesBeforeMove[j] + capVolAfterTrialMove[i][j] - capVolBeforeTrialMove[i][j]; 

                   capVolumesAfterMove[j] = newCapVolJ; // Alan: store newCapVolJ in this array and later update capVolumesBeforeMove[j] if the trial move is accepted

                   newSWR = a[j] * Math.pow(((volOfMicrogel[j] - newCapVolJ) / volOfMicrogel[j]), 1.0/3.0);
                   mixF = nMon * ((newSWR * newSWR * newSWR - 1) * Math.log(1 - 1 / (newSWR * newSWR * newSWR)) + chi * (1 - 1 / (newSWR * newSWR * newSWR)));
                   elasticF = 1.5 * nMon * xLinkFrac * (newSWR * newSWR - Math.log(newSWR) - 1);

                   floryFRAfterMove[j] = mixF + elasticF;
                   
                   // Calculate change in energy for microgel j
                   dFloryForJ += floryFRAfterMove[j] - oldFreeEnergyJ;
 
                }
                else{//if they are not overlapping
                   newPairEnergy[i][j] = 0; newPairEnergy[j][i] = 0; // Alan: need to set pair energy to zero if no overlap
                   capVolAfterTrialMove[i][j] = 0;
                   capVolAfterTrialMove[j][i] = 0;
 
                   // Calculate FR free energy for microgel j after the move
                   double newCapVolJ = capVolumesBeforeMove[j] - capVolBeforeTrialMove[i][j]; // Alan
                   capVolumesAfterMove[j] = newCapVolJ; // Alan: store newCapVolJ in this array and later update capVolumesBeforeMove[j] if the trial move is accepted

                  //  double newCapVolJ = capVolumesBeforeMove[j];
                  // if (capVolBeforeTrialMove[i][j] > 0) {
                  //    newCapVolJ -= capVolBeforeTrialMove[i][j];
                  // }

                   if(newCapVolJ < -1.e-10) { // Alan: check for negative cap volumes
                      System.out.println("STOP: newCapVolJ < 0!!!"); 
                      System.out.println("newCapVolJ = " + newCapVolJ + ", capVolumesBeforeMove[j] = " + capVolumesBeforeMove[j] + ", capVolAfterTrialMove[i][j] = " + capVolAfterTrialMove[i][j] + ", capVolBeforeTrialMove[i][j] = " + capVolBeforeTrialMove[i][j]);
                      try {
                         System.in.read(); // pause and wait for key to be pressed
                      } catch (IOException e) {
                      e.printStackTrace();
                      }  
                   }  

                   newSWR = a[j] * Math.pow(((volOfMicrogel[j] - newCapVolJ) / volOfMicrogel[j]), 1.0 / 3.0);
                   mixF = nMon * ((newSWR * newSWR * newSWR - 1) * Math.log(1 - 1 / (newSWR * newSWR * newSWR)) + chi * (1 - 1 / (newSWR * newSWR * newSWR)));
                   elasticF = 1.5 * nMon * xLinkFrac * (newSWR * newSWR - Math.log(newSWR) - 1);
                   floryFRAfterMove[j] = mixF + elasticF;
                   
                   // Calculate change in energy for microgel j
                   dFloryForJ += floryFRAfterMove[j] - oldFreeEnergyJ;
                   oldFloryForJ += oldFreeEnergyJ; // Alan
                   newFloryForJ += floryFRAfterMove[j]; // Alan
                }
          }
          // Update cap volumes after the trial move
          //capVolSumAfterMoveArray[i] = capVolSumAfterMoveI; // Alan: capVolSumAfterMoveArray is never used
          capVolumesAfterMove[i] = capVolSumAfterMoveI; // Alan: need to store this cap volume for updates
          
          // Flory-Rehner free energy for microgel i after the move
          newSWR = a[i] * Math.pow(((volOfMicrogel[i] - capVolSumAfterMoveI)/volOfMicrogel[i]), 1.0/3.0);
          mixF = nMon*((newSWR*newSWR*newSWR-1)*Math.log(1-1/(newSWR*newSWR*newSWR))+chi*(1-1/(newSWR*newSWR*newSWR)));
          elasticF = 1.5*nMon*xLinkFrac*(newSWR*newSWR-Math.log(newSWR) - 1);
          floryFRAfterMove[i] = mixF + elasticF;

          // Calculate change in energy for microgel i
          double dFloryForI = floryFRAfterMove[i] - oldFreeEnergyI; // Use updated old free energy

          double oldFlory = oldFreeEnergyI + oldFloryForJ; // Alan
          double newFlory = floryFRAfterMove[i] + newFloryForJ; // Alan
 
          // Calculate total change in FR energy for the system
          dFlory = dFloryForI + dFloryForJ;
   
          de = lambda*(pairEnergySum-energy[i])+(1-lambda)*(pairPotentialCorrection)+dFlory; // change in total energy due to trial move

         if(Math.exp(-de) < Math.random()){ // Metropolis algorithm
            //System.out.println("Move rejected!");
            x[i] -= dxtrial; // reject move
            y[i] -= dytrial;
            z[i] -= dztrial;
	         a[i] -= datrial;
         } 
          
         else { //accept move and update energy
             //System.out.println("Move accepted!");
             dxOverN += dxtrial / (double) N;
             dyOverN += dytrial / (double) N;
             dzOverN += dztrial / (double) N;
  
             energy[i] = pairEnergySum;
  
             for (int k = 0; k < 12; k++) {
               j = neighborList[i][k];
               if (j == i) continue;
                   energy[j] += newPairEnergy[i][j] - pairEnergy[i][j];
                   pairEnergy[i][j] = newPairEnergy[i][j];
                   pairEnergy[j][i] = newPairEnergy[i][j];
                   // Update capVolBeforeMove with the new accepted values
                   capVolBeforeTrialMove[i][j] = capVolAfterTrialMove[i][j]; // Alan uncommented 
                   capVolBeforeTrialMove[j][i] = capVolAfterTrialMove[j][i]; // Alan uncommented
                   // Update old free energy for microgel i before the next trial move
                   floryFRBeforeMove[j] = floryFRAfterMove[j]; // Alan: moved here, so that j!=i
                   capVolumesBeforeMove[j] = capVolumesAfterMove[j]; // Alan: update cap volume before next trial move
                   if (capVolumesBeforeMove[j] < -1e-10) { System.out.println("STOP: capVolumesBeforeMove[j] < 0 !!!");} // Alan: check for negative cap volumes
             }
             // Update capVolSumBeforeMoveArray with the new value for microgel i
             capVolumesBeforeMove[i] = capVolSumAfterMoveI; // Alan 
             if (capVolumesBeforeMove[i] < -1e-10) { System.out.println("STOP: capVolumesBeforeMove[i] < 0 !!!");} // Alan: check for negative cap volumes
 
          }
       }
        calculateTotalEnergy(lambda); // new total energy
    }

    public void calculateTotalEnergy(double lambda) {
       totalEnergy = 0;  //when lambda = 1
       totalVirial = 0;
       totalPairEnergy = 0;  //when lambda = 1
       totalFreeEnergy = 0;
       springEnergySum = 0; // total spring energy of the system
       squaredDisplacementSum = 0;
       double virialSum;
       double xij, yij, zij, r, r2, sigma, dxi, dyi, dzi;
       double derivPart, HertzEnergy;
       double mixF, elasticF, totalF;
       double fOverR, fx, fy, fz;
       for(int i = 0; i < N; i++) {
          pairEnergySum = 0; //for the lamba = 1 system
          virialSum = 0;
          
          // the Euclidean distance between the current position and the initial position of each particle
          dxi = x[i]-x0[i]-dxOverN;
          dyi = y[i]-y0[i]-dyOverN;
          dzi = z[i]-z0[i]-dzOverN;
         
          squaredDisplacement = Math.pow(dxi, 2)+Math.pow(dyi, 2)+Math.pow(dzi, 2); // the square displacement 
          squaredDisplacementSum+=squaredDisplacement;
          springEnergy = springConstant*squaredDisplacement; //since it is a single particle (i) component
          springEnergySum += springEnergy; 
          for(int k = 0; k < 12; k++) {
            int j = neighborList[i][k];
            if (j == i) continue;
                xij = PBC.separation(x[i]-x[j], side);
                yij = PBC.separation(y[i]-y[j], side);
                zij = PBC.separation(z[i]-z[j], side);
                r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
 
                sigma = a[i]+ a[j];
                if (r2 < sigma*sigma){
                   r = Math.sqrt(r2);
                   // Hertz pair potential amplitude for lambda = 1 system
                   B = scale*Young*nChains*Math.pow(sigma, 2.)*Math.sqrt(a[i]*a[j])/(Math.pow(a[i], 3.)+Math.pow(a[j], 3.)); 
                   derivPart = B*Math.pow(1-r/sigma, 1.5);
                   HertzEnergy = derivPart*(1-r/sigma);
                   pairEnergySum += HertzEnergy;
                   pairEnergy[i][j] = HertzEnergy; // pair energy of particles i and j
                   fOverR = (2.5/sigma)*derivPart/r;
                   fx = fOverR*xij; // pair force in x-direction
                   fy = fOverR*yij; // pair force in y-direction
                   fz = fOverR*zij; // pair force in z-direction
                   virialSum += xij*fx + yij*fy + zij*fz;
                }
                else { // Alan: set nonoverlapping pair energies to zero
                   pairEnergy[i][j] = 0; pairEnergy[j][i] = 0;// Alan
                }   
          }
 
          // Flory-Rehner single-particle free energy (associated with swelling)
          mixF = nMon*((a[i]*a[i]*a[i]-1)*Math.log(1-1/a[i]/a[i]/a[i])+chi*(1-1/a[i]/a[i]/a[i]));
          elasticF = 1.5*nMon*xLinkFrac*(a[i]*a[i]-Math.log(a[i])-1);
          
          totalF = elasticF + mixF; // Flory-Rehner free energy after the trial move
 
          energy[i] = pairEnergySum; // total energy of particle i for the lambda system // because the change in energy is computed in the step method and the change has been used in the step method
 
          mixFRSR = nMon*((reservoirSR*reservoirSR*reservoirSR-1)*Math.log(1-1/reservoirSR/reservoirSR/reservoirSR)+chi*(1-1/reservoirSR/reservoirSR/reservoirSR));
          elasticFRSR = 1.5*nChains*(reservoirSR*reservoirSR-Math.log(reservoirSR)-1);
 
          totalFRSR = mixFRSR + elasticFRSR; // total free energy associated with the reservoir state
          totalF = totalF - totalFRSR;
       
          totalPairEnergy += pairEnergySum; //when lambda = 1
          totalFreeEnergy += totalF;
          totalVirial += virialSum;
       }
          totalPairEnergy *= 0.5; // correct for double counting pairs
    
          if (steps > delay){
             if ((steps-delay)%snapshotInterval==0){ // include configurations only after even number of intervals
                numberOfConfigurations++;
                totalVirial *= 0.5; // correct for double counting pairs
                totalEnergy = totalPairEnergy + totalFreeEnergy;
                pairEnergyAccumulator += totalPairEnergy; // when lambda = 1
                energyAccumulator += totalEnergy; // running totals
                freeEnergyAccumulator += totalFreeEnergy;
                virialAccumulator += totalVirial;
                squaredDisplacementAccumulator += squaredDisplacementSum; // accumulates the mean sqaure displacement
                springEnergyAccumulator += springEnergySum; // accumulates all the spring energies
                // Accumulate the volume fraction
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
 
    public double meanSpringEnergy(){ // mean springEnergyAccumulator
      return springEnergyAccumulator/N/numberOfConfigurations;
   }

    public double meanSquareDisplacement(){
       return squaredDisplacementAccumulator/N/numberOfConfigurations; // quantity <r>^2
    }
 
    public double meanFreeEnergy() { // mean free energy per particle
       return freeEnergyAccumulator/N/numberOfConfigurations; // quantity <F>/N
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
       //HertzSpheresNonLocalFacet.FloryRehnerPressure f = new FloryRehnerPressure(nMon, nChains, chi);
       HertzSpheresNonLocalFacetV2.FloryRehnerPressure f = new FloryRehnerPressure(nMon, nChains, chi);
       double xleft = 1.1;
       double xright = 15.;
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
