   /*
   * HertzianSpheresModel.java
   *
   * Monte Carlo model of fixed-size Hertzian spheres.
   *
   * This class provides the common physical model used for both:
   *
   *   1. Fluid-phase equation-of-state and free-energy calculations.
   *   2. Solid-phase Frenkel-Ladd free-energy calculations.
   *
   * Pair potential:
   *
   *                 {
   *                 | betaEpsilon * (1 - r/sigma)^(5/2),  r < sigma
   *   beta*u(r)  =  |
   *                 | 0,                                  r >= sigma
   *                 }
   *
   * The particle diameter sigma and interaction strength betaEpsilon are fixed
   * and are identical in the fluid and solid simulations.
   *
   * This version was derived from HertzSpheresNonLocalFacet.java. The following
   * compressible-microgel features have been removed for the fixed-size
   * Hertzian-sphere benchmark:
   *
   *   - particle-radius fluctuations;
   *   - Flory-Rehner free energy;
   *   - facet and spherical-cap calculations;
   *   - cross-link-dependent Hertz amplitude;
   *   - Young's-modulus-dependent Hertz amplitude;
   *   - interpenetration corrections.
   *
   * Fluid Monte Carlo moves are accepted using only the change in Hertzian
   * pair energy. Einstein springs and centre-of-mass corrections are used only
   * during the Frenkel-Ladd solid free-energy calculation.
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


   public class HertzianSpheresModel {
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
      public Random random = new Random(12345); // Random number generator for Monte Carlo moves
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
      public double modifiedSwellingRatioAccumulator;
      public double sigma = 1.0;
      public double betaEpsilon = 12000.0;

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
         a = new double[N];

         // // Initialize radii with the fully swollen reservoirSR
         // for (int i = 0; i < N; i++) {
         // a[i] = reservoirSR; // initial particle radii (fully swollen)
         // }

         // Fixed particle radius in units where sigma = 1.
         for (int i = 0; i < N; i++) {
            a[i] = sigma / 2.0;
         }

         // density is rho* = N*sigma^3/V.
         side = sigma * Math.cbrt(N / density);
         totalVol = side * side * side;

         energy = new double[N];
         pairEnergy = new double[N][N];
         newPairEnergy = new double[N][N];

         steps = 0;
         energyAccumulator = 0.0;
         pairEnergyAccumulator = 0.0;
         freeEnergyAccumulator = 0.0;
         virialAccumulator = 0.0;
         springEnergyAccumulator = 0.0;
         numberOfConfigurations = 0.0;
         squaredDisplacementAccumulator = 0.0;
         volFracAccumulator = 0.0;

         dxOverN = 0.0;
         dyOverN = 0.0;
         dzOverN = 0.0;
         volFrac = 0; // counter for system volume fraction

         numberBins = (int) (maxRadius/sizeBinWidth);

         sizeDist = new double[numberBins];

         for(int i=0; i<numberBins; i++){ // initialize size histogram
            sizeDist[i] = 0;
         }

         // initialize positions
         if (configuration.toUpperCase().equals("SC")) {
            setSCpositions();
         }
         if (configuration.toUpperCase().equals("FCC")) {
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
                        x[i] = PBC.position(ix * d/2. + (random.nextDouble()-0.5) * d, side);
                        y[i] = PBC.position(iy * d/2. + (random.nextDouble()-0.5) * d, side);
                        z[i] = PBC.position(iz * d/2. + (random.nextDouble()-0.5) * d, side);
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
                        x[i] = PBC.position(ix * d/2. + (random.nextDouble()-0.5) * d, side);
                        y[i] = PBC.position(iy * d/2. + (random.nextDouble()-0.5) * d, side);
                        z[i] = PBC.position(iz * d/2. + (random.nextDouble()-0.5) * d, side);
                        i++;
                     }
                  }
               }
            }
         }
         calculateTotalEnergy(lambda); // initial energy
      }

      /**
       * Performs one Monte Carlo sweep of the Hertzian fluid.
      */
      public void step() {
         steps++;

         for (int i = 0; i < N; i++) {
            double dxtrial =
                     tolerance * 2.0 * (random.nextDouble() - 0.5);
            double dytrial =
                     tolerance * 2.0 * (random.nextDouble() - 0.5);
            double dztrial =
                     tolerance * 2.0 * (random.nextDouble() - 0.5);

            // Apply the trial displacement.
            x[i] += dxtrial;
            y[i] += dytrial;
            z[i] += dztrial;

            double newEnergyI = 0.0;

            // Calculate the new Hertzian energy involving particle i.
            for (int j = 0; j < N; j++) {
                  if (j == i) {
                     continue;
                  }

                  double xij = PBC.separation(x[i] - x[j], side);
                  double yij = PBC.separation(y[i] - y[j], side);
                  double zij = PBC.separation(z[i] - z[j], side);

                  double r = Math.sqrt(
                        xij*xij + yij*yij + zij*zij
                  );

                  // This also sets the pair energy to zero when r >= sigma.
                  newPairEnergy[i][j] = hertzEnergy(r);
                  newEnergyI += newPairEnergy[i][j];
            }

            double deltaEnergy = newEnergyI - energy[i];

            boolean accept =
                     deltaEnergy <= 0.0
                     || random.nextDouble() < Math.exp(-deltaEnergy);

            if (accept) {
                  energy[i] = newEnergyI;

                  for (int j = 0; j < N; j++) {
                     if (j == i) {
                        continue;
                     }

                     energy[j] +=
                              newPairEnergy[i][j] - pairEnergy[i][j];

                     pairEnergy[i][j] = newPairEnergy[i][j];
                     pairEnergy[j][i] = newPairEnergy[i][j];
                  }
            } else {
                  // Restore the old position.
                  x[i] -= dxtrial;
                  y[i] -= dytrial;
                  z[i] -= dztrial;
            }
         }

         calculateTotalEnergy(lambda);
      }

      /**
       * Returns the dimensionless Hertzian pair energy beta*u(r).
      */
      private double hertzEnergy(double r) {
         if (r >= sigma) {
            return 0.0;
         }

         double overlap = 1.0 - r / sigma;
         return betaEpsilon * Math.pow(overlap, 2.5);
      }

      public void changeDensity(double newDensity) {
         double newSide = sigma * Math.cbrt(N / newDensity);
         double scaleFactor = newSide / side;

         for (int i = 0; i < N; i++) {
            x[i] *= scaleFactor;
            y[i] *= scaleFactor;
            z[i] *= scaleFactor;
         }

         density = newDensity;
         side = newSide;
         totalVol = side * side * side;

         steps = 0;
         numberOfConfigurations = 0.0;
         energyAccumulator = 0.0;
         pairEnergyAccumulator = 0.0;
         freeEnergyAccumulator = 0.0;
         virialAccumulator = 0.0;

         calculateTotalEnergy(lambda);
      }

      /**
       * Returns the dimensionless pair virial contribution -beta*r*du/dr.
      */
      private double hertzVirial(double r) {
         if (r >= sigma) {
            return 0.0;
         }

         double overlap = 1.0 - r / sigma;

         return 2.5 * betaEpsilon
                  * (r / sigma)
                  * Math.pow(overlap, 1.5);
      }

      public void calculateTotalEnergy(double lambda) {
         totalPairEnergy = 0.0;
         totalVirial = 0.0;
         totalFreeEnergy = 0.0;
         totalEnergy = 0.0;
         springEnergySum = 0.0;
         squaredDisplacementSum = 0.0;

         for (int i = 0; i < N; i++) {
            energy[i] = 0.0;

            for (int j = 0; j < N; j++) {
                  pairEnergy[i][j] = 0.0;
            }
         }

         // Calculate Hertzian energy and virial using each pair once.
         for (int i = 0; i < N - 1; i++) {
            for (int j = i + 1; j < N; j++) {
                  double xij = PBC.separation(x[i] - x[j], side);
                  double yij = PBC.separation(y[i] - y[j], side);
                  double zij = PBC.separation(z[i] - z[j], side);

                  double r = Math.sqrt(
                        xij*xij + yij*yij + zij*zij
                  );

                  double pairU = hertzEnergy(r);
                  double pairW = hertzVirial(r);

                  pairEnergy[i][j] = pairU;
                  pairEnergy[j][i] = pairU;

                  energy[i] += pairU;
                  energy[j] += pairU;

                  totalPairEnergy += pairU;
                  totalVirial += pairW;
            }
         }

         // Einstein energy is needed only for solid thermodynamic integration.
         for (int i = 0; i < N; i++) {
            double dxi = x[i] - x0[i] - dxOverN;
            double dyi = y[i] - y0[i] - dyOverN;
            double dzi = z[i] - z0[i] - dzOverN;

            double displacementSquared =
                     dxi*dxi + dyi*dyi + dzi*dzi;

            squaredDisplacementSum += displacementSquared;
            springEnergySum += springConstant * displacementSquared;
         }

         // The physical energy is Hertzian only.
         totalEnergy = totalPairEnergy;

         if (steps > delay && (steps - delay) % snapshotInterval == 0) {
            numberOfConfigurations++;

            pairEnergyAccumulator += totalPairEnergy;
            energyAccumulator += totalEnergy;
            virialAccumulator += totalVirial;
            squaredDisplacementAccumulator += squaredDisplacementSum;
            springEnergyAccumulator += springEnergySum;
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

      /**
       * Compute hard-sphere second virial coefficient for comparison.
       * For hard spheres: B2 = (2π/3) * σ^3
       */
      public double hardSphereB2(double sigma) {
         return (2.0 * Math.PI / 3.0) * Math.pow(sigma, 3.0);
      }

      /**
       * Computes the second virial coefficient of the fixed-size
      * Hertzian pair potential.
      */
      public double secondVirialCoefficient() {
         /*
          * B2 = 2*pi*integral_0^sigma r^2[1-exp(-beta*u(r))]dr.
          * Composite Simpson integration is used because the former
          * 10-point Gauss-Legendre rule loses about 1.8% at
          * betaEpsilon = 12000.
          */
         final int intervals = 2000; // even, deterministic, and well converged
         final double h = sigma / intervals;
         double sum = 0.0;

         for (int i = 0; i <= intervals; i++) {
            double r = i * h;
            double integrand =
               r * r * (1.0 - Math.exp(-hertzEnergy(r)));

            if (i == 0 || i == intervals) {
               sum += integrand;
            } else if (i % 2 == 0) {
               sum += 2.0 * integrand;
            } else {
               sum += 4.0 * integrand;
            }
         }

         return 2.0 * Math.PI * h * sum / 3.0;
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

      public double meanModifiedRadius() {
         if (numberOfConfigurations == 0) return 0;
         return modifiedSwellingRatioAccumulator / N / numberOfConfigurations;
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
