package org.opensourcephysics.sip.Hertz;

import org.opensourcephysics.numerics.*;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

/**
 * HertzSpheres performs a Monte Carlo simulation of nonionic microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @authors Alan Denton and Matt Urich
 * @version 1.2 25-05-2021
 * 
 */

public class HertzSpheresSolidPhaseWithFaceting {
   
    public int N; // number of particles
    public String initConfig; // initial configuration of particles
    public int nx; // number of columns and rows in initial crystal lattice
    public double side, totalVol; // side length and volume of cubic simulation box 
    // Note: lengths are in units of dry microgel radius, energies in thermal (kT) units
    public double d; // distance between neighboring lattice sites 
    public double x[], y[], z[], a[]; // coordinates (x, y, z) and radius (a) of particles 
    public double x0[], y0[], z0[];//to store the initial positions
    public double oldFloryFRParticleJ[], oldFloryFRParticleI[];
    public double energy[], energy1[]; // total energies of particles (Hertzian plus Flory-Rehner)
    public double newNMon[];
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
    public double boltzmannFactorAccumulator, boltzmannFactor, squaredDisplacementAccumulator;
    public double dxOverN, dyOverN, dzOverN;
    public String fileExtension; 
    public double numberOfConfigurations;
    public double density, nnDistance;
    public double squaredDisplacement, squaredDisplacementSum;
    public double newFloryFR, mixF, elasticF, FloryFR;
    public double volCap[], volumeOfMicrogelI[], newFlory[], oldFloryFR[];
    // Create a new Random object
    public Random random = new Random();

    public boolean initialized = false; //boolean to use to track initializing the first element

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

      oldFloryFRParticleJ = new double[N]; // array to store the Flory-Rehner free energy for each particle before the trial move
      oldFloryFRParticleI = new double[N];
      volumeOfMicrogelI = new double[N];
      newNMon = new double[N];
      newFlory = new double[N];
      oldFloryFR = new double[N];

      a = new double[N]; // particle radii (units of dry radius)

      monRadius = 0.3; // estimate of monomer radius [nm]
      nMon = 0.63*Math.pow(dryR, 3.)/Math.pow(monRadius, 3.); // number of monomers in a microgel
      nChains = xLinkFrac*nMon; // number of chains in a microgel 
      reservoirSR = reservoirSwellingRatio(nMon, nChains, chi);

      reservoirVolFrac = dryVolFrac*Math.pow(reservoirSR, 3.);

      for (int i=0; i<N; i++){
	      a[i] = reservoirSR; // initial particle radii (fully swollen)
      }

      energy = new double[N]; // total energies of particles
      volCap = new double[N]; // the volume of the microgel cap after faceting
      
      // Set only the first element to 2
      if(!initialized){
         volCap[0] = 2;
         initialized = true;
      }

      pairEnergy = new double[N][N]; // pair energies of particles
      newPairEnergy = new double[N][N]; // for lambda = 1
      steps = 0;
      energyAccumulator = 0;
      pairEnergyAccumulator = 0; //for the lamba = 1 system
      freeEnergyAccumulator = 0;
      virialAccumulator = 0;

      for(int i=0; i<numberBins; i++){ // initialize size histogram
         sizeDist[i] = 0;
      }

      // initialize positions
      if (configuration.toUpperCase().equals("SC")) {
         scale = 1; // no scaling of Hertz interactions
         setSCpositions();
      }
      if (configuration.toUpperCase().equals("FCC")) {
         scale = 1; // no scaling of HertzmicrogelVol
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
      double xij, yij, zij, de, r, r2, sigma, da2;
      double mixFParticleI, elasticFParticleI, totalF, mixFParticleJ, elasticFParticleJ;
      double amin, overlapVol, facetVolPerGel;
      double totalFR = 0;
      double facetVolAccumulator = 0, overlapVolSum;
      boolean isFirstIteration = true; // Initialize a boolean flag

      for (int i = 0; i < N; i++) { // attempt a trial move (displacement and size change) 
         //volume of microgel before the trial move
         volumeOfMicrogelI[i] = (4/3)*Math.PI*Math.pow(a[i], 3);

         dxtrial = tolerance*2.*(random.nextDouble()-0.5);         
         dytrial = tolerance*2.*(random.nextDouble()-0.5);
         dztrial = tolerance*2.*(random.nextDouble()-0.5);
         datrial = atolerance*2.*(random.nextDouble()-0.5);

         // Euclidean distance of trial displacements
         trialDisplacementDistance = (dxtrial*dxtrial)+(dytrial*dytrial)+(dztrial*dztrial);
         
         //the Euclidean distance between the current position and the initial position of each particle and the accumulated shift in the center of mass
         double dxi = x[i]-x0[i]-dxOverN;
         double dyi = y[i]-y0[i]-dyOverN;
         double dzi = z[i]-z0[i]-dzOverN;
         // displacement = (dxi*dxi)+(dyi*dyi)+(dzi*dzi);
         
         if (isFirstIteration){
            //the new number of monomers
            newNMon[i] = nMon*(volumeOfMicrogelI[i]-volCap[0])/volumeOfMicrogelI[i]; //number of Monomers to compute the new FR free energy

            //Flory-Rehner Free energy before size and displacement changes for particle i the first time thorugh
            mixF = newNMon[i]*((a[i]*a[i]*a[i]-1)*Math.log(1-1/a[i]/a[i]/a[i])+chi*(1-1/a[i]/a[i]/a[i]));
            elasticF = 1.5*newNMon[i]*xLinkFrac*(a[i]*a[i]-Math.log(a[i])-1);
            oldFloryFR[i] = mixF + elasticF; // total FR before the trial move for particle i
            isFirstIteration = false; // Set the flag to false after the first iteration

         }

         x[i] += dxtrial; // trial displacement for particle i
         y[i] += dytrial;                       
         z[i] += dztrial;                 
         a[i] += datrial; // trial radius change for particle i

         pairPotentialCorrection = 2*springConstant*Math.sqrt((dxtrial*dxi+dytrial*dyi+dztrial*dzi)*(dxtrial*dxi+dytrial*dyi+dztrial*dzi))+(1-1/(double)N)*trialDisplacementDistance; // change in spring energy
         pairEnergySum = 0;
         overlapVolSum=0;

         for (int j = 0; j < N; j++) {
            
            if(j != i) { // consider interactions with other particles

               xij = PBC.separation(x[i]-x[j], side);
               yij = PBC.separation(y[i]-y[j], side);
               zij = PBC.separation(z[i]-z[j], side);

               r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
	            sigma = a[i] + a[j];
               da2 = Math.pow(a[i]-a[j], 2); // difference of radii squared

	            if (r2 < sigma*sigma){
                  r = Math.sqrt(r2);
                  B = scale*Young*nChains*Math.pow(sigma, 2.)*Math.sqrt(a[i]*a[j])/(Math.pow(a[i], 3.)+Math.pow(a[j], 3.)); 
		            // Hertz pair potential amplitude (scaled by factor of scale*Young)
                  HertzEnergy = B*Math.pow(1-r/sigma, 2.5);
                  newPairEnergy[i][j] = HertzEnergy;
                  pairEnergySum+=HertzEnergy;

                  overlapVol = 0.5*(Math.PI*(sigma-r)*(sigma-r)*(r2+2*r*sigma-3*da2)/12./r); // volume of lens-shaped overlap region
                  overlapVolSum+=overlapVol;
               // facetVolPerGel = 0.5*overlapVol;

               }
            }
         } 
         volCap[i] = overlapVolSum;
         newNMon[i] = nMon*(volumeOfMicrogelI[i]-volCap[i])/volumeOfMicrogelI[i]; //number of Monomers to compute the new FR free energy
         
         // Flory-Rehner single-particle free energy (associated with swelling) after trial changes
         mixF = newNMon[i]*((a[i]*a[i]*a[i]-1)*Math.log(1-1/a[i]/a[i]/a[i])+chi*(1-1/a[i]/a[i]/a[i]));
         elasticF = 1.5*newNMon[i]*xLinkFrac*(a[i]*a[i]-Math.log(a[i])-1);
         newFloryFR = elasticF + mixF; // newFloryFR after faceting

         // the change in energy
         de = lambda*(pairEnergySum-energy[i])+(1-lambda)*(pairPotentialCorrection)+(newFloryFR-oldFloryFR[i]); // change in total energy due to trial move
        
         if(Math.exp(-de) < random.nextDouble()){ // Metropolis algorithm
            System.out.println("Move rejected");
            x[i] -= dxtrial; // reject move
            y[i] -= dytrial;
            z[i] -= dztrial;
	         a[i] -= datrial;
         } 
         else { // accept move and update energies
            System.out.println("Move accepted");
            // change in displacement per particle trial move: delta r/N
            dxOverN += dxtrial/(double)N;
            dyOverN += dytrial/(double)N;
            dzOverN += dztrial/(double)N;

            energy[i] = pairEnergySum; // update energy of moved particle

            // Update Flory-Rehner Free Energy based on acceptance
            oldFloryFR[i] = newFloryFR;

            for (int j = 0; j < N; ++j) { // update energies of other particles 
               if(j != i) { // ensures that the particle i is not interacting with itself
                  energy[j] += newPairEnergy[i][j]-pairEnergy[i][j]; // For each pair of particles (i, j), it updates the energy of particle j based on the change in pair energy
                  pairEnergy[i][j] = newPairEnergy[i][j]; // updates the pair energy for the pair (i, j) with the new value calculated in the trial move 
                  pairEnergy[j][i] = newPairEnergy[i][j];
               }
            }
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
      double xij, yij, zij, r, r2, sigma, dxi, dyi, dzi, da2;
      double derivPart, HertzEnergy;
      double mixF, elasticF, totalF;
      double fOverR, fx, fy, fz;
      double amin, overlapVol, facetVolPerGel;
      double facetVolAccumulator = 0, overlapVolSum;

      for(int i = 0; i < N; ++i) {
         volumeOfMicrogelI[i] = (4/3)*Math.PI*Math.pow(a[i], 3);

         pairEnergySum = 0; //for the lamba = 1 system
         overlapVolSum=0;
         virialSum = 0;
         
         // the Euclidean distance between the current position and the initial position of each particle
         dxi = x[i]-x0[i]-dxOverN;
         dyi = y[i]-y0[i]-dyOverN;
         dzi = z[i]-z0[i]-dzOverN;
        
         squaredDisplacement = Math.pow(dxi, 2)+Math.pow(dyi, 2)+Math.pow(dzi, 2); // the square displacement 
         squaredDisplacementSum+=squaredDisplacement;
         springEnergy = springConstant*squaredDisplacement; //since it is a single particle (i) component
         springEnergySum += springEnergy; 
         for(int j = 0; j < N; ++j) {
            if(j != i) { // consider interactions with other particles
               xij = PBC.separation(x[i]-x[j], side);
               yij = PBC.separation(y[i]-y[j], side);
               zij = PBC.separation(z[i]-z[j], side);
               r2 = xij*xij + yij*yij + zij*zij; // particle separation squared

	            sigma = a[i]+ a[j];
               da2 = Math.pow(a[i]-a[j], 2); // difference of radii squared
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

                  overlapVol = 0.5*(Math.PI*(sigma-r)*(sigma-r)*(r2+2*r*sigma-3*da2)/12./r); // volume of lens-shaped overlap region
                  overlapVolSum+=overlapVol;
               // facetVolPerGel = 0.5*overlapVol;

	            }
            }
         }
         volCap[i] = overlapVolSum;
         newNMon[i] = nMon*(volumeOfMicrogelI[i]-volCap[i])/volumeOfMicrogelI[i]; //number of Monomers to compute the new FR free energy

         // Flory-Rehner single-particle free energy (associated with swelling)
	      mixF = newNMon[i]*((a[i]*a[i]*a[i]-1)*Math.log(1-1/a[i]/a[i]/a[i])+chi*(1-1/a[i]/a[i]/a[i]));
	      elasticF = 1.5*newNMon[i]*xLinkFrac*(a[i]*a[i]-Math.log(a[i])-1);
         //elasticF = 1.5*nMon*xLinkFrac*(a[i]*a[i]-Math.log(a[i])-1);
         
	      totalF = elasticF + mixF; // Flory-Rehner free energy after the trial move
         
         energy[i] = pairEnergySum; // total energy of particle i for the lambda system // because the change in energy is computed in the step method and the change has been used in the step method

         mixFRSR =  newNMon[i]*((reservoirSR*reservoirSR*reservoirSR-1)*Math.log(1-1/reservoirSR/reservoirSR/reservoirSR)+chi*(1-1/reservoirSR/reservoirSR/reservoirSR));
         elasticFRSR = 1.5*nChains*(reservoirSR*reservoirSR-Math.log(reservoirSR)-1);

         totalFRSR = mixFRSR + elasticFRSR; // total free energy associated with the reservoir state
         //totalF = totalF - totalFRSR;
         newFloryFR = totalF - totalFRSR;
         totalF = newFloryFR;

         oldFloryFR[i] = newFloryFR; 

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
               boltzmannFactorAccumulator += boltzmannFactor;
               springEnergyAccumulator += springEnergySum; // accumulates all the spring energies
               squaredDisplacementAccumulator += squaredDisplacementSum; // accumulates the mean sqaure displacement
               
            }
         }

   }

   // mean energy per particle [kT units]
   public double meanEnergy() {
      return energyAccumulator/N/numberOfConfigurations; // quantity <E>/N
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
    public void calculateVolumeFraction() { // instantaneous volume fraction 
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
        volFrac += microgelVol/totalVol;
    }

    /* reservoir swelling ratio as root of Flory-Rehner pressure */
    public double reservoirSwellingRatio(double nMon, double nChains, double chi){
      Function f = new FloryRehnerPressure(nMon, nChains, chi);
      double xleft = 1.1;
      double xright = 6.;
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
