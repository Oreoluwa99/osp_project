package org.opensourcephysics.sip.Hertz;

import org.opensourcephysics.numerics.*;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * HertzSpheres performs a Monte Carlo simulation of nonionic microgels interacting via the 
 * Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @authors Alan Denton and Matt Urich
 * @version 1.2 25-05-2021
 * 
 */

public class HertzSpheresLiquidLambdaEdited {
	
    public int N; // number of particles
    public String initConfig; // initial configuration of particles
    public int nx; // number of columns and rows in initial crystal lattice
    public double side, totalVol; // side length and volume of cubic simulation box 
    // Note: lengths are in units of dry microgel radius, energies in thermal (kT) units
    public double d; // distance between neighboring lattice sites 
    public double x[], y[], z[], a[]; // coordinates (x, y, z) and radius (a) of particles 
    public double energy[]; // total energies of particles (Hertzian plus Flory-Rehner)
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
    public double dryVolFrac, volFrac, reservoirVolFrac; 
    public int steps; // number of Monte Carlo (MC) steps
    public double delay, stop; // MC steps after which statistics are collected and not collected
    public int snapshotInterval; // interval by which successive samples are separated
    public double sizeDist[], sizeBinWidth; // particle radius histogram and bin width
    public int numberBins; // number of histogram bins 
    public double grBinWidth, maxRadius, deltaK; // bin widths and range
    public double meanR; // mean particle radius
    public String fileExtension; 
    public double lambda, dlambda, dphi;
    public double mixFRSR, elasticFRSR, totalFRSR;
    public double numberOfConfigurations; // number of configurations
    public double meanVirialFreeEnergy;


   /**
    * Initialize the model.
    * 
    * @param configuration
    *  Initial lattice structure 
    */
   public void initialize(String configuration) {
      x = new double[N]; // particle coordinates
      y = new double[N];
      z = new double[N];
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
      pairEnergy = new double[N][N]; // pair energies of particles
      newPairEnergy = new double[N][N];
      steps = 1;
      energyAccumulator = 0;
      pairEnergyAccumulator = 0;
      freeEnergyAccumulator = 0;
      virialAccumulator = 0;
      numberOfConfigurations = 0;
      volFrac = 0; // counter for system volume fraction
      side = Math.cbrt(4.*Math.PI*N/dryVolFrac/3.); // side length of cubic simulation box
      totalVol = side*side*side; // box volume [units of dry radius cubed]
      //System.out.println("total volume: " + totalVol);

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
      double newEnergy, HertzEnergy;
      double xij, yij, zij, de, r, r2, sigma;
      double mixF, elasticF, totalF;
      
      for (int i = 0; i < N; i++) { // attempt a trial move (displacement and size change)
	  
         dxtrial = tolerance*2.*(Math.random()-0.5);
         dytrial = tolerance*2.*(Math.random()-0.5);
         dztrial = tolerance*2.*(Math.random()-0.5);
         datrial = atolerance*2.*(Math.random()-0.5);
	 
         x[i] = PBC.position(x[i]+dxtrial, side); // trial displacement (periodic boundary conditions)
         y[i] = PBC.position(y[i]+dytrial, side);
         z[i] = PBC.position(z[i]+dztrial, side);
	      a[i] += datrial; // trial radius change
	 
         newEnergy = 0; // keep track of the change in energy after trial move
         
         for (int j = 0; j < N; j++) {
            if(j != i) { // consider interactions with other particles
               xij = PBC.separation(x[i]-x[j], side);
               yij = PBC.separation(y[i]-y[j], side);
               zij = PBC.separation(z[i]-z[j], side);
               r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
	            sigma = a[i] + a[j];
	            if (r2 < sigma*sigma){
                  r = Math.sqrt(r2);
                  // Hertz pair potential amplitude (scaled by factor of scale*Young)
                  B = lambda*scale*Young*nChains*Math.pow(sigma, 2.)*Math.sqrt(a[i]*a[j])/(Math.pow(a[i], 3.)+Math.pow(a[j], 3.)); 
		            HertzEnergy = B*Math.pow(1-r/sigma, 2.5);
                  newPairEnergy[i][j] = HertzEnergy;
		            newEnergy += HertzEnergy;
	            }
            }
            
         } 
	 
         // Flory-Rehner single-particle free energy (associated with swelling)
	      mixF = nMon*((a[i]*a[i]*a[i]-1)*Math.log(1-1/a[i]/a[i]/a[i])+chi*(1-1/a[i]/a[i]/a[i]));
	      elasticF = 1.5*nChains*(a[i]*a[i]-Math.log(a[i])-1);
	      totalF = elasticF + mixF;

         // Flory-Rehner single particle free energy for the fully swollen microgel
         mixFRSR = nMon*((reservoirSR*reservoirSR*reservoirSR-1)*Math.log(1-1/reservoirSR/reservoirSR/reservoirSR)+chi*(1-1/reservoirSR/reservoirSR/reservoirSR));
         elasticFRSR = 1.5*nChains*(reservoirSR*reservoirSR-Math.log(reservoirSR)-1);

         totalFRSR = mixFRSR + elasticFRSR;
         totalF = totalF - totalFRSR;
	 
	      newEnergy += totalF; // total energy after trial move
         de = newEnergy-energy[i]; // change in total energy due to trial move
	 
         if(Math.exp(-de) < Math.random()) { // Metropolis algorithm
            x[i] = PBC.position(x[i]-dxtrial, side); // reject move
            y[i] = PBC.position(y[i]-dytrial, side);
            z[i] = PBC.position(z[i]-dztrial, side);
	         a[i] -= datrial;
         }
         else { // accept move and update energies
            energy[i] = newEnergy; // update energy of moved particle
            for (int j = 0; j < N; ++j) { // update energies of other particles 
               if(j != i) {
                  energy[j] += newPairEnergy[i][j]-pairEnergy[i][j];
                  pairEnergy[i][j] = newPairEnergy[i][j];
                  pairEnergy[j][i] = newPairEnergy[i][j];
               }
            }
         }
      }
      calculateTotalEnergy(lambda); // new total energy
	  
   }

   // total energy and virial (for pressure)
   public void calculateTotalEnergy(double lambda) {
      totalEnergy = 0;
      totalVirial = 0;
      totalPairEnergy = 0;
      totalFreeEnergy = 0;
      double pairEnergySum, virialSum;
      double xij, yij, zij, r, r2, sigma;
      double derivPart, HertzEnergy;
      double mixF, elasticF, totalF;
      double fOverR, fx, fy, fz;
      for(int i = 0; i < N; ++i) {
         pairEnergySum = 0;
         virialSum = 0;
         for(int j = 0; j < N; ++j) {
            if(j != i) { // consider interactions with other particles
               xij = PBC.separation(x[i]-x[j], side);
               yij = PBC.separation(y[i]-y[j], side);
               zij = PBC.separation(z[i]-z[j], side);
               r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
	            sigma = a[i]+ a[j];
	            if (r2 < sigma*sigma){
                  r = Math.sqrt(r2);
                  // Hertz pair potential amplitude (scaled by a factor of scale*Young)
                  B = lambda*scale*Young*nChains*Math.pow(sigma, 2.)*Math.sqrt(a[i]*a[j])/(Math.pow(a[i], 3.)+Math.pow(a[j], 3.)); 
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
            }
         }

         // Flory-Rehner single-particle free energy (associated with swelling)
	      mixF = nMon*((a[i]*a[i]*a[i]-1)*Math.log(1-1/a[i]/a[i]/a[i])+chi*(1-1/a[i]/a[i]/a[i]));
	      elasticF = 1.5*nChains*(a[i]*a[i]-Math.log(a[i])-1);
	      if (a[i] < 0){ // particle radius should never be negative
	         System.out.println("Negative radius: " + a[i]);
	         System.out.println("Elastic free energy: " + elasticF);
	      }
	      totalF = elasticF + mixF; // Flory-Rehner free energy

         mixFRSR = nMon*((reservoirSR*reservoirSR*reservoirSR-1)*Math.log(1-1/reservoirSR/reservoirSR/reservoirSR)+chi*(1-1/reservoirSR/reservoirSR/reservoirSR));
         elasticFRSR = 1.5*nChains*(reservoirSR*reservoirSR-Math.log(reservoirSR)-1);

         totalFRSR = mixFRSR + elasticFRSR;
         totalF = totalF - totalFRSR;

         energy[i] = pairEnergySum + totalF; // total energy of particle i
	      totalPairEnergy += pairEnergySum;
         totalFreeEnergy += totalF;
	      totalVirial += virialSum;

      }

      if (steps > delay){
         if ((steps-delay)%snapshotInterval==0){ // include configurations only after even number of intervals
            numberOfConfigurations++;
            totalVirial *= 0.5; // correct for double counting pairs
            totalPairEnergy *= 0.5; // correct for double counting pairs
            totalEnergy = totalPairEnergy + totalFreeEnergy;
            pairEnergyAccumulator += totalPairEnergy;
            energyAccumulator += totalEnergy; // running totals
            freeEnergyAccumulator += totalFreeEnergy;
            virialAccumulator += totalVirial; 
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
       }
       else
	      meanR = 0;
       
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
