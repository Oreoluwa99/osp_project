package org.opensourcephysics.sip.HertzMixture;

import org.opensourcephysics.numerics.*;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * HertzMixture performs a Monte Carlo simulation of a binary mixture of microgels interacting 
 * via the Hertz elastic pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @author Alan Denton 
 * @version 1.2 07.09.2021
 * 
 */

public class HertzMixture {
	
    public int N1, N2; // numbers of particles of the two species
    public int nx; // number of columns and rows in initial crystal lattice
    public double side, totalVol; // side length and volume of cubic simulation box 
    // NOTE: lengths are in units of dry radius of microgel species 1, energies in thermal (kT) units 
    public double d; // center-to-center distance between lattice points
    public int lattice[][][]; // lattice occupancy array
    public double x1[], y1[], z1[], a1[], x2[], y2[], z2[], a2[]; // coordinates and radii of particles 
    public double energy1[], energy2[]; // total energies of particles (Hertzian plus Flory-Rehner) 
    public double pairEnergy11[][], pairEnergy12[][], pairEnergy22[][]; // pair energies of particles
    public double newPairEnergy11[][], newPairEnergy12[][], newPairEnergy22[][];
    public double totalEnergy, totalPairEnergy, totalFreeEnergy; // total energy of system
    public double totalVirial; // total virial of system (for pressure calculation)
    // accumulators for computing various thermodynamic properties:
    public double energyAccumulator, pairEnergyAccumulator, freeEnergyAccumulator, virialAccumulator;
    public double tolerance, atolerance; // tolerances for trial displacements and radius changes
    public double monRadius, nMon1, nMon2, nChains1, nChains2, chi1, chi2; // microgel parameters
    public double xLinkFrac1, xLinkDensityRatio; // microgel parameters
    public double Young, B11, B12, B22; // Young's modulus calibration, prefactors for Hertz potentials
    public double scale; // scale factor (0 or 1) for prefactor of Hertz potentials
    public double dryRadius1, dryRadius2, dryRadiusRatio, reservoirSR1, reservoirSR2;
    public double R1, R2, meanR1, meanR2, meanSR1, meanSR2; // particle radii and swelling ratios
    public double dryVolFrac1, dryVolFrac2, reservoirVolFrac, volFrac;
    public int steps; // number of Monte Carlo steps
    public double delay, stop; // MC steps after which statistics are collected and stopping point
    public int snapshotInterval; // interval by which successive samples are separated
    public double sizeDist1[], sizeDist2[], sizeBinWidth; // particle radii histograms and bin width
    public int numberBins; // number of histogram bins 
    public double grBinWidth, deltaK, maxSwellingRatio; // bin widths and range
    public String fileExtension; // for file naming

   /**
   * Initialize the model.
   * 
   * @param configuration
   *  Initial lattice structure 
   */
   public void initialize(String configuration) {
      x1 = new double[N1]; // particle coordinates
      y1 = new double[N1];
      z1 = new double[N1];
      a1 = new double[N1]; // particle radii (units of dry radius of species 1)
      x2 = new double[N2];
      y2 = new double[N2];
      z2 = new double[N2];
      a2 = new double[N2];

      monRadius = 0.3; // monomer radius [nm]
      nMon1 = 0.63*Math.pow(dryRadius1/monRadius, 3.); // number of monomers in microgel 1
      nChains1 = xLinkFrac1*nMon1; // number of chains in microgel 1
      dryRadiusRatio = dryRadius2/dryRadius1;
      nMon2 = nMon1*Math.pow(dryRadiusRatio, 3.); // number of monomers in microgel 2
      nChains2 = nChains1*xLinkDensityRatio*Math.pow(dryRadiusRatio, 3.);
      dryVolFrac2 = dryVolFrac1*Math.pow(dryRadiusRatio, 3.)*(double)N2/(double)N1;
      reservoirSR1 = reservoirSwellingRatio(nMon1, nChains1, chi1);
      reservoirSR2 = reservoirSwellingRatio(nMon2, nChains2, chi2);
      reservoirVolFrac = dryVolFrac1*Math.pow(reservoirSR1, 3.) + dryVolFrac2*Math.pow(reservoirSR2, 3.);

      for (int i = 0;i<N1;i++){
	  a1[i] = reservoirSR1; // swelling ratio of microgel species 1 (fully swollen)
      }
      for (int i = 0;i<N2;i++){
	  a2[i] = reservoirSR2; // swelling ratio of microgel species 2 (fully swollen)
      }
      energy1 = new double[N1]; // total energies of particles
      energy2 = new double[N2];
      pairEnergy11 = new double[N1][N1]; // pair energies of particles
      pairEnergy12 = new double[N1][N2];
      pairEnergy22 = new double[N2][N2];
      newPairEnergy11 = new double[N1][N1];
      newPairEnergy12 = new double[N1][N2];
      newPairEnergy22 = new double[N2][N2];
      steps = 1;
      energyAccumulator = 0;
      pairEnergyAccumulator = 0;
      freeEnergyAccumulator = 0;
      virialAccumulator = 0;
      volFrac = 0; // counter for system volume fraction
      side = Math.cbrt(4.*Math.PI*N1/dryVolFrac1/3.); // side length of cubic simulation box
      totalVol = side*side*side; // box volume

      numberBins = (int) (maxSwellingRatio/sizeBinWidth);
      sizeDist1 = new double[(int) (numberBins)];
      sizeDist2 = new double[(int) (numberBins)];

      for(int i=0; i<numberBins; i++){ // initialize size histograms
         sizeDist1[i] = 0;
         sizeDist2[i] = 0;
      }

      // initialize positions
      if (configuration.toUpperCase().equals("SC")) {
         scale = 1; // no scaling of Hertz interactions
         setSCPositions();
      } 
      if (configuration.toUpperCase().equals("FCC")) {
         scale = 1; // no scaling of Hertz interactions
         setFCCPositions();
      } 
      if (configuration.toUpperCase().equals("DISORDERED-FCC")) {
         scale = 1; // no scaling of Hertz interactions
         setDisorderedFCCPositions();
      } 
      if (configuration.toUpperCase().equals("RANDOM-FCC")) {
         scale = 0; // B = 0 initially to allow particle positions to randomize
         setRandomFCCPositions();
      } 
      if (configuration.toUpperCase().equals("CENTER")) {
         scale = 1; // no scaling of Hertz interactions
         setCenterPositions();
      } 
   }

   /**
   * Place particles on sites of a simple cubic lattice.
   */
   public void setSCPositions() {
      System.out.println("SC");
      int ix, iy, iz;
      double dnx = Math.cbrt(N1+N2);
      d = side / dnx; // center-center distance between neighboring lattice sites
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // N1 is not a perfect cube
      }

      int i = 0;
      for (iy = 0; iy < nx; iy++) { // loop through particles in a column
         for (ix = 0; ix < nx; ix++) { // loop through particles in a row
            for (iz = 0; iz < nx; iz++) {
               if (i < N1) { // check for remaining particles
                  x1[i] = PBC.position(ix * d, side);
                  y1[i] = PBC.position(iy * d, side);
                  z1[i] = PBC.position(iz * d, side);
                  if (i < N2) {
                     x2[i] = PBC.position(x1[i] + d/2, side);
                     y2[i] = PBC.position(y1[i] + d/2, side);
                     z2[i] = PBC.position(z1[i] + d/2, side);
                  } 
                  i++;
               }
            }
         }
      }
      calculateTotalEnergy(); // initial energy
   }

   /**
   * Place particles on sites of an FCC lattice.
   */
   public void setFCCPositions() {
      System.out.println("FCC");
      int ix, iy, iz;
      double dnx = Math.cbrt((N1+N2)/4.);
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // (N1+N2)/4 is not a perfect cube
      }

      int i = 0;
      for (ix = 0; ix < 2*nx; ix++) { // loop through particles in a row
         for (iy = 0; iy < 2*nx; iy++) { // loop through particles in a column
            for (iz = 0; iz < 2*nx; iz++) { // loop through particles in a layer
               if (i < N1) { // check for remaining particles
                  if ((ix+iy+iz)%2 == 0) { // check for remaining particles
                     x1[i] = PBC.position(ix * d/2., side);
                     y1[i] = PBC.position(iy * d/2., side);
                     z1[i] = PBC.position(iz * d/2., side);
                     if (i < N2) {
                        x2[i] = PBC.position(x1[i] + d/2., side);
                        y2[i] = PBC.position(y1[i] + d/2., side);
                        z2[i] = PBC.position(z1[i] + d/2., side);
                     } 
                     i++;
                  }
               }
            }
         }
      }
      calculateTotalEnergy(); // initial energy
   }

   /**
   * Place particles at positions randomly displaced from sites of an FCC lattice.
   */
   public void setRandomFCCPositions() {
      System.out.println("random FCC");
      int ix, iy, iz;
      double dnx = Math.cbrt((N1+N2)/4.);
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // (N1+N2)/4 is not a perfect cube
      }

      int i = 0;
      for (ix = 0; ix < 2*nx; ix++) { // loop through particles in a row
         for (iy = 0; iy < 2*nx; iy++) { // loop through particles in a column
            for (iz = 0; iz < 2*nx; iz++) { // loop through particles in a layer
               if (i < N1) { // check for remaining particles
                  if ((ix+iy+iz)%2 == 0) { // check for remaining particles
                     x1[i] = PBC.position(ix * d/2. + (Math.random()-0.5) * d, side);
                     y1[i] = PBC.position(iy * d/2. + (Math.random()-0.5) * d, side);
                     z1[i] = PBC.position(iz * d/2. + (Math.random()-0.5) * d, side);
                     if (i < N2) {
                        x2[i] = PBC.position(x1[i] + d/2. + (Math.random()-0.5) * d, side);
                        y2[i] = PBC.position(y1[i] + d/2. + (Math.random()-0.5) * d, side);
                        z2[i] = PBC.position(z1[i] + d/2. + (Math.random()-0.5) * d, side);
                     } 
                     i++;
                  }
               }
            }
         }
      }
      calculateTotalEnergy(); // initial energy
   }

   /**
   * Place particles on sites of a disordered FCC lattice.
   */
   public void setDisorderedFCCPositions() {
      System.out.println("disordered FCC");
      int i, ix, iy, iz;
      boolean iplace;
      double dnx = Math.cbrt((N1+N2)/4.); // number of unit cells
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // (N1+N2)/4 is not a perfect cube
      }
      lattice = new int[2*nx][2*nx][2*nx];

      for (ix = 0; ix < 2*nx; ix++) { // initialize lattice occupancy array
         for (iy = 0; iy < 2*nx; iy++) { 
            for (iz = 0; iz < 2*nx; iz++) { 
               lattice[ix][iy][iz] = 0;
            }
         }
      }

      for (i = 0; i < N2; i++) { // randomly place particles of (minority, big) species 2 
         iplace = false;
         while (iplace == false) {
            ix = (int) (2*nx*Math.random());
            iy = (int) (2*nx*Math.random());
            iz = (int) (2*nx*Math.random());
            if ((ix+iy+iz)%2 == 0) { 
               if (lattice[ix][iy][iz] == 0) {
                  x2[i] = PBC.position(ix * d/2., side);
                  y2[i] = PBC.position(iy * d/2., side);
                  z2[i] = PBC.position(iz * d/2., side);
                  lattice[ix][iy][iz] = 1;
                  iplace = true;
               }
            }
         }
      }

      i = 0;
      for (ix = 0; ix < 2*nx; ix++) { // place particles of species 1 on empty sites
         for (iy = 0; iy < 2*nx; iy++) { 
            for (iz = 0; iz < 2*nx; iz++) { 
               if ((ix+iy+iz)%2 == 0) { // check for remaining particles
                  if (lattice[ix][iy][iz] == 0) {
                     x1[i] = PBC.position(ix * d/2., side);
                     y1[i] = PBC.position(iy * d/2., side);
                     z1[i] = PBC.position(iz * d/2., side);
                     i++;
                  }
               }
            }
         }
      }
      calculateTotalEnergy(); // initial energy
   }

   /**
   * Place all particles at the center of the cell.
   */
   public void setCenterPositions() {
      int i;
      for(i = 0; i < N1; i++) {
          x1[i] = side/2.;
          y1[i] = side/2.;
          z1[i] = side/2.;
      }
      for(i = 0; i < N2; i++) {
          x2[i] = side/2.;
          y2[i] = side/2.;
          z2[i] = side/2.;
      }
      calculateTotalEnergy(); // initial energy
   }

   /**
   * Do a Monte Carlo simulation step.
   */
   public void step() {
      steps++;
      double dxtrial, dytrial, dztrial, atrial;
      double newEnergy, HertzEnergy;
      double dx, dy, dz, de, r, r2, sigma;
      double mixF, elasticF, totalF;

      for (int i = 0; i < N1; ++i) { // attempt a trial move (displacement and size change) of species 1
	  
         dxtrial = tolerance*2.*(Math.random()-0.5);
         dytrial = tolerance*2.*(Math.random()-0.5);
         dztrial = tolerance*2.*(Math.random()-0.5);
         atrial = atolerance*2.*(Math.random()-0.5);
	 
         x1[i] = PBC.position(x1[i]+dxtrial, side); // trial displacement (periodic boundary conditions)
         y1[i] = PBC.position(y1[i]+dytrial, side);
         z1[i] = PBC.position(z1[i]+dztrial, side);
	 a1[i] += atrial; // trial change in swelling ratio (not radius) 
	 
         newEnergy = 0; // keep track of the change in energy after trial move
	 
         // consider interactions with other particles
         for (int j = 0; j < N1; ++j) { // 11 pairs
            if(j != i) {
               dx = PBC.separation(x1[i]-x1[j], side);
               dy = PBC.separation(y1[i]-y1[j], side);
               dz = PBC.separation(z1[i]-z1[j], side);
               r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
	       sigma = a1[i] + a1[j]; // sum of particle radii
	       if (r2 < sigma*sigma){
                  r = Math.sqrt(r2);
                  // Hertz pair potential amplitude (scaled by a factor of scale*Young)
                  B11 = scale*Young*nChains1*sigma*sigma*Math.sqrt(a1[i]*a1[j])/(Math.pow(a1[i], 3.)+Math.pow(a1[j], 3.));
		  HertzEnergy = B11*Math.pow(1-r/sigma, 2.5);
		  newEnergy += HertzEnergy;
                  newPairEnergy11[i][j] = HertzEnergy;
	       }
            }
         }

         for (int j = 0; j < N2; ++j) { // 12 pairs
            dx = PBC.separation(x1[i]-x2[j], side);
            dy = PBC.separation(y1[i]-y2[j], side);
            dz = PBC.separation(z1[i]-z2[j], side);
            r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
            R1 = a1[i]; // radius of particle 1
            R2 = a2[j]*dryRadiusRatio; // radius of particle 2
	    sigma = R1 + R2; // sum of particle radii
            B12 = scale*Young*sigma*sigma*Math.sqrt(R1*R2)*(nChains1*nChains2/(nChains1*Math.pow(R2, 3.)+nChains2*Math.pow(R1, 3.)));
	    if (r2 < sigma*sigma){
                r = Math.sqrt(r2);
	        HertzEnergy = B12*Math.pow(1-r/sigma, 2.5);
	        newEnergy += HertzEnergy;
                newPairEnergy12[i][j] = HertzEnergy;
	    }
         }
	 
         // Flory-Rehner single-particle free energy (associated with swelling)
	 mixF = nMon1*((a1[i]*a1[i]*a1[i]-1)*Math.log(1-1/a1[i]/a1[i]/a1[i])+chi1*(1-1/a1[i]/a1[i]/a1[i]));
	 elasticF = 1.5*nChains1*(a1[i]*a1[i]-Math.log(a1[i])-1);
	 totalF = elasticF + mixF;
	 
	 newEnergy += totalF;
         de = newEnergy-energy1[i];
	 
         if(Math.exp(-de) < Math.random()) {
            x1[i] = PBC.position(x1[i]-dxtrial, side); // reject move
            y1[i] = PBC.position(y1[i]-dytrial, side);
            z1[i] = PBC.position(z1[i]-dztrial, side);
	    a1[i] -= atrial;
         }
         else { // accept move and update energies
            energy1[i] = newEnergy;
            for (int j = 0; j < N1; ++j) { // update energies of other particles 
               if(j != i) {
                  energy1[j] += newPairEnergy11[i][j]-pairEnergy11[i][j];
                  pairEnergy11[i][j] = newPairEnergy11[i][j];
                  pairEnergy11[j][i] = newPairEnergy11[i][j];
               }
            }
            for (int j = 0; j < N2; ++j) { // update energies of other particles 
               energy2[j] += newPairEnergy12[i][j]-pairEnergy12[i][j];
               pairEnergy12[i][j] = newPairEnergy12[i][j];
            }
         }
      }

      for (int i = 0; i < N2; ++i) { // attempt a trial move (displacement and size change) of species 2
	  
         dxtrial = tolerance*2.*(Math.random()-0.5);
         dytrial = tolerance*2.*(Math.random()-0.5);
         dztrial = tolerance*2.*(Math.random()-0.5);
         atrial = atolerance*2.*(Math.random()-0.5);
	 
         x2[i] = PBC.position(x2[i]+dxtrial, side);
         y2[i] = PBC.position(y2[i]+dytrial, side);
         z2[i] = PBC.position(z2[i]+dztrial, side);
	 a2[i] += atrial; //trial swelling ratio change
	 
         newEnergy = 0;
	 
         for (int j = 0; j < N2; ++j) { // 22 pairs
            if(j != i) {
               dx = PBC.separation(x2[i]-x2[j], side);
               dy = PBC.separation(y2[i]-y2[j], side);
               dz = PBC.separation(z2[i]-z2[j], side);
               r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
	       sigma = (a2[i] + a2[j])*dryRadiusRatio; // sum of particle radii
	       if (r2 < sigma*sigma){
                  r = Math.sqrt(r2);
                  // Hertz pair potential amplitude (scaled by a factor of scale*Young)
                  B22 = scale*Young*nChains2*sigma*sigma*Math.sqrt(a2[i]*a2[j])/(Math.pow(a2[i], 3.)+Math.pow(a2[j], 3.));
		  HertzEnergy = B22*Math.pow(1-r/sigma, 2.5);
		  newEnergy += HertzEnergy;
                  newPairEnergy22[i][j] = HertzEnergy;
	       }
            }
         }

         for (int j = 0; j < N1; ++j) { // 21 pairs
            dx = PBC.separation(x2[i]-x1[j], side);
            dy = PBC.separation(y2[i]-y1[j], side);
            dz = PBC.separation(z2[i]-z1[j], side);
            r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
            R1 = a1[j]; // radius of particle 1
            R2 = a2[i]*dryRadiusRatio; // radius of particle 2
	    sigma = R1 + R2; // sum of particle radii
            B12 = scale*Young*sigma*sigma*Math.sqrt(R1*R2)*(nChains1*nChains2/(nChains1*Math.pow(R2, 3.)+nChains2*Math.pow(R1, 3.)));
	    if (r2 < sigma*sigma){
                r = Math.sqrt(r2);
	        HertzEnergy = B12*Math.pow(1-r/sigma, 2.5);
	        newEnergy += HertzEnergy;
                newPairEnergy12[j][i] = HertzEnergy;
	    }
         }
	 
         // Flory-Rehner single-particle free energy (associated with swelling)
	 mixF = nMon2*((a2[i]*a2[i]*a2[i]-1)*Math.log(1-1/a2[i]/a2[i]/a2[i])+chi2*(1-1/a2[i]/a2[i]/a2[i]));
	 elasticF = 1.5*nChains2*(a2[i]*a2[i]-Math.log(a2[i])-1);
	 totalF = elasticF + mixF;
	 
	 newEnergy += totalF;
         de = newEnergy-energy2[i];
	 
         if(Math.exp(-de) < Math.random()) {
            x2[i] = PBC.position(x2[i]-dxtrial, side); // reject move
            y2[i] = PBC.position(y2[i]-dytrial, side);
            z2[i] = PBC.position(z2[i]-dztrial, side);
	    a2[i] -= atrial;
         }
         else { // accept move and update energies
            energy2[i] = newEnergy;
            for (int j = 0; j < N1; ++j) { // update energies of other particles 
               energy1[j] += newPairEnergy12[j][i]-pairEnergy12[j][i];
               pairEnergy12[j][i] = newPairEnergy12[j][i];
            }
            for (int j = 0; j < N2; ++j) { // update energies of other particles 
               if(j != i) {
                  energy2[j] += newPairEnergy22[i][j]-pairEnergy22[i][j];
                  pairEnergy22[i][j] = newPairEnergy22[i][j];
                  pairEnergy22[j][i] = newPairEnergy22[i][j];
               }
            }
         }
      }
      calculateTotalEnergy(); // new total energy
   }


   // total energy and virial
   public void calculateTotalEnergy() {
      totalEnergy = 0;
      totalPairEnergy = 0;
      totalVirial = 0;
      double pairEnergySum, virialSum;
      double dx, dy, dz, r, r2, sigma;
      double derivPart, HertzEnergy;
      double mixF, elasticF, totalF;
      double fOverR, fx, fy, fz;
      int i, j;

// loop over species 1

      for(i = 0; i < N1; ++i) { 
         pairEnergySum = 0;
         virialSum = 0;
         // consider interactions with other particles
         for(j = 0; j < N1; ++j) { // loop over 11 pairs
            if(j != i) {
               dx = PBC.separation(x1[i]-x1[j], side);
               dy = PBC.separation(y1[i]-y1[j], side);
               dz = PBC.separation(z1[i]-z1[j], side);
               r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
	       sigma = a1[i]+ a1[j]; // sum of particle radii
	       if (r2 < sigma*sigma){
                  r = Math.sqrt(r2);
                  // Hertz pair potential amplitude (scaled by a factor of scale*Young)
                  B11 = scale*Young*nChains1*sigma*sigma*Math.sqrt(a1[i]*a1[j])/(Math.pow(a1[i], 3.)+Math.pow(a1[j], 3.));
		  derivPart = B11*Math.pow(1-r/sigma, 1.5);
		  HertzEnergy = derivPart*(1-r/sigma);
		  pairEnergySum += HertzEnergy;
                  pairEnergy11[i][j] = HertzEnergy; // pair energy of particles i and j
/*
		  fOverR = (2.5/sigma)*derivPart/r;
		  fx = fOverR*dx; // force in x-direction
		  fy = fOverR*dy; // force in y-direction
		  fz = fOverR*dz; // force in z-direction
		  virialSum += dx*fx+dy*fy+dz*fz;
*/
	       }
            }
         }

         for(j = 0; j < N2; ++j) { // loop over 12 pairs
            dx = PBC.separation(x1[i]-x2[j], side);
            dy = PBC.separation(y1[i]-y2[j], side);
            dz = PBC.separation(z1[i]-z2[j], side);
            r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
            R1 = a1[i]; // radius of particle 1
            R2 = a2[j]*dryRadiusRatio; // radius of particle 2
	    sigma = R1 + R2; // sum of particle radii
            B12 = scale*Young*sigma*sigma*Math.sqrt(R1*R2)*(nChains1*nChains2/(nChains1*Math.pow(R2, 3.)+nChains2*Math.pow(R1, 3.)));
	    if (r2 < sigma*sigma){
               r = Math.sqrt(r2);
	       derivPart = B12*Math.pow(1-r/sigma, 1.5);
	       HertzEnergy = derivPart*(1-r/sigma);
	       pairEnergySum += HertzEnergy;
               pairEnergy12[i][j] = HertzEnergy; // pair energy of particles i and j
/*
	       fOverR = (2.5/sigma)*derivPart/r;
	       fx = fOverR*dx; // force in x-direction
	       fy = fOverR*dy; // force in y-direction
	       fz = fOverR*dz; // force in z-direction
	       virialSum += dx*fx+dy*fy+dz*fz;
*/
	    }
         }

         // Flory-Rehner single-particle free energy (associated with swelling)
	 mixF = nMon1*((a1[i]*a1[i]*a1[i]-1)*Math.log(1-1/a1[i]/a1[i]/a1[i])+chi1*(1-1/a1[i]/a1[i]/a1[i]));
	 elasticF = 1.5*nChains1*(a1[i]*a1[i]-Math.log(a1[i])-1);
/*
	 if (a1[i] < 0){
	     System.out.println("Negative Swelling Ratio: " + a1[i]);
	     System.out.println("Elastic Free Energy: " + elasticF);
	 }
*/
	 totalF = elasticF + mixF;

	 energy1[i] = pairEnergySum + totalF;
         totalPairEnergy += pairEnergySum;
         totalFreeEnergy += totalF;
	 totalVirial += virialSum;
	 
      }

// loop over species 2

      // consider interactions with other particles
      for(i = 0; i < N2; ++i) { 
         pairEnergySum = 0;
         virialSum = 0;
         for(j = 0; j < N2; ++j) { // loop over 22 pairs
            if(j != i) {
               dx = PBC.separation(x2[i]-x2[j], side);
               dy = PBC.separation(y2[i]-y2[j], side);
               dz = PBC.separation(z2[i]-z2[j], side);
               r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
	       sigma = (a2[i] + a2[j])*dryRadiusRatio; // sum of particle radii
	       if (r2 < sigma*sigma){
                  r = Math.sqrt(r2);
                  // Hertz pair potential amplitude (scaled by a factor of scale*Young)
                  B22 = scale*Young*nChains2*sigma*sigma*Math.sqrt(a2[i]*a2[j])/(Math.pow(a2[i], 3.)+Math.pow(a2[j], 3.));
		  derivPart = B22*Math.pow(1-r/sigma, 1.5);
		  HertzEnergy = derivPart*(1-r/sigma);
		  pairEnergySum += HertzEnergy;
                  pairEnergy22[i][j] = HertzEnergy; // pair energy of particles i and j
/*
		  fOverR = (2.5/sigma)*derivPart/r;
		  fx = fOverR*dx; // force in x-direction
		  fy = fOverR*dy; // force in y-direction
		  fz = fOverR*dz; // force in z-direction
		  virialSum += dx*fx+dy*fy+dz*fz;
*/
	       }
            }
         }

         for(j = 0; j < N1; ++j) { // loop over 21 pairs
            dx = PBC.separation(x2[i]-x1[j], side);
            dy = PBC.separation(y2[i]-y1[j], side);
            dz = PBC.separation(z2[i]-z1[j], side);
            r2 = dx*dx+dy*dy+dz*dz; // particle separation squared
            R1 = a1[j]; // radius of particle 1
            R2 = a2[i]*dryRadiusRatio; // radius of particle 2
	    sigma = R1 + R2; // sum of particle radii
            B12 = scale*Young*sigma*sigma*Math.sqrt(R1*R2)*(nChains1*nChains2/(nChains1*Math.pow(R2, 3.)+nChains2*Math.pow(R1, 3.)));
	    if (r2 < sigma*sigma){
               r = Math.sqrt(r2);
	       derivPart = B12*Math.pow(1-r/sigma, 1.5);
	       HertzEnergy = derivPart*(1-r/sigma);
	       pairEnergySum += HertzEnergy;
               pairEnergy12[j][i] = HertzEnergy; // pair energy of particles i and j
/*
	       fOverR = (2.5/sigma)*derivPart/r;
	       fx = fOverR*dx; // force in x-direction
	       fy = fOverR*dy; // force in y-direction
	       fz = fOverR*dz; // force in z-direction
	       virialSum += dx*fx+dy*fy+dz*fz;
*/
	    }
         }

         // Flory-Rehner single-particle free energy (associated with swelling)
	 mixF = nMon2*((a2[i]*a2[i]*a2[i]-1)*Math.log(1-1/a2[i]/a2[i]/a2[i])+chi2*(1-1/a2[i]/a2[i]/a2[i]));
	 elasticF = 1.5*nChains2*(a2[i]*a2[i]-Math.log(a2[i])-1);
/*
	 if (a2[i] < 0){
	     System.out.println("Negative Swelling Ratio: " + a2[i]);
	     System.out.println("Elastic Free Energy: " + elasticF);
	 }
*/
	 totalF = elasticF + mixF;

	 energy2[i] = pairEnergySum + totalF;
         totalPairEnergy += pairEnergySum;
         totalFreeEnergy += totalF;
	 totalVirial += virialSum;

      }

      totalVirial *= 0.5; // correct for double counting
      totalPairEnergy *= 0.5;
      totalEnergy = totalPairEnergy + totalFreeEnergy;

      virialAccumulator += totalVirial; // running total
      energyAccumulator += totalEnergy;
      pairEnergyAccumulator += totalPairEnergy;
      freeEnergyAccumulator += totalFreeEnergy;
   }

   // mean energy per particle (kT units)
   public double meanEnergy() {
       return energyAccumulator/(N1+N2)/(double)(steps); // quantity <E>/N
       //return energyAccumulator/(N1+N2)/(double)(steps-delay); // quantity <E>/N
   }

   // mean pair energy per particle (kT units)
   public double meanPairEnergy() {
      return pairEnergyAccumulator/(N1+N2)/(double)(steps); // quantity <Epair>/N
   }

   public double meanFreeEnergy() { // mean free energy per particle (kT units)
      return freeEnergyAccumulator/(N1+N2)/(double)(steps); // quantity <F>/N
   }

   public double meanPressure() { // mean pressure (dimensionless)
      double meanVirial;
      meanVirial = virialAccumulator/(double)(steps);
      //meanVirial = virialAccumulator/(double)(steps-delay);
      return 1+(1./3.)*meanVirial/(N1+N2); // quantity PV/NkT
   }

   public void sizeDistributions() {
       int i, bin;
       for (i=0; i<N1; i++){
	   bin = (int) (a1[i]/sizeBinWidth);
	   sizeDist1[bin]++; // radius 1 [units of dry radius 1]
       }
       for (i=0; i<N2; i++){
	   bin = (int) (a2[i]*dryRadiusRatio/sizeBinWidth);
	   sizeDist2[bin]++; // radius 2 [units of dry radius 1]
       }
   }

   public double meanSwellingRatio1() {
       double sum = 0;
       for (int i=0;i<N1;i++){
           sum += a1[i];
        }
        meanSR1 = sum/N1;
       return meanSR1;
    }
    
   public double meanSwellingRatio2() {
       double sum = 0;
       for (int i=0;i<N2;i++){
           sum += a2[i];
        }
        meanSR2 = sum/N2;
       return meanSR2;
    }

   public double meanRadius1() { // mean radius 1 [units of dry radius 1]
       double sum;
       sum=0;
       if (steps > delay){
	   for (int i=0;i<maxSwellingRatio/sizeBinWidth;i++){
	       sum += (i*sizeBinWidth)*sizeDist1[i];
	   }
	   meanR1 = sum/N1;
       }
       else
	   meanR1 = 0;
       
       return meanR1;
   }

   public double meanRadius2() { // mean radius 2 [units of dry radius 1]
       double sum;
       sum=0;
       if (steps > delay){
	   for (int i=0;i<maxSwellingRatio/sizeBinWidth;i++){
	       sum += (i*sizeBinWidth)*sizeDist2[i];
	   }
	   meanR2 = sum/N2;
       }
       else
	   meanR2 = 0;
       
       return meanR2;
   }

    /* compute mean volume fraction after stopping */
    public void calculateVolumeFraction() { // instantaneous volume fraction 
        double microgelVol, overlapVol;
        double xij, yij, zij, r, r2, sigma, da2, amin;

        microgelVol = 0;
        for (int i=0; i<N1; i++) {
            microgelVol += (4./3.)*Math.PI*a1[i]*a1[i]*a1[i]; // sum volumes of spheres
            for (int j=i+1; j<N1; j++){
                xij = PBC.separation(x1[i]-x1[j], side);
                yij = PBC.separation(y1[i]-y1[j], side);
                zij = PBC.separation(z1[i]-z1[j], side);
                r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
                sigma = a1[i]+a1[j]; // sum of radii of two particles (units of dry radius)
                da2 = Math.pow(a1[i]-a1[j], 2); // difference of radii squared 
                if (r2 < sigma*sigma){ // particles are overlapping 
                    if (r2 < da2){ // one particle is entirely inside the other
                        amin = Math.min(a1[i],a1[j]); // minimum of radii
                        overlapVol = (4./3.)*Math.PI*amin*amin*amin; // volume of smaller particle
                    }
                    else { // one particle is NOT entirely inside the other
                        r = Math.sqrt(r2);
                        overlapVol = Math.PI*(sigma-r)*(sigma-r)*(r2+2*r*sigma-3*da2)/12./r; // volume of lens-shaped overlap region
                    }
                    microgelVol -= overlapVol; // subtract overlap volume so we don't double-count
                }
            }
            for (int j=1; j<N2; j++){
                xij = PBC.separation(x1[i]-x2[j], side);
                yij = PBC.separation(y1[i]-y2[j], side);
                zij = PBC.separation(z1[i]-z2[j], side);
                r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
                sigma = a1[i]+a2[j]; // sum of radii of two particles (units of dry radius)
                da2 = Math.pow(a1[i]-a2[j], 2); // difference of radii squared 
                if (r2 < sigma*sigma){ // particles are overlapping 
                    if (r2 < da2){ // one particle is entirely inside the other
                        amin = Math.min(a1[i],a2[j]); // minimum of radii
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

        for (int i=0; i<N2; i++) {
            microgelVol += (4./3.)*Math.PI*a2[i]*a2[i]*a2[i]; // sum volumes of spheres
            for (int j=i+1; j<N2; j++){
                xij = PBC.separation(x2[i]-x2[j], side);
                yij = PBC.separation(y2[i]-y2[j], side);
                zij = PBC.separation(z2[i]-z2[j], side);
                r2 = xij*xij + yij*yij + zij*zij; // particle separation squared
                sigma = a2[i]+a2[j]; // sum of radii of two particles (units of dry radius)
                da2 = Math.pow(a2[i]-a2[j], 2); // difference of radii squared 
                if (r2 < sigma*sigma){ // particles are overlapping 
                    if (r2 < da2){ // one particle is entirely inside the other
                        amin = Math.min(a2[i],a2[j]); // minimum of radii
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


/*
    public void calculateVolumeFraction() { // accumulator (NOT actual volume fraction until stopping)
	double R2, totalVol, microgelVol;       // at stopping, correct volFrac in writeData 

	microgelVol=0;
	totalVol=side*side*side;
	for (int i=0; i<N1; i++) {
	    microgelVol += (4./3.)*Math.PI*a1[i]*a1[i]*a1[i];
	}
	for (int i=0; i<N2; i++) {
            R2 = a2[i]*dryRadiusRatio;
	    microgelVol += (4./3.)*Math.PI*Math.pow(R2, 3.);
	}

	volFrac += microgelVol/totalVol;

    }
*/

/*
    public void snapshot(){
       if (steps>delay){
          if ((steps-delay) % snapshotInterval == 0){
	     sizeDistributions();                
	     rdf.update();
	     calculateVolumeFraction();
	    }
	}
   }
*/

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
