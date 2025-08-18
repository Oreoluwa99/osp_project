package org.opensourcephysics.sip.Hertz;

import org.opensourcephysics.numerics.*;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Performs a Monte Carlo simulation of volume-charged ionic microgels interacting via a
 * Hertz-Yukawa effective pair potential and swelling according to the Flory-Rehner free energy.
 * 
 * @author Based on MC code ch15/LJParticles, modified for microgels by Matt Urich and Alan Denton.
 * @version 2.2 05.05.2020
 * 
 */

public class VolumeChargedMicrogels {
	
// lengths are in units of the dry microgel radius
// energies are in units of k_B*T (k_B is Boltzmann constant)
    public int N; // number of particles
    public int nx; // number of columns and rows
    public double side; // box side length
    public double d; // center-to-center distance between particles 
    public double x[], y[], z[], a[]; // x, y, z coordinates of particles and radius array
    public double energy[]; // energies of particles
    public double pairEnergy[][], newPairEnergy[][]; // pair energies of particles
    public double totalEnergy;
    public double energyAccumulator;
    public double totalVirial;
    public double virialAccumulator;
    public double tolerance, atolerance; // tolerances for trial displacements and radius changes
    public double B; // amplitude of Hertz potential
    public double nChains, nMon, chi;
    public double dryR, reservoirR; // dry radius, reservoir swollen radius
    public double meanR, meanSR; // mean radius, mean swelling ratio
    public double volFrac, dryVolFrac, reservoirVolFrac;
    public double sizeDist[], sizeBinWidth, maxSwellingRatio; // radius histogram, bin width, range 
    public double radBinWidth,deltaK;

    public double lBjerrum, gamma; // Bjerrum length, electrostatic coupling constant
    public double Z; // valence (counterions per microgel)
    public double csalt, Nsalt; // salt concentration and number of salt ion pairs per microgel
    public double Avogadro; // Avogadro's number
    public double k; // screening constant
    public int steps; // number of Monte Carlo steps
    public double delay, stop;
    public int snapshotDelay;
    public String fileExtension; // allows unique output file names 


   /**
   * Initialize the model.
   * 
   * @param configuration
   *  Initial lattice structure 
   */
   public void initialize(String configuration) {
      x = new double[N];
      y = new double[N];
      z = new double[N];
      a = new double[N];

      reservoirR = reservoirSwellingRatio(nMon, nChains, chi);
      reservoirVolFrac = dryVolFrac*Math.pow(reservoirR, 3.);

      for (int i = 0;i<N;i++){
	  a[i]=reservoirR;
      }
      energy = new double[N];
      pairEnergy = new double[N][N];
      newPairEnergy = new double[N][N];
      steps = 1;
      energyAccumulator = 0;
      virialAccumulator = 0;
      volFrac = 0;

      gamma = lBjerrum/dryR; // electrostatic coupling constant
      side = Math.cbrt(4*Math.PI*N/dryVolFrac/3.); // side length of simulation box

      Avogadro = 6.0221413e-07; // Avogadro's number [micromol/nm^3]
      // number of salt ion pairs per microgel
      Nsalt = csalt*Avogadro*Math.pow(side*dryR,3)/N;
      // Debye screening constant (counterions only; units of inverse dry radius)
      k = Math.sqrt(3*(Z+2*Nsalt)*dryVolFrac*gamma); 

      sizeDist = new double[(int) (maxSwellingRatio/sizeBinWidth)];

      // initialize positions
      if (configuration.toUpperCase().equals("CUBE")) {
         setCubePositions();
      } 
      if (configuration.toUpperCase().equals("FCC")) {
         setFCCPositions();
      } 
   }

   /**
   * Place particles on sites of a cubic lattice.
   */
   public void setCubePositions() {
      int ix, iy, iz;
      double dnx = Math.cbrt(N);
      d = side / dnx; // center-center distance between two particles
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
      calculateTotalEnergy(); // initial energy
   }

   /**
   * Place particles on sites of an fcc lattice.
   */
   public void setFCCPositions() {
      int ix, iy, iz;
      double dnx = Math.cbrt(N/4.);
      d = side / dnx; // lattice constant
      nx = (int) dnx;
      if (dnx - nx > 0.00001) {
         nx++; // N is not a perfect cube
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
      calculateTotalEnergy(); // initial energy
   }

   /**
   * Do a Monte Carlo simulation step.
   */
   public void step() {
      steps++;
      double dxtrial, dytrial, dztrial, atrial;
      double newEnergy, HertzEnergy;
      double dx, dy, dz, de, r, r2, sigma, aMean, R;
      double vbare, vind, voverlap, A, YukawaEnergy;
      double a3, mixF, elasticF, elecF, totalF, selfU, volU, K, ka;
      
      // trial move
      for (int i = 0; i < N; ++i) {
	  
         dxtrial = tolerance*2.*(Math.random()-0.5);
         dytrial = tolerance*2.*(Math.random()-0.5);
         dztrial = tolerance*2.*(Math.random()-0.5);
         atrial = atolerance*2.*(Math.random()-0.5);
	 
         x[i] = PBC.position(x[i]+dxtrial, side); // trial displacement
         y[i] = PBC.position(y[i]+dytrial, side);
         z[i] = PBC.position(z[i]+dztrial, side);
	 a[i] += atrial; // trial radius change
	 
/*
	 if (a[i]<0){
	     System.out.println("Negative radius: " + a[i]);
	 }
*/

         newEnergy = 0;
	 
         for (int j = 0; j < N; ++j) {
            if(j != i) {
               dx = PBC.separation(x[i]-x[j], side);
               dy = PBC.separation(y[i]-y[j], side);
               dz = PBC.separation(z[i]-z[j], side);
               r2 = dx*dx+dy*dy+dz*dz;
	       r = Math.sqrt(r2);
	       sigma = a[i] + a[j];
	       aMean = 0.5*sigma;
	       R = r/aMean;
	       ka = k*aMean;
	       if (r < sigma){
	           vbare = (gamma/aMean)*Z*Z*(6./5.-0.5*R*R+(3./16.)*Math.pow(R,3.)-Math.pow(R,5)/160);
                   vind = -0.5*(gamma/r)*Math.pow(3*Z/ka/ka,2)*(Math.pow(1+1/ka,2)
                      *Math.exp(-2*ka)*Math.sinh(k*r)+(1-1/(ka*ka))*(1-Math.exp(-k*r)
                      +Math.pow(k*r,2)/2+Math.pow(k*r,4)/24)-(2./3.)*ka*ka*(1-0.4*ka*ka)*R
                      -Math.pow(k,4)*Math.pow(r,3)*aMean/9-Math.pow(k,4)*Math.pow(r,6)/aMean/aMean/720);

                   voverlap = vbare + vind;

		   HertzEnergy = B*Math.pow(1-r/sigma,2.5);
		   newEnergy += voverlap + HertzEnergy; // pair energy
                   newPairEnergy[i][j] = voverlap + HertzEnergy;
	       }
	       else {
	           A = gamma*(3*Z*(Math.cosh(k*a[i])-Math.sinh(k*a[i])/(k*a[i]))/Math.pow(k*a[i],2))
                       *(3*Z*(Math.cosh(k*a[j])-Math.sinh(k*a[j])/(k*a[j]))/Math.pow(k*a[j],2));
	           YukawaEnergy = A*Math.exp(-k*r)/r;
	           newEnergy += YukawaEnergy; // pair energy
                   newPairEnergy[i][j] = YukawaEnergy;
	       }
	    }
         }
         a3 = Math.pow(a[i],3);
	 mixF = nMon*((a3-1)*Math.log(1-1/a3)+chi*(1-1/a3));
	 elasticF = 1.5*nChains*(a[i]*a[i]-Math.log(a[i])-1);
	 
         K = k*a[i]; // screening constant (units of inverse swollen radius)
         selfU = 0.6*Z*Z*gamma/a[i]; // self energy of a single ionic microgel
         // microion volume energy (neglecting constants independent of a[i])
         volU = -(3*Z*Z*gamma/a[i])*(0.2-0.5/K/K+(0.75/Math.pow(K,3))
                *(1-1/K/K+(1+2/K+1/K/K)*Math.exp(-2*K))); 
         elecF = selfU + volU; // electrostatic free energy

	 totalF = mixF + elasticF + elecF; // single-particle free energy
	 newEnergy += totalF;
	 
         de = newEnergy-energy[i];
	 
         if(Math.exp(-de) < Math.random() || a[i] < 1) {
            x[i] = PBC.position(x[i]-dxtrial, side); // reject move
            y[i] = PBC.position(y[i]-dytrial, side);
            z[i] = PBC.position(z[i]-dztrial, side);
	    a[i] -= atrial;
         }
         else { // accept move and update energies
            energy[i] = newEnergy;
            for (int j = 0; j < N; ++j) { // update energies of other particles 
               if(j != i) {
                  energy[j] += newPairEnergy[i][j]-pairEnergy[i][j];
                  pairEnergy[i][j] = newPairEnergy[i][j];
                  pairEnergy[j][i] = newPairEnergy[i][j];
               }
            }
         }
      }
      calculateTotalEnergy(); // new total energy
   }

   // total energy and pressure
   public void calculateTotalEnergy() {
      double dx, dy, dz, r, r2, sigma, aMean, R;
      double vbare, vind, voverlap, A, YukawaEnergy;
      double a3, mixF, elasticF, elecF, totalF, selfU, volU, volP, K, ka;
      double fOverR, fx, fy, fz;
      double fbare, find, foverlap;
      double derivPart, HertzEnergy;
      double pairEnergySum, energySum, virial1, virial2, virialSum;

      totalEnergy = 0;
      totalVirial = 0;

      for(int i = 0; i < N; ++i) { // sum over particles
         pairEnergySum = 0;
         energySum = 0;
         virialSum = 0;
         for(int j = 0; j < N; ++j) { // sum over pairs of particles
            if(j != i) {
               dx = PBC.separation(x[i]-x[j], side);
               dy = PBC.separation(y[i]-y[j], side);
               dz = PBC.separation(z[i]-z[j], side);
               r2 = dx*dx+dy*dy+dz*dz;
               r = Math.sqrt(r2);
	       sigma = a[i]+ a[j];
	       aMean = 0.5*sigma;
	       R = r/aMean;
	       ka = k*aMean;
	       if (r < sigma){
	           vbare = (gamma/aMean)*Z*Z*(6./5.-0.5*R*R+(3./16.)*Math.pow(R,3.)-Math.pow(R,5)/160);
                   vind = -0.5*(gamma/r)*Math.pow(3*Z/ka/ka,2)*(Math.pow(1+1/ka,2)
                      *Math.exp(-2*ka)*Math.sinh(k*r)+(1-1/(ka*ka))*(1-Math.exp(-k*r)
                      +Math.pow(k*r,2)/2+Math.pow(k*r,4)/24)-(2./3.)*ka*ka*(1-0.4*ka*ka)*R
                      -Math.pow(k,4)*Math.pow(r,3)*aMean/9-Math.pow(k,4)*Math.pow(r,6)/aMean/aMean/720);

                   voverlap = vbare + vind;

 	           fbare = (gamma/aMean/aMean)*Z*Z*(R-(9./16.)*R*R+Math.pow(R,4)/32.);

                   find = vind/r+(0.5*k*gamma/r)*Math.pow(3.*Z/ka/ka,2)
	                  *(Math.pow(1.+1./ka,2)*Math.exp(-2.*ka)*Math.cosh(k*r)
                          +(1.-1./ka/ka)*(Math.exp(-k*r)+k*r+Math.pow(k*r,3)/6.)
                          -(2./3.)*ka*(1.-0.4*ka*ka)-(1./3.)*Math.pow(k,3)*r*r*aMean
                          -(1./120.)*Math.pow(k,3)*Math.pow(r,5)/aMean/aMean);

	           foverlap = fbare + find;

	           derivPart = B*Math.pow(1-r/sigma,1.5);
	           HertzEnergy = derivPart*(1-r/sigma);
                   pairEnergy[i][j] = voverlap + HertzEnergy; // pair energy of particles i and j
	           pairEnergySum += voverlap + HertzEnergy; // total pair energy of particle i
	           fOverR = ((2.5/sigma)*derivPart+foverlap)/r;
	           fx = fOverR*dx; // force in x-direction
	           fy = fOverR*dy; // force in y-direction
	           fz = fOverR*dz; // force in z-direction
                   virial1 = (dx*fx+dy*fy+dz*fz)/3.; // standard virial contribution to pressure
                   // contribution from density-dependent effective electrostatic pair potential
                   virial2 = 0.5*Math.pow(3*Z/ka/ka,2)*(gamma/r)*(-2+R*R+3/ka/ka-0.5*Math.pow(k*r,2)
                             +(2./3.)*k*k*r*aMean+(1./24.)*k*k*Math.pow(r,4)/aMean/aMean
                             +(2+k*r/2.-3./ka/ka-0.5*R/ka)*Math.exp(-k*r)
                             -(4+ka+6/ka+3/ka/ka)*Math.exp(-2*ka)*Math.sinh(k*r)
                             +0.5*Math.pow(1+1/ka,2)*Math.exp(-2*ka)*k*r*Math.cosh(k*r));
                   virialSum += virial1 + virial2;
	       }
	       else {
	           A = gamma*(3*Z*(Math.cosh(k*a[i])-Math.sinh(k*a[i])/(k*a[i]))/Math.pow(k*a[i],2))
                       *(3*Z*(Math.cosh(k*a[j])-Math.sinh(k*a[j])/(k*a[j]))/Math.pow(k*a[j],2));
                   YukawaEnergy = A*Math.exp(-k*r)/r;
                   pairEnergy[i][j] = YukawaEnergy; // pair energy of particles i and j
                   pairEnergySum += YukawaEnergy; // pair energy
                   fOverR = (YukawaEnergy/r)*(k+1/r);
                   fx = fOverR*dx; // force in x-direction
                   fy = fOverR*dy; // force in y-direction
                   fz = fOverR*dz; // force in z-direction
                   virial1 = (dx*fx+dy*fy+dz*fz)/3.; // standard virial contribution to pressure
                   // contribution from density-dependent effective electrostatic pair potential
                   virial2 = (ka*ka*Math.sinh(ka)/(ka*Math.cosh(ka)-Math.sinh(ka))-3-k*r/2.)*YukawaEnergy;
                   virialSum += virial1 + virial2;
	       }
            }
	 } // end of j loop over particles

         virialSum *= 0.5; // correct for double-counting pairs

         a3 = Math.pow(a[i],3);
	 mixF = nMon*((a3-1)*Math.log(1-1/a3)+chi*(1-1/a3));
	 elasticF = 1.5*nChains*(a[i]*a[i]-Math.log(a[i])-1);
	 
         K = k*a[i]; // screening constant (units of inverse swollen radius)
         selfU = 0.6*Z*Z*gamma/a[i]; // self energy of a single ionic microgel
         // microion volume energy (neglecting constants independent of a[i])
         volU = -(3*Z*Z*gamma/a[i])*(0.2-0.5/K/K+(0.75/Math.pow(K,3))
                *(1-1/K/K+(1+2/K+1/K/K)*Math.exp(-2*K))); 
         elecF = selfU + volU; // electrostatic free energy

	 totalF = mixF + elasticF + elecF; // single-particle free energy

	 energy[i] = pairEnergySum + totalF;
         energySum += 0.5*pairEnergySum + totalF; // correcting for double-counting pairs

         // microion volume pressure 
         volP = Z + 2*Nsalt+ 1.5*Z*Z*(gamma/a[i])*(-1/K/K+(9./4.)*Math.pow(K,-3)-3.75*Math.pow(K,-5)
                +(1.5/K/K+5.25*Math.pow(K,-3)+7.5*Math.pow(K,-4)+3.75*Math.pow(K,-5))*Math.exp(-2*K));

         virialSum += volP; 

	 if (steps>delay){
	     totalEnergy += energySum; 
	     totalVirial += virialSum; 
	 }
	 
      } // end of i loop over particles

      energyAccumulator += totalEnergy; // running total
      virialAccumulator += totalVirial; // running total
   }

   // mean pair energy per particle
   public double meanEnergy() {
       //return energyAccumulator/N/(double)(steps); // quantity <E_pair>/N
       return energyAccumulator/N/(double)(steps-delay); // quantity <E_pair>/N
   }


   // mean pressure (no tail correction)
   public double meanPressure() {
      double meanVirial;
      //meanVirial = virialAccumulator/(double)(steps);
      meanVirial = virialAccumulator/(double)(steps-delay);
      return 1+meanVirial/N; // quantity PV/NkT
      //return 1+meanVirial/N-(100./60.)*2*Nsalt; // quantity (P-Pr)V/NkT
   }

   public void sizeDistributions() {
       int bin;
       for (int i=0;i<N;i++){
	   bin = (int) (a[i]/sizeBinWidth);
	   sizeDist[bin]++;
       }
	
   }
    
   public double meanSwellingRatio() {
       double sum = 0;
       for (int i=0;i<N;i++){
	   sum += a[i];
	}
	meanSR = sum/N;
       return meanSR;
    }

   public double meanRadius() {
       double sum = 0;
       for (int i=0;i<maxSwellingRatio/sizeBinWidth;i++){
           sum += (i*sizeBinWidth)*sizeDist[i];
       }
       meanR = sum/N;
       return meanR;
    }

    public void calculateVolumeFraction() { // accumulator for instantaneous volume fraction 
	double totalVol, microgelVol;       // after stopping, we correct volFrac in writeData method

	microgelVol=0;
	totalVol=side*side*side;
	for (int i=0;i<N;i++) {
	    microgelVol += (4./3.)*Math.PI*a[i]*a[i]*a[i];
	}
	volFrac += microgelVol/totalVol;
    }

/*
    public void snapshot(){
	if (steps>delay){
	    if ((steps-delay) % snapshotDelay == 0){
		sizeDistributions();                
		rdf.update();
		calculateVolumeFraction();
	    }
	}
    }
*/

    /* finds reservoir swelling ratio as root of Flory-Rehner pressure */
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
         * constructs the FloryRehnerPressure function with the given parameters
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
          //double pMixing = 3.*nMon*(x*x*Math.log(1.-1./x/x/x)+1./x+chi/Math.pow(x,-4.));
          //double pElastic = 1.5*nChains*(2.*x-1./x);
          double pMixing = nMon*(x*x*x*Math.log(1.-1./x/x/x)+1.+chi/x/x/x);
          double pElastic = nChains*(x*x-0.5);
          double pressure = pMixing + pElastic;
          return pressure;
        }

    }

}
