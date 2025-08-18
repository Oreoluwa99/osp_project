package org.opensourcephysics.sip.ch03;
import org.opensourcephysics.controls.*;

/**
 * Charge3DApp solves and displays in a 3D box the time evolution of a 
 * charged particle in electric and magnetic fields (GTC, Problem 3.17)
 *
 * @author Alan Denton 02/10/10
 */
public class Charge3DApp extends AbstractSimulation {
  Charge3D charge = new Charge3D();

  /**
   * Initializes the simulation.
   */
  public void initialize() {
    double dt = control.getDouble("dt");
    double x = 0;
    double y = 0;
    double z = 0;
    double vx = control.getDouble("initial vx");
    double vy = control.getDouble("initial vy");
    double vz = control.getDouble("initial vz");
    double ex = control.getDouble("Ex");  // x-component of electric field 
    double ey = control.getDouble("Ey");  // y-component of electric field 
    double ez = control.getDouble("Ez");  // z-component of electric field 
    double bx = control.getDouble("Bx");  // x-component of magnetic field 
    double by = control.getDouble("By");  // y-component of magnetic field 
    double bz = control.getDouble("Bz");  // z-component of magnetic field 
    charge.setState(x, vx, y, vy, z, vz, ex, ey, ez, bx, by, bz);
    charge.setStepSize(dt);
  }

  /**
   * Does a time step.
   */
  public void doStep() {
    charge.step();  // advance the state by one time step
  }

  /**
   * Resets the simulation.
   */
  public void reset() {
    control.setValue("initial vx", 10);
    control.setValue("initial vy", 0);
    control.setValue("initial vz", 0);
    control.setValue("Ex", 0);
    control.setValue("Ey", 10);
    control.setValue("Ez", 0);
    control.setValue("Bx", 0);
    control.setValue("By", 0);
    control.setValue("Bz", 10);
    control.setValue("dt", 0.01);
    enableStepsPerDisplay(true);
  }

  /**
   * Starts the Java application.
   * @param args  command line parameters
   */
  public static void main(String[] args) {
    SimulationControl.createApp(new Charge3DApp());
  }
}
