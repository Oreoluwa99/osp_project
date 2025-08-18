package org.opensourcephysics.sip.ch03;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;

/**
 * ChargeApp solves and displays the time evolution of a charged particle
 * in electric and magnetic fields (GTC, Problem 3.17)
 *
 * @author Alan Denton 02/10/10
 */
public class ChargeApp extends AbstractSimulation {
  PlotFrame frame = new PlotFrame("x", "y", "Charge Trajectory");
  Charge charge = new Charge();

  /** 
   * Constructs a ChargeApp
   */
  public ChargeApp() {
    frame.addDrawable(charge);
    frame.setPreferredMinMax(-10, 10, -10, 10);
    frame.setSquareAspect(true);
  }

  /**
   * Initializes the simulation.
   */
  public void initialize() {
    double dt = control.getDouble("dt");
    double x = control.getDouble("initial x");
    double y = control.getDouble("initial y");
    double z = 0;
    double v = control.getDouble("initial v");
    double theta = control.getDouble("initial angle (degrees)");
    theta = theta*Math.PI/180;  // convert from radians to degrees;
    double vx = v*Math.cos(theta);
    double vy = v*Math.sin(theta);
    double b = control.getDouble("Bz");  // magnetic field in z direction
    double e = control.getDouble("Ey");  // electric field in y direction
    charge.setState(x, vx, y, vy, b, e);
    charge.setStepSize(dt);
  }

  /**
   * Does a time step.
   */
  public void doStep() {
    charge.step();  // advance the state by one time step
    frame.setMessage("t = "+decimalFormat.format(charge.state[6]));
  }

  /**
   * Resets the simulation.
   */
  public void reset() {
    control.setValue("initial x", 0);
    control.setValue("initial y", 0);
    control.setValue("initial v", 40);
    control.setValue("initial angle (degrees)", 0);
    control.setValue("Bz", 10);
    control.setValue("Ey", 10);
    control.setValue("dt", 0.01);
    enableStepsPerDisplay(true);
  }

  /**
   * Starts the Java application.
   * @param args  command line parameters
   */
  public static void main(String[] args) {
    SimulationControl.createApp(new ChargeApp());
  }
}
