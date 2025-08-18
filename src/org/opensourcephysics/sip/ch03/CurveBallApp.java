package org.opensourcephysics.sip.ch03;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;

/**
 * CurveBallApp solves and displays the time evolution of a spinning ball by stepping a projectile model.
 *
 * @author Alan Denton 09/18/08
 */
public class CurveBallApp extends AbstractSimulation {
  PlotFrame plotFrame = new PlotFrame("t", "y,z", "Position versus time");
  CurveBall ball = new CurveBall();
  PlotFrame animationFrame = new PlotFrame("x", "y,z", "Trajectory");

  public CurveBallApp() {
    animationFrame.addDrawable(ball);
    plotFrame.setXYColumnNames(0, "t", "x");
    plotFrame.setXYColumnNames(1, "t", "y");
    plotFrame.setXYColumnNames(2, "t", "z");
  }

  /**
   * Initializes the simulation.
   */
  public void initialize() {
    double dt = control.getDouble("dt");
    double x = control.getDouble("initial x");
    double y = 0;
    double z = control.getDouble("initial z");
    double v = control.getDouble("initial v");
    double theta = control.getDouble("initial angle (degrees)");
    theta = theta*Math.PI/180; // convert from radians to degrees;
    double vx = v*Math.cos(theta);  
    double vy = 0; // initial velocity in x-z plane
    double vz = v*Math.sin(theta);
    double wx = 0; // no spin about x-axis
    double wy = control.getDouble("initial wy");
    double wz = control.getDouble("initial wz");
    double cd = control.getDouble("drag coefficient cd"); 
    double cm = control.getDouble("Magnus coefficient cm"); 
    ball.setState(x, vx, wx, y, vy, wy, z, vz, wz, cd, cm);
    ball.setStepSize(dt);
//  double size = (vx*vx+vy*vy)/10; // estimate of size needed for display
//  animationFrame.setPreferredMinMax(-1, size, -1, size);
    double size = (vx*vx+vy*vy)/60; // estimate of size needed for display
    animationFrame.setPreferredMinMax(-1, size, -1, 0.1*size);
  }

  /**
   * Does a time step.
   */
  public void doStep() {
    plotFrame.append(1, ball.state[9], ball.state[3]); // y vs time data added
    plotFrame.append(2, ball.state[9], ball.state[6]); // z vs time data added
    animationFrame.append(0, ball.state[0], ball.state[6]); // x-z trajectory data added
    animationFrame.append(1, ball.state[0], ball.state[3]); // x-y trajectory data added
    ball.step(); // advance the state by one time step
  }

  /**
   * Resets the simulation.
   */
  public void reset() {
    control.setValue("initial x", 0);
    control.setValue("initial z", 1.8);
    control.setValue("initial v", 40);
    control.setValue("initial angle (degrees)", 0);
    control.setValue("initial wy", 0);
    control.setValue("initial wz", 200);
    control.setValue("dt", 0.01);
    control.setValue("drag coefficient cd", 0.006); 
    control.setValue("Magnus coefficient cm", 0.0004); 
    enableStepsPerDisplay(true);
  }

  /**
   * Starts the Java application.
   * @param args  command line parameters
   */
  public static void main(String[] args) {
    SimulationControl.createApp(new CurveBallApp());
  }
}
