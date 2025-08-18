package org.opensourcephysics.sip.ch04;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;

/**
 * DDOApp solves and displays the time evolution of a driven, damped harmonic oscillator by stepping a SHO model.
 *
 * @author Alan Denton
 * @version 1.0
 */
public class DDOApp extends AbstractSimulation {
  PlotFrame plotFrame1 = new PlotFrame("Time", "Displacement", "Displacement vs. time");
  PlotFrame plotFrame2 = new PlotFrame("Time", "Energy", "Energy vs. time");
  PlotFrame plotFrame3 = new PlotFrame("Position", "Momentum", "Phase Space");
  DDO ddo;
  EulerRichardson odeSolver;
  DisplayFrame displayFrame = new DisplayFrame("Simple Harmonic Oscillator");

  /**
   * Constructs the DDOApp and intializes the display.
   */
  public DDOApp() {
    ddo = new DDO(0, 0, 0, 0, 0, 0);
    odeSolver = new EulerRichardson(ddo); // create numerical method
    displayFrame.addDrawable(ddo);
    displayFrame.setPreferredMinMax(-.5, .5, -.5, .5);
  }

  /*
   * Initializes the animation using the values in the control
   */
  public void initialize() {
    ddo.x0 = control.getDouble("x0");
    ddo.state[0] = ddo.x0;
    ddo.v0 = control.getDouble("v0");
    ddo.state[1] = ddo.v0;
    ddo.omega0 = control.getDouble("omega0");
    ddo.gamma = control.getDouble("gamma");
    ddo.omega = control.getDouble("omega");
    ddo.a0 = control.getDouble("A0");
    odeSolver.setStepSize(control.getDouble("dt"));
  }

  /**
   * Does an animation step
   */
  public void doStep() {
    plotFrame1.append(0, ddo.state[2], ddo.state[0]); // displacement vs time
    plotFrame2.append(0, ddo.state[2], ddo.energy);   // energy vs time
    plotFrame3.append(0, ddo.state[0], ddo.state[1]); // velocity vs position
    odeSolver.step();   // advance state by current step size
  }

  /**
   * Resets animation to a predefined state
   */
  public void reset() {
    control.setValue("x0", 0);     
    control.setValue("v0", 0);     
    control.setValue("omega0", 3);
    control.setValue("gamma", 0.5);
    control.setValue("omega", 2); 
    control.setValue("A0", 1);   
    control.setValue("dt", 0.01);     
    enableStepsPerDisplay(true);
/*
Note: To change number of steps per display, edit 2 occurrences of "stepsPerDisplay = 1" in ~/numerics/AbstractSimulation.java
*/
    odeSolver.setStepSize(0.01);
    ddo.state[0] = 0; // initial displacement
    ddo.state[1] = 0; // initial velocity
    ddo.state[2] = 0; // initial time
    initialize();
  }

  /**
   * Start Java application
   * @param args  command line parameters
   */
  public static void main(String[] args) {
    Control control = SimulationControl.createApp(new DDOApp());
  }
}
